classdef OuterLQR < handle
%======================================================================
% OUTERLQR
%======================================================================
% Slow rigid-body controller for coupled flexible-wing aircraft.
%
% Rigid-body state:
%   x_rb = [u; v; w; phi; theta; psi; p; q; r]
%
% Outer-loop input:
%   u_rb = [delta_e; delta_a; delta_r; thrust]
%
% Control law:
%   u_rb = u_trim - K*(x_rb - x_ref)
%
% The command is then allocated into:
%   1) uFlexOuter : low-frequency shared wing-surface vector
%   2) rbCmd      : elevator/rudder/thrust command for rigid-body plant
%
% This object does NOT replace nMHE/nMPC. It wraps around them as the
% slower rigid-body loop.
%======================================================================

    properties
        cfg
        trim

        Ts double
        useDiscrete logical = true

        nx double = 9
        nuRB double = 4

        A double
        B double
        Ad double
        Bd double

        Q double
        R double

        K double

        xTrim double
        uTrimRB double

        nSurf double
        varPer double
        nuFlex double

        aileronSurfaceMap double
        elevatorSurfaceMap double

        useAileronOnWingSurfaces logical = true
        useElevatorOnWingSurfaces logical = false

        structuredOuter struct = struct()
        sharedWingDaisyChain struct = struct()
        simplePdRigidPath struct = struct()
        scheduledEquilibriumFeedforward struct = struct()

        rbParams struct

        eigCL
        controllabilityRank double
        debug logical = true

        last struct
    end

    methods
        %--------------------------------------------------------------
        function obj = OuterLQR(cfg, trim)
        % Constructor. Builds/linearizes rigid-body model and solves LQR.
        %--------------------------------------------------------------
            obj.cfg  = cfg;
            obj.trim = trim;

            if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'outerTs')
                obj.Ts = cfg.ctrl.outerTs;
            elseif isfield(cfg,'ctrl') && isfield(cfg.ctrl,'outerEvery')
                obj.Ts = cfg.ctrl.outerEvery * cfg.ctrl.Ts;
            else
                obj.Ts = 4 * cfg.ctrl.Ts;
            end

            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'useDiscrete')
                obj.useDiscrete = logical(cfg.outerLQR.useDiscrete);
            else
                obj.useDiscrete = true;
            end
            % Full rigid-body LQR model dimensions.
            % x_rb = [u; v; w; phi; theta; psi; p; q; r]
            % u_rb = [delta_e; delta_a; delta_r; thrust]
            obj.nx   = 9;
            obj.nuRB = 4;
            obj.nSurf = cfg.ctrl.n_surf;
            obj.varPer = cfg.ctrl.var_per;
            obj.nuFlex = obj.nSurf * obj.varPer;

            obj.rbParams = obj.getRigidParams(cfg);

            obj.xTrim   = obj.buildTrimState(cfg, trim);
            obj.uTrimRB = obj.buildTrimInput(cfg, trim);

            obj.aileronSurfaceMap  = obj.buildSurfaceMap(cfg, 'aileronSurfaceMap',  [1; -1]);
            obj.elevatorSurfaceMap = obj.buildSurfaceMap(cfg, 'elevatorSurfaceMap', [1;  1]);

            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'useAileronOnWingSurfaces')
                obj.useAileronOnWingSurfaces = logical(cfg.outerLQR.useAileronOnWingSurfaces);
            end

            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'useElevatorOnWingSurfaces')
                obj.useElevatorOnWingSurfaces = logical(cfg.outerLQR.useElevatorOnWingSurfaces);
            end
            obj.structuredOuter = obj.loadStructuredOuterConfig(cfg);
            obj.sharedWingDaisyChain = obj.loadSharedWingDaisyChainConfig(cfg);
            obj.simplePdRigidPath = obj.loadSimplePdRigidPathConfig(cfg);
            obj.scheduledEquilibriumFeedforward = ...
                obj.loadScheduledEquilibriumFeedforwardConfig(cfg);

            [obj.A,obj.B] = obj.buildLinearModel(cfg);
            % ------------------------------------------------------------------
            % Validate rigid-body linear model before attempting LQR.
            % ------------------------------------------------------------------
            if isempty(obj.A) || isempty(obj.B)
                error(['OuterLQR:EmptyLinearModel\n', ...
                       'OuterLQR received empty A/B matrices. Either provide valid ', ...
                       'cfg.outerLQR.A and cfg.outerLQR.B, or allow buildLinearModel() ', ...
                       'to finite-difference the fallback rigid-body model.']);
            end
            
            if ~isequal(size(obj.A),[obj.nx,obj.nx])
                error('OuterLQR:BadA', ...
                      'A must be %d x %d. Got %d x %d.', ...
                      obj.nx,obj.nx,size(obj.A,1),size(obj.A,2));
            end
            
            if ~isequal(size(obj.B),[obj.nx,obj.nuRB])
                error('OuterLQR:BadB', ...
                      'B must be %d x %d. Got %d x %d.', ...
                      obj.nx,obj.nuRB,size(obj.B,1),size(obj.B,2));
            end
            [obj.Q,obj.R] = obj.buildWeights(cfg);

            obj.solveGain();

            obj.last = struct();

            if obj.debug
                obj.printDesignReport();
            end
        end

        %--------------------------------------------------------------
        function [uFlexOuter, rbCmd, info] = computeControl(obj, rbMeas, ref, t, realizedWingState)
        % Compute low-frequency outer-loop command.
        %
        % Inputs:
        %   rbMeas : struct with fields:
        %       v_B     3x1
        %       euler   3x1
        %       omega_B 3x1
        %
        %   ref : optional struct with any of:
        %       v_B_ref, euler_ref, omega_B_ref
        %
        % Outputs:
        %   uFlexOuter : nuFlex x 1 command for shared wing surfaces
        %   rbCmd      : struct with elevator/rudder/thrust/aileron
        %   info       : diagnostics
        %--------------------------------------------------------------
            if nargin < 4
                t = nan;
            end

            if nargin < 3 || isempty(ref)
                ref = struct();
            end

            if nargin < 5
                realizedWingState = struct();
            end

            x = obj.rbToState(rbMeas);
            xRef = obj.refToState(ref);

            dx = obj.wrapStateError(x - xRef);

            sharedWingInfo = struct('enabled',false,'symmetricIncrement',0);
            simplePdInfo = struct('enabled',false,'thetaGain',0,'qGain',0);
            speedGuidanceInfo = struct('enabled',false);
            equilibriumInfo = struct('enabled',false,'owner',"fixed_initial_trim", ...
                'uTrimRB',obj.uTrimRB(:),'uTrimWing',zeros(obj.nuFlex,1));
            if obj.simplePdRigidPath.enable
                [uFlexOuter,uRB,du,simplePdInfo] = obj.computeSimplePdCommand(dx);
                [~,rbCmd] = obj.allocateCommand(uRB);
            elseif obj.sharedWingDaisyChain.enable
                [uFlexOuter,uRB,du,sharedWingInfo] = ...
                    obj.computeSharedWingDaisyChainCommand(dx,realizedWingState);
                [~,rbCmd] = obj.allocateCommand(uRB);
            else
                [uTrimActive,uTrimWing,equilibriumInfo] = ...
                    obj.resolveScheduledEquilibriumFeedforward(realizedWingState);
                if equilibriumInfo.enabled && ...
                        obj.scheduledEquilibriumFeedforward.speedGuidance.enable
                    [du,speedGuidanceInfo] = obj.computeScheduledSpeedGuidance( ...
                        dx,x,xRef,ref);
                else
                    du = -obj.K * dx;
                end
                uRB = obj.saturateOuterRBInput(uTrimActive + du,uTrimActive);
                [uFlexOuter, rbCmd] = obj.allocateCommand(uRB);
                if equilibriumInfo.enabled
                    assert(norm(uFlexOuter,inf)<=1e-14, ...
                        'OuterLQR:ScheduledEquilibriumWingConflict', ...
                        ['The scheduled-equilibrium discriminator cannot ', ...
                         'coexist with a second direct wing allocation.']);
                    uFlexOuter = uTrimWing;
                end
            end

            info = struct();
            info.t = t;
            info.x = x;
            info.xRef = xRef;
            info.dx = dx;
            info.du = du;
            info.uRB = uRB;
            info.uFlexOuter = uFlexOuter;
            info.rbCmd = rbCmd;
            info.sharedWingDaisyChain = sharedWingInfo;
            info.simplePdRigidPath = simplePdInfo;
            info.scheduledEquilibriumFeedforward = equilibriumInfo;
            info.speedGuidance = speedGuidanceInfo;

            obj.last = info;
        end

        %--------------------------------------------------------------
        function [uFlexOuter, info] = computeStructuredWingIncrement(obj, xhat, rbMeas, t)
        %COMPUTESTRUCTUREDWINGINCREMENT Evaluate the approved audit-only map.
        %
        % xhat is the flexible nMHE state. rbMeas supplies body velocity,
        % Euler angles, and body rates; inertial position is intentionally
        % excluded. Invalid contract data fails closed to zero wing command.
        %--------------------------------------------------------------
            if nargin < 4
                t = nan;
            end

            uFlexOuter = zeros(obj.nuFlex,1);
            info = struct('t',t,'enabled',false,'accepted',false, ...
                'status',"DISABLED",'symmetricWingIncrement',0, ...
                'descriptorState',[],'outerOutput',[]);
            if ~isfield(obj.structuredOuter,'enable') || ...
                    ~obj.structuredOuter.enable
                return
            end
            info.enabled = true;

            try
                design = obj.structuredOuter.design;
                xhat = xhat(:);
                rbState = obj.rbToState(rbMeas);
                descriptorTrim = [obj.trim.states(:); obj.xTrim(:)];
                xDescriptor = [xhat; rbState(:)];
                assert(numel(xhat) == design.flexibleStateCount, ...
                    'OuterLQR:StructuredFlexibleState', ...
                    'The nMHE state length does not match the certified map.');
                assert(numel(xDescriptor) == design.runtimeStateCount && ...
                    numel(descriptorTrim) == design.runtimeStateCount, ...
                    'OuterLQR:StructuredState', ...
                    'The hybrid descriptor state has an unexpected length.');
                assert(all(isfinite(xDescriptor)) && all(isfinite(descriptorTrim)), ...
                    'OuterLQR:StructuredFinite', ...
                    'The hybrid descriptor state must be finite.');

                outputMap = design.outputMap;
                outputScales = design.outerOutputScales(:);
                gain = design.symmetricWingGain(:).';
                inputScales = design.outerInputScales(:);
                assert(isequal(size(outputMap), ...
                    [numel(outputScales), design.runtimeStateCount]), ...
                    'OuterLQR:StructuredOutputMap', ...
                    'The certified physical-output map has an unexpected size.');
                assert(numel(gain) == numel(outputScales) && ...
                    numel(inputScales) >= 1 && all(outputScales > 0), ...
                    'OuterLQR:StructuredGain', ...
                    'The certified structured-controller scaling is invalid.');

                yOuter = outputMap * (xDescriptor - descriptorTrim);
                uSymmetric = -inputScales(1) * gain * (yOuter ./ outputScales);
                assert(isfinite(uSymmetric), 'OuterLQR:StructuredCommand', ...
                    'The structured wing increment is nonfinite.');
                surfaceMap = design.symmetricWingSurfaceMap(:);
                assert(numel(surfaceMap) == obj.nSurf && ...
                    all(isfinite(surfaceMap)), 'OuterLQR:StructuredSurfaceMap', ...
                    'The symmetric wing-surface map is incompatible with cfg.ctrl.n_surf.');

                deflection = surfaceMap * uSymmetric;
                if obj.varPer == 1
                    uFlexOuter = deflection;
                elseif obj.varPer == 2
                    uFlexOuter = [deflection; zeros(obj.nSurf,1)];
                else
                    error('OuterLQR:StructuredVarPer', ...
                        'Unsupported cfg.ctrl.var_per = %d.', obj.varPer);
                end
                assert(all(isfinite(uFlexOuter)), 'OuterLQR:StructuredCommand', ...
                    'The structured wing command is nonfinite.');
                info.accepted = true;
                info.status = "ACCEPTED";
                info.symmetricWingIncrement = uSymmetric;
                info.descriptorState = xDescriptor;
                info.outerOutput = yOuter;
            catch err
                uFlexOuter = zeros(obj.nuFlex,1);
                info.status = "FAIL_CLOSED_" + string(err.identifier);
            end
        end
    end

    methods (Access = private)
        %--------------------------------------------------------------
        function request = loadScheduledEquilibriumFeedforwardConfig(~,cfg)
        %LOADSCHEDULEDEQUILIBRIUMFEEDFORWARDCONFIG Validate Case-B audit path.
            request = struct('enable',false);
            if ~isfield(cfg,'outerLQR') || ...
                    ~isfield(cfg.outerLQR,'scheduledEquilibriumFeedforward')
                return
            end
            request = cfg.outerLQR.scheduledEquilibriumFeedforward;
            if ~isfield(request,'enable') || ~logical(request.enable)
                request = struct('enable',false);
                return
            end
            assert(isfield(request,'auditOnly') && ...
                logical(request.auditOnly) && isfield(request,'source') && ...
                string(request.source)== ...
                "phase18c-v17a-caseb-scheduled-virtual-longitudinal-guidance-allocation-audit-v1", ...
                'OuterLQR:ScheduledEquilibriumFeedforwardScope', ...
                ['Scheduled equilibrium feedforward must remain bound to ', ...
                 'the approved default-inactive Case-B audit.']);
            if ~isfield(request,'speedGuidance')
                request.speedGuidance = struct('enable',false);
            end
            guidance = request.speedGuidance;
            if ~isfield(guidance,'enable') || ~logical(guidance.enable)
                request.speedGuidance = struct('enable',false);
                return
            end
            required = {'timeConstantSeconds','maximumAccelerationMps2'};
            assert(isstruct(guidance) && all(isfield(guidance,required)) && ...
                isscalar(guidance.timeConstantSeconds) && ...
                isfinite(guidance.timeConstantSeconds) && ...
                guidance.timeConstantSeconds>0 && ...
                isscalar(guidance.maximumAccelerationMps2) && ...
                isfinite(guidance.maximumAccelerationMps2) && ...
                guidance.maximumAccelerationMps2>0, ...
                'OuterLQR:ScheduledSpeedGuidanceConfig', ...
                'The Case-B speed-guidance configuration is invalid.');
            if isfield(guidance,'finiteHorizonAccelerationPerNewton')
                assert(isscalar(guidance.finiteHorizonAccelerationPerNewton) && ...
                    isfinite(guidance.finiteHorizonAccelerationPerNewton) && ...
                    guidance.finiteHorizonAccelerationPerNewton>0 && ...
                    isfield(guidance,'finiteHorizonSeconds') && ...
                    isscalar(guidance.finiteHorizonSeconds) && ...
                    isfinite(guidance.finiteHorizonSeconds) && ...
                    guidance.finiteHorizonSeconds>0 && ...
                    isfield(guidance,'effectivenessSource') && ...
                    strlength(string(guidance.effectivenessSource))>0, ...
                    'OuterLQR:ScheduledSpeedGuidanceEffectiveness', ...
                    ['The finite-horizon thrust effectiveness must be ', ...
                     'positive, finite, and source identified.']);
            end
        end

        %--------------------------------------------------------------
        function [du,info] = computeScheduledSpeedGuidance(obj,dx,x,xRef,ref)
        %COMPUTESCHEDULEDSPEEDGUIDANCE Separate axial and pitch demands.
            guidance = obj.scheduledEquilibriumFeedforward.speedGuidance;
            assert(isstruct(ref) && ...
                isfield(ref,'speedReferenceRateMps2') && ...
                isscalar(ref.speedReferenceRateMps2) && ...
                isfinite(ref.speedReferenceRateMps2), ...
                'OuterLQR:ScheduledSpeedGuidanceReference', ...
                'The speed-guidance path requires a finite reference rate.');
            referenceVelocity = xRef(1:3);
            measuredVelocity = x(1:3);
            referenceSpeed = norm(referenceVelocity);
            measuredSpeed = norm(measuredVelocity);
            assert(referenceSpeed>0 && isfinite(measuredSpeed), ...
                'OuterLQR:ScheduledSpeedGuidanceState', ...
                'The speed-guidance velocity state is invalid.');

            tangent = referenceVelocity/referenceSpeed;
            velocityError = measuredVelocity-referenceVelocity;
            axialError = tangent.'*velocityError;
            feedbackError = dx;
            feedbackError(1:3) = velocityError-tangent*axialError;
            du = -obj.K*feedbackError;

            speedError = referenceSpeed-measuredSpeed;
            accelerationUnbounded = double(ref.speedReferenceRateMps2) + ...
                speedError/double(guidance.timeConstantSeconds);
            accelerationDemand = min(max(accelerationUnbounded, ...
                -double(guidance.maximumAccelerationMps2)), ...
                double(guidance.maximumAccelerationMps2));
            effectiveness = 1/double(obj.rbParams.mass);
            effectivenessOwner = "rigid_mass_inverse";
            effectivenessHorizon = 0;
            if isfield(guidance,'finiteHorizonAccelerationPerNewton')
                effectiveness = double( ...
                    guidance.finiteHorizonAccelerationPerNewton);
                effectivenessOwner = string(guidance.effectivenessSource);
                effectivenessHorizon = double(guidance.finiteHorizonSeconds);
            end
            thrustIncrement = accelerationDemand/effectiveness;
            du(4) = thrustIncrement;
            info = struct('enabled',true, ...
                'owner',"reference_acceleration_plus_axial_speed_error", ...
                'referenceSpeedMps',referenceSpeed, ...
                'measuredSpeedMps',measuredSpeed, ...
                'speedErrorMps',speedError, ...
                'referenceAccelerationMps2', ...
                    double(ref.speedReferenceRateMps2), ...
                'accelerationUnboundedMps2',accelerationUnbounded, ...
                'accelerationDemandMps2',accelerationDemand, ...
                'thrustAccelerationEffectivenessMps2PerNewton', ...
                    effectiveness, ...
                'thrustEffectivenessOwner',effectivenessOwner, ...
                'thrustEffectivenessHorizonSeconds', ...
                    effectivenessHorizon, ...
                'thrustIncrementNewton',thrustIncrement, ...
                'axialVelocityErrorMps',axialError, ...
                'axialErrorRemovedFromPitchFeedback',true);
        end

        %--------------------------------------------------------------
        function [uTrimRB,uTrimWing,info] = ...
                resolveScheduledEquilibriumFeedforward(obj,runtimeState)
        %RESOLVESCHEDULEDEQUILIBRIUMFEEDFORWARD Read the atomic packet trim.
            uTrimRB = obj.uTrimRB(:);
            uTrimWing = zeros(obj.nuFlex,1);
            info = struct('enabled',false,'owner',"fixed_initial_trim", ...
                'uTrimRB',uTrimRB,'uTrimWing',uTrimWing);
            if ~obj.scheduledEquilibriumFeedforward.enable
                return
            end
            assert(isstruct(runtimeState) && ...
                isfield(runtimeState,'scheduledQueryTrim'), ...
                'OuterLQR:ScheduledEquilibriumFeedforwardMissing', ...
                'The active scheduled packet did not supply its query trim.');
            queryTrim = runtimeState.scheduledQueryTrim;
            required = {'wing','elevator','thrust'};
            assert(isstruct(queryTrim) && all(isfield(queryTrim,required)), ...
                'OuterLQR:ScheduledEquilibriumFeedforwardContract', ...
                'The scheduled query-trim contract is incomplete.');
            commandTrim = queryTrim;
            owner = "active_scheduled_reciprocal_packet_query_trim";
            if obj.scheduledEquilibriumFeedforward.speedGuidance.enable
                assert(isfield(runtimeState,'referenceQueryTrim'), ...
                    'OuterLQR:ScheduledSpeedGuidanceTrim', ...
                    ['The virtual-speed guidance path requires the ', ...
                     'source-interpolated reference equilibrium.']);
                commandTrim = runtimeState.referenceQueryTrim;
                owner = "reference_scheduled_reciprocal_packet_query_trim";
            end
            wing = double(commandTrim.wing(:));
            elevator = double(commandTrim.elevator);
            thrust = double(commandTrim.thrust);
            assert(numel(wing)==obj.nuFlex && isscalar(elevator) && ...
                isscalar(thrust) && thrust>=0 && ...
                all(isfinite([wing;elevator;thrust])), ...
                'OuterLQR:ScheduledEquilibriumFeedforwardValue', ...
                'The scheduled query trim is nonfinite or violates T >= 0.');
            uTrimRB = [elevator;0;0;thrust];
            uTrimWing = wing;
            info = struct('enabled',true, ...
                'owner',owner,'uTrimRB',uTrimRB,'uTrimWing',uTrimWing, ...
                'activeScheduledQueryTrim',queryTrim);
        end

        %--------------------------------------------------------------
        function simplePd = loadSimplePdRigidPathConfig(~,cfg)
        %LOADSIMPLEPDRIGIDPATHCONFIG Validate the disabled-by-default audit.
            simplePd = struct('enable',false,'auditOnly',true);
            if ~isfield(cfg,'outerLQR') || ...
                    ~isfield(cfg.outerLQR,'simplePdRigidPath')
                return
            end
            simplePd = cfg.outerLQR.simplePdRigidPath;
            if ~isfield(simplePd,'enable') || ~simplePd.enable
                simplePd.enable = false;
                return
            end
            assert(isfield(simplePd,'auditOnly') && simplePd.auditOnly && ...
                isfield(simplePd,'source') && string(simplePd.source)== ...
                "phase18c-v17a-casea-simple-pd-rigid-path-discriminator-v1", ...
                'OuterLQR:SimplePdAudit', ...
                'The simple PD path must remain bound to its audit contract.');
        end

        %--------------------------------------------------------------
        function [uFlexOuter,uRB,du,info] = computeSimplePdCommand(obj,dx)
        %COMPUTESIMPLEPDCOMMAND LQR-row-derived pitch/rate diagnostic.
            thetaIndex = 5;
            qIndex = 8;
            elevatorIndex = 1;
            thetaGain = obj.K(elevatorIndex,thetaIndex);
            qGain = obj.K(elevatorIndex,qIndex);
            assert(all(isfinite([thetaGain,qGain])), ...
                'OuterLQR:SimplePdGain', ...
                'The source-derived PD gains must be finite.');
            du = zeros(obj.nuRB,1);
            du(elevatorIndex) = -thetaGain*dx(thetaIndex) - qGain*dx(qIndex);
            uRB = obj.saturateOuterRBInput(obj.uTrimRB + du);
            uFlexOuter = zeros(obj.nuFlex,1);
            info = struct('enabled',true,'thetaGain',thetaGain,'qGain',qGain, ...
                'thetaError',dx(thetaIndex),'qError',dx(qIndex), ...
                'elevatorIncrement',du(elevatorIndex));
        end

        %--------------------------------------------------------------
        function structured = loadStructuredOuterConfig(~, cfg)
        %LOADSTRUCTUREDOUTERCONFIG Validate the disabled-by-default audit map.
        %--------------------------------------------------------------
            structured = struct('enable',false);
            if ~isfield(cfg,'outerLQR') || ...
                    ~isfield(cfg.outerLQR,'structuredOuter')
                return
            end
            structured = cfg.outerLQR.structuredOuter;
            if ~isfield(structured,'enable') || ~structured.enable
                structured.enable = false;
                return
            end
            required = {'auditOnly','design','artifactSha256', ...
                'expectedArtifactSha256'};
            for field = required
                assert(isfield(structured,field{1}), ...
                    'OuterLQR:StructuredConfig', ...
                    'The structured outer configuration is incomplete.');
            end
            assert(structured.auditOnly && ...
                string(structured.artifactSha256) == ...
                string(structured.expectedArtifactSha256), ...
                'OuterLQR:StructuredArtifactHash', ...
                'The structured outer artifact is not the approved hash.');
        end

        %--------------------------------------------------------------
        function chain = loadSharedWingDaisyChainConfig(~,cfg)
        %LOADSHAREDWINGDAISYCHAINCONFIG Validate an exact audit-only design.
            chain = struct('enable',false);
            if ~isfield(cfg,'outerLQR') || ...
                    ~isfield(cfg.outerLQR,'sharedWingDaisyChain')
                return
            end
            chain = cfg.outerLQR.sharedWingDaisyChain;
            if ~isfield(chain,'enable') || ~logical(chain.enable)
                chain = struct('enable',false);
                return
            end
            required = {'auditOnly','design','artifactSha256','sourceId'};
            for field = required
                assert(isfield(chain,field{1}), ...
                    'OuterLQR:SharedWingDaisyChainConfig', ...
                    'The shared-wing daisy-chain configuration lacks %s.',field{1});
            end
            assert(logical(chain.auditOnly), ...
                'OuterLQR:SharedWingDaisyChainScope', ...
                'The shared-wing daisy-chain configuration must be audit-only.');
            design = chain.design;
            requiredDesign = {'gain','rigidStateIndices', ...
                'wingIncrementLimitRadians','elevatorIncrementLimitRadians', ...
                'thrustIncrementLimitNewtons','sampleTimeSeconds'};
            for field = requiredDesign
                assert(isfield(design,field{1}), ...
                    'OuterLQR:SharedWingDaisyChainDesign', ...
                    'The shared-wing daisy-chain design lacks %s.',field{1});
            end
            realized = isfield(design,'designType') && ...
                string(design.designType)=="realized_multirate_fusion_actuator";
            expectedGainSize = [3,4];
            if realized, expectedGainSize = [3,6]; end
            assert(isequal(size(design.gain),expectedGainSize) && ...
                numel(design.rigidStateIndices)==4 && ...
                all(isfinite(design.gain),'all') && ...
                all(isfinite(design.rigidStateIndices)) && ...
                isfinite(design.wingIncrementLimitRadians) && ...
                design.wingIncrementLimitRadians>0 && ...
                abs(design.sampleTimeSeconds-cfg.ctrl.outerTs)<= ...
                    10*eps(design.sampleTimeSeconds), ...
                'OuterLQR:SharedWingDaisyChainDesign', ...
                'The shared-wing daisy-chain design is incompatible with runtime.');
            if realized
                assert(isfield(design,'stateOrder') && ...
                    isequal(string(design.stateOrder(:)).', ...
                    ["u","w","theta","q","outerLpSym","actuatorDeltaSym"]), ...
                    'OuterLQR:SharedWingRealizedStateOrder', ...
                    'The realised multirate state order is incompatible with OuterLQR.');
            end
        end

        %--------------------------------------------------------------
        function [uFlex,uRB,duRigid,info] = ...
                computeSharedWingDaisyChainCommand(obj,dx,realizedWingState)
        %COMPUTESHAREDWINGDAISYCHAINCOMMAND Evaluate coordinated outer inputs.
            design = obj.sharedWingDaisyChain.design;
            stateIndices = double(design.rigidStateIndices(:));
            assert(isequal(stateIndices,[1;3;5;8]), ...
                'OuterLQR:SharedWingDaisyChainStateOrder', ...
                'The shared-wing design state order is incompatible with OuterLQR.');
            if nargin < 3
                realizedWingState = struct();
            end
            realized = isfield(design,'designType') && ...
                string(design.designType)=="realized_multirate_fusion_actuator";
            if realized
                required = {'fusionOuterLP','actuatorDelta'};
                for field = required
                    assert(isfield(realizedWingState,field{1}), ...
                        'OuterLQR:SharedWingRealizedState', ...
                        ['The realised multirate shared-wing design requires ', ...
                         'the current %s state.'],field{1});
                end
                outerLp = realizedWingState.fusionOuterLP(:);
                actuatorDelta = realizedWingState.actuatorDelta(:);
                assert(~isempty(outerLp) && ~isempty(actuatorDelta) && ...
                    isfinite(outerLp(1)) && isfinite(actuatorDelta(1)), ...
                    'OuterLQR:SharedWingRealizedState', ...
                    'The realised multirate shared-wing state must be finite.');
                commandState = [dx(stateIndices);outerLp(1);actuatorDelta(1)];
            else
                commandState = dx(stateIndices);
            end
            virtualCommand = -double(design.gain) * commandState;
            wingLimit = double(design.wingIncrementLimitRadians);
            wingIncrement = min(max(virtualCommand(1),-wingLimit),wingLimit);
            uRB = obj.uTrimRB;
            uRB(1) = uRB(1) + virtualCommand(2);
            uRB(4) = uRB(4) + virtualCommand(3);
            uRB = obj.saturateOuterRBInput(uRB);
            duRigid = uRB-obj.uTrimRB;
            if obj.varPer==1
                uFlex = wingIncrement*ones(obj.nSurf,1);
            elseif obj.varPer==2
                uFlex = [wingIncrement*ones(obj.nSurf,1);zeros(obj.nSurf,1)];
            else
                error('OuterLQR:SharedWingDaisyChainVarPer', ...
                    'Unsupported cfg.ctrl.var_per = %d.',obj.varPer);
            end
            info = struct('enabled',true,'symmetricIncrement',wingIncrement, ...
                'virtualCommand',virtualCommand,'rigidIncrement',duRigid, ...
                'realizedMultirate',realized,'commandState',commandState, ...
                'wingSaturated',abs(wingIncrement-virtualCommand(1))> ...
                    10*eps(max(1,wingLimit)));
        end

        %--------------------------------------------------------------
        function [A,B] = buildLinearModel(obj,cfg)
        %BUILDLINEARMODEL Build rigid-body linear model for outer LQR.
        %
        % Priority:
        %   1. Use cfg.outerLQR.A/B only if both exist, are nonempty, and have
        %      correct dimensions.
        %   2. Otherwise finite-difference the internal fallback rigid-body model.
        %
        % Full rigid-body state:
        %   x_rb = [u; v; w; phi; theta; psi; p; q; r]
        %
        % Full outer input:
        %   u_rb = [delta_e; delta_a; delta_r; thrust]
        
            % --------------------------------------------------------------
            % 1. Try user-supplied A/B, but only if valid and nonempty.
            % --------------------------------------------------------------
            hasUserA = isfield(cfg,'outerLQR') && ...
                       isfield(cfg.outerLQR,'A') && ...
                       ~isempty(cfg.outerLQR.A);
        
            hasUserB = isfield(cfg,'outerLQR') && ...
                       isfield(cfg.outerLQR,'B') && ...
                       ~isempty(cfg.outerLQR.B);
        
            if hasUserA || hasUserB
                if ~(hasUserA && hasUserB)
                    error(['OuterLQR:PartialLinearModel\n', ...
                           'cfg.outerLQR.A and cfg.outerLQR.B must either both be ', ...
                           'provided and nonempty, or both omitted/empty.']);
                end
        
                A = cfg.outerLQR.A;
                B = cfg.outerLQR.B;
        
                if ~isequal(size(A),[obj.nx,obj.nx])
                    error('OuterLQR:BadUserA', ...
                          'cfg.outerLQR.A must be %d x %d. Got %d x %d.', ...
                          obj.nx,obj.nx,size(A,1),size(A,2));
                end
        
                if ~isequal(size(B),[obj.nx,obj.nuRB])
                    error('OuterLQR:BadUserB', ...
                          'cfg.outerLQR.B must be %d x %d. Got %d x %d.', ...
                          obj.nx,obj.nuRB,size(B,1),size(B,2));
                end
        
                fprintf('[OuterLQR] Using user-supplied cfg.outerLQR.A/B.\n');
                return
            end
        
            % --------------------------------------------------------------
            % 2. Finite-difference fallback rigid-body model.
            % --------------------------------------------------------------
            fprintf('[OuterLQR] No valid cfg.outerLQR.A/B supplied. Building A/B by finite difference.\n');
        
            x0 = obj.xTrim(:);
            u0 = obj.uTrimRB(:);
        
            if numel(x0) ~= obj.nx
                error('OuterLQR:BadTrimState', ...
                      'obj.xTrim must have length %d. Got %d.', obj.nx,numel(x0));
            end
        
            if numel(u0) ~= obj.nuRB
                error('OuterLQR:BadTrimInput', ...
                      'obj.uTrimRB must have length %d. Got %d.', obj.nuRB,numel(u0));
            end
        
            f0 = obj.rbDynamics(x0,u0);
        
            if isempty(f0) || numel(f0) ~= obj.nx
                error('OuterLQR:BadDynamicsOutput', ...
                      'rbDynamics(xTrim,uTrim) must return length %d. Got length %d.', ...
                      obj.nx,numel(f0));
            end
        
            if any(~isfinite(f0))
                error('OuterLQR:BadDynamicsOutput', ...
                      'rbDynamics(xTrim,uTrim) returned NaN/Inf.');
            end
        
            A = zeros(obj.nx,obj.nx);
            B = zeros(obj.nx,obj.nuRB);
        
            % State perturbation floors:
            %   velocity states: m/s
            %   angle states: rad
            %   rate states: rad/s
            hX_floor = [ ...
                1e-4; 1e-4; 1e-4; ...
                1e-6; 1e-6; 1e-6; ...
                1e-6; 1e-6; 1e-6];
        
            for i = 1:obj.nx
                h = 1e-6 * max(1,abs(x0(i)));
                h = max(h,hX_floor(i));
        
                xp = x0;
                xm = x0;
        
                xp(i) = xp(i) + h;
                xm(i) = xm(i) - h;
        
                fp = obj.rbDynamics(xp,u0);
                fm = obj.rbDynamics(xm,u0);
        
                A(:,i) = (fp(:) - fm(:))/(2*h);
            end
        
            % Input perturbation floors:
            %   surfaces: rad
            %   thrust: N
            hU_floor = [1e-5; 1e-5; 1e-5; 1e-4];
        
            for j = 1:obj.nuRB
                h = 1e-6 * max(1,abs(u0(j)));
                h = max(h,hU_floor(j));
        
                up = u0;
                um = u0;
        
                up(j) = up(j) + h;
                um(j) = um(j) - h;
        
                fp = obj.rbDynamics(x0,up);
                fm = obj.rbDynamics(x0,um);
        
                B(:,j) = (fp(:) - fm(:))/(2*h);
            end
        
            % --------------------------------------------------------------
            % 3. Final sanity checks.
            % --------------------------------------------------------------
            if any(~isfinite(A(:))) || any(~isfinite(B(:)))
                error('OuterLQR:FiniteDifferenceFailed', ...
                      'Finite-differenced A/B contain NaN or Inf.');
            end
        
            fprintf('[OuterLQR] Built finite-difference A/B.\n');
            fprintf('           ||A||_F = %.3e, ||B||_F = %.3e\n', norm(A,'fro'), norm(B,'fro'));
            fprintf('           B column norms: ');
            fprintf('%.3e ',vecnorm(B,2,1));
            fprintf('\n');
        end

        %--------------------------------------------------------------
        function xdot = rbDynamics(obj,x,uRB)
        % Rigid-body fallback dynamics for LQR linearization.
        %
        % x   = [u;v;w;phi;theta;psi;p;q;r]
        % uRB= [delta_e;delta_a;delta_r;thrust]
        %--------------------------------------------------------------
            vB     = x(1:3);
            euler  = x(4:6);
            omegaB = x(7:9);

            delta_e = uRB(1);
            delta_a = uRB(2);
            delta_r = uRB(3);
            thrust  = uRB(4);

            m = obj.rbParams.mass;
            I = obj.rbParams.I_B;

            U = norm(vB);

            if U < 1e-8
                U = max(1e-8,obj.safeGetFlightSpeed(obj.cfg));
                alpha = obj.safeGetAlpha(obj.cfg);
                beta = 0;
            else
                alpha = atan2(vB(3),vB(1));
                beta  = asin(max(-1,min(1,vB(2)/U)));
            end

            % ----------------------------------------------------------
            % Frozen wing clamp force/moment at trim.
            % ----------------------------------------------------------
            if isfield(obj.trim,'Clamp6')
                Clamp6 = obj.trim.Clamp6(:);
            elseif isfield(obj.trim,'clamp6')
                Clamp6 = obj.trim.clamp6(:);
            else
                Clamp6 = zeros(6,1);
            end

            if numel(Clamp6) ~= 6
                Clamp6 = zeros(6,1);
            end

            Fwing = Clamp6(1:3);
            Mwing = Clamp6(4:6);

            % ----------------------------------------------------------
            % Tail / fin models if available.
            % ----------------------------------------------------------
            if exist('tailAeroForceMoment','file') == 2
                [Ftail,Mtail] = tailAeroForceMoment(U,alpha,delta_e,obj.cfg);
            else
                Ftail = zeros(3,1);
                Mtail = zeros(3,1);
            end

            if exist('finAeroForceMoment','file') == 2
                [Ffin,Mfin] = finAeroForceMoment(U,beta,delta_r,obj.cfg);
            else
                Ffin = zeros(3,1);
                Mfin = zeros(3,1);
            end

            Ftail = Ftail(:); Mtail = Mtail(:);
            Ffin  = Ffin(:);  Mfin  = Mfin(:);

            % ----------------------------------------------------------
            % Simple aileron roll moment and angular damping fallback.
            % ----------------------------------------------------------
            L_da = obj.getOuterScalar('aileronRollMomentGain',0);
            Lp   = obj.getOuterScalar('Lp',0);
            Mq   = obj.getOuterScalar('Mq',0);
            Nr   = obj.getOuterScalar('Nr',0);

            M_ail = [L_da*delta_a + Lp*omegaB(1);
                     Mq*omegaB(2);
                     Nr*omegaB(3)];

            % ----------------------------------------------------------
            % Gravity and thrust.
            % ----------------------------------------------------------
            C_BI = obj.dcmBodyToInertial(euler);
            C_IB = C_BI.';

            Fgrav_I = [0;0;m*obj.rbParams.g];
            Fgrav_B = C_IB * Fgrav_I;

            Fthrust = [thrust;0;0];

            if isfield(obj.rbParams,'rThrust_B')
                rT = obj.rbParams.rThrust_B(:);
            else
                rT = zeros(3,1);
            end

            Mthrust = cross(rT,Fthrust);

            Ftot = Fwing + Ftail + Ffin + Fgrav_B + Fthrust;
            Mtot = Mwing + Mtail + Mfin + M_ail + Mthrust;

            % ----------------------------------------------------------
            % Newton-Euler equations.
            % ----------------------------------------------------------
            vBdot = Ftot/m - cross(omegaB,vB);
            eulDot = obj.eulerRates321(euler,omegaB);
            omegaDot = I \ (Mtot - cross(omegaB,I*omegaB));

            xdot = [vBdot; eulDot; omegaDot];
        end

        %--------------------------------------------------------------
        function [Q,R] = buildWeights(obj,cfg)
        % Bryson-rule or user-provided LQR weights.
        %--------------------------------------------------------------
            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'Q') && ~isempty(cfg.outerLQR.Q)
                Q = cfg.outerLQR.Q;
            else
                xMax = cfg.outerLQR.xMax(:);
                assert(numel(xMax)==obj.nx, 'cfg.outerLQR.xMax must have length 9.');
                Q = diag(1./max(xMax,eps).^2);
            end

            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'R') && ~isempty(cfg.outerLQR.R)
                R = cfg.outerLQR.R;
            else
                uMax = cfg.outerLQR.uMax(:);
                assert(numel(uMax)==obj.nuRB, 'cfg.outerLQR.uMax must have length 4.');
                R = diag(1./max(uMax,eps).^2);
            end

            Q = 0.5*(Q+Q.');
            R = 0.5*(R+R.');
        end

        %--------------------------------------------------------------
        function solveGain(obj)
        % Solve continuous or discrete LQR.
        %--------------------------------------------------------------
            Ctrb = ctrb(obj.A,obj.B);
            obj.controllabilityRank = rank(Ctrb,1e-8);
            
            if obj.useDiscrete
                sysc = ss(obj.A,obj.B,eye(obj.nx),zeros(obj.nx,obj.nuRB));
                sysd = c2d(sysc,obj.Ts,'zoh');

                obj.Ad = sysd.A;
                obj.Bd = sysd.B;
            %     fprintf('\n');
            % fprintf('======================================================================\n');
            % fprintf(' OUTER LQR DIAGNOSTICS BEFORE DLQR\n');
            % fprintf('======================================================================\n');
            % 
            % fprintf('  size(A) = %d x %d\n', size(obj.Ad,1), size(obj.Ad,2));
            % fprintf('  size(B) = %d x %d\n', size(obj.Bd,1), size(obj.Bd,2));
            % 
            % fprintf('  rank(ctrb(Ad,Bd)) = %d / %d\n', ...
            %         rank(ctrb(obj.Ad,obj.Bd),1e-8), size(obj.Ad,1));
            % 
            % fprintf('  input column norms of Bd:\n');
            % disp(vecnorm(obj.Bd,2,1));
            % 
            % fprintf('  eig(R):\n');
            % disp(eig(0.5*(obj.R+obj.R.')).');
            % 
            % fprintf('  eig(Ad):\n');
            % disp(eig(obj.Ad).');
            % 
            % tolPBH = 1e-8;
            % lam = eig(obj.Ad);
            % badModes = [];
            % 
            % for ii = 1:numel(lam)
            %     if abs(lam(ii)) >= 1 - 1e-8
            %         PBH = [lam(ii)*eye(size(obj.Ad)) - obj.Ad, obj.Bd];
            %         rPBH = rank(PBH,tolPBH);
            % 
            %         if rPBH < size(obj.Ad,1)
            %             badModes(end+1,1) = lam(ii); %#ok<AGROW>
            %             fprintf('  UNCONTROLLABLE UNIT-CIRCLE MODE: lambda = %.6g%+.6gi, PBH rank %d/%d\n', ...
            %                     real(lam(ii)), imag(lam(ii)), rPBH, size(obj.Ad,1));
            %         end
            %     end
            % end
            % 
            % fprintf('======================================================================\n\n');
                % [obj.K,~,~] = dlqr(obj.Ad,obj.Bd,obj.Q,obj.R);
                % [obj.K,~,~] = dlqr(sysd,obj.Q,obj.R);
                ix = obj.cfg.outerLQR.designStateIdx;
                iu = obj.cfg.outerLQR.designInputIdx;
                
                Ause = obj.Ad(ix,ix);
                Buse = obj.Bd(ix,iu);
                Quse = obj.Q(ix,ix);
                Ruse = obj.R(iu,iu);
                % obj.Ad = Ause;
                % obj.Bd = Buse;
                % obj.Q = Quse;
                % obj.R = Ruse;
                Kred = dlqr(Ause,Buse,Quse,Ruse);
                
                Kfull = zeros(4,9);
                Kfull(iu,ix) = Kred;
                obj.K = Kfull;
                obj.eigCL = eig(obj.Ad - obj.Bd*obj.K);
            else
                [obj.K,~,~] = lqr(obj.A,obj.B,obj.Q,obj.R);
                obj.eigCL = eig(obj.A - obj.B*obj.K); 
            end
        end

        %--------------------------------------------------------------
        function x = rbToState(~,rb)
        % Convert plant rigid-body sensor/state struct to LQR state.
        %--------------------------------------------------------------
            if isstruct(rb)
                if isfield(rb,'v_B')
                    vB = rb.v_B(:);
                elseif isfield(rb,'vB')
                    vB = rb.vB(:);
                else
                    error('OuterLQR:rbToState','rb must contain v_B.');
                end

                if isfield(rb,'euler')
                    euler = rb.euler(:);
                elseif isfield(rb,'angles')
                    euler = rb.angles(:);
                else
                    error('OuterLQR:rbToState','rb must contain euler.');
                end

                if isfield(rb,'omega_B')
                    omegaB = rb.omega_B(:);
                elseif isfield(rb,'omegaB')
                    omegaB = rb.omegaB(:);
                elseif isfield(rb,'rates')
                    omegaB = rb.rates(:);
                else
                    error('OuterLQR:rbToState','rb must contain omega_B.');
                end
            else
                rb = rb(:);
                if numel(rb) == 9
                    vB = rb(1:3);
                    euler = rb(4:6);
                    omegaB = rb(7:9);
                elseif numel(rb) == 12
                    vB = rb(4:6);
                    euler = rb(7:9);
                    omegaB = rb(10:12);
                else
                    error('OuterLQR:rbToState', ...
                          'rb vector must have length 9 or 12.');
                end
            end

            x = [vB(:); euler(:); omegaB(:)];
        end

        %--------------------------------------------------------------
        function xRef = refToState(obj,ref)
        % Build reference state. Defaults to trim.
        %--------------------------------------------------------------
            xRef = obj.xTrim;

            if isempty(ref)
                return
            end

            if isstruct(ref)
                if isfield(ref,'v_B_ref')
                    xRef(1:3) = ref.v_B_ref(:);
                elseif isfield(ref,'v_B')
                    xRef(1:3) = ref.v_B(:);
                end

                if isfield(ref,'euler_ref')
                    xRef(4:6) = ref.euler_ref(:);
                elseif isfield(ref,'euler')
                    xRef(4:6) = ref.euler(:);
                end

                if isfield(ref,'omega_B_ref')
                    xRef(7:9) = ref.omega_B_ref(:);
                elseif isfield(ref,'omega_B')
                    xRef(7:9) = ref.omega_B(:);
                end
            else
                ref = ref(:);
                if numel(ref) == 9
                    xRef = ref;
                else
                    error('OuterLQR:refToState','Reference vector must have length 9.');
                end
            end
        end

        %--------------------------------------------------------------
        function dx = wrapStateError(~,dx)
        % Wrap Euler angle errors.
        %--------------------------------------------------------------
            dx(4:6) = atan2(sin(dx(4:6)),cos(dx(4:6)));
        end

        %--------------------------------------------------------------
        function uRB = saturateOuterRBInput(obj,uRB,uCenter)
        % Limit rigid-body input perturbations about the trim command.
        %--------------------------------------------------------------
            if nargin < 3
                uCenter = obj.uTrimRB(:);
            end
            uCenter = uCenter(:);
            assert(numel(uCenter)==obj.nuRB && all(isfinite(uCenter)) && ...
                uCenter(4)>=0, ...
                'OuterLQR:InputLimitCenter', ...
                'The outer-input limit center must be finite with T >= 0.');
            uMax = obj.cfg.outerLQR.uMax(:);
            assert(numel(uMax) == obj.nuRB, ...
                'cfg.outerLQR.uMax must have length 4.');
            assert(all(isfinite(uMax)) && all(uMax > 0), ...
                'cfg.outerLQR.uMax must contain finite positive limits.');

            du = uRB(:) - uCenter;
            du = min(max(du, -uMax), uMax);
            uRB = uCenter + du;

            % The perturbation limit is trim-relative; the physical thrust
            % floor applies to the reconstructed total command.
            uRB(4) = max(0,uRB(4));
        end

        %--------------------------------------------------------------
        function [uFlexOuter, rbCmd] = allocateCommand(obj,uRB)
        % Map [delta_e; delta_a; delta_r; T] into wing-surface vector and
        % rigid-body command struct.
        %--------------------------------------------------------------
            delta_e = uRB(1);
            delta_a = uRB(2);
            delta_r = uRB(3);
            thrust  = uRB(4);

            defOuter = zeros(obj.nSurf,1);

            if obj.useElevatorOnWingSurfaces
                defOuter = defOuter + obj.elevatorSurfaceMap * delta_e;
            end

            if obj.useAileronOnWingSurfaces
                defOuter = defOuter + obj.aileronSurfaceMap * delta_a;
            end

            if obj.varPer == 1
                uFlexOuter = defOuter;
            elseif obj.varPer == 2
                % Rate channels are actuator states/commands, not directly
                % generated by LQR. Let actuator chain create rates.
                uFlexOuter = [defOuter; zeros(obj.nSurf,1)];
            else
                error('OuterLQR:allocateCommand', ...
                      'Unsupported cfg.ctrl.var_per = %d.', obj.varPer);
            end

            rbCmd = struct();
            rbCmd.delta_e = delta_e;
            rbCmd.delta_a = delta_a;
            rbCmd.delta_r = delta_r;
            rbCmd.thrust  = thrust;
        end

        %--------------------------------------------------------------
        function xTrim = buildTrimState(obj,cfg,trim)
        % Construct trim x_rb = [u;v;w;phi;theta;psi;p;q;r].
        %--------------------------------------------------------------
            U0 = obj.safeGetFlightSpeed(cfg);

            if isfield(trim,'alphaDeg')
                alpha0 = deg2rad(trim.alphaDeg);
            else
                alpha0 = obj.safeGetAlpha(cfg);
            end

            vB0 = [U0*cos(alpha0); 0; U0*sin(alpha0)];

            if isfield(trim,'rb') && isfield(trim.rb,'euler')
                euler0 = trim.rb.euler(:);
            elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'euler0')
                euler0 = cfg.rigidEOMset.euler0(:);
            else
                euler0 = [0; alpha0; 0];
            end
            assert(numel(euler0) == 3 && all(isfinite(euler0)), ...
                'OuterLQR:InvalidTrimEuler', ...
                'The rigid trim Euler reference must contain three finite values.');

            omega0 = zeros(3,1);

            xTrim = [vB0; euler0; omega0];
        end

        %--------------------------------------------------------------
        function uTrim = buildTrimInput(~,cfg,trim)
        % u_rb = [delta_e; delta_a; delta_r; thrust].
        %--------------------------------------------------------------
            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,'uTrim')
                uTrim = cfg.outerLQR.uTrim(:);
                assert(numel(uTrim)==4,'cfg.outerLQR.uTrim must have length 4.');
                return
            end

            delta_e = 0;
            delta_a = 0;
            delta_r = 0;

            if isfield(trim,'deltaDeg')
                d = deg2rad(trim.deltaDeg);
                if isscalar(d)
                    delta_e = d;
                elseif numel(d) >= 1
                    delta_e = d(1);
                end
            end

            if isfield(trim,'thrust')
                T = trim.thrust;
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'thrust')
                T = cfg.trim.thrust;
            else
                T = 0;
            end

            uTrim = [delta_e; delta_a; delta_r; T];
        end

        %--------------------------------------------------------------
        function map = buildSurfaceMap(obj,cfg,fieldName,defaultMap)
        % Build nSurf x 1 allocation map.
        %--------------------------------------------------------------
            if isfield(cfg,'outerLQR') && isfield(cfg.outerLQR,fieldName)
                map = cfg.outerLQR.(fieldName)(:);
            else
                map = defaultMap(:);
            end

            if numel(map) ~= obj.nSurf
                if isscalar(map)
                    map = repmat(map,obj.nSurf,1);
                else
                    error('OuterLQR:SurfaceMap', ...
                          '%s must be scalar or length n_surf.', fieldName);
                end
            end
        end

        %--------------------------------------------------------------
        function p = getRigidParams(~,cfg)
        % Rigid mass/inertia.
        %--------------------------------------------------------------
            p = struct();

            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'mass')
                p.mass = cfg.rigidEOMset.mass;
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'mass')
                p.mass = cfg.flight.mass;
            else
                rParams = RigidBody.methods.paramsRigid_PazyUAV();
                p.mass = rParams.mass;
            end

            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I_B')
                p.I_B = cfg.rigidEOMset.I_B;
            elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I')
                p.I_B = cfg.rigidEOMset.I;
            else
                rParams = RigidBody.methods.paramsRigid_PazyUAV();
                if isfield(rParams,'I_B')
                    p.I_B = rParams.I_B;
                elseif isfield(rParams,'I')
                    p.I_B = rParams.I;
                elseif all(isfield(rParams,{'Ixx','Iyy','Izz'}))
                    p.I_B = diag([rParams.Ixx,rParams.Iyy,rParams.Izz]);
                else
                    error('OuterLQR:RigidParams','Could not infer rigid-body inertia.');
                end
            end

            p.I_B = 0.5*(p.I_B+p.I_B.');
            if isfield(cfg,'rigidEOMset') && ...
                    isfield(cfg.rigidEOMset,'massOwnership')
                ownership = cfg.rigidEOMset.massOwnership;
                assert(isstruct(ownership) && isfield(ownership,'owner'), ...
                    'OuterLQR:RigidMassOwnership', ...
                    'The rigid mass-ownership contract is incomplete.');
                if string(ownership.owner)== ...
                        "nonwing_source_owned_wing_reaction"
                    expected = RigidBody.methods.paramsRigid_PazyUAV();
                    expected = expected.nonwing;
                    assert(abs(p.mass-expected.mass)<=1e-14 && ...
                        norm(p.I_B-expected.J,'fro')<=1e-14, ...
                        'OuterLQR:MixedRigidMassOwnership', ...
                        ['Outer LQR must use the same non-wing mass and ', ...
                         'inertia as the nonlinear runtime.']);
                end
                p.massOwnership = ownership;
            end

            if isfield(cfg,'flight') && isfield(cfg.flight,'g')
                p.g = cfg.flight.g;
            else
                p.g = 9.807;
            end

            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'rThrust_B')
                p.rThrust_B = cfg.rigidEOMset.rThrust_B(:);
            else
                p.rThrust_B = zeros(3,1);
            end
        end

        %--------------------------------------------------------------
        function val = getOuterScalar(obj,name,defaultVal)
            if isfield(obj.cfg,'outerLQR') && isfield(obj.cfg.outerLQR,name)
                val = obj.cfg.outerLQR.(name);
            else
                val = defaultVal;
            end
        end

        %--------------------------------------------------------------
        function U = safeGetFlightSpeed(~,cfg)
            if isfield(cfg,'flight') && isfield(cfg.flight,'U_inf')
                U = cfg.flight.U_inf;
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'U_inf')
                U = cfg.trim.U_inf;
            else
                U = 0;
            end
        end

        %--------------------------------------------------------------
        function alpha = safeGetAlpha(~,cfg)
            if isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
                alpha = deg2rad(cfg.flight.aoa_deg);
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'alphaDeg')
                alpha = deg2rad(cfg.trim.alphaDeg);
            else
                alpha = 0;
            end
        end

        %--------------------------------------------------------------
        function C_BI = dcmBodyToInertial(~,euler)
            phi = euler(1); theta = euler(2); psi = euler(3);

            cphi = cos(phi);   sphi = sin(phi);
            cth  = cos(theta); sth  = sin(theta);
            cpsi = cos(psi);   spsi = sin(psi);

            C_BI = [ ...
                cth*cpsi, sphi*sth*cpsi - cphi*spsi, cphi*sth*cpsi + sphi*spsi;
                cth*spsi, sphi*sth*spsi + cphi*cpsi, cphi*sth*spsi - sphi*cpsi;
                -sth,     sphi*cth,                  cphi*cth ];
        end

        %--------------------------------------------------------------
        function eulerDot = eulerRates321(~,euler,omegaB)
            phi = euler(1);
            theta = euler(2);

            p = omegaB(1);
            q = omegaB(2);
            r = omegaB(3);

            ctheta = cos(theta);

            if abs(ctheta) < 1e-8
                ctheta = sign(ctheta + eps)*1e-8;
            end

            T = [ ...
                1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
                0, cos(phi),            -sin(phi);
                0, sin(phi)/ctheta,      cos(phi)/ctheta ];

            eulerDot = T*[p;q;r];
        end

        %--------------------------------------------------------------
        function printDesignReport(obj)
            fprintf('\n');
            fprintf('======================================================================\n');
            fprintf(' OUTER LQR DESIGN REPORT\n');
            fprintf('======================================================================\n');
            fprintf('  Ts_outer             : %.6g s\n', obj.Ts);
            fprintf('  useDiscrete          : %d\n', obj.useDiscrete);
            fprintf('  nx, nuRB             : %d, %d\n', obj.nx, obj.nuRB);
            fprintf('  controllability rank : %d / %d\n', obj.controllabilityRank, obj.nx);

            if obj.useDiscrete
                fprintf('  max |eig(Ad-BdK)|    : %.6g\n', max(abs(obj.eigCL)));
            else
                fprintf('  max real eig(A-BK)   : %.6g\n', max(real(obj.eigCL)));
            end

            fprintf('  K size               : %d x %d\n', size(obj.K,1), size(obj.K,2));
            fprintf('======================================================================\n\n');
        end
    end
end
