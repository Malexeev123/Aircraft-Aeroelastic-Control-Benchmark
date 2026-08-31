classdef PlantRunTime < handle
    properties
        % ---------------- configuration / model handles --------------------
    cfg
    beam
    aero
    base
    idx
    trim

    model                   % AeroFlex.sim.ROMIntegrator

    % ---------------- flexible ROM state -------------------------------
    xFlex                   % current flexible/aeroelastic ROM state
    xFlex0                  % initial flexible/aeroelastic ROM state

    nx double               % flexible ROM state dimension
    nu double               % control input dimension
    nw double               % disturbance dimension

    % ---------------- timing -------------------------------------------
    dt double               % internal plant time step
    t double = 0            % current plant time
    k double = 0            % current plant step index

    % ---------------- body-case / coupled simulation --------------------
    bodyCase string = "wingOnly"
    useRigidBody logical = false

    % ---------------- rigid-body state and parameters -------------------
    rb                      % struct containing rigid-body state
    rbParams                % mass/inertia/geometry/etc.

    % ---------------- flight-condition cache ----------------------------
    flightCond              % struct: Uinf, alpha, beta, euler, omega_B

    % ---------------- last force/moment bookkeeping ---------------------
    last                    % struct storing Clamp6, Ftot_B, Mtot_B, etc.

    % ---------------- trim / fixed rigid inputs -------------------------
    trimThrust double = 0
    trimTailDelta double = 0
    trimWingControl
    lastWingControl
    trimMu

    % ---------------- debug --------------------------------------------
    % debugCoupled logical = false
    debugCoupled logical = true
         
    sensor 
    x 
   
    % ---------------- Scheduling ----------------------------------------
    ROMlib
    sched
    scheduleEnabled logical = false
    controlSampleScheduleHoldStepsRemaining double = 0
    controlSampleIndex double = 0
    reciprocalTangentCorrection
    
    end
    methods
        function obj = PlantRunTime(cfg, beam, aero, base, x0, trim)
        %PLANTRUNTIME Runtime plant wrapper for wing-only or rigid-flex simulation.
        %
        % Existing call site:
        %   plant = AeroFlex.sim.PlantRunTime(cfg, beam, aero, base, x0)
        %
        % Supports:
        %   cfg.sim.bodyCase = "wing_only"
        %   cfg.sim.bodyCase = "fully_coupled"
        
            obj.cfg  = cfg;
            obj.beam = beam;
            obj.aero = aero;
            obj.base = base;
            obj.trim = trim;
            obj.dt = cfg.sim.dt;
            obj.t  = 0;
            obj.k  = 0;
            % --------------------------------------------------------------
            % Body case
            % --------------------------------------------------------------
            if isfield(cfg,'sim') && isfield(cfg.sim,'bodyCase')
                obj.bodyCase = lower(string(cfg.sim.bodyCase));
            elseif isfield(cfg,'sim') && isfield(cfg.sim,'body_case')
                obj.bodyCase = lower(string(cfg.sim.body_case));
            else
                obj.bodyCase = "wingonly";
            end
            
            obj.bodyCase = erase(obj.bodyCase,"_");
            obj.useRigidBody = obj.bodyCase == "coupledfull";

            obj.idx = AeroFlex.core.buildIndexStruct(beam.Nm, aero.Na);

            obj.model  = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);

            if obj.bodyCase == "wingonly"
                obj.model.parConst.RateProject = struct( ...
                'projSet', false, ...
                'Pz', obj.beam.Pz);
            else
                obj.model.parConst.RateProject = struct( ...
                'projSet', true, ...
                'Pz', obj.beam.Pz);
            end

            if ismethod(obj.model,'rebuildConstrainedOperator')
                obj.model = obj.model.rebuildConstrainedOperator();
            end
            obj.idx = obj.model.idx;
            obj.sensor = AeroFlex.sensor.WingVelSensor(beam,cfg);
            obj.nx = size(obj.model.L,1);

            % Failure Handling:
            if isfield(cfg,'nu')
                obj.nu = cfg.nu;
            elseif isfield(cfg,'ctrl') && isfield(cfg.ctrl,'n_surf') && isfield(cfg.ctrl,'var_per')
                obj.nu = cfg.ctrl.n_surf * cfg.ctrl.var_per;
            else
                obj.nu = size(aero.forceMap.B_delta,2) + size(aero.forceMap.B_ddelta,2);
            end
        
            if isfield(cfg,'nw')
                obj.nw = cfg.nw;
            elseif isfield(aero,'forceMap') && isfield(aero.forceMap,'Bw')
                obj.nw = size(aero.forceMap.Bw,2);
            else
                obj.nw = 1;
            end
            % --------------------------------------------------------------
            % Flexible ROM initial state
            % --------------------------------------------------------------
            validateattributes(x0, {'numeric'}, {'real', 'vector', 'finite'}, ...
                mfilename, 'x0');
            x0 = x0(:);
        
            if numel(x0) < obj.nx
                error('PlantRunTime:InitialCondition', ...
                      'x0 has length %d, but flexible ROM needs at least nx = %d.', ...
                      numel(x0), obj.nx);
            end
        
            % x0 is an explicit runtime input.  Its optional appended rigid
            % tail is consumed by initRigidState below; the propagated
            % flexible state uses its leading nx entries.
            obj.x = x0;
            obj.xFlex0 = x0(1:obj.nx);
            obj.xFlex  = obj.xFlex0;
        
            % --------------------------------------------------------------
            % Rigid-body initialization
            % --------------------------------------------------------------
            obj.rbParams = obj.getRigidParams(cfg);
            obj.rb       = obj.initRigidState(cfg, x0);
        
            % --------------------------------------------------------------
            % Fixed trim commands
            % --------------------------------------------------------------
            
            obj.trimThrust = trim.thrust;    

            % if isfield(cfg,'trim') && isfield(cfg.trim,'deltaDeg')
            %     obj.trimTailDelta = deg2rad(cfg.trim.deltaDeg);
            %     if numel(obj.trimTailDelta) > 1
            %         obj.trimTailDelta = obj.trimTailDelta(1);
            %     end
            % else
            %     obj.trimTailDelta = 0;
            % end

            % -------------------------------------------------------------------------
            % Rear elevator/tail trim.
            % Prefer the solved trim struct. Fall back to cfg fields.
            % Stored internally in radians.
            % -------------------------------------------------------------------------
            if nargin >= 6 && isstruct(trim) && isfield(trim,'deltaElev')
                obj.trimTailDelta = trim.deltaElev;
            else
                obj.trimTailDelta = 0;
            end
            
            if numel(obj.trimTailDelta) > 1
                obj.trimTailDelta = obj.trimTailDelta(1);
            end
            
            obj.trimWingControl = obj.getTrimWingControl(trim);
            obj.lastWingControl = obj.trimWingControl;

            % --------------------------------------------------------------
            % Flight-condition cache
            % --------------------------------------------------------------
            obj.flightCond = struct();
            obj.flightCond.Uinf    = obj.safeGetFlightSpeed(cfg);
            obj.flightCond.alpha   = obj.safeGetAlpha(cfg);
            obj.flightCond.beta    = 0;
            obj.flightCond.euler   = obj.rb.euler;
            obj.flightCond.omega_B = obj.rb.omega_B;
            
            obj.trimMu = [obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha)];
            % --------------------------------------------------------------
            % Last-load cache
            % --------------------------------------------------------------
            obj.last = struct();
            obj.last.Clamp6 = zeros(6,1);
            obj.last.Fwing_B = zeros(3,1);
            obj.last.Mwing_B = zeros(3,1);
            obj.last.Ftail_B = zeros(3,1);
            obj.last.Mtail_B = zeros(3,1);
            obj.last.Ffin_B  = zeros(3,1);
            obj.last.Mfin_B  = zeros(3,1);
            obj.last.Fgrav_B = zeros(3,1);
            obj.last.Fthrust_B = zeros(3,1);
            obj.last.Mthrust_B = zeros(3,1);
            obj.last.Freciprocal_B = zeros(3,1);
            obj.last.Mreciprocal_B = zeros(3,1);
            obj.last.Ftot_B = zeros(3,1);
            obj.last.Mtot_B = zeros(3,1);
            obj.last.qRatio     = 1;
            obj.last.uWingTotal = obj.trimWingControl;
            obj.last.gust       = zeros(obj.nw,1);
            obj.last.sched_mu   = obj.trimMu;
            obj.last.rigidAttitudeChiOwner = struct( ...
                'enabled',false,'thetaAbsoluteRad',obj.rb.euler(2), ...
                'thetaEquilibriumRad',deg2rad(obj.trimMu(2)), ...
                'chiRad',zeros(3,1));
            obj.reciprocalTangentCorrection = struct( ...
                'enabled',false,'hiddenState',zeros(0,1), ...
                'lastVisibleCorrection',zeros(0,1), ...
                'lastInput',zeros(0,1));

            if isfield(cfg,'debug') && isfield(cfg.debug,'coupledPlant')
                obj.debugCoupled = logical(cfg.debug.coupledPlant);
            else
                obj.debugCoupled = true;
                % obj.debugCoupled = false;
            end
        end
        %--------------------------------------------------------------
        function installReciprocalTangentCorrection(obj,candidate)
        %INSTALLRECIPROCALTANGENTCORRECTION Install exact Case-A audit model.
        %   The correction replaces the accepted runtime's local tangent by
        %   adding source-minus-runtime increments after the unchanged
        %   nonlinear step. It is audit-only, default inactive, and valid
        %   only for the exact coupledFull Case-A source.

            assert(obj.bodyCase=="coupledfull", ...
                'PlantRunTime:ReciprocalTangentBodyCase', ...
                'The reciprocal tangent is valid only for coupledFull.');
            assert(isstruct(candidate) && isfield(candidate,'auditOnly') && ...
                logical(candidate.auditOnly) && ...
                isfield(candidate,'defaultActive') && ...
                ~logical(candidate.defaultActive) && ...
                isfield(candidate,'runtimeEligible') && ...
                logical(candidate.runtimeEligible), ...
                'PlantRunTime:ReciprocalTangentEligibility', ...
                'The reciprocal tangent must be qualified audit-only evidence.');
            assert(isfield(candidate,'schemaVersion') && ...
                ismember(string(candidate.schemaVersion),[ ...
                "phase18c-v17a-exact-casea-reciprocal-tangent-v2", ...
                "phase18c-v17a-exact-casea-reciprocal-tangent-v3", ...
                "phase18c-v17a-exact-casea-reciprocal-tangent-v4", ...
                "phase18c-v17a-scheduled-reciprocal-tangent-v1", ...
                "phase18c-v17a-scheduled-reciprocal-tangent-symmetric-audit-v1"]), ...
                'PlantRunTime:ReciprocalTangentSchema', ...
                'The reciprocal tangent schema is not supported.');
            assert(abs(double(candidate.sampleTimeSeconds)-obj.dt)<= ...
                10*eps(max(obj.dt,double(candidate.sampleTimeSeconds))), ...
                'PlantRunTime:ReciprocalTangentSampleTime', ...
                'The reciprocal tangent and runtime sample times differ.');
            visibleCount = double(candidate.visible.count);
            hiddenCount = double(candidate.hidden.count);
            inputCount = double(candidate.input.count);
            completeRuntimeSchema = ismember(string(candidate.schemaVersion), ...
                ["phase18c-v17a-exact-casea-reciprocal-tangent-v3", ...
                 "phase18c-v17a-exact-casea-reciprocal-tangent-v4", ...
                 "phase18c-v17a-scheduled-reciprocal-tangent-v1", ...
                 "phase18c-v17a-scheduled-reciprocal-tangent-symmetric-audit-v1"]);
            if completeRuntimeSchema
                assert(isfield(candidate,'runtimeState') && ...
                    isfield(candidate.runtimeState,'count') && ...
                    isfield(candidate.correction,'deltaAFullState') && ...
                    isfield(candidate.partition,'hiddenFromFullState') && ...
                    isfield(candidate.base,'correctionStateEquilibrium'), ...
                    'PlantRunTime:ReciprocalTangentFullStateContract', ...
                    'The complete runtime-state subtraction is incomplete.');
                correctionStateCount = double(candidate.runtimeState.count);
                visibleCorrectionMap = ...
                    candidate.correction.deltaAFullState;
                hiddenStateMap = candidate.partition.hiddenFromFullState;
                correctionStateEquilibrium = ...
                    candidate.base.correctionStateEquilibrium(:);
            else
                correctionStateCount = visibleCount;
                visibleCorrectionMap = candidate.correction.deltaAvv;
                hiddenStateMap = candidate.partition.Ahv;
                correctionStateEquilibrium = ...
                    candidate.base.visibleEquilibrium(:);
            end
            assert(visibleCount==69 && hiddenCount>0 && inputCount==7 && ...
                correctionStateCount>=visibleCount && ...
                isequal(size(visibleCorrectionMap), ...
                    [visibleCount,correctionStateCount]) && ...
                isequal(size(candidate.partition.Avh), ...
                    [visibleCount,hiddenCount]) && ...
                isequal(size(hiddenStateMap), ...
                    [hiddenCount,correctionStateCount]) && ...
                isequal(size(candidate.partition.Ahh), ...
                    [hiddenCount,hiddenCount]) && ...
                isequal(size(candidate.correction.deltaBv), ...
                    [visibleCount,inputCount]) && ...
                isequal(size(candidate.partition.Bh), ...
                    [hiddenCount,inputCount]), ...
                'PlantRunTime:ReciprocalTangentDimensions', ...
                'The reciprocal tangent dimensions or state order changed.');
            values = [visibleCorrectionMap(:); ...
                candidate.partition.Avh(:);hiddenStateMap(:); ...
                candidate.partition.Ahh(:); ...
                candidate.correction.deltaBv(:);candidate.partition.Bh(:); ...
                correctionStateEquilibrium; ...
                candidate.base.inputEquilibrium(:)];
            assert(all(isfinite(values)) && ...
                numel(correctionStateEquilibrium)==correctionStateCount && ...
                numel(candidate.base.inputEquilibrium)==inputCount && ...
                candidate.checks.elevatorChannelPresent && ...
                candidate.checks.thrustChannelPresent, ...
                'PlantRunTime:ReciprocalTangentContract', ...
                'The reciprocal tangent contains invalid values or inputs.');

            replacesFixedRootIncrement = string(candidate.schemaVersion)== ...
                "phase18c-v17a-exact-casea-reciprocal-tangent-v4" || ...
                startsWith(string(candidate.schemaVersion), ...
                "phase18c-v17a-scheduled-reciprocal-tangent");
            fixedRootWingEquilibrium = zeros(6,1);
            sourceRootOutput = struct();
            if replacesFixedRootIncrement
                assert(isfield(candidate,'rigidReactionOwner') && ...
                    string(candidate.rigidReactionOwner)== ...
                    "source_reciprocal_replaces_fixed_root_increment" && ...
                    isfield(candidate.base,'fixedRootWingEquilibrium') && ...
                    numel(candidate.base.fixedRootWingEquilibrium)==6 && ...
                    all(isfinite(candidate.base.fixedRootWingEquilibrium(:))), ...
                    'PlantRunTime:ReciprocalTangentReactionOwner', ...
                    ['The V4 reciprocal tangent must explicitly own the ', ...
                     'incremental rigid-wing reaction and preserve a finite ', ...
                     'six-component equilibrium wrench.']);
                fixedRootWingEquilibrium = ...
                    candidate.base.fixedRootWingEquilibrium(:);
                assert(isfield(candidate,'rootWrenchOutput') && ...
                    isfield(candidate.rootWrenchOutput,'runtimeStateMap') && ...
                    isfield(candidate.rootWrenchOutput,'runtimeInputMap') && ...
                    isfield(candidate.rootWrenchOutput,'equilibrium') && ...
                    isequal(size(candidate.rootWrenchOutput.runtimeStateMap), ...
                        [6,correctionStateCount+hiddenCount]) && ...
                    isequal(size(candidate.rootWrenchOutput.runtimeInputMap), ...
                        [6,inputCount]) && ...
                    numel(candidate.rootWrenchOutput.equilibrium)==6 && ...
                    all(isfinite([ ...
                        candidate.rootWrenchOutput.runtimeStateMap(:); ...
                        candidate.rootWrenchOutput.runtimeInputMap(:); ...
                        candidate.rootWrenchOutput.equilibrium(:)])), ...
                    'PlantRunTime:ReciprocalTangentRootOutput', ...
                    ['The V4 reciprocal tangent requires a finite ', ...
                     'source-owned six-component root-wrench output.']);
                sourceRootOutput = candidate.rootWrenchOutput;
            end

            wrenchMapMethod = "central_difference";
            if startsWith(string(candidate.schemaVersion), ...
                    "phase18c-v17a-scheduled-reciprocal-tangent")
                wrenchMapMethod = "analytic_rk4_audit";
            end
            [rigidWrenchMap,wrenchMapDetails] = ...
                obj.buildReciprocalRigidWrenchMap(wrenchMapMethod);
            correctionOperator = [visibleCorrectionMap, ...
                candidate.partition.Avh,candidate.correction.deltaBv];
            desiredRigidOperator = correctionOperator( ...
                candidate.visible.rigid,:);
            dynamicWrenchMap = rigidWrenchMap(1:6,:);
            wrenchCommandOperator = dynamicWrenchMap\ ...
                desiredRigidOperator(1:6,:);
            rigidProjectionResidual = norm(dynamicWrenchMap* ...
                wrenchCommandOperator-desiredRigidOperator(1:6,:),'fro')/ ...
                max(1,norm(desiredRigidOperator(1:6,:),'fro'));
            derivativeCertified = ...
                wrenchMapDetails.halfStepConvergence<=1e-7;
            if wrenchMapMethod=="analytic_rk4_audit"
                derivativeCertified = ...
                    wrenchMapDetails.complexStepRelativeError<=1e-10 && ...
                    wrenchMapDetails.stepParityRelativeError<=1e-12;
            end
            assert(wrenchMapDetails.dynamicRank==6 && ...
                derivativeCertified && ...
                wrenchMapDetails.conditionNumber<=1e3 && ...
                rigidProjectionResidual<=1e-4, ...
                'PlantRunTime:ReciprocalTangentRigidWrenchMap', ...
                ['The reciprocal rigid correction is not realizable by ', ...
                 'the existing six-component body-wrench path: method=%s, ', ...
                 'rank=%d, condition=%.6e, derivative=%.6e, ', ...
                 'step-parity=%.6e, projection=%.6e.'], ...
                wrenchMapMethod, ...
                wrenchMapDetails.rank,wrenchMapDetails.conditionNumber, ...
                wrenchMapDetails.derivativeCertificationResidual, ...
                wrenchMapDetails.stepParityRelativeError, ...
                rigidProjectionResidual);

            refinementUpdateLimit = 2;
            if isfield(obj.cfg,'ctrl') && isfield(obj.cfg.ctrl, ...
                    'caseBRigidWrenchRefinementAudit')
                request = obj.cfg.ctrl.caseBRigidWrenchRefinementAudit;
                required = {'enabled','auditOnly','changeId','caseId', ...
                    'maximumCorrectionUpdates'};
                assert(isstruct(request) && isscalar(request) && ...
                    all(isfield(request,required)), ...
                    'PlantRunTime:RigidWrenchRefinementRequest', ...
                    'The Case-B rigid-wrench refinement request is incomplete.');
                if logical(request.enabled)
                    assert(logical(request.auditOnly) && ...
                        string(request.changeId)== ...
                        "phase18c-v17a-caseb-rigid-wrench-refinement-extension-audit-v1" && ...
                        string(request.caseId)=="formal_case_b" && ...
                        double(request.maximumCorrectionUpdates)==3 && ...
                        startsWith(string(candidate.schemaVersion), ...
                            "phase18c-v17a-scheduled-reciprocal-tangent"), ...
                        'PlantRunTime:RigidWrenchRefinementApproval', ...
                        ['The additional rigid-wrench correction is ', ...
                         'authorized only for the approved scheduled ', ...
                         'Case-B audit.']);
                    refinementUpdateLimit = 3;
                end
            end

            obj.reciprocalTangentCorrection = struct( ...
                'enabled',true,'candidate',candidate, ...
                'correctionStateCount',correctionStateCount, ...
                'correctionStateEquilibrium',correctionStateEquilibrium, ...
                'visibleCorrectionMap',visibleCorrectionMap, ...
                'hiddenStateMap',hiddenStateMap, ...
                'hiddenState',zeros(candidate.hidden.count,1), ...
                'lastVisibleCorrection',zeros(candidate.visible.count,1), ...
                'lastInput',zeros(candidate.input.count,1), ...
                'rigidWrenchMap',rigidWrenchMap, ...
                'rigidWrenchCommandOperator',wrenchCommandOperator, ...
                'rigidWrenchMapDetails',wrenchMapDetails, ...
                'rigidProjectionResidual',rigidProjectionResidual, ...
                'lastWrenchCorrection',zeros(6,1), ...
                'lastRigidChartCorrection',zeros(3,1), ...
                'lastRigidRealizationResidual',0, ...
                'lastWrenchRefinementCount',0, ...
                'rigidWrenchRefinementUpdateLimit', ...
                    refinementUpdateLimit, ...
                'lastStepAccepted',true, ...
                'lastDomainRejected',false, ...
                'lastDomainReason',"", ...
                'lastDomainDetails',struct([]), ...
                'replacesFixedRootIncrement',replacesFixedRootIncrement, ...
                'fixedRootWingEquilibrium',fixedRootWingEquilibrium, ...
                'lastFixedRootIncrement',zeros(6,1), ...
                'lastBodyAppliedWingWrench',fixedRootWingEquilibrium, ...
                'sourceRootOutput',sourceRootOutput, ...
                'lastSourceRootWrench',zeros(6,1));
        end
        %--------------------------------------------------------------
        function applyScheduledReciprocalPacket(obj,packet)
        %APPLYSCHEDULEDRECIPROCALPACKET Refresh common maps, retain memory.
            required = {'schemaVersion','plantCandidate','queryTrim', ...
                'queryNativeEquilibrium','sourceIds','schedule'};
            assert(isstruct(packet) && isscalar(packet) && ...
                all(isfield(packet,required)) && ...
                string(packet.schemaVersion)== ...
                    "phase18c-v17a-scheduled-reciprocal-packet-v1" && ...
                numel(packet.queryNativeEquilibrium)==obj.nx && ...
                all(isfinite(packet.queryNativeEquilibrium(:))), ...
                'PlantRunTime:ScheduledReciprocalPacket', ...
                'The scheduled reciprocal packet is incomplete.');
            candidate = packet.plantCandidate;
            hiddenState = zeros(candidate.hidden.count,1);
            if obj.reciprocalTangentEnabled() && ...
                    isfield(obj.reciprocalTangentCorrection,'candidate') && ...
                    startsWith(string(obj.reciprocalTangentCorrection. ...
                        candidate.schemaVersion), ...
                        "phase18c-v17a-scheduled-reciprocal-tangent") && ...
                    numel(obj.reciprocalTangentCorrection.hiddenState)== ...
                    candidate.hidden.count
                hiddenState = obj.reciprocalTangentCorrection.hiddenState(:);
                if isfield(packet,'memoryCoordinateAudit') && ...
                        isfield(packet.memoryCoordinateAudit,'enabled') && ...
                        logical(packet.memoryCoordinateAudit.enabled)
                    hiddenState = obj.transportScheduledStableMemory( ...
                        packet,hiddenState);
                end
            end

            rigidSave = obj.rb;
            cleanup = onCleanup(@() obj.restoreRigidSnapshot(rigidSave));
            queryRigid = packet.queryTrim.rigid(:);
            assert(numel(queryRigid)==9 && ...
                all(isfinite(queryRigid)) && ...
                packet.queryTrim.thrust>=0, ...
                'PlantRunTime:ScheduledReciprocalTrim', ...
                'The scheduled reciprocal trim is invalid.');
            obj.rb.v_B = queryRigid(1:3);
            obj.rb.omega_B = queryRigid(4:6);
            obj.rb.euler = queryRigid(7:9);
            rigidCommand = struct('delta_e',packet.queryTrim.elevator, ...
                'delta_a',0,'delta_r',0,'thrust',packet.queryTrim.thrust);
            loads = obj.computeCoupledLoads( ...
                packet.queryNativeEquilibrium(:),packet.queryTrim.wing(:), ...
                zeros(obj.nw,1),rigidCommand);
            candidate.base.fixedRootWingEquilibrium = ...
                [loads.Fwing_B(:);loads.Mwing_B(:)];
            clear cleanup
            obj.restoreRigidSnapshot(rigidSave);

            obj.installReciprocalTangentCorrection(candidate);
            obj.reciprocalTangentCorrection.hiddenState = hiddenState;
            obj.reciprocalTangentCorrection.scheduledQueryTrim = ...
                packet.queryTrim;
            obj.reciprocalTangentCorrection.scheduledSourceIds = ...
                string(packet.sourceIds(:));
            obj.reciprocalTangentCorrection.scheduledWeights = ...
                double(packet.schedule.weights(:));
            if isfield(packet,'members')
                obj.reciprocalTangentCorrection.scheduledMembers = ...
                    packet.members;
            end
        end
        %--------------------------------------------------------------
        function hiddenState = transportScheduledStableMemory( ...
                obj,packet,hiddenState)
        %TRANSPORTSCHEDULEDSTABLEMEMORY Preserve physical latent state.
        %   The scheduled packet stores independent source-local stable
        %   charts z = h - F*v.  Across a schedule refresh, reconstruct h
        %   in the old source chart and express it in the new source chart.

            runtime = obj.reciprocalTangentCorrection;
            assert(isfield(runtime,'scheduledMembers') && ...
                isfield(packet,'members') && ...
                numel(runtime.scheduledMembers)==numel(packet.members), ...
                'PlantRunTime:ScheduledStableMemoryTransport', ...
                'Stable reciprocal memory cannot be transported safely.');
            oldMembers = runtime.scheduledMembers;
            newMembers = packet.members;
            assert(isequal(string({oldMembers.sourceId})', ...
                string({newMembers.sourceId})'), ...
                'PlantRunTime:ScheduledStableMemorySourceOrder', ...
                'The stable reciprocal source order changed at refresh.');
            state = obj.currentReciprocalRuntimeState();
            hiddenCounts = arrayfun(@(member) ...
                double(member.candidate.hidden.count),newMembers);
            offsets = [0;cumsum(hiddenCounts(:))];
            assert(numel(hiddenState)==offsets(end), ...
                'PlantRunTime:ScheduledStableMemoryCount', ...
                'The stable reciprocal hidden-memory count changed.');
            for memberIndex = 1:numel(newMembers)
                localHidden = offsets(memberIndex)+(1:hiddenCounts(memberIndex));
                oldCandidate = oldMembers(memberIndex).candidate;
                newCandidate = newMembers(memberIndex).candidate;
                assert(isfield(oldCandidate,'memoryCoordinateAudit') && ...
                    isfield(newCandidate,'memoryCoordinateAudit') && ...
                    isfield(oldCandidate.memoryCoordinateAudit, ...
                    'outputInjection') && ...
                    isfield(newCandidate.memoryCoordinateAudit, ...
                    'outputInjection') && ...
                    isfield(oldMembers(memberIndex),'runtimeFromQueryState') && ...
                    isfield(newMembers(memberIndex),'runtimeFromQueryState'), ...
                    'PlantRunTime:ScheduledStableMemoryContract', ...
                    'The stable reciprocal memory contract is incomplete.');
                oldDeparture = oldMembers(memberIndex). ...
                    runtimeFromQueryState*state- ...
                    oldCandidate.base.correctionStateEquilibrium(:);
                newDeparture = newMembers(memberIndex). ...
                    runtimeFromQueryState*state- ...
                    newCandidate.base.correctionStateEquilibrium(:);
                oldVisible = oldDeparture( ...
                    oldCandidate.runtimeState.sourceVisible);
                newVisible = newDeparture( ...
                    newCandidate.runtimeState.sourceVisible);
                oldF = oldCandidate.memoryCoordinateAudit.outputInjection;
                newF = newCandidate.memoryCoordinateAudit.outputInjection;
                rawHidden = hiddenState(localHidden)+oldF*oldVisible;
                hiddenState(localHidden) = rawHidden-newF*newVisible;
            end
            assert(all(isfinite(hiddenState)), ...
                'PlantRunTime:ScheduledStableMemoryFinite', ...
                'The stable reciprocal memory transport was nonfinite.');
        end
        %--------------------------------------------------------------
        function state = currentReciprocalRuntimeState(obj)
            state = [obj.xFlex(obj.idx.q1);obj.xFlex(obj.idx.q2); ...
                obj.xFlex(obj.idx.qGam);obj.rb.v_B(:); ...
                obj.rb.omega_B(:);obj.rb.euler(:); ...
                obj.xFlex(obj.idx.qxi);obj.xFlex(obj.idx.chi)];
        end
        %--------------------------------------------------------------
        function disableReciprocalTangentCorrection(obj)
        %DISABLERECIPROCALTANGENTCORRECTION Restore accepted runtime behavior.
            obj.reciprocalTangentCorrection = struct( ...
                'enabled',false,'hiddenState',zeros(0,1), ...
                'lastVisibleCorrection',zeros(0,1), ...
                'lastInput',zeros(0,1));
        end
        %--------------------------------------------------------------
        function configureCaseBRigidWrenchRefinementAudit(obj,request)
        %CONFIGURECASEBRIGIDWRENCHREFINEMENTAUDIT Bind the approved S269 audit.
            required = {'enabled','auditOnly','changeId','caseId', ...
                'maximumCorrectionUpdates'};
            assert(isstruct(request) && isscalar(request) && ...
                all(isfield(request,required)) && ...
                logical(request.enabled) && logical(request.auditOnly) && ...
                string(request.changeId)== ...
                "phase18c-v17a-caseb-rigid-wrench-refinement-extension-audit-v1" && ...
                string(request.caseId)=="formal_case_b" && ...
                double(request.maximumCorrectionUpdates)==3 && ...
                obj.reciprocalTangentEnabled() && ...
                isfield(obj.reciprocalTangentCorrection, ...
                    'scheduledMembers'), ...
                'PlantRunTime:RigidWrenchRefinementApproval', ...
                ['The additional rigid-wrench correction is authorized ', ...
                 'only for the approved scheduled Case-B audit.']);
            obj.cfg.ctrl.caseBRigidWrenchRefinementAudit = request;
            obj.reciprocalTangentCorrection. ...
                rigidWrenchRefinementUpdateLimit = 3;
        end
        %--------------------------------------------------------------
        function [z_k,t_k] = readSensors(obj, cfg, log)
            t_k = obj.t;
            % y = log( round(t_k/cfg.sim.dt)+1, :).' ;
            if cfg.forceRealSense
                z_k = obj.sensor.measure(obj.x(obj.model.idx.q1));
            else
                y = log( round(t_k/cfg.sim.dt)+1, :).' ;

                z_k = y;
            end
            % t_k = obj.t;
        end
        %--------------------------------------------------------------
        function x = step(obj,u_k,g_k,mode)
            % [obj.x,~] = obj.model.step(obj.x,u_k,g_k,[],[]);
            switch mode
                case "ROM"          % internal model
                    Nint = round(obj.cfg.ctrl.Ts / obj.cfg.sim.dt);
                    % % if ~exists(S)
                    % 
                    % S = [eye(144), zeros(144, 5)];
                    % % end
                    % x1 = obj.x;
                    for m = 1:Nint
                        [obj.x,~] = obj.model.step(obj.x, u_k, g_k, [], []);
                        % S = [eye(144), zeros(144, 5)];
                        % [obj.x,S] = obj.model.step(obj.x, u_k, g_k,  S, true);
                    end
                    % x2 = obj.x;
                    % % xS = x1+S*[x1; g_k; u_k];
                    % % xS = S*[x1; g_k; u_k];
                    % xS = S(1:length(x1),1:length(x1))*x1;
                    % SensErr = xS-x2;
                    % disp(SensErr);
                    % disp(norm(SensErr));
                    % [obj.x,~] = obj.model.step(obj.x,u_k,g_k,[],[]);
                case "SHARPy"       % read from UDP
                    obj.x = receiveStateFromUDP();              %# your code
                case "openLoop"     % cached log
                    obj.x = getCachedState(obj.t);              %# your code
            end
            % obj.t = obj.t + obj.cfg.sim.dt;
            obj.t = obj.t + obj.cfg.ctrl.Ts;      % ONE sample advance
            x = obj.x;
        end
        function xNext = stepCoupled(obj, u_cmd, gust, rbCmd)
        %STEPCOUPLED Advance one coupled rigid-flexible plant step.
        %
        % Commands are increments about trim.  The total wing, elevator, and thrust
        % commands are assembled internally.
        %
        % Runtime scheduling is handled as an admissible model switch:
        %   1. compute current flight condition,
        %   2. compose the total commands used by the plant,
        %   3. test/apply the candidate scheduled ROM at the current x,u,w,
        %   4. evaluate loads and propagate.
        
            if nargin < 4 || isempty(rbCmd)
                rbCmd = struct();
            end
        
            u_cmd = obj.sanitizeControl(u_cmd);
            gust  = obj.sanitizeDisturbance(gust);
        
            % Current rigid-body-derived scheduling variables.
            obj.flightCond = obj.computeFlightConditionFromRB();
        
            obj.updateWingFlightCondition();

            % Total controls used by both the load calculation and propagation.
            % These must be known before a scheduled model is tested/applied.
            uWing = obj.composeWingControl(u_cmd);
            rbTot = obj.composeRigidCommand(rbCmd);
        
            obj.lastWingControl = uWing;
            obj.last.uWingTotal = uWing;
            obj.last.gust       = gust;

            reciprocalVisible = [];
            reciprocalInput = [];
            if obj.reciprocalTangentEnabled()
                [reciprocalVisible,reciprocalInput] = ...
                    obj.reciprocalTangentInputs(u_cmd,gust,rbCmd);
            end
        
            % Apply/hold scheduled ROM using the same x,u,w that will be propagated.
            obj.updateScheduledROM(gust, uWing, rbTot);
            obj.applyPackageRelativeRigidAttitudeChi();

            % Loads are evaluated at the beginning of the step.  This keeps the first
            % runtime load check directly comparable with the trim residual.
            loads = obj.computeCoupledLoads(obj.xFlex, uWing, gust, rbTot);
            reciprocalStep = struct();
            if obj.reciprocalTangentEnabled()
                [loads,reciprocalStep] = ...
                    obj.prepareReciprocalTangentCorrection( ...
                        loads,reciprocalVisible,reciprocalInput);
            end
        
            % Optional first-step fixed-point check.  Keep it out of the main loop
            % because it costs an extra ROM step.
            if obj.debugCoupled && obj.k == 0
                if obj.scheduleEnabled
                    [xProbe, ~] = obj.model.step(obj.xFlex, uWing, gust, [], false);
                    qRatioPrint = 1;
                else
                    qRatioPrint = obj.last.qRatio;
                    [xProbe, ~] = obj.model.step( ...
                        obj.xFlex, uWing, gust, [], false, qRatioPrint);
                end
        
                dxProbe = xProbe - obj.xFlex;
        
                fprintf(['[flex fixed-point check] |dx|=%.3e |dx/dt|=%.3e ', ...
                         '|dq1|=%.3e |dq2|=%.3e |dqGam|=%.3e\n'], ...
                        norm(dxProbe), norm(dxProbe)/obj.model.dt, ...
                        norm(dxProbe(obj.idx.q1)), ...
                        norm(dxProbe(obj.idx.q2)), ...
                        norm(dxProbe(obj.idx.qGam)));
        
                fprintf(['[wing input] delta=[%+.6e %+.6e] ', ...
                         'deltaDot=[%+.6e %+.6e] qRatio=%.9f\n'], ...
                        uWing(1), uWing(2), uWing(3), uWing(4), qRatioPrint);
            end
        
            % Propagate flexible ROM with the active scheduled/frozen model.
            obj.xFlex = obj.stepWingROM(uWing, gust);
        
            % Propagate rigid body with beginning-of-step loads.
            obj.rb = obj.stepRigidBody(obj.rb, loads.Ftot_B, loads.Mtot_B);

            if obj.reciprocalTangentEnabled()
                obj.applyReciprocalTangentCorrection(reciprocalStep);
            end
        
            obj.k = obj.k + 1;
            obj.t = obj.t + obj.dt;
        
            obj.last.Clamp6     = loads.Clamp6;
            obj.last.Fwing_B    = loads.Fwing_B;
            obj.last.Mwing_B    = loads.Mwing_B;
            obj.last.Ftail_B    = loads.Ftail_B;
            obj.last.Mtail_B    = loads.Mtail_B;
            obj.last.Ffin_B     = loads.Ffin_B;
            obj.last.Mfin_B     = loads.Mfin_B;
            obj.last.Fgrav_B    = loads.Fgrav_B;
            obj.last.Fthrust_B  = loads.Fthrust_B;
            obj.last.Mthrust_B  = loads.Mthrust_B;
            obj.last.Freciprocal_B = loads.Freciprocal_B;
            obj.last.Mreciprocal_B = loads.Mreciprocal_B;
            obj.last.Ftot_B     = loads.Ftot_B;
            obj.last.Mtot_B     = loads.Mtot_B;
            if obj.reciprocalTangentEnabled()
                obj.last.reciprocalTangentCorrection = struct( ...
                    'enabled',true, ...
                    'visibleIncrement', ...
                        obj.reciprocalTangentCorrection.lastVisibleCorrection, ...
                    'hiddenState', ...
                        obj.reciprocalTangentCorrection.hiddenState, ...
                    'inputIncrement', ...
                        obj.reciprocalTangentCorrection.lastInput, ...
                    'netBodyWrenchCorrection', ...
                        obj.reciprocalTangentCorrection. ...
                        lastWrenchCorrection, ...
                    'EulerChartCorrection', ...
                        obj.reciprocalTangentCorrection. ...
                        lastRigidChartCorrection, ...
                    'rigidRealizationResidual', ...
                        obj.reciprocalTangentCorrection. ...
                        lastRigidRealizationResidual, ...
                    'rigidWrenchRefinementCount', ...
                        obj.reciprocalTangentCorrection. ...
                        lastWrenchRefinementCount, ...
                    'domainRejected', ...
                        obj.reciprocalTangentCorrection.lastDomainRejected, ...
                    'domainReason', ...
                        obj.reciprocalTangentCorrection.lastDomainReason, ...
                    'domainDetails', ...
                        obj.reciprocalTangentCorrection.lastDomainDetails, ...
                    'stepAccepted',obj.reciprocalTangentCorrection. ...
                        lastStepAccepted, ...
                    'replacesFixedRootIncrement', ...
                        obj.reciprocalTangentCorrection. ...
                        replacesFixedRootIncrement, ...
                    'fixedRootIncrement', ...
                        obj.reciprocalTangentCorrection. ...
                        lastFixedRootIncrement, ...
                    'bodyAppliedWingWrench', ...
                        obj.reciprocalTangentCorrection. ...
                        lastBodyAppliedWingWrench, ...
                    'sourceRootWrench', ...
                        obj.reciprocalTangentCorrection. ...
                        lastSourceRootWrench, ...
                    'loadOutputQualified', ...
                        obj.reciprocalTangentCorrection. ...
                        replacesFixedRootIncrement && ...
                        all(structfun(@logical, ...
                        obj.reciprocalTangentCorrection. ...
                        sourceRootOutput.checks)), ...
                    'clampOwner', ...
                        'fixed-root structural diagnostic only', ...
                    'rootLoadOwner', ...
                        'source wing dynamic equilibrium at root cut', ...
                    'bodyCorrectionOwner', ...
                        'source reciprocal incremental rigid-wing path');
            else
                obj.last.reciprocalTangentCorrection = struct( ...
                    'enabled',false,'loadOutputQualified',false);
            end
        
            if obj.debugCoupled
                fprintf(['[stepCoupled] t=%.4f | U=%.3f | alpha=%+.3f deg | ', ...
                         '|Ftot|=%.3e | |Mtot|=%.3e | Fx=%+.3e Fz=%+.3e My=%+.3e\n'], ...
                        obj.t, obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha), ...
                        norm(loads.Ftot_B), norm(loads.Mtot_B), ...
                        loads.Ftot_B(1), loads.Ftot_B(3), loads.Mtot_B(2));
            end
        
            xNext = obj.packCoupledState();
        end
        % function xNext = stepCoupled(obj, u_cmd, gust, rbCmd)
        % %STEPCOUPLED Advance one coupled rigid-flexible plant step.
        % %
        % % Inputs are increments about trim. Trim wing deflection, elevator, and
        % % thrust are added internally.
        % 
        %     if nargin < 4 || isempty(rbCmd)
        %         rbCmd = struct();
        %     end
        % 
        %     u_cmd = obj.sanitizeControl(u_cmd);
        %     gust  = obj.sanitizeDisturbance(gust);
        % 
        %     obj.flightCond = obj.computeFlightConditionFromRB();
        % 
        %     % Scheduled model update. Use frozen trim mode until t=0 trim replay passes.
        %     obj.updateScheduledROM(gust);
        % 
        %     % Total controls applied to the flexible wing and rigid aircraft.
        %     uWing = obj.composeWingControl(u_cmd);
        %     rbTot = obj.composeRigidCommand(rbCmd);
        % 
        %     obj.lastWingControl = uWing;
        %     obj.last.uWingTotal = uWing;
        %     obj.last.gust       = gust;
        % 
        %     % Loads are evaluated at the beginning of the step. This makes the
        %     % first-step residual directly comparable with the trim residual.
        %     loads = obj.computeCoupledLoads(obj.xFlex, uWing, gust, rbTot);
        % 
        %     xBefore = obj.xFlex;
        %     uBefore = uWing(:);
        % 
        %     if ~isfield(obj.last,'qRatio') || isempty(obj.last.qRatio)
        %         qRatio = 1;
        %     else
        %         qRatio = obj.last.qRatio;
        %     end
        % 
        %     %%% Not sure if we use this anymore since scheduling:
        %     % This might work better? was used in legacy but not sure if
        %     % needed now.
        %     if obj.cfg.library.updateMode ~= "perPlantStep"
        %         obj.updateROMScalingFromRigidBody();
        %         obj.updateWingFlightCondition();
        %     end
        %     %%%
        %     % obj.model.parConst.RateProject = struct( ...
        %     % 'projSet', true, ...
        %     % 'Pz', obj.beam.Pz,'Pr', obj.beam.Pr,'red', obj.beam.red);
        %     [xProbe, ~] = obj.model.step( ...
        %         xBefore, uBefore, gust, [], false, qRatio);
        % 
        %     dxProbe = xProbe - xBefore;
        % 
        %     fprintf(['[flex fixed-point check] |dx|=%.3e |dx/dt|=%.3e ', ...
        %              '|dq1|=%.3e |dq2|=%.3e |dqGam|=%.3e\n'], ...
        %             norm(dxProbe), norm(dxProbe)/obj.model.dt, ...
        %             norm(dxProbe(obj.idx.q1)), ...
        %             norm(dxProbe(obj.idx.q2)), ...
        %             norm(dxProbe(obj.idx.qGam)));
        % 
        %     fprintf(['[wing input] delta=[%+.6e %+.6e] ', ...
        %              'deltaDot=[%+.6e %+.6e] qRatio=%.9f\n'], ...
        %             uBefore(1), uBefore(2), uBefore(3), uBefore(4), qRatio);
        %     % Propagate flexible ROM with the same total wing command used in loads.
        %     obj.xFlex = obj.stepWingROM(uWing, gust);
        % 
        %     % Propagate rigid body with beginning-of-step loads.
        %     obj.rb = obj.stepRigidBody(obj.rb, loads.Ftot_B, loads.Mtot_B);
        % 
        %     obj.k = obj.k + 1;
        %     obj.t = obj.t + obj.dt;
        % 
        %     obj.last.Clamp6     = loads.Clamp6;
        %     obj.last.Fwing_B    = loads.Fwing_B;
        %     obj.last.Mwing_B    = loads.Mwing_B;
        %     obj.last.Ftail_B    = loads.Ftail_B;
        %     obj.last.Mtail_B    = loads.Mtail_B;
        %     obj.last.Ffin_B     = loads.Ffin_B;
        %     obj.last.Mfin_B     = loads.Mfin_B;
        %     obj.last.Fgrav_B    = loads.Fgrav_B;
        %     obj.last.Fthrust_B  = loads.Fthrust_B;
        %     obj.last.Mthrust_B  = loads.Mthrust_B;
        %     obj.last.Ftot_B     = loads.Ftot_B;
        %     obj.last.Mtot_B     = loads.Mtot_B;
        % 
        %     if obj.debugCoupled
        %         fprintf(['[stepCoupled] t=%.4f | U=%.3f | alpha=%+.3f deg | ', ...
        %                  '|Ftot|=%.3e | |Mtot|=%.3e | Fx=%+.3e Fz=%+.3e My=%+.3e\n'], ...
        %                 obj.t, obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha), ...
        %                 norm(loads.Ftot_B), norm(loads.Mtot_B), ...
        %                 loads.Ftot_B(1), loads.Ftot_B(3), loads.Mtot_B(2));
        %     end
        % 
        %     xNext = obj.packCoupledState();
        % end
        % function xNext = stepCoupled(obj, u_cmd, gust, rbCmd)
        % %STEPCOUPLED One explicit rigid-flexible coupled plant step.
        % %
        % % Sequence:
        % %   1. Read rigid-body flight condition.
        % %   2. Update wing flight-condition cache.
        % %   3. Propagate flexible wing ROM by one internal dt.
        % %   4. Compute wing clamp/root reaction.
        % %   5. Compute tail/fin/thrust/gravity loads.
        % %   6. Sum forces/moments.
        % %   7. Propagate rigid-body equations.
        % %   8. Return packed coupled state.
        % %
        % % This is intentionally explicit coupling. The coupling becomes stiff,
        % %% Note to wrap steps 3--7 in a fixed-point subiteration later.
        %      if nargin < 4 || isempty(rbCmd)
        %         % rbCmd = struct('delta_e',0,'delta_a',0,'delta_r',0,'thrust',obj.trimThrust);
        %         rbCmd = [];
        %      end
        %     % -------------------------------------------------------------------------
        %     % Rigid-body command defaults.
        %     % rbCmd fields are interpreted as increments about trim unless otherwise
        %     % stated.
        %     % -------------------------------------------------------------------------
        %     if isempty(rbCmd)
        %         rbCmd = struct();
        %     end
        % 
        %     if ~isfield(rbCmd,'delta_e') || isempty(rbCmd.delta_e)
        %         rbCmd.delta_e = 0;
        %     end
        % 
        %     if ~isfield(rbCmd,'delta_a') || isempty(rbCmd.delta_a)
        %         rbCmd.delta_a = 0;
        %     end
        % 
        %     if ~isfield(rbCmd,'delta_r') || isempty(rbCmd.delta_r)
        %         rbCmd.delta_r = 0;
        %     end
        % 
        %     if ~isfield(rbCmd,'thrust') || isempty(rbCmd.thrust)
        %         rbCmd.thrust = 0;
        %     end
        % 
        %     % Total rear elevator and thrust applied to rigid aircraft.
        %     delta_e_total = obj.trimTailDelta + rbCmd.delta_e;
        %     thrust_total  = obj.trimThrust    + rbCmd.thrust;
        %     rbCmd.delta_e =delta_e_total;
        %     rbCmd.thrust = thrust_total;
        %     % --------------------------------------------------------------
        %     % 0. Validate/sanitize inputs
        %     % --------------------------------------------------------------
        %     u_cmd = obj.sanitizeControl(u_cmd);
        %     gust  = obj.sanitizeDisturbance(gust);
        % 
        %     % --------------------------------------------------------------
        %     % 1. Current rigid-body-derived flight condition
        %     % --------------------------------------------------------------
        %     obj.flightCond = obj.computeFlightConditionFromRB();
        % 
        %     % --------------------------------------------------------------
        %     % 2. Push flight condition into flexible ROM, if requested
        %     % --------------------------------------------------------------
        %     % gustEff = obj.mapRigidMotionToWingInput(gust);
        % 
        %     % Scheduled update based on condition
        %     obj.updateScheduledROM(gust);
        % 
        %     % Old debug scaling
        %     % obj.updateROMScalingFromRigidBody();
        %     % obj.updateWingFlightCondition();
        %     obj.last.qRatio = 1;
        %     % --------------------------------------------------------------
        %     % 3. Propagate flexible wing ROM
        %     % --------------------------------------------------------------
        %     obj.xFlex = obj.stepWingROM(u_cmd, gust);
        % 
        %     % --------------------------------------------------------------
        %     % 4. Compute wing clamp/root reaction
        %     % --------------------------------------------------------------
        %     Clamp6 = obj.computeWingClampReaction();
        %     Clamp6 = obj.mirrorWingClamp(Clamp6);
        % 
        %     Fwing_B = Clamp6(1:3);
        %     Mwing_B = Clamp6(4:6);
        % 
        %     % --------------------------------------------------------------
        %     % 5. Rigid-body aerodynamic / gravity / thrust loads
        %     % --------------------------------------------------------------
        %     if isempty(rbCmd )
        %         [Ftail_B, Mtail_B] = obj.computeTailLoads(u_cmd);
        %         [Ffin_B,  Mfin_B ] = obj.computeFinLoads(u_cmd);
        % 
        %         Fgrav_B = obj.computeGravityBody();
        % 
        %         [Fthrust_B, Mthrust_B] = obj.computeThrustLoads(u_cmd);
        %     else
        %         % Rigid-body tail/fin/thrust see rbCmd:
        %         [Ftail_B, Mtail_B] = obj.computeTailLoadsFromCmd(rbCmd);
        %         [Ffin_B,  Mfin_B ] = obj.computeFinLoadsFromCmd(rbCmd);
        %         Fgrav_B = obj.computeGravityBody();
        % 
        %         [Fthrust_B, Mthrust_B] = obj.computeThrustLoadsFromCmd(rbCmd);
        % 
        %     end
        %     % --------------------------------------------------------------
        %     % 6. Sum total aircraft force/moment in body axes
        %     % --------------------------------------------------------------
        %     Ftot_B = Fwing_B + Ftail_B + Ffin_B + Fgrav_B + Fthrust_B;
        %     Mtot_B = Mwing_B + Mtail_B + Mfin_B + Mthrust_B;
        % 
        %     % --------------------------------------------------------------
        %     % 7. Propagate rigid-body EOM
        %     % --------------------------------------------------------------
        %     obj.rb = obj.stepRigidBody(obj.rb, Ftot_B, Mtot_B);
        % 
        %     % --------------------------------------------------------------
        %     % 8. Advance time and store diagnostics
        %     % --------------------------------------------------------------
        %     obj.k = obj.k + 1;
        %     obj.t = obj.t + obj.dt;
        % 
        %     obj.last.Clamp6 = Clamp6;
        %     obj.last.Fwing_B = Fwing_B;
        %     obj.last.Mwing_B = Mwing_B;
        %     obj.last.Ftail_B = Ftail_B;
        %     obj.last.Mtail_B = Mtail_B;
        %     obj.last.Ffin_B  = Ffin_B;
        %     obj.last.Mfin_B  = Mfin_B;
        %     obj.last.Fgrav_B = Fgrav_B;
        %     obj.last.Fthrust_B = Fthrust_B;
        %     obj.last.Mthrust_B = Mthrust_B;
        %     obj.last.Ftot_B = Ftot_B;
        %     obj.last.Mtot_B = Mtot_B;
        % 
        %     if obj.debugCoupled
        %         fprintf(['[stepCoupled] t=%.4f | U=%.3f | alpha=%.3f deg | ', ...
        %                  '|Fwing|=%.3e | |Ftot|=%.3e | |Mtot|=%.3e\n'], ...
        %                 obj.t, obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha), ...
        %                 norm(Fwing_B), norm(Ftot_B), norm(Mtot_B));
        %     end
        % 
        %     % Return packed coupled state.
        %     xNext = obj.packCoupledState();
        % end

        function updateROMScalingFromRigidBody(obj)
        %UPDATEROMSCALINGFROMRIGIDBODY Apply cheap U-dependent scaling.
        %
        % This is only useful if the ROMIntegrator actually uses these fields during
        % model.step. If L already contains the scaling, this does not update L.
        
            mode = "fixed";
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'rigidToWingMode')
                mode = lower(string(obj.cfg.sim.rigidToWingMode));
            end
        
            if mode ~= "qbar" && mode ~= "alpha_gust"
                return
            end
        
            U0 = obj.cfg.flight.U_inf;
            U  = obj.flightCond.Uinf;
        
            if U0 <= 0
                return
            end
        
            qRatio = (U/U0)^2;
        
            % Store for diagnostics.
            obj.last.qRatio = qRatio;
        
           
        end
        % function attachROMLibrary(obj, ROMlib)
        %     obj.ROMlib = AeroFlex.sched.loadLibrary(ROMlib);
        %     obj.scheduleEnabled = true;
        %     mu0 = [obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha)];
        %     obj.sched = AeroFlex.sched.evalLibrary(obj.ROMlib, mu0, obj.cfg.library);
        %     obj = AeroFlex.sched.applyToPlant(obj, obj.sched);
        % end
        % function attachROMLibrary(obj, ROMlibIn)
        % %ATTACHROMLIBRARY Attach interpolated ROM library.
        % 
        %     if isstruct(ROMlibIn)
        %         obj.ROMlib = ROMlibIn;
        %     else
        %         obj.ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
        %     end
        % 
        %     obj.scheduleEnabled = true;
        % 
        %     obj.trimMu = [obj.safeGetFlightSpeed(obj.cfg), rad2deg(obj.safeGetAlpha(obj.cfg))];
        % 
        %     obj.sched = AeroFlex.sched.evalLibrary(obj.ROMlib, obj.trimMu, obj.cfg.library);
        %     obj.model = AeroFlex.sim.ROMIntegrator(obj.cfg, obj.beam, obj.aero, obj.sched.base);
        %     obj = AeroFlex.sched.applyToPlant(obj, obj.sched);
        % 
        %     % obj.model.parConst.RateProject = struct( ...
        %     % 'projSet', true, ...
        %     % 'Pz', obj.beam.Pz,'Pr', obj.beam.Pr,'red', obj.beam.red);
        % 
        %     % 'Pz', obj.beam.Pz);
        %     % 'Pz', obj.beam.Pz,'Pr', obj.beam.Pr,'red', obj.beam.red);
        %     % obj.model = AeroFlex.sched.applyToROMIntegrator(obj.model, obj.sched, obj.cfg);
        % 
        % 
        %     obj.last.sched_mu = obj.trimMu;
        % 
        %     % if obj.debugCoupled
        %     %     uWing0 = obj.trimWingControl;
        %     %     rb0 = obj.composeRigidCommand(struct());
        %     %     loads0 = obj.computeCoupledLoads(obj.xFlex, uWing0, zeros(obj.nw,1), rb0);
        %     % 
        %     %     fprintf(['[PlantRunTime:initTrimCheck] mu=[%.3f %.3f] | ', ...
        %     %              'Fx=%+.3e Fz=%+.3e My=%+.3e | |F|=%.3e | |M|=%.3e\n'], ...
        %     %             obj.trimMu(1), obj.trimMu(2), ...
        %     %             loads0.Ftot_B(1), loads0.Ftot_B(3), loads0.Mtot_B(2), ...
        %     %             norm(loads0.Ftot_B), norm(loads0.Mtot_B));
        %     % end
        % 
        %     if obj.debugCoupled
        %         uWing0 = obj.trimWingControl;
        %         rb0 = obj.composeRigidCommand(struct());
        %         loads0 = obj.computeCoupledLoads(obj.xFlex, uWing0, zeros(obj.nw,1), rb0);
        % 
        %         fprintf(['[PlantRunTime:initTrimCheck] mu=[%.3f %.3f] | ', ...
        %                  'Fx=%+.3e Fz=%+.3e My=%+.3e | |F|=%.3e | |M|=%.3e\n'], ...
        %                 obj.trimMu(1), obj.trimMu(2), ...
        %                 loads0.Ftot_B(1), loads0.Ftot_B(3), loads0.Mtot_B(2), ...
        %                 norm(loads0.Ftot_B), norm(loads0.Mtot_B));
        % 
        %         obj.printTrimReplayComparison(loads0);
        %     end
        % end
        function attachROMLibrary(obj, ROMlibIn)
        %ATTACHROMLIBRARY Attach a Stage-2 compatible ROM library.
        %
        % The library is applied as a complete scheduled package.  Do not rebuild a
        % nominal ROM with only sched.base; that mixes scheduled and unscheduled data.
        
            if isstruct(ROMlibIn)
                obj.ROMlib = ROMlibIn;
            else
                obj.ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
            end
        
            if ~isfield(obj.ROMlib,'compatibleCoordinates') || ~obj.ROMlib.compatibleCoordinates
                error('PlantRunTime:RawLibrary', ...
                    'PlantRunTime requires the Stage-2 compatible library.');
            end
        
            obj.scheduleEnabled = true;
        
            physicalTrimMu = [obj.safeGetFlightSpeed(obj.cfg), ...
                              rad2deg(obj.safeGetAlpha(obj.cfg))];
            obj.trimMu = physicalTrimMu;
            freezePackage = isfield(obj.cfg,'library') && ...
                isfield(obj.cfg.library,'updateMode') && ...
                strcmpi(string(obj.cfg.library.updateMode), "frozenTrim");
            if freezePackage && isfield(obj.cfg.library,'schedulerQueryPoint') && ...
                    numel(obj.cfg.library.schedulerQueryPoint) == 2
                obj.trimMu = double(obj.cfg.library.schedulerQueryPoint(:).');
            end
        
            obj.sched = AeroFlex.sched.evalLibrary( ...
                obj.ROMlib, obj.trimMu, obj.cfg.library);

            requireExact = isfield(obj.cfg.library,'requireExactNode') && ...
                logical(obj.cfg.library.requireExactNode);
            if requireExact
                tol = 1e-10;
                if isfield(obj.cfg.library,'interpTol')
                    tol = obj.cfg.library.interpTol;
                end
                isExact = numel(obj.sched.pointIds) == 1 && ...
                    numel(obj.sched.weights) == 1 && ...
                    abs(obj.sched.weights(1) - 1) <= tol && ...
                    norm(obj.sched.pointMu(1,:) - obj.trimMu, inf) <= tol;
                if ~isExact
                    error('PlantRunTime:ExactNodeRequired', ...
                        ['Frozen runtime requires one exact source node at ', ...
                         'U=%.12g m/s, alpha=%.12g deg.'], ...
                        obj.trimMu(1), obj.trimMu(2));
                end
            end
        
            % Apply the complete scheduled ROM.  Leave xFlex untouched.
            obj = AeroFlex.sched.applyToPlant(obj, obj.sched);
        
            % The sensor basis follows the active/reference beam basis.
            obj.sensor = AeroFlex.sensor.WingVelSensor(obj.beam, obj.cfg);
        
            obj.last.sched_mu = obj.trimMu;
            obj.last.qRatio   = 1;
            obj.applyPackageRelativeRigidAttitudeChi();
        
            if obj.debugCoupled
                uWing0 = obj.trimWingControl;
                rb0 = obj.composeRigidCommand(struct());
                loads0 = obj.computeCoupledLoads(obj.xFlex, uWing0, zeros(obj.nw,1), rb0);
        
                fprintf(['[PlantRunTime:initTrimCheck] mu=[%.3f %.3f] | ', ...
                         'Fx=%+.3e Fz=%+.3e My=%+.3e | |F|=%.3e | |M|=%.3e\n'], ...
                        obj.trimMu(1), obj.trimMu(2), ...
                        loads0.Ftot_B(1), loads0.Ftot_B(3), loads0.Mtot_B(2), ...
                        norm(loads0.Ftot_B), norm(loads0.Mtot_B));
        
                obj.printTrimReplayComparison(loads0);
            end
        end

        function [package,info] = prepareControlSample( ...
                obj,gust,uPreviousInterval,rbCmd,plantStepCount)
        %PREPARECONTROLSAMPLE Select and hold one model for a control sample.
        %
        % The scheduler uses the current plant state and the command applied
        % over the preceding interval. The accepted package is then held for
        % plantStepCount calls to stepCoupled so the plant, estimator, and
        % controller use one interval map.

            if nargin < 2 || isempty(gust)
                gust = zeros(obj.nw,1);
            end
            if nargin < 3 || isempty(uPreviousInterval)
                uPreviousInterval = zeros(obj.nu,1);
            end
            if nargin < 4 || isempty(rbCmd)
                rbCmd = struct();
            end
            if nargin < 5 || isempty(plantStepCount)
                plantStepCount = 1;
            end

            assert(isscalar(plantStepCount) && isfinite(plantStepCount) && ...
                plantStepCount >= 1 && plantStepCount == round(plantStepCount), ...
                'PlantRunTime:ControlSampleStepCount', ...
                'plantStepCount must be a positive integer.');
            assert(obj.controlSampleScheduleHoldStepsRemaining == 0, ...
                'PlantRunTime:ControlSampleOverlap', ...
                ['A new control sample was prepared before %d held plant ', ...
                 'steps from the preceding sample were consumed.'], ...
                obj.controlSampleScheduleHoldStepsRemaining);

            gust = obj.sanitizeDisturbance(gust);
            uPreviousInterval = obj.sanitizeControl(uPreviousInterval);
            obj.flightCond = obj.computeFlightConditionFromRB();
            obj.updateWingFlightCondition();
            uWing = obj.composeWingControl(uPreviousInterval);
            rbTot = obj.composeRigidCommand(rbCmd);

            muBefore = obj.last.sched_mu;
            obj.last.scheduleStateTransport = struct('enabled',false);
            obj.updateScheduledROM(gust,uWing,rbTot);
            obj.applyPackageRelativeRigidAttitudeChi();
            package = obj.activePredictionPackage();

            if obj.useRigidBody && obj.scheduleEnabled
                obj.controlSampleScheduleHoldStepsRemaining = plantStepCount;
            end
            obj.controlSampleIndex = obj.controlSampleIndex + 1;

            info = struct();
            info.controlSampleIndex = obj.controlSampleIndex;
            info.plantStepCount = plantStepCount;
            info.activeMu = package.mu(:).';
            info.requestedScheduleUInfMps = nan;
            if isfield(rbCmd,'scheduleUInfMps') && ...
                    ~isempty(rbCmd.scheduleUInfMps)
                info.requestedScheduleUInfMps = ...
                    double(rbCmd.scheduleUInfMps);
            end
            info.composedScheduleUInfMps = nan;
            if isfield(rbTot,'scheduleUInfMps') && ...
                    ~isempty(rbTot.scheduleUInfMps)
                info.composedScheduleUInfMps = ...
                    double(rbTot.scheduleUInfMps);
            end
            info.scheduleChanged = numel(muBefore) ~= numel(info.activeMu) || ...
                norm(muBefore(:)-info.activeMu(:),inf) > 0;
            info.scheduleMuRaw = nan(1,2);
            if isfield(obj.last,'sched_mu_raw') && ...
                    numel(obj.last.sched_mu_raw) == 2
                info.scheduleMuRaw = double(obj.last.sched_mu_raw(:).');
            end
            info.scheduleMuFiltered = nan(1,2);
            if isfield(obj.last,'sched_mu_filter') && ...
                    numel(obj.last.sched_mu_filter) == 2
                info.scheduleMuFiltered = ...
                    double(obj.last.sched_mu_filter(:).');
            end
            info.schedulerGate = struct();
            if isfield(obj.last,'schedGate') && isstruct(obj.last.schedGate)
                info.schedulerGate = obj.last.schedGate;
            end
            info.stateTransport = obj.last.scheduleStateTransport;
            info.holdStepsRemaining = ...
                obj.controlSampleScheduleHoldStepsRemaining;
        end

        function package = activePredictionPackage(obj)
        %ACTIVEPREDICTIONPACKAGE Return the exact active runtime model bundle.

            if isstruct(obj.sched) && ~isempty(obj.sched)
                package = obj.sched;
            else
                package = struct();
            end

            package.L = obj.model.L;
            package.idx = obj.idx;
            package.parConst = obj.model.parConst;
            package.beam = obj.beam;
            package.aero = obj.aero;
            package.base = obj.base;
            package.mu = obj.last.sched_mu(:).';
            if isprop(obj.model,'internalCoupledCoordinate')
                package.internalCoupledCoordinate = ...
                    obj.model.internalCoupledCoordinate;
            end
        end

        function [xNext,interval] = stepCoupledIntervalAudit( ...
                obj,uCommand,gust,rbCommand,packet)
        %STEPCOUPLEDINTERVALAUDIT Commit one qualified packed P1 interval.
        %   The pure packed core completes before this handle is mutated.
        %   Omission of the guarded request leaves stepCoupled unchanged.

            request = obj.getNestedMember(obj.cfg, ...
                {'sim','scheduledReciprocalPlantIntervalAudit'},struct());
            assert(isstruct(request) && isscalar(request) && ...
                isfield(request,'enabled') && logical(request.enabled) && ...
                isfield(request,'auditOnly') && logical(request.auditOnly) && ...
                isfield(request,'changeId') && string(request.changeId)== ...
                "phase18c-v17a-casebc-p1-plant-interval-runtime-integration-v1" && ...
                isfield(request,'caseId') && ismember(string(request.caseId), ...
                ["formal_case_b","formal_case_c"]), ...
                'PlantRunTime:P1Request', ...
                'The packed plant interval requires its approved Case-B/C request.');
            assert(isstruct(packet) && isscalar(packet) && ...
                isfield(packet,'schemaVersion') && string(packet.schemaVersion)== ...
                "phase18c-v17a-casebc-p1-plant-interval-packet-v1" && ...
                isfield(packet,'changeId') && ...
                string(packet.changeId)==string(request.changeId), ...
                'PlantRunTime:P1Packet', ...
                'The packed plant interval packet is absent or incompatible.');

            uCommand = obj.sanitizeControl(uCommand);
            gust = obj.sanitizeDisturbance(gust);
            if nargin<4 || isempty(rbCommand)
                rbCommand = struct();
            end
            rbTotal = obj.composeRigidCommand(rbCommand);
            deltaAileron = 0;
            deltaRudder = 0;
            elevatorIncrement = 0;
            thrustIncrement = 0;
            if isfield(rbCommand,'delta_a') && ~isempty(rbCommand.delta_a)
                deltaAileron = double(rbCommand.delta_a);
            end
            if isfield(rbCommand,'delta_r') && ~isempty(rbCommand.delta_r)
                deltaRudder = double(rbCommand.delta_r);
            end
            if isfield(rbCommand,'delta_e') && ~isempty(rbCommand.delta_e)
                elevatorIncrement = double(rbCommand.delta_e);
            end
            if isfield(rbCommand,'thrust') && ~isempty(rbCommand.thrust)
                thrustIncrement = double(rbCommand.thrust);
            end
            assert(numel(gust)==1 && numel(uCommand)==4 && ...
                isscalar(deltaAileron) && isfinite(deltaAileron) && ...
                deltaAileron==0 && isscalar(deltaRudder) && ...
                isfinite(deltaRudder) && deltaRudder==0 && ...
                isscalar(elevatorIncrement) && isfinite(elevatorIncrement) && ...
                isscalar(thrustIncrement) && isfinite(thrustIncrement) && ...
                rbTotal.thrust>=0, ...
                'PlantRunTime:P1SymmetricScope', ...
                ['P1 requires scalar gust, four symmetric-wing channels, ', ...
                 'zero lateral commands, and nonnegative total thrust.']);
            substepCount = double(packet.substepCount);
            assert(substepCount==14 && ...
                obj.controlSampleScheduleHoldStepsRemaining==substepCount && ...
                obj.reciprocalTangentEnabled() && ...
                numel(obj.reciprocalTangentCorrection.hiddenState)==546, ...
                'PlantRunTime:P1RuntimeOwner', ...
                'The qualified fourteen-substep scheduled owner is not active.');

            rigidInitial = [obj.rb.r_I(:);obj.rb.v_B(:); ...
                obj.rb.euler(:);obj.rb.omega_B(:)];
            rigidIncrement = [elevatorIncrement;thrustIncrement];
            interval = AeroFlex.sim. ...
                scheduledReciprocalPlantIntervalCoreAudit( ...
                obj.xFlex(:),rigidInitial, ...
                obj.reciprocalTangentCorrection.hiddenState(:), ...
                uCommand,gust,rigidIncrement,packet);

            expectedStateSize = [632,substepCount+1];
            assert(isequal(size(interval.stateHistory),expectedStateSize) && ...
                isequal(size(interval.rigidStateBeforeHistory), ...
                    [12,substepCount]) && ...
                isequal(size(interval.sourceRatios),[14,4,substepCount]) && ...
                numel(interval.domainDetails)==substepCount && ...
                all(isfinite(interval.stateHistory),'all') && ...
                all(isfinite(interval.totalWrenchHistory),'all') && ...
                all(isfinite(interval.sourceRootHistory),'all') && ...
                all(isfinite(interval.sourceRatios),'all') && ...
                interval.thrustTotal>=0 && ...
                norm(interval.stateHistory(:,1)- ...
                    [obj.xFlex(:);rigidInitial; ...
                    obj.reciprocalTangentCorrection.hiddenState(:)],inf)==0, ...
                'PlantRunTime:P1Output', ...
                'The packed plant interval failed its precommit contract.');

            oldStep = obj.k;
            oldTime = obj.t;
            finalState = interval.stateHistory(:,end);
            obj.xFlex = finalState(1:74);
            obj.rb.r_I = finalState(75:77);
            obj.rb.v_B = finalState(78:80);
            obj.rb.euler = finalState(81:83);
            obj.rb.omega_B = finalState(84:86);
            obj.reciprocalTangentCorrection.hiddenState = finalState(87:end);
            if numel(obj.x)>=obj.nx
                obj.x(1:obj.nx) = obj.xFlex;
            end
            obj.k = oldStep+substepCount;
            obj.t = oldTime+substepCount*obj.dt;
            obj.controlSampleScheduleHoldStepsRemaining = ...
                obj.controlSampleScheduleHoldStepsRemaining-substepCount;

            finalCondition = interval.flightConditionHistory(:,end);
            obj.flightCond = struct('Uinf',finalCondition(1), ...
                'alpha',finalCondition(2),'beta',finalCondition(3), ...
                'euler',finalCondition(4:6),'omega_B',finalCondition(7:9));
            obj.lastWingControl = interval.wingTotal(:);
            obj.last.uWingTotal = interval.wingTotal(:);
            obj.last.gust = gust(:);

            finalIndex = substepCount;
            obj.last.Clamp6 = interval.clampHistory(:,finalIndex);
            obj.last.Fwing_B = interval.wingWrenchHistory(1:3,finalIndex);
            obj.last.Mwing_B = interval.wingWrenchHistory(4:6,finalIndex);
            obj.last.Ftail_B = interval.tailWrenchHistory(1:3,finalIndex);
            obj.last.Mtail_B = interval.tailWrenchHistory(4:6,finalIndex);
            obj.last.Ffin_B = interval.finWrenchHistory(1:3,finalIndex);
            obj.last.Mfin_B = interval.finWrenchHistory(4:6,finalIndex);
            obj.last.Fgrav_B = interval.gravityWrenchHistory(1:3,finalIndex);
            obj.last.Fthrust_B = interval.thrustWrenchHistory(1:3,finalIndex);
            obj.last.Mthrust_B = interval.thrustWrenchHistory(4:6,finalIndex);
            obj.last.Freciprocal_B = ...
                interval.reciprocalWrenchHistory(1:3,finalIndex);
            obj.last.Mreciprocal_B = ...
                interval.reciprocalWrenchHistory(4:6,finalIndex);
            obj.last.Ftot_B = interval.totalWrenchHistory(1:3,finalIndex);
            obj.last.Mtot_B = interval.totalWrenchHistory(4:6,finalIndex);

            runtime = obj.reciprocalTangentCorrection;
            runtime.lastVisibleCorrection = ...
                interval.visibleCorrectionHistory(:,finalIndex);
            runtime.lastInput = interval.inputDeparture(:);
            runtime.lastWrenchCorrection = ...
                interval.wrenchCorrectionHistory(:,finalIndex);
            runtime.lastRigidChartCorrection = ...
                interval.rigidChartCorrectionHistory(:,finalIndex);
            runtime.lastRigidRealizationResidual = ...
                interval.rigidResidual(finalIndex);
            runtime.lastStepAccepted = logical( ...
                interval.stepAccepted(finalIndex));
            runtime.lastDomainRejected = logical( ...
                interval.domainRejected(finalIndex));
            runtime.lastDomainReason = string( ...
                interval.domainReason(finalIndex));
            runtime.lastDomainDetails = interval.domainDetails{finalIndex};
            runtime.lastWrenchRefinementCount = ...
                interval.refinementCount(finalIndex);
            runtime.lastFixedRootIncrement = ...
                interval.fixedRootIncrementHistory(:,finalIndex);
            runtime.lastBodyAppliedWingWrench = ...
                interval.bodyAppliedWingWrenchHistory(:,finalIndex);
            runtime.lastSourceRootWrench = ...
                interval.sourceRootHistory(:,finalIndex);
            obj.reciprocalTangentCorrection = runtime;
            obj.last.reciprocalTangentCorrection = struct( ...
                'enabled',true, ...
                'visibleIncrement',runtime.lastVisibleCorrection, ...
                'hiddenState',runtime.hiddenState, ...
                'inputIncrement',runtime.lastInput, ...
                'netBodyWrenchCorrection',runtime.lastWrenchCorrection, ...
                'EulerChartCorrection',runtime.lastRigidChartCorrection, ...
                'rigidRealizationResidual', ...
                    runtime.lastRigidRealizationResidual, ...
                'rigidWrenchRefinementCount', ...
                    runtime.lastWrenchRefinementCount, ...
                'domainRejected',runtime.lastDomainRejected, ...
                'domainReason',runtime.lastDomainReason, ...
                'domainDetails',runtime.lastDomainDetails, ...
                'stepAccepted',runtime.lastStepAccepted, ...
                'replacesFixedRootIncrement', ...
                    runtime.replacesFixedRootIncrement, ...
                'fixedRootIncrement',runtime.lastFixedRootIncrement, ...
                'bodyAppliedWingWrench',runtime.lastBodyAppliedWingWrench, ...
                'sourceRootWrench',runtime.lastSourceRootWrench, ...
                'loadOutputQualified',runtime.replacesFixedRootIncrement && ...
                    all(structfun(@logical,runtime.sourceRootOutput.checks)), ...
                'clampOwner','fixed-root structural diagnostic only', ...
                'rootLoadOwner', ...
                    'source wing dynamic equilibrium at root cut', ...
                'bodyCorrectionOwner', ...
                    'source reciprocal incremental rigid-wing path');

            traceEnabled = isfield(obj.cfg,'debug') && ...
                isfield(obj.cfg.debug,'schedule') && ...
                logical(obj.cfg.debug.schedule);
            if traceEnabled
                obj.last.scheduleUpdateTrace = struct( ...
                    'reason',"control_sample_hold", ...
                    'step',oldStep+substepCount-1, ...
                    'timeSeconds',oldTime+(substepCount-1)*obj.dt, ...
                    'holdStepsRemaining',1);
            end
            if obj.debugCoupled
                fprintf(['[stepCoupledIntervalAudit] t=%.4f | U=%.3f | ', ...
                    'alpha=%+.3f deg | |Ftot|=%.3e | |Mtot|=%.3e\n'], ...
                    obj.t,obj.flightCond.Uinf,rad2deg(obj.flightCond.alpha), ...
                    norm(obj.last.Ftot_B),norm(obj.last.Mtot_B));
            end
            xNext = obj.packCoupledState();
        end


        function diag = diagnoseTrimLocalDerivatives(obj, opts)
        %DIAGNOSETRIMLOCALDERIVATIVES Finite-difference load derivatives at trim.
        %
        % This routine is intentionally read-only with respect to time integration:
        %   - it does not call stepCoupled
        %   - it does not call stepWingROM
        %   - it does not update scheduling
        %   - it restores rb, flightCond, and last before returning
        %
        % Perturbations:
        %   alpha    : changes body velocity direction at fixed Euler attitude
        %   q        : changes body pitch rate only
        %   U        : changes speed at fixed alpha
        %   elevator : changes rear elevator command increment about trim
        %   wing     : changes both wing control surfaces equally about trim
        
        if nargin < 2 || isempty(opts)
            opts = struct();
        end
        
        % if obj.cfg.sim.body_case == "wingOnly" || obj.cfg.sim.bodyCase == "wingOnly"
        %     return;
        % end
        hAlpha = getOpt(opts, 'hAlpha', deg2rad(0.05));
        hQ     = getOpt(opts, 'hQ',     0.01);
        hU     = getOpt(opts, 'hU',     0.05);
        hElev  = getOpt(opts, 'hElev',  deg2rad(0.05));
        hWing  = getOpt(opts, 'hWing',  deg2rad(0.05));
        
        scales = [1.0, 0.5, 0.25];
        
        rbSave     = obj.rb;
        fcSave     = obj.flightCond;
        lastSave   = obj.last;
        xFlexSave  = obj.xFlex;
        
        rb0 = obj.rb;
        
        U0 = norm(rb0.v_B(:));
        if U0 < 1e-9
            U0 = obj.safeGetFlightSpeed(obj.cfg);
        end
        
        alpha0 = atan2(rb0.v_B(3), rb0.v_B(1));
        
        q0 = 0;
        if isfield(rb0,'omega_B') && numel(rb0.omega_B) >= 2
            q0 = rb0.omega_B(2);
        end
        
        gust0 = zeros(obj.nw,1);
        
        diag = struct();
        diag.base.U = U0;
        diag.base.alpha = alpha0;
        diag.base.q = q0;
        diag.steps = struct('hAlpha',hAlpha,'hQ',hQ,'hU',hU, ...
                            'hElev',hElev,'hWing',hWing);
        
        labels = { ...
            'MyWing'; ...
            'MyTail'; ...
            'MyTotal'; ...
            'FzWing'; ...
            'FzTail'; ...
            'FzTotal'; ...
            'FxTotal'};
        
        vars = {'alpha','q','U','elevator','wing'};
        units = {'per rad','per rad/s','per m/s','per rad','per rad'};
        
        D = nan(numel(labels), numel(vars), numel(scales));
        
        try
            baseLoads = evalLoads(U0, alpha0, q0, 0, 0);
            diag.base.loads = baseLoads;
        
            for is = 1:numel(scales)
                s = scales(is);
        
                % alpha derivative: perturb velocity direction at fixed attitude.
                D(:,1,is) = centralDiff('alpha', s*hAlpha);
        
                % pitch-rate derivative: perturb q only.
                D(:,2,is) = centralDiff('q', s*hQ);
        
                % speed derivative: perturb U at fixed alpha.
                D(:,3,is) = centralDiff('U', s*hU);
        
                % elevator derivative: perturb rear elevator command increment.
                D(:,4,is) = centralDiff('elevator', s*hElev);
        
                % wing derivative: perturb both wing surfaces symmetrically.
                D(:,5,is) = centralDiff('wing', s*hWing);
            end
        
        catch ME
            obj.rb = rbSave;
            obj.flightCond = fcSave;
            obj.last = lastSave;
            obj.xFlex = xFlexSave;
            rethrow(ME);
        end
        
        obj.rb = rbSave;
        obj.flightCond = fcSave;
        obj.last = lastSave;
        obj.xFlex = xFlexSave;
        
        diag.labels = labels;
        diag.vars = vars;
        diag.units = units;
        diag.D = D;
        
        fprintf('\n[PlantRunTime: local trim load derivatives]\n');
        fprintf('  Base: U = %.6f m/s, alpha = %.6f deg, q = %.6e rad/s\n', ...
            U0, rad2deg(alpha0), q0);
        fprintf('  Convention check: for the current z-down/body-q convention, a statically restoring aircraft should normally have dMyTotal/dalpha < 0.\n');
        fprintf('                    pitch damping should normally give dMyTotal/dq < 0.\n\n');
        
        for iv = 1:numel(vars)
            fprintf('  Derivatives with respect to %s [%s]\n', vars{iv}, units{iv});
            fprintf('    %-10s %14s %14s %14s\n', 'quantity', 'h', 'h/2', 'h/4');
        
            for il = 1:numel(labels)
                fprintf('    %-10s %+14.6e %+14.6e %+14.6e\n', ...
                    labels{il}, D(il,iv,1), D(il,iv,2), D(il,iv,3));
            end
            fprintf('\n');
        end
        
        diag.summary = struct();
        diag.summary.dMyWing_dalpha  = D(1,1,end);
        diag.summary.dMyTail_dalpha  = D(2,1,end);
        diag.summary.dMyTotal_dalpha = D(3,1,end);
        diag.summary.dMyWing_dq      = D(1,2,end);
        diag.summary.dMyTail_dq      = D(2,2,end);
        diag.summary.dMyTotal_dq     = D(3,2,end);
        
        fprintf('  Summary using h/4:\n');
        fprintf('    dMyWing/dalpha  = %+ .6e N*m/rad\n', diag.summary.dMyWing_dalpha);
        fprintf('    dMyTail/dalpha  = %+ .6e N*m/rad\n', diag.summary.dMyTail_dalpha);
        fprintf('    dMyTotal/dalpha = %+ .6e N*m/rad\n', diag.summary.dMyTotal_dalpha);
        fprintf('    dMyWing/dq      = %+ .6e N*m/(rad/s)\n', diag.summary.dMyWing_dq);
        fprintf('    dMyTail/dq      = %+ .6e N*m/(rad/s)\n', diag.summary.dMyTail_dq);
        fprintf('    dMyTotal/dq     = %+ .6e N*m/(rad/s)\n\n', diag.summary.dMyTotal_dq);
        
            function d = centralDiff(name, h)
                switch name
                    case 'alpha'
                        lp = evalLoads(U0, alpha0 + h, q0, 0, 0);
                        lm = evalLoads(U0, alpha0 - h, q0, 0, 0);
        
                    case 'q'
                        lp = evalLoads(U0, alpha0, q0 + h, 0, 0);
                        lm = evalLoads(U0, alpha0, q0 - h, 0, 0);
        
                    case 'U'
                        lp = evalLoads(U0 + h, alpha0, q0, 0, 0);
                        lm = evalLoads(U0 - h, alpha0, q0, 0, 0);
        
                    case 'elevator'
                        lp = evalLoads(U0, alpha0, q0, +h, 0);
                        lm = evalLoads(U0, alpha0, q0, -h, 0);
        
                    case 'wing'
                        lp = evalLoads(U0, alpha0, q0, 0, +h);
                        lm = evalLoads(U0, alpha0, q0, 0, -h);
        
                    otherwise
                        error('Unknown derivative variable "%s".', name);
                end
        
                d = (extractLoadVector(lp) - extractLoadVector(lm))/(2*h);
            end
        
            function loads = evalLoads(U, alpha, qPitch, elevInc, wingInc)
                rb = rb0;
        
                rb.v_B = [U*cos(alpha); 0; U*sin(alpha)];
        
                if ~isfield(rb,'omega_B') || isempty(rb.omega_B)
                    rb.omega_B = zeros(3,1);
                end
                rb.omega_B = rb.omega_B(:);
                rb.omega_B(2) = qPitch;
        
                obj.rb = rb;
                obj.flightCond = obj.computeFlightConditionFromRB();
        
                uInc = zeros(obj.nu,1);
                if obj.nu >= 1
                    uInc(1) = wingInc;
                end
                if obj.nu >= 2
                    uInc(2) = wingInc;
                end
        
                uWing = obj.composeWingControl(uInc);
        
                rbCmd = struct();
                rbCmd.delta_e = elevInc;
                rbCmd.delta_a = 0;
                rbCmd.delta_r = 0;
                rbCmd.thrust  = 0;
        
                rbTot = obj.composeRigidCommand(rbCmd);
        
                loads = obj.computeCoupledLoads(xFlexSave, uWing, gust0, rbTot);
            end
        
            function y = extractLoadVector(loads)
                y = [ ...
                    loads.Mwing_B(2); ...
                    loads.Mtail_B(2); ...
                    loads.Mtot_B(2); ...
                    loads.Fwing_B(3); ...
                    loads.Ftail_B(3); ...
                    loads.Ftot_B(3); ...
                    loads.Ftot_B(1)];
            end
        
            function val = getOpt(S, name, defaultVal)
                if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
                    val = S.(name);
                else
                    val = defaultVal;
                end
            end
        end
        function debugScheduleConsistency(obj, tag, sched)
        %DEBUGSCHEDULECONSISTENCY Print schedule/update consistency checks.
        
            if nargin < 3 || isempty(sched)
                sched = obj.sched;
            end
        
            fprintf('\n[sched:%s]\n', tag);
        
            if isfield(sched,'mu')
                fprintf('  mu = [%.6f %.6f]\n', sched.mu(1), sched.mu(2));
            end
        
            if isfield(sched,'weights')
                fprintf('  ids = [%s], w = [%s]\n', ...
                    num2str(sched.pointIds(:).'), num2str(sched.weights(:).',' %.4f'));
            end
            if isfield(sched,'x_eq') && ~isempty(sched.x_eq) && numel(sched.x_eq) == numel(obj.xFlex)
                fprintf('  |xFlex - sched.x_eq| = %.3e\n', norm(obj.xFlex - sched.x_eq(:)));
            end
            
            % if isfield(sched,'u_eq') && ~isempty(sched.u_eq)
            %     fprintf('Size obj.lastWingControl(:): %.3f. Size sched.u_eq(:): %.3f. Using first two.' ...
            %         , size(obj.lastWingControl(:), size(sched.u_eq(:))));
            %     fprintf('  |uWing - sched.u_eq| = %.3e\n', norm(obj.lastWingControl(1:2) - sched.u_eq(1:2)));
            % end
            if isfield(sched,'u_eq') && ~isempty(sched.u_eq)
                nU = min(numel(obj.lastWingControl), numel(sched.u_eq));
                fprintf('  size(uWing)=%d, size(u_eq)=%d\n', ...
                    numel(obj.lastWingControl), numel(sched.u_eq));
                fprintf('  |uWing - sched.u_eq| = %.3e\n', ...
                    norm(obj.lastWingControl(1:nU) - sched.u_eq(1:nU)));
            end
            if isfield(sched,'info') && isfield(sched.info,'mode')
                fprintf('  mode = %s\n', string(sched.info.mode));
            end
        
            if ~isempty(obj.model) && isfield(sched,'L')
                fprintf('  ||model.L - sched.L||_F = %.3e\n', ...
                    norm(obj.model.L - sched.L,'fro'));
            end
        
            rpOn = false;
            if isfield(obj.model.parConst,'RateProject') && ...
                    isfield(obj.model.parConst.RateProject,'projSet')
                rpOn = logical(obj.model.parConst.RateProject.projSet);
            end
        
            fprintf('  RateProject = %d\n', rpOn);
        
            if rpOn && isfield(obj.model.parConst.RateProject,'Pz') && ...
                    ~isempty(obj.model.parConst.RateProject.Pz) && ...
                    ~isempty(obj.model.Ldyn)
        
                q1 = obj.idx.q1(:);
                Pz = obj.model.parConst.RateProject.Pz;
        
                eProj = norm(obj.model.Ldyn(q1,:) - Pz*obj.model.L(q1,:),'fro');
                fprintf('  ||Ldyn(q1,:) - Pz*L(q1,:)||_F = %.3e\n', eProj);
            elseif ~isempty(obj.model.Ldyn)
                fprintf('  ||Ldyn - L||_F = %.3e\n', norm(obj.model.Ldyn - obj.model.L,'fro'));
            end
        
            if isfield(sched,'beam')
                if isfield(sched.beam,'Pz') && isfield(obj.beam,'Pz')
                    fprintf('  ||sched.Pz - plant.Pz||_F = %.3e\n', ...
                        norm(sched.beam.Pz - obj.beam.Pz,'fro'));
                end
        
                if isfield(sched.beam,'Pr') && isfield(obj.beam,'Pr')
                    fprintf('  ||sched.Pr - plant.Pr||_F = %.3e\n', ...
                        norm(sched.beam.Pr - obj.beam.Pr,'fro'));
                end
        
                if isfield(sched.beam,'red') && isfield(sched.beam.red,'phi1_sA') && ...
                        isfield(obj.beam,'red') && isfield(obj.beam.red,'phi1_sA')
                    fprintf('  ||sched.phi1_sA - plant.phi1_sA||_F = %.3e\n', ...
                        norm(sched.beam.red.phi1_sA - obj.beam.red.phi1_sA,'fro'));
                end
            end
        
            try
                pc = obj.model.parConst;
                pc.u_ctrl = obj.lastWingControl(:);
                pc.gust = zeros(obj.nw,1);
                pc.N_Thrust = zeros(obj.beam.Nm,1);
        
                raw = AeroFlex.sim.nonlinear_terms(obj.xFlex,pc,obj.idx) + obj.model.L*obj.xFlex;
                q1raw = raw(obj.idx.q1);
        
                Pz = obj.beam.Pz;
                Pr = obj.beam.Pr;
        
                fprintf('  |raw xdot| = %.3e\n', norm(raw));
                fprintf('  |q1raw|    = %.3e\n', norm(q1raw));
                fprintf('  |Pz*q1raw| = %.3e\n', norm(Pz*q1raw));
                fprintf('  |Pr*q1raw| = %.3e\n', norm(Pr*q1raw));
        
                Clamp6 = obj.beam.red.phi1_sA*(Pr*q1raw);
                fprintf('  |Clamp6 raw| = %.3e, Fx=%+.3e, Fz=%+.3e, My=%+.3e\n', ...
                    norm(Clamp6), Clamp6(1), Clamp6(3), Clamp6(5));
            catch ME
                fprintf('  residual check skipped: %s\n', ME.message);
            end
        
            fprintf('\n');
        end
    end
    
    methods (Access = private)
        function enabled = reciprocalTangentEnabled(obj)
            enabled = isstruct(obj.reciprocalTangentCorrection) && ...
                isfield(obj.reciprocalTangentCorrection,'enabled') && ...
                logical(obj.reciprocalTangentCorrection.enabled);
        end

        function [visible,input] = reciprocalTangentInputs( ...
                obj,uCommand,gust,rbCommand)
            candidate = obj.reciprocalTangentCorrection.candidate;
            assert(numel(gust)==1 && numel(uCommand)==4, ...
                'PlantRunTime:ReciprocalTangentInputOrder', ...
                ['Expected gust plus left/right deflection and ', ...
                 'left/right rate.']);
            elevatorIncrement = 0;
            thrustIncrement = 0;
            if isfield(rbCommand,'delta_e') && ~isempty(rbCommand.delta_e)
                elevatorIncrement = double(rbCommand.delta_e);
            end
            if isfield(rbCommand,'thrust') && ~isempty(rbCommand.thrust)
                thrustIncrement = double(rbCommand.thrust);
            end
            assert(isscalar(elevatorIncrement) && ...
                isfinite(elevatorIncrement) && isscalar(thrustIncrement) && ...
                isfinite(thrustIncrement) && ...
                obj.trimThrust+thrustIncrement>=0, ...
                'PlantRunTime:ReciprocalTangentRigidInput', ...
                ['Elevator and thrust increments must be finite scalars, ', ...
                 'and total thrust must remain nonnegative.']);
            sourceVisible = [obj.xFlex(obj.idx.q1);obj.xFlex(obj.idx.q2); ...
                obj.xFlex(obj.idx.qGam);obj.rb.v_B(:); ...
                obj.rb.omega_B(:);obj.rb.euler(:)];
            visible = sourceVisible;
            if obj.reciprocalTangentCorrection.correctionStateCount> ...
                    candidate.visible.count
                visible = [sourceVisible;obj.xFlex(obj.idx.qxi); ...
                    obj.xFlex(obj.idx.chi)];
            end
            if isfield(obj.reciprocalTangentCorrection, ...
                    'scheduledQueryTrim')
                queryTrim = obj.reciprocalTangentCorrection. ...
                    scheduledQueryTrim;
                actualWing = obj.composeWingControl(uCommand);
                input = [gust(:);actualWing-queryTrim.wing(:); ...
                    obj.trimTailDelta+elevatorIncrement- ...
                        queryTrim.elevator; ...
                    obj.trimThrust+thrustIncrement-queryTrim.thrust];
            else
                input = [gust(:);uCommand(:);elevatorIncrement; ...
                    thrustIncrement];
            end
            if isfield(candidate,'symmetricLongitudinalAudit') && ...
                    logical(candidate.symmetricLongitudinalAudit.enabled)
                symmetryResidual = max(abs([input(2)-input(3); ...
                    input(4)-input(5)]));
                symmetryTolerance = 5e-13+1e-10*max(1,norm(input(2:5),inf));
                assert(symmetryResidual<=symmetryTolerance, ...
                    'PlantRunTime:SymmetricReciprocalInput', ...
                    ['The symmetric-longitudinal reciprocal runtime ', ...
                     'received a differential wing command (%.3e).'], ...
                    symmetryResidual);
            end
            assert(numel(visible)==obj.reciprocalTangentCorrection. ...
                correctionStateCount && ...
                numel(input)==candidate.input.count, ...
                'PlantRunTime:ReciprocalTangentRuntimeOrder', ...
                'The runtime visible-state or input order changed.');
        end

        function [loads,step] = prepareReciprocalTangentCorrection( ...
                obj,loads,visibleBefore,input)
            runtime = obj.reciprocalTangentCorrection;
            candidate = runtime.candidate;
            visibleDeparture = visibleBefore(:)- ...
                runtime.correctionStateEquilibrium;
            inputDeparture = input(:)-candidate.base.inputEquilibrium(:);
            hiddenBefore = runtime.hiddenState(:);

            fixedRootWingWrench = [loads.Fwing_B(:);loads.Mwing_B(:)];
            fixedRootIncrement = zeros(6,1);
            if runtime.replacesFixedRootIncrement
                fixedRootIncrement = fixedRootWingWrench- ...
                    runtime.fixedRootWingEquilibrium;
                % The retained-rigid source owns the incremental moving-base
                % wing/body exchange. Preserve the qualified steady trim
                % wrench, but remove the fixed-root incremental reaction
                % before applying the source-owned reciprocal increment.
                loads.Ftot_B = loads.Ftot_B-fixedRootIncrement(1:3);
                loads.Mtot_B = loads.Mtot_B-fixedRootIncrement(4:6);
            end

            visibleCorrection = ...
                runtime.visibleCorrectionMap*visibleDeparture+ ...
                candidate.partition.Avh*hiddenBefore+ ...
                candidate.correction.deltaBv*inputDeparture;
            hiddenNext = runtime.hiddenStateMap*visibleDeparture+ ...
                candidate.partition.Ahh*hiddenBefore+ ...
                candidate.partition.Bh*inputDeparture;
            operatorInput = [visibleDeparture;hiddenBefore;inputDeparture];
            sourceRootWrench = zeros(6,1);
            if runtime.replacesFixedRootIncrement
                sourceRootWrench = runtime.sourceRootOutput.equilibrium(:)+ ...
                    runtime.sourceRootOutput.runtimeStateMap* ...
                        [visibleDeparture;hiddenBefore]+ ...
                    runtime.sourceRootOutput.runtimeInputMap*inputDeparture;
            end
            finiteCandidate = all(isfinite([visibleCorrection;hiddenNext; ...
                sourceRootWrench]));
            if ~finiteCandidate
                if runtime.replacesFixedRootIncrement
                    loads.Ftot_B = loads.Ftot_B+fixedRootIncrement(1:3);
                    loads.Mtot_B = loads.Mtot_B+fixedRootIncrement(4:6);
                end
                step = struct('visibleCorrection', ...
                    zeros(candidate.visible.count,1), ...
                    'hiddenNext',hiddenBefore, ...
                    'inputDeparture',inputDeparture, ...
                    'wrenchCorrection',zeros(6,1), ...
                    'rigidChartCorrection',zeros(3,1), ...
                    'rigidRealizationResidual',0, ...
                    'rigidWrenchRefinementCount',0, ...
                    'fixedRootIncrement',fixedRootIncrement, ...
                    'bodyAppliedWingWrench',fixedRootWingWrench, ...
                    'sourceRootWrench',zeros(6,1), ...
                    'domainRejected',true, ...
                    'domainReason',"nonfinite reciprocal prediction", ...
                    'accepted',false);
                return
            end
            domain = obj.scheduledReciprocalSourceDomainCheck( ...
                visibleBefore,inputDeparture,hiddenBefore);
            if ~domain.accepted
                if runtime.replacesFixedRootIncrement
                    loads.Ftot_B = loads.Ftot_B+fixedRootIncrement(1:3);
                    loads.Mtot_B = loads.Mtot_B+fixedRootIncrement(4:6);
                end
                step = struct('visibleCorrection', ...
                    zeros(candidate.visible.count,1), ...
                    'hiddenNext',hiddenBefore, ...
                    'inputDeparture',inputDeparture, ...
                    'wrenchCorrection',zeros(6,1), ...
                    'rigidChartCorrection',zeros(3,1), ...
                    'rigidRealizationResidual',0, ...
                    'rigidWrenchRefinementCount',0, ...
                    'fixedRootIncrement',fixedRootIncrement, ...
                    'bodyAppliedWingWrench',fixedRootWingWrench, ...
                    'sourceRootWrench',zeros(6,1), ...
                    'domainRejected',true, ...
                    'domainReason',domain.reason, ...
                    'domainDetails',domain.details, ...
                    'accepted',false);
                return
            end
            wrenchCorrection = runtime.rigidWrenchCommandOperator* ...
                operatorInput;
            desiredRigidCorrection = visibleCorrection( ...
                candidate.visible.rigid);
            traceReciprocal = isfield(obj.cfg,'debug') && ...
                isfield(obj.cfg.debug,'reciprocalTangent') && ...
                logical(obj.cfg.debug.reciprocalTangent);
            traceStride = max(1,round(obj.cfg.ctrl.Ts/obj.dt));
            traceWrenchNorm = norm(wrenchCorrection,inf);
            if traceReciprocal && (mod(obj.k,traceStride)==0 || ...
                    traceWrenchNorm >= 1e2)
                fprintf(['[reciprocal tangent checkpoint] k=%d t=%.6f ', ...
                    '|visible|=%.3e |hidden|=%.3e |input|=%.3e ', ...
                    '|fixedRoot|=%.3e |sourceRoot|=%.3e |w0|=%.3e\n'], ...
                    obj.k,obj.t,norm(visibleDeparture,inf), ...
                    norm(hiddenBefore,inf),norm(inputDeparture,inf), ...
                    norm(fixedRootIncrement,inf),norm(sourceRootWrench,inf), ...
                    traceWrenchNorm);
            end
            baselineNext = obj.stepRigidBody(obj.rb, ...
                loads.Ftot_B,loads.Mtot_B);
            dynamicWrenchMap = runtime.rigidWrenchMap(1:6,:);
            realizationTolerance = 5e-13+1e-8* ...
                max(norm(desiredRigidCorrection(1:6),inf),1e-10);
            refinementCount = 0;
            refinementUpdateLimit = 2;
            if isfield(runtime,'rigidWrenchRefinementUpdateLimit')
                refinementUpdateLimit = double( ...
                    runtime.rigidWrenchRefinementUpdateLimit);
            end
            assert(ismember(refinementUpdateLimit,[2,3]), ...
                'PlantRunTime:RigidWrenchRefinementLimit', ...
                'The rigid-wrench refinement update limit is invalid.');
            for refinement = 0:refinementUpdateLimit
                correctedNext = obj.stepRigidBody(obj.rb, ...
                    loads.Ftot_B+wrenchCorrection(1:3), ...
                    loads.Mtot_B+wrenchCorrection(4:6));
                realizedRigidCorrection = obj.reciprocalRigidState( ...
                    correctedNext)-obj.reciprocalRigidState(baselineNext);
                dynamicResidual = desiredRigidCorrection(1:6)- ...
                    realizedRigidCorrection(1:6);
                dynamicRealizationResidual = norm(dynamicResidual,inf);
                if dynamicRealizationResidual<=realizationTolerance || ...
                        refinement==refinementUpdateLimit
                    break
                end
                wrenchCorrection = wrenchCorrection+ ...
                    dynamicWrenchMap\dynamicResidual;
                refinementCount = refinementCount+1;
            end
            accepted = dynamicRealizationResidual<=realizationTolerance;
            rigidChartCorrection = desiredRigidCorrection(7:9)- ...
                realizedRigidCorrection(7:9);
            if traceReciprocal && (mod(obj.k,traceStride)==0 || ...
                    traceWrenchNorm >= 1e2 || ...
                    norm(wrenchCorrection,inf) >= 1e2 || ~accepted)
                fprintf(['[reciprocal tangent realization] k=%d t=%.6f ', ...
                    'refine=%d accepted=%d |w|=%.3e |residual|=%.3e ', ...
                    'tol=%.3e\n'],obj.k,obj.t,refinementCount,accepted, ...
                    norm(wrenchCorrection,inf),dynamicRealizationResidual, ...
                    realizationTolerance);
            end
            finiteCorrection = all(isfinite([visibleCorrection;hiddenNext; ...
                wrenchCorrection;realizedRigidCorrection]));
            if ~finiteCorrection
                wrenchMapRcond = nan;
                if size(dynamicWrenchMap,1)==size(dynamicWrenchMap,2)
                    wrenchMapRcond = rcond(dynamicWrenchMap);
                end
                nonfiniteTrace = struct( ...
                    'step',obj.k,'timeSeconds',obj.t, ...
                    'visibleFinite',all(isfinite(visibleCorrection)), ...
                    'hiddenFinite',all(isfinite(hiddenNext)), ...
                    'wrenchFinite',all(isfinite(wrenchCorrection)), ...
                    'realizedRigidFinite',all(isfinite(realizedRigidCorrection)), ...
                    'residualFinite',all(isfinite(dynamicResidual)), ...
                    'wrenchMapSize',size(dynamicWrenchMap), ...
                    'wrenchMapRcond',wrenchMapRcond, ...
                    'wrenchNorm',norm(wrenchCorrection,inf), ...
                    'residualNorm',norm(dynamicResidual,inf), ...
                    'realizationTolerance',realizationTolerance);
                obj.last.reciprocalTangentNonfinite = nonfiniteTrace;
                if isfield(obj.cfg,'debug') && isfield(obj.cfg.debug, ...
                        'reciprocalTangent') && logical(obj.cfg.debug.reciprocalTangent)
                    fprintf(['[reciprocal tangent nonfinite] k=%d t=%.6f ', ...
                        'visible=%d hidden=%d wrench=%d realized=%d residual=%d ', ...
                        'rcond=%.3e |wrench|=%.3e |residual|=%.3e tol=%.3e\n'], ...
                        obj.k,obj.t,nonfiniteTrace.visibleFinite, ...
                        nonfiniteTrace.hiddenFinite,nonfiniteTrace.wrenchFinite, ...
                        nonfiniteTrace.realizedRigidFinite, ...
                        nonfiniteTrace.residualFinite, ...
                        nonfiniteTrace.wrenchMapRcond,nonfiniteTrace.wrenchNorm, ...
                        nonfiniteTrace.residualNorm,realizationTolerance);
                end
            end
            assert(finiteCorrection, ...
                'PlantRunTime:ReciprocalTangentFinite', ...
                'The reciprocal tangent produced a nonfinite correction.');

            if accepted
                loads.Freciprocal_B = wrenchCorrection(1:3);
                loads.Mreciprocal_B = wrenchCorrection(4:6);
                loads.Ftot_B = loads.Ftot_B+loads.Freciprocal_B;
                loads.Mtot_B = loads.Mtot_B+loads.Mreciprocal_B;
            else
                if runtime.replacesFixedRootIncrement
                    % Fail closed to the preceding qualified fixed-root path.
                    loads.Ftot_B = loads.Ftot_B+fixedRootIncrement(1:3);
                    loads.Mtot_B = loads.Mtot_B+fixedRootIncrement(4:6);
                end
                visibleCorrection = zeros(candidate.visible.count,1);
                hiddenNext = hiddenBefore;
                wrenchCorrection = zeros(6,1);
                rigidChartCorrection = zeros(3,1);
            end
            if runtime.replacesFixedRootIncrement && accepted
                bodyAppliedWingWrench = ...
                    runtime.fixedRootWingEquilibrium+wrenchCorrection;
            else
                bodyAppliedWingWrench = ...
                    fixedRootWingWrench+wrenchCorrection;
            end
            step = struct('visibleCorrection',visibleCorrection, ...
                'hiddenNext',hiddenNext,'inputDeparture',inputDeparture, ...
                'wrenchCorrection',wrenchCorrection, ...
                'rigidChartCorrection',rigidChartCorrection, ...
                'rigidRealizationResidual',dynamicRealizationResidual, ...
                'rigidWrenchRefinementCount',refinementCount, ...
                'fixedRootIncrement',fixedRootIncrement, ...
                'bodyAppliedWingWrench',bodyAppliedWingWrench, ...
                'sourceRootWrench',sourceRootWrench, ...
                'domainRejected',false,'domainReason',"", ...
                'domainDetails',domain.details, ...
                'accepted',accepted);
        end

        function domain = scheduledReciprocalSourceDomainCheck( ...
                obj,visibleBefore,inputDeparture,hiddenBefore)
        %SCHEDULEDRECIPROCALSOURCEDOMAINCHECK Reject only audit-bound paths.
            domain = struct('accepted',true,'reason',"", ...
                'details',repmat(struct('sourceId',"", ...
                'inputRatio',0,'visibleRatio',0,'hiddenRatio',0, ...
                'pulseHiddenRatio',0,'memoryOwner',"one_step_pulse", ...
                'rootWrenchRatio',0,'accepted',true),0,1));
            runtime = obj.reciprocalTangentCorrection;
            if ~isfield(runtime,'scheduledMembers') || ...
                    isempty(runtime.scheduledMembers) || ...
                    ~isfield(runtime.scheduledMembers(1).candidate, ...
                    'sourceDomain')
                return
            end
            state = obj.currentReciprocalRuntimeState();
            members = runtime.scheduledMembers;
            weights = double(runtime.scheduledWeights(:));
            assert(numel(weights)==numel(members) && ...
                all(isfinite(weights)) && all(weights>=-1e-14) && ...
                abs(sum(weights)-1)<=1e-12, ...
                'PlantRunTime:ScheduledReciprocalSourceDomainWeights', ...
                'The scheduled reciprocal source-domain weights are invalid.');
            activeSource = weights>0;
            assert(any(activeSource), ...
                'PlantRunTime:ScheduledReciprocalSourceDomainEmpty', ...
                'The scheduled reciprocal source-domain stencil is empty.');
            hiddenCounts = arrayfun(@(member) ...
                double(member.candidate.hidden.count),members);
            offsets = [0;cumsum(hiddenCounts(:))];
            assert(numel(hiddenBefore)==offsets(end), ...
                'PlantRunTime:ScheduledReciprocalSourceDomainCount', ...
                'The source-domain hidden-memory partition changed.');
            details = repmat(struct('sourceId',"", ...
                'inputRatio',0,'visibleRatio',0,'hiddenRatio',0, ...
                'pulseHiddenRatio',0,'memoryOwner',"one_step_pulse", ...
                'rootWrenchRatio',0,'accepted',true),numel(members),1);
            for memberIndex = 1:numel(members)
                if ~activeSource(memberIndex)
                    continue
                end
                member = members(memberIndex);
                candidate = member.candidate;
                sourceDomain = candidate.sourceDomain;
                localHidden = offsets(memberIndex)+(1:hiddenCounts(memberIndex));
                departure = member.runtimeFromQueryState*state- ...
                    candidate.base.correctionStateEquilibrium(:);
                assert(isfield(member,'queryRuntimeDeparture') && ...
                    isfield(member,'queryInputDeparture') && ...
                    isfield(member,'queryHiddenEquilibrium') && ...
                    isfield(member,'queryRootWrench') && ...
                    isequal(size(member.queryRuntimeDeparture),size(departure)) && ...
                    isequal(size(member.queryInputDeparture),size(inputDeparture)) && ...
                    isequal(size(member.queryHiddenEquilibrium), ...
                    [hiddenCounts(memberIndex),1]) && ...
                    isequal(size(member.queryRootWrench),[6,1]) && ...
                    all(isfinite([member.queryRuntimeDeparture(:); ...
                    member.queryInputDeparture(:); ...
                    member.queryHiddenEquilibrium(:); ...
                    member.queryRootWrench(:)])), ...
                    'PlantRunTime:ScheduledReciprocalSourceDomainCenter', ...
                    'The scheduled reciprocal source-domain center is invalid.');
                runtimeDeviation = departure-member.queryRuntimeDeparture(:);
                queryInput = inputDeparture+ ...
                    runtime.candidate.base.inputEquilibrium(:);
                actualInput = [queryInput(1); ...
                    runtime.scheduledQueryTrim.wing(:)+queryInput(2:5); ...
                    runtime.scheduledQueryTrim.elevator+queryInput(6); ...
                    runtime.scheduledQueryTrim.thrust+queryInput(7)];
                sourceInputDeparture = [actualInput(1); ...
                    actualInput(2:5)-member.trim.wing(:); ...
                    actualInput(6)-member.trim.elevator; ...
                    actualInput(7)-member.trim.thrust]- ...
                    candidate.base.inputEquilibrium(:);
                inputDeviation = sourceInputDeparture- ...
                    member.queryInputDeparture(:);
                hidden = hiddenBefore(localHidden)- ...
                    member.queryHiddenEquilibrium(:);
                visible = runtimeDeviation(candidate.runtimeState.sourceVisible);
                rootWrench = candidate.rootWrenchOutput.equilibrium(:)+ ...
                    candidate.rootWrenchOutput.runtimeStateMap* ...
                    [departure;hidden+member.queryHiddenEquilibrium(:)]+ ...
                    candidate.rootWrenchOutput.runtimeInputMap* ...
                    sourceInputDeparture;
                rootWrenchDeviation = rootWrench-member.queryRootWrench(:);
                % Source-domain pulses certify the symmetric Case-B chart,
                % not independent left/right commands.  Reject an asymmetric
                % input rather than silently averaging it into the envelope.
                symmetricInput = [inputDeviation(1); ...
                    mean(inputDeviation(2:3)); ...
                    mean(inputDeviation(4:5)); ...
                    inputDeviation(6:7)];
                inputSymmetry = norm([inputDeviation(2)-inputDeviation(3); ...
                    inputDeviation(4)-inputDeviation(5)],inf);
                symmetryTolerance = 100*eps(max(1,norm(inputDeviation,inf)));
                inputRatio = norm(symmetricInput./sourceDomain.inputBounds,inf);
                visibleRatio = norm(visible,inf)/ ...
                    sourceDomain.visibleDepartureInfinity;
                pulseHiddenRatio = norm(hidden,inf)/ ...
                    sourceDomain.hiddenDepartureInfinity;
                hiddenRatio = pulseHiddenRatio;
                memoryOwner = "one_step_pulse";
                if isfield(sourceDomain,'recurrentMemoryCertificate') && ...
                        logical(sourceDomain. ...
                        recurrentMemoryCertificate.enabled)
                    certificate = sourceDomain.recurrentMemoryCertificate;
                    factor = certificate.lyapunovFactorLower;
                    assert(isequal(size(factor), ...
                        [hiddenCounts(memberIndex),hiddenCounts(memberIndex)]) && ...
                        isfinite(certificate.invariantRadius) && ...
                        certificate.invariantRadius>0 && ...
                        certificate.gateLimit==1, ...
                        'PlantRunTime:ScheduledReciprocalMemoryCertificate', ...
                        'The recurrent-memory certificate is invalid for %s.', ...
                        member.sourceId);
                    hiddenRatio = norm(factor'*hidden,2)/ ...
                        certificate.invariantRadius;
                    memoryOwner = "discrete_lyapunov_invariant";
                end
                rootWrenchRatio = norm(rootWrenchDeviation,inf)/ ...
                    sourceDomain.rootWrenchDepartureInfinity;
                tolerance = 100*eps(max(1,max([inputRatio,visibleRatio, ...
                    hiddenRatio,rootWrenchRatio])));
                accepted = inputSymmetry<=symmetryTolerance && ...
                    max([inputRatio,visibleRatio,hiddenRatio, ...
                    rootWrenchRatio])<=1+tolerance;
                details(memberIndex) = struct('sourceId', ...
                    string(member.sourceId),'inputRatio',inputRatio, ...
                    'visibleRatio',visibleRatio,'hiddenRatio',hiddenRatio, ...
                    'pulseHiddenRatio',pulseHiddenRatio, ...
                    'memoryOwner',memoryOwner, ...
                    'rootWrenchRatio',rootWrenchRatio,'accepted',accepted);
                if ~accepted
                    domain.accepted = false;
                    domain.reason = "source_domain_"+string(member.sourceId);
                    break
                end
            end
            domain.details = details;
        end

        function applyReciprocalTangentCorrection(obj,step)
            candidate = obj.reciprocalTangentCorrection.candidate;
            visibleCorrection = step.visibleCorrection;

            obj.xFlex(obj.idx.q1) = obj.xFlex(obj.idx.q1)+ ...
                visibleCorrection(candidate.visible.q1);
            obj.xFlex(obj.idx.q2) = obj.xFlex(obj.idx.q2)+ ...
                visibleCorrection(candidate.visible.q2);
            obj.xFlex(obj.idx.qGam) = obj.xFlex(obj.idx.qGam)+ ...
                visibleCorrection(candidate.visible.qGam);
            obj.rb.euler = obj.rb.euler+step.rigidChartCorrection;
            obj.x(1:obj.nx) = obj.xFlex;

            obj.reciprocalTangentCorrection.hiddenState = step.hiddenNext;
            obj.reciprocalTangentCorrection.lastVisibleCorrection = ...
                visibleCorrection;
            obj.reciprocalTangentCorrection.lastInput = step.inputDeparture;
            obj.reciprocalTangentCorrection.lastWrenchCorrection = ...
                step.wrenchCorrection;
            obj.reciprocalTangentCorrection.lastRigidChartCorrection = ...
                step.rigidChartCorrection;
            obj.reciprocalTangentCorrection. ...
                lastRigidRealizationResidual = ...
                step.rigidRealizationResidual;
            obj.reciprocalTangentCorrection.lastStepAccepted = step.accepted;
            obj.reciprocalTangentCorrection.lastDomainRejected = ...
                logical(step.domainRejected);
            obj.reciprocalTangentCorrection.lastDomainReason = ...
                string(step.domainReason);
            if isfield(step,'domainDetails')
                obj.reciprocalTangentCorrection.lastDomainDetails = ...
                    step.domainDetails;
            else
                obj.reciprocalTangentCorrection.lastDomainDetails = struct([]);
            end
            obj.reciprocalTangentCorrection.lastWrenchRefinementCount = ...
                step.rigidWrenchRefinementCount;
            obj.reciprocalTangentCorrection.lastFixedRootIncrement = ...
                step.fixedRootIncrement;
            obj.reciprocalTangentCorrection.lastBodyAppliedWingWrench = ...
                step.bodyAppliedWingWrench;
            obj.reciprocalTangentCorrection.lastSourceRootWrench = ...
                step.sourceRootWrench;
        end

        function [map,details] = buildReciprocalRigidWrenchMap(obj,method)
            if nargin<2
                method = "central_difference";
            end
            method = string(method);
            assert(isscalar(method) && any(method==[ ...
                "central_difference","analytic_rk4_audit"]), ...
                'PlantRunTime:ReciprocalRigidWrenchMapMethod', ...
                'The reciprocal rigid-wrench map method is invalid.');
            rbCommand = obj.composeRigidCommand(struct());
            loads = obj.computeCoupledLoads(obj.xFlex,obj.lastWingControl, ...
                zeros(obj.nw,1),rbCommand);
            fullSteps = [1e-4*ones(3,1);1e-5*ones(3,1)];
            halfStepConvergence = nan;
            complexStepRelativeError = nan;
            stepParityRelativeError = nan;
            if method=="central_difference"
                map = localCentralMap(fullSteps);
                halfMap = localCentralMap(fullSteps/2);
                halfStepConvergence = norm(map-halfMap,'fro')/ ...
                    max(norm(halfMap,'fro'),eps);
                derivativeCertificationResidual = halfStepConvergence;
            else
                state = [obj.rb.r_I(:);obj.rb.v_B(:); ...
                    obj.rb.euler(:);obj.rb.omega_B(:)];
                wrench = [loads.Ftot_B(:);loads.Mtot_B(:)];
                [analyticNext,~,inputTangent] = ...
                    AeroFlex.sim.rigidRk4Tangent( ...
                        state,wrench,obj.rbParams,obj.dt);
                reciprocalRows = [4:6,10:12,7:9];
                map = inputTangent(reciprocalRows,:);
                complexMap = localComplexMap(wrench);
                complexStepRelativeError = norm(map-complexMap,'fro')/ ...
                    max(1,norm(complexMap,'fro'));
                nonlinearNext = obj.stepRigidBody( ...
                    obj.rb,loads.Ftot_B,loads.Mtot_B);
                nonlinearState = [nonlinearNext.r_I(:); ...
                    nonlinearNext.v_B(:);nonlinearNext.euler(:); ...
                    nonlinearNext.omega_B(:)];
                stepParityRelativeError = ...
                    norm(analyticNext-nonlinearState,inf)/ ...
                    max(1,norm(nonlinearState,inf));
                derivativeCertificationResidual = ...
                    complexStepRelativeError;
            end
            singularValues = svd(map);
            dynamicSingularValues = svd(map(1:6,:));
            details = struct('method',method, ...
                'rank',rank(map,1e-10*singularValues(1)), ...
                'dynamicRank',rank(map(1:6,:), ...
                    1e-10*dynamicSingularValues(1)), ...
                'conditionNumber',dynamicSingularValues(1)/ ...
                    dynamicSingularValues(end), ...
                'halfStepConvergence',halfStepConvergence, ...
                'halfStepConvergenceTolerance',1e-7, ...
                'complexStepRelativeError',complexStepRelativeError, ...
                'complexStepTolerance',1e-10, ...
                'stepParityRelativeError',stepParityRelativeError, ...
                'stepParityTolerance',1e-12, ...
                'derivativeCertificationResidual', ...
                    derivativeCertificationResidual, ...
                'forceStepsNewton',fullSteps(1:3), ...
                'momentStepsNewtonMeter',fullSteps(4:6));

            function result = localCentralMap(steps)
                result = zeros(9,6);
                for column = 1:6
                    plus = zeros(6,1); minus = zeros(6,1);
                    plus(column) = steps(column);
                    minus(column) = -steps(column);
                    plusNext = obj.stepRigidBody(obj.rb, ...
                        loads.Ftot_B+plus(1:3), ...
                        loads.Mtot_B+plus(4:6));
                    minusNext = obj.stepRigidBody(obj.rb, ...
                        loads.Ftot_B+minus(1:3), ...
                        loads.Mtot_B+minus(4:6));
                    result(:,column) = (obj.reciprocalRigidState(plusNext)- ...
                        obj.reciprocalRigidState(minusNext))/(2*steps(column));
                end
            end

            function result = localComplexMap(wrench)
                complexStep = 1e-30;
                result = zeros(9,6);
                for column = 1:6
                    perturbed = complex(wrench);
                    perturbed(column) = ...
                        perturbed(column)+1i*complexStep;
                    next = obj.stepRigidBody(obj.rb, ...
                        perturbed(1:3),perturbed(4:6));
                    result(:,column) = imag( ...
                        obj.reciprocalRigidState(next))/complexStep;
                end
            end
        end

        function restoreRigidSnapshot(obj,rigid)
        %RESTORERIGIDSNAPSHOT Restore a temporary rigid-state evaluation.
            obj.rb = rigid;
        end

        function state = reciprocalRigidState(~,rb)
            state = [rb.v_B(:);rb.omega_B(:);rb.euler(:)];
        end

        function xNext = stepWingROM(obj, uWing, gust)
        %STEPWINGROM Advance the flexible/aeroelastic ROM one plant step.
        %
        % With a scheduled ROM, the flight-condition scaling is already inside the
        % scheduled L/parConst.  Do not pass qRatio in that case.
        
            try
                if obj.scheduleEnabled
                    [xNext, ~] = obj.model.step(obj.xFlex, uWing, gust, [], false);
                else
                    if ~isfield(obj.last,'qRatio') || isempty(obj.last.qRatio)
                        obj.last.qRatio = 1;
                    end
        
                    [xNext, ~] = obj.model.step( ...
                        obj.xFlex, uWing, gust, [], false, obj.last.qRatio);
                end
            catch ME
                error('PlantRunTime:stepWingROM', ...
                    'Could not advance ROMIntegrator.step(). Original error:\n%s', ...
                    ME.message);
            end
        end
        % function xNext = stepWingROM(obj, uWing, gust)
        % %STEPWINGROM Advance flexible/aeroelastic ROM one step.
        % 
        %     if ~isfield(obj.last,'qRatio') || isempty(obj.last.qRatio)
        %         obj.last.qRatio = 1;
        %     end
        % 
        %     try
        %         [xNext, ~] = obj.model.step(obj.xFlex, uWing, gust, [], false, obj.last.qRatio);
        %     catch ME
        %         error('PlantRunTime:stepWingROM', ...
        %               'Could not advance ROMIntegrator.step(). Original error:\n%s', ...
        %               ME.message);
        %     end
        % end
        % function xNext = stepWingROM(obj, u_cmd, gust)
        % %STEPWINGROM Advance flexible/aeroelastic ROM by one plant time step.
        % 
        %     Sdummy = [eye(obj.nx), zeros(obj.nx, obj.nw + obj.nu)];
        % 
        %     try
        %         [xNext, ~] = obj.model.step(obj.xFlex, u_cmd, gust, Sdummy, false, obj.last.qRatio);
        %     catch
        %         % Try robustly because weird errors in some commits.
        %         try
        %             [xNext, ~] = obj.model.step(obj.xFlex, u_cmd, gust, Sdummy, true, obj.last.qRatio);
        %         catch ME
        %             error('PlantRunTime:stepWingROM', ...
        %                   'Could not advance ROMIntegrator.step(). Original error:\n%s', ...
        %                   ME.message);
        %         end
        %     end
        % end
        function [Clamp6, Fwing_B, Mwing_B] = computeWingClampReaction(obj, xFlex, uWing, gust)
        %COMPUTEWINGCLAMPREACTION Recover root resultant from the raw ROM RHS.
        %
        % Propagation uses the projected operator Ldyn.  Clamp load recovery must use
        % raw model.L so that Pr*q1dot contains the reaction needed to hold the root.
        
            if nargin < 2 || isempty(xFlex)
                xFlex = obj.xFlex;
            end
        
            if nargin < 3 || isempty(uWing)
                uWing = obj.lastWingControl;
            end
        
            if nargin < 4 || isempty(gust)
                gust = zeros(obj.nw,1);
            end
        
            pc = obj.model.parConst;
            pc.u_ctrl = uWing(:);
            pc.gust   = gust(:);
        
            % Thrust is handled by the rigid-body load path unless deliberately
            % injected into the flexible model elsewhere.
            pc.N_Thrust = zeros(numel(obj.idx.q1),1);
        
            xdot = AeroFlex.sim.nonlinear_terms(xFlex, pc, obj.idx) + ...
                   obj.model.L*xFlex;
        
            q1dot = xdot(obj.idx.q1);
        
            if isstruct(obj.sched) && isfield(obj.sched,'equilibriumCentered') && ...
                    isfield(obj.sched.equilibriumCentered,'enabled') && ...
                    obj.sched.equilibriumCentered.enabled
                Clamp6 = AeroFlex.sched.recoverCenteredRootWrench( ...
                    obj.sched,xFlex,uWing,obj.centeredRecoveryGust(gust));
            else
                Clamp6 = AeroFlex.beam.recoverRootWrench( ...
                    obj.beam.red.phi1_sA,obj.beam.Pr,obj.beam.Pr*q1dot);
            end
            Clamp6 = Clamp6(:);
        
            if obj.shouldMirrorWingClamp()
                Clamp6 = obj.mirrorWingClamp(Clamp6);
            end
        
            if any(~isfinite(Clamp6))
                error('PlantRunTime:ClampNaN', ...
                    'Non-finite wing clamp load.');
            end
        
            [Fwing_B, Mwing_B] = obj.mapClampToRigidLoads(Clamp6);
        end
        % function [Clamp6, Fwing_B, Mwing_B] = computeWingClampReaction(obj, xFlex, uWing, gust)
        % %COMPUTEWINGCLAMPREACTION Compute wing resultant from the same RHS used by trim.
        % 
        %     if nargin < 2 || isempty(xFlex)
        %         xFlex = obj.xFlex;
        %     end
        % 
        %     if nargin < 3 || isempty(uWing)
        %         uWing = obj.lastWingControl;
        %     end
        % 
        %     if nargin < 4
        %         gust = [];
        %     end
        % 
        %     pc = obj.model.parConst;
        %     pc.u_ctrl = uWing(:);
        %     pc.gust   = gust;
        % 
        %     % Fuselage thrust is handled by rigid loads, not flexible modal forcing.
        %     pc.N_Thrust = zeros(obj.beam.Nm,1);
        % 
        %     xdot = AeroFlex.sim.nonlinear_terms(xFlex, pc, obj.idx) + obj.model.L*xFlex;
        % 
        %     q1dot = xdot(obj.idx.q1);
        % 
        %     Clamp6 = obj.beam.red.phi1_sA * (obj.beam.Pr*q1dot);
        %     Clamp6 = Clamp6(:);
        % 
        %     if obj.shouldMirrorWingClamp()
        %         Clamp6 = obj.mirrorWingClamp(Clamp6);
        %     end
        % 
        %     [Fwing_B, Mwing_B] = obj.mapClampToRigidLoads(Clamp6);
        % 
        %     if any(~isfinite(Clamp6))
        %         error('PlantRunTime:ClampNaN', 'Non-finite wing clamp load.');
        %     end
        % end
        % function Clamp6 = computeWingClampReaction(obj)
        % %COMPUTEWINGCLAMPREACTION Compute root reaction [Fx;Fy;Fz;Mx;My;Mz].
        % %
        % % Mirrors trim force path:
        % %   xdot = nonlinear_terms(x) + L*x
        % %   q1dot = xdot(idx.q1)
        % %   Clamp6 = phi1_sA * (Pr*q1dot)
        % 
        %     % if ~isfield(obj.beam,'Pr')
        %     %     error('PlantRunTime:ClampReaction', ...
        %     %           'beam.Pr is required for clamp reaction extraction.');
        %     % end
        % 
        %     % if ~isfield(obj.beam,'red') || ~isfield(obj.beam.red,'phi1_sA')
        %     %     error('PlantRunTime:ClampReaction', ...
        %     %           'beam.red.phi1_sA is required to map clamp modal reactions to 6D resultant.');
        %     % end
        % 
        %     xdot = AeroFlex.sim.nonlinear_terms(obj.xFlex, obj.model.parConst, obj.idx) ...
        %          + obj.model.L * obj.xFlex;
        % 
        %     q1dot = xdot(obj.idx.q1);
        % 
        %     Clamp6 = obj.beam.red.phi1_sA * (obj.beam.Pr * q1dot);
        %     Clamp6 = obj.mirrorWingClamp(Clamp6);
        % 
        % 
        %     Clamp6 = Clamp6(:);
        % 
        %     if numel(Clamp6) ~= 6
        %         error('PlantRunTime:ClampReaction', ...
        %               'Clamp6 must have length 6. Got length %d.', numel(Clamp6));
        %     end
        % 
        %     % symmetry treatment if modeling one semi-wing but applying
        %     % symmetric pair loads to the aircraft.
        %     if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'symmetricWingPair') ...
        %             && obj.cfg.sim.symmetricWingPair
        % 
        %         % For a symmetric left/right pair:
        %         %   Fx, Fz, My double.
        %         %   Fy, Mx, Mz cancel under exact symmetry.
        %         Clamp6 = [ ...
        %             2*Clamp6(1);
        %             0;
        %             2*Clamp6(3);
        %             0;
        %             2*Clamp6(5);
        %             0 ];
        %     elseif isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'wingMultiplicity')
        %         Clamp6 = obj.cfg.sim.wingMultiplicity * Clamp6;
        %         % Clamp6 = 2 * Clamp6;
        %     end
        % end

        function fc = computeFlightConditionFromRB(obj)
        %COMPUTEFLIGHTCONDITIONFROMRB Compute Uinf, alpha, beta from body velocity.
        
            vB = obj.rb.v_B(:);
        
            U = norm(vB);
        
            if U < 1e-9
                U = obj.safeGetFlightSpeed(obj.cfg);
                alpha = obj.safeGetAlpha(obj.cfg);
                beta = 0;
            else
                % Body axes convention assumed:
                %   x forward, y right, z down.
                %
                % alpha = atan2(w,u)
                % beta  = asin(v/U)
                alpha = atan2(vB(3), vB(1));
                beta  = asin(max(-1,min(1,vB(2)/U)));
            end
        
            fc = struct();
            fc.Uinf    = U;
            fc.alpha   = alpha;
            fc.beta    = beta;
            fc.euler   = obj.rb.euler(:);
            fc.omega_B = obj.rb.omega_B(:);
        end
        
        function updateWingFlightCondition(obj)
        %UPDATEWINGFLIGHTCONDITION Optionally inject rigid-body attitude/AoA into ROM.
        %
        % By default this is a no-op except for storing obj.flightCond.
        %
        % If cfg.sim.coupleRigidAttitudeIntoWingChi = true and idx.chi exists,
        % then the chi portion of the flexible ROM state is overwritten with
        % a rigid attitude/AoA-compatible value.
        
            if ~isfield(obj.cfg,'sim') || ~isfield(obj.cfg.sim,'coupleRigidAttitudeIntoWingChi')
                return
            end
        
            if ~obj.cfg.sim.coupleRigidAttitudeIntoWingChi
                return
            end
        
            if ~isfield(obj.idx,'chi')
                return
            end
        
            aoaRef = obj.safeGetAlpha(obj.cfg);
        
            % Conservative choice consistent with trim initialization:
            %   chi = [0; alpha - alpha_ref; 0]
            %
            % In works:
            %   chi = obj.flightCond.euler;
            chi = [0;
                   obj.flightCond.alpha - aoaRef;
                   0];
        
            obj.xFlex(obj.idx.chi) = chi;
        end

        function chi = applyPackageRelativeRigidAttitudeChi(obj)
        %APPLYPACKAGERELATIVERIGIDATTITUDECHI Bind ROM chi to rigid attitude.
        %   This default-inactive Case-B/C audit contract preserves
        %   thetaAbsolute = thetaEquilibrium(activePackage) + chi(2).
        %   qxi remains the local flexible orientation field. Lateral chi
        %   entries remain zero in the symmetric-longitudinal scope.

            chi = [];
            if ~isfield(obj.cfg,'sim') || ...
                    ~isfield(obj.cfg.sim,'packageRelativeRigidAttitudeChi')
                return
            end
            request = obj.cfg.sim.packageRelativeRigidAttitudeChi;
            assert(isstruct(request) && isfield(request,'enable') && ...
                isscalar(request.enable), ...
                'PlantRunTime:RigidAttitudeChiRequest', ...
                ['cfg.sim.packageRelativeRigidAttitudeChi must be a ', ...
                 'scalar request with an enable field.']);
            if ~logical(request.enable)
                return
            end
            assert(obj.useRigidBody, ...
                'PlantRunTime:RigidAttitudeChiBodyCase', ...
                ['Package-relative rigid-attitude chi is restricted to ', ...
                 'the coupled-full runtime.']);
            assert(isfield(request,'auditOnly') && ...
                isscalar(request.auditOnly) && logical(request.auditOnly), ...
                'PlantRunTime:RigidAttitudeChiAuditOnly', ...
                'The package-relative chi request must remain audit-only.');
            assert(isfield(request,'symmetricLongitudinal') && ...
                isscalar(request.symmetricLongitudinal) && ...
                logical(request.symmetricLongitudinal), ...
                'PlantRunTime:RigidAttitudeChiScope', ...
                ['The approved package-relative chi owner is restricted ', ...
                 'to symmetric-longitudinal operation.']);
            assert(isfield(obj.idx,'chi') && numel(obj.idx.chi) == 3, ...
                'PlantRunTime:RigidAttitudeChiIndex', ...
                'The active ROM must expose exactly three chi states.');
            assert(isfield(obj.last,'sched_mu') && ...
                isnumeric(obj.last.sched_mu) && numel(obj.last.sched_mu) == 2 && ...
                all(isfinite(obj.last.sched_mu(:))), ...
                'PlantRunTime:RigidAttitudeChiPackage', ...
                'The active package must expose finite [U, alphaDeg].');
            assert(isstruct(obj.rb) && isfield(obj.rb,'euler') && ...
                isnumeric(obj.rb.euler) && numel(obj.rb.euler) == 3 && ...
                all(isfinite(obj.rb.euler(:))), ...
                'PlantRunTime:RigidAttitudeChiMeasurement', ...
                'The current measured/fused rigid Euler attitude is invalid.');

            thetaAbsolute = double(obj.rb.euler(2));
            thetaEquilibrium = deg2rad(double(obj.last.sched_mu(2)));
            thetaDeparture = atan2(sin(thetaAbsolute-thetaEquilibrium), ...
                cos(thetaAbsolute-thetaEquilibrium));
            chi = [0;thetaDeparture;0];
            obj.xFlex(obj.idx.chi) = chi;
            if numel(obj.x) >= obj.nx
                obj.x(obj.idx.chi) = chi;
            end
            obj.last.rigidAttitudeChiOwner = struct( ...
                'enabled',true, ...
                'owner','measured_rigid_euler_minus_active_package_equilibrium', ...
                'thetaAbsoluteRad',thetaAbsolute, ...
                'thetaEquilibriumRad',thetaEquilibrium, ...
                'chiRad',chi, ...
                'physicalReconstructionErrorRad', ...
                atan2(sin(thetaEquilibrium+chi(2)-thetaAbsolute), ...
                    cos(thetaEquilibrium+chi(2)-thetaAbsolute)));
        end

        function [Ftail_B, Mtail_B] = computeTailLoads(obj, u_cmd)
        %COMPUTETAILLOADS Tail force/moment in body axes.
        
            U     = obj.flightCond.Uinf;
            alpha = obj.flightCond.alpha;
        
            deltaTail = obj.trimTailDelta;
        
            % Later map a surface command to elevator, to do it here.
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'tailDeltaIndex')
                ii = obj.cfg.sim.tailDeltaIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    deltaTail = u_cmd(ii);
                end
            end
        
            if exist('tailAeroForceMoment','file') == 2
                [Ftail_B, Mtail_B] = tailAeroForceMoment(U, alpha, deltaTail, obj.cfg);
            else
                Ftail_B = zeros(3,1);
                Mtail_B = zeros(3,1);
            end
        
            Ftail_B = Ftail_B(:);
            Mtail_B = Mtail_B(:);
        end

        function [Ffin_B, Mfin_B] = computeFinLoads(obj, u_cmd)
        %COMPUTEFINLOADS Vertical fin force/moment in body axes.
        
            U    = obj.flightCond.Uinf;
            beta = obj.flightCond.beta;
        
            deltaRudder = 0;
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'rudderIndex')
                ii = obj.cfg.sim.rudderIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    deltaRudder = u_cmd(ii);
                end
            end
        
            if exist('finAeroForceMoment','file') == 2
                [Ffin_B, Mfin_B] = finAeroForceMoment(U, beta, deltaRudder, obj.cfg);
            else
                Ffin_B = zeros(3,1);
                Mfin_B = zeros(3,1);
            end
        
            Ffin_B = Ffin_B(:);
            Mfin_B = Mfin_B(:);
        end
    
        function Fgrav_B = computeGravityBody(obj)
        %COMPUTEGRAVITYBODY Gravity force in body axes.
        %
        % Body/NED convention:
        %   inertial gravity vector is [0;0;m*g] in NED.
        %   body axes x forward, y right, z down.
        
            m = obj.rbParams.mass;
            g = obj.rbParams.g;
        
            eul = obj.rb.euler(:);
            C_BI = obj.dcmBodyToInertial(eul);
            C_IB = C_BI.';
        
            Fgrav_I = [0; 0; m*g];
            Fgrav_B = C_IB * Fgrav_I;
        end

        function [Fthrust_B, Mthrust_B] = computeThrustLoads(obj, u_cmd)
        %COMPUTETHRUSTLOADS Thrust force/moment in body axes.
        %
        % Default: constant trim thrust along +x body through CG.
        % TO DO at a later times:
        %   cfg.sim.thrustIndex can map one entry of u_cmd to thrust.
        
            T = obj.trimThrust;
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'thrustIndex')
                ii = obj.cfg.sim.thrustIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    T = u_cmd(ii);
                end
            end
        
            Fthrust_B = [T; 0; 0];
        
            if isfield(obj.rbParams,'rThrust_B')
                rT = obj.rbParams.rThrust_B(:);
            else
                rT = zeros(3,1);
            end
        
            Mthrust_B = cross(rT, Fthrust_B);
        end

        function rbNext = stepRigidBody(obj, rb, F_B, M_B)
        %STEPRIGIDBODY One explicit RK4 step of rigid-body 6-DoF equations.
        
            dt = obj.dt;
        
            y0 = [rb.r_I(:);
                  rb.v_B(:);
                  rb.euler(:);
                  rb.omega_B(:)];
        
            f = @(y) obj.rigidDerivatives(y, F_B, M_B);
        
            k1 = f(y0);
            k2 = f(y0 + 0.5*dt*k1);
            k3 = f(y0 + 0.5*dt*k2);
            k4 = f(y0 + dt*k3);
        
            y1 = y0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        
            rbNext = rb;
            rbNext.r_I     = y1(1:3);
            rbNext.v_B     = y1(4:6);
            rbNext.euler   = y1(7:9);
            rbNext.omega_B = y1(10:12);
        end

        function ydot = rigidDerivatives(obj, y, F_B, M_B)
        %RIGIDDERIVATIVES Newton-Euler rigid-body equations.
        
            r_I     = y(1:3);
            v_B     = y(4:6);
            euler   = y(7:9);
            omega_B = y(10:12);
        
            r_I_unused = r_I;
        
            m = obj.rbParams.mass;
            I = obj.rbParams.I_B;
        
            C_BI = obj.dcmBodyToInertial(euler);
        
            rdot_I = C_BI * v_B;
        
            % Body-frame translational dynamics:
            %   m*(dv_B/dt + omega x v_B) = F_B
            vdot_B = F_B(:)/m - cross(omega_B, v_B);
        
            % Rotational dynamics:
            %   I*omega_dot + omega x (I*omega) = M_B
            omegaDot_B = I \ (M_B(:) - cross(omega_B, I*omega_B));
        
            eulerDot = obj.eulerRates321(euler, omega_B);
        
            ydot = [rdot_I;
                    vdot_B;
                    eulerDot;
                    omegaDot_B];
        end
        function Clamp6_total = mirrorWingClamp(obj,Clamp6_left)
            S = diag([1,-1,1, -1,1,-1]);
            Clamp6_total = Clamp6_left + S*Clamp6_left;
            % Clamp6_total = Clamp6_left;
        end
        function rb = initRigidState(obj, cfg, x0)
        %INITRIGIDSTATE Initialize rigid-body state.
        
            rb = struct();
        
            rb.r_I = zeros(3,1);
        
            U0 = obj.safeGetFlightSpeed(cfg);
            alpha0 = obj.safeGetAlpha(cfg);
        
            % Body-axis velocity. With z down:
            %   u = U*cos(alpha), w = U*sin(alpha)
            % rb.v_B = [U0*cos(alpha0);
            %           0;
            %           U0*sin(alpha0)];
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'vB0')
                rb.v_B = cfg.rigidEOMset.vB0(:);
            else
                rb.v_B = [U0*cos(alpha0);
                          0;
                          U0*sin(alpha0)];
            end

            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'euler0')
                rb.euler = cfg.rigidEOMset.euler0(:);
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'alphaDeg')
                rb.euler = [0; deg2rad(cfg.trim.alphaDeg); 0];
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
                rb.euler = [0; deg2rad(cfg.flight.aoa_deg); 0];
            else
                rb.euler = [0; alpha0; 0];
            end
        
            rb.omega_B = zeros(3,1);
        
            % Optional: parse appended rigid states if x0 is longer than flexible nx.
            % Expected appended format:
            %   [r_I(3); v_B(3); euler(3); omega_B(3)]
            if numel(x0) >= obj.nx + 12
                rbTail = x0(obj.nx+1:obj.nx+12);
                rb.r_I     = rbTail(1:3);
                rb.v_B     = rbTail(4:6);
                rb.euler   = rbTail(7:9);
                rb.omega_B = rbTail(10:12);
            end
        end

        function p = getRigidParams(obj, cfg)
        %GETRIGIDPARAMS Retrieve rigid-body mass/inertia parameters.
        
            p = struct();
        
            % ---------------- mass ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'mass')
                p.mass = cfg.rigidEOMset.mass;
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'mass')
                p.mass = cfg.flight.mass;
            else
                % uses this will implement other later.
                try
                    rParams = RigidBody.methods.paramsRigid_PazyUAV();
                    p.mass = rParams.mass;
                catch
                    error('PlantRunTime:RigidParams', ...
                          ['Rigid-body mass not found. Set cfg.rigidEOMset.mass ', ...
                           'or cfg.flight.mass, or make RigidBody.methods.paramsRigid_PazyUAV() available.']);
                end
            end
        
            % ---------------- inertia ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I_B')
                p.I_B = cfg.rigidEOMset.I_B;
            elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I')
                p.I_B = cfg.rigidEOMset.I;
            else
                try
                    rParams = RigidBody.methods.paramsRigid_PazyUAV();
        
                    if isfield(rParams,'I_B')
                        p.I_B = rParams.I_B;
                    elseif isfield(rParams,'I')
                        p.I_B = rParams.I;
                    elseif all(isfield(rParams,{'Ixx','Iyy','Izz'}))
                        p.I_B = diag([rParams.Ixx, rParams.Iyy, rParams.Izz]);
                    else
                        error('No inertia field found in paramsRigid_PazyUAV.');
                    end
                catch
                    error('PlantRunTime:RigidParams', ...
                          ['Rigid-body inertia not found. Set cfg.rigidEOMset.I_B ', ...
                           'or cfg.rigidEOMset.I.']);
                end
            end
        
            p.I_B = 0.5*(p.I_B + p.I_B.');
        
            if rcond(p.I_B) < 1e-12
                error('PlantRunTime:RigidParams', ...
                      'Rigid-body inertia matrix is singular or ill-conditioned.');
            end
        
            % ---------------- gravity ----------------
            if isfield(cfg,'flight') && isfield(cfg.flight,'g')
                p.g = cfg.flight.g;
            else
                p.g = 9.807;
            end
        
            % ---------------- thrust offset ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'rThrust_B')
                p.rThrust_B = cfg.rigidEOMset.rThrust_B(:);
            else
                p.rThrust_B = zeros(3,1);
            end

            if isfield(cfg,'rigidEOMset') && ...
                    isfield(cfg.rigidEOMset,'massOwnership')
                ownership = cfg.rigidEOMset.massOwnership;
                assert(isstruct(ownership) && isscalar(ownership) && ...
                    isfield(ownership,'owner'), ...
                    'PlantRunTime:RigidMassOwnership', ...
                    'The rigid mass-ownership contract is incomplete.');
                if string(ownership.owner) == ...
                        "nonwing_source_owned_wing_reaction"
                    expected = RigidBody.methods.paramsRigid_PazyUAV();
                    expected = expected.nonwing;
                    assert(abs(p.mass-expected.mass)<=1e-14 && ...
                        norm(p.I_B-expected.J,'fro')<=1e-14 && ...
                        isfield(cfg.rigidEOMset,'rWingRoot_B') && ...
                        norm(cfg.rigidEOMset.rWingRoot_B(:)- ...
                            expected.wingRootFromCM,inf)<=1e-14 && ...
                        isfield(cfg,'tail') && isfield(cfg.tail,'r_B') && ...
                        norm(cfg.tail.r_B(:)- ...
                            expected.tailArmFromCM,inf)<=1e-14, ...
                        'PlantRunTime:MixedRigidMassOwnership', ...
                        ['The V17A non-wing mass, inertia, wing-root arm, ', ...
                         'and tail arm must come from one owner.']);
                    p.rWingRoot_B = cfg.rigidEOMset.rWingRoot_B(:);
                    assert(numel(p.rWingRoot_B)==3 && ...
                        all(isfinite(p.rWingRoot_B)), ...
                        'PlantRunTime:RigidWingRootArm', ...
                        'The non-wing wing-root moment arm must be finite 3x1.');
                end
                p.massOwnership = ownership;
            end
        end

        function u = sanitizeControl(obj, u)
        %SANITIZECONTROL Ensure control vector has length nu.
        
            if isempty(u)
                u = zeros(obj.nu,1);
            end
        
            u = u(:);
        
            if isscalar(u) && obj.nu > 1
                u = repmat(u,obj.nu,1);
            end
        
            if numel(u) ~= obj.nu
                error('PlantRunTime:ControlDimension', ...
                      'Control vector has length %d, expected nu = %d.', ...
                      numel(u), obj.nu);
            end
        end

        function w = sanitizeDisturbance(obj, w)
        %SANITIZEDISTURBANCE Ensure disturbance vector has length nw.
        
            if isempty(w)
                w = zeros(obj.nw,1);
            end
        
            w = w(:);
        
            if isscalar(w) && obj.nw > 1
                w = repmat(w,obj.nw,1);
            end
        
            if numel(w) ~= obj.nw
                error('PlantRunTime:DisturbanceDimension', ...
                      'Disturbance vector has length %d, expected nw = %d.', ...
                      numel(w), obj.nw);
            end
        end

        function xPack = packCoupledState(obj)
        %PACKCOUPLEDSTATE Return current plant state.
        %
        % For wing_only:
        %   xPack = xFlex
        %
        % For fully_coupled:
        %   xPack = [xFlex; r_I; v_B; euler; omega_B]
        
            if obj.useRigidBody
                xPack = [obj.xFlex(:);
                         obj.rb.r_I(:);
                         obj.rb.v_B(:);
                         obj.rb.euler(:);
                         obj.rb.omega_B(:)];
            else
                xPack = obj.xFlex(:);
            end
        end

        function C_BI = dcmBodyToInertial(~, euler)
        %DCMBODYTOINERTIAL 3-2-1 body-to-inertial DCM.
        %
        % euler = [phi; theta; psi] = roll, pitch, yaw.
        %
        % Assumes NED-style inertial axes and body axes x-forward, y-right, z-down.
        
            phi   = euler(1);
            theta = euler(2);
            psi   = euler(3);
        
            cphi = cos(phi);   sphi = sin(phi);
            cth  = cos(theta); sth  = sin(theta);
            cpsi = cos(psi);   spsi = sin(psi);
        
            C_BI = [ ...
                cth*cpsi, ...
                sphi*sth*cpsi - cphi*spsi, ...
                cphi*sth*cpsi + sphi*spsi;
        
                cth*spsi, ...
                sphi*sth*spsi + cphi*cpsi, ...
                cphi*sth*spsi - sphi*cpsi;
        
                -sth, ...
                sphi*cth, ...
                cphi*cth ];
        end


        function eulerDot = eulerRates321(~, euler, omega_B)
        %EULERRATES321 3-2-1 Euler kinematics.
        %
        % omega_B = [p; q; r].
        
            phi   = euler(1);
            theta = euler(2);
        
            p = omega_B(1);
            q = omega_B(2);
            r = omega_B(3);
        
            ctheta = cos(theta);
        
            if abs(ctheta) < 1e-8
                warning('PlantRunTime:EulerSingularity', ...
                        'Pitch angle near Euler singularity.');
                ctheta = sign(ctheta)*1e-8;
            end
        
            T = [ ...
                1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
                0, cos(phi),            -sin(phi);
                0, sin(phi)/ctheta,      cos(phi)/ctheta ];
        
            eulerDot = T * [p; q; r];
        end

        function U = safeGetFlightSpeed(~, cfg)
        %SAFEGETFLIGHTSPEED Retrieve nominal freestream speed.
        
            if isfield(cfg,'flight') && isfield(cfg.flight,'U_inf')
                U = cfg.flight.U_inf;
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'U_inf')
                U = cfg.trim.U_inf;
            else
                U = 0;
            end
        end
        
        function alpha = safeGetAlpha(obj, cfg)
        %SAFEGETALPHA Retrieve trim/reference angle of attack in radians.
        
            if isstruct(obj.trim) && isfield(obj.trim,'alphaDeg') && ...
                    ~isempty(obj.trim.alphaDeg) && isfinite(obj.trim.alphaDeg)
                alpha = deg2rad(obj.trim.alphaDeg);
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'alphaDeg') && ...
                    ~isempty(cfg.trim.alphaDeg) && isfinite(cfg.trim.alphaDeg)
                alpha = deg2rad(cfg.trim.alphaDeg);
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
                alpha = deg2rad(cfg.flight.aoa_deg);
            else
                alpha = 0;
            end
        end
        % function alpha = safeGetAlpha(obj, cfg)
        % %SAFEGETALPHA Retrieve nominal/reference angle of attack in radians.
        % 
        % 
        %     if isprop(obj,'trim') && isfield(obj.trim,'alphaDeg')
        %         alpha = deg2rad(obj.trim.alphaDeg);
        %     elseif isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
        %         alpha = deg2rad(cfg.flight.aoa_deg);
        %     else
        %         alpha = 0;
        %     end
        % end
        % function [Ftail_B, Mtail_B] = computeTailLoadsFromCmd(obj,rbCmd)
        %     U = obj.flightCond.Uinf;
        %     alpha = obj.flightCond.alpha;
        % 
        %     if isfield(rbCmd,'delta_e')
        %         delta_e = rbCmd.delta_e;
        %     else
        %         delta_e = 0;
        %     end
        % 
        %     % if exist('tailAeroForceMoment','file') == 2
        %         % [Ftail_B,Mtail_B] = tailAeroForceMoment(U,alpha,delta_e,obj.cfg);
        %         [Ftail_B,Mtail_B] = RigidBody.localTailAeroForceMoment(U,alpha,delta_e,obj.cfg);
        %     % else
        %         % Ftail_B = zeros(3,1);
        %         % Mtail_B = zeros(3,1);
        %     % end
        % end
        function [Ftail_B, Mtail_B] = computeTailLoadsFromCmd(obj, rbCmd)
            %COMPUTETAILLOADSFROMCMD Horizontal-tail loads in body axes.
            %
            % rbCmd is interpreted as the total rigid command after composeRigidCommand,
            % so rbCmd.delta_e already includes trim elevator plus any increment.
            
            U = obj.flightCond.Uinf;
            alpha = obj.flightCond.alpha;
            
            if isfield(rbCmd,'delta_e') && ~isempty(rbCmd.delta_e)
                delta_e = rbCmd.delta_e;
            else
                delta_e = obj.trimTailDelta;
            end
            
            q_pitch = 0;
            if isfield(obj.flightCond,'omega_B') && numel(obj.flightCond.omega_B) >= 2
                q_pitch = obj.flightCond.omega_B(2);
            elseif isfield(obj.rb,'omega_B') && numel(obj.rb.omega_B) >= 2
                q_pitch = obj.rb.omega_B(2);
            end
            
            [Ftail_B, Mtail_B] = RigidBody.localTailAeroForceMoment( ...
                U, alpha, delta_e, obj.cfg, q_pitch);
            
            Ftail_B = Ftail_B(:);
            Mtail_B = Mtail_B(:);
            
            if any(~isfinite(Ftail_B)) || any(~isfinite(Mtail_B))
                error('PlantRunTime:TailLoadNaN', ...
                    'Tail load calculation returned a non-finite force or moment.');
            end
        end

        function [Ffin_B, Mfin_B] = computeFinLoadsFromCmd(obj,rbCmd)
            U = obj.flightCond.Uinf;
            beta = obj.flightCond.beta;
        
            if isfield(rbCmd,'delta_r')
                delta_r = rbCmd.delta_r;
            else
                delta_r = 0;
            end
        
            if exist('finAeroForceMoment','file') == 2
                [Ffin_B,Mfin_B] = finAeroForceMoment(U,beta,delta_r,obj.cfg);
            else
                Ffin_B = zeros(3,1);
                Mfin_B = zeros(3,1);
            end
        end
        function [Fthrust_B, Mthrust_B] = computeThrustLoadsFromCmd(obj,rbCmd)
            if isfield(rbCmd,'thrust')
                T = rbCmd.thrust;
            else
                T = obj.trimThrust;
            end
        
            Fthrust_B = [T;0;0];
        
            if isfield(obj.rbParams,'rThrust_B')
                rT = obj.rbParams.rThrust_B(:);
            else
                rT = zeros(3,1);
            end
        
            Mthrust_B = cross(rT,Fthrust_B);
        end

        
        
        % function updateScheduledROM(obj, gust)
        %     if ~isprop(obj,'scheduleEnabled') || ~obj.scheduleEnabled || isempty(obj.ROMlib)
        %         return
        %     end
        %     mu = AeroFlex.sched.computeMuFromRB(obj.rb, obj.cfg, gust);
        %     schedNew = AeroFlex.sched.evalLibrary(obj.ROMlib, mu, obj.cfg.library);
        %     obj.sched = schedNew;
        %     obj = AeroFlex.sched.applyToPlant(obj, schedNew);
        % end
        % function updateScheduledROM(obj, gust)
        % %UPDATESCHEDULEDROM Apply scheduled ROM according to cfg.library.updateMode.
        % 
        %     if ~obj.scheduleEnabled || isempty(obj.ROMlib)
        %         return
        %     end
        % 
        %     mode = "perplantstep";
        %     if isfield(obj.cfg,'library') && isfield(obj.cfg.library,'updateMode')
        %         mode = lower(string(obj.cfg.library.updateMode));
        %     end
        % 
        %     switch mode
        %         case {"none","off","fixed","frozentrim"}
        %             return
        % 
        %         case {"perplantstep","dynamic"}
        %             % mu = [obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha)];
        %             mu = AeroFlex.sched.computeMuFromRB(obj.rb, obj.cfg, gust);
        % 
        %         otherwise
        %             error('PlantRunTime:ScheduleMode', ...
        %                   'Unknown cfg.library.updateMode = "%s".', mode);
        %     end
        % 
        %     if any(~isfinite(mu))
        %         error('PlantRunTime:ScheduleMu', ...
        %               'Non-finite scheduled point: [%s].', num2str(mu));
        %     end
        % 
        %     try
        %         schedNew = AeroFlex.sched.evalLibrary(obj.ROMlib, mu, obj.cfg.library);
        %     catch ME
        %         if isfield(obj.cfg.library,'holdLastOnExtrapolate') && obj.cfg.library.holdLastOnExtrapolate
        %             warning('PlantRunTime:ScheduleHold', ...
        %                 'Schedule update failed at mu=[%.4f %.4f]. Holding last schedule.', mu(1), mu(2));
        %             return
        %         else
        %             rethrow(ME);
        %         end
        %     end
        % 
        %     obj.sched = schedNew;
        %     % obj = AeroFlex.sched.applyToPlant(obj, schedNew);
        % 
        %     % 'Pz', obj.beam.Pz,'Pr', obj.beam.Pr,'red', obj.beam.red);
        %     % obj.model = AeroFlex.sched.applyToROMIntegrator(obj.model, obj.sched, obj.cfg);
        %     % obj.model.parConst.RateProject = struct( ...
        %     % 'projSet', true, ...
        %     % 'Pz', obj.beam.Pz,'Pr', obj.beam.Pr,'red', obj.beam.red);
        % 
        %     % 'Pz', obj.beam.Pz);
        %     obj = AeroFlex.sched.applyToPlant(obj, schedNew);
        % 
        %     obj.last.qRatio = 1;
        %     obj.last.sched_mu = mu;
        % 
        %     if isfield(obj.cfg,'debug') && isfield(obj.cfg.debug,'schedule') && obj.cfg.debug.schedule
        %         obj.debugScheduleConsistency('perPlantStep',schedNew);
        %     end
        % end
        function updateScheduledROM(obj, gust, uWing, rbTot)
        %UPDATESCHEDULEDROM Apply scheduled ROM according to cfg.library.updateMode.
        %
        % Stage-2 common coordinates mean the flexible state is not transformed when
        % a schedule is accepted.  However, changing the ROM is a hybrid model switch.
        % A candidate schedule is accepted only if it is locally admissible at the
        % current x,u,w.
        
            traceEnabled = isfield(obj.cfg,'debug') && ...
                isfield(obj.cfg.debug,'schedule') && ...
                logical(obj.cfg.debug.schedule);
            if traceEnabled
                obj.last.scheduleUpdateTrace = struct( ...
                    'reason',"entered",'step',obj.k,'timeSeconds',obj.t, ...
                    'holdStepsRemaining', ...
                    obj.controlSampleScheduleHoldStepsRemaining);
            end

            if obj.controlSampleScheduleHoldStepsRemaining > 0
                if traceEnabled
                    obj.last.scheduleUpdateTrace.reason = "control_sample_hold";
                end
                obj.controlSampleScheduleHoldStepsRemaining = ...
                    obj.controlSampleScheduleHoldStepsRemaining - 1;
                return
            end

            if ~obj.scheduleEnabled || isempty(obj.ROMlib)
                if traceEnabled
                    obj.last.scheduleUpdateTrace.reason = "scheduler_disabled";
                end
                return
            end
        
            if nargin < 2 || isempty(gust)
                gust = zeros(obj.nw,1);
            end
            gust = obj.sanitizeDisturbance(gust);
        
            if nargin < 3 || isempty(uWing)
                uWing = obj.lastWingControl;
            end
            uWing = uWing(:);
        
            if nargin < 4 || isempty(rbTot)
                rbTot = obj.composeRigidCommand(struct());
            end
        
            mode = "perplantstep";
            if isfield(obj.cfg,'library') && isfield(obj.cfg.library,'updateMode')
                mode = lower(string(obj.cfg.library.updateMode));
            end
        
            switch mode
                case {"none","off","fixed","frozentrim"}
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.reason = "fixed_mode";
                    end
                    return
        
                case {"perplantstep","dynamic"}
                    % Use the flight-condition cache computed immediately before this
                    % call. This avoids hidden convention differences in helper routines.
                    muRaw = [obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha)];
                    mu = obj.filteredSchedulePoint(muRaw, rbTot);
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.muRaw = muRaw;
                        obj.last.scheduleUpdateTrace.muCandidate = mu;
                        obj.last.scheduleUpdateTrace.hasScheduleUInfMps = ...
                            isstruct(rbTot) && isfield(rbTot,'scheduleUInfMps');
                    end
                    % mu = obj.filteredSchedulePoint(muRaw);

                otherwise
                    error('PlantRunTime:ScheduleMode', ...
                        'Unknown cfg.library.updateMode = "%s".', mode);
            end
        
            if any(~isfinite(mu))
                error('PlantRunTime:ScheduleMu', ...
                    'Non-finite scheduled point: [%s].', num2str(mu));
            end
        
            % Meaningful deadband.  Tiny changes of 1e-5 deg were causing repeated
            % model switches in the validated no-gust open-loop case.
            if isfield(obj.last,'sched_mu') && ~isempty(obj.last.sched_mu)
                dU = abs(mu(1) - obj.last.sched_mu(1));
                da = abs(mu(2) - obj.last.sched_mu(2));
        
                dUmin = obj.libOpt('updateDeadbandU', 0.05);
                dAmin = obj.libOpt('updateDeadbandAlphaDeg', 0.02);
        
                if dU < dUmin && da < dAmin
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.reason = "deadband_hold";
                    end

                    if ~isfield(obj.last,'schedDeadbandHolds') || isempty(obj.last.schedDeadbandHolds)
                        obj.last.schedDeadbandHolds = 0;
                    end
                
                    obj.last.schedDeadbandHolds = obj.last.schedDeadbandHolds + 1;
                
                    if obj.libOpt('printDeadbandHolds', false)
                        everyHold = obj.libOpt('deadbandPrintEvery', 100);
                
                        if obj.last.schedDeadbandHolds <= 3 || ...
                                mod(obj.last.schedDeadbandHolds, everyHold) == 0
                
                            fprintf(['[sched:deadband] step=%d t=%.6f ', ...
                                     'dU=%.3e < %.3e | da=%.3e < %.3e deg | holds=%d\n'], ...
                                    obj.k, obj.t, dU, dUmin, da, dAmin, ...
                                    obj.last.schedDeadbandHolds);
                        end
                    end
                
                    return
                end
            end
        
            if isfield(obj.last,'schedRejectUntilStep') && obj.k < obj.last.schedRejectUntilStep && ...
                    isfield(obj.last,'schedRejectMu') && ~isempty(obj.last.schedRejectMu)
            
                retryScale = [ ...
                    obj.libOpt('updateDeadbandU', 0.05), ...
                    obj.libOpt('updateDeadbandAlphaDeg', 0.02)];
            
                retryFactor = obj.libOpt('rejectRetryFactor', 1.0);
                dmuReject = abs(mu(:).' - obj.last.schedRejectMu(:).');
            
                if all(dmuReject < retryFactor*retryScale)
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.reason = "reject_backoff";
                    end
                    return
                end
            end

            % Optional minimum spacing between accepted schedule switches.
            minSteps = obj.libOpt('minUpdateIntervalSteps', 0);
            if minSteps > 0 && isfield(obj.last,'sched_step') && ~isempty(obj.last.sched_step)
                if (obj.k - obj.last.sched_step) < minSteps
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.reason = "minimum_spacing";
                    end
                    return
                end
            end
        
            [runtimeLibrary,runtimeLibraryCfg,caseViewTrace] = ...
                obj.resolveRuntimeCaseView(mu);
            if isfield(caseViewTrace,'holdLastAcceptedPackage') && ...
                    logical(caseViewTrace.holdLastAcceptedPackage)
                if ~isfield(obj.last,'schedOutOfCorridorHolds') || ...
                        isempty(obj.last.schedOutOfCorridorHolds)
                    obj.last.schedOutOfCorridorHolds = 0;
                end
                obj.last.schedOutOfCorridorHolds = ...
                    obj.last.schedOutOfCorridorHolds + 1;
                obj.last.schedLastHold = caseViewTrace;
                if traceEnabled
                    obj.last.scheduleUpdateTrace.reason = ...
                        "runtime_case_view_hold";
                    obj.last.scheduleUpdateTrace.caseView = caseViewTrace;
                end
                return
            end
            try
                [schedNew,packageReuse] = obj.resolveScheduledQueryPackage( ...
                    runtimeLibrary,mu,runtimeLibraryCfg);
                if traceEnabled
                    obj.last.scheduleUpdateTrace.exactPackageReuse = packageReuse;
                end
            catch ME
                if traceEnabled
                    obj.last.scheduleUpdateTrace.errorIdentifier = ...
                        string(ME.identifier);
                    obj.last.scheduleUpdateTrace.errorMessage = ...
                        string(ME.message);
                end
                if isfield(obj.cfg,'library') && ...
                        isfield(obj.cfg.library,'holdLastOnExtrapolate') && ...
                        obj.cfg.library.holdLastOnExtrapolate
        
                    warning('PlantRunTime:ScheduleHold', ...
                        'Schedule update failed at mu=[%.4f %.4f]. Holding last schedule.', ...
                        mu(1), mu(2));
                    if traceEnabled
                        obj.last.scheduleUpdateTrace.reason = "library_hold";
                    end
                    return
                else
                    rethrow(ME);
                end
            end
        
            oldSched = obj.currentScheduleSnapshot();
            if traceEnabled
                obj.last.scheduleUpdateTrace.caseView = caseViewTrace;
            end

            % forces_0 is an affine trim/preload term.  It should not be switched
            % aggressively with gust-induced alpha.  Apply the runtime force-offset
            % policy before admissibility checks.
            schedRaw = schedNew;
            schedNew = obj.applyScheduleForceOffsetPolicy(oldSched, schedNew);
            transition = struct('enabled',false);
            xCandidate = obj.xFlex;
            oldDiag = obj.evaluateScheduleAtState( ...
                oldSched,obj.xFlex,uWing,gust,rbTot);
            if obj.atomicStateTransportEnabled()
                transition = obj.buildAtomicScheduleTransition( ...
                    oldSched,schedNew,schedRaw,mu,uWing,gust,rbTot,oldDiag);
                xCandidate = transition.newState;
            end
            newDiag = obj.evaluateScheduleAtState( ...
                schedNew,xCandidate,uWing,gust,rbTot);
            
            % Optional raw-candidate diagnostic.  This keeps visibility into what the
            % unmodified scheduled forces_0 would have done.
            rawDiag = [];
            if obj.libOpt('forceOffsetRawDiagnostic', false)
                rawDiag = obj.evaluateScheduleAtState(schedRaw, obj.xFlex, uWing, gust, rbTot);
            end
            
            blockDiag = [];
            if obj.libOpt('blockJumpDiagnostic', false)
                blockDiag = obj.evaluateScheduleBlockJumps( ...
                    oldSched, schedNew, obj.xFlex, uWing, gust, rbTot, oldDiag);
            end

            % Diagnostic only: many scheduled ROMs store u_eq. If the force maps are
            % local increments about u_eq, then evaluating a candidate with absolute
            % uWing can create an artificial moment jump. This branch tests that
            % hypothesis without changing the committed dynamics.
            relUDiag = [];
            if obj.libOpt('inputOffsetDiagnostic', true) && isfield(schedNew,'u_eq') && ...
                    ~isempty(schedNew.u_eq)
            
                uWingRelNew = obj.scheduleInputIncrement(schedNew, uWing);
                relUDiag = obj.evaluateScheduleAtState(schedNew, obj.xFlex, uWingRelNew, gust, rbTot);
                relUDiag.uUsed = uWingRelNew;
            end
            
            affDiag = [];
            if obj.libOpt('affineStateDiagnostic', true) && ...
                    isfield(oldSched,'x_eq') && isfield(schedNew,'x_eq') && ...
                    ~isempty(oldSched.x_eq) && ~isempty(schedNew.x_eq) && ...
                    numel(oldSched.x_eq) == numel(obj.xFlex) && ...
                    numel(schedNew.x_eq) == numel(obj.xFlex)
        
                xAff = obj.xFlex(:) + oldSched.x_eq(:) - schedNew.x_eq(:);
                affDiag = obj.evaluateScheduleAtState(schedNew, xAff, uWing, gust, rbTot);
            end
            if obj.libOpt('affineStateDiagnostic', true) && isempty(affDiag)
                if ~isfield(obj.last,'affineDiagSkipped') || isempty(obj.last.affineDiagSkipped)
                    obj.last.affineDiagSkipped = 0;
                end
            
                obj.last.affineDiagSkipped = obj.last.affineDiagSkipped + 1;
            
                if obj.last.affineDiagSkipped <= 3
                    fprintf(['[sched:affine diagnostic skipped] ', ...
                             'oldHasXeq=%d newHasXeq=%d | nx=%d\n'], ...
                            isfield(oldSched,'x_eq'), isfield(schedNew,'x_eq'), ...
                            numel(obj.xFlex));
                end
            end
        
            % gate = obj.assessScheduleCandidate(oldSched, schedNew, oldDiag, newDiag, affDiag);
            gate = obj.assessScheduleCandidate(oldSched, schedNew, oldDiag, ...
                newDiag, affDiag, relUDiag, transition);
            if ~isempty(rawDiag)
                gate.rawForces0 = struct();
                gate.rawForces0.rPz = rawDiag.rPz;
                gate.rawForces0.dF_B = rawDiag.Ftot_B - oldDiag.Ftot_B;
                gate.rawForces0.dM_B = rawDiag.Mtot_B - oldDiag.Mtot_B;
            end
            
            gate.blockDiag = blockDiag;

            % The default scheduler remains a single fail-closed gate.  An
            % explicitly enabled audit may recover a coordinate-continuity
            % rejection by testing progressively smaller moves toward the
            % same raw target.
            % Every candidate is subjected to the unchanged interpolation,
            % state-transport, residual, and centered-load gates before the
            % active package or state is modified.
            backtrackingEligible = ...
                obj.scheduleGateBacktrackingReasonEligible(gate.reason);
            if ~gate.accept && backtrackingEligible && ...
                    obj.scheduleGateBacktrackingEnabled()
                candidate = obj.backtrackScheduleGateCandidate( ...
                    mu,oldSched,oldDiag,uWing,gust,rbTot,gate);
                if candidate.accepted
                    mu = candidate.mu;
                    schedNew = candidate.schedNew;
                    transition = candidate.transition;
                    newDiag = candidate.newDiag;
                    blockDiag = candidate.blockDiag;
                    affDiag = candidate.affDiag;
                    gate = candidate.gate;
                    gate.rawForces0 = candidate.rawForces0;
                    gate.blockDiag = blockDiag;
                else
                    gate.scheduleBacktracking = candidate.info;
                end
            end

            % temp
            % gate.accept = true;
            if ~gate.accept
                if ~isfield(obj.last,'schedRejectedCount') || isempty(obj.last.schedRejectedCount)
                    obj.last.schedRejectedCount = 0;
                end
                obj.last.schedRejectedCount = obj.last.schedRejectedCount + 1;
                obj.last.schedLastReject = gate;
        
                if obj.shouldPrintScheduleGate(false)
                    obj.printScheduleGate('rejected', mu, oldSched, schedNew, oldDiag, newDiag, affDiag, gate);
                end
                
                rejectBackoffSteps = obj.libOpt('rejectBackoffSteps', 0);
                if rejectBackoffSteps > 0
                    obj.last.schedRejectUntilStep = obj.k + rejectBackoffSteps;
                    obj.last.schedRejectMu = mu;
                end

                return
            end
        
            % Candidate is admissible. Commit it.
            if transition.enabled
                obj.xFlex = transition.newState;
                obj.x(1:obj.nx) = transition.newState;
                obj.last.scheduleStateTransport = transition;
            end
            obj.sched = schedNew;
            obj = AeroFlex.sched.applyToPlant(obj, schedNew);
        
            obj.last.qRatio        = 1;
            obj.last.sched_mu      = mu;
            obj.last.sched_step    = obj.k;
            obj.last.sched_ids     = schedNew.pointIds;
            obj.last.sched_weights = schedNew.weights;
            obj.last.schedGate     = gate;
        
            if obj.shouldPrintScheduleGate(true)
                obj.printScheduleGate('accepted', mu, oldSched, schedNew, oldDiag, newDiag, affDiag, gate);
            end
        end

        function [queryLibrary,queryCfg,trace] = resolveRuntimeCaseView(obj,mu)
        %RESOLVERUNTIMECASEVIEW Bind an approved profile-owned runtime stencil.
        % The default path deliberately returns the complete library unchanged.

            queryLibrary = obj.ROMlib;
            queryCfg = obj.cfg.library;
            trace = struct('enabled',false);
            if ~isfield(queryCfg,'runtimeCaseViewOwner') || ...
                    ~isstruct(queryCfg.runtimeCaseViewOwner) || ...
                    ~isfield(queryCfg.runtimeCaseViewOwner,'enabled') || ...
                    ~logical(queryCfg.runtimeCaseViewOwner.enabled)
                return
            end

            owner = queryCfg.runtimeCaseViewOwner;
            assert(isfield(owner,'auditOnly') && logical(owner.auditOnly), ...
                'PlantRunTime:RuntimeCaseViewAuditOnly', ...
                'The runtime case-view owner must remain explicitly audit-only.');
            required = {'speedRangeMps','alphaOffsetDeg','tolerance', ...
                'sourceSetKeys','caseManifestSha256'};
            missing = required(~isfield(owner,required));
            assert(isempty(missing),'PlantRunTime:RuntimeCaseViewConfig', ...
                'The runtime case-view owner is missing: %s.',strjoin(missing,', '));
            speedRange = double(owner.speedRangeMps(:).');
            alphaOffset = double(owner.alphaOffsetDeg);
            tolerance = double(owner.tolerance);
            assert(isequal(size(speedRange),[1,2]) && ...
                all(isfinite(speedRange)) && speedRange(1) <= speedRange(2) && ...
                isscalar(alphaOffset) && isfinite(alphaOffset) && ...
                isscalar(tolerance) && isfinite(tolerance) && tolerance >= 0, ...
                'PlantRunTime:RuntimeCaseViewConfig', ...
                'The runtime case-view centerline bounds are invalid.');
            inCorridor = mu(1) >= speedRange(1)-tolerance && ...
                mu(1) <= speedRange(2)+tolerance && ...
                abs(mu(2)-(alphaOffset-mu(1))) <= tolerance;
            if ~inCorridor
                holdLast = isfield(queryCfg,'holdLastOnExtrapolate') && ...
                    logical(queryCfg.holdLastOnExtrapolate);
                if holdLast
                    trace = struct('enabled',true, ...
                        'query',double(mu(:).'), ...
                        'holdLastAcceptedPackage',true, ...
                        'reason',"out_of_corridor", ...
                        'speedRangeMps',speedRange, ...
                        'alphaOffsetDeg',alphaOffset, ...
                        'tolerance',tolerance);
                    return
                end
                error('PlantRunTime:RuntimeCaseViewOutOfCorridor', ...
                    ['The Case-B runtime query [%.12g %.12g] is outside ', ...
                     'the approved centerline corridor.'],mu(1),mu(2));
            end

            if isfield(queryCfg,'caseView')
                queryCfg = rmfield(queryCfg,'caseView');
            end
            library = AeroFlex.sched.loadLibrary(obj.ROMlib);
            pointCount = numel(library.points);
            if pointCount == 1
                pointMu = double(library.points(1).mu(:).');
                assert(norm(mu(:).'-pointMu,inf) <= tolerance, ...
                    'PlantRunTime:RuntimeCaseViewExactMiss', ...
                    'The single-source Case-B runtime view missed its exact node.');
                ids = 1;
                weightInfo = struct('mode',"exact");
            elseif pointCount == 2
                pointMu = reshape([library.points.mu],numel(mu),[]).';
                direction = pointMu(2,:)-pointMu(1,:);
                fraction = dot(mu(:).'-pointMu(1,:),direction) / ...
                    dot(direction,direction);
                projection = pointMu(1,:)+fraction*direction;
                assert(all(isfinite(direction)) && norm(direction,2) > 0 && ...
                    norm(projection-mu(:).',inf) <= tolerance && ...
                    fraction >= -tolerance && fraction <= 1+tolerance, ...
                    'PlantRunTime:RuntimeCaseViewSegment', ...
                    'The two-source Case-B runtime view does not contain its query.');
                ids = [1,2];
                weightInfo = struct('mode',"case_linear2d");
            else
                [~,ids,weightInfo] = AeroFlex.sched.interpWeights( ...
                    library,mu,queryCfg);
            end
            ids = ids(:).';
            queryLibrary = library;
            queryLibrary.points = library.points(ids);
            % `points`, `mu`, and `sourcePaths` are parallel library fields.
            % Keep the governed view self-consistent before `evalLibrary`
            % validates its case-view geometry.
            queryLibrary.mu = reshape([queryLibrary.points.mu], ...
                numel(mu),[]).';
            if isfield(library,'sourcePaths') && ...
                    numel(library.sourcePaths) == pointCount
                queryLibrary.sourcePaths = library.sourcePaths(ids);
            end
            sourceIds = string({queryLibrary.points.name});
            sourceKey = strjoin(sort(sourceIds),'|');
            approvedKeys = string(owner.sourceSetKeys(:));
            assert(nnz(approvedKeys == sourceKey) == 1, ...
                'PlantRunTime:RuntimeCaseViewStencil', ...
                ['The Case-B runtime query [%.12g %.12g] resolved an ', ...
                 'unapproved source set {%s}.'],mu(1),mu(2),sourceKey);
            if numel(sourceIds) == 1
                mode = "exact";
            elseif numel(sourceIds) == 2
                mode = "case_linear2d";
            elseif numel(sourceIds) == 3
                mode = "case_barycentric2d";
            else
                error('PlantRunTime:RuntimeCaseViewStencil', ...
                    'The Case-B runtime source set has unsupported cardinality %d.', ...
                    numel(sourceIds));
            end
            queryCfg.caseView = struct( ...
                'testId',"caseb_runtime_" + replace(sourceKey,"|","_"), ...
                'query',double(mu(:).'),'heldOutSourceIds',strings(0,1), ...
                'permittedSourceIds',sourceIds, ...
                'expectedInterpolationMode',mode, ...
                'caseManifestSha256',string(owner.caseManifestSha256));
            trace = struct('enabled',true,'query',double(mu(:).'), ...
                'sourceIds',sourceIds,'sourceSetKey',sourceKey, ...
                'interpolationMode',mode,'rawMode',string(weightInfo.mode));
        end

        function S = currentScheduleSnapshot(obj)
        %CURRENTSCHEDULESNAPSHOT Build a lightweight snapshot of the active ROM.
        %
        % This is used only for scheduling diagnostics/admissibility checks.  The
        % active model data are always taken from the plant/model objects.  Optional
        % scheduling metadata such as x_eq, u_eq, pointIds, and weights are recovered
        % from obj.sched or obj.last when available.
        
            S = struct();
        
            % Active committed dynamics and load-recovery data.
            S.L        = obj.model.L;
            S.idx      = obj.idx;
            S.parConst = obj.model.parConst;
            S.beam     = obj.beam;
        
            % Scheduling point.
            if isfield(obj.last,'sched_mu') && ~isempty(obj.last.sched_mu)
                S.mu = obj.last.sched_mu;
            else
                S.mu = [obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha)];
            end
        
            % obj is a class object, not a struct.  Use isprop here, not isfield.
            if isprop(obj,'sched') && isstruct(obj.sched) && ~isempty(obj.sched)
        
                if isfield(obj.sched,'pointIds')
                    S.pointIds = obj.sched.pointIds;
                end
        
                if isfield(obj.sched,'weights')
                    S.weights = obj.sched.weights;
                end
        
                if isfield(obj.sched,'info')
                    S.info = obj.sched.info;
                end
        
                if isfield(obj.sched,'x_eq') && ~isempty(obj.sched.x_eq)
                    S.x_eq = obj.sched.x_eq(:);
                end
        
                if isfield(obj.sched,'u_eq') && ~isempty(obj.sched.u_eq)
                    S.u_eq = obj.sched.u_eq(:);
                end

                if isfield(obj.sched,'scheduleStateCoordinate')
                    S.scheduleStateCoordinate = ...
                        obj.sched.scheduleStateCoordinate;
                end

                % Preserve the active package-owned centered-recovery
                % contract.  Omitting these fields makes the old side of a
                % scheduler transition fall back to legacy raw recovery
                % while the candidate side uses centered recovery, creating
                % a false load jump at the exact-node/interior boundary.
                centeredFields = {'equilibriumCentered', ...
                    'physicalRecoveryReferenceState', ...
                    'physicalRecoveryReferenceControl'};
                for centeredFieldIndex = 1:numel(centeredFields)
                    centeredField = centeredFields{centeredFieldIndex};
                    if isfield(obj.sched,centeredField)
                        S.(centeredField) = obj.sched.(centeredField);
                    end
                end
            end
        
            % Fallbacks from plant.last.  These are useful if the active sched struct
            % was not preserved by an older applyToPlant implementation.
            if ~isfield(S,'pointIds') && isfield(obj.last,'sched_pointIds')
                S.pointIds = obj.last.sched_pointIds;
            end
        
            if ~isfield(S,'weights') && isfield(obj.last,'sched_weights')
                S.weights = obj.last.sched_weights;
            end
        
            if ~isfield(S,'info') && isfield(obj.last,'sched_info')
                S.info = obj.last.sched_info;
            end
        
            if ~isfield(S,'x_eq') && isfield(obj.last,'sched_x_eq') && ...
                    ~isempty(obj.last.sched_x_eq)
                S.x_eq = obj.last.sched_x_eq(:);
            end
        
            if ~isfield(S,'u_eq') && isfield(obj.last,'sched_u_eq') && ...
                    ~isempty(obj.last.sched_u_eq)
                S.u_eq = obj.last.sched_u_eq(:);
            end
        end

        function tf = hasMember(~, S, name)
        %HASMEMBER True for either struct fields or object properties.
        
            tf = false;
        
            if isstruct(S)
                tf = isfield(S, name);
            elseif isobject(S)
                tf = isprop(S, name);
            end
        end
        
        function val = getMember(obj, S, name, defaultVal)
        %GETMEMBER Read a struct field or object property.
        
            if nargin < 4
                defaultVal = [];
            end
        
            if obj.hasMember(S, name)
                val = S.(name);
            else
                val = defaultVal;
            end
        end
        
        function val = getNestedMember(obj, S, path, defaultVal)
        %GETNESTEDMEMBER Safe nested access for mixed structs/objects.
        
            if nargin < 4
                defaultVal = [];
            end
        
            val = S;
        
            for k = 1:numel(path)
                if isempty(val) || ~obj.hasMember(val, path{k})
                    val = defaultVal;
                    return
                end
        
                val = obj.getMember(val, path{k}, defaultVal);
            end
        end
        function D = evaluateScheduleAtState(obj, schedUse, xEval, uWing, gust, rbTot)
        %EVALUATESCHEDULEATSTATE Evaluate RHS and loads for a candidate schedule.
        %
        % This routine is read-only. It does not modify obj.model, obj.beam, obj.idx,
        % or obj.xFlex. It is the core admissibility test for runtime scheduling.
        
            if nargin < 3 || isempty(xEval)
                xEval = obj.xFlex;
            end
            if nargin < 4 || isempty(uWing)
                uWing = obj.lastWingControl;
            end
            if nargin < 5 || isempty(gust)
                gust = zeros(obj.nw,1);
            end
            if nargin < 6 || isempty(rbTot)
                rbTot = obj.composeRigidCommand(struct());
            end
        
            xEval = xEval(:);
            uWing = uWing(:);
            gust  = gust(:);
        
            idxUse  = schedUse.idx;
            LUse    = schedUse.L;
            pcUse   = schedUse.parConst;
            beamUse = schedUse.beam;
        
            pcUse.u_ctrl = uWing;
            pcUse.gust   = gust;
        
            % Thrust is a rigid-body force for the coupled Pazy setup unless
            % explicitly modeled as a modal wing force.
            pcUse.N_Thrust = zeros(numel(idxUse.q1),1);
        
            raw = AeroFlex.sim.nonlinear_terms(xEval, pcUse, idxUse) + LUse*xEval;
            q1raw = raw(idxUse.q1);
        
            % Projection/reaction maps may live either in parConst.RateProject or in the
            % scheduled beam object.  The beam may be a struct or a BeamModel object, so
            % use mixed struct/object-safe access.
            Pz = obj.getNestedMember(pcUse, {'RateProject','Pz'}, []);
            
            if isempty(Pz)
                Pz = obj.getNestedMember(beamUse, {'Pz'}, []);
            end
            
            Pr = obj.getNestedMember(beamUse, {'Pr'}, []);
            phi1_sA = obj.getNestedMember(beamUse, {'red','phi1_sA'}, []);
            
            if isempty(Pz)
                error('PlantRunTime:ScheduleDiagnostic', ...
                    ['Could not find Pz in candidate schedule. ', ...
                     'This usually means a diagnostic mixed schedule replaced parConst ', ...
                     'with an interpolated bundle whose RateProject.Pz is empty, while ', ...
                     'beam.Pz could not be read.']);
            end
            
            if isempty(Pr)
                error('PlantRunTime:ScheduleDiagnostic', ...
                    'Could not find Pr in candidate schedule.');
            end
            
            if isempty(phi1_sA)
                error('PlantRunTime:ScheduleDiagnostic', ...
                    'Could not find beam.red.phi1_sA in candidate schedule.');
            end
        
            if isfield(schedUse,'equilibriumCentered') && ...
                    isfield(schedUse.equilibriumCentered,'enabled') && ...
                    logical(schedUse.equilibriumCentered.enabled)
                % Use the same package-owned absolute reaction contract as
                % computeWingClampReaction.  The scheduler gate must certify
                % the load path that will actually be applied after commit.
                Clamp6 = AeroFlex.sched.recoverCenteredRootWrench( ...
                    schedUse,xEval,uWing,obj.centeredRecoveryGust(gust));
            else
                Clamp6 = phi1_sA*(Pr*q1raw);
            end
            Clamp6 = Clamp6(:);
        
            if obj.shouldMirrorWingClamp()
                Clamp6 = obj.mirrorWingClamp(Clamp6);
            end
        
            [Fwing_B, Mwing_B] = obj.mapClampToRigidLoads(Clamp6);
        
            [Ftail_B, Mtail_B] = obj.computeTailLoadsFromCmd(rbTot);
            [Ffin_B,  Mfin_B ] = obj.computeFinLoadsFromCmd(rbTot);
            Fgrav_B = obj.computeGravityBody();
            [Fthrust_B, Mthrust_B] = obj.computeThrustLoadsFromCmd(rbTot);
        
            Ftot_B = Fwing_B + Ftail_B(:) + Ffin_B(:) + Fgrav_B(:) + Fthrust_B(:);
            Mtot_B = Mwing_B + Mtail_B(:) + Mfin_B(:) + Mthrust_B(:);
        
            D = struct();
            D.raw       = raw;
            D.q1raw     = q1raw;
            D.rRaw      = norm(raw);
            D.rQ1       = norm(q1raw);
            D.rPz       = norm(Pz*q1raw);
            D.rPr       = norm(Pr*q1raw);
            D.Clamp6    = Clamp6;
            D.Fwing_B   = Fwing_B(:);
            D.Mwing_B   = Mwing_B(:);
            D.Ftot_B    = Ftot_B(:);
            D.Mtot_B    = Mtot_B(:);
            D.Ftail_B   = Ftail_B(:);
            D.Mtail_B   = Mtail_B(:);
            D.Fgrav_B   = Fgrav_B(:);
            D.Fthrust_B = Fthrust_B(:);
        end

        function enabled = atomicStateTransportEnabled(obj)
        %ATOMICSTATETRANSPORTENABLED Read the default-disabled Case-B flag.
            request = obj.getNestedMember( ...
                obj.cfg,{'library','atomicStateTransport'},struct());
            enabled = isstruct(request) && isscalar(request) && ...
                isfield(request,'enabled') && ...
                isscalar(request.enabled) && logical(request.enabled);
        end

        function enabled = scheduleGateBacktrackingEnabled(obj)
        %SCHEDULEGATEBACKTRACKINGENABLED Read the audit-only scheduler flag.
            request = obj.getNestedMember( ...
                obj.cfg,{'library','scheduleGateBacktracking'},struct());
            enabled = isstruct(request) && isscalar(request) && ...
                isfield(request,'enabled') && isscalar(request.enabled) && ...
                logical(request.enabled);
            if enabled
                assert(isfield(request,'auditOnly') && ...
                    isscalar(request.auditOnly) && logical(request.auditOnly), ...
                    'PlantRunTime:ScheduleBacktrackingAuditOnly', ...
                    ['Schedule-gate backtracking is default-inactive and ', ...
                    'requires an explicit audit-only request.']);
            end
        end

        function eligible = scheduleGateBacktrackingReasonEligible(~,reason)
        %SCHEDULEGATEBACKTRACKINGREASONELIGIBLE Bound recoverable gate reasons.
            eligibleReasons = ["state_transport_invariant", ...
                "transported_vector_field_curvature","load_jump"];
            eligible = any(string(reason)==eligibleReasons);
        end

        function candidate = backtrackScheduleGateCandidate( ...
                obj,rawMu,oldSched,oldDiag,uWing,gust,rbTot,rawGate)
        %BACKTRACKSCHEDULEGATECANDIDATE Find the largest fully admissible move.
        %
        % The active package and state are not modified here.  A candidate is
        % returned for atomic commit only after every existing gate accepts it.
            request = obj.getNestedMember( ...
                obj.cfg,{'library','scheduleGateBacktracking'},struct());
            contraction = double(obj.getNestedMember( ...
                request,{'contraction'},0.5));
            maximumTrials = double(obj.getNestedMember( ...
                request,{'maximumTrials'},12));
            assert(isscalar(contraction) && isfinite(contraction) && ...
                contraction > 0 && contraction < 1, ...
                'PlantRunTime:ScheduleBacktrackingContraction', ...
                'The backtracking contraction must lie strictly between zero and one.');
            assert(isscalar(maximumTrials) && isfinite(maximumTrials) && ...
                maximumTrials >= 2 && maximumTrials == floor(maximumTrials), ...
                'PlantRunTime:ScheduleBacktrackingTrials', ...
                'The backtracking trial limit must be an integer of at least two.');
            assert(isfield(obj.last,'sched_mu') && ...
                all(isfinite(obj.last.sched_mu(:))), ...
                'PlantRunTime:ScheduleBacktrackingOrigin', ...
                'Backtracking requires a finite preceding accepted coordinate.');

            originMu = double(obj.last.sched_mu(:).');
            rawMu = double(rawMu(:).');
            candidate = struct('accepted',false,'mu',originMu, ...
                'schedNew',[],'schedRaw',[], ...
                'transition',struct('enabled',false), ...
                'xCandidate',obj.xFlex,'newDiag',[], ...
                'rawDiag',[],'rawForces0',[], ...
                'blockDiag',[],'affDiag',[],'relUDiag',[], ...
                'gate',rawGate,'info',struct());
            trialScale = nan(maximumTrials,1);
            trialMu = nan(maximumTrials,2);
            trialReason = strings(maximumTrials,1);
            trialMaximumLoadRatio = nan(maximumTrials,1);
            trialAccepted = false(maximumTrials,1);
            trialScale(1) = 1;
            trialMu(1,:) = rawMu;
            trialReason(1) = string(rawGate.reason);
            trialMaximumLoadRatio(1) = obj.scheduleGateLoadRatio(rawGate);

            completedTrials = 1;
            for trialIndex = 2:maximumTrials
                scale = contraction^(trialIndex-1);
                muTry = originMu + scale*(rawMu-originMu);
                completedTrials = trialIndex;
                trialScale(trialIndex) = scale;
                trialMu(trialIndex,:) = muTry;
                try
                    schedRawTry = AeroFlex.sched.evalLibrary( ...
                        obj.ROMlib,muTry,obj.cfg.library);
                    schedTry = obj.applyScheduleForceOffsetPolicy( ...
                        oldSched,schedRawTry);
                    transitionTry = struct('enabled',false);
                    xTry = obj.xFlex;
                    if obj.atomicStateTransportEnabled()
                        transitionTry = obj.buildAtomicScheduleTransition( ...
                            oldSched,schedTry,schedRawTry,muTry,uWing,gust, ...
                            rbTot,oldDiag);
                        xTry = transitionTry.newState;
                    end
                    newDiagTry = obj.evaluateScheduleAtState( ...
                        schedTry,xTry,uWing,gust,rbTot);
                    rawDiagTry = [];
                    rawForces0Try = [];
                    if obj.libOpt('forceOffsetRawDiagnostic',false)
                        rawDiagTry = obj.evaluateScheduleAtState( ...
                            schedRawTry,obj.xFlex,uWing,gust,rbTot);
                        rawForces0Try = struct( ...
                            'rPz',rawDiagTry.rPz, ...
                            'dF_B',rawDiagTry.Ftot_B-oldDiag.Ftot_B, ...
                            'dM_B',rawDiagTry.Mtot_B-oldDiag.Mtot_B);
                    end
                    blockDiagTry = [];
                    if obj.libOpt('blockJumpDiagnostic',false)
                        blockDiagTry = obj.evaluateScheduleBlockJumps( ...
                            oldSched,schedTry,obj.xFlex,uWing,gust,rbTot,oldDiag);
                    end
                    relUDiagTry = [];
                    if obj.libOpt('inputOffsetDiagnostic',true) && ...
                            isfield(schedTry,'u_eq') && ~isempty(schedTry.u_eq)
                        uWingRelNew = obj.scheduleInputIncrement(schedTry,uWing);
                        relUDiagTry = obj.evaluateScheduleAtState( ...
                            schedTry,obj.xFlex,uWingRelNew,gust,rbTot);
                        relUDiagTry.uUsed = uWingRelNew;
                    end
                    affDiagTry = [];
                    if obj.libOpt('affineStateDiagnostic',true) && ...
                            isfield(oldSched,'x_eq') && isfield(schedTry,'x_eq') && ...
                            ~isempty(oldSched.x_eq) && ~isempty(schedTry.x_eq) && ...
                            numel(oldSched.x_eq) == numel(obj.xFlex) && ...
                            numel(schedTry.x_eq) == numel(obj.xFlex)
                        xAff = obj.xFlex(:)+oldSched.x_eq(:)-schedTry.x_eq(:);
                        affDiagTry = obj.evaluateScheduleAtState( ...
                            schedTry,xAff,uWing,gust,rbTot);
                    end
                    gateTry = obj.assessScheduleCandidate( ...
                        oldSched,schedTry,oldDiag,newDiagTry,affDiagTry, ...
                        relUDiagTry,transitionTry);
                    gateTry.blockDiag = blockDiagTry;
                    if ~isempty(rawForces0Try)
                        gateTry.rawForces0 = rawForces0Try;
                    end
                    trialReason(trialIndex) = string(gateTry.reason);
                    trialMaximumLoadRatio(trialIndex) = ...
                        obj.scheduleGateLoadRatio(gateTry);
                    trialAccepted(trialIndex) = logical(gateTry.accept);
                    if gateTry.accept
                        candidate.accepted = true;
                        candidate.mu = muTry;
                        candidate.schedNew = schedTry;
                        candidate.schedRaw = schedRawTry;
                        candidate.transition = transitionTry;
                        candidate.xCandidate = xTry;
                        candidate.newDiag = newDiagTry;
                        candidate.rawDiag = rawDiagTry;
                        candidate.rawForces0 = rawForces0Try;
                        candidate.blockDiag = blockDiagTry;
                        candidate.affDiag = affDiagTry;
                        candidate.relUDiag = relUDiagTry;
                        candidate.gate = gateTry;
                        break
                    end
                catch exception
                    trialReason(trialIndex) = ...
                        "construction_failed:"+string(exception.identifier);
                end
            end
            keep = 1:completedTrials;
            candidate.info = struct( ...
                'enabled',true,'rawTarget',rawMu,'origin',originMu, ...
                'rawReason',string(rawGate.reason), ...
                'acceptedCoordinate',candidate.mu, ...
                'acceptedScale',double(candidate.accepted)* ...
                    trialScale(completedTrials), ...
                'completedTrials',completedTrials, ...
                'maximumTrials',maximumTrials,'contraction',contraction, ...
                'trialScale',trialScale(keep), ...
                'trialMu',trialMu(keep,:), ...
                'trialReason',trialReason(keep), ...
                'trialMaximumLoadRatio',trialMaximumLoadRatio(keep), ...
                'trialAccepted',trialAccepted(keep), ...
                'coordinateLag',rawMu-candidate.mu);
            candidate.gate.scheduleBacktracking = candidate.info;
        end

        function ratio = scheduleGateLoadRatio(~,gate)
        %SCHEDULEGATELOADRATIO Maximum normalized protected-load jump.
            ratio = NaN;
            if isfield(gate,'loadJumpLimit') && ...
                    numel(gate.loadJumpLimit) >= 3 && ...
                    isfield(gate,'dF_B') && numel(gate.dF_B) >= 3 && ...
                    isfield(gate,'dM_B') && numel(gate.dM_B) >= 2
                jump = abs([gate.dF_B(1);gate.dF_B(3);gate.dM_B(2)]);
                ratio = max(jump./gate.loadJumpLimit(:));
            end
        end

        function transition = buildAtomicScheduleTransition( ...
                obj,oldSched,newSched,rawNewSched,newMu,uWing,gust,rbTot,oldDiag)
        %BUILDATOMICSCHEDULETRANSITION Qualify one old-to-new chart change.
            transition = struct('enabled',true,'valid',false, ...
                'reason',"incomplete",'transform',[], ...
                'newState',obj.xFlex,'condition',inf, ...
                'physicalPreservationRelativeError',inf, ...
                'roundTripRelativeError',inf, ...
                'equilibriumResidual',inf, ...
                'physicalDerivativeFirstOrderVariation',inf, ...
                'physicalDerivativeCurvature',inf, ...
                'loadFirstOrderVariation',[inf;inf;inf], ...
                'loadCurvature',[inf;inf;inf]);
            try
                [toNew,newMap] = obj.scheduleStateMap(oldSched,newSched);
                xNew = toNew*obj.xFlex(:);
                midMu = 0.5*(double(oldSched.mu(:).')+double(newMu(:).'));
                schedMid = AeroFlex.sched.evalLibrary( ...
                    obj.ROMlib,midMu,obj.cfg.library);
                schedMid = obj.applyScheduleForceOffsetPolicy( ...
                    oldSched,schedMid);
                [toMid,midMap] = obj.scheduleStateMap(oldSched,schedMid);
                xMid = toMid*obj.xFlex(:);
                midDiag = obj.evaluateScheduleAtState( ...
                    schedMid,xMid,uWing,gust,rbTot);
                newDiag = obj.evaluateScheduleAtState( ...
                    newSched,xNew,uWing,gust,rbTot);
                equilibriumDiag = obj.evaluateScheduleAtState( ...
                    rawNewSched,rawNewSched.x_eq(:),rawNewSched.u_eq(:), ...
                    zeros(obj.nw,1),rbTot);
                oldPhysical = obj.schedulePhysicalDerivative(oldSched,oldDiag);
                midPhysical = obj.schedulePhysicalDerivative(schedMid,midDiag);
                newPhysical = obj.schedulePhysicalDerivative(newSched,newDiag);
                oldLoad = [oldDiag.Ftot_B(1);oldDiag.Ftot_B(3);oldDiag.Mtot_B(2)];
                midLoad = [midDiag.Ftot_B(1);midDiag.Ftot_B(3);midDiag.Mtot_B(2)];
                newLoad = [newDiag.Ftot_B(1);newDiag.Ftot_B(3);newDiag.Mtot_B(2)];
                transition.transform = toNew;
                transition.newState = xNew;
                transition.condition = max(newMap.condition,midMap.condition);
                transition.physicalPreservationRelativeError = max( ...
                    newMap.physicalPreservationRelativeError, ...
                    midMap.physicalPreservationRelativeError);
                transition.roundTripRelativeError = max( ...
                    newMap.roundTripRelativeError,midMap.roundTripRelativeError);
                transition.equilibriumResidual = equilibriumDiag.rPz;
                transition.physicalDerivativeFirstOrderVariation = ...
                    norm(newPhysical-oldPhysical,inf);
                transition.physicalDerivativeCurvature = ...
                    norm(newPhysical-2*midPhysical+oldPhysical,inf);
                transition.loadFirstOrderVariation = abs(newLoad-oldLoad);
                transition.loadCurvature = abs(newLoad-2*midLoad+oldLoad);
                conditionLimit = obj.libOpt( ...
                    'stateTransportConditionLimit',1e6);
                physicalTolerance = obj.libOpt( ...
                    'stateTransportPhysicalTolerance',5e-5);
                roundTripTolerance = obj.libOpt( ...
                    'stateTransportRoundTripTolerance',5e-4);
                transition.valid = transition.condition <= conditionLimit && ...
                    transition.physicalPreservationRelativeError <= ...
                        physicalTolerance && ...
                    transition.roundTripRelativeError <= roundTripTolerance && ...
                    all(isfinite(xNew));
                if transition.valid
                    transition.reason = "qualified";
                else
                    transition.reason = "mapping_tolerance";
                end
            catch exception
                transition.reason = "construction_failed:"+string(exception.identifier);
                transition.errorMessage = string(exception.message);
            end
        end

        function [transform,diagnostic] = scheduleStateMap(obj,from,to)
        %SCHEDULESTATEMAP Map a reduced state through physical/full owners.
            required = {'q1ToPhysical','q2ToPhysical','qxiToPhysical', ...
                'qGammaToFull','qGammaFromFull'};
            assert(isfield(from,'scheduleStateCoordinate') && ...
                isfield(to,'scheduleStateCoordinate'), ...
                'PlantRunTime:ScheduleStateCoordinate', ...
                'Both schedule packages must provide state-coordinate maps.');
            fromMap = from.scheduleStateCoordinate;
            toMap = to.scheduleStateCoordinate;
            assert(all(isfield(fromMap,required)) && all(isfield(toMap,required)), ...
                'PlantRunTime:ScheduleStateCoordinate', ...
                'A schedule state-coordinate contract is incomplete.');
            stateIndex = from.idx;
            assert(isequal(stateIndex,to.idx), ...
                'PlantRunTime:ScheduleStateOrdering', ...
                'Scheduled packages changed the flexible-state ordering.');
            transform = speye(obj.nx);
            blocks = cell(4,1);
            blocks{1} = toMap.q1ToPhysical\fromMap.q1ToPhysical;
            blocks{2} = toMap.q2ToPhysical\fromMap.q2ToPhysical;
            blocks{3} = toMap.qxiToPhysical\fromMap.qxiToPhysical;
            blocks{4} = toMap.qGammaFromFull*fromMap.qGammaToFull;
            groups = {stateIndex.q1,stateIndex.q2, ...
                stateIndex.qxi,stateIndex.qGam};
            for groupIndex = 1:numel(groups)
                transform(groups{groupIndex},groups{groupIndex}) = ...
                    blocks{groupIndex};
            end
            reverse = speye(obj.nx);
            reverseBlocks = { ...
                fromMap.q1ToPhysical\toMap.q1ToPhysical, ...
                fromMap.q2ToPhysical\toMap.q2ToPhysical, ...
                fromMap.qxiToPhysical\toMap.qxiToPhysical, ...
                fromMap.qGammaFromFull*toMap.qGammaToFull};
            for groupIndex = 1:numel(groups)
                reverse(groups{groupIndex},groups{groupIndex}) = ...
                    reverseBlocks{groupIndex};
            end
            state = obj.xFlex(:);
            mapped = transform*state;
            errors = [ ...
                norm(toMap.q1ToPhysical*mapped(stateIndex.q1)- ...
                    fromMap.q1ToPhysical*state(stateIndex.q1),inf), ...
                norm(toMap.q2ToPhysical*mapped(stateIndex.q2)- ...
                    fromMap.q2ToPhysical*state(stateIndex.q2),inf), ...
                norm(toMap.qxiToPhysical*mapped(stateIndex.qxi)- ...
                    fromMap.qxiToPhysical*state(stateIndex.qxi),inf), ...
                norm(toMap.qGammaToFull*mapped(stateIndex.qGam)- ...
                    fromMap.qGammaToFull*state(stateIndex.qGam),inf)];
            scales = [ ...
                norm(fromMap.q1ToPhysical*state(stateIndex.q1),inf), ...
                norm(fromMap.q2ToPhysical*state(stateIndex.q2),inf), ...
                norm(fromMap.qxiToPhysical*state(stateIndex.qxi),inf), ...
                norm(fromMap.qGammaToFull*state(stateIndex.qGam),inf)];
            blockConditions = cellfun(@(value) cond(full(value)),blocks);
            diagnostic = struct( ...
                'condition',max(blockConditions), ...
                'physicalPreservationRelativeError', ...
                    max(errors./max(1,scales)), ...
                'roundTripRelativeError', ...
                    norm(reverse*mapped-state,inf)/max(1,norm(state,inf)));
        end

        function value = schedulePhysicalDerivative(~,sched,diagnostic)
        %SCHEDULEPHYSICALDERIVATIVE Recover the constrained q1 vector field.
            Pz = sched.beam.Pz;
            H1 = sched.scheduleStateCoordinate.q1ToPhysical;
            value = H1*(Pz*diagnostic.q1raw);
        end

        % function gate = assessScheduleCandidate(obj, oldSched, newSched, oldD, newD, affD)
        function gate = assessScheduleCandidate(obj, oldSched, newSched, ...
                oldD, newD, affD, relUD, transition)
        %ASSESSSCHEDULECANDIDATE Decide whether a runtime schedule switch is safe.
        
            gate = struct();
            gate.accept = true;
            gate.reason = "accepted";
        
            dF = newD.Ftot_B - oldD.Ftot_B;
            dM = newD.Mtot_B - oldD.Mtot_B;
        
            gate.dF_B = dF;
            gate.dM_B = dM;
            gate.rPzOld = oldD.rPz;
            gate.rPzNew = newD.rPz;
            gate.rPzGrowth = newD.rPz/max(oldD.rPz, eps);
        
            if ~isempty(affD)
                gate.rPzNewAffine = affD.rPz;
                gate.affineBetterRatio = affD.rPz/max(newD.rPz, eps);
            else
                gate.rPzNewAffine = NaN;
                gate.affineBetterRatio = NaN;
            end
        
            if nargin < 7
                relUD = [];
            end
            if nargin < 8 || isempty(transition)
                transition = struct('enabled',false);
            end
            
            if ~isempty(relUD)
                gate.rPzNewRelU = relUD.rPz;
                gate.relUBetterRatio = relUD.rPz/max(newD.rPz, eps);
                gate.dFRelU_B = relUD.Ftot_B - oldD.Ftot_B;
                gate.dMRelU_B = relUD.Mtot_B - oldD.Mtot_B;
            else
                gate.rPzNewRelU = NaN;
                gate.relUBetterRatio = NaN;
                gate.dFRelU_B = [NaN; NaN; NaN];
                gate.dMRelU_B = [NaN; NaN; NaN];
            end

            gate.dL = norm(oldSched.L - newSched.L,'fro');
        
            gate.dPz  = NaN;
            gate.dPr  = NaN;
            gate.dPhi = NaN;
            
            PzOld = obj.getNestedMember(oldSched, {'beam','Pz'}, []);
            PzNew = obj.getNestedMember(newSched, {'beam','Pz'}, []);
            
            PrOld = obj.getNestedMember(oldSched, {'beam','Pr'}, []);
            PrNew = obj.getNestedMember(newSched, {'beam','Pr'}, []);
            
            PhiOld = obj.getNestedMember(oldSched, {'beam','red','phi1_sA'}, []);
            PhiNew = obj.getNestedMember(newSched, {'beam','red','phi1_sA'}, []);
            
            if ~isempty(PzOld) && ~isempty(PzNew) && isequal(size(PzOld), size(PzNew))
                gate.dPz = norm(PzOld - PzNew, 'fro');
            end
            
            if ~isempty(PrOld) && ~isempty(PrNew) && isequal(size(PrOld), size(PrNew))
                gate.dPr = norm(PrOld - PrNew, 'fro');
            end
            
            if ~isempty(PhiOld) && ~isempty(PhiNew) && isequal(size(PhiOld), size(PhiNew))
                gate.dPhi = norm(PhiOld - PhiNew, 'fro');
            end
        
            gate.dForces0 = NaN;
            if isfield(oldSched.parConst,'forces_0') && isfield(newSched.parConst,'forces_0')
                gate.dForces0 = norm(oldSched.parConst.forces_0 - newSched.parConst.forces_0);
            end
        
            % gate.dBdel = NaN;
            % gate.dDdel = NaN;
            % if isfield(oldSched.parConst,'Bdel') && isfield(newSched.parConst,'Bdel')
            %     gate.dBdel = norm(oldSched.parConst.Bdel - newSched.parConst.Bdel,'fro');
            % end
            % if isfield(oldSched.parConst,'Ddel') && isfield(newSched.parConst,'Ddel')
            %     gate.dDdel = norm(oldSched.parConst.Ddel - newSched.parConst.Ddel,'fro');
            % end
            gate.dBdel = NaN;
            gate.dDdel = NaN;
            gate.rBdel = NaN;
            gate.rDdel = NaN;
            
            if isfield(oldSched.parConst,'Bdel') && isfield(newSched.parConst,'Bdel')
                dB = oldSched.parConst.Bdel - newSched.parConst.Bdel;
                gate.dBdel = norm(dB,'fro');
                gate.rBdel = gate.dBdel / max(norm(oldSched.parConst.Bdel,'fro'), eps);
            end
            
            if isfield(oldSched.parConst,'Ddel') && isfield(newSched.parConst,'Ddel')
                dD = oldSched.parConst.Ddel - newSched.parConst.Ddel;
                gate.dDdel = norm(dD,'fro');
                gate.rDdel = gate.dDdel / max(norm(oldSched.parConst.Ddel,'fro'), eps);
            end
        
            % pzAbsTol = obj.libOpt('admissiblePzAbsTol', 5e-5);
            % pzGrowth = obj.libOpt('admissiblePzGrowth', 5.0);
            pzAbsTol = obj.libOpt('admissiblePzAbsTol', 5e-3);
            pzGrowth = obj.libOpt('admissiblePzGrowth', 10.0);

            pzLimit = max(pzAbsTol, pzGrowth*max(oldD.rPz, eps));
            gate.pzLimit = pzLimit;

            if transition.enabled
                gate.stateTransport = rmfield(transition,intersect( ...
                    fieldnames(transition),{'transform','newState'}));
                gate.rPzEquilibrium = transition.equilibriumResidual;
                gate.physicalDerivativeCurvature = ...
                    transition.physicalDerivativeCurvature;
                if ~transition.valid
                    gate.accept = false;
                    gate.reason = "state_transport_invariant";
                    return
                end
                if transition.equilibriumResidual > pzAbsTol || ...
                        transition.physicalDerivativeCurvature > pzAbsTol
                    gate.accept = false;
                    gate.reason = "transported_vector_field_curvature";
                    return
                end
            elseif newD.rPz > pzLimit
                gate.accept = false;
                gate.reason = "projected_residual_jump";
                return
            end
        
            % Schedule switches are judged by the jump they introduce in the current
            % dynamic loads.  A trim-level absolute tolerance is too strict once the
            % aircraft is already responding to a gust, so use an absolute floor plus a
            % relative allowance based on the current load scale.
            
            loadJumpAbsTol = obj.libOpt('loadJumpTol', [1e-3; 1e-3; 5e-4]);
            loadJumpAbsTol = loadJumpAbsTol(:);
            
            if numel(loadJumpAbsTol) < 3
                loadJumpAbsTol = [1e-3; 1e-3; 5e-4];
            end
            
            loadJumpRelTol = obj.libOpt('loadJumpRelTol', [0.02; 0.02; 0.02]);
            loadJumpRelTol = loadJumpRelTol(:);
            
            if numel(loadJumpRelTol) < 3
                loadJumpRelTol = [0.02; 0.02; 0.02];
            end
            
            loadJumpRefFloor = obj.libOpt('loadJumpRefFloor', [1.0; 1.0; 0.10]);
            loadJumpRefFloor = loadJumpRefFloor(:);
            
            if numel(loadJumpRefFloor) < 3
                loadJumpRefFloor = [1.0; 1.0; 0.10];
            end
            
            jump = [abs(dF(1)); abs(dF(3)); abs(dM(2))];
            if transition.enabled
                % Atomic chart transport does not waive the applied-load
                % jump gate.  Enforce the larger of the physical old/new
                % change and the midpoint curvature using the same locked
                % component limits.
                jump = max(jump,transition.loadCurvature);
            end
            
            loadRef = max(abs([oldD.Ftot_B(1); oldD.Ftot_B(3); oldD.Mtot_B(2)]), ...
                          loadJumpRefFloor(1:3));
            
            loadJumpLimit = loadJumpAbsTol(1:3) + loadJumpRelTol(1:3).*loadRef;
            
            gate.loadJumpTol = loadJumpAbsTol(1:3);
            gate.loadJumpRelTol = loadJumpRelTol(1:3);
            gate.loadJumpLimit = loadJumpLimit;
            
            if any(jump > loadJumpLimit)
                gate.accept = false;
                gate.reason = "load_jump";
                return
            end
        
            % Simplex changes are not forbidden, but they are useful to print.
            gate.simplexChanged = false;
            if isfield(oldSched,'pointIds') && isfield(newSched,'pointIds')
                gate.simplexChanged = ~isequal(oldSched.pointIds(:), newSched.pointIds(:));
            end
        end
        
        function tf = shouldPrintScheduleGate(obj, accepted)
        %SHOULDPRINTSCHEDULEGATE Rate-limit schedule diagnostics.
        
            if ~isfield(obj.cfg,'debug') || ~isfield(obj.cfg.debug,'schedule') || ...
                    ~logical(obj.cfg.debug.schedule)
                tf = false;
                return
            end
        
            if accepted
                every = obj.libOpt('scheduleAcceptedPrintEvery', 25);
                tf = every > 0 && (mod(obj.k, every) == 0 || obj.k == 0);
            else
                maxRejectPrint = obj.libOpt('scheduleRejectPrintMax', 12);
                if ~isfield(obj.last,'schedRejectedCount') || isempty(obj.last.schedRejectedCount)
                    tf = true;
                else
                    tf = obj.last.schedRejectedCount <= maxRejectPrint || ...
                         mod(obj.last.schedRejectedCount, 25) == 0;
                end
            end
        end
        
        function printScheduleGate(obj, decision, mu, oldSched, newSched, oldD, newD, affD, gate)
        % function printScheduleGate(obj, decision, mu, oldSched, newSched, oldD, newD, affD, gate)
        %PRINTSCHEDULEGATE Compact diagnostics for schedule acceptance/rejection.
        
            fprintf('\n[sched:%s]\n', decision);
            fprintf('  step=%d t=%.6f mu=[%.6f %.6f]\n', obj.k, obj.t, mu(1), mu(2));
        
            if isfield(oldSched,'mu') && ~isempty(oldSched.mu)
                fprintf('  dmu=[%+.3e %+.3e]\n', mu(1)-oldSched.mu(1), mu(2)-oldSched.mu(2));
            end
        
            if isfield(newSched,'pointIds') && isfield(newSched,'weights')
                fprintf('  ids=[%s] w=[%s]\n', ...
                    num2str(newSched.pointIds(:).'), ...
                    num2str(newSched.weights(:).',' %.4f'));
            end
        
            fprintf('  reason=%s\n', string(gate.reason));
            fprintf('  dL=%.3e dPz=%.3e dPr=%.3e dPhi=%.3e\n', ...
                gate.dL, gate.dPz, gate.dPr, gate.dPhi);
            fprintf('  dForces0=%.3e dBdel=%.3e dDdel=%.3e\n', ...
                gate.dForces0, gate.dBdel, gate.dDdel);
        
            fprintf('  |Pz*q1raw| old=%.3e new=%.3e limit=%.3e growth=%.3e\n', ...
                gate.rPzOld, gate.rPzNew, gate.pzLimit, gate.rPzGrowth);
        
            if ~isempty(affD)
                fprintf('  affine diagnostic: |Pz*q1raw|=%.3e ratio=%.3e\n', ...
                    gate.rPzNewAffine, gate.affineBetterRatio);
            end
            if isfield(gate,'rPzNewRelU') && isfinite(gate.rPzNewRelU)
                fprintf('  input-offset diagnostic: |Pz*q1raw|=%.3e ratio=%.3e\n', ...
                    gate.rPzNewRelU, gate.relUBetterRatio);
                fprintf('  input-offset load jump: dFx=%+.3e dFz=%+.3e dMy=%+.3e\n', ...
                    gate.dFRelU_B(1), gate.dFRelU_B(3), gate.dMRelU_B(2));
            end
        
            fprintf('  load jump: dFx=%+.3e dFz=%+.3e dMy=%+.3e\n', ...
                gate.dF_B(1), gate.dF_B(3), gate.dM_B(2));
            
            if isfield(gate,'loadJumpLimit') && ~isempty(gate.loadJumpLimit)
                fprintf('  load limits: Fx=%.3e Fz=%.3e My=%.3e\n', ...
                    gate.loadJumpLimit(1), gate.loadJumpLimit(2), gate.loadJumpLimit(3));
            end

            if isfield(gate,'rawForces0') && ~isempty(gate.rawForces0)
                fprintf('  raw forces_0 candidate: |Pz*q1raw|=%.3e dFx=%+.3e dFz=%+.3e dMy=%+.3e\n', ...
                    gate.rawForces0.rPz, ...
                    gate.rawForces0.dF_B(1), ...
                    gate.rawForces0.dF_B(3), ...
                    gate.rawForces0.dM_B(2));
            end

            if isfield(gate,'blockDiag') && ~isempty(gate.blockDiag)
                obj.printScheduleBlockJumps(gate.blockDiag);
            end

            fprintf('  old loads: Fx=%+.3e Fz=%+.3e My=%+.3e\n', ...
                oldD.Ftot_B(1), oldD.Ftot_B(3), oldD.Mtot_B(2));
        
            fprintf('  new loads: Fx=%+.3e Fz=%+.3e My=%+.3e\n\n', ...
                newD.Ftot_B(1), newD.Ftot_B(3), newD.Mtot_B(2));
        end
        
        function val = libOpt(obj, name, defaultVal)
        %LIBOPT Read cfg.library option with default.
        
            val = defaultVal;
        
            if isfield(obj.cfg,'library') && isstruct(obj.cfg.library) && ...
                    isfield(obj.cfg.library,name) && ~isempty(obj.cfg.library.(name))
                val = obj.cfg.library.(name);
            end
        end
        function uTotal = composeWingControl(obj, uIncrement)
        %COMPOSEWINGCONTROL Add trim wing command to control increments.
        
            uIncrement = obj.sanitizeControl(uIncrement);
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'commandsAreTrimIncrements') ...
                    && ~obj.cfg.sim.commandsAreTrimIncrements
                uTotal = uIncrement;
                return
            end
        
            uTotal = obj.trimWingControl(:) + uIncrement(:);
        end
        
        function rbTot = composeRigidCommand(obj, rbCmd)
        %COMPOSERIGIDCOMMAND Add trim elevator and thrust to rigid increments.
        
            if nargin < 2 || isempty(rbCmd)
                rbCmd = struct();
            end
        
            if ~isfield(rbCmd,'delta_e') || isempty(rbCmd.delta_e)
                rbCmd.delta_e = 0;
            end
        
            if ~isfield(rbCmd,'delta_a') || isempty(rbCmd.delta_a)
                rbCmd.delta_a = 0;
            end
        
            if ~isfield(rbCmd,'delta_r') || isempty(rbCmd.delta_r)
                rbCmd.delta_r = 0;
            end
        
            if ~isfield(rbCmd,'thrust') || isempty(rbCmd.thrust)
                rbCmd.thrust = 0;
            end
        
            rbTot = rbCmd;
            % For future trajecotry pass something like:
            % rbCmd.alphaScheduleDeg = alpha_ref_deg;
            % for
            % xk = plant.stepCoupled(uk, gk, rbCmd);

            rbTot.delta_e = obj.trimTailDelta + rbCmd.delta_e;
            rbTot.thrust  = obj.trimThrust    + rbCmd.thrust;
        end
        
        function loads = computeCoupledLoads(obj, xFlex, uWing, gust, rbCmd)
        %COMPUTECOUPLEDLOADS Compute all rigid-body force/moment terms.
            
            [Clamp6, Fwing_B, Mwing_B] = obj.computeWingClampReaction(xFlex, uWing, gust);
        
            [Ftail_B, Mtail_B] = obj.computeTailLoadsFromCmd(rbCmd);
            [Ffin_B,  Mfin_B ] = obj.computeFinLoadsFromCmd(rbCmd);
            Fgrav_B = obj.computeGravityBody();
            [Fthrust_B, Mthrust_B] = obj.computeThrustLoadsFromCmd(rbCmd);
        
            loads = struct();
            loads.Clamp6    = Clamp6;
            loads.Fwing_B   = Fwing_B;
            loads.Mwing_B   = Mwing_B;
            loads.Ftail_B   = Ftail_B(:);
            loads.Mtail_B   = Mtail_B(:);
            loads.Ffin_B    = Ffin_B(:);
            loads.Mfin_B    = Mfin_B(:);
            loads.Fgrav_B   = Fgrav_B(:);
            loads.Fthrust_B = Fthrust_B(:);
            loads.Mthrust_B = Mthrust_B(:);
            loads.Freciprocal_B = zeros(3,1);
            loads.Mreciprocal_B = zeros(3,1);
        
            loads.Ftot_B = loads.Fwing_B + loads.Ftail_B + loads.Ffin_B + ...
                           loads.Fgrav_B + loads.Fthrust_B;
        
            loads.Mtot_B = loads.Mwing_B + loads.Mtail_B + loads.Mfin_B + ...
                           loads.Mthrust_B;
        end
        
        function uTrim = getTrimWingControl(obj, trim)
        %GETTRIMWINGCONTROL Return [delta_1; delta_2; ddelta_1; ddelta_2] at trim.
        %
        % trim.deltaDeg is historical naming. Values are treated as radians.
        
            uTrim = zeros(obj.nu,1);
        
            if nargin < 2 || ~isstruct(trim)
                return
            end
        
            if isfield(trim,'deltaWing')
                d = trim.deltaWing(:);
            elseif isfield(trim,'deltaDeg')
                d = trim.deltaDeg(:);
            else
                return
            end
        
            n = min(numel(uTrim), numel(d));
            uTrim(1:n) = d(1:n);
        end
        
        function [Fwing_B, Mwing_B] = mapClampToRigidLoads(obj, Clamp6)
        %MAPCLAMPTORIGIDLOADS Convert clamp resultant to rigid-body force/moment.
        
            Fwing_B = Clamp6(1:3);
            Mwing_B = Clamp6(4:6);
        
            useShift = false;
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'shiftWingMomentToCG')
                useShift = logical(obj.cfg.sim.shiftWingMomentToCG);
            elseif isfield(obj.cfg,'trim') && isfield(obj.cfg.trim,'shiftWingMomentToCG')
                useShift = logical(obj.cfg.trim.shiftWingMomentToCG);
            end
        
            if useShift
                if isfield(obj.rbParams,'rWingRoot_B')
                    % The non-wing owner defines this vector from the rigid
                    % center of mass to the wing-root wrench application point.
                    Mwing_B = Mwing_B + cross( ...
                        obj.rbParams.rWingRoot_B(:),Fwing_B);
                elseif isfield(obj.cfg,'geom') && isfield(obj.cfg.geom,'rWingClamp_B') ...
                        && isfield(obj.cfg.geom,'rCG_B')
                    rA_CG_B = obj.cfg.geom.rWingClamp_B(:) - obj.cfg.geom.rCG_B(:);
                    Mwing_B = Mwing_B + cross(rA_CG_B, Fwing_B);
                elseif isfield(obj.rbParams,'rWingClamp_B')
                    Mwing_B = Mwing_B + cross(obj.rbParams.rWingClamp_B(:), Fwing_B);
                end
            end
        end

        function gust = centeredRecoveryGust(obj,gust)
        %CENTEREDRECOVERYGUST Preserve default recovery pending audit evidence.
            enabled = isfield(obj.cfg,'sim') && ...
                isfield(obj.cfg.sim,'centeredRootWrenchGustRecoveryAudit') && ...
                isstruct(obj.cfg.sim.centeredRootWrenchGustRecoveryAudit) && ...
                isfield(obj.cfg.sim.centeredRootWrenchGustRecoveryAudit,'enabled') && ...
                logical(obj.cfg.sim.centeredRootWrenchGustRecoveryAudit.enabled);
            if ~enabled
                gust = [];
                return
            end
            request = obj.cfg.sim.centeredRootWrenchGustRecoveryAudit;
            assert(isfield(request,'auditOnly') && logical(request.auditOnly) && ...
                isfield(request,'source') && string(request.source) == ...
                "phase18c-v17a-centered-root-wrench-gust-recovery-audit-v1", ...
                'PlantRunTime:CenteredRecoveryGustAudit', ...
                'Centered root-wrench gust recovery requires its audit-only request.');
            gust = obj.sanitizeDisturbance(gust);
        end
        
        function tf = shouldMirrorWingClamp(obj)
        %SHOULDMIRRORWINGCLAMP Return true for half-wing-to-full-aircraft mapping.
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'mirrorWingClamp')
                tf = logical(obj.cfg.sim.mirrorWingClamp);
            elseif isfield(obj.cfg,'trim') && isfield(obj.cfg.trim,'mirrorWingClamp')
                tf = logical(obj.cfg.trim.mirrorWingClamp);
            else
                tf = false;
            end
        end
        function printTrimReplayComparison(obj, loads0)
        %PRINTTRIMREPLAYCOMPARISON Compare trim and runtime load decomposition.
        
            if ~isstruct(obj.trim) || ~isfield(obj.trim,'debug') || ...
                    ~isfield(obj.trim.debug,'loads')
                return
            end
        
            Ltr = obj.trim.debug.loads;
        
            fprintf('\n[trim replay comparison]\n');
        
            printVec('Fwing',   loads0.Fwing_B,   Ltr.Fwing_B);
            printVec('Ftail',   loads0.Ftail_B,   Ltr.Ftail_B);
            % printVec('Ffin',    loads0.Ffin_B,    Ltr.Ffin_B);
            printVec('Fgrav',   loads0.Fgrav_B,   Ltr.Fgrav_B);
            printVec('Fthrust', loads0.Fthrust_B, Ltr.Fthrust_B);
            printVec('Ftot',    loads0.Ftot_B,    Ltr.Ftot_B);
        
            printVec('Mwing',   loads0.Mwing_B,   Ltr.Mwing_B);
            printVec('Mtail',   loads0.Mtail_B,   Ltr.Mtail_B);
            % printVec('Mfin',    loads0.Mfin_B,    Ltr.Mfin_B);
            printVec('Mthrust', loads0.Mthrust_B, Ltr.Mthrust_B);
            printVec('Mtot',    loads0.Mtot_B,    Ltr.Mtot_B);
        
            fprintf('\n');
        
            function printVec(name, runVal, trimVal)
                runVal = runVal(:);
                trimVal = trimVal(:);
                d = runVal - trimVal;
        
                fprintf('  %-8s run=[%+.3e %+.3e %+.3e] trim=[%+.3e %+.3e %+.3e] diff=[%+.3e %+.3e %+.3e]\n', ...
                    name, runVal(1),runVal(2),runVal(3), ...
                    trimVal(1),trimVal(2),trimVal(3), ...
                    d(1),d(2),d(3));
            end
        end

        function uLocal = scheduleInputIncrement(~, schedUse, uWing)
        %SCHEDULEINPUTINCREMENT Candidate local control increment.
        %
        % The committed dynamics are not changed by this helper. It is used first
        % as a diagnostic to determine whether scheduled force maps expect control
        % increments about sched.u_eq.
        
            uLocal = uWing(:);
        
            if ~isfield(schedUse,'u_eq') || isempty(schedUse.u_eq)
                return
            end
        
            uEq = schedUse.u_eq(:);
            n = min(numel(uLocal), numel(uEq));
        
            uLocal(1:n) = uLocal(1:n) - uEq(1:n);
        end

        function mu = filteredSchedulePoint(obj, muRaw, rbTot)
        %FILTEREDSCHEDULEPOINT Runtime scheduling point for the ROM library.
        %
        % The library alpha is an operating-point coordinate, not necessarily the
        % instantaneous disturbed body alpha.  For commanded trajectories, use the
        % reference schedule.  For gust-only validation, use trimHold.  For slow
        % open-loop manoeuvres, filtered can be used cautiously.
        
            if nargin < 3 || isempty(rbTot)
                rbTot = struct();
            end
        
            muRaw = double(muRaw(:).');
        
            if numel(muRaw) ~= 2
                error('PlantRunTime:BadSchedulePoint', ...
                    'Expected muRaw=[U_inf, alpha_deg]. Got %d entries.', numel(muRaw));
            end
        
            if ~isfield(obj.last,'sched_mu_anchor') || isempty(obj.last.sched_mu_anchor)
                if isfield(obj.last,'sched_mu') && ~isempty(obj.last.sched_mu)
                    obj.last.sched_mu_anchor = obj.last.sched_mu(:).';
                else
                    obj.last.sched_mu_anchor = muRaw;
                end
            end
        
            if ~isfield(obj.last,'sched_mu_filter') || isempty(obj.last.sched_mu_filter)
                obj.last.sched_mu_filter = obj.last.sched_mu_anchor;
            end
        
            muFilt = obj.last.sched_mu_filter;
        
            % U scheduling.  An explicit audit reference is used only when
            % the caller supplies it through the rigid-command metadata.
            % The normal path continues to use the current flight condition.
            scheduleU = muRaw(1);
            if isstruct(rbTot) && isfield(rbTot,'scheduleUInfMps') && ...
                    ~isempty(rbTot.scheduleUInfMps)
                scheduleU = double(rbTot.scheduleUInfMps);
                assert(isscalar(scheduleU) && isfinite(scheduleU) && ...
                    scheduleU > 0, 'PlantRunTime:ScheduleReferenceSpeed', ...
                    'scheduleUInfMps must be a positive finite scalar.');
            end

            % tauU=0 means use the selected scheduling speed directly.
            tauU = obj.libOpt('scheduleTauU', 0.0);
            if tauU > 0
                bU = obj.dt/(tauU + obj.dt);
                muFilt(1) = muFilt(1) + bU*(scheduleU - muFilt(1));
            else
                muFilt(1) = scheduleU;
            end
            
            %   Alpha scheduling modes:
            %   trimHold : hold alpha at the trim/initial schedule point
            %   filtered : low-pass filter rigid-body alpha
            %   instant  : use instantaneous rigid-body alpha

            mode = lower(string(obj.libOpt('scheduleAlphaMode', 'reference')));
            
            switch mode
        
                case {"reference","ref","trajectory","commanded"}
                    alphaRefDeg = obj.getScheduleAlphaReference(rbTot);
        
                    tauA = obj.libOpt('scheduleTauAlpha', 0.10);
                    if tauA > 0
                        bA = obj.dt/(tauA + obj.dt);
                        muFilt(2) = muFilt(2) + bA*(alphaRefDeg - muFilt(2));
                    else
                        muFilt(2) = alphaRefDeg;
                    end
        
                case {"trimhold","trim","frozenalpha","hold"}
                    muFilt(2) = obj.last.sched_mu_anchor(2);
        
                case {"filtered","slow"}
                    tauA = obj.libOpt('scheduleTauAlpha', 0.10);
                    if tauA > 0
                        bA = obj.dt/(tauA + obj.dt);
                        muFilt(2) = muFilt(2) + bA*(muRaw(2) - muFilt(2));
                    else
                        muFilt(2) = muRaw(2);
                    end
        
                case {"instant","raw"}
                    muFilt(2) = muRaw(2);
        
                otherwise
                    error('PlantRunTime:BadScheduleAlphaMode', ...
                        'Unknown cfg.library.scheduleAlphaMode="%s".', mode);
            end
        
            obj.last.sched_mu_raw = muRaw;
            obj.last.sched_mu_filter = muFilt;
        
            mu = muFilt;
        end

        % function mu = filteredSchedulePoint(obj, muRaw)
        % %FILTEREDSCHEDULEPOINT Runtime scheduling point for the ROM library.
        % %
        % % The alpha coordinate in the ROM library is a quasi-steady operating-point
        % % parameter.  It should not normally chase the fast gust-induced rigid-body
        % % alpha.  For gust/open-loop studies, hold alpha at the trim schedule and
        % % allow only U to vary.  For slow manoeuvre studies, use 'filtered'.
        % 
        %     muRaw = double(muRaw(:).');
        % 
        %     if numel(muRaw) ~= 2
        %         error('PlantRunTime:BadSchedulePoint', ...
        %             'Expected muRaw=[U_inf, alpha_deg]. Got %d entries.', numel(muRaw));
        %     end
        % 
        %     % Anchor the schedule at the initial trim schedule.
        %     if ~isfield(obj.last,'sched_mu_anchor') || isempty(obj.last.sched_mu_anchor)
        %         if isfield(obj.last,'sched_mu') && ~isempty(obj.last.sched_mu)
        %             obj.last.sched_mu_anchor = obj.last.sched_mu(:).';
        %         else
        %             obj.last.sched_mu_anchor = muRaw;
        %         end
        %     end
        % 
        %     if ~isfield(obj.last,'sched_mu_filter') || isempty(obj.last.sched_mu_filter)
        %         obj.last.sched_mu_filter = obj.last.sched_mu_anchor;
        %     end
        % 
        %     muFilt = obj.last.sched_mu_filter;
        % 
        %     % U scheduling may remain active.  If tauU=0, U follows the measured value.
        %     tauU = obj.libOpt('scheduleTauU', 0.0);
        %     if tauU > 0
        %         bU = obj.dt/(tauU + obj.dt);
        %         muFilt(1) = muFilt(1) + bU*(muRaw(1) - muFilt(1));
        %     else
        %         muFilt(1) = muRaw(1);
        %     end
        % 
        %     % Alpha scheduling modes:
        %     %   trimHold : hold alpha at the trim/initial schedule point
        %     %   filtered : low-pass filter rigid-body alpha
        %     %   instant  : use instantaneous rigid-body alpha
        %     %
        %     % For gust runs, use trimHold.
        %     mode = char(obj.libOpt('scheduleAlphaMode', 'trimHold'));
        % 
        %     switch lower(mode)
        %         case {'trimhold','trim','frozenalpha','hold'}
        %             muFilt(2) = obj.last.sched_mu_anchor(2);
        % 
        %         case {'filtered','slow'}
        %             tauA = obj.libOpt('scheduleTauAlpha', 0.10);
        %             if tauA > 0
        %                 bA = obj.dt/(tauA + obj.dt);
        %                 muFilt(2) = muFilt(2) + bA*(muRaw(2) - muFilt(2));
        %             else
        %                 muFilt(2) = muRaw(2);
        %             end
        % 
        %         case {'instant','raw'}
        %             muFilt(2) = muRaw(2);
        % 
        %         otherwise
        %             error('PlantRunTime:BadScheduleAlphaMode', ...
        %                 'Unknown cfg.library.scheduleAlphaMode="%s".', mode);
        %     end
        % 
        %     obj.last.sched_mu_raw = muRaw;
        %     obj.last.sched_mu_filter = muFilt;
        % 
        %     mu = muFilt;
        % end

        function blockDiag = evaluateScheduleBlockJumps(obj, oldSched, newSched, x, uWing, gust, rbTot, oldD)
        %EVALUATESCHEDULEBLOCKJUMPS Isolate which scheduled block causes a jump.
        %
        % Each diagnostic evaluates the current state with one group of scheduled
        % fields replaced by the candidate schedule.  This does not change the
        % committed model.
        
            blockDiag = struct( ...
                'name', {}, ...
                'rPz', {}, ...
                'dFx', {}, ...
                'dFz', {}, ...
                'dMy', {});
        
            if nargin < 8 || isempty(oldD)
                oldD = obj.evaluateScheduleAtState(oldSched, x, uWing, gust, rbTot);
            end
        
            % Linear operator only.
            mix = oldSched;
            if isfield(newSched,'L') && isequal(size(newSched.L), size(oldSched.L))
                mix.L = newSched.L;
                mix = obj.attachProjectionForDiagnostic(mix);
                blockDiag(end+1) = obj.evalScheduleBlockCase('L only', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
            end
            
            % All parConst fields, but keep old L.
            mix = oldSched;
            if isfield(newSched,'parConst')
                mix.parConst = newSched.parConst;
                mix = obj.attachProjectionForDiagnostic(mix);
                blockDiag(end+1) = obj.evalScheduleBlockCase('parConst all', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
            end
        
            % Steady force offset.
            mix = oldSched;
            mix = obj.copyParConstFields(mix, newSched, {'forces_0'});
            mix = obj.attachProjectionForDiagnostic(mix);
            blockDiag(end+1) = obj.evalScheduleBlockCase('forces_0', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
            
            % Nonlinear structural/gravity tensors.
            mix = oldSched;
            mix = obj.copyParConstFields(mix, newSched, {'Gamma1','Gamma2','Gamma_g','Gamma_xi'});
            mix = obj.attachProjectionForDiagnostic(mix);
            blockDiag(end+1) = obj.evalScheduleBlockCase('Gamma tensors', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
        
            % Control-surface force maps.
            mix = oldSched;
            mix = obj.copyParConstFields(mix, newSched, {'Bdel','Ddel','Bddel','Dddel'});
            mix = obj.attachProjectionForDiagnostic(mix);
            blockDiag(end+1) = obj.evalScheduleBlockCase('control maps', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
        
            % Gust input maps.
            mix = oldSched;
            mix = obj.copyParConstFields(mix, newSched, {'Bw','Dw','gust_input'});
            mix = obj.attachProjectionForDiagnostic(mix);
            blockDiag(end+1) = obj.evalScheduleBlockCase('gust maps', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
        
            % Scaling/time constants.
            mix = oldSched;
            mix = obj.copyParConstFields(mix, newSched, {'scaleAero','scaleAero1','Fscale','scaleA','t_inf','dt','Na'});
            mix = obj.attachProjectionForDiagnostic(mix);
            blockDiag(end+1) = obj.evalScheduleBlockCase('scales', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
        
            % Full candidate for comparison.
            mix = obj.attachProjectionForDiagnostic(newSched);
            blockDiag(end+1) = obj.evalScheduleBlockCase('full candidate', mix, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
            % blockDiag(end+1) = obj.evalScheduleBlockCase('full candidate', newSched, x, uWing, gust, rbTot, oldD); %#ok<AGROW>
        end
        
        function entry = evalScheduleBlockCase(obj, name, schedUse, x, uWing, gust, rbTot, oldD)
        %EVALSCHEDULEBLOCKCASE Evaluate one block-swap case.
        
            D = obj.evaluateScheduleAtState(schedUse, x, uWing, gust, rbTot);
        
            entry = struct();
            entry.name = name;
            entry.rPz = D.rPz;
            entry.dFx = D.Ftot_B(1) - oldD.Ftot_B(1);
            entry.dFz = D.Ftot_B(3) - oldD.Ftot_B(3);
            entry.dMy = D.Mtot_B(2) - oldD.Mtot_B(2);
        end
        
        function schedOut = copyParConstFields(~, schedOut, schedIn, fields)
        %COPYPARCONSTFIELDS Copy compatible parConst fields from schedIn to schedOut.
        
            if ~isfield(schedOut,'parConst') || ~isfield(schedIn,'parConst')
                return
            end
        
            for k = 1:numel(fields)
                f = fields{k};
        
                if ~isfield(schedIn.parConst, f)
                    continue
                end
        
                if isfield(schedOut.parConst, f)
                    a = schedOut.parConst.(f);
                    b = schedIn.parConst.(f);
        
                    if isnumeric(a) && isnumeric(b) && ~isequal(size(a), size(b))
                        continue
                    end
                end
        
                schedOut.parConst.(f) = schedIn.parConst.(f);
            end
        end
        
        function printScheduleBlockJumps(obj, blockDiag)
        %PRINTSCHEDULEBLOCKJUMPS Print compact block-jump diagnostics.
        
            if isempty(blockDiag)
                return
            end
        
            nMax = obj.libOpt('blockJumpPrintMax', 6);
        
            dMyAbs = abs([blockDiag.dMy]);
            [~,ord] = sort(dMyAbs, 'desc');
            ord = ord(1:min(nMax, numel(ord)));
        
            fprintf('  block-jump diagnostic, sorted by |dMy|:\n');
            fprintf('    %-16s  |Pz*q1raw|       dFx          dFz          dMy\n', 'block');
        
            for ii = 1:numel(ord)
                k = ord(ii);
                fprintf('    %-16s  %.3e  %+.3e  %+.3e  %+.3e\n', ...
                    blockDiag(k).name, ...
                    blockDiag(k).rPz, ...
                    blockDiag(k).dFx, ...
                    blockDiag(k).dFz, ...
                    blockDiag(k).dMy);
            end
        end

        function schedOut = attachProjectionForDiagnostic(obj, schedOut)
        %ATTACHPROJECTIONFORDIAGNOSTIC Ensure diagnostic schedules carry Pz.
        %
        % interpParConst intentionally does not interpolate RateProject.  Diagnostic
        % block swaps can therefore create a schedule whose parConst.RateProject.Pz is
        % empty.  For diagnostic evaluation, recover Pz from the active beam map.
        
            if ~isfield(schedOut,'parConst') || isempty(schedOut.parConst)
                return
            end
        
            Pz = obj.getNestedMember(schedOut, {'parConst','RateProject','Pz'}, []);
        
            if isempty(Pz)
                Pz = obj.getNestedMember(schedOut, {'beam','Pz'}, []);
            end
        
            if isempty(Pz)
                return
            end
        
            schedOut.parConst.RateProject = struct( ...
                'projSet', true, ...
                'Pz', Pz);
        end

        function [schedNew,info] = resolveScheduledQueryPackage( ...
                obj,runtimeLibrary,mu,runtimeLibraryCfg)
        %RESOLVESCHEDULEDQUERYPACKAGE Reuse an exact audit-local package only.
            info = struct('enabled',false,'hit',false,'key',"", ...
                'lookupSeconds',0,'buildSeconds',0,'hits',0,'misses',0);
            if ~isfield(obj.cfg,'library') || ...
                    ~isfield(obj.cfg.library,'exactPackageReuseAudit') || ...
                    isempty(obj.cfg.library.exactPackageReuseAudit)
                schedNew = AeroFlex.sched.evalLibrary( ...
                    runtimeLibrary,mu,runtimeLibraryCfg);
                return
            end
            request = obj.cfg.library.exactPackageReuseAudit;
            assert(isstruct(request) && isscalar(request) && ...
                isfield(request,'enabled') && logical(request.enabled) && ...
                isfield(request,'auditOnly') && logical(request.auditOnly) && ...
                isfield(request,'source') && string(request.source)== ...
                "phase18c-v17a-caseb-exact-package-reuse-output-decoupling-audit-v1", ...
                'PlantRunTime:ExactPackageReuseScope', ...
                'Exact package reuse is restricted to the approved Case-B audit.');
            maximumEntries = 64;
            if isfield(request,'maximumEntries') && ~isempty(request.maximumEntries)
                maximumEntries = double(request.maximumEntries);
            end
            assert(isscalar(maximumEntries) && isfinite(maximumEntries) && ...
                maximumEntries >= 1 && maximumEntries == round(maximumEntries), ...
                'PlantRunTime:ExactPackageReuseCapacity', ...
                'The exact package-reuse cache capacity must be a positive integer.');
            candidate = runtimeLibraryCfg.fullCoordinateRuntimeCandidate;
            assert(isstruct(candidate) && isfield(candidate,'registryPath') && ...
                isfield(candidate,'fieldRoot'), ...
                'PlantRunTime:ExactPackageReuseCandidate', ...
                'Exact package reuse requires the active full-coordinate candidate.');
            lookupTimer = tic;
            [weights,ids] = AeroFlex.sched.interpWeights( ...
                runtimeLibrary,mu,runtimeLibraryCfg);
            sourceIds = string({runtimeLibrary.points(ids).name});
            registryHash = localExactPackageFileHash(candidate.registryPath);
            assert(strlength(registryHash)>0, ...
                'PlantRunTime:ExactPackageReuseRegistry', ...
                'The full-coordinate registry hash is unavailable.');
            key = localExactPackageReuseKey(mu,ids,weights,sourceIds, ...
                candidate,registryHash);
            cache = localExactPackageCache(obj.last,maximumEntries);
            hitIndex = find(string({cache.entries.key})==key,1,'first');
            info.enabled = true;
            info.key = key;
            info.lookupSeconds = toc(lookupTimer);
            if ~isempty(hitIndex)
                schedNew = cache.entries(hitIndex).package;
                info.hit = true;
                cache.hits = cache.hits+1;
            else
                buildTimer = tic;
                schedNew = AeroFlex.sched.evalLibrary( ...
                    runtimeLibrary,mu,runtimeLibraryCfg);
                info.buildSeconds = toc(buildTimer);
                assert(isequal(double(schedNew.mu(:)),double(mu(:))) && ...
                    isequal(double(schedNew.weights(:)),double(weights(:))) && ...
                    isequal(double(schedNew.pointIds(:)),double(ids(:))), ...
                    'PlantRunTime:ExactPackageReuseIdentity', ...
                    'The cached-query resolver disagreed with evalLibrary.');
                entry = struct('key',key,'package',schedNew);
                if numel(cache.entries) >= maximumEntries
                    cache.entries(1) = [];
                end
                cache.entries(end+1,1) = entry;
                cache.misses = cache.misses+1;
                cache.buildSeconds = cache.buildSeconds+info.buildSeconds;
            end
            cache.lookupSeconds = cache.lookupSeconds+info.lookupSeconds;
            cache.keys = string({cache.entries.key});
            info.hits = cache.hits;
            info.misses = cache.misses;
            obj.last.exactPackageReuseAudit = cache;
        end

        function schedOut = applyScheduleForceOffsetPolicy(obj, oldSched, schedIn)
        %APPLYSCHEDULEFORCEOFFSETPOLICY Runtime policy for scheduled forces_0.
        %
        % forces_0 is an affine preload/trim term.  During gust response, scheduling
        % it with instantaneous alpha injects a spurious generalized force jump.  The
        % default policy for gust validation is to hold the last accepted offset.
        
            schedOut = schedIn;
        
            if ~isfield(schedOut,'parConst') || ~isfield(schedOut.parConst,'forces_0')
                return
            end
        
            mode = char(obj.libOpt('forceOffsetMode', 'holdLast'));
        
            switch lower(mode)
        
                case {'scheduled','raw','none'}
                    % Use the interpolated library value as-is.
                    return
        
                case {'holdlast','hold','trimhold','frozen'}
                    if isfield(oldSched,'parConst') && isfield(oldSched.parConst,'forces_0') && ...
                            isequal(size(oldSched.parConst.forces_0), size(schedOut.parConst.forces_0))
        
                        schedOut.parConst.forces_0 = oldSched.parConst.forces_0;
                    end
        
                case {'blend','slow'}
                    if isfield(oldSched,'parConst') && isfield(oldSched.parConst,'forces_0') && ...
                            isequal(size(oldSched.parConst.forces_0), size(schedOut.parConst.forces_0))
        
                        beta = obj.libOpt('forceOffsetBlend', 0.02);
                        beta = min(max(beta,0),1);
        
                        schedOut.parConst.forces_0 = ...
                            (1-beta)*oldSched.parConst.forces_0 + ...
                             beta *schedOut.parConst.forces_0;
                    end
        
                otherwise
                    error('PlantRunTime:BadForceOffsetMode', ...
                        'Unknown cfg.library.forceOffsetMode = "%s".', mode);
            end
        
            if ~isfield(schedOut,'policy') || ~isstruct(schedOut.policy)
                schedOut.policy = struct();
            end
            schedOut.policy.forceOffsetMode = mode;
        end

        function alphaRefDeg = getScheduleAlphaReference(obj, rbTot)
        %GETSCHEDULEALPHAREFERENCE Reference alpha used by the ROM scheduler.
        %
        % rbTot may carry a trajectory/outer-loop scheduling command.  If no
        % reference is supplied, fall back to the trim schedule.  This keeps gust-only
        % tests from scheduling on the gust-induced alpha response.
        
            alphaRefDeg = obj.last.sched_mu_anchor(2);
        
            if ~isstruct(rbTot)
                return
            end
        
            if isfield(rbTot,'alphaScheduleDeg') && ~isempty(rbTot.alphaScheduleDeg)
                alphaRefDeg = rbTot.alphaScheduleDeg;
                return
            end
        
            if isfield(rbTot,'alpha_ref_deg') && ~isempty(rbTot.alpha_ref_deg)
                alphaRefDeg = rbTot.alpha_ref_deg;
                return
            end
        
            if isfield(rbTot,'alphaCmdDeg') && ~isempty(rbTot.alphaCmdDeg)
                alphaRefDeg = rbTot.alphaCmdDeg;
                return
            end
        
            if isfield(rbTot,'alphaScheduleRad') && ~isempty(rbTot.alphaScheduleRad)
                alphaRefDeg = rad2deg(rbTot.alphaScheduleRad);
                return
            end
        
            if isfield(rbTot,'alpha_ref_rad') && ~isempty(rbTot.alpha_ref_rad)
                alphaRefDeg = rad2deg(rbTot.alpha_ref_rad);
                return
            end
        end
    end

end

function cache = localExactPackageCache(last,maximumEntries)
%LOCALEXACTPACKAGECACHE Recover bounded run-local exact-query cache state.
if isfield(last,'exactPackageReuseAudit') && ...
        isstruct(last.exactPackageReuseAudit) && ...
        isfield(last.exactPackageReuseAudit,'entries')
    cache = last.exactPackageReuseAudit;
else
    cache = struct('entries',struct('key',{},'package',{}), ...
        'hits',0,'misses',0,'buildSeconds',0,'lookupSeconds',0);
end
required = {'entries','hits','misses','buildSeconds','lookupSeconds'};
if ~all(isfield(cache,required)) || numel(cache.entries)>maximumEntries
    cache = struct('entries',struct('key',{},'package',{}), ...
        'hits',0,'misses',0,'buildSeconds',0,'lookupSeconds',0);
end
end

function key = localExactPackageReuseKey(mu,ids,weights,sourceIds, ...
        candidate,registryHash)
%LOCALEXACTPACKAGEREUSEKEY Construct an exact binary schedule-package key.
queryHex = string(num2hex(double(mu(:))));
weightHex = string(num2hex(double(weights(:))));
idText = compose("%d",double(ids(:)));
fieldRoot = char(java.io.File(char(string(candidate.fieldRoot))).getCanonicalPath());
architecture = "full_coordinate_atomic_lift_interpolate_project";
key = strjoin(["query";queryHex;"ids";idText;"weights";weightHex; ...
    "sources";sourceIds(:);"architecture";architecture; ...
    "registry";string(registryHash);"fieldRoot";string(fieldRoot)],"|");
end

function digest = localExactPackageFileHash(path)
%LOCALEXACTPACKAGEFILEHASH Return the registry SHA-256 bound into a cache key.
path = char(string(path));
if ~isfile(path)
    digest = "";
    return
end
fileId = fopen(path,'rb');
if fileId < 0
    digest = "";
    return
end
cleanup = onCleanup(@()fclose(fileId));
data = fread(fileId,Inf,'*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(data,'int8'));
bytes = typecast(engine.digest(),'uint8');
digest = lower(string(reshape(dec2hex(bytes,2).',1,[])));
end
