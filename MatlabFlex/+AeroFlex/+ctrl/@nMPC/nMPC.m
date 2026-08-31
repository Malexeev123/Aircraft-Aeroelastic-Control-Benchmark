classdef nMPC < AeroFlex.ctrl.ControllerBase
%========================================================
%  NON-LINEAR MODEL-PREDICTIVE CONTROLLER (multiple-shoot)
%========================================================
%  Multiple-shooting NMPC transcription:
%
%      z = [X_0; X_1; ...; X_Nc; U_0; ...; U_{Nc-1}; s_T]
%
%  subject to
%
%      X_0 = xhat
%      X_{k+1} - Phi_k(X_k,U_k,W_k) = 0,   k = 0,...,Nc-1
%
%  with cost
%
%      J = sum_{k=0}^{Nc-1} 1/2 (X_k-xTrim)' Qc (X_k-xTrim)
%        + sum_{k=0}^{Nc-1} 1/2 (U_k-uTrim)' Rc (U_k-uTrim)
%        + 1/2 (X_Nc-xTrim)' Pc (X_Nc-xTrim)
%
%  Optional terminal set:
%
%      1/2 (X_Nc-xTrim)' Pc (X_Nc-xTrim) <= aT
%
%  An explicitly enabled terminal-viability mode replaces this inequality
%  by 1/2*(X_Nc-xTrim)'Pc*(X_Nc-xTrim)-aT-s_T <= 0, s_T >= 0,
%  with a configured quadratic penalty. It is default-disabled and a
%  nonzero s_T is reported as terminal-set non-satisfaction.
%
%  Notes:
%  * Disturbance is NOT an optimization variable here.
%    It is predicted from the MHE gust estimate and treated as exogenous.
%  * U contains trim-relative surface deflection and rate increments.
%    The locked total trim is added only at the prediction-model boundary.
%  * Only the first optimized control U_0 is applied.
%  * Warm start is shifted in receding-horizon fashion and then repaired
%    by forward prediction from the newest xhat.
%========================================================

    properties
        % ---------- timing & horizon ----------------------------------
        dt double
        Ts double
        Nc double

        % ---------- dimensions ----------------------------------------
        nx double
        nu double
        nw double = 1

        % ---------- weights / bounds ----------------------------------
        Qc double
        Rc double
        Pc double

        xL double
        xU double
        uL double
        uU double

        % Optional fields retained for compatibility. These are folded
        % into uL/uU by buildControlBounds() if var_per == 2.
        urL double
        urU double
        wL double
        wU double

        % ---------- internal model ------------------------------------
        model
        Sprev
        nativeStateCount (1,1) double
        reciprocalProviderEnabled logical = false
        reciprocalProvider
        reciprocalProviderFormulation (1,1) string = "disabled"
        reciprocalUnboundedAerodynamicLagStates logical = false
        reciprocalLatentIndex double
        reciprocalLatentTrim double
        reciprocalInitialLatentState double
        reciprocalPredictedLatentHorizon double
        reciprocalFutureContextHorizon cell
        identicalFutureContextCondensationEnabled logical = false
        identicalFutureContextPreparationCount (1,1) double = 0
        identicalFutureContextReuseCount (1,1) double = 0
        identicalFutureContextFallbackCount (1,1) double = 0
        identicalFutureContextLastAction (1,1) string = "disabled"

        % ---------- stored horizons / diagnostics ----------------------
        Xhist
        Uhist
        Whist
        uPrev
        k double = 0
        n_surf
        % ---------- optimizer state -----------------------------------
        z0 double
        H double
        solverOpts

        % ---------- disturbance prediction -----------------------------
        wHat double
        Nd double

        % ---------- terminal constraint --------------------------------
        aT double
        terminalViabilityEnabled logical = false
        terminalViabilityPenalty double = 1e6

        % ---------- trim/reference -------------------------------------
        xTrim double
        uTrim double
        uModelTrim double
        actuatorDeflectionAlpha double = 0.5

        % ---------- debug ----------------------------------------------
        debug logical = false
        gradientChecksEnabled logical = false
        dbg struct
        
        solverName string = "fmincon"
        sqpSolver
        prioritySqpSolver
        sqpCheckDone logical = false
        prioritySqpCheckDone logical = false

        % Default-disabled physical-output audit.  The production objective
        % remains exactly Qc/Rc/Pc unless the explicit audit field is enabled.
        wingtipOutputCostEnabled logical = false
        wingtipOutputCostStageWeight double = 0
        wingtipOutputCostTerminalWeight double = 0
        wingtipOutputGradient double
        wingtipOutputOwnerPolicy string = "legacy_casea_closed_checkpoint"

        % Default-inactive formal Case-A subspace restriction.  The general
        % two-surface controller remains independent so later lateral work
        % can retain differential authority.
        symmetricSurfaceSubspaceEnabled logical = false
        symmetricSurfaceSubspaceCaseId (1,1) string = "disabled"
        symmetricSurfaceSubspaceMap double = [1;1]
        symmetricSurfaceSubspaceChangeId (1,1) string = ""

        % Default-inactive exact rate-chart/full-condensing RTI seed. The
        % configured full NLP solver remains the correction and fallback.
        realtimeRtiEnabled logical = false
        realtimeRtiChangeId (1,1) string = ""
        realtimeRtiIterationCount (1,1) double = 2
        realtimeRtiChartMode (1,1) string = "disabled"
        realtimeRtiSolver
        realtimeRtiSolutionOwnerEnabled logical = false
        nativeReducedHorizonRtiRequested logical = false
        nativeReducedHorizonRtiActive logical = false
        nativeReducedHorizonRtiKernel function_handle = function_handle.empty
        nativeReducedHorizonRtiKernelIdentity (1,1) string = "disabled"
        nativeReducedHorizonRtiLastFallback (1,1) string = ""
        nativeValueHorizonRequested logical = false
        nativeValueHorizonActive logical = false
        nativeValueHorizonKernel function_handle = function_handle.empty
        nativeValueHorizonKernelIdentity (1,1) string = "disabled"
        nativeValueHorizonLastFallback (1,1) string = ""
        nativeCausalRolloutRequested logical = false
        nativeCausalRolloutActive logical = false
        nativeCausalRolloutKernel function_handle = function_handle.empty
        nativeCausalRolloutKernelIdentity (1,1) string = "disabled"
        nativeCausalRolloutLastFallback (1,1) string = ""
        preparedHorizonDataReuseRequested logical = false
        preparedHorizonDataReuseActive logical = false
        acceptedReplayReuseRequested logical = false
        acceptedReplayReuseActive logical = false
        acceptedLatentHorizonCondensationRequested logical = false
        acceptedLatentHorizonCondensationActive logical = false

        % Default-inactive scheduled source-domain replay gate for the
        % approved Case-B reciprocal RTI audit.
        scheduledSourceDomainRtiConstraintEnabled logical = false
        scheduledSourceDomainRtiConstraintChangeId (1,1) string = ""

    end

    %==================================================================
    methods
    %------------------------------------------------------------------
    function obj = nMPC(cfg,beam,aero,base,trim)
    % Constructor
    %------------------------------------------------------------------
        obj@AeroFlex.ctrl.ControllerBase(cfg,trim);

        % ---------- timing --------------------------------------------
        obj.Ts = cfg.ctrl.Ts;
        obj.Nc = cfg.ctrl.Nc;

        % ---------- model and dimensions ------------------------------
        obj.model = cfg.modelHandle(cfg,beam,aero,base);
        obj.dt    = obj.model.dt;
        obj.nativeStateCount = size(obj.model.L,1);
        obj.nx    = obj.nativeStateCount;
        obj.nu    = cfg.ctrl.n_surf * cfg.ctrl.var_per;
        obj.n_surf = cfg.ctrl.n_surf;
        if isfield(cfg,'nw')
            obj.nw = cfg.nw;
        else
            obj.nw = 1;
        end

        assert(cfg.ctrl.var_per == 2 && obj.nu == 2*obj.n_surf, ...
            'nMPC:ControlOrdering', ...
            'nMPC requires [deflection; rate] channels for each surface set.');
        obj.configureSymmetricSurfaceSubspace(cfg);
        nSubsteps = round(obj.Ts/obj.dt);
        timeTolerance = 100*eps(max([1,obj.Ts,obj.dt]));
        assert(nSubsteps >= 1 && ...
            abs(nSubsteps*obj.dt-obj.Ts) <= timeTolerance, ...
            'nMPC:SampleAlignment', ...
            'Controller Ts must be an integer multiple of prediction-model dt.');

        % ---------- weights -------------------------------------------
        obj.Qc = cfg.Qc;
        obj.Rc = cfg.Rc;
        obj.Pc = cfg.Pc;

        % ---------- trim/reference ------------------------------------
        obj.xTrim = trim.states(:);
        obj.uModelTrim = obj.buildModelTrim(trim);
        obj.uTrim = zeros(obj.nu,1);

        assert(numel(obj.xTrim)==obj.nx && all(isfinite(obj.xTrim)), ...
            'nMPC:TrimState', ...
            'trim.states must contain %d finite prediction states.',obj.nx);
        obj.configureReciprocalProvider(cfg);

        if isfield(cfg.ctrl,'actuatorDeflectionAlpha')
            obj.actuatorDeflectionAlpha = cfg.ctrl.actuatorDeflectionAlpha;
        end
        assert(isscalar(obj.actuatorDeflectionAlpha) && ...
            isfinite(obj.actuatorDeflectionAlpha) && ...
            obj.actuatorDeflectionAlpha >= 0 && ...
            obj.actuatorDeflectionAlpha <= 1, ...
            'nMPC:ActuatorDeflectionAlpha', ...
            'cfg.ctrl.actuatorDeflectionAlpha must be in [0,1].');

        if isfield(cfg,'debug') && isstruct(cfg.debug) && ...
                isfield(cfg.debug,'level') && cfg.debug.level >= 3
            obj.debug = true;
            if isfield(cfg.ctrl,'sqp') && ...
                    isfield(cfg.ctrl.sqp,'CheckGradientsOnce')
                obj.gradientChecksEnabled = ...
                    logical(cfg.ctrl.sqp.CheckGradientsOnce);
            end
        end

        % ---------- state bounds --------------------------------------
        % Same style as your nMHE. These are broad safety/solver bounds.
        margin  = 10;
        absFloor = 1;

        delta = max(abs(margin .* obj.xTrim), absFloor);
        xLcfg = obj.xTrim - delta;
        xUcfg = obj.xTrim + delta;

        obj.xL = obj.expandToLength(xLcfg,obj.nx,'xL');
        obj.xU = obj.expandToLength(xUcfg,obj.nx,'xU');
        if obj.reciprocalUnboundedAerodynamicLagStates
            aerodynamicLagIndex = obj.model.idx.qGam(:);
            assert(all(aerodynamicLagIndex>=1 & ...
                aerodynamicLagIndex<=obj.nativeStateCount), ...
                'nMPC:ReciprocalAerodynamicLagIndex', ...
                'The reciprocal aerodynamic-lag state index is invalid.');
            obj.xL(aerodynamicLagIndex) = -inf;
            obj.xU(aerodynamicLagIndex) = inf;
        end

        % ---------- control bounds ------------------------------------
        [obj.uL,obj.uU,obj.urL,obj.urU] = obj.buildControlBounds(cfg);

        % Optional disturbance bounds retained for compatibility.
        if isfield(cfg,'wL')
            obj.wL = cfg.wL;
        else
            obj.wL = -inf(obj.nw,1);
        end

        if isfield(cfg,'wU')
            obj.wU = cfg.wU;
        else
            obj.wU = inf(obj.nw,1);
        end

        % ---------- terminal radius -----------------------------------
        if isfield(cfg.ctrl,'aT')
            obj.aT = cfg.ctrl.aT;
        else
            obj.aT = inf;   % disables terminal inequality
        end

        outputCost = struct();
        hasProfileOutputCost = isfield(cfg.ctrl,'nmpcWingtipOutput') && ...
            ~isempty(cfg.ctrl.nmpcWingtipOutput);
        hasAuditOutputCost = isfield(cfg.ctrl,'nmpcWingtipOutputAudit') && ...
            ~isempty(cfg.ctrl.nmpcWingtipOutputAudit);
        assert(~(hasProfileOutputCost && hasAuditOutputCost), ...
            'nMPC:WingtipOutputConfigConflict', ...
            ['Select either nmpcWingtipOutput or the retained audit alias, ', ...
             'not both.']);
        if hasProfileOutputCost
            outputCost = cfg.ctrl.nmpcWingtipOutput;
        elseif hasAuditOutputCost
            outputCost = cfg.ctrl.nmpcWingtipOutputAudit;
        end
        if ~isempty(fieldnames(outputCost))
            assert(isstruct(outputCost) && isfield(outputCost,'enabled') && ...
                isfield(outputCost,'stageWeight') && ...
                isfield(outputCost,'terminalWeight'), ...
                'nMPC:WingtipOutputConfig', ...
                ['The wingtip output configuration requires enabled, ', ...
                 'stageWeight, and terminalWeight.']);
            obj.wingtipOutputCostEnabled = logical(outputCost.enabled);
            obj.wingtipOutputCostStageWeight = outputCost.stageWeight;
            obj.wingtipOutputCostTerminalWeight = outputCost.terminalWeight;
            if isfield(outputCost,'ownerPolicy')
                obj.wingtipOutputOwnerPolicy = string(outputCost.ownerPolicy);
            end
        end
        assert(isscalar(obj.wingtipOutputCostEnabled) && ...
            isscalar(obj.wingtipOutputCostStageWeight) && ...
            isscalar(obj.wingtipOutputCostTerminalWeight) && ...
            isfinite(obj.wingtipOutputCostStageWeight) && ...
            isfinite(obj.wingtipOutputCostTerminalWeight) && ...
            obj.wingtipOutputCostStageWeight >= 0 && ...
            obj.wingtipOutputCostTerminalWeight >= 0, ...
            'nMPC:WingtipOutputConfig', ...
            'Wingtip output weights must be finite nonnegative scalars.');
        if obj.wingtipOutputCostEnabled
            if isfield(outputCost,'gradient')
                gradient = double(outputCost.gradient(:).');
                assert(numel(gradient) == obj.nx && ...
                    all(isfinite(gradient)), ...
                    'nMPC:WingtipOutputConfig', ...
                    'The package-owned wingtip gradient is invalid.');
                obj.wingtipOutputGradient = gradient;
            else
                obj.wingtipOutputGradient = ...
                    obj.buildWingtipOutputGradient(beam,base);
            end
        else
            obj.wingtipOutputGradient = zeros(1,obj.nx);
        end
        if isfield(cfg.ctrl,'terminalViability') && ...
                ~isempty(cfg.ctrl.terminalViability)
            terminalViability = cfg.ctrl.terminalViability;
            assert(isstruct(terminalViability) && ...
                isfield(terminalViability,'enabled') && ...
                isfield(terminalViability,'penalty'), ...
                'nMPC:TerminalViabilityConfig', ...
                ['cfg.ctrl.terminalViability requires logical enabled and ', ...
                 'positive finite penalty fields.']);
            obj.terminalViabilityEnabled = logical(terminalViability.enabled);
            obj.terminalViabilityPenalty = terminalViability.penalty;
        end
        assert(isscalar(obj.terminalViabilityEnabled) && ...
            isscalar(obj.terminalViabilityPenalty) && ...
            isfinite(obj.terminalViabilityPenalty) && ...
            obj.terminalViabilityPenalty > 0, ...
            'nMPC:TerminalViabilityConfig', ...
            'Terminal-viability penalty must be positive and finite.');
        assert(~obj.terminalViabilityEnabled || isfinite(obj.aT), ...
            'nMPC:TerminalViabilityConfig', ...
            'Terminal viability requires a finite terminal-set radius aT.');

        % ---------- gust prediction length -----------------------------
        obj.Nd   = ceil(obj.Nc/2);
        obj.wHat = zeros(obj.nw,1);

        % ---------- histories / stored horizons ------------------------
        obj.Xhist = repmat(obj.xTrim,1,obj.Nc+1);
        obj.Uhist = repmat(obj.uTrim,1,obj.Nc);
        obj.Whist = zeros(obj.nw,obj.Nc);

        obj.uPrev = obj.uTrim;

        % ---------- decision-vector initial guess ----------------------
        idx = obj.buildIndexMaps();
        nVar = obj.decisionVariableCount();

        obj.z0 = zeros(nVar,1);

        for j = 1:obj.Nc+1
            obj.z0(idx.x{j}) = obj.xTrim;
        end

        for j = 1:obj.Nc
            obj.z0(idx.u{j}) = obj.uTrim;
        end

        obj.H = speye(nVar);

        % Initial STM placeholder. Not required by this direct multiple-
        % shooting transcription, but retained for interface compatibility.
        obj.Sprev = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];

        % ---------- solver options ------------------------------------
        obj.solverOpts = optimoptions('fmincon', ...
            'Algorithm','sqp', ...      % interior-point | sqp
            'SpecifyObjectiveGradient',true, ...
            'SpecifyConstraintGradient',true, ...
            'Display','none', ...
            'OptimalityTolerance',1e-4, ...
            'StepTolerance',1e-7, ...
            'ConstraintTolerance',1e-6, ...
            'MaxIterations',150);

        % Useful while debugging:
        % obj.solverOpts.FiniteDifferenceStepSize = 1e-8;

        % obj.solverOpts.CheckGradients = true;
        % obj.solverOpts.FiniteDifferenceType = 'central';

        obj.dbg = struct('t',[],'U',[],'cont',[]);

        %======================================================================
        % Solver selection: fmincon or custom SQP
        %======================================================================
        if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'mpcSolver')
            obj.solverName = lower(string(cfg.ctrl.mpcSolver));
        else
            obj.solverName = "fmincon";
        end
        
        switch obj.solverName
        
            case "fmincon"
                obj.sqpSolver = [];
        
            case "custom_sqp"
                sqpOpts = AeroFlex.optim.SQPSolver.defaultOptions();
        
                if isfield(cfg.ctrl,'sqp') && ~isempty(cfg.ctrl.sqp)
                    f = fieldnames(cfg.ctrl.sqp);
        
                    for ii = 1:numel(f)
                        sqpOpts.(f{ii}) = cfg.ctrl.sqp.(f{ii});
                    end
                end
        
                obj.sqpSolver = AeroFlex.optim.SQPSolver(sqpOpts);
                % The lexicographic objective has independent curvature.
                % Do not transfer the primary BFGS approximation into it.
                obj.prioritySqpSolver = AeroFlex.optim.SQPSolver(sqpOpts);
        
            otherwise
                error('nMPC:Solver', ...
                      'Unknown cfg.ctrl.mpcSolver = "%s". Use "fmincon" or "custom_sqp".', ...
                      obj.solverName);
        end
        
        obj.sqpCheckDone = ~obj.gradientChecksEnabled;
        obj.prioritySqpCheckDone = ~obj.gradientChecksEnabled;
        obj.configureRealtimeRtiAudit(cfg);
        obj.configureNativeReducedHorizonRtiAudit(cfg);
        obj.configureHighLeverageRuntimeAudit(cfg);
        obj.configureIdenticalFutureContextCondensationAudit(cfg);
        obj.configureScheduledSourceDomainRtiConstraintAudit(cfg);

        % The configured model provider owns the scheduled rate-projection
        % policy. Do not overwrite it inside the controller.

    end

    %------------------------------------------------------------------
    function synchronizeTimingFromModel(obj)
    %SYNCHRONIZETIMINGFROMMODEL Adopt the active scheduled-model sample step.
        assert(isprop(obj.model,'dt') && isscalar(obj.model.dt) && ...
            isfinite(obj.model.dt) && obj.model.dt > 0, ...
            'nMPC:ModelStep', ...
            'The active controller prediction model must define a positive scalar dt.');
        obj.dt = obj.model.dt;
        nSubsteps = round(obj.Ts/obj.dt);
        timeTolerance = 100*eps(max([1,obj.Ts,obj.dt]));
        assert(nSubsteps >= 1 && ...
            abs(nSubsteps*obj.dt-obj.Ts) <= timeTolerance, ...
            'nMPC:SampleAlignment', ...
            'Controller Ts must be an integer multiple of the active prediction-model dt.');
    end

    %------------------------------------------------------------------
    function info = transportScheduledState( ...
            obj,transform,newReference,~,~,outputContract)
    %TRANSPORTSCHEDULEDSTATE Move persistent controller state to a new chart.
        transform = full(transform);
        newReference = newReference(:);
        assert(isequal(size(transform),[obj.nx,obj.nx]) && ...
            all(isfinite(transform),'all') && numel(newReference) == obj.nx && ...
            all(isfinite(newReference)), ...
            'nMPC:ScheduledStateTransport', ...
            'The scheduled controller state map is invalid.');
        inverse = transform\eye(obj.nx);
        index = obj.buildIndexMaps();
        newXhist = transform*obj.Xhist;
        newZ0 = obj.z0;
        for node = 1:obj.Nc+1
            newZ0(index.x{node}) = transform*obj.z0(index.x{node});
        end
        newQc = inverse.'*obj.Qc*inverse;
        newPc = inverse.'*obj.Pc*inverse;
        newQc = 0.5*(newQc+newQc.');
        newPc = 0.5*(newPc+newPc.');
        delta = max(abs(10*newReference),1);
        obj.Xhist = newXhist;
        obj.z0 = newZ0;
        obj.Qc = newQc;
        obj.Pc = newPc;
        obj.xL = newReference-delta;
        obj.xU = newReference+delta;
        if obj.reciprocalUnboundedAerodynamicLagStates
            qGamIndex = obj.model.idx.qGam(:);
            obj.xL(qGamIndex) = -inf;
            obj.xU(qGamIndex) = inf;
        end
        obj.xTrim = newReference;
        if obj.wingtipOutputCostEnabled
            if nargin >= 6 && ~isempty(outputContract)
                assert(isstruct(outputContract) && ...
                    isfield(outputContract,'symmetricGradient') && ...
                    isfield(outputContract,'ownerPolicy'), ...
                    'nMPC:ScheduledWingtipOutput', ...
                    'The scheduled physical-output contract is invalid.');
                gradient = double(outputContract.symmetricGradient(:).');
                assert(numel(gradient) == obj.nx && ...
                    all(isfinite(gradient)), ...
                    'nMPC:ScheduledWingtipOutput', ...
                    'The scheduled physical-output gradient is invalid.');
                obj.wingtipOutputGradient = gradient;
                obj.wingtipOutputOwnerPolicy = ...
                    string(outputContract.ownerPolicy);
            else
                assert(obj.wingtipOutputOwnerPolicy ~= ...
                    "package_owned_signed_mirrored_tip", ...
                    'nMPC:ScheduledWingtipOutput', ...
                    ['A package-owned wingtip objective requires an atomic ', ...
                     'output-contract rebuild at every schedule change.']);
                % Preserve the closed legacy checkpoint behavior.
                obj.wingtipOutputGradient = ...
                    obj.wingtipOutputGradient*inverse;
            end
        end
        obj.H = speye(numel(obj.z0));
        obj.Sprev = [eye(obj.nx),zeros(obj.nx,obj.nw+obj.nu)];
        info = struct('condition',cond(transform), ...
            'stageWeightSymmetryError',norm(obj.Qc-obj.Qc.','fro'), ...
            'terminalWeightSymmetryError',norm(obj.Pc-obj.Pc.','fro'), ...
            'accepted',true);
    end

    %------------------------------------------------------------------
    function [uk,Uinfo] = computeControl(obj,xhat,whatEst,uPrev,t_k,uBaseHorizon,priorityReference)
    % Real-time NMPC call.
    %
    % Inputs:
    %   xhat    : current MHE state estimate, nx x 1
    %   whatEst : current or horizon disturbance estimate from nMHE
    %             If this is a full MHE horizon vector, the last nw entries
    %             are used as the current disturbance estimate.
    %   uPrev   : actuator command actually applied at previous sample
    %   t_k     : current time, for debug plots
    %   uBaseHorizon : optional known applied-command base, nu x Nc. When
    %                  supplied, the decision is the applied command and
    %                  the existing input penalty is evaluated on its
    %                  residual from this base. This audit-only interface
    %                  leaves the prediction, bounds, rate geometry, and
    %                  final nonlinear acceptance unchanged.
    %   priorityReference : optional audit-only lexicographic wing-reference
    %                  selection. The existing nMPC problem remains primary;
    %                  the reference is considered only by a second solve
    %                  constrained to a primary-cost equivalence guard.
    %
    % Output:
    %   uk      : first optimized control move
    %   Uinfo   : diagnostics
    %------------------------------------------------------------------
        obj.synchronizeTimingFromModel();
        xhat = xhat(:);
        assert(numel(xhat)==obj.nx, ...
            'nMPC:Dimension','xhat must have length nx.');

        % Use the actual trim-relative plant command if supplied.
        if nargin >= 4 && ~isempty(uPrev)
            obj.uPrev = obj.expandToLength(uPrev,obj.nu,'uPrev');
        end
        if obj.symmetricSurfaceSubspaceEnabled
            priorResidual = obj.surfaceSubspaceResidual(obj.uPrev);
            assert(norm(priorResidual,inf) <= 1e-12, ...
                'nMPC:SymmetricSurfacePriorEndpoint', ...
                ['The formal Case-A symmetric-surface policy cannot start ', ...
                 'from an asymmetric prior endpoint (residual %.3e).'], ...
                norm(priorResidual,inf));
        end

        if nargin < 6 || isempty(uBaseHorizon)
            uBaseHorizon = repmat(obj.uTrim,1,obj.Nc);
        else
            assert(isnumeric(uBaseHorizon) && isequal( ...
                size(uBaseHorizon),[obj.nu,obj.Nc]) && ...
                all(isfinite(uBaseHorizon),'all'), ...
                'nMPC:KnownBaseHorizon', ...
                'uBaseHorizon must be a finite nu-by-Nc command matrix.');
        end
        if nargin < 7 || isempty(priorityReference)
            priorityReference = struct('enabled',false);
        end
        priority = obj.parsePriorityReference(priorityReference);
        primaryBaseHorizon = uBaseHorizon;
        if priority.enabled
            % A priority reference is not an independently applied command.
            % Retain the unmodified trim-relative input objective in level one.
            primaryBaseHorizon = repmat(obj.uTrim,1,obj.Nc);
        end

        % Parse the newest gust estimate from MHE.
        obj.wHat = obj.parseDisturbanceEstimate(whatEst);

        % Store latest measured/estimated quantities for diagnostics.
        obj.k = obj.k + 1;

        obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
        obj.Uhist(:,1:end-1) = obj.Uhist(:,2:end);
        obj.Whist(:,1:end-1) = obj.Whist(:,2:end);

        obj.Xhist(:,end) = xhat;
        obj.Uhist(:,end) = obj.uPrev;
        obj.Whist(:,end) = obj.wHat;

        % Build exogenous future disturbance profile.
        wHorz = obj.buildPredictedGust();

        % Repair the shifted warm start so that:
          % X_0 = xhat
        % and the guessed state horizon is dynamically consistent under
        % the current shifted control guess and predicted disturbance.

        zControlSeed = obj.z0;
        if obj.realtimeRtiEnabled
            zControlSeed = obj.repairRealtimeRtiControlGuess(zControlSeed);
        end
        sharedValueHorizonData = struct();
        if obj.preparedHorizonDataReuseActive
            sharedValueHorizonData = ...
                obj.buildNativeControllerValueHorizonData();
        end
        zSolve0 = obj.repairStateGuess( ...
            zControlSeed,xhat,wHorz,sharedValueHorizonData);
        % Assemble and solve the full horizon NLP.
        %======================================================================
        % Assemble and solve NMPC multiple-shooting NLP
        %======================================================================
        nlp = obj.assembleWindow( ...
            xhat,wHorz,primaryBaseHorizon,sharedValueHorizonData);
        optimizerInitial = zSolve0;
        rtiCandidate = zeros(0,1);
        rtiInfo = obj.emptyRealtimeRtiInfo();
        if obj.realtimeRtiEnabled
            try
                [rtiCandidate,rtiInfo] = obj.buildRealtimeRtiSeed( ...
                    nlp,zSolve0,wHorz);
                rtiInfo = obj.qualifyScheduledSourceDomainRtiCandidate( ...
                    rtiInfo,rtiCandidate,xhat,whatEst);
                if rtiInfo.qualified
                    optimizerInitial = rtiCandidate;
                else
                    rtiInfo.fallbackToFullInitial = true;
                end
            catch rtiException
                rtiInfo = obj.emptyRealtimeRtiInfo();
                rtiInfo.attempted = true;
                rtiInfo.fallbackToFullInitial = true;
                rtiInfo.message = "RTI preparation failed closed: " + ...
                    string(rtiException.message);
                rtiInfo.identifier = string(rtiException.identifier);
            end
        end
        if obj.gradientChecksEnabled && ~obj.sqpCheckDone && ...
                obj.solverName == "custom_sqp"
            obj.localCheckNMPCEqualityGradient(nlp,zSolve0,nlp.lb,nlp.ub);
            obj.localCheckNMPCEqualityGradientBlocks( ...
                nlp,zSolve0,nlp.lb,nlp.ub);
            obj.sqpSolver.checkGradients( ...
                nlp.cost,nlp.nonl,zSolve0,nlp.lb,nlp.ub);
            obj.sqpCheckDone = true;
        end
        if obj.realtimeRtiSolutionOwnerEnabled
            if rtiInfo.qualified
                zOpt = rtiCandidate;
                fval = rtiInfo.objective;
                exitflag = 1;
                output = struct( ...
                    'message',"Qualified condensed RTI owned the action.", ...
                    'iterations',rtiInfo.completedIterations, ...
                    'constrviolation',rtiInfo.constraintViolationInf, ...
                    'firstorderopt',obj.realtimeRtiFirstOrderMetric(rtiInfo));
                lambda = obj.emptySolverMultipliers();
            else
                zOpt = rtiCandidate;
                fval = inf;
                exitflag = -2;
                output = struct('message', ...
                    "Condensed RTI was rejected; holding the prior endpoint.");
                lambda = obj.emptySolverMultipliers();
            end
        else
        switch obj.solverName
                
            case "fmincon"
        
                [zOpt,fval,exitflag,output,lambda] = fmincon( ...
                    nlp.cost,optimizerInitial, ...
                    [],[],[],[], ...
                    nlp.lb,nlp.ub, ...
                    nlp.nonl,obj.solverOpts);
        
            case "custom_sqp"
        
                [zOpt,fval,exitflag,output,lambda] = obj.sqpSolver.solve( ...
                    nlp.cost,optimizerInitial,nlp.lb,nlp.ub,nlp.nonl);
        
            otherwise
                error('nMPC:Solver','Unhandled solverName = "%s".', obj.solverName);
        end
        end

        priorityInfo = struct('enabled',priority.enabled, ...
            'attempted',false,'accepted',false,'fallbackToPrimary',false, ...
            'primaryCost',nan,'primaryCostGuard',nan,'secondaryCost',nan, ...
            'selectedPrimaryCost',nan,'guardViolation',nan, ...
            'guardSatisfied',false,'exitflag',nan, ...
            'referenceHorizon',priority.referenceHorizon);
        if priority.enabled && obj.realtimeRtiSolutionOwnerEnabled
            % The online owner solves the unchanged primary safety problem.
            % Do not silently invoke the full-NLP secondary solve.
            priorityInfo.fallbackToPrimary = true;
        elseif priority.enabled && exitflag > 0 && ~isempty(zOpt) && ...
                numel(zOpt) == numel(zSolve0) && all(isfinite(zOpt))
            [primaryCost,~] = nlp.cost(zOpt);
            primaryCostGuard = primaryCost + priority.primaryCostTolerance;
            priorityNlp = obj.buildPriorityNlp( ...
                nlp,priority.referenceHorizon,primaryCostGuard);
            priorityInfo.attempted = true;
            priorityInfo.primaryCost = primaryCost;
            priorityInfo.primaryCostGuard = primaryCostGuard;
            referenceIsNontrivial = norm(priority.referenceHorizon( ...
                1:obj.n_surf,:) - obj.uTrim(1:obj.n_surf),inf) > 0;
            if obj.gradientChecksEnabled && referenceIsNontrivial && ...
                    ~obj.prioritySqpCheckDone && ...
                    obj.solverName == "custom_sqp"
                obj.prioritySqpSolver.checkGradients( ...
                    priorityNlp.cost,priorityNlp.nonl,zOpt, ...
                    priorityNlp.lb,priorityNlp.ub);
                obj.prioritySqpCheckDone = true;
            end
            switch obj.solverName
                case "fmincon"
                    [zPriority,fSecondary,flagPriority,outputPriority,lambdaPriority] = ...
                        fmincon(priorityNlp.cost,zOpt,[],[],[],[], ...
                        priorityNlp.lb,priorityNlp.ub, ...
                        priorityNlp.nonl,obj.solverOpts);
                case "custom_sqp"
                    [zPriority,fSecondary,flagPriority,outputPriority,lambdaPriority] = ...
                        obj.prioritySqpSolver.solve(priorityNlp.cost,zOpt, ...
                        priorityNlp.lb,priorityNlp.ub,priorityNlp.nonl);
            end
            priorityInfo.exitflag = flagPriority;
            priorityInfo.secondaryCost = fSecondary;
            if flagPriority > 0 && ~isempty(zPriority) && ...
                    numel(zPriority) == numel(zSolve0) && all(isfinite(zPriority))
                [selectedPrimaryCost,~] = nlp.cost(zPriority);
                priorityInfo.selectedPrimaryCost = selectedPrimaryCost;
                priorityInfo.guardViolation = selectedPrimaryCost - primaryCostGuard;
                priorityInfo.guardSatisfied = selectedPrimaryCost <= ...
                    primaryCostGuard + 100*eps(max(1,abs(primaryCostGuard)));
                if priorityInfo.guardSatisfied
                    zOpt = zPriority;
                    fval = selectedPrimaryCost;
                    exitflag = flagPriority;
                    output = outputPriority;
                    lambda = lambdaPriority;
                    priorityInfo.accepted = true;
                else
                    % Solver feasibility tolerance does not authorize a
                    % primary-cost relaxation. Reject this secondary result
                    % and retain the accepted unchanged primary command.
                    priorityInfo.fallbackToPrimary = true;
                end
            else
                priorityInfo.fallbackToPrimary = true;
            end
        end

        idx = obj.buildIndexMaps();

        % Diagnostics even if the solve fails.
        if isempty(zOpt) || numel(zOpt) ~= numel(zSolve0) || ...
                any(~isfinite(zOpt))
            continuityNorm = inf;
            constraintViolationInf = inf;
        else
            [~,ceq] = nlp.nonl(zOpt);
            continuityNorm = norm(ceq);
            if isfield(output,'constrviolation') && ...
                    isscalar(output.constrviolation) && ...
                    isfinite(output.constrviolation)
                constraintViolationInf = output.constrviolation;
            else
                [c,ceq] = nlp.nonl(zOpt);
                constraintViolationInf = max([0;c(:);abs(ceq(:)); ...
                    nlp.lb(:)-zOpt(:);zOpt(:)-nlp.ub(:)]);
            end
        end

        if exitflag <= 0
            warning('nMPC:SolveFailure', ...
            'nMPC solver "%s" stopped with exitflag %d. Holding previous input. Message: %s', ...
            obj.solverName, exitflag, output.message);

            if obj.scheduledSourceDomainRtiConstraintEnabled && ...
                    isfield(rtiInfo,'sourceDomainRejected') && ...
                    logical(rtiInfo.sourceDomainRejected)
                uk = obj.uTrim;
                Ufallback = repmat(obj.uTrim,1,obj.Nc);
                obj.uPrev = uk;
                obj.z0 = obj.resetRealtimeRtiControlHorizon(obj.z0);
            else
                uk = obj.uPrev;
                Ufallback = repmat(obj.uPrev,1,obj.Nc);
            end

            Uinfo.cost       = fval;
            Uinfo.exitflag   = exitflag;
            Uinfo.uHorizon   = Ufallback;
            Uinfo.candidateUHorizon = nan(obj.nu,obj.Nc);
            Uinfo.candidateFirstCommand = nan(obj.nu,1);
            Uinfo.wHorizon   = wHorz;
            Uinfo.uBaseHorizon = primaryBaseHorizon;
            Uinfo.priority = priorityInfo;
            Uinfo.continuity = continuityNorm;
            Uinfo.continuityTwoNorm = continuityNorm;
            Uinfo.constraintViolationInf = constraintViolationInf;
            Uinfo.output     = output;
            Uinfo.lambda     = lambda;
            Uinfo.terminalMode = "hard";
            Uinfo.terminalSlack = nan;
            Uinfo.terminalSetValue = nan;
            Uinfo.terminalSetSatisfied = false;
            Uinfo.surfaceSubspace = obj.surfaceSubspaceDiagnostics(Ufallback);
            Uinfo.realtimeRti = rtiInfo;
            if ~isempty(zOpt) && numel(zOpt) == numel(zSolve0) && ...
                    all(isfinite(zOpt))
                candidateUHorizon = reshape( ...
                    zOpt(idx.u{1}(1):idx.u{end}(end)),obj.nu,obj.Nc);
                Uinfo.candidateUHorizon = candidateUHorizon;
                Uinfo.candidateFirstCommand = candidateUHorizon(:,1);
                terminalInfo = obj.evaluateTerminal(zOpt);
                Uinfo.terminalMode = terminalInfo.mode;
                Uinfo.terminalSlack = terminalInfo.slack;
                Uinfo.terminalSetValue = terminalInfo.value;
                Uinfo.terminalSetSatisfied = terminalInfo.satisfied;
            end

            if obj.debug
                obj.debugPlots(t_k,Uinfo);
            end

            return
        end

        % Extract solution.
        Xopt = reshape(zOpt(1:(obj.Nc+1)*obj.nx),obj.nx,obj.Nc+1);
        Uopt = reshape(zOpt(idx.u{1}(1):idx.u{end}(end)),obj.nu,obj.Nc);
        terminalInfo = obj.evaluateTerminal(zOpt);
        if obj.reciprocalProviderEnabled
            acceptedReplayEndpoints = zeros(0,0);
            acceptedReplayInfo = struct('enabled', ...
                obj.acceptedReplayReuseActive,'cacheHit',false, ...
                'valueReplayCacheHits',0,'valueReplayCacheMisses',0);
            if obj.acceptedReplayReuseActive
                [acceptedReplayEndpoints,acceptedReplayInfo] = ...
                    nlp.getAcceptedReplay(zOpt);
                assert(acceptedReplayInfo.cacheHit, ...
                    'nMPC:AcceptedReplayReuseMiss', ...
                    ['The accepted controller decision was not the exact ', ...
                     'decision from the final value replay.']);
            end
            obj.reciprocalPredictedLatentHorizon = ...
                obj.rebuildReciprocalLatentHorizon( ...
                    Xopt,Uopt,wHorz,sharedValueHorizonData, ...
                    acceptedReplayEndpoints);
            rtiInfo.acceptedLatentHorizonCondensationApplied = ...
                obj.acceptedLatentHorizonCondensationActive;
            rtiInfo.acceptedReplayReuseApplied = ...
                acceptedReplayInfo.cacheHit;
            rtiInfo.valueReplayCacheHits = ...
                acceptedReplayInfo.valueReplayCacheHits;
            rtiInfo.valueReplayCacheMisses = ...
                acceptedReplayInfo.valueReplayCacheMisses;
        end

        % Apply only the first control move.
        % uk = zOpt(idx.u{2});
        uk = zOpt(idx.u{1});

        % Store optimized prediction horizons.
        obj.Xhist = Xopt;
        obj.Uhist = Uopt;
        obj.Whist = wHorz;

        obj.uPrev = uk;

        % Receding horizon warm start for next sample.
        obj.z0 = obj.shiftGuess(zOpt);

        Uinfo.cost       = fval;
        Uinfo.exitflag   = exitflag;
        Uinfo.output     = output;
        Uinfo.lambda     = lambda;
        Uinfo.uHorizon   = Uopt;
        Uinfo.candidateUHorizon = Uopt;
        Uinfo.candidateFirstCommand = uk;
        Uinfo.uBaseHorizon = primaryBaseHorizon;
        Uinfo.priority = priorityInfo;
        Uinfo.wHorizon   = wHorz;
        Uinfo.continuity = continuityNorm;
        Uinfo.continuityTwoNorm = continuityNorm;
        Uinfo.constraintViolationInf = constraintViolationInf;
        Uinfo.terminalMode = terminalInfo.mode;
        Uinfo.terminalSlack = terminalInfo.slack;
        Uinfo.terminalSetValue = terminalInfo.value;
        Uinfo.terminalSetSatisfied = terminalInfo.satisfied;
        Uinfo.surfaceSubspace = obj.surfaceSubspaceDiagnostics(Uopt);
        Uinfo.realtimeRti = rtiInfo;
        
        if isfield(output,'constrviolation')
            Uinfo.constrviolation = output.constrviolation;
        end
        
        if isfield(output,'firstorderopt')
            Uinfo.firstorderopt = output.firstorderopt;
        end
        
        if isfield(output,'stepsize')
            Uinfo.stepsize = output.stepsize;
        end
        
        if isfield(output,'qpExitflag')
            Uinfo.qpExitflag = output.qpExitflag;
        end
        
        if isfield(output,'slackEqInf')
            Uinfo.slackEqInf = output.slackEqInf;
        end
        
        if isfield(output,'slackIneqInf')
            Uinfo.slackIneqInf = output.slackIneqInf;
        end

        if obj.debug
            obj.debugPlots(t_k,Uinfo);
        end
    end

    %------------------------------------------------------------------
    function [nlp,zSolve0,wHorz,idx] = assembleAuditWindow( ...
            obj,xhat,whatEst,auditOnly)
    %ASSEMBLEAUDITWINDOW Return the current NLP without changing runtime state.
        assert(isscalar(auditOnly) && islogical(auditOnly) && auditOnly, ...
            'nMPC:AuditWindowGuard', ...
            'assembleAuditWindow requires explicit logical auditOnly=true.');
        xhat = xhat(:);
        assert(numel(xhat) == obj.nx && all(isfinite(xhat)), ...
            'nMPC:AuditWindowState', ...
            'The audit state must contain nx finite entries.');
        w0 = obj.parseDisturbanceEstimate(whatEst);
        wHorz = zeros(obj.nw,obj.Nc);
        if norm(w0) > 0
            count = max(1,min(obj.Nd,obj.Nc));
            wHorz(:,1:count) = w0*linspace(1,0,count);
        end
        zControlSeed = obj.z0;
        if obj.realtimeRtiEnabled
            zControlSeed = obj.repairRealtimeRtiControlGuess(zControlSeed);
        end
        zSolve0 = obj.repairStateGuess(zControlSeed,xhat,wHorz);
        baseHorizon = repmat(obj.uTrim,1,obj.Nc);
        nlp = obj.assembleWindow(xhat,wHorz,baseHorizon);
        idx = obj.buildIndexMaps();
    end

    %------------------------------------------------------------------
    function setReciprocalPredictionContext(obj,latentState,contextHorizon)
    %SETRECIPROCALPREDICTIONCONTEXT Bind causal memory and future context.
        assert(obj.reciprocalProviderEnabled, ...
            'nMPC:ReciprocalProviderDisabled', ...
            'The reciprocal prediction context requires its audit selector.');
        latentState = latentState(:);
        assert(numel(latentState)== ...
            obj.reciprocalProvider.hiddenStateCount && ...
            all(isfinite(latentState)) && iscell(contextHorizon) && ...
            numel(contextHorizon)==obj.Nc, ...
            'nMPC:ReciprocalPredictionContext', ...
            ['The reciprocal predictor requires one finite latent endpoint ', ...
             'and one future context per control interval.']);
        validationState = [obj.xTrim;latentState];
        referenceRawContext = struct();
        referencePreparedContext = struct();
        for interval = 1:obj.Nc
            assert(isstruct(contextHorizon{interval}) && ...
                isscalar(contextHorizon{interval}), ...
                'nMPC:ReciprocalPredictionContext', ...
                'Future reciprocal context %d is invalid.',interval);
            rawContext = contextHorizon{interval};
            reusePreparedModel = ...
                obj.identicalFutureContextCondensationEnabled && ...
                interval>1 && ...
                ~isempty(fieldnames(referenceRawContext)) && ...
                ~isfield(rawContext,'scheduledModel') && ...
                obj.identicalFutureContextBindingEqual( ...
                    rawContext,referenceRawContext);
            if reusePreparedModel
                rawContext.scheduledModel = ...
                    referencePreparedContext.scheduledModel;
                contextHorizon{interval} = rawContext;
                obj.identicalFutureContextReuseCount = ...
                    obj.identicalFutureContextReuseCount+1;
                obj.identicalFutureContextLastAction = "exact_reuse";
            else
                contextHorizon{interval} = ...
                    obj.prepareScheduledReciprocalContext(rawContext);
                obj.identicalFutureContextPreparationCount = ...
                    obj.identicalFutureContextPreparationCount+1;
                if obj.identicalFutureContextCondensationEnabled
                    if interval==1
                        referenceRawContext = rawContext;
                        referencePreparedContext = ...
                            contextHorizon{interval};
                        obj.identicalFutureContextLastAction = ...
                            "reference_prepared";
                    else
                        obj.identicalFutureContextFallbackCount = ...
                            obj.identicalFutureContextFallbackCount+1;
                        obj.identicalFutureContextLastAction = ...
                            "nonidentical_prepared";
                    end
                end
            end
            obj.reciprocalProvider.propagateControlInterval( ...
                validationState,0,zeros(obj.nu,1), ...
                contextHorizon{interval},false);
        end
        obj.reciprocalInitialLatentState = latentState;
        obj.reciprocalFutureContextHorizon = ...
            reshape(contextHorizon,1,obj.Nc);
    end

    %------------------------------------------------------------------
    function snapshot = identicalFutureContextCondensationAuditSnapshot(obj)
    %IDENTICALFUTURECONTEXTCONDENSATIONAUDITSNAPSHOT Report exact reuse.
        snapshot = struct( ...
            'enabled',obj.identicalFutureContextCondensationEnabled, ...
            'preparationCount',obj.identicalFutureContextPreparationCount, ...
            'reuseCount',obj.identicalFutureContextReuseCount, ...
            'fallbackCount',obj.identicalFutureContextFallbackCount, ...
            'lastAction',obj.identicalFutureContextLastAction);
    end

    %------------------------------------------------------------------
    function applyScheduledReciprocalPacket(obj,packet)
    %APPLYSCHEDULEDRECIPROCALPACKET Refresh query maps without memory reset.
        assert(obj.reciprocalProviderEnabled && ...
            isa(obj.reciprocalProvider, ...
                'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider'), ...
            'nMPC:ScheduledReciprocalProviderDisabled', ...
            'The scheduled reciprocal provider is not active.');
        obj.reciprocalProvider.applyScheduledPacket(obj.model,packet);
    end

    %------------------------------------------------------------------
    function audit = evaluateReciprocalProviderConstraintAudit( ...
            obj,xhat,whatEst,zOverride)
    %EVALUATERECIPROCALPROVIDERCONSTRAINTAUDIT Read-only nMPC NLP hook.
        assert(obj.reciprocalProviderEnabled, ...
            'nMPC:ReciprocalProviderAuditDisabled', ...
            'The reciprocal nMPC audit requires its approved selector.');
        [nlp,z,~,~] = obj.assembleAuditWindow(xhat,whatEst,true);
        if nargin >= 4 && ~isempty(zOverride)
            z = zOverride(:);
            assert(numel(z)==numel(obj.z0) && all(isfinite(z)), ...
                'nMPC:ReciprocalProviderAuditDecision', ...
                'The reciprocal nMPC audit decision is invalid.');
        end
        [cost,gradient,hessian] = nlp.cost(z);
        [inequality,equality,gradientInequality,gradientEquality] = ...
            nlp.nonl(z);
        audit = struct('z',z,'cost',cost,'gradient',gradient, ...
            'hessian',hessian,'inequality',inequality, ...
            'equality',equality, ...
            'gradientInequality',gradientInequality, ...
            'gradientEquality',gradientEquality, ...
            'lowerBound',nlp.lb,'upperBound',nlp.ub, ...
            'nativeStateCount',obj.nativeStateCount, ...
            'latentStateCount', ...
                obj.reciprocalProvider.hiddenStateCount, ...
            'decisionCount',numel(z));
    end

    %------------------------------------------------------------------
    function audit = evaluateReciprocalProviderZeroControlAudit( ...
            obj,xhat,whatEst)
    %EVALUATERECIPROCALPROVIDERZEROCONTROLAUDIT Read-only trim invariance hook.
        assert(obj.reciprocalProviderEnabled, ...
            'nMPC:ReciprocalProviderAuditDisabled', ...
            'The reciprocal nMPC audit requires its approved selector.');
        [~,z,wHorz,idx] = obj.assembleAuditWindow(xhat,whatEst,true);
        for interval = 1:obj.Nc
            z(idx.u{interval}) = 0;
        end
        z = obj.repairStateGuess(z,xhat,wHorz);
        stateHorizon = reshape(z(1:(obj.Nc+1)*obj.nx), ...
            obj.nx,obj.Nc+1);
        audit = struct('stateHorizon',stateHorizon, ...
            'firstDeparture',stateHorizon(:,2)-obj.xTrim, ...
            'terminalDeparture',stateHorizon(:,end)-obj.xTrim, ...
            'initialLatentState',obj.reciprocalInitialLatentState, ...
            'trimLatentState',obj.reciprocalLatentTrim, ...
            'disturbanceHorizon',wHorz);
    end

    %------------------------------------------------------------------
    function audit = evaluateReciprocalProviderSourceDomainAudit( ...
            obj,xhat,whatEst,zOverride)
    %EVALUATERECIPROCALPROVIDERSOURCEDOMAINAUDIT Inspect a proposed horizon.
    %   The returned ratios are source-query relative and match the
    %   PlantRunTime source-domain owner.  This is read-only; command
    %   selection remains unchanged until the audit gate is qualified.
        assert(obj.reciprocalProviderEnabled && isa( ...
            obj.reciprocalProvider, ...
            'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider'), ...
            'nMPC:ScheduledReciprocalProviderAuditDisabled', ...
            'The source-domain audit requires the scheduled provider.');
        useAcceptedHorizon = nargin >= 4 && isempty(zOverride) && ...
            ~isempty(obj.Xhist) && ~isempty(obj.Uhist) && ...
            ~isempty(obj.Whist);
        if useAcceptedHorizon
            X = obj.Xhist;
            U = obj.Uhist;
            wHorz = obj.Whist;
        else
            [~,z,wHorz,idx] = obj.assembleAuditWindow(xhat,whatEst,true);
            X = reshape(z(1:(obj.Nc+1)*obj.nx),obj.nx,obj.Nc+1);
            U = reshape(z(idx.u{1}(1):idx.u{end}(end)),obj.nu,obj.Nc);
        end
        if nargin >= 4 && ~isempty(zOverride)
            z = zOverride(:);
            assert(numel(z)==obj.decisionVariableCount() && ...
                all(isfinite(z)), ...
                'nMPC:SourceDomainAuditDecision', ...
                'The source-domain audit decision is invalid.');
            X = reshape(z(1:(obj.Nc+1)*obj.nx),obj.nx,obj.Nc+1);
            U = reshape(z(idx.u{1}(1):idx.u{end}(end)),obj.nu,obj.Nc);
        end
        assert(isequal(size(X),[obj.nx,obj.Nc+1]) && ...
            isequal(size(U),[obj.nu,obj.Nc]) && ...
            isequal(size(wHorz),[obj.nw,obj.Nc]) && ...
            all(isfinite([X(:);U(:);wHorz(:)])), ...
            'nMPC:SourceDomainAuditHorizon', ...
            'The source-domain audit horizon is invalid.');
        latentHorizon = obj.rebuildReciprocalLatentHorizon(X,U,wHorz);
        interval = cell(1,obj.Nc);
        accepted = true(1,obj.Nc);
        for intervalIndex = 1:obj.Nc
            if intervalIndex == 1
                uPrevious = obj.uPrev;
            else
                uPrevious = U(:,intervalIndex-1);
            end
            uModel = obj.buildIntervalModelControl( ...
                uPrevious,U(:,intervalIndex));
            interval{intervalIndex} = obj.reciprocalProvider. ...
                evaluateControlIntervalSourceDomain( ...
                [X(:,intervalIndex);latentHorizon(:,intervalIndex)], ...
                wHorz(:,intervalIndex),uModel-obj.uModelTrim, ...
                obj.reciprocalFutureContextHorizon{intervalIndex});
            accepted(intervalIndex) = interval{intervalIndex}.accepted;
        end
        audit = struct('accepted',all(accepted), ...
            'firstRejectedInterval',find(~accepted,1,'first'), ...
            'interval',{interval},'stateHorizon',X, ...
            'latentHorizon',latentHorizon,'controlHorizon',U, ...
            'disturbanceHorizon',wHorz);
    end

    %------------------------------------------------------------------
    function configureNativeReducedHorizonCheckpointAudit(obj,request)
    %CONFIGURENATIVEREDUCEDHORIZONCHECKPOINTAUDIT Bind a saved Case-B owner.
        cfg = struct('ctrl',struct( ...
            'nativeReducedHorizonRtiAudit',request));
        obj.configureNativeReducedHorizonRtiAudit(cfg);
    end

    %------------------------------------------------------------------
    function configureHighLeverageRuntimeCheckpointAudit(obj,request)
    %CONFIGUREHIGHLEVERAGERUNTIMECHECKPOINTAUDIT Bind an approved saved E3 owner.
        cfg = struct('ctrl',struct('highLeverageRuntimeAudit',request));
        obj.configureHighLeverageRuntimeAudit(cfg);
    end

    %------------------------------------------------------------------
    function configurePreparedHorizonDataReuseCheckpointAudit(obj,request)
    %CONFIGUREPREPAREDHORIZONDATAREUSECHECKPOINTAUDIT Extend saved V73 exactly.
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled') && logical(request.enabled) && ...
            isfield(request,'auditOnly') && logical(request.auditOnly) && ...
            isfield(request,'changeId') && string(request.changeId)== ...
                "phase18c-v17a-casebc-high-leverage-runtime-audit-v1" && ...
            isfield(request,'scope') && string(request.scope)== ...
                "saved_state_discriminator" && ...
            obj.nativeValueHorizonActive && ...
            obj.nativeCausalRolloutActive && ...
            obj.acceptedLatentHorizonCondensationActive && ...
            obj.realtimeRtiSolutionOwnerEnabled && ...
            obj.realtimeRtiIterationCount==1, ...
            'nMPC:PreparedHorizonDataReuseCheckpointRequest', ...
            'Exact data reuse requires the qualified saved V73 controller.');
        obj.preparedHorizonDataReuseRequested = true;
        obj.preparedHorizonDataReuseActive = true;
    end

    %------------------------------------------------------------------
    function configureAcceptedReplayReuseCheckpointAudit(obj,request)
    %CONFIGUREACCEPTEDREPLAYREUSECHECKPOINTAUDIT Reuse the exact final E3 replay.
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled') && logical(request.enabled) && ...
            isfield(request,'auditOnly') && logical(request.auditOnly) && ...
            isfield(request,'changeId') && string(request.changeId)== ...
                "phase18c-v17a-casebc-high-leverage-runtime-audit-v1" && ...
            isfield(request,'scope') && string(request.scope)== ...
                "saved_state_discriminator" && ...
            obj.nativeValueHorizonActive && ...
            obj.acceptedLatentHorizonCondensationActive && ...
            obj.realtimeRtiSolutionOwnerEnabled && ...
            obj.realtimeRtiIterationCount==1, ...
            'nMPC:AcceptedReplayReuseCheckpointRequest', ...
            'Exact replay reuse requires the qualified saved V73 controller.');
        obj.acceptedReplayReuseRequested = true;
        obj.acceptedReplayReuseActive = true;
    end

    %------------------------------------------------------------------
    function configureIdenticalFutureContextCondensationCheckpointAudit( ...
            obj,request)
    %CONFIGUREIDENTICALFUTURECONTEXTCONDENSATIONCHECKPOINTAUDIT Bind C1.
        cfg = struct('ctrl',struct( ...
            'identicalFutureContextRealizationCondensationAudit',request));
        obj.configureIdenticalFutureContextCondensationAudit(cfg);
    end

    %------------------------------------------------------------------
    function audit = evaluateRealtimeRtiSeedAudit(obj,xhat,whatEst)
    %EVALUATEREALTIMERTISEEDAUDIT Read-only controller seed evaluation.
        assert(obj.realtimeRtiEnabled, ...
            'nMPC:RealtimeRtiAuditDisabled', ...
            'The condensed nMPC RTI audit requires its approved selector.');
        [nlp,zInitial,wHorizon] = obj.assembleAuditWindow( ...
            xhat,whatEst,true);
        [candidate,info] = obj.buildRealtimeRtiSeed( ...
            nlp,zInitial,wHorizon);
        audit = struct('initial',zInitial,'candidate',candidate, ...
            'disturbanceHorizon',wHorizon,'info',info);
    end

    %------------------------------------------------------------------
    end % public methods

    %==================================================================
    methods(Access=private)
    %------------------------------------------------------------------
    function configureRealtimeRtiAudit(obj,cfg)
    %CONFIGUREREALTIMERTIAUDIT Install the approved default-inactive seed.
        legacyField = 'realtimeRateChartRtiAudit';
        ownerField = 'condensedRtiRuntimeOwner';
        if ~isfield(cfg,'ctrl')
            return
        end
        legacyPresent = isfield(cfg.ctrl,legacyField) && ...
            ~isempty(cfg.ctrl.(legacyField));
        ownerPresent = isfield(cfg.ctrl,ownerField) && ...
            ~isempty(cfg.ctrl.(ownerField));
        if ~legacyPresent && ~ownerPresent
            return
        end
        assert(~(legacyPresent && ownerPresent), ...
            'nMPC:RealtimeRtiRequestConflict', ...
            'Legacy RTI seeding and RTI solution ownership are exclusive.');
        if ownerPresent
            request = cfg.ctrl.(ownerField);
        else
            request = cfg.ctrl.(legacyField);
        end
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled'), 'nMPC:RealtimeRtiRequest', ...
            'The real-time RTI audit request requires enabled.');
        enabled = request.enabled;
        assert(isscalar(enabled) && (islogical(enabled) || ...
            (isnumeric(enabled) && isfinite(enabled) && ...
            ismember(enabled,[0 1]))), 'nMPC:RealtimeRtiRequest', ...
            'The real-time RTI enabled selector must be logical.');
        if ~logical(enabled)
            return
        end
        if ownerPresent
            required = {'auditOnly','changeId','chartMode','condensingMode', ...
                'solutionOwner','nmpcIterationCount','fallbackPolicy'};
            assert(all(isfield(request,required)) && ...
                ~logical(request.auditOnly) && ...
                string(request.changeId)== ...
                    "phase18c-v17a-condensed-rti-compiled-runtime-owner-promotion-v1" && ...
                string(request.chartMode)=="auto_rate" && ...
                string(request.condensingMode)=="full_state" && ...
                string(request.solutionOwner)=="qualified_rti" && ...
                string(request.fallbackPolicy)=="hold_previous", ...
                'nMPC:RealtimeRtiRequest', ...
                'The enabled nMPC RTI solution-owner request is unauthorized.');
            iterationCount = double(request.nmpcIterationCount);
            obj.realtimeRtiSolutionOwnerEnabled = true;
        else
            required = {'auditOnly','changeId','chartMode','condensingMode', ...
                'iterationCount','correctionSolver'};
            assert(all(isfield(request,required)) && ...
                logical(request.auditOnly) && ...
                string(request.changeId)== ...
                    "phase18c-v17a-realtime-rate-chart-partial-condensing-rti-audit-v1" && ...
                string(request.chartMode)=="auto_rate" && ...
                string(request.condensingMode)=="full_state" && ...
                string(request.correctionSolver)=="configured_full_nlp", ...
                'nMPC:RealtimeRtiRequest', ...
                'The enabled real-time RTI audit request is not authorized.');
            iterationCount = double(request.iterationCount);
        end
        assert(isscalar(iterationCount) && isfinite(iterationCount) && ...
            iterationCount==fix(iterationCount) && ...
            iterationCount>=1 && iterationCount<=3, ...
            'nMPC:RealtimeRtiIterationCount', ...
            'The real-time RTI audit permits one to three iterations.');
        assert(obj.n_surf==2 && obj.nu==4, ...
            'nMPC:RealtimeRtiControlDimension', ...
            'The approved rate chart requires two endpoint/rate surface pairs.');

        rtiOptions = AeroFlex.optim.SQPSolver.defaultOptions();
        rtiOptions.ElasticMode = false;
        rtiOptions.CheckGradientsOnce = false;
        rtiOptions.QPOptions = optimoptions('quadprog', ...
            'Algorithm','interior-point-convex','Display','off', ...
            'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9, ...
            'StepTolerance',1e-12,'MaxIterations',200);
        if isfield(request,'qpMode')
            qpMode = string(request.qpMode);
            assert(isscalar(qpMode) && ismember(qpMode,[ ...
                "configured_quadprog","direct_then_active_set"]), ...
                'nMPC:RealtimeRtiQpMode', ...
                'The condensed nMPC RTI QP mode is unsupported.');
            rtiOptions.CondensedRtiQpMode = qpMode;
        end
        obj.realtimeRtiSolver = AeroFlex.optim.SQPSolver(rtiOptions);
        obj.realtimeRtiEnabled = true;
        obj.realtimeRtiChangeId = string(request.changeId);
        obj.realtimeRtiIterationCount = iterationCount;
        obj.realtimeRtiChartMode = string(request.chartMode);
    end

    %------------------------------------------------------------------
    function configureNativeReducedHorizonRtiAudit(obj,cfg)
    %CONFIGURENATIVEREDUCEDHORIZONRTIAUDIT Bind scheduled-only H2 dispatch.
        fieldName = 'nativeReducedHorizonRtiAudit';
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl,fieldName) || ...
                isempty(cfg.ctrl.(fieldName))
            return
        end
        request = cfg.ctrl.(fieldName);
        required = {'enabled','auditOnly','changeId','controllerEnabled', ...
            'controllerBinarySha256'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)) && logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1", ...
            'nMPC:NativeReducedHorizonRequest', ...
            'The native reduced-horizon controller request is unauthorized.');
        obj.nativeReducedHorizonRtiRequested = ...
            logical(request.enabled) && logical(request.controllerEnabled);
        if ~obj.nativeReducedHorizonRtiRequested
            return
        end
        if ~obj.realtimeRtiEnabled || obj.Nc~=10 || obj.nx~=74 || ...
                obj.nw~=1 || ~obj.symmetricSurfaceSubspaceEnabled || ...
                ~obj.reciprocalProviderEnabled || ...
                ~isa(obj.reciprocalProvider, ...
                    'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            obj.nativeReducedHorizonRtiLastFallback = ...
                "unsupported controller architecture; exact RTI retained";
            return
        end
        functionName = ...
            'AeroFlex_ctrl_scheduledReciprocalControllerReducedTangentHorizonAudit_mex';
        kernelPath = string(which(functionName));
        expectedHash = lower(string(request.controllerBinarySha256));
        if kernelPath=="" || strlength(expectedHash)~=64 || ...
                obj.localFileHash(kernelPath)~=expectedHash
            obj.nativeReducedHorizonRtiLastFallback = ...
                "controller horizon MEX unavailable or hash mismatch";
            return
        end
        obj.nativeReducedHorizonRtiKernel = str2func(functionName);
        obj.nativeReducedHorizonRtiKernelIdentity = functionName+"|"+ ...
            expectedHash+"|"+string(computer('arch'))+"|"+string(mexext);
        obj.nativeReducedHorizonRtiActive = true;
    end

    %------------------------------------------------------------------
    function configureHighLeverageRuntimeAudit(obj,cfg)
    %CONFIGUREHIGHLEVERAGERUNTIMEAUDIT Bind the scheduled-only E3 value MEX.
        fieldName = 'highLeverageRuntimeAudit';
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl,fieldName) || ...
                isempty(cfg.ctrl.(fieldName))
            return
        end
        request = cfg.ctrl.(fieldName);
        required = {'enabled','auditOnly','changeId','e3Enabled', ...
            'controllerValueBinarySha256'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)) && logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casebc-high-leverage-runtime-audit-v1", ...
            'nMPC:HighLeverageRuntimeRequest', ...
            'The high-leverage controller request is unauthorized.');
        obj.acceptedLatentHorizonCondensationRequested = false;
        obj.acceptedLatentHorizonCondensationActive = false;
        obj.nativeCausalRolloutRequested = false;
        obj.nativeCausalRolloutActive = false;
        obj.preparedHorizonDataReuseRequested = false;
        obj.preparedHorizonDataReuseActive = false;
        obj.acceptedReplayReuseRequested = false;
        obj.acceptedReplayReuseActive = false;
        obj.nativeValueHorizonRequested = ...
            logical(request.enabled) && logical(request.e3Enabled);
        if ~obj.nativeValueHorizonRequested
            return
        end
        if ~obj.nativeReducedHorizonRtiActive || obj.Nc~=10 || ...
                obj.nx~=74 || obj.nw~=1 || ...
                ~obj.symmetricSurfaceSubspaceEnabled || ...
                ~obj.reciprocalProviderEnabled || ...
                ~isa(obj.reciprocalProvider, ...
                    'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            obj.nativeValueHorizonLastFallback = ...
                "unsupported controller architecture; H2 replay retained";
            return
        end
        functionName = ...
            'AeroFlex_ctrl_scheduledReciprocalControllerValueHorizonAudit_mex';
        kernelPath = string(which(functionName));
        expectedHash = lower(string(request.controllerValueBinarySha256));
        if kernelPath=="" || strlength(expectedHash)~=64 || ...
                obj.localFileHash(kernelPath)~=expectedHash
            obj.nativeValueHorizonLastFallback = ...
                "controller value-horizon MEX unavailable or hash mismatch";
            return
        end
        obj.nativeValueHorizonKernel = str2func(functionName);
        obj.nativeValueHorizonKernelIdentity = functionName+"|"+ ...
            expectedHash+"|"+string(computer('arch'))+"|"+string(mexext);
        obj.nativeValueHorizonActive = true;
        if isfield(request,'controllerCausalRolloutEnabled') && ...
                logical(request.controllerCausalRolloutEnabled)
            assert(isfield(request,'controllerCausalRolloutBinarySha256'), ...
                'nMPC:NativeCausalRolloutRequest', ...
                'The controller causal-rollout binary hash is required.');
            obj.nativeCausalRolloutRequested = true;
            causalFunctionName = ...
                'AeroFlex_ctrl_scheduledReciprocalControllerCausalRolloutAudit_mex';
            causalKernelPath = string(which(causalFunctionName));
            causalExpectedHash = lower(string( ...
                request.controllerCausalRolloutBinarySha256));
            if causalKernelPath=="" || strlength(causalExpectedHash)~=64 || ...
                    obj.localFileHash(causalKernelPath)~=causalExpectedHash
                obj.nativeCausalRolloutLastFallback = ...
                    "controller causal-rollout MEX unavailable or hash mismatch";
            else
                obj.nativeCausalRolloutKernel = str2func(causalFunctionName);
                obj.nativeCausalRolloutKernelIdentity = ...
                    causalFunctionName+"|"+causalExpectedHash+"|"+ ...
                    string(computer('arch'))+"|"+string(mexext);
                obj.nativeCausalRolloutActive = true;
            end
        end
        if isfield(request,'controllerLatentHorizonCondensationEnabled') && ...
                logical(request.controllerLatentHorizonCondensationEnabled)
            assert(obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled, ...
                'nMPC:AcceptedLatentHorizonCondensationRequest', ...
                ['Exact accepted latent-horizon condensation requires ', ...
                 'the qualified online condensed-RTI solution owner.']);
            obj.acceptedLatentHorizonCondensationRequested = true;
            obj.acceptedLatentHorizonCondensationActive = true;
        end
        if isfield(request,'preparedHorizonDataReuseEnabled') && ...
                logical(request.preparedHorizonDataReuseEnabled)
            assert(isfield(request,'preparedHorizonDataReuseScope') && ...
                string(request.preparedHorizonDataReuseScope)== ...
                    "saved_state_discriminator" && ...
                obj.nativeCausalRolloutActive && ...
                obj.acceptedLatentHorizonCondensationActive && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                isfield(request,'r1Enabled') && ...
                logical(request.r1Enabled) && ...
                isfield(request,'nmpcIterationCount') && ...
                double(request.nmpcIterationCount)==1, ...
                'nMPC:PreparedHorizonDataReuseRequest', ...
                ['Exact prepared-horizon data reuse requires the ', ...
                 'qualified V73 controller owner.']);
            obj.preparedHorizonDataReuseRequested = true;
            obj.preparedHorizonDataReuseActive = true;
        end
        if isfield(request,'acceptedReplayReuseEnabled') && ...
                logical(request.acceptedReplayReuseEnabled)
            assert(isfield(request,'acceptedReplayReuseScope') && ...
                string(request.acceptedReplayReuseScope)== ...
                    "saved_state_discriminator" && ...
                obj.acceptedLatentHorizonCondensationActive && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                isfield(request,'r1Enabled') && ...
                logical(request.r1Enabled) && ...
                isfield(request,'nmpcIterationCount') && ...
                double(request.nmpcIterationCount)==1, ...
                'nMPC:AcceptedReplayReuseRequest', ...
                ['Exact accepted-replay reuse requires the qualified ', ...
                 'V73 controller owner.']);
            obj.acceptedReplayReuseRequested = true;
            obj.acceptedReplayReuseActive = true;
        end
        if isfield(request,'r1Enabled') && logical(request.r1Enabled)
            assert(isfield(request,'nmpcIterationCount') && ...
                isfield(request,'r1Scope') && ...
                string(request.r1Scope)=="saved_state_discriminator" && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                obj.realtimeRtiIterationCount==2 && ...
                isscalar(request.nmpcIterationCount) && ...
                double(request.nmpcIterationCount)==1, ...
                'nMPC:HighLeverageR1Request', ...
                ['The one-correction controller discriminator requires ', ...
                 'the qualified two-correction H2/E3 owner.']);
            obj.realtimeRtiIterationCount = 1;
        end
    end

    %------------------------------------------------------------------
    function configureIdenticalFutureContextCondensationAudit(obj,cfg)
    %CONFIGUREIDENTICALFUTURECONTEXTCONDENSATIONAUDIT Bind exact C1 reuse.
        fieldName = ...
            'identicalFutureContextRealizationCondensationAudit';
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl,fieldName) || ...
                isempty(cfg.ctrl.(fieldName))
            return
        end
        request = cfg.ctrl.(fieldName);
        required = {'enabled','auditOnly','changeId','caseId', ...
            'contextPolicy','expectedContextCount'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)), ...
            'nMPC:IdenticalFutureContextRequest', ...
            'The identical-future-context request is incomplete.');
        obj.identicalFutureContextCondensationEnabled = false;
        obj.identicalFutureContextLastAction = "disabled";
        if ~logical(request.enabled)
            return
        end
        assert(logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casebc-identical-future-context-realization-condensation-v1" && ...
            string(request.caseId)=="formal_case_b" && ...
            string(request.contextPolicy)== ...
                "held_measured_rigid_and_held_outer_tail_thrust" && ...
            double(request.expectedContextCount)==obj.Nc && obj.Nc==10 && ...
            obj.reciprocalProviderEnabled && ...
            isa(obj.reciprocalProvider, ...
                'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider') && ...
            obj.nativeReducedHorizonRtiActive && ...
            obj.nativeValueHorizonActive && ...
            obj.realtimeRtiSolutionOwnerEnabled && ...
            obj.realtimeRtiIterationCount==1, ...
            'nMPC:IdenticalFutureContextApproval', ...
            ['Exact future-context condensation requires the qualified ', ...
             'scheduled H2/E3/R1 ten-interval Case-B owner.']);
        obj.identicalFutureContextCondensationEnabled = true;
        obj.identicalFutureContextPreparationCount = 0;
        obj.identicalFutureContextReuseCount = 0;
        obj.identicalFutureContextFallbackCount = 0;
        obj.identicalFutureContextLastAction = "configured";
    end

    %------------------------------------------------------------------
    function configureScheduledSourceDomainRtiConstraintAudit(obj,cfg)
    %CONFIGURESCHEDULEDSOURCEDOMAINRTICONSTRAINTAUDIT Install the Case-B gate.
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl, ...
                'scheduledSourceDomainRtiConstraintAudit') || ...
                isempty(cfg.ctrl.scheduledSourceDomainRtiConstraintAudit)
            return
        end
        request = cfg.ctrl.scheduledSourceDomainRtiConstraintAudit;
        required = {'enabled','auditOnly','changeId','caseId', ...
            'replayAllIntervals','fallbackPolicy'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)), ...
            'nMPC:ScheduledSourceDomainRtiRequest', ...
            'The scheduled source-domain RTI request is incomplete.');
        if ~logical(request.enabled)
            return
        end
        assert(obj.realtimeRtiSolutionOwnerEnabled && ...
            obj.reciprocalProviderEnabled && ...
            isa(obj.reciprocalProvider, ...
            'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider') && ...
            logical(request.auditOnly) && ...
            string(request.changeId)== ...
            "phase18c-v17a-caseb-scheduled-source-domain-rti-constraint-audit-v1" && ...
            string(request.caseId)=="formal_case_b" && ...
            logical(request.replayAllIntervals) && ...
            string(request.fallbackPolicy)=="trim_relative", ...
            'nMPC:ScheduledSourceDomainRtiRequest', ...
            ['The scheduled source-domain gate is authorized only for the ', ...
             'approved Case-B reciprocal RTI audit contract.']);
        obj.scheduledSourceDomainRtiConstraintEnabled = true;
        obj.scheduledSourceDomainRtiConstraintChangeId = ...
            string(request.changeId);
    end

    %------------------------------------------------------------------
    function info = qualifyScheduledSourceDomainRtiCandidate( ...
            obj,info,candidate,xhat,whatEst)
    %QUALIFYSCHEDULEDSOURCEDOMAINRTICANDIDATE Reject invalid replay horizons.
        info.sourceDomainGateEnabled = ...
            obj.scheduledSourceDomainRtiConstraintEnabled;
        info.sourceDomainAccepted = true;
        info.sourceDomainRejected = false;
        info.sourceDomainFirstRejectedInterval = nan;
        if ~obj.scheduledSourceDomainRtiConstraintEnabled || ...
                ~isfield(info,'qualified') || ~logical(info.qualified)
            return
        end
        try
            audit = obj.evaluateReciprocalProviderSourceDomainAudit( ...
                xhat,whatEst,candidate);
            info.sourceDomainMovingBaseline = ...
                obj.evaluateSourceDomainMovingBaseline( ...
                candidate,xhat,whatEst,audit);
            info.sourceDomainAccepted = logical(audit.accepted);
            info.sourceDomainFirstRejectedInterval = ...
                double(audit.firstRejectedInterval);
            if ~audit.accepted
                info.qualified = false;
                info.sourceDomainRejected = true;
                info.message = "Condensed RTI source-domain replay rejected.";
            end
        catch replayException
            info.qualified = false;
            info.sourceDomainAccepted = false;
            info.sourceDomainRejected = true;
            info.sourceDomainFirstRejectedInterval = nan;
            info.message = "Condensed RTI source-domain replay failed closed: " + ...
                string(replayException.message);
            info.identifier = string(replayException.identifier);
        end
    end

    %------------------------------------------------------------------
    function result = evaluateSourceDomainMovingBaseline( ...
            obj,candidate,xhat,whatEst,candidateAudit)
    %EVALUATESOURCEDOMAINMOVINGBASELINE Compare a candidate with zero inner action.
    %   This audit does not alter source-domain acceptance.  It exposes
    %   whether a static-envelope rejection is already present in the
    %   scheduled zero-inner continuation at a moving Case-B condition.
        [~,~,wHorz,~] = obj.assembleAuditWindow(xhat,whatEst,true);
        baseline = obj.resetRealtimeRtiControlHorizon(candidate);
        baseline = obj.repairStateGuess(baseline,xhat,wHorz);
        baselineAudit = obj.evaluateReciprocalProviderSourceDomainAudit( ...
            xhat,whatEst,baseline);
        result = struct( ...
            'candidate',obj.sourceDomainAuditSummary(candidateAudit), ...
            'zeroInnerBaseline',obj.sourceDomainAuditSummary(baselineAudit), ...
            'stateDifferenceInfinity',max(abs( ...
                candidateAudit.stateHorizon-baselineAudit.stateHorizon), ...
                [],'all'), ...
            'latentDifferenceInfinity',max(abs( ...
                candidateAudit.latentHorizon-baselineAudit.latentHorizon), ...
                [],'all'), ...
            'controlDifferenceInfinity',max(abs( ...
                candidateAudit.controlHorizon-baselineAudit.controlHorizon), ...
                [],'all'));
    end

    %------------------------------------------------------------------
    function summary = sourceDomainAuditSummary(~,audit)
    %SOURCEDOMAINAUDITSUMMARY Retain scalar source-domain diagnostics.
        maximumInputRatio = 0;
        maximumVisibleRatio = 0;
        maximumHiddenRatio = 0;
        maximumRootWrenchRatio = 0;
        for intervalIndex = 1:numel(audit.interval)
            interval = audit.interval{intervalIndex};
            for substepIndex = 1:numel(interval.substeps)
                details = interval.substeps(substepIndex).details;
                if isempty(details)
                    continue
                end
                maximumInputRatio = max(maximumInputRatio, ...
                    max([details.inputRatio]));
                maximumVisibleRatio = max(maximumVisibleRatio, ...
                    max([details.visibleRatio]));
                maximumHiddenRatio = max(maximumHiddenRatio, ...
                    max([details.hiddenRatio]));
                maximumRootWrenchRatio = max(maximumRootWrenchRatio, ...
                    max([details.rootWrenchRatio]));
            end
        end
        summary = struct('accepted',logical(audit.accepted), ...
            'firstRejectedInterval',double(audit.firstRejectedInterval), ...
            'maximumInputRatio',maximumInputRatio, ...
            'maximumVisibleRatio',maximumVisibleRatio, ...
            'maximumHiddenRatio',maximumHiddenRatio, ...
            'maximumRootWrenchRatio',maximumRootWrenchRatio);
    end

    %------------------------------------------------------------------
    function metric = realtimeRtiFirstOrderMetric(~,info)
        metric = nan;
        if isfield(info,'iterations') && ~isempty(info.iterations)
            values = [info.iterations.reducedKktInf];
            values = values(isfinite(values));
            if ~isempty(values)
                metric = values(end);
            end
        end
    end

    %------------------------------------------------------------------
    function lambda = emptySolverMultipliers(~)
        lambda = struct('eqnonlin',zeros(0,1), ...
            'ineqnonlin',zeros(0,1),'lower',zeros(0,1), ...
            'upper',zeros(0,1));
    end

    %------------------------------------------------------------------
    function configureSymmetricSurfaceSubspace(obj,cfg)
    %CONFIGURESYMMETRICSURFACESUBSPACE Install a case-qualified audit policy.
        if ~isfield(cfg,'ctrl')
            return
        end
        hasCaseA = isfield(cfg.ctrl,'formalCaseASymmetricSurfaceAudit') && ...
            ~isempty(cfg.ctrl.formalCaseASymmetricSurfaceAudit);
        hasCaseB = isfield(cfg.ctrl,'formalCaseBSymmetricSurfaceAudit') && ...
            ~isempty(cfg.ctrl.formalCaseBSymmetricSurfaceAudit);
        if ~(hasCaseA || hasCaseB)
            return
        end
        if hasCaseA && hasCaseB
            requestA = cfg.ctrl.formalCaseASymmetricSurfaceAudit;
            requestB = cfg.ctrl.formalCaseBSymmetricSurfaceAudit;
            assert(isstruct(requestA) && isscalar(requestA) && ...
                isfield(requestA,'enabled') && isstruct(requestB) && ...
                isscalar(requestB) && isfield(requestB,'enabled') && ...
                ~(logical(requestA.enabled) && logical(requestB.enabled)), ...
                'nMPC:SymmetricSurfaceRequest', ...
                'Only one case-qualified symmetric-surface audit may be active.');
        end
        if hasCaseB && logical(cfg.ctrl.formalCaseBSymmetricSurfaceAudit.enabled)
            request = cfg.ctrl.formalCaseBSymmetricSurfaceAudit;
            expectedChangeId = ...
                "phase18c-v17a-caseb-nmpc-symmetric-surface-subspace-extension-v1";
            expectedCaseId = "formal_case_b";
        else
            request = cfg.ctrl.formalCaseASymmetricSurfaceAudit;
            expectedChangeId = ...
                "phase18c-v17a-casea-nmpc-symmetric-surface-subspace-audit-v1";
            expectedCaseId = "formal_case_a";
        end
        required = {'enabled','auditOnly','changeId','caseId', ...
            'surfaceMode','sourceMap'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)), ...
            'nMPC:SymmetricSurfaceRequest', ...
            'The formal Case-A symmetric-surface request is incomplete.');
        enabled = request.enabled;
        auditOnly = request.auditOnly;
        assert(isscalar(enabled) && (islogical(enabled) || ...
            (isnumeric(enabled) && isfinite(enabled) && ...
            ismember(enabled,[0 1]))) && ...
            isscalar(auditOnly) && (islogical(auditOnly) || ...
            (isnumeric(auditOnly) && isfinite(auditOnly) && ...
            ismember(auditOnly,[0 1]))), ...
            'nMPC:SymmetricSurfaceRequest', ...
            'The enabled and auditOnly selectors must be logical scalars.');
        if ~logical(enabled)
            return
        end
        sourceMap = double(request.sourceMap(:));
        changeId = string(request.changeId);
        caseId = lower(string(request.caseId));
        surfaceMode = lower(string(request.surfaceMode));
        assert(logical(auditOnly) && obj.n_surf==2 && ...
            changeId==expectedChangeId && caseId==expectedCaseId && ...
            surfaceMode=="symmetric_longitudinal" && ...
            isequal(sourceMap,[1;1]), ...
            'nMPC:SymmetricSurfaceRequest', ...
            ['The symmetric-surface policy is authorized only for the ', ...
             'approved case-qualified [1;1] audit contract.']);
        obj.symmetricSurfaceSubspaceEnabled = true;
        obj.symmetricSurfaceSubspaceCaseId = caseId;
        obj.symmetricSurfaceSubspaceMap = sourceMap;
        obj.symmetricSurfaceSubspaceChangeId = changeId;
    end

    %------------------------------------------------------------------
    function configureReciprocalProvider(obj,cfg)
    %CONFIGURERECIPROCALPROVIDER Install the default-inactive audit model.
        hasScheduledRequest = isfield(cfg,'ctrl') && isfield(cfg.ctrl, ...
            'scheduledReciprocalControllerModelProviderAudit') && ...
            ~isempty(cfg.ctrl.scheduledReciprocalControllerModelProviderAudit) && ...
            isfield(cfg.ctrl.scheduledReciprocalControllerModelProviderAudit, ...
                'enabled') && logical(cfg.ctrl. ...
                scheduledReciprocalControllerModelProviderAudit.enabled);
        hasExactRequest = isfield(cfg,'ctrl') && isfield(cfg.ctrl, ...
            'reciprocalControllerModelProviderAudit') && ...
            ~isempty(cfg.ctrl.reciprocalControllerModelProviderAudit) && ...
            isfield(cfg.ctrl.reciprocalControllerModelProviderAudit, ...
                'enabled') && ...
            logical(cfg.ctrl.reciprocalControllerModelProviderAudit.enabled);
        assert(~(hasScheduledRequest && hasExactRequest), ...
            'nMPC:ReciprocalProviderConflict', ...
            ['Select either the exact-source or scheduled reciprocal ', ...
             'provider, not both.']);
        if hasScheduledRequest
            request = cfg.ctrl. ...
                scheduledReciprocalControllerModelProviderAudit;
            required = {'enabled','auditOnly','changeId','initialPacket', ...
                'initialContext','formulation'};
            assert(isstruct(request) && isscalar(request), ...
                'nMPC:ScheduledReciprocalProviderRequest', ...
                'The scheduled reciprocal provider request must be scalar.');
            missing = required(~isfield(request,required));
            assert(isempty(missing), ...
                'nMPC:ScheduledReciprocalProviderRequest', ...
                'The scheduled reciprocal provider request lacks: %s.', ...
                strjoin(missing,', '));
            assert(logical(request.auditOnly) && string(request.changeId)== ...
                "phase18c-v17a-casebc-scheduled-reciprocal-runtime-integration-v1", ...
                'nMPC:ScheduledReciprocalProviderApproval', ...
                'The scheduled reciprocal provider approval binding changed.');
            obj.reciprocalProviderEnabled = logical(request.enabled);
            if ~obj.reciprocalProviderEnabled
                return
            end
            formulation = lower(string(request.formulation));
            assert(isscalar(formulation) && ...
                formulation=="condensed_internal_memory", ...
                'nMPC:ScheduledReciprocalProviderFormulation', ...
                ['The scheduled nMPC provider requires condensed ', ...
                 'internal memory.']);
            packet = request.initialPacket;
            obj.reciprocalProvider = ...
                AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider( ...
                    obj.model,packet.members,packet.queryTrim);
            if localScheduledAggregateRequested(request)
                obj.reciprocalProvider.configureScheduledAggregateRuntime( ...
                    struct('enabled',true,'auditOnly',true, ...
                        'changeId',string(request.changeId)));
                obj.reciprocalProvider.applyScheduledPacket( ...
                    obj.model,packet);
            end
            if isfield(cfg.ctrl,'compiledReciprocalIntervalProvider') && ...
                    ~isempty(cfg.ctrl.compiledReciprocalIntervalProvider)
                compiledRequest = ...
                    cfg.ctrl.compiledReciprocalIntervalProvider;
                if isfield(compiledRequest,'controllerEnabled')
                    compiledRequest.enabled = logical( ...
                        compiledRequest.enabled) && logical( ...
                        compiledRequest.controllerEnabled);
                end
                obj.reciprocalProvider.configureFixedIntervalKernelAudit( ...
                    compiledRequest);
            end
            if isfield(request,'packetCacheAudit')
                obj.reciprocalProvider. ...
                    configureValidatedPacketCacheAudit( ...
                        request.packetCacheAudit);
            end
            obj.reciprocalProviderFormulation = formulation;
            obj.reciprocalLatentIndex = obj.nativeStateCount+ ...
                (1:obj.reciprocalProvider.hiddenStateCount);
            augmentedTrim = obj.reciprocalProvider.initialize(obj.xTrim);
            obj.reciprocalLatentTrim = ...
                augmentedTrim(obj.reciprocalLatentIndex);
            obj.reciprocalInitialLatentState = obj.reciprocalLatentTrim;
            obj.reciprocalPredictedLatentHorizon = repmat( ...
                obj.reciprocalLatentTrim,1,obj.Nc+1);
            obj.reciprocalFutureContextHorizon = ...
                repmat({request.initialContext},1,obj.Nc);
            if isfield(request,'unboundedAerodynamicLagStates') && ...
                    ~isempty(request.unboundedAerodynamicLagStates)
                selector = request.unboundedAerodynamicLagStates;
                assert(isscalar(selector) && (islogical(selector) || ...
                    (isnumeric(selector) && isfinite(selector) && ...
                    ismember(selector,[0 1]))), ...
                    'nMPC:ReciprocalAerodynamicLagBounds', ...
                    ['The reciprocal aerodynamic-lag bound selector must ', ...
                     'be a logical scalar.']);
                obj.reciprocalUnboundedAerodynamicLagStates = ...
                    logical(selector);
            end
            return
        end
        if ~hasExactRequest
            return
        end
        request = cfg.ctrl.reciprocalControllerModelProviderAudit;
        required = {'enabled','auditOnly','changeId','candidatePath', ...
            'candidateSha256','initialContext'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)) && logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casea-reciprocal-controller-model-provider-audit-v1", ...
            'nMPC:ReciprocalProviderRequest', ...
            'The reciprocal provider request is incomplete or unauthorized.');
        obj.reciprocalProviderEnabled = logical(request.enabled);
        if ~obj.reciprocalProviderEnabled
            return
        end
        if isfield(request,'formulation') && ~isempty(request.formulation)
            formulation = lower(string(request.formulation));
        else
            formulation = "augmented_reference";
        end
        assert(isscalar(formulation) && ...
            formulation=="condensed_internal_memory", ...
            'nMPC:ReciprocalProviderFormulation', ...
            ['The nMPC reciprocal predictor requires the condensed ', ...
             'internal-memory formulation.']);
        candidatePath = char(string(request.candidatePath));
        expectedHash = lower(string(request.candidateSha256));
        assert(isfile(candidatePath) && ...
            obj.localFileHash(candidatePath)==expectedHash, ...
            'nMPC:ReciprocalProviderCandidateHash', ...
            'The reciprocal provider candidate is unavailable or changed.');
        loaded = load(candidatePath,'candidate');
        assert(isfield(loaded,'candidate'), ...
            'nMPC:ReciprocalProviderCandidate', ...
            'The reciprocal provider artifact does not contain candidate.');
        obj.reciprocalProvider = ...
            AeroFlex.ctrl.ReciprocalControllerModelProvider( ...
                obj.model,loaded.candidate,obj.uModelTrim);
        if isfield(cfg.ctrl,'compiledReciprocalIntervalProvider') && ...
                ~isempty(cfg.ctrl.compiledReciprocalIntervalProvider)
            obj.reciprocalProvider.configureFixedIntervalKernelAudit( ...
                cfg.ctrl.compiledReciprocalIntervalProvider);
        end
        obj.reciprocalProviderFormulation = formulation;
        obj.reciprocalLatentIndex = obj.nativeStateCount+ ...
            (1:obj.reciprocalProvider.hiddenStateCount);
        augmentedTrim = obj.reciprocalProvider.initialize(obj.xTrim);
        obj.reciprocalLatentTrim = ...
            augmentedTrim(obj.reciprocalLatentIndex);
        obj.reciprocalInitialLatentState = obj.reciprocalLatentTrim;
        obj.reciprocalPredictedLatentHorizon = repmat( ...
            obj.reciprocalLatentTrim,1,obj.Nc+1);
        obj.reciprocalFutureContextHorizon = ...
            repmat({request.initialContext},1,obj.Nc);
        if isfield(request,'unboundedAerodynamicLagStates') && ...
                ~isempty(request.unboundedAerodynamicLagStates)
            selector = request.unboundedAerodynamicLagStates;
            assert(isscalar(selector) && (islogical(selector) || ...
                (isnumeric(selector) && isfinite(selector) && ...
                ismember(selector,[0 1]))), ...
                'nMPC:ReciprocalAerodynamicLagBounds', ...
                ['The reciprocal aerodynamic-lag bound selector must be ', ...
                 'a logical scalar.']);
            obj.reciprocalUnboundedAerodynamicLagStates = ...
                logical(selector);
        end
    end

    function context = prepareScheduledReciprocalContext(obj,context)
    %PREPARESCHEDULEDRECIPROCALCONTEXT Retain one future native model.
        if ~isa(obj.reciprocalProvider, ...
                'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            return
        end
        required = {'scheduledPacket','scheduledPackage'};
        assert(isstruct(context) && isscalar(context) && ...
            all(isfield(context,required)), ...
            'nMPC:ScheduledReciprocalContext', ...
            ['Each future scheduled interval requires its packet and ', ...
             'native package.']);
        package = context.scheduledPackage;
        assert(isstruct(package) && isscalar(package) && ...
            all(isfield(package,{'beam','aero','base','L','idx'})), ...
            'nMPC:ScheduledReciprocalPackage', ...
            'The future scheduled package is incomplete.');
        integrator = obj.cfg.modelHandle( ...
            obj.cfg,package.beam,package.aero,package.base);
        context.scheduledModel = AeroFlex.sched.applyToROMIntegrator( ...
            integrator,package,obj.cfg);
    end

    function tf = identicalFutureContextBindingEqual(obj,left,right)
    %IDENTICALFUTURECONTEXTBINDINGEQUAL Compare the exact compact owner key.
        tf = false;
        contextFields = {'elevatorIncrement','thrustIncrement', ...
            'rigidState','ownerPolicy','scheduledReciprocal', ...
            'scheduledPacket','scheduledPackage'};
        if ~isstruct(left) || ~isscalar(left) || ...
                ~isstruct(right) || ~isscalar(right) || ...
                ~all(isfield(left,contextFields)) || ...
                ~all(isfield(right,contextFields))
            return
        end
        leftPacket = left.scheduledPacket;
        rightPacket = right.scheduledPacket;
        packetFields = {'schemaVersion','queryPackageHash', ...
            'bankManifestSha256','sourceIds','schedule','queryTrim', ...
            'hiddenStateCount','memoryCoordinateAudit', ...
            'sourceDomainAudit','symmetricLongitudinalReciprocalAudit'};
        if ~isstruct(leftPacket) || ~isscalar(leftPacket) || ...
                ~isstruct(rightPacket) || ~isscalar(rightPacket) || ...
                ~all(isfield(leftPacket,packetFields)) || ...
                ~all(isfield(rightPacket,packetFields))
            return
        end
        leftPackage = left.scheduledPackage;
        rightPackage = right.scheduledPackage;
        packageFields = {'mu','pointIds','weights','sourceContractId', ...
            'idx','parConst','internalCoupledCoordinate'};
        if ~isstruct(leftPackage) || ~isscalar(leftPackage) || ...
                ~isstruct(rightPackage) || ~isscalar(rightPackage) || ...
                ~all(isfield(leftPackage,packageFields)) || ...
                ~all(isfield(rightPackage,packageFields))
            return
        end
        tf = isequaln(left.elevatorIncrement,right.elevatorIncrement) && ...
            isequaln(left.thrustIncrement,right.thrustIncrement) && ...
            isequaln(left.rigidState,right.rigidState) && ...
            isequaln(left.ownerPolicy,right.ownerPolicy) && ...
            isequaln(left.scheduledReciprocal, ...
                right.scheduledReciprocal) && ...
            isequaln(leftPacket.schemaVersion, ...
                rightPacket.schemaVersion) && ...
            isequaln(leftPacket.queryPackageHash, ...
                rightPacket.queryPackageHash) && ...
            isequaln(leftPacket.bankManifestSha256, ...
                rightPacket.bankManifestSha256) && ...
            isequaln(leftPacket.sourceIds,rightPacket.sourceIds) && ...
            isequaln(leftPacket.schedule,rightPacket.schedule) && ...
            isequaln(leftPacket.queryTrim,rightPacket.queryTrim) && ...
            isequaln(leftPacket.hiddenStateCount, ...
                rightPacket.hiddenStateCount) && ...
            isequaln(leftPacket.memoryCoordinateAudit, ...
                rightPacket.memoryCoordinateAudit) && ...
            isequaln(leftPacket.sourceDomainAudit, ...
                rightPacket.sourceDomainAudit) && ...
            isequaln(leftPacket.symmetricLongitudinalReciprocalAudit, ...
                rightPacket.symmetricLongitudinalReciprocalAudit) && ...
            isequaln(leftPackage.mu,rightPackage.mu) && ...
            isequaln(leftPackage.pointIds,rightPackage.pointIds) && ...
            isequaln(leftPackage.weights,rightPackage.weights) && ...
            isequaln(leftPackage.sourceContractId, ...
                rightPackage.sourceContractId) && ...
            isequaln(leftPackage.idx,rightPackage.idx) && ...
            isequaln(leftPackage.parConst,rightPackage.parConst) && ...
            isequaln(leftPackage.internalCoupledCoordinate, ...
                rightPackage.internalCoupledCoordinate) && ...
            obj.Nc==10;
    end

    %------------------------------------------------------------------
    function idx = buildIndexMaps(obj)
        nx = obj.nx;
        nu = obj.nu;
        N  = obj.Nc;

        idx.x = cell(N+1,1);
        for k = 1:N+1
            idx.x{k} = (k-1)*nx + (1:nx);
        end

        startU = (N+1)*nx;

        idx.u = cell(N,1);
        for k = 1:N
            idx.u{k} = startU + (k-1)*nu + (1:nu);
        end
        if obj.terminalViabilityEnabled
            idx.terminalSlack = startU + N*nu + 1;
        else
            idx.terminalSlack = [];
        end
    end

    %------------------------------------------------------------------
    function zShift = shiftGuess(obj,z)
    % Shift optimal solution one sample forward:
    %
    %   X_0^{guess,new} <- X_1^*
    %   ...
    %   X_{N-1}^{guess,new} <- X_N^*
    %   X_N^{guess,new} <- X_N^*
    %
    %   U_0^{guess,new} <- U_1^*
    %   ...
    %   U_{N-2}^{guess,new} <- U_{N-1}^*
    %   U_{N-1}^{guess,new} <- U_{N-1}^*
    %
    % The next computeControl() call then overwrites X_0 with the newest
    % xhat and repairs the state guess by forward prediction.
    %------------------------------------------------------------------
        idx = obj.buildIndexMaps();
        N   = obj.Nc;

        zShift = zeros(size(z));

        for k = 1:N
            zShift(idx.x{k}) = z(idx.x{k+1});
        end
        zShift(idx.x{N+1}) = z(idx.x{N+1});

        for k = 1:N-1
            zShift(idx.u{k}) = z(idx.u{k+1});
        end
        zShift(idx.u{N}) = z(idx.u{N});
        if obj.terminalViabilityEnabled
            zShift(idx.terminalSlack) = z(idx.terminalSlack);
        end
    end

    %------------------------------------------------------------------
    function z = repairRealtimeRtiControlGuess(obj,z)
    %REPAIRREALTIMERTICONTROLGUESS Make the shifted endpoint/rate seed exact.
        assert(obj.realtimeRtiEnabled && numel(z)==numel(obj.z0), ...
            'nMPC:RealtimeRtiControlRepair', ...
            'The RTI control repair requires its enabled audit contract.');
        idx = obj.buildIndexMaps();
        control = reshape(z(idx.u{1}(1):idx.u{end}(end)), ...
            obj.nu,obj.Nc);
        endpoints = control(1:obj.n_surf,:);
        priorEndpoint = obj.uPrev(1:obj.n_surf);
        positionLower = obj.uL(1:obj.n_surf);
        positionUpper = obj.uU(1:obj.n_surf);
        rateLower = obj.uL(obj.n_surf+(1:obj.n_surf));
        rateUpper = obj.uU(obj.n_surf+(1:obj.n_surf));
        for interval = 1:obj.Nc
            reachableLower = max(positionLower, ...
                priorEndpoint+obj.Ts*rateLower);
            reachableUpper = min(positionUpper, ...
                priorEndpoint+obj.Ts*rateUpper);
            assert(all(reachableLower<=reachableUpper), ...
                'nMPC:RealtimeRtiControlRepairReachability', ...
                ['The applied endpoint has no position- and rate-feasible ', ...
                 'successor at interval %d.'],interval);
            if obj.symmetricSurfaceSubspaceEnabled
                symmetricLower = max(reachableLower);
                symmetricUpper = min(reachableUpper);
                assert(symmetricLower<=symmetricUpper, ...
                    'nMPC:RealtimeRtiControlRepairSymmetricReachability', ...
                    ['The applied endpoint has no symmetric position- and ', ...
                     'rate-feasible successor at interval %d.'],interval);
                symmetricEndpoint = mean(endpoints(:,interval));
                symmetricEndpoint = min(max(symmetricEndpoint, ...
                    symmetricLower),symmetricUpper);
                endpoints(:,interval) = symmetricEndpoint;
            else
                endpoints(:,interval) = min(max( ...
                    endpoints(:,interval),reachableLower),reachableUpper);
            end
            priorEndpoint = endpoints(:,interval);
        end
        priorEndpoint = obj.uPrev(1:obj.n_surf);
        rates = (endpoints-[priorEndpoint,endpoints(:,1:end-1)])/obj.Ts;
        control = [endpoints;rates];
        lower = repmat(obj.uL,1,obj.Nc);
        upper = repmat(obj.uU,1,obj.Nc);
        boundViolation = max([lower(:)-control(:); ...
            control(:)-upper(:);0]);
        assert(boundViolation<=1e-12, ...
            'nMPC:RealtimeRtiControlRepairBounds', ...
            ['The repaired endpoint/rate warm start exceeds a physical ', ...
             'bound by %.3e.'],boundViolation);
        rateResidual = endpoints-[priorEndpoint,endpoints(:,1:end-1)]- ...
            obj.Ts*rates;
        assert(norm(rateResidual(:),inf)<=1e-14, ...
            'nMPC:RealtimeRtiControlRepairEquality', ...
            'The repaired endpoint/rate warm start is not exact.');
        if obj.symmetricSurfaceSubspaceEnabled
            assert(norm(control(1,:)-control(2,:),inf)<=1e-14 && ...
                norm(control(3,:)-control(4,:),inf)<=1e-14, ...
                'nMPC:RealtimeRtiControlRepairSymmetry', ...
                'The Case-A RTI warm start is not symmetric.');
        end
        z(idx.u{1}(1):idx.u{end}(end)) = control(:);
    end

    %------------------------------------------------------------------
    function z = resetRealtimeRtiControlHorizon(obj,z)
    %RESETREALTIMERTICONTROLHORIZON Retain a finite trim-relative fallback.
        idx = obj.buildIndexMaps();
        for interval = 1:obj.Nc
            z(idx.u{interval}) = obj.uTrim;
        end
    end

    %------------------------------------------------------------------
    function [candidate,info] = buildRealtimeRtiSeed( ...
            obj,nlp,zInitial,wHorz)
    %BUILDREALTIMERTISEED Construct an exact condensed RTI full-NLP seed.
        assert(obj.realtimeRtiEnabled && ~isempty(obj.realtimeRtiSolver), ...
            'nMPC:RealtimeRtiDisabled', ...
            'The condensed RTI seed requires its approved audit selector.');
        if obj.nativeReducedHorizonRtiActive
            [candidate,info] = obj.buildNativeReducedHorizonRtiSeed( ...
                nlp,zInitial,wHorz);
            return
        end
        transformation = obj.buildRealtimeRtiTransformation();
        stateCoordinateCount = obj.nx*(obj.Nc+1);
        request = struct('auditOnly',true, ...
            'transformation',transformation, ...
            'eliminatedColumns',1:stateCoordinateCount, ...
            'eliminatedRows',1:stateCoordinateCount, ...
            'iterationCount',obj.realtimeRtiIterationCount, ...
            'constraintTolerance',1e-8, ...
            'inequalityTolerance',1e-10);
        [candidate,solverInfo] = obj.realtimeRtiSolver.solveCondensedRti( ...
            nlp.cost,zInitial,nlp.lb,nlp.ub,nlp.nonl,nlp.nonl,request);
        info = solverInfo;
        info.enabled = true;
        info.attempted = true;
        info.fallbackToFullInitial = false;
        info.changeId = obj.realtimeRtiChangeId;
        info.chartMode = obj.realtimeRtiChartMode;
        if obj.symmetricSurfaceSubspaceEnabled
            info.activeChart = obj.symmetricSurfaceSubspaceCaseId + ...
                "_symmetric_rate";
        else
            info.activeChart = "general_symmetric_differential_rate";
        end
        info.condensingMode = "full_state";
        info.correctionSolver = obj.solverName;
        info.nativeReducedHorizon = false;
        info.nativeKernelIdentity = "disabled";
        info.horizonPreparationSeconds = 0;
        info.causalWarmStartCondensationApplied = ...
            obj.nativeCausalRolloutActive;
    end

    %------------------------------------------------------------------
    function [candidate,info] = buildNativeReducedHorizonRtiSeed( ...
            obj,nlp,zInitial,wHorz)
    %BUILDNATIVEREDUCEDHORIZONRTISEED Apply exact causal Nc=10 condensing.
        totalTimer = tic;
        preparationTimer = tic;
        if obj.preparedHorizonDataReuseActive
            assert(isfield(nlp,'valueHorizonData'), ...
                'nMPC:PreparedHorizonDataReuseMissing', ...
                'The shared controller horizon data is unavailable.');
            data = nlp.valueHorizonData;
        else
            data = obj.buildNativeControllerValueHorizonData();
        end
        wingExpansion = obj.nativeReducedControllerWingExpansion();
        data.wingExpansion = wingExpansion;
        data.initialStateTarget = zInitial(1:obj.nx);
        data.disturbance = wHorz;
        horizonPreparationSeconds = toc(preparationTimer);
        packetFunction = @(z)obj.nativeReducedControllerPacket(z,data);
        request = struct('auditOnly',true, ...
            'changeId', ...
                "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1", ...
            'iterationCount',obj.realtimeRtiIterationCount, ...
            'constraintTolerance',1e-8, ...
            'inequalityTolerance',1e-10, ...
            'nonlinearBacktracking',struct('enabled',false));
        [candidate,solverInfo] = ...
            obj.realtimeRtiSolver.solvePreparedReducedRti( ...
                nlp.cost,zInitial,nlp.lb,nlp.ub,nlp.nonl, ...
                packetFunction,request);
        solverInfo.totalSeconds = toc(totalTimer);
        solverInfo.horizonPreparationSeconds = horizonPreparationSeconds;
        info = solverInfo;
        info.enabled = true;
        info.attempted = true;
        info.fallbackToFullInitial = false;
        info.changeId = ...
            "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1";
        info.chartMode = obj.realtimeRtiChartMode;
        info.activeChart = ...
            obj.symmetricSurfaceSubspaceCaseId+"_symmetric_rate";
        info.condensingMode = "native_reduced_horizon";
        info.correctionSolver = obj.solverName;
        info.nativeReducedHorizon = true;
        info.nativeKernelIdentity = ...
            obj.nativeReducedHorizonRtiKernelIdentity;
        info.causalWarmStartCondensationApplied = ...
            obj.nativeCausalRolloutActive;
        info.preparedHorizonDataReuseApplied = ...
            obj.preparedHorizonDataReuseActive;
        info.horizonDataBuildsPerSample = ...
            double(obj.preparedHorizonDataReuseActive);
    end

    %------------------------------------------------------------------
    function packet = nativeReducedControllerPacket(obj,z,data)
    %NATIVEREDUCEDCONTROLLERPACKET Assemble the exact symmetric rate chart.
        idx = obj.buildIndexMaps();
        nativeNodes = reshape(z(1:obj.nx*(obj.Nc+1)), ...
            obj.nx,obj.Nc+1);
        wingTotal = zeros(4,14,obj.Nc);
        wingIncrement = zeros(4,14,obj.Nc);
        for interval = 1:obj.Nc
            control = z(idx.u{interval});
            if interval==1
                previous = obj.uPrev(1:obj.n_surf);
            else
                previousControl = z(idx.u{interval-1});
                previous = previousControl(1:obj.n_surf);
            end
            uModel = obj.buildIntervalModelControl( ...
                [previous;zeros(obj.n_surf,1)],control);
            increment = uModel-obj.uModelTrim;
            wingIncrement(:,:,interval) = repmat(increment,1,14);
            wingTotal(:,:,interval) = repmat( ...
                data.packets(interval).queryTrimWing+increment,1,14);
            rateResidual = control(1:obj.n_surf)-previous- ...
                obj.Ts*control(obj.n_surf+(1:obj.n_surf));
            assert(norm(rateResidual,inf)<=1e-12 && ...
                abs(control(1)-control(2))<=1e-12, ...
                'nMPC:NativeReducedHorizonRateChart', ...
                'The native reduced horizon requires an exact symmetric rate seed.');
        end
        disturbance = reshape(data.disturbance,1,obj.Nc);
        [~,~,stateExpansion,stateOffset] = ...
            obj.nativeReducedHorizonRtiKernel( ...
                nativeNodes,data.latentInitial,data.initialStateTarget, ...
                disturbance,wingTotal,wingIncrement,data.elevator, ...
                data.thrust,data.rigid,data.wingExpansion,data.packets);

        terminalCount = double(obj.terminalViabilityEnabled);
        reducedCount = obj.Nc+terminalCount;
        expansion = spalloc(numel(z),reducedCount, ...
            obj.nx*obj.Nc*(obj.Nc+1)+5*obj.Nc+terminalCount);
        offset = zeros(numel(z),1);
        for node = 1:obj.Nc+1
            expansion(idx.x{node},1:obj.Nc) = stateExpansion(:,:,node);
            offset(idx.x{node}) = stateOffset(:,node);
        end
        surfaceMap = obj.symmetricSurfaceSubspaceMap;
        for interval = 1:obj.Nc
            rows = idx.u{interval};
            for rate = 1:interval
                expansion(rows(1:obj.n_surf),rate) = obj.Ts*surfaceMap;
            end
            expansion(rows(obj.n_surf+(1:obj.n_surf)),interval) = ...
                surfaceMap;
        end
        if terminalCount>0
            expansion(idx.terminalSlack,end) = 1;
        end

        if isfinite(obj.aT)
            terminalState = nativeNodes(:,end);
            terminalError = terminalState-obj.xTrim;
            terminalValue = 0.5*terminalError.'*obj.Pc*terminalError;
            inequality = terminalValue-obj.aT;
            gradient = zeros(numel(z),1);
            gradient(idx.x{end}) = obj.Pc*terminalError;
            if terminalCount>0
                inequality = inequality-z(idx.terminalSlack);
                gradient(idx.terminalSlack) = -1;
            end
            A = gradient.'*expansion;
            b = -inequality-gradient.'*offset;
        else
            A = sparse(0,reducedCount);
            b = zeros(0,1);
        end
        packet = struct('expansion',expansion,'offset',offset, ...
            'A',sparse(A),'b',b,'Aeq',sparse(0,reducedCount), ...
            'beq',zeros(0,1));
    end

    %------------------------------------------------------------------
    function expansion = nativeReducedControllerWingExpansion(obj)
    %NATIVEREDUCEDCONTROLLERWINGEXPANSION Map symmetric rates to ROM input.
        expansion = zeros(4,obj.Nc,obj.Nc);
        surfaceMap = obj.symmetricSurfaceSubspaceMap;
        for interval = 1:obj.Nc
            for rate = 1:interval-1
                expansion(1:obj.n_surf,rate,interval) = ...
                    obj.Ts*surfaceMap;
            end
            expansion(1:obj.n_surf,interval,interval) = ...
                obj.actuatorDeflectionAlpha*obj.Ts*surfaceMap;
            expansion(obj.n_surf+(1:obj.n_surf),interval,interval) = ...
                surfaceMap;
        end
    end

    %------------------------------------------------------------------
    function transformation = buildRealtimeRtiTransformation(obj)
    %BUILDREALTIMERTITRANSFORMATION Map reduced rate steps to full decisions.
        stateCount = obj.nx*(obj.Nc+1);
        terminalCount = double(obj.terminalViabilityEnabled);
        if obj.symmetricSurfaceSubspaceEnabled
            coordinateCount = 1;
            endpointMap = [1;1];
        else
            coordinateCount = 2;
            endpointMap = [1 1;1 -1];
        end
        reducedCount = stateCount+coordinateCount*obj.Nc+terminalCount;
        transformation = spalloc(numel(obj.z0),reducedCount, ...
            stateCount+8*coordinateCount*obj.Nc+terminalCount);
        transformation(1:stateCount,1:stateCount) = speye(stateCount);
        for interval = 1:obj.Nc
            fullRows = stateCount+(interval-1)*obj.nu+(1:obj.nu);
            current = stateCount+(interval-1)*coordinateCount+ ...
                (1:coordinateCount);
            transformation(fullRows(obj.n_surf+(1:obj.n_surf)),current) = ...
                endpointMap;
            for prior = 1:interval
                rateColumn = stateCount+(prior-1)*coordinateCount+ ...
                    (1:coordinateCount);
                transformation(fullRows(1:obj.n_surf),rateColumn) = ...
                    obj.Ts*endpointMap;
            end
        end
        if terminalCount>0
            transformation(end,end) = 1;
        end
        assert(size(transformation,1)==numel(obj.z0) && ...
            size(transformation,2)==reducedCount && ...
            all(sum(abs(transformation),1)>0), ...
            'nMPC:RealtimeRtiTransformationDimension', ...
            'The condensed RTI rate transformation is incomplete.');
    end

    %------------------------------------------------------------------
    function info = emptyRealtimeRtiInfo(obj)
        info = struct('enabled',obj.realtimeRtiEnabled, ...
            'attempted',false,'qualified',false, ...
            'fallbackToFullInitial',false,'message',"not attempted", ...
            'identifier',"",'changeId',obj.realtimeRtiChangeId, ...
            'chartMode',obj.realtimeRtiChartMode, ...
            'activeChart',"disabled",'condensingMode',"disabled", ...
            'correctionSolver',obj.solverName,'iterationCount',0, ...
            'completedIterations',0,'iterations',struct([]), ...
            'nonlinearEqualityInf',nan,'nonlinearInequalityMax',nan, ...
            'boundViolationMax',nan,'constraintViolationInf',nan, ...
            'objective',nan,'valueReplaySeconds',0,'totalSeconds',0, ...
            'nativeReducedHorizon',false, ...
            'nativeKernelIdentity',"disabled", ...
            'horizonPreparationSeconds',0, ...
            'acceptedLatentHorizonCondensationApplied',false, ...
            'causalWarmStartCondensationApplied', ...
                obj.nativeCausalRolloutActive);
        info.preparedHorizonDataReuseApplied = ...
            obj.preparedHorizonDataReuseActive;
        info.horizonDataBuildsPerSample = ...
            double(obj.preparedHorizonDataReuseActive);
    end

    %------------------------------------------------------------------
    function z = repairStateGuess(obj,z,x0,wHorz,valueHorizonData)
    % Forward-predict state nodes from the newest measured/estimated x0
    % using the same interval-input map as the continuity constraints.
    %------------------------------------------------------------------
        idx  = obj.buildIndexMaps();
        Nint = round(obj.Ts/obj.dt);

        x = x0(:);
        z(idx.x{1}) = x;
        if obj.reciprocalProviderEnabled
            assert(numel(obj.reciprocalInitialLatentState)== ...
                obj.reciprocalProvider.hiddenStateCount && ...
                numel(obj.reciprocalFutureContextHorizon)==obj.Nc, ...
                'nMPC:ReciprocalWarmStartContext', ...
                'The reciprocal nMPC warm start lacks causal context.');
            augmentedState = [x;obj.reciprocalInitialLatentState];
        end

        if obj.reciprocalProviderEnabled && obj.nativeCausalRolloutActive
            if nargin < 5 || isempty(fieldnames(valueHorizonData))
                data = obj.buildNativeControllerValueHorizonData();
            else
                data = valueHorizonData;
            end
            endpoints = obj.evaluateNativeControllerCausalRollout( ...
                z,wHorz,data,x);
            for k = 1:obj.Nc
                z(idx.x{k+1}) = endpoints(1:obj.nativeStateCount,k);
            end
            return
        end

        for k = 1:obj.Nc
            uCurrent = z(idx.u{k});
            if k == 1
                uPrevious = obj.uPrev;
            else
                uPrevious = z(idx.u{k-1});
            end
            uModel = obj.buildIntervalModelControl(uPrevious,uCurrent);
            wk = wHorz(:,k);

            if obj.reciprocalProviderEnabled
                wingIncrement = uModel-obj.uModelTrim;
                try
                    augmentedState = ...
                        obj.reciprocalProvider.propagateControlInterval( ...
                            augmentedState,wk,wingIncrement, ...
                            obj.reciprocalFutureContextHorizon{k},false);
                catch propagationException
                    kernel = obj.reciprocalProvider. ...
                        fixedIntervalKernelAuditSnapshot();
                    error('nMPC:ReciprocalWarmStartPropagation', ...
                        ['Reciprocal warm-start propagation failed at ', ...
                         'horizon interval %d (state inf %.6g, latent inf ', ...
                         '%.6g, wing increment inf %.6g, disturbance inf ', ...
                         '%.6g). Compiled-kernel fallback: %s'], ...
                        k,norm(augmentedState,inf), ...
                        norm(augmentedState(obj.reciprocalLatentIndex),inf), ...
                        norm(wingIncrement,inf),norm(wk,inf), ...
                        char(kernel.lastFallback));
                end
                x = augmentedState(1:obj.nativeStateCount);
            else
                for m = 1:Nint
                    [x,~] = obj.model.step(x,uModel,wk,[],false);
                end
            end

            z(idx.x{k+1}) = x;
        end
    end

    %------------------------------------------------------------------
    function latentHorizon = rebuildReciprocalLatentHorizon( ...
            obj,X,U,W,valueHorizonData,acceptedReplayEndpoints)
    %REBUILDRECIPROCALLATENTHORIZON Recover accepted prediction memory.
        assert(obj.reciprocalProviderEnabled && ...
            isequal(size(X),[obj.nativeStateCount,obj.Nc+1]) && ...
            isequal(size(U),[obj.nu,obj.Nc]) && ...
            isequal(size(W),[obj.nw,obj.Nc]), ...
            'nMPC:ReciprocalLatentRebuild', ...
            'The accepted reciprocal nMPC horizon is incompatible.');
        latentHorizon = zeros( ...
            obj.reciprocalProvider.hiddenStateCount,obj.Nc+1);
        latentHorizon(:,1) = obj.reciprocalInitialLatentState;
        if obj.acceptedLatentHorizonCondensationActive
            if nargin >= 6 && ~isempty(acceptedReplayEndpoints)
                endpoints = acceptedReplayEndpoints;
            else
                if nargin < 5 || isempty(fieldnames(valueHorizonData))
                    data = obj.buildNativeControllerValueHorizonData();
                else
                    data = valueHorizonData;
                end
                idx = obj.buildIndexMaps();
                decision = zeros(obj.decisionVariableCount(),1);
                decision(1:(obj.Nc+1)*obj.nx) = X(:);
                for interval = 1:obj.Nc
                    decision(idx.u{interval}) = U(:,interval);
                end
                endpoints = obj.evaluateNativeControllerValueHorizon( ...
                    decision,W,data);
            end
            assert(isequal(size(endpoints), ...
                [obj.reciprocalProvider.stateCount,obj.Nc]) && ...
                all(isfinite(endpoints),'all'), ...
                'nMPC:AcceptedReplayReuseOutput', ...
                'The accepted controller replay endpoints are invalid.');
            latentHorizon(:,2:obj.Nc+1) = ...
                endpoints(obj.reciprocalLatentIndex,:);
            assert(all(isfinite(latentHorizon),'all'), ...
                'nMPC:ReciprocalLatentRebuild', ...
                'The accepted reciprocal latent prediction is nonfinite.');
            return
        end
        for interval = 1:obj.Nc
            if interval == 1
                uPrevious = obj.uPrev;
            else
                uPrevious = U(:,interval-1);
            end
            uModel = obj.buildIntervalModelControl( ...
                uPrevious,U(:,interval));
            augmented = ...
                obj.reciprocalProvider.propagateControlInterval( ...
                    [X(:,interval);latentHorizon(:,interval)], ...
                    W(:,interval),uModel-obj.uModelTrim, ...
                    obj.reciprocalFutureContextHorizon{interval},false);
            latentHorizon(:,interval+1) = ...
                augmented(obj.reciprocalLatentIndex);
        end
        assert(all(isfinite(latentHorizon),'all'), ...
            'nMPC:ReciprocalLatentRebuild', ...
            'The accepted reciprocal latent prediction is nonfinite.');
    end

    %------------------------------------------------------------------
    function uModel = buildIntervalModelControl(obj,uPrevious,uCurrent)
    % Convert trim-relative endpoint decisions to the total interval input.
        uPrevious = obj.expandToLength(uPrevious,obj.nu,'uPrevious');
        uCurrent = obj.expandToLength(uCurrent,obj.nu,'uCurrent');
        nSurf = obj.n_surf;

        deltaModel = (1-obj.actuatorDeflectionAlpha) * ...
            uPrevious(1:nSurf) + obj.actuatorDeflectionAlpha * ...
            uCurrent(1:nSurf);
        rateModel = uCurrent(nSurf+1:2*nSurf);
        uModel = obj.uModelTrim + [deltaModel; rateModel];
    end

    %------------------------------------------------------------------
    function wPred = buildPredictedGust(obj)
    % Disturbance prediction used by MPC:
    %   - current MHE gust estimate at the first MPC interval,
    %   - linearly decays over the first Nd intervals,
    %   - zero afterward.
    %
    % Output:
    %   wPred: nw x Nc
    %------------------------------------------------------------------
        wPred = zeros(obj.nw,obj.Nc);

        w0 = obj.wHat(:);

        if isempty(w0) || norm(w0) == 0
            return
        end

        NdEff = max(1,min(obj.Nd,obj.Nc));

        % alpha = [1 ... 0] over the decay portion.
        alpha = linspace(1,0,NdEff);

        wPred(:,1:NdEff) = w0 * alpha;
    end

    %------------------------------------------------------------------
    function w = parseDisturbanceEstimate(obj,whatEst)
    % Accepts either:
    %   1) current disturbance estimate, nw x 1, or
    %   2) full MHE horizon vector, nw*Ne x 1.
    %
    % Uses the last nw entries as the current disturbance estimate.
    %------------------------------------------------------------------
        if isempty(whatEst)
            w = zeros(obj.nw,1);
            return
        end

        wVec = whatEst(:);

        if numel(wVec) < obj.nw
            error('nMPC:Dimension', ...
                'whatEst has fewer entries than nw.');
        end

        w = wVec(end-obj.nw+1:end);
    end

    %------------------------------------------------------------------
    function priority = parsePriorityReference(obj,value)
    %PARSEPRIORITYREFERENCE Validate the disabled-by-default audit contract.
        assert(isstruct(value) && isfield(value,'enabled') && ...
            isscalar(value.enabled), 'nMPC:PriorityReference', ...
            'priorityReference requires scalar logical enabled.');
        priority = struct('enabled',logical(value.enabled), ...
            'referenceHorizon',repmat(obj.uTrim,1,obj.Nc), ...
            'primaryCostTolerance',0);
        if ~priority.enabled
            return
        end
        required = {'auditOnly','referenceHorizon','primaryCostTolerance'};
        for field = required
            assert(isfield(value,field{1}), 'nMPC:PriorityReference', ...
                'Enabled priorityReference requires %s.',field{1});
        end
        assert(logical(value.auditOnly), 'nMPC:PriorityReference', ...
            'The priority reference is audit-only.');
        reference = value.referenceHorizon;
        assert(isnumeric(reference) && isequal(size(reference),[obj.nu,obj.Nc]) && ...
            all(isfinite(reference),'all'), 'nMPC:PriorityReference', ...
            'referenceHorizon must be finite nu-by-Nc.');
        tolerance = value.primaryCostTolerance;
        assert(isscalar(tolerance) && isfinite(tolerance) && tolerance >= 0, ...
            'nMPC:PriorityReference', ...
            'primaryCostTolerance must be finite and nonnegative.');
        priority.referenceHorizon = reference;
        priority.primaryCostTolerance = tolerance;
    end

    %------------------------------------------------------------------
    function nlp = buildPriorityNlp(obj,primaryNlp,referenceHorizon,primaryCostGuard)
    %BUILDPRIORITYNLP Select an outer-reference-near primary-equivalent path.
        assert(isequal(size(referenceHorizon),[obj.nu,obj.Nc]) && ...
            all(isfinite(referenceHorizon),'all') && ...
            isscalar(primaryCostGuard) && isfinite(primaryCostGuard), ...
            'nMPC:PriorityReference', ...
            'Priority NLP inputs are invalid.');
        idx = obj.buildIndexMaps();
        nVar = obj.decisionVariableCount();
        nSurf = obj.n_surf;
        deflectionScale = max(abs([obj.uL(1:nSurf),obj.uU(1:nSurf)]),[],2);
        assert(all(isfinite(deflectionScale)) && all(deflectionScale > 0), ...
            'nMPC:PriorityReference', ...
            'The priority-reference deflection normalization is invalid.');
        Nc = obj.Nc;
        nlp = primaryNlp;
        nlp.cost = @secondaryCost;
        nlp.nonl = @priorityNonlinear;

        function [J,g] = secondaryCost(z)
            J = 0;
            g = zeros(nVar,1);
            for k = 1:Nc
                du = z(idx.u{k}(1:nSurf)) - referenceHorizon(1:nSurf,k);
                normalizedError = du ./ deflectionScale;
                J = J + 0.5*(normalizedError.'*normalizedError);
                g(idx.u{k}(1:nSurf)) = g(idx.u{k}(1:nSurf)) + ...
                    du ./ (deflectionScale.^2);
            end
        end

        function [c,ceq,gradc,gradceq] = priorityNonlinear(z)
            [cPrimary,ceq,gradPrimary,gradceq] = primaryNlp.nonl(z);
            [primaryCost,primaryGradient] = primaryNlp.cost(z);
            c = [cPrimary; primaryCost-primaryCostGuard];
            gradc = [gradPrimary,primaryGradient];
        end
    end

    %------------------------------------------------------------------
    function data = buildNativeControllerValueHorizonData(obj)
    %BUILDNATIVECONTROLLERVALUEHORIZONDATA Pack invariant E3 interval data.
        contexts = reshape(obj.reciprocalFutureContextHorizon,1,obj.Nc);
        obj.reciprocalProvider.configureHorizonSensitivityStencilAudit( ...
            contexts);
        packets = AeroFlex.ctrl. ...
            buildScheduledReciprocalHorizonPacketAudit(contexts);
        elevator = zeros(1,14,obj.Nc);
        thrust = zeros(1,14,obj.Nc);
        rigid = zeros(9,14,obj.Nc);
        for interval = 1:obj.Nc
            context = contexts{interval};
            elevator(:,:,interval) = context.elevatorIncrement;
            thrust(:,:,interval) = context.thrustIncrement;
            rigid(:,:,interval) = context.rigidState;
        end
        data = struct('packets',packets,'elevator',elevator, ...
            'thrust',thrust,'rigid',rigid, ...
            'latentInitial',obj.reciprocalInitialLatentState);
    end

    %------------------------------------------------------------------
    function [endpoints,rateResidual,symmetryResidual] = ...
            evaluateNativeControllerValueHorizon(obj,z,wHorz,data)
    %EVALUATENATIVECONTROLLERVALUEHORIZON Execute exact E3 value replay.
        idx = obj.buildIndexMaps();
        X = reshape(z(1:obj.nx*(obj.Nc+1)),obj.nx,obj.Nc+1);
        U = reshape(z([idx.u{:}]),obj.nu,obj.Nc);
        wingTotal = zeros(4,14,obj.Nc);
        wingIncrement = zeros(4,14,obj.Nc);
        rateResidual = zeros(obj.n_surf,obj.Nc);
        symmetryResidual = zeros(obj.Nc,1);
        for interval = 1:obj.Nc
            if interval==1
                previous = obj.uPrev;
            else
                previous = U(:,interval-1);
            end
            current = U(:,interval);
            uModel = obj.buildIntervalModelControl(previous,current);
            increment = uModel-obj.uModelTrim;
            wingIncrement(:,:,interval) = repmat(increment,1,14);
            wingTotal(:,:,interval) = repmat( ...
                data.packets(interval).queryTrimWing+increment,1,14);
            rateResidual(:,interval) = current(1:obj.n_surf)- ...
                previous(1:obj.n_surf)- ...
                obj.Ts*current(obj.n_surf+(1:obj.n_surf));
            symmetryResidual(interval) = current(1)-current(2);
        end
        endpoints = obj.nativeValueHorizonKernel( ...
            X(:,1:obj.Nc),data.latentInitial,wHorz,wingTotal, ...
            wingIncrement,data.elevator,data.thrust,data.rigid,data.packets);
        assert(isequal(size(endpoints),[obj.reciprocalProvider.stateCount, ...
            obj.Nc]) && all(isfinite(endpoints),'all'), ...
            'nMPC:NativeValueHorizonOutput', ...
            'The controller E3 value-horizon output is invalid.');
        rateResidual = rateResidual(:);
        if ~obj.symmetricSurfaceSubspaceEnabled
            symmetryResidual = zeros(0,1);
        end
    end

    %------------------------------------------------------------------
    function endpoints = evaluateNativeControllerCausalRollout( ...
            obj,z,wHorz,data,x0)
    %EVALUATENATIVECONTROLLERCAUSALROLLOUT Execute exact feasible rollout.
        assert(obj.nativeCausalRolloutActive && numel(x0)==obj.nx && ...
            isequal(size(wHorz),[obj.nw,obj.Nc]), ...
            'nMPC:NativeCausalRolloutInput', ...
            'The controller causal-rollout input is incompatible.');
        idx = obj.buildIndexMaps();
        U = reshape(z([idx.u{:}]),obj.nu,obj.Nc);
        wingTotal = zeros(4,14,obj.Nc);
        wingIncrement = zeros(4,14,obj.Nc);
        for interval = 1:obj.Nc
            if interval==1
                previous = obj.uPrev;
            else
                previous = U(:,interval-1);
            end
            uModel = obj.buildIntervalModelControl(previous,U(:,interval));
            increment = uModel-obj.uModelTrim;
            wingIncrement(:,:,interval) = repmat(increment,1,14);
            wingTotal(:,:,interval) = repmat( ...
                data.packets(interval).queryTrimWing+increment,1,14);
        end
        endpoints = obj.nativeCausalRolloutKernel( ...
            x0,data.latentInitial,wHorz,wingTotal,wingIncrement, ...
            data.elevator,data.thrust,data.rigid,data.packets);
        assert(isequal(size(endpoints),[obj.reciprocalProvider.stateCount, ...
            obj.Nc]) && all(isfinite(endpoints),'all'), ...
            'nMPC:NativeCausalRolloutOutput', ...
            'The controller causal-rollout output is invalid.');
    end

    %------------------------------------------------------------------
    function nlp = assembleWindow( ...
            obj,x0Current,wHorz,uBaseHorizon,valueHorizonDataOverride)
    % Build cost, nonlinear constraints, and bounds for fmincon.
    %------------------------------------------------------------------
        nx  = obj.nx;
        nu  = obj.nu;
        nw  = obj.nw;
        Nc  = obj.Nc;
        Nint = round(obj.Ts/obj.dt);

        idx = obj.buildIndexMaps();

        nVar = obj.decisionVariableCount();

        Q = obj.Qc; R = obj.Rc; P = obj.Pc;
        xRef = obj.xTrim;
        assert(isequal(size(uBaseHorizon),[nu,Nc]) && ...
            all(isfinite(uBaseHorizon),'all'), 'nMPC:KnownBaseHorizon', ...
            'uBaseHorizon must be a finite nu-by-Nc command matrix.');
        valueHorizonData = struct();
        if obj.nativeValueHorizonActive
            if nargin >= 5 && ~isempty(fieldnames(valueHorizonDataOverride))
                valueHorizonData = valueHorizonDataOverride;
            else
                valueHorizonData = ...
                    obj.buildNativeControllerValueHorizonData();
            end
        end
        valueReplayCacheDecision = zeros(0,1);
        valueReplayCacheInequality = zeros(0,1);
        valueReplayCacheEquality = zeros(0,1);
        valueReplayCacheEndpoints = zeros(0,0);
        valueReplayCacheHits = 0;
        valueReplayCacheMisses = 0;

        % --------------------- COST ----------------------------------
        function [J,g,H] = costFun(z)
            J = 0;
            g = zeros(nVar,1);
            if nargout > 2
                H = sparse(nVar,nVar);
            end

            % Stage cost: k = 0,...,Nc-1
            for k = 0:Nc-1
                xk = z(idx.x{k+1});
                uk = z(idx.u{k+1});

                dX = xk - xRef;
                dU = uk - uBaseHorizon(:,k+1);

                J = J + 0.5*dX.'*Q*dX + 0.5*dU.'*R*dU;

                g(idx.x{k+1}) = g(idx.x{k+1}) + Q*dX;
                g(idx.u{k+1}) = g(idx.u{k+1}) + R*dU;
                if nargout > 2
                    H(idx.x{k+1},idx.x{k+1}) = ...
                        H(idx.x{k+1},idx.x{k+1})+Q;
                    H(idx.u{k+1},idx.u{k+1}) = ...
                        H(idx.u{k+1},idx.u{k+1})+R;
                end
                if obj.wingtipOutputCostEnabled
                    tipDeviation = obj.wingtipOutputGradient*dX;
                    J = J + 0.5*obj.wingtipOutputCostStageWeight* ...
                        tipDeviation^2;
                    g(idx.x{k+1}) = g(idx.x{k+1}) + ...
                        obj.wingtipOutputCostStageWeight* ...
                        obj.wingtipOutputGradient.'*tipDeviation;
                    if nargout > 2
                        H(idx.x{k+1},idx.x{k+1}) = ...
                            H(idx.x{k+1},idx.x{k+1})+ ...
                            obj.wingtipOutputCostStageWeight* ...
                            (obj.wingtipOutputGradient.'* ...
                            obj.wingtipOutputGradient);
                    end
                end
            end

            % Terminal cost: X_Nc only
            xN = z(idx.x{Nc+1});
            dXT = xN - xRef;

            J = J + 0.5*dXT.'*P*dXT;
            g(idx.x{Nc+1}) = g(idx.x{Nc+1}) + P*dXT;
            if nargout > 2
                H(idx.x{Nc+1},idx.x{Nc+1}) = ...
                    H(idx.x{Nc+1},idx.x{Nc+1})+P;
            end
            if obj.wingtipOutputCostEnabled
                terminalTipDeviation = obj.wingtipOutputGradient*dXT;
                J = J + 0.5*obj.wingtipOutputCostTerminalWeight* ...
                    terminalTipDeviation^2;
                g(idx.x{Nc+1}) = g(idx.x{Nc+1}) + ...
                    obj.wingtipOutputCostTerminalWeight* ...
                    obj.wingtipOutputGradient.'*terminalTipDeviation;
                if nargout > 2
                    H(idx.x{Nc+1},idx.x{Nc+1}) = ...
                        H(idx.x{Nc+1},idx.x{Nc+1})+ ...
                        obj.wingtipOutputCostTerminalWeight* ...
                        (obj.wingtipOutputGradient.'* ...
                        obj.wingtipOutputGradient);
                end
            end
            if obj.terminalViabilityEnabled
                terminalSlack = z(idx.terminalSlack);
                J = J + 0.5*obj.terminalViabilityPenalty*terminalSlack^2;
                g(idx.terminalSlack) = ...
                    obj.terminalViabilityPenalty*terminalSlack;
                if nargout > 2
                    H(idx.terminalSlack,idx.terminalSlack) = ...
                        obj.terminalViabilityPenalty;
                end
            end
        end

        % ----------------- NONLINEAR CONSTRAINTS ----------------------
        function [c,ceq,gradc,gradceq] = nonlFun(z)
            % Equality constraints enforce the initial state, shooting
            % continuity, and endpoint deflection/rate consistency.
            includeGradients = nargout > 2;
            currentReplayEndpoints = zeros(0,0);
            if obj.acceptedReplayReuseActive && ~includeGradients
                if isequal(z,valueReplayCacheDecision)
                    c = valueReplayCacheInequality;
                    ceq = valueReplayCacheEquality;
                    gradc = [];
                    gradceq = [];
                    valueReplayCacheHits = valueReplayCacheHits+1;
                    return
                end
                valueReplayCacheMisses = valueReplayCacheMisses+1;
            end
            nSurf = obj.n_surf;
            alphaDefl = obj.actuatorDeflectionAlpha;

            ceqRate = zeros(nSurf*Nc,1);
            if includeGradients
                JRateEq = spalloc(nSurf*Nc,nVar,3*nSurf*Nc);
            else
                JRateEq = sparse(0,nVar);
            end
            if obj.symmetricSurfaceSubspaceEnabled
                % Equal endpoint deflections plus the existing two physical
                % endpoint/rate equations imply equal rates when the prior
                % endpoint is symmetric.  Adding a separate rate-difference
                % row would make one equality dependent per interval.
                ceqSymmetry = zeros(Nc,1);
                if includeGradients
                    JSymmetry = spalloc(Nc,nVar,2*Nc);
                else
                    JSymmetry = sparse(0,nVar);
                end
            else
                ceqSymmetry = zeros(0,1);
                JSymmetry = sparse(0,nVar);
            end

            ceq = zeros(nx*(Nc+1),1);
            if includeGradients
                Jceq = spalloc(nx*(Nc+1),nVar, ...
                    nx + Nc*(nx + nx*nx + nx*nu + nx*nSurf));
            else
                Jceq = sparse(0,nVar);
            end

            % ---- initial state equality: X_0 - xhat = 0 --------------
            ceq(1:nx) = z(idx.x{1}) - x0Current;
            if includeGradients
                Jceq(1:nx,idx.x{1}) = eye(nx);
            end

            % ---- continuity constraints ------------------------------
            row0 = nx;
            row0rate = 0;
            if obj.reciprocalProviderEnabled
                assert(numel(obj.reciprocalInitialLatentState)== ...
                    obj.reciprocalProvider.hiddenStateCount && ...
                    numel(obj.reciprocalFutureContextHorizon)==Nc, ...
                    'nMPC:ReciprocalConstraintContext', ...
                    'The reciprocal nMPC constraint lacks future context.');
                if isa(obj.reciprocalProvider, ...
                        'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
                    obj.reciprocalProvider. ...
                        configureHorizonSensitivityStencilAudit( ...
                            obj.reciprocalFutureContextHorizon);
                end
                latent = obj.reciprocalInitialLatentState;
                if includeGradients
                    latentDerivative = sparse( ...
                        obj.reciprocalProvider.hiddenStateCount,nVar);
                else
                    latentDerivative = sparse(0,nVar);
                end
            end

            useNativeValueHorizon = obj.nativeValueHorizonActive && ...
                obj.reciprocalProviderEnabled && ~includeGradients;
            if useNativeValueHorizon
                [endpoints,ceqRate,ceqSymmetry] = ...
                    obj.evaluateNativeControllerValueHorizon( ...
                        z,wHorz,valueHorizonData);
                currentReplayEndpoints = endpoints;
                stateNodes = reshape(z(1:nx*(Nc+1)),nx,Nc+1);
                ceq(nx+(1:nx*Nc)) = reshape( ...
                    stateNodes(:,2:Nc+1)-endpoints(1:nx,:),[],1);
            else
            for k = 0:Nc-1
                xk = z(idx.x{k+1});
                uk = z(idx.u{k+1});

                wk = wHorz(:,k+1);

                x = xk;
                delta_k = uk(1:nSurf);
                rate_k  = uk(nSurf+1:2*nSurf);
            
                if k == 0
                    delta_prev = obj.uPrev(1:nSurf);
                    prevIsDecision = false;
                else
                    uPrevDecision = z(idx.u{k});
                    delta_prev = uPrevDecision(1:nSurf);
                    prevIsDecision = true;
                end
            
                rowsRate = row0rate + (1:nSurf);
                ceqRate(rowsRate) = delta_k - delta_prev - obj.Ts*rate_k;

                if includeGradients
                    % d/d delta_k
                    JRateEq(rowsRate,idx.u{k+1}(1:nSurf)) = ...
                        speye(nSurf);

                    % d/d rate_k
                    JRateEq(rowsRate, ...
                        idx.u{k+1}(nSurf+1:2*nSurf)) = ...
                        -obj.Ts*speye(nSurf);
                end

                if obj.symmetricSurfaceSubspaceEnabled
                    rowSymmetry = k+1;
                    ceqSymmetry(rowSymmetry) = delta_k(1)-delta_k(2);
                    if includeGradients
                        JSymmetry(rowSymmetry,idx.u{k+1}(1:2)) = [1 -1];
                    end
                end

                % d/d delta_{k-1}
                if includeGradients && prevIsDecision
                    JRateEq(rowsRate, idx.u{k}(1:nSurf)) = -speye(nSurf);
                end
            
                row0rate = row0rate + nSurf;

                 % Deflection used by the aerodynamic/ROM model over interval k.
                %
                % alphaDefl = 0.5 gives midpoint deflection:
                %
                %   delta_model = 0.5*(delta_start + delta_end)
                uModel = obj.buildIntervalModelControl( ...
                    [delta_prev; zeros(nSurf,1)],uk);
                if obj.reciprocalProviderEnabled
                    [augmentedEnd,localSensitivity] = ...
                        obj.reciprocalProvider.propagateControlInterval( ...
                            [x;latent],wk,uModel-obj.uModelTrim, ...
                            obj.reciprocalFutureContextHorizon{k+1}, ...
                            includeGradients);
                    xEnd = augmentedEnd(1:obj.nativeStateCount);
                    rows = row0+(1:nx);
                    ceq(rows) = z(idx.x{k+2})-xEnd;

                    if includeGradients
                        augmentedDerivative = sparse( ...
                            obj.reciprocalProvider.stateCount,nVar);
                        augmentedDerivative(1:nx,idx.x{k+1}) = speye(nx);
                        augmentedDerivative(obj.reciprocalLatentIndex,:) = ...
                            latentDerivative;
                        stateSensitivity = sparse(localSensitivity(:, ...
                            1:obj.reciprocalProvider.stateCount));
                        nextDerivative = ...
                            stateSensitivity*augmentedDerivative;
                        controlSensitivity = sparse(localSensitivity(:, ...
                            obj.reciprocalProvider.stateCount+obj.nw+ ...
                            (1:obj.nu)));
                        nextDerivative(:,idx.u{k+1}(1:nSurf)) = ...
                            nextDerivative(:,idx.u{k+1}(1:nSurf))+ ...
                            alphaDefl*controlSensitivity(:,1:nSurf);
                        nextDerivative(:, ...
                            idx.u{k+1}(nSurf+1:2*nSurf)) = ...
                            nextDerivative(:, ...
                            idx.u{k+1}(nSurf+1:2*nSurf))+ ...
                            controlSensitivity(:,nSurf+1:2*nSurf);
                        if prevIsDecision
                            nextDerivative(:,idx.u{k}(1:nSurf)) = ...
                                nextDerivative(:,idx.u{k}(1:nSurf))+ ...
                                (1-alphaDefl)*controlSensitivity(:,1:nSurf);
                        end
                        Jceq(rows,:) = -nextDerivative(1:nx,:);
                        Jceq(rows,idx.x{k+2}) = ...
                            Jceq(rows,idx.x{k+2})+speye(nx);
                        latentDerivative = nextDerivative( ...
                            obj.reciprocalLatentIndex,:);
                    end
                    latent = augmentedEnd(obj.reciprocalLatentIndex);
                    row0 = row0+nx;
                    continue
                end
                if includeGradients
                    % S = d x_current / d [xk; wk; uk].
                    S = [eye(nx), zeros(nx,nw+nu)];
                    for m = 1:Nint
                        [x,S] = obj.model.step(x,uModel,wk,S,true);
                    end
                    Sx = S(:,1:nx);
                    Su = S(:,nx+nw+1:end);
                    Su_delta = Su(:,1:nSurf);
                    Su_rate  = Su(:,nSurf+1:2*nSurf);
                else
                    for m = 1:Nint
                        x = obj.model.step(x,uModel,wk,[],false);
                    end
                end

                xEnd = x;

                rows = row0 + (1:nx);

                % Defect convention:
                %   X_{k+1} - Phi(X_k,U_k,W_k) = 0
                ceq(rows) = z(idx.x{k+2}) - xEnd;

                if includeGradients
                    Jceq(rows,idx.x{k+2}) =  speye(nx);
                    Jceq(rows,idx.x{k+1}) = -(Sx);
                    Jceq(rows,idx.u{k+1}(1:nSurf)) = ...
                        Jceq(rows,idx.u{k+1}(1:nSurf))- ...
                        alphaDefl*Su_delta;

                    Jceq(rows,idx.u{k+1}(nSurf+1:2*nSurf)) = ...
                        Jceq(rows,idx.u{k+1}(nSurf+1:2*nSurf))-Su_rate;

                    if prevIsDecision
                        Jceq(rows,idx.u{k}(1:nSurf)) = ...
                            Jceq(rows,idx.u{k}(1:nSurf))- ...
                            (1-alphaDefl)*Su_delta;
                    end
                end

                row0 = row0 + nx;
            end
            end
            ceq = [ceq; ceqRate; ceqSymmetry];

            if includeGradients
                % fmincon wants gradceq as nVar x nEq.
                gradceq = [Jceq; JRateEq; JSymmetry].';
            else
                gradceq = [];
            end

            % ---- terminal set inequality -----------------------------
            xN = z(idx.x{Nc+1});
            eN = xN - xRef;
            terminalValue = 0.5 * eN.' * P * eN;
            if obj.terminalViabilityEnabled
                terminalSlack = z(idx.terminalSlack);
                c = terminalValue - obj.aT - terminalSlack;
            else
                c = terminalValue - obj.aT;
            end
            if includeGradients
                gradc = zeros(nVar,1);
                gradc(idx.x{Nc+1}) = P*eN;
                if obj.terminalViabilityEnabled
                    gradc(idx.terminalSlack) = -1;
                end
            else
                gradc = [];
                if obj.acceptedReplayReuseActive
                    valueReplayCacheDecision = z;
                    valueReplayCacheInequality = c;
                    valueReplayCacheEquality = ceq;
                    valueReplayCacheEndpoints = currentReplayEndpoints;
                end
            end
        end

        function [endpoints,cacheInfo] = getAcceptedReplay(decision)
            cacheHit = obj.acceptedReplayReuseActive && ...
                isequal(decision,valueReplayCacheDecision) && ...
                ~isempty(valueReplayCacheEndpoints);
            if cacheHit
                endpoints = valueReplayCacheEndpoints;
            else
                endpoints = zeros(0,0);
            end
            cacheInfo = struct('enabled',obj.acceptedReplayReuseActive, ...
                'cacheHit',cacheHit, ...
                'valueReplayCacheHits',valueReplayCacheHits, ...
                'valueReplayCacheMisses',valueReplayCacheMisses);
        end

        % --------------------- BOUNDS --------------------------------
        xLrep = repmat(obj.xL,Nc+1,1);
        xUrep = repmat(obj.xU,Nc+1,1);

        uLrep = repmat(obj.uL,Nc,1);
        uUrep = repmat(obj.uU,Nc,1);

        nlp.cost = @costFun;
        nlp.nonl = @nonlFun;
        nlp.getAcceptedReplay = @getAcceptedReplay;
        if obj.preparedHorizonDataReuseActive
            nlp.valueHorizonData = valueHorizonData;
        end
        if obj.terminalViabilityEnabled
            nlp.lb = [xLrep; uLrep; 0];
            nlp.ub = [xUrep; uUrep; inf];
        else
            nlp.lb = [xLrep; uLrep];
            nlp.ub = [xUrep; uUrep];
        end
    end

    %------------------------------------------------------------------
    function nVar = decisionVariableCount(obj)
        nVar = (obj.Nc+1)*obj.nx + obj.Nc*obj.nu + ...
            double(obj.terminalViabilityEnabled);
    end

    %------------------------------------------------------------------
    function residual = surfaceSubspaceResidual(obj,control)
    %SURFACESUBSPACERESIDUAL Deflection/rate difference for two surfaces.
        assert(obj.n_surf==2 && size(control,1)==obj.nu, ...
            'nMPC:SymmetricSurfaceDimension', ...
            'The formal Case-A subspace requires two surface pairs.');
        residual = [control(1,:)-control(2,:); ...
            control(3,:)-control(4,:)];
    end

    %------------------------------------------------------------------
    function diagnostics = surfaceSubspaceDiagnostics(obj,control)
    %SURFACESUBSPACEDIAGNOSTICS Report without changing unrestricted use.
        diagnostics = struct('enabled', ...
            obj.symmetricSurfaceSubspaceEnabled, ...
            'caseId',obj.symmetricSurfaceSubspaceCaseId, ...
            'changeId',obj.symmetricSurfaceSubspaceChangeId, ...
            'maximumResidualInf',NaN);
        if obj.n_surf==2 && size(control,1)==obj.nu
            residual = obj.surfaceSubspaceResidual(control);
            diagnostics.maximumResidualInf = norm(residual(:),inf);
        end
    end

    %------------------------------------------------------------------
    function tipGradient = buildWingtipOutputGradient(obj,beam,base)
    %BUILDWINGTIPOUTPUTGRADIENT Central-difference physical-output map.
    % The map is fixed at the current trim solely for the audit objective.
        x0 = obj.xTrim;
        h = 1e-6*max(1,abs(x0));
        tipGradient = zeros(1,obj.nx);
        for column = 1:obj.nx
            xp = x0; xm = x0;
            xp(column) = xp(column) + h(column);
            xm(column) = xm(column) - h(column);
            tipGradient(column) = (obj.wingtipHeave(xp,beam,base) - ...
                obj.wingtipHeave(xm,beam,base))/(2*h(column));
        end
        assert(all(isfinite(tipGradient)), 'nMPC:WingtipOutputAuditMap', ...
            'The trim-linear wingtip-heave output gradient is non-finite.');
    end

    %------------------------------------------------------------------
    function tip = wingtipHeave(obj,state,beam,base)
    %WINGTIPHEAVE Recover physical A-frame vertical deflection at the tip.
        q2 = state(obj.model.idx.q2);
        qxi = state(obj.model.idx.qxi);
        q0 = -beam.Omega \ q2;
        x0 = beam.phi0*q0;
        xiFull = base.phi_xi_modes*qxi;
        nNode = double(beam.fem.num_node)-1;
        nDof = double(beam.nFlex)/nNode;
        selectTip = double((nNode/2+1)*(nDof+1)) + (1:3);
        transform = sparse(7*nNode,6*nNode);
        quaternion = reshape(xiFull,4,[]).';
        for node = 1:nNode
            i6 = (node-1)*6+(1:6);
            i7 = (node-1)*7+(1:7);
            transform(i7,i6) = blkdiag( ...
                AeroFlex.core.T_phi_quat(quaternion(node,:)), ...
                AeroFlex.core.halfTangentialOperator(quaternion(node,:)));
        end
        x0A = transform*x0;
        tip = x0A(selectTip(3));
        assert(isscalar(tip) && isfinite(tip), 'nMPC:WingtipOutputAuditMap', ...
            'The recovered wingtip heave is non-finite.');
    end

    %------------------------------------------------------------------
    function terminalInfo = evaluateTerminal(obj,z)
        idx = obj.buildIndexMaps();
        terminalError = z(idx.x{obj.Nc+1}) - obj.xTrim;
        terminalInfo.value = 0.5*terminalError.'*obj.Pc*terminalError - obj.aT;
        terminalInfo.mode = "hard";
        terminalInfo.slack = 0;
        if obj.terminalViabilityEnabled
            terminalInfo.mode = "viability_slack";
            terminalInfo.slack = z(idx.terminalSlack);
        end
        terminalInfo.satisfied = terminalInfo.value <= 0;
    end

    %------------------------------------------------------------------
    function [uL,uU,urL,urU] = buildControlBounds(obj,cfg)
    % Builds full nu x 1 lower/upper bounds for the control decision.
    %
    % Supported cases:
    %
    % 1) cfg.uL and cfg.uU are already nu x 1.
    %
    % 2) var_per == 2 and nu == 2*n_surf:
    %      U = [delta_1; ...; delta_nsurf; rate_1; ...; rate_nsurf]
    %
    %    Then scalar/vector cfg.uL/cfg.uU are used for deflection bounds,
    %    and cfg.urateL/cfg.urateU are used for rate-channel bounds.
    %
    % If your ordering is instead [delta_1; rate_1; delta_2; rate_2],
    % change this helper accordingly.
    %------------------------------------------------------------------
        nSurf = cfg.ctrl.n_surf;
        varPer = cfg.ctrl.var_per;

        uLraw = cfg.uL;
        uUraw = cfg.uU;

        hasRateBounds = isfield(cfg,'urateL') && isfield(cfg,'urateU');

        if numel(uLraw)==obj.nu && numel(uUraw)==obj.nu
            uL = uLraw(:);
            uU = uUraw(:);

            if hasRateBounds
                urL = cfg.urateL(:);
                urU = cfg.urateU(:);
            else
                urL = [];
                urU = [];
            end
        elseif hasRateBounds && varPer==2 && obj.nu==2*nSurf
            defL  = obj.expandToLength(cfg.uL,nSurf,'cfg.uL');
            defU  = obj.expandToLength(cfg.uU,nSurf,'cfg.uU');
            rateL = obj.expandToLength(cfg.urateL,nSurf,'cfg.urateL');
            rateU = obj.expandToLength(cfg.urateU,nSurf,'cfg.urateU');

            uL = [defL; rateL];
            uU = [defU; rateU];

            urL = rateL;
            urU = rateU;
        else
            uL = obj.expandToLength(uLraw,obj.nu,'cfg.uL');
            uU = obj.expandToLength(uUraw,obj.nu,'cfg.uU');

            if hasRateBounds
                urL = cfg.urateL(:);
                urU = cfg.urateU(:);
            else
                urL = [];
                urU = [];
            end
        end

        % Configuration limits are physical total surface limits. Convert
        % only the deflection channels to the controller's increment frame.
        uL(1:nSurf) = uL(1:nSurf) - obj.uModelTrim(1:nSurf);
        uU(1:nSurf) = uU(1:nSurf) - obj.uModelTrim(1:nSurf);
        assert(all(isfinite(uL) | isinf(uL)) && ...
            all(isfinite(uU) | isinf(uU)) && all(uL <= uU), ...
            'nMPC:ControlBounds','Control bounds are invalid after trim shift.');
    end

    %------------------------------------------------------------------
    function uModelTrim = buildModelTrim(obj,trim)
    % Resolve the total model trim input in radians and radians per second.
        nSurf = obj.n_surf;
        uModelTrim = zeros(obj.nu,1);

        if isfield(trim,'u_ctrl') && ~isempty(trim.u_ctrl)
            uModelTrim = obj.expandToLength( ...
                trim.u_ctrl,obj.nu,'trim.u_ctrl');
        elseif isfield(trim,'deltaWing') && ~isempty(trim.deltaWing)
            delta = obj.expandToLength( ...
                trim.deltaWing,nSurf,'trim.deltaWing');
            uModelTrim(1:nSurf) = delta;
        elseif isfield(trim,'deltaWingDeg') && ~isempty(trim.deltaWingDeg)
            delta = obj.expandToLength( ...
                deg2rad(trim.deltaWingDeg),nSurf,'trim.deltaWingDeg');
            uModelTrim(1:nSurf) = delta;
        elseif isfield(trim,'deltaDeg') && ~isempty(trim.deltaDeg)
            deltaRaw = trim.deltaDeg(:);
            if numel(deltaRaw) == obj.nu
                deltaRaw = deltaRaw(1:nSurf);
            end
            delta = obj.expandToLength( ...
                deg2rad(deltaRaw),nSurf,'trim.deltaDeg');
            uModelTrim(1:nSurf) = delta;
        else
            error('nMPC:TrimControl', ...
                'A total wing trim control is required.');
        end

        assert(all(isfinite(uModelTrim)), 'nMPC:TrimControl', ...
            'The total model trim control must be finite.');
        assert(norm(uModelTrim(nSurf+1:end),inf) <= 100*eps, ...
            'nMPC:TrimRate', ...
            'The locked trim must have zero surface-rate channels.');
    end

    %------------------------------------------------------------------
    function v = expandToLength(~,v,n,name)
    % Converts scalar/vector input to n x 1 column with validation.
    %------------------------------------------------------------------
        v = v(:);

        if isempty(v)
            error('nMPC:Config','%s is empty.',name);
        end

        if isscalar(v)
            v = repmat(v,n,1);
            return
        end

        if numel(v) ~= n
            error('nMPC:Dimension', ...
                '%s must be scalar or length %d. Got length %d.', ...
                name,n,numel(v));
        end
    end

    %------------------------------------------------------------------
    function digest = localFileHash(~,path)
        engine = javaMethod( ...
            'getInstance','java.security.MessageDigest','SHA-256');
        fileId = fopen(path,'rb');
        assert(fileId>=0,'nMPC:ReciprocalProviderHashOpen', ...
            'Cannot open %s.',path);
        cleanup = onCleanup(@() fclose(fileId));
        while ~feof(fileId)
            bytes = fread(fileId,1024*1024,'*uint8');
            if isempty(bytes), break, end
            engine.update(typecast(bytes(:),'int8'));
        end
        digest = lower(string(reshape(dec2hex( ...
            typecast(engine.digest(),'uint8'),2).',1,[])));
        clear cleanup
    end

    %------------------------------------------------------------------
    function debugPlots(obj,t_k,info)
    % Shows how the optimized control horizon evolves over time.
    %------------------------------------------------------------------
        uH = info.uHorizon;   % expected size: nu x Nc

        if isempty(obj.dbg) || ~isfield(obj.dbg,'t')
            obj.dbg.t    = [];
            obj.dbg.U    = [];
            obj.dbg.cont = [];
        end

        obj.dbg.t(end+1)      = t_k;
        obj.dbg.U(end+1,:)    = uH(:).';
        obj.dbg.cont(end+1)   = info.continuity;

        figure(5001); clf

        for ii = 1:obj.nu
            subplot(obj.nu,1,ii); hold on

            cols = ii:obj.nu:size(obj.dbg.U,2);

            for s = 1:obj.Nc
                stairs(obj.dbg.t,obj.dbg.U(:,cols(s)),'LineWidth',1.0);
            end

            hold off
            grid on
            ylabel(sprintf('u_%d',ii));

            if ii==1
                title('NMPC optimized control horizons');
            end
        end

        xlabel('time [s]');

        figure(5004); clf
        plot(obj.dbg.t,obj.dbg.cont,'LineWidth',1.0);
        grid on
        xlabel('time [s]');
        ylabel('||ceq||_2');
        title('NMPC multiple-shooting feasibility');
    end
        function localCheckNMPCEqualityGradient(obj,nlp,z,lb,ub)
        %LOCALCHECKNMPCEQUALITYGRADIENT Split equality-gradient check into
        % initial-condition rows and dynamic-continuity rows.
        
            nx = obj.nx;
            Nc = obj.Nc;
        
            z  = z(:);
            lb = lb(:);
            ub = ub(:);
            n  = numel(z);
        
            z = min(max(z,lb),ub);
        
            [~,ceq0,~,Gceq] = nlp.nonl(z);
        
            if size(Gceq,1) ~= n
                Gceq = Gceq.';
            end
        
            rng(11);
            d = randn(n,1);
            d = d / max(norm(d),eps);
        
            h = 1e-6;
        
            active = abs(d) > 0;
            hBound = inf;
        
            idxU = active & isfinite(ub);
            if any(idxU)
                hBound = min(hBound,min((ub(idxU)-z(idxU))./abs(d(idxU))));
            end
        
            idxL = active & isfinite(lb);
            if any(idxL)
                hBound = min(hBound,min((z(idxL)-lb(idxL))./abs(d(idxL))));
            end
        
            if isfinite(hBound)
                h = min(h,0.49*hBound);
            end
        
            if h <= eps
                fprintf('[nMPC grad check] Could not find bound-safe FD step.\n');
                return
            end
        
            [~,ceqp] = nlp.nonl(z + h*d);
            [~,ceqm] = nlp.nonl(z - h*d);
        
            ceqFD = (ceqp(:)-ceqm(:))/(2*h);
            ceqAN = Gceq.'*d;
        
            % Expected layout:
            %   ceq = [X0 - xhat; defects]
            if numel(ceq0) >= nx*(Nc+1)
                rowsInit = 1:nx;
                rowsDyn  = nx+1:numel(ceq0);
            else
                rowsInit = [];
                rowsDyn  = 1:numel(ceq0);
            end
        
            denAll = max([1,norm(ceqFD,inf),norm(ceqAN,inf)]);
            errAll = norm(ceqFD-ceqAN,inf)/denAll;
        
            if ~isempty(rowsInit)
                denInit = max([1,norm(ceqFD(rowsInit),inf),norm(ceqAN(rowsInit),inf)]);
                errInit = norm(ceqFD(rowsInit)-ceqAN(rowsInit),inf)/denInit;
            else
                errInit = nan;
            end
        
            denDyn = max([1,norm(ceqFD(rowsDyn),inf),norm(ceqAN(rowsDyn),inf)]);
            errDyn = norm(ceqFD(rowsDyn)-ceqAN(rowsDyn),inf)/denDyn;
        
            fprintf('\n');
            fprintf('======================================================================\n');
            fprintf(' NMPC EQUALITY-GRADIENT SPLIT CHECK\n');
            fprintf('======================================================================\n');
            fprintf('  h                         : %.3e\n', h);
            fprintf('  total equality rel error  : %.3e\n', errAll);
            fprintf('  initial-state rel error   : %.3e\n', errInit);
            fprintf('  dynamic-defect rel error  : %.3e\n', errDyn);
            fprintf('  ||ceqFD||inf              : %.3e\n', norm(ceqFD,inf));
            fprintf('  ||ceqAN||inf              : %.3e\n', norm(ceqAN,inf));
            fprintf('  ||ceqFD-ceqAN||inf        : %.3e\n', norm(ceqFD-ceqAN,inf));
            fprintf('======================================================================\n\n');
        end
        function localCheckNMPCEqualityGradientBlocks(obj,nlp,z,lb,ub)
            %LOCALCHECKNMPCEQUALITYGRADIENTBLOCKS Check dynamic equality Jacobian by
            % perturbing only X variables and only U variables separately.
            
                idx = obj.buildIndexMaps();
            
                nx = obj.nx;
                nu = obj.nu;
                Nc = obj.Nc;
            
                z  = z(:);
                lb = lb(:);
                ub = ub(:);
                n  = numel(z);
            
                z = min(max(z,lb),ub);
            
                [~,ceq0,~,Gceq] = nlp.nonl(z);
            
                if size(Gceq,1) ~= n
                    Gceq = Gceq.';
                end
            
                rowsDyn = nx+1:numel(ceq0);  % assumes ceq = [X0-xhat; dynamics]
            
                h = 1e-6;
            
                % --------------------------------------------------------------
                % Build state-only direction
                % --------------------------------------------------------------
                dX = zeros(n,1);
            
                rng(101);
                for k = 1:Nc+1
                    dX(idx.x{k}) = randn(nx,1);
                end
            
                dX = dX / max(norm(dX),eps);
            
                % --------------------------------------------------------------
                % Build control-only direction
                % --------------------------------------------------------------
                dU = zeros(n,1);
            
                rng(202);
                for k = 1:Nc
                    dU(idx.u{k}) = randn(nu,1);
                end
            
                dU = dU / max(norm(dU),eps);
            
                % --------------------------------------------------------------
                % Evaluate errors
                % --------------------------------------------------------------
                errX = obj.localDirectionalEqualityError(nlp,z,dX,h,lb,ub,Gceq,rowsDyn);
                errU = obj.localDirectionalEqualityError(nlp,z,dU,h,lb,ub,Gceq,rowsDyn);
            
                fprintf('\n');
                fprintf('======================================================================\n');
                fprintf(' NMPC DYNAMIC-EQUALITY BLOCK GRADIENT CHECK\n');
                fprintf('======================================================================\n');
                fprintf('  dynamic defect rel error, X-only direction : %.3e\n', errX);
                fprintf('  dynamic defect rel error, U-only direction : %.3e\n', errU);
                fprintf('======================================================================\n\n');
        end
        function err = localDirectionalEqualityError(obj,nlp,z,d,h,lb,ub,Gceq,rowsDyn)
            %LOCALDIRECTIONALEQUALITYERROR Directional FD check for selected equality rows.
            %
            % This is used by localCheckNMPCEqualityGradientBlocks to isolate whether
            % the dynamic-defect Jacobian error comes from X sensitivities or U
            % sensitivities.
                       
                z  = z(:);
                d  = d(:);
                lb = lb(:);
                ub = ub(:);
            
                active = abs(d) > 0;
                hBound = inf;
            
                idxU = active & isfinite(ub);
                if any(idxU)
                    hBound = min(hBound,min((ub(idxU)-z(idxU))./abs(d(idxU))));
                end
            
                idxL = active & isfinite(lb);
                if any(idxL)
                    hBound = min(hBound,min((z(idxL)-lb(idxL))./abs(d(idxL))));
                end
            
                if isfinite(hBound)
                    h = min(h,0.49*hBound);
                end
            
                if h <= eps
                    err = nan;
                    return
                end
            
                [~,ceqp] = nlp.nonl(z + h*d);
                [~,ceqm] = nlp.nonl(z - h*d);
            
                ceqFD = (ceqp(:)-ceqm(:))/(2*h);
                ceqAN = Gceq.'*d;
            
                rowsDyn = rowsDyn(:);
            
                den = max([1,norm(ceqFD(rowsDyn),inf),norm(ceqAN(rowsDyn),inf)]);
            
                err = norm(ceqFD(rowsDyn)-ceqAN(rowsDyn),inf)/den;
            end
    %------------------------------------------------------------------
    end % private methods
end

function enabled = localScheduledAggregateRequested(request)
enabled = false;
if ~isfield(request,'runtimeAcceleration') || ...
        ~isstruct(request.runtimeAcceleration) || ...
        ~isscalar(request.runtimeAcceleration) || ...
        ~isfield(request.runtimeAcceleration,'enabled') || ...
        ~logical(request.runtimeAcceleration.enabled) || ...
        ~isfield(request.runtimeAcceleration,'scheduledAggregate')
    return
end
selector = request.runtimeAcceleration.scheduledAggregate;
assert(isscalar(selector) && (islogical(selector) || ...
    (isnumeric(selector) && isfinite(selector) && ismember(selector,[0 1]))), ...
    'nMPC:ScheduledAggregateSelector', ...
    'The scheduled aggregate selector must be logical.');
enabled = logical(selector);
end
