classdef nMHEv2 < AeroFlex.ctrl.EstimatorBase
%========================================================
%  NON-LINEAR  MOVING–HORIZON ESTIMATOR  (multiple-shoot)
%========================================================
%  * buffer-free – stores its own X,U,Y,W histories
%  * Diehl-style warm start (shift & predict)
%  * Exact Jacobian of continuity constraints
%========================================================

    properties
        % ---- scalar config
        dt  double
        Ts  double
        Ne  double
        method (1,1) string {mustBeMember(method,["multiple","single"])} = "multiple"
        % ---- dimensions
        nx  double
        nu  double
        ny  double
        nw  double
        % ---- weights / bounds
        Qe  double;  Pe double;  Re double
        xL  double;  xU double
        wL  double;  wU double
        uModelTrim double
        stateModelTrim double
        % ---- model & buffers
        model    % ROMIntegrator
        sensor
        measurementContract struct
        Sprev
        Xhist
        Uhist
        Yhist
        measurementContextHistory double
        knownChiContextEnabled logical = false
        knownChiHistory double
        knownChiIndex double
        Whist
        k   double                         % sample counter
        % ---- optimiser state
        z0  double
        H   double
        solverOpts
        firstFullSolve
        % ---- outputs
        % xhat double
        % what double
        % ---- debug
        debug logical = false
        dbg   struct
        xlast

        solverName string = "fmincon"
        coldStartStrategy (1,1) string = "trim"
        sqpSolver
        sqpCheckDone logical = false
        fminconIterationDiagnostics logical = false
        fminconIterationTrace struct = struct()
        fminconIterationPlot struct = struct()
        candidateQualityDiagnostics logical = false
        preparationFeedbackAudit logical = false
        preparedFeedback struct
        scheduledPackageHistoryEnabled logical = false
        scheduledPackageHistory cell
        scheduledPackageHistoryValid logical = true
        scheduledPackageHistoryReason string = "disabled"
        incomingPackageStateToActive double
        nativeStateCount (1,1) double
        reciprocalProviderEnabled logical = false
        reciprocalProvider
        reciprocalProviderFormulation (1,1) string = "disabled"
        reciprocalCondensedEnabled logical = false
        reciprocalDisturbanceWarmStart (1,1) string = "hold_last"
        reciprocalUnboundedAerodynamicLagStates logical = false
        reciprocalLatentIndex double
        reciprocalLatentTrim double
        reciprocalLatentHistory double
        reciprocalContextHistory cell
        reciprocalPendingContext struct = struct()
        reciprocalPendingContextValid logical = false
        nativeMeasurementContract struct
        % Default-inactive exact full-state-condensing RTI seed. The full
        % multiple-shooting NLP remains the correction and fallback owner.
        realtimeRtiEnabled logical = false
        realtimeRtiChangeId (1,1) string = ""
        realtimeRtiIterationCount (1,1) double = 3
        realtimeRtiSolver
        realtimeRtiSolutionOwnerEnabled logical = false
        realtimeRtiNonlinearBacktracking struct = struct('enabled',false)
        realtimeRtiFullCorrectionFallback struct = struct('enabled',false)
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
        terminalStmCondensationRequested logical = false
        terminalStmCondensationActive logical = false
        latentHistoryCondensationRequested logical = false
        latentHistoryCondensationActive logical = false
        preparedHorizonDataReuseRequested logical = false
        preparedHorizonDataReuseActive logical = false
        acceptedReplayReuseRequested logical = false
        acceptedReplayReuseActive logical = false
    end
%======================================================================
methods
%----------------------------------------------------------------------
function obj = nMHEv2(cfg,beam,aero,base, trim)
% Constructor
%----------------------------------------------------------------------
    obj@AeroFlex.ctrl.EstimatorBase(cfg, trim);

    % -------- 1. basic numbers ---------------------------------------
    obj.dt   = cfg.sim.dt;
    obj.Ts   = cfg.ctrl.Ts;
    obj.Ne   = cfg.ctrl.Ne;
    obj.method = lower(string(cfg.ctrl.mhe_method));
    if isfield(cfg.ctrl,'mheColdStartStrategy') && ...
            ~isempty(cfg.ctrl.mheColdStartStrategy)
        coldStartStrategy = lower(string(cfg.ctrl.mheColdStartStrategy));
        assert(isscalar(coldStartStrategy) && ...
            ismember(coldStartStrategy,["trim","linear_batch"]), ...
            'nMHEv2:ColdStartStrategy', ...
            ['cfg.ctrl.mheColdStartStrategy must be "trim" or ', ...
             '"linear_batch".']);
        obj.coldStartStrategy = coldStartStrategy;
    end

    obj.model  = cfg.modelHandle(cfg,beam,aero,base);
    obj.sensor = cfg.sensorHandle(beam,cfg);

    if isfield(cfg,'ctrl') && isfield(cfg.ctrl, ...
            'mheScheduledPackageHistoryAudit') && ...
            ~isempty(cfg.ctrl.mheScheduledPackageHistoryAudit)
        request = cfg.ctrl.mheScheduledPackageHistoryAudit;
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled') && isfield(request,'auditOnly') && ...
            isscalar(request.enabled) && isscalar(request.auditOnly) && ...
            logical(request.auditOnly), ...
            'nMHEv2:ScheduledPackageHistoryRequest', ...
            ['Scheduled package history is available only through an ', ...
             'explicit audit-only request.']);
        obj.scheduledPackageHistoryEnabled = logical(request.enabled);
    end

    obj.nativeStateCount = size(obj.model.L,1);
    obj.nx = obj.nativeStateCount;
    obj.nu = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    obj.nw = cfg.nw;
    if isfield(cfg.ctrl,'mheKnownRigidAttitudeChi') && ...
            ~isempty(cfg.ctrl.mheKnownRigidAttitudeChi)
        request = cfg.ctrl.mheKnownRigidAttitudeChi;
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enable') && isfield(request,'auditOnly') && ...
            isscalar(request.enable) && isscalar(request.auditOnly) && ...
            logical(request.auditOnly), ...
            'nMHEv2:KnownChiRequest', ...
            ['Known rigid-attitude chi context is available only through ', ...
             'an explicit audit-only request.']);
        obj.knownChiContextEnabled = logical(request.enable);
    end
    if obj.knownChiContextEnabled
        assert(isfield(obj.model.idx,'chi') && ...
            numel(obj.model.idx.chi) == 3, ...
            'nMHEv2:KnownChiIndex', ...
            'The known rigid-attitude context requires exactly three chi states.');
        obj.knownChiIndex = obj.model.idx.chi(:);
    else
        obj.knownChiIndex = zeros(0,1);
    end
    nativeTrimState = trim.states(:);
    assert(numel(nativeTrimState)==obj.nativeStateCount && ...
        all(isfinite(nativeTrimState)), ...
        'nMHEv2:NativeTrimState', ...
        'The native trim state is incompatible with the prediction model.');
    obj.measurementContract = ...
        obj.buildMeasurementContract(cfg,nativeTrimState);
    obj.nativeMeasurementContract = obj.measurementContract;
    obj.uModelTrim = obj.buildModelTrim(trim);
    xTrim = obj.configureReciprocalProvider(cfg,nativeTrimState);
    obj.ny = obj.measurementContract.ny;

    assert(obj.nx > 0 && obj.nu > 0 && obj.ny > 0 && obj.nw > 0, ...
        'nMHEv2:Dimensions', ...
        'Estimator state, input, measurement, and disturbance dimensions must be positive.');
    nSubsteps = round(obj.Ts/obj.dt);
    sampleTolerance = 100*eps(max([1,obj.Ts,obj.dt]));
    assert(obj.Ts > 0 && obj.dt > 0 && nSubsteps >= 1 && ...
        abs(nSubsteps*obj.dt-obj.Ts) <= sampleTolerance, ...
        'nMHEv2:SampleAlignment', ...
        'Estimator Ts must be an integer multiple of prediction-model dt.');
    assert(cfg.ny == obj.ny, 'nMHEv2:MeasurementMap', ...
        'cfg.ny must match the selected measurement contract.');
    % -------- 2. weights / bounds ------------------------------------
    obj.Qe = cfg.Qe;     obj.Pe = cfg.Pe;     obj.Re = cfg.Re;
    if obj.reciprocalProviderEnabled && ~obj.reciprocalCondensedEnabled
        assert(isequal(size(obj.Pe), ...
            [obj.nativeStateCount,obj.nativeStateCount]), ...
            'nMHEv2:ReciprocalNativeArrivalWeight', ...
            ['The reciprocal audit must retain the unchanged native ', ...
             'arrival-weight dimensions.']);
        latentCount = numel(obj.reciprocalLatentIndex);
        obj.Pe = blkdiag(obj.Pe,zeros(latentCount,latentCount));
    end
    assert(numel(xTrim) == obj.nx && all(isfinite(xTrim)), ...
        'nMHEv2:TrimState', ...
        'trim.states must contain %d finite prediction states.',obj.nx);
    obj.stateModelTrim = xTrim;
    assert(isequal(size(obj.Qe),[obj.ny,obj.ny]) && ...
        all(isfinite(obj.Qe),'all'), ...
        'nMHEv2:MeasurementWeight', ...
        'cfg.Qe must be a finite ny-by-ny matrix.');
    assert(isequal(size(obj.Pe),[obj.nx,obj.nx]) && ...
        all(isfinite(obj.Pe),'all'), ...
        'nMHEv2:ArrivalWeight', ...
        'cfg.Pe must be a finite nx-by-nx matrix.');
    assert(isequal(size(obj.Re),[obj.nw,obj.nw]) && ...
        all(isfinite(obj.Re),'all'), ...
        'nMHEv2:DisturbanceWeight', ...
        'cfg.Re must be a finite nw-by-nw matrix.');
    if isfield(cfg,'RdWe') && ~isempty(cfg.RdWe)
        assert(isequal(size(cfg.RdWe),[obj.nw,obj.nw]) && ...
            all(isfinite(cfg.RdWe),'all'), ...
            'nMHEv2:DisturbanceRateWeight', ...
            'cfg.RdWe must be a finite nw-by-nw matrix.');
    end
    if isfield(cfg,'mhe') && isfield(cfg.mhe,'RwTerminal') && ...
            ~isempty(cfg.mhe.RwTerminal)
        assert(isequal(size(cfg.mhe.RwTerminal),[obj.nw,obj.nw]) && ...
            all(isfinite(cfg.mhe.RwTerminal),'all'), ...
            'nMHEv2:TerminalDisturbanceWeight', ...
            'cfg.mhe.RwTerminal must be a finite nw-by-nw matrix.');
    end
    % margin = 1.60;                            % 30 % either side
    margin = 10;                            % 30 % either side
    % cfg.xL = (1 - margin) .* xTrim;
    % cfg.xU = (1 + margin) .* xTrim;

    % absFloor = 1e-1;  
    absFloor = 1;  

    delta   = max(abs(margin .* xTrim) , absFloor);   % n×1 non‑negative
    cfg.xL      = xTrim - delta;
    cfg.xU      = xTrim + delta;
    obj.xL = cfg.xL;     obj.xU = cfg.xU;
    if obj.reciprocalUnboundedAerodynamicLagStates
        aerodynamicLagIndex = obj.model.idx.qGam(:);
        assert(all(aerodynamicLagIndex>=1 & ...
            aerodynamicLagIndex<=obj.nativeStateCount), ...
            'nMHEv2:ReciprocalAerodynamicLagIndex', ...
            'The reciprocal aerodynamic-lag state index is invalid.');
        obj.xL(aerodynamicLagIndex) = -inf;
        obj.xU(aerodynamicLagIndex) = inf;
    end
    obj.wL = obj.expandToLength(cfg.wL,obj.nw,'cfg.wL');
    obj.wU = obj.expandToLength(cfg.wU,obj.nw,'cfg.wU');
    finiteStateBound = isfinite(obj.xL) & isfinite(obj.xU);
    if obj.reciprocalUnboundedAerodynamicLagStates
        requiredFinite = true(obj.nx,1);
        requiredFinite(obj.model.idx.qGam(:)) = false;
        validUnbounded = all(isinf(obj.xL(~requiredFinite)) & ...
            obj.xL(~requiredFinite)<0 & isinf(obj.xU(~requiredFinite)) & ...
            obj.xU(~requiredFinite)>0);
    else
        requiredFinite = true(obj.nx,1);
        validUnbounded = true;
    end
    assert(all(finiteStateBound(requiredFinite)) && validUnbounded && ...
        all(~isnan(obj.xL)) && all(~isnan(obj.xU)) && ...
        all(obj.xL <= obj.xU), ...
        'nMHEv2:StateBounds','The estimator state bounds are invalid.');
    assert(all(isfinite(obj.wL) | isinf(obj.wL)) && ...
        all(isfinite(obj.wU) | isinf(obj.wU)) && all(obj.wL <= obj.wU), ...
        'nMHEv2:DisturbanceBounds', ...
        'The estimator disturbance bounds are invalid.');

    if isfield(cfg,'debug') && isstruct(cfg.debug) && ...
            isfield(cfg.debug,'level') && ~isempty(cfg.debug.level)
        obj.debug = cfg.debug.level >= 3;
    end
    if isfield(cfg,'ctrl') && isfield(cfg.ctrl, ...
            'mheFminconIterationDiagnostics') && ...
            ~isempty(cfg.ctrl.mheFminconIterationDiagnostics)
        request = cfg.ctrl.mheFminconIterationDiagnostics;
        assert(isscalar(request) && (islogical(request) || ...
            (isnumeric(request) && isfinite(request) && ...
            ismember(request,[0 1]))), ...
            'nMHEv2:FminconIterationDiagnostics', ...
            'cfg.ctrl.mheFminconIterationDiagnostics must be a logical scalar.');
        obj.fminconIterationDiagnostics = logical(request);
    end
    if isfield(cfg,'ctrl') && isfield(cfg.ctrl, ...
            'mheCandidateQualityDiagnostics') && ...
            ~isempty(cfg.ctrl.mheCandidateQualityDiagnostics)
        request = cfg.ctrl.mheCandidateQualityDiagnostics;
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled') && isfield(request,'auditOnly') && ...
            isscalar(request.enabled) && isscalar(request.auditOnly) && ...
            logical(request.auditOnly), ...
            'nMHEv2:CandidateQualityDiagnostics', ...
            ['Candidate-quality diagnostics are available only through ', ...
             'an explicit audit-only request.']);
        obj.candidateQualityDiagnostics = logical(request.enabled);
    end

    % -------- 3. optimiser options -----------------------------------
    obj.solverOpts = optimoptions('fmincon',...
        'Algorithm','sqp',...     % interior-point | sqp
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true,...
        'HessianApproximation','lbfgs',...
        'Display','none',...
        'OptimalityTolerance',1e-4, ...
        'StepTolerance',1e-7, ...
        'ConstraintTolerance',1e-6);

    % obj.solverOpts.FiniteDifferenceType ='central';
    % obj.solverOpts.CheckGradients = true;
            % 'FiniteDifferenceStepSize',1e-8,...

    %======================================================================
    % nMHE solver selection
    %======================================================================
    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'mheSolver')
        obj.solverName = lower(string(cfg.ctrl.mheSolver));
    else
        obj.solverName = "fmincon";
    end
    
    switch obj.solverName
    
        case "fmincon"
            obj.sqpSolver = [];
    
            % Make sure gradient checking is not accidentally active online.
            if isprop(obj.solverOpts,'CheckGradients')
                obj.solverOpts.CheckGradients = false;
            end
            if obj.fminconIterationDiagnostics
                obj.solverOpts.OutputFcn = @(z,values,state) ...
                    obj.fminconOutputFunction(z,values,state);
            end
    
        case "custom_sqp"
            sqpOpts = AeroFlex.optim.SQPSolver.defaultOptions();
    
            if isfield(cfg.ctrl,'mheSqp') && ~isempty(cfg.ctrl.mheSqp)
                f = fieldnames(cfg.ctrl.mheSqp);
                for ii = 1:numel(f)
                    sqpOpts.(f{ii}) = cfg.ctrl.mheSqp.(f{ii});
                end
            end
    
            obj.sqpSolver = AeroFlex.optim.SQPSolver(sqpOpts);
    
        otherwise
            error('nMHEv2:Solver', ...
                  'Unknown cfg.ctrl.mheSolver = "%s". Use "fmincon" or "custom_sqp".', ...
                  obj.solverName);
    end

    obj.configureRealtimeRtiAudit(cfg);
    obj.configureNativeReducedHorizonRtiAudit(cfg);
    obj.configureHighLeverageRuntimeAudit(cfg);

    obj.sqpCheckDone = ~obj.debug;
    obj.preparedFeedback = struct('valid',false);
    if isfield(cfg.ctrl,'mhePreparationFeedbackAudit') && ...
            ~isempty(cfg.ctrl.mhePreparationFeedbackAudit)
        request = cfg.ctrl.mhePreparationFeedbackAudit;
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled') && isfield(request,'auditOnly') && ...
            isscalar(request.enabled) && isscalar(request.auditOnly) && ...
            logical(request.auditOnly), ...
            'nMHEv2:PreparationFeedbackAudit', ...
            'Preparation/feedback operation is available only by audit-only request.');
        obj.preparationFeedbackAudit = logical(request.enabled);
        assert(~obj.preparationFeedbackAudit || obj.solverName == "custom_sqp", ...
            'nMHEv2:PreparationFeedbackSolver', ...
            'Preparation/feedback audit requires the custom SQP estimator.');
    end

    % -------- 4. rolling buffers  ------------------------------------
    obj.Xhist = zeros(obj.nx,obj.Ne+1);
    obj.Uhist = zeros(obj.nu,obj.Ne);
    obj.Yhist = zeros(obj.ny,obj.Ne+1);
    obj.measurementContextHistory = repmat( ...
        obj.measurementContract.linearizationAtTrim.context(:),1,obj.Ne+1);
    if obj.knownChiContextEnabled
        obj.knownChiHistory = repmat( ...
            xTrim(obj.knownChiIndex),1,obj.Ne+1);
    else
        obj.knownChiHistory = zeros(0,obj.Ne+1);
    end
    obj.Whist = zeros(obj.nw,obj.Ne);
    obj.scheduledPackageHistory = cell(1,obj.Ne);
    obj.incomingPackageStateToActive = eye(obj.nx);
    if obj.scheduledPackageHistoryEnabled
        obj.scheduledPackageHistoryValid = false;
        obj.scheduledPackageHistoryReason = "awaiting_completed_intervals";
    end
    obj.k     = 0;

    % -------- 5. initial guess ---------------------------------------
    % xTrim already contains the configured native or audit-augmented
    % equilibrium state. Do not replace it with the native trim vector here.
    obj.xhat        = xTrim;
    obj.Xhist(:,:)  = repmat(xTrim,1,obj.Ne+1);
    obj.Whist(:,:)  = 0;                      % zero-gust history
    obj.what = zeros(obj.Ne*obj.nw,1);

    idx   = obj.buildIndexMaps();
    obj.z0           = zeros((obj.Ne+1)*obj.nx + obj.Ne*obj.nw,1);
    for j = 1:obj.Ne+1
        obj.z0(idx.x{j}) = xTrim;
    end
    % obj.z0(idx.x{:}) = xTrim;                 % arrival state
    % obj.z0(idx.x{end}) = xTrim;                 % arrival state
    obj.H            = speye(numel(obj.z0));  % LBFGS initial Hessian

    % initial STM (∂x/∂x0  and  ∂x/∂w) – identity / zeros
    obj.Sprev = [eye(obj.nx) , zeros(obj.nx,obj.nw+obj.nu)];
    obj.dbg = struct('t', [], 'W', [], 'cont', 0);
    obj.firstFullSolve = true;
end
%----------------------------------------------------------------------
function synchronizeTimingFromModel(obj)
%SYNCHRONIZETIMINGFROMMODEL Adopt the active scheduled-model sample step.
    assert(isprop(obj.model,'dt') && isscalar(obj.model.dt) && ...
        isfinite(obj.model.dt) && obj.model.dt > 0, ...
        'nMHEv2:ModelStep', ...
        'The active estimator prediction model must define a positive scalar dt.');
    obj.dt = obj.model.dt;
    nSubsteps = round(obj.Ts/obj.dt);
    sampleTolerance = 100*eps(max([1,obj.Ts,obj.dt]));
    assert(nSubsteps >= 1 && abs(nSubsteps*obj.dt-obj.Ts) <= sampleTolerance, ...
        'nMHEv2:SampleAlignment', ...
        'Estimator Ts must be an integer multiple of the active prediction-model dt.');
end
%----------------------------------------------------------------------
function info = transportScheduledState(obj,transform,newReference,cfg)
%TRANSPORTSCHEDULEDSTATE Move persistent estimator state to a new chart.
    transform = full(transform);
    newReference = newReference(:);
    assert(isequal(size(transform),[obj.nx,obj.nx]) && ...
        all(isfinite(transform),'all') && numel(newReference) == obj.nx && ...
        all(isfinite(newReference)), ...
        'nMHEv2:ScheduledStateTransport', ...
        'The scheduled estimator state map is invalid.');
    inverse = transform\eye(obj.nx);
    newContract = obj.buildMeasurementContract(cfg,newReference);
    assert(newContract.ny == obj.ny && isequal(size(obj.Qe),[obj.ny,obj.ny]), ...
        'nMHEv2:ScheduledMeasurementDimension', ...
        'A scheduled measurement-map rebuild changed the measurement dimension.');
    index = obj.buildIndexMaps();
    newXhist = transform*obj.Xhist;
    if obj.knownChiContextEnabled
        obj.knownChiHistory = newXhist(obj.knownChiIndex,:);
    end
    newZ0 = obj.z0;
    for node = 1:obj.Ne+1
        newZ0(index.x{node}) = transform*obj.z0(index.x{node});
    end
    newPe = inverse.'*obj.Pe*inverse;
    newPe = 0.5*(newPe+newPe.');
    delta = max(abs(10*newReference),1);
    newLower = newReference-delta;
    newUpper = newReference+delta;
    obj.xhat = transform*obj.xhat;
    if ~isempty(obj.xlast)
        obj.xlast = transform*obj.xlast;
    end
    obj.Xhist = newXhist;
    obj.z0 = newZ0;
    obj.Pe = newPe;
    obj.xL = newLower;
    obj.xU = newUpper;
    if obj.reciprocalUnboundedAerodynamicLagStates
        qGamIndex = obj.model.idx.qGam(:);
        obj.xL(qGamIndex) = -inf;
        obj.xU(qGamIndex) = inf;
    end
    obj.stateModelTrim = newReference;
    obj.measurementContract = newContract;
    obj.H = speye(numel(obj.z0));
    obj.Sprev = [eye(obj.nx),zeros(obj.nx,obj.nw+obj.nu)];
    obj.preparedFeedback = struct('valid',false);
    if obj.scheduledPackageHistoryEnabled
        for interval = 1:numel(obj.scheduledPackageHistory)
            entry = obj.scheduledPackageHistory{interval};
            if isempty(entry)
                continue
            end
            assert(isstruct(entry) && isfield(entry,'stateToActive') && ...
                isequal(size(entry.stateToActive),[obj.nx,obj.nx]), ...
                'nMHEv2:ScheduledPackageHistoryTransport', ...
                'The retained package has no valid active-chart map.');
            entry.stateToActive = transform*entry.stateToActive;
            entry.stateFromActive = entry.stateToActive\eye(obj.nx);
            assert(all(isfinite(entry.stateToActive),'all') && ...
                all(isfinite(entry.stateFromActive),'all'), ...
                'nMHEv2:ScheduledPackageHistoryTransport', ...
                'A scheduled package history map became nonfinite.');
            obj.scheduledPackageHistory{interval} = entry;
        end
        obj.scheduledPackageHistoryReason = "transported_to_active_chart";
        obj.incomingPackageStateToActive = ...
            transform*obj.incomingPackageStateToActive;
    end
    info = struct('condition',cond(transform), ...
        'arrivalWeightSymmetryError',norm(obj.Pe-obj.Pe.','fro'), ...
        'measurementRows',obj.ny,'accepted',true);
end
%----------------------------------------------------------------------
function info = prepareNextEstimate(obj,u_expected)
%PREPARENEXTESTIMATE Audit-only continuity preparation for the next sample.
% The future measurement is deliberately excluded; feedback retains the
% actual measurement objective and uses this cache only after exact checks.
    info = struct('enabled',obj.preparationFeedbackAudit, ...
        'prepared',false,'discarded',false,'seconds',0, ...
        'reason',"disabled");
    obj.preparedFeedback = struct('valid',false);
    if ~obj.preparationFeedbackAudit
        return
    end
    u_expected = u_expected(:);
    assert(numel(u_expected) == obj.nu && all(isfinite(u_expected)), ...
        'nMHEv2:PreparationInput', ...
        'Prepared applied input must be a finite trim-relative command.');
    if obj.firstFullSolve || obj.k < obj.Ne + 1
        info.reason = "warmup";
        return
    end
    obj.synchronizeTimingFromModel();
    uForecast = obj.Uhist;
    if obj.Ne > 1
        uForecast(:,1:end-1) = uForecast(:,2:end);
    end
    uForecast(:,end) = u_expected;
    idx = obj.buildIndexMaps();
    zForecast = obj.z0;
    if obj.Ne > 1
        zForecast(idx.w{end}) = zForecast(idx.w{end-1});
    end
    nlp = obj.assembleWindow(uForecast);
    timer = tic;
    prepared = obj.sqpSolver.prepareInitialConstraints( ...
        zForecast,nlp.lb,nlp.ub,nlp.nonl);
    info.seconds = toc(timer);
    obj.preparedFeedback = struct('valid',true, ...
        'expectedInput',u_expected,'z',zForecast, ...
        'modelL',obj.model.L,'modelParConst',obj.model.parConst, ...
        'prepared',prepared, ...
        'preparationSeconds',info.seconds);
    info.prepared = true;
    info.reason = "prepared";
end
%----------------------------------------------------------------------
function [xhat,what,info] = estimate(obj,z_k,u_prev,t_k,varargin)
% Real-time entry point (called every Ts)
%----------------------------------------------------------------------
    obj.synchronizeTimingFromModel();
    z_k = z_k(:);
    u_prev = u_prev(:);
    [measurementContext,knownChi,intervalPackage] = ...
        obj.parseEstimateContexts(varargin{:});
    assert(numel(z_k) == obj.ny && all(isfinite(z_k)), ...
        'nMHEv2:Measurement', ...
        'z_k must contain %d finite measurements.',obj.ny);
    assert(numel(u_prev) == obj.nu && all(isfinite(u_prev)), ...
        'nMHEv2:AppliedInput', ...
        'u_prev must contain %d finite trim-relative applied inputs.',obj.nu);
    assert(isscalar(t_k) && isfinite(t_k), ...
        'nMHEv2:Time','t_k must be a finite scalar.');
    obj.xlast = obj.xhat;
    % ---- 1. roll/append histories every call ----------------------------
    obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
    obj.Yhist(:,1:end-1) = obj.Yhist(:,2:end);
    obj.Whist(:,1:end-1) = obj.Whist(:,2:end);
    obj.Uhist(:,1:end-1) = obj.Uhist(:,2:end);
    obj.measurementContextHistory(:,1:end-1) = ...
        obj.measurementContextHistory(:,2:end);
    obj.knownChiHistory(:,1:end-1) = obj.knownChiHistory(:,2:end);
    if obj.reciprocalCondensedEnabled
        obj.reciprocalLatentHistory(:,1:end-1) = ...
            obj.reciprocalLatentHistory(:,2:end);
        obj.reciprocalLatentHistory(:,end) = ...
            obj.reciprocalLatentHistory(:,end-1);
    end
    if obj.scheduledPackageHistoryEnabled
        obj.scheduledPackageHistory(1:end-1) = ...
            obj.scheduledPackageHistory(2:end);
    end
    if obj.reciprocalProviderEnabled
        assert(obj.reciprocalPendingContextValid, ...
            'nMHEv2:ReciprocalCompletedIntervalMissing', ...
            ['Every estimator sample requires the exact completed-interval ', ...
             'reciprocal context before history is advanced.']);
        obj.reciprocalContextHistory(1:end-1) = ...
            obj.reciprocalContextHistory(2:end);
        obj.reciprocalContextHistory{end} = ...
            obj.reciprocalPendingContext;
        obj.reciprocalPendingContextValid = false;
    end
    
    obj.Yhist(:,end) = z_k;
    obj.Uhist(:,end) = u_prev;
    obj.measurementContextHistory(:,end) = measurementContext;
    obj.knownChiHistory(:,end) = knownChi;
    if obj.scheduledPackageHistoryEnabled
        obj.scheduledPackageHistory{end} = ...
            obj.makeScheduledPackageHistoryEntry(intervalPackage);
        obj.incomingPackageStateToActive = eye(obj.nx);
        obj.scheduledPackageHistoryValid = all( ...
            ~cellfun(@isempty,obj.scheduledPackageHistory));
        if obj.scheduledPackageHistoryValid
            obj.scheduledPackageHistoryReason = "complete";
        else
            obj.scheduledPackageHistoryReason = "awaiting_completed_intervals";
        end
    end
    
    % obj.k = min(obj.k + 1, obj.Ne);
    obj.k = min(obj.k + 1, obj.Ne+1);
    
    % ---- 2. until horizon is sufficiently populated ---------------------
    % if obj.k < obj.Ne
    if obj.k < obj.Ne+1
        if obj.knownChiContextEnabled
            obj.xhat(obj.knownChiIndex) = knownChi;
        end
        obj.Xhist(:,end) = obj.xhat;
        if obj.reciprocalCondensedEnabled
            % Keep prescribed reciprocal memory causal while the native
            % shooting horizon is being populated. At exact trim this is
            % the zero-departure invariant; away from trim it prevents the
            % first optimization from restarting the latent dynamics.
            augmentedWarmup = obj.reciprocalProvider.propagateInterval( ...
                [obj.Xhist(:,end-1); ...
                 obj.reciprocalLatentHistory(:,end-1)], ...
                obj.Whist(:,end), ...
                obj.reciprocalContextForInterval(obj.Ne),false);
            obj.reciprocalLatentHistory(:,end) = ...
                augmentedWarmup(obj.reciprocalLatentIndex);
        end
    
        xhat = obj.outputEstimate();
        what = obj.what;
        idx = obj.buildIndexMaps();

        % The optimizer warm start retains the complete internal state. The
        % public estimate may intentionally expose only the native 74 states.
        obj.z0(idx.x{obj.k}) = obj.xhat;
        info = struct( ...
            'cost',0, ...
            'continuity',0, ...
            'exitflag',0, ...
            'wHorizon',zeros(obj.Ne,obj.nw), ...
            'solveAttempted',false, ...
            'accepted',false, ...
            'fallback',false, ...
            'status',"warmup", ...
            'initializer',obj.coldStartInfo("warmup"), ...
            'realtimeRti',obj.emptyRealtimeRtiInfo());
    
        return
    end

    idx = obj.buildIndexMaps();
    if obj.scheduledPackageHistoryEnabled
        assert(obj.scheduledPackageHistoryValid, ...
            'nMHEv2:ScheduledPackageHistoryIncomplete', ...
            ['A complete source-bound package history is required before ', ...
             'a scheduled-horizon nMHE solve.']);
    end
    % 
    % I think this stays commented out
    % So do not overwrite an internal shooting node before the solve.
    % The shifted solution is already stored in obj.z0 from the previous call.
    % obj.z0(idx.x{obj.k}) = obj.xhat;

    % % current arrival state -------------------------
    % obj.z0(idx.x{1}) = obj.xhat;
    % Nint=round(obj.Ts/obj.dt);
    % 
    % % forward-predict the horizon with zero gust ----
    % x_tmp = obj.xhat;  S_tmp = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
    %     for j = 1:obj.Ne
    %         uj = obj.Uhist(:,j);  wj = obj.z0(idx.w{j});
    %         for m = 1:Nint
    %             [x_tmp,S_tmp] = obj.model.step(x_tmp, uj, wj, S_tmp, true);
    %         end
    %         obj.z0(idx.x{j+1}) = x_tmp;
    %     end
    % obj.Sprev = S_tmp;      

    % Build a solve-local candidate. Persistent z0 is updated only after
    % an accepted solve.
    zSolve0 = obj.z0;
    zSolve0 = obj.bindKnownChiNodes(zSolve0);
    initializerInfo = obj.coldStartInfo("shifted_warm_start");

    % Initial disturbance guess ---------------------
    if obj.firstFullSolve
        switch obj.coldStartStrategy
            case "trim"
                zSolve0([idx.w{:}]) = 0;
                initializerInfo = obj.coldStartInfo("trim");
                initializerInfo.accepted = true;
            case "linear_batch"
                initializerTimer = tic;
                [zLinear,initializerInfo] = ...
                    obj.buildLinearBatchColdStart();
                initializerInfo.seconds = toc(initializerTimer);
                if initializerInfo.accepted
                    zSolve0 = zLinear;
                else
                    zSolve0([idx.w{:}]) = 0;
                    initializerInfo.fallbackToTrim = true;
                    warning('nMHEv2:ColdStartFailure', ...
                        ['The linear-batch cold start was rejected ', ...
                         '(exitflag %g, state-bound violation %.3e, ', ...
                         'disturbance-bound violation %.3e). ', ...
                         'Using the legacy trim horizon. %s'], ...
                        initializerInfo.exitflag, ...
                        initializerInfo.stateBoundViolationInf, ...
                        initializerInfo.disturbanceBoundViolationInf, ...
                        initializerInfo.message);
                end
            otherwise
                error('nMHEv2:ColdStartStrategy', ...
                    'Unhandled cold-start strategy "%s".', ...
                    obj.coldStartStrategy);
        end
    else
        % Extend the shifted disturbance history using the selected
        % source-bound audit policy. The established default remains a
        % held terminal disturbance.
        zSolve0 = obj.applyDisturbanceWarmStart(zSolve0,idx);
    end
    sharedValueHorizonData = struct();
    if obj.preparedHorizonDataReuseActive
        sharedValueHorizonData = ...
            obj.buildNativeEstimatorValueHorizonData();
    end
    if obj.reciprocalProviderEnabled
        % Rebuild every shooting node from the retained arrival and the
        % actual completed interval histories. This changes only the solver
        % initial guess; measurements and decision variables remain free.
        zSolve0 = obj.buildReciprocalFeasibleWarmStart( ...
            zSolve0,sharedValueHorizonData);
        initializerInfo.strategy = "reciprocal_feasible_rollout";
        initializerInfo.attempted = true;
        initializerInfo.accepted = true;
    end
    zSolve0 = obj.bindKnownChiNodes(zSolve0);
    % % ---- 4. assemble NLP & solve ------------------------------------
    % nlp = obj.assembleWindow();
    % [p,fval,exitflag] = fmincon(nlp.cost,obj.z0,...
    %                             [],[],[],[],nlp.lb,nlp.ub,...
    %                             nlp.nonl,obj.solverOpts);
    %======================================================================
    % Assemble and solve nMHE multiple-shooting NLP
    %======================================================================
    nlp = obj.assembleWindow([],sharedValueHorizonData);
    optimizerInitial = zSolve0;
    realtimeRtiCandidate = zeros(0,1);
    realtimeRtiInfo = obj.emptyRealtimeRtiInfo();
    if obj.realtimeRtiEnabled
        try
            [realtimeRtiCandidate,realtimeRtiInfo] = ...
                obj.buildRealtimeRtiSeed(nlp,zSolve0);
            if realtimeRtiInfo.qualified
                optimizerInitial = realtimeRtiCandidate;
            else
                realtimeRtiInfo.fallbackToFullInitial = true;
            end
        catch realtimeRtiException
            realtimeRtiInfo = obj.emptyRealtimeRtiInfo();
            realtimeRtiInfo.attempted = true;
            realtimeRtiInfo.fallbackToFullInitial = true;
            realtimeRtiInfo.message = ...
                "RTI preparation failed closed: " + ...
                string(realtimeRtiException.message);
            realtimeRtiInfo.identifier = ...
                string(realtimeRtiException.identifier);
        end
    end
    [preparedInitial,preparationInfo] = ...
        obj.consumePreparedFeedback(u_prev,zSolve0);
    candidateQuality = struct('enabled',obj.candidateQualityDiagnostics, ...
        'initialCost',NaN,'initialConstraintViolationInf',NaN, ...
        'initialFeasible',false,'candidateToInitialCostRatio',NaN, ...
        'candidateCostIncrease',NaN);
    if obj.candidateQualityDiagnostics
        [candidateQuality.initialCost,~] = nlp.cost(zSolve0);
        [cInitial,ceqInitial] = nlp.nonl(zSolve0);
        candidateQuality.initialConstraintViolationInf = ...
            obj.constraintViolation( ...
                cInitial,ceqInitial,zSolve0,nlp.lb,nlp.ub);
        candidateQuality.initialFeasible = ...
            isfinite(candidateQuality.initialCost) && ...
            candidateQuality.initialConstraintViolationInf <= ...
            obj.constraintTolerance();
    end
    
    if obj.realtimeRtiSolutionOwnerEnabled
        if realtimeRtiInfo.qualified
            p = realtimeRtiCandidate;
            fval = realtimeRtiInfo.objective;
            exitflag = 1;
            output = struct( ...
                'message',"Qualified condensed RTI owned the estimate.", ...
                'iterations',realtimeRtiInfo.completedIterations, ...
                'constrviolation',realtimeRtiInfo.constraintViolationInf, ...
                'firstorderopt',obj.realtimeRtiFirstOrderMetric( ...
                    realtimeRtiInfo));
            lambda = obj.emptySolverMultipliers();
        elseif obj.realtimeRtiFullCorrectionFallback.enabled
            correction = obj.solveFullNlpAudit(nlp,zSolve0);
            p = correction.decision;
            fval = correction.cost;
            exitflag = correction.exitflag;
            output = correction.output;
            lambda = correction.lambda;
            realtimeRtiInfo.fullCorrection = correction;
        else
            p = realtimeRtiCandidate;
            fval = inf;
            exitflag = -2;
            output = struct('message', ...
                "Condensed RTI was rejected; holding the last estimate.");
            lambda = obj.emptySolverMultipliers();
        end
    else
    switch obj.solverName
    
        case "fmincon"
    
            [p,fval,exitflag,output,lambda] = fmincon( ...
                nlp.cost,optimizerInitial, ...
                [],[],[],[], ...
                nlp.lb,nlp.ub, ...
                nlp.nonl,obj.solverOpts);
            if obj.fminconIterationDiagnostics
                output.iterationTrace = obj.fminconIterationTrace;
            end
    
        case "custom_sqp"
    
            if obj.debug && ~obj.sqpCheckDone && ...
                    obj.sqpSolver.options.CheckGradientsOnce

                % Global directional check.
                obj.sqpSolver.checkGradients( ...
                    nlp.cost,nlp.nonl,zSolve0,nlp.lb,nlp.ub);

                % MHE-specific block check: separates X and W sensitivity errors.
                obj.localCheckMHEEqualityGradientBlocks( ...
                    nlp,zSolve0,nlp.lb,nlp.ub);

                obj.sqpCheckDone = true;
            end
    
            [p,fval,exitflag,output,lambda] = obj.sqpSolver.solve( ...
                nlp.cost,optimizerInitial,nlp.lb,nlp.ub,nlp.nonl,preparedInitial);
    
        otherwise
            error('nMHEv2:Solver','Unhandled solverName = "%s".', obj.solverName);
    end
    end

    % ---- 5. classify candidate before persistent-state updates -------
    candidateFinite = ~isempty(p) && numel(p) == numel(zSolve0) && ...
        all(isfinite(p)) && isfinite(fval);
    continuityNorm = inf;
    constraintViolationInf = inf;
    useTerminalStmCondensation = ...
        obj.terminalStmCondensationActive && ...
        obj.nativeValueHorizonActive && ...
        obj.reciprocalCondensedEnabled && ...
        obj.realtimeRtiSolutionOwnerEnabled;
    if candidateFinite
        if useTerminalStmCondensation
            % E3 supplies the exact horizon values. The identical terminal
            % STM is refreshed once after acceptance instead of carrying
            % unused Jacobians through all eight intervals here.
            [cCandidate,ceqCandidate] = nlp.nonl(p);
        else
            % Refresh the exact final sensitivity/STM at the classified
            % candidate. This remains required when merit trials use the
            % audit-only state-propagation path.
            [cCandidate,ceqCandidate,~,~] = nlp.nonl(p);
        end
        candidateFinite = all(isfinite(cCandidate)) && ...
            all(isfinite(ceqCandidate));
        if candidateFinite
            continuityNorm = norm(ceqCandidate,2);
            equalityViolation = norm(ceqCandidate,inf);
            if isempty(cCandidate)
                inequalityViolation = 0;
            else
                inequalityViolation = max(max(cCandidate),0);
            end
            constraintViolationInf = max( ...
                equalityViolation,inequalityViolation);
            if obj.reciprocalProviderEnabled
                constraintViolationInf = obj.constraintViolation( ...
                    cCandidate,ceqCandidate,p,nlp.lb,nlp.ub);
            end
        end
    end
    candidateAccepted = exitflag > 0 && candidateFinite && ...
        constraintViolationInf <= obj.constraintTolerance();

    info.cost = fval;
    if obj.realtimeRtiSolutionOwnerEnabled
        info.solverUsed = "condensed_rti";
    else
        info.solverUsed = obj.solverName;
    end
    info.exitflag = exitflag;
    info.output = output;
    info.lambda = lambda;
    info.continuity = continuityNorm;
    info.constraintViolationInf = constraintViolationInf;
    info.solveAttempted = true;
    info.accepted = candidateAccepted;
    info.initializer = initializerInfo;
    info.preparation = preparationInfo;
    realtimeRtiInfo.candidateValueOnlyValidation = ...
        useTerminalStmCondensation;
    realtimeRtiInfo.terminalSensitivityRefreshIntervals = 0;
    realtimeRtiInfo.latentHistoryCondensationApplied = false;
    info.realtimeRti = realtimeRtiInfo;
    if obj.candidateQualityDiagnostics
        candidateQuality.candidateToInitialCostRatio = fval/max( ...
            abs(candidateQuality.initialCost),eps);
        candidateQuality.candidateCostIncrease = ...
            fval-candidateQuality.initialCost;
    end
    info.candidateQuality = candidateQuality;

    if ~candidateAccepted
        xhat = obj.outputEstimate();
        what = obj.what;
        info.wHorizon = reshape(obj.what,obj.nw,[])';
        info.fallback = true;
        info.status = "fallback_last_accepted";
        warning('nMHEv2:SolveFailure', ...
            ['nMHE solver "%s" rejected exitflag %d with constraint ', ...
             'violation %.3e. Holding the last accepted estimate. %s'], ...
            obj.solverName,exitflag,constraintViolationInf, ...
            obj.outputMessage(output));
        if obj.debug
            obj.debugPlots(t_k,info);
        end
        return
    end

    % ---- 6. accept and shift the finite successful solution ----------
    Xsol = reshape(p([idx.x{:}]),obj.nx,obj.Ne+1);
    Wsol = reshape(p([idx.w{:}]),obj.nw,obj.Ne);

    if obj.reciprocalCondensedEnabled
        acceptedReplayEndpoints = zeros(0,0);
        acceptedReplayInfo = struct('enabled',obj.acceptedReplayReuseActive, ...
            'cacheHit',false,'valueReplayCacheHits',0, ...
            'valueReplayCacheMisses',0);
        if obj.acceptedReplayReuseActive
            [acceptedReplayEndpoints,acceptedReplayInfo] = ...
                nlp.getAcceptedReplay(p);
            assert(acceptedReplayInfo.cacheHit, ...
                'nMHEv2:AcceptedReplayReuseMiss', ...
                ['The accepted estimator decision was not the exact ', ...
                 'decision from the final value replay.']);
        end
        obj.reciprocalLatentHistory = ...
            obj.rebuildReciprocalLatentHistory( ...
                Xsol,Wsol,sharedValueHorizonData,acceptedReplayEndpoints);
        info.realtimeRti.latentHistoryCondensationApplied = ...
            obj.latentHistoryCondensationActive;
        info.realtimeRti.acceptedReplayReuseApplied = ...
            acceptedReplayInfo.cacheHit;
        info.realtimeRti.valueReplayCacheHits = ...
            acceptedReplayInfo.valueReplayCacheHits;
        info.realtimeRti.valueReplayCacheMisses = ...
            acceptedReplayInfo.valueReplayCacheMisses;
    end

    obj.Xhist = Xsol;
    obj.Whist = Wsol;
    obj.xhat = Xsol(:,end);
    obj.what = Wsol(:);
    obj.xlast = Xsol(:,1);
    obj.z0 = obj.shiftGuess(p);
    if useTerminalStmCondensation
        obj.Sprev = obj.refreshReciprocalTerminalStmExact( ...
            Xsol,Wsol,obj.reciprocalLatentHistory);
        info.realtimeRti.terminalSensitivityRefreshIntervals = 1;
    else
        obj.Sprev = nlp.getSTM();
    end
    obj.firstFullSolve = false;

    xhat = obj.outputEstimate();
    what = obj.what;
    info.wHorizon = Wsol.';
    info.fallback = false;
    info.status = "accepted";
    
    if isstruct(output) && isfield(output,'constrviolation')
        info.constrviolation = output.constrviolation;
    end
    
    if isstruct(output) && isfield(output,'firstorderopt')
        info.firstorderopt = output.firstorderopt;
    end
    
    if isstruct(output) && isfield(output,'stepsize')
        info.stepsize = output.stepsize;
    end
    
    if isstruct(output) && isfield(output,'qpExitflag')
        info.qpExitflag = output.qpExitflag;
    end
    if obj.debug  
        obj.debugPlots(t_k,info);  
        % obj.debugPlots2(t_k,info);  
    
    end
    % 
    % wH = info.wHorizon(:);
    % 
    % fprintf('[MHE w check] t = %.4f | min(w)=%.4e | max(w)=%.4e | wL=%.4e | wU=%.4e | hitL=%d | hitU=%d\n', ...
    % t_k, min(wH), max(wH), obj.wL(1), obj.wU(1), ...
    % any(wH <= obj.wL(1) + 1e-6), ...
    % any(wH >= obj.wU(1) - 1e-6));
end
%----------------------------------------------------------------------
function audit = evaluateKnownChiConstraintAudit(obj,z)
%EVALUATEKNOWNCHICONSTRAINTAUDIT Read-only derivative qualification hook.
    assert(obj.knownChiContextEnabled, ...
        'nMHEv2:KnownChiAuditDisabled', ...
        'The known-chi derivative audit requires its explicit request.');
    if nargin < 2 || isempty(z)
        z = obj.bindKnownChiNodes(obj.z0);
    else
        z = obj.bindKnownChiNodes(z(:));
    end
    nlp = obj.assembleWindow();
    [cost,gradient,hessian] = nlp.cost(z);
    [inequality,equality,gradientInequality,gradientEquality] = ...
        nlp.nonl(z);
    audit = struct('z',z,'cost',cost,'gradient',gradient, ...
        'hessian',hessian,'inequality',inequality, ...
        'equality',equality, ...
        'gradientInequality',gradientInequality, ...
        'gradientEquality',gradientEquality, ...
        'lowerBound',nlp.lb,'upperBound',nlp.ub, ...
        'knownChiHistory',obj.knownChiHistory);
end
%----------------------------------------------------------------------
function setReciprocalCompletedIntervalContext(obj,context)
%SETRECIPROCALCOMPLETEDINTERVALCONTEXT Supply causal plant-substep history.
    assert(obj.reciprocalProviderEnabled, ...
        'nMHEv2:ReciprocalProviderDisabled', ...
        'The reciprocal interval context requires its approved audit selector.');
    context = obj.prepareScheduledReciprocalContext(context);
    % The provider owns the complete validation contract. A zero-disturbance
    % dry propagation validates dimensions, finiteness, and command ownership
    % without reading or changing estimator state.
    obj.reciprocalProvider.propagateInterval( ...
        obj.reciprocalValidationState(),0,context,false);
    obj.reciprocalPendingContext = context;
    obj.reciprocalPendingContextValid = true;
end
%----------------------------------------------------------------------
function audit = evaluateReciprocalProviderConstraintAudit(obj,z)
%EVALUATERECIPROCALPROVIDERCONSTRAINTAUDIT Read-only NLP derivative hook.
    assert(obj.reciprocalProviderEnabled, ...
        'nMHEv2:ReciprocalProviderAuditDisabled', ...
        'The reciprocal constraint audit requires its approved selector.');
    if nargin < 2 || isempty(z)
        z = obj.z0;
    else
        z = z(:);
    end
    nlp = obj.assembleWindow();
    [cost,gradient,hessian] = nlp.cost(z);
    [inequality,equality,gradientInequality,gradientEquality] = ...
        nlp.nonl(z);
    audit = struct('z',z,'cost',cost,'gradient',gradient, ...
        'hessian',hessian,'inequality',inequality, ...
        'equality',equality, ...
        'gradientInequality',gradientInequality, ...
        'gradientEquality',gradientEquality, ...
        'lowerBound',nlp.lb,'upperBound',nlp.ub, ...
        'reciprocalContextHistory',{obj.reciprocalContextHistory}, ...
        'nativeStateCount',obj.nativeStateCount, ...
        'latentStateCount',numel(obj.reciprocalLatentIndex));
end
%----------------------------------------------------------------------
function audit = evaluateReciprocalDisturbanceWarmStartAudit(obj,z)
%EVALUATERECIPROCALDISTURBANCEWARMSTARTAUDIT Inspect the causal seed only.
    assert(obj.reciprocalProviderEnabled, ...
        'nMHEv2:ReciprocalProviderAuditDisabled', ...
        'The reciprocal warm-start audit requires its approved selector.');
    index = obj.buildIndexMaps();
    z = z(:);
    assert(numel(z)==numel(obj.z0) && all(isfinite(z)), ...
        'nMHEv2:ReciprocalWarmStartAuditDecision', ...
        'The reciprocal warm-start audit decision is invalid.');
    if obj.Ne > 2
        predecessor = z(index.w{end-2});
    else
        predecessor = z(index.w{end-1});
    end
    updated = obj.applyDisturbanceWarmStart(z,index);
    audit = struct( ...
        'strategy',obj.reciprocalDisturbanceWarmStart, ...
        'previousDisturbance',z(index.w{end-1}), ...
        'predecessorDisturbance',predecessor, ...
        'terminalDisturbance',updated(index.w{end}), ...
        'lowerBound',obj.wL,'upperBound',obj.wU);
end
%----------------------------------------------------------------------
function state = reciprocalInternalMemoryEndpoint(obj)
%RECIPROCALINTERNALMEMORYENDPOINT Return accepted controller-model memory.
    assert(obj.reciprocalProviderEnabled && ...
        obj.reciprocalCondensedEnabled && ...
        isequal(size(obj.reciprocalLatentHistory), ...
            [obj.reciprocalProvider.hiddenStateCount,obj.Ne+1]) && ...
        all(isfinite(obj.reciprocalLatentHistory),'all'), ...
        'nMHEv2:ReciprocalInternalMemoryEndpoint', ...
        ['The reciprocal internal-memory endpoint requires an enabled, ', ...
         'valid condensed estimator.']);
    state = obj.reciprocalLatentHistory(:,end);
end
%----------------------------------------------------------------------
function audit = reciprocalInternalMemoryAudit(obj)
%RECIPROCALINTERNALMEMORYAUDIT Return compact causal-memory diagnostics.
    assert(obj.reciprocalProviderEnabled && ...
        obj.reciprocalCondensedEnabled, ...
        'nMHEv2:ReciprocalInternalMemoryAudit', ...
        ['The reciprocal-memory audit requires the enabled condensed ', ...
         'provider formulation.']);
    context = obj.reciprocalContextForInterval(obj.Ne);
    audit = struct( ...
        'endpoint',obj.reciprocalLatentHistory(:,end), ...
        'trim',obj.reciprocalLatentTrim(:), ...
        'modelTrim',obj.uModelTrim(:), ...
        'wingTotal',context.wingTotal, ...
        'wingIncrement',context.wingIncrement, ...
        'elevatorIncrement',context.elevatorIncrement, ...
        'thrustIncrement',context.thrustIncrement, ...
        'rigidState',context.rigidState, ...
        'scheduledReciprocal',{context.scheduledReciprocal});
end
%----------------------------------------------------------------------
function applyScheduledReciprocalPacket(obj,packet)
%APPLYSCHEDULEDRECIPROCALPACKET Refresh query maps without memory reset.
    assert(obj.reciprocalProviderEnabled && ...
        isa(obj.reciprocalProvider, ...
            'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider'), ...
        'nMHEv2:ScheduledReciprocalProviderDisabled', ...
        'The scheduled reciprocal provider is not active.');
    obj.reciprocalProvider.applyScheduledPacket(obj.model,packet);
end
%----------------------------------------------------------------------
function audit = evaluateRealtimeRtiSeedAudit(obj,z,runFullCorrection)
%EVALUATEREALTIMERTISEEDAUDIT Read-only condensed-estimator seed hook.
    assert(obj.realtimeRtiEnabled, ...
        'nMHEv2:RealtimeRtiAuditDisabled', ...
        'The condensed nMHE RTI audit requires its approved selector.');
    if nargin < 2 || isempty(z)
        z = obj.realtimeRtiAuditInitial();
    else
        z = z(:);
    end
    if nargin < 3 || isempty(runFullCorrection)
        runFullCorrection = false;
    end
    assert(isscalar(runFullCorrection) && ...
        (islogical(runFullCorrection) || ...
        (isnumeric(runFullCorrection) && isfinite(runFullCorrection) && ...
        ismember(runFullCorrection,[0 1]))), ...
        'nMHEv2:RealtimeRtiAuditCorrection', ...
        'The full-correction selector must be logical.');
    assert(numel(z)==numel(obj.z0) && all(isfinite(z)), ...
        'nMHEv2:RealtimeRtiAuditDecision', ...
        'The condensed nMHE RTI audit decision is invalid.');
    nlp = obj.assembleWindow();
    [candidate,info] = obj.buildRealtimeRtiSeed(nlp,z);
    audit = struct('initial',z,'candidate',candidate,'info',info);
    if logical(runFullCorrection)
        audit.fullCorrection = struct( ...
            'native',obj.solveFullNlpAudit(nlp,z), ...
            'rtiWarm',obj.solveFullNlpAudit(nlp,candidate));
    end
end
%----------------------------------------------------------------------
function configureNativeReducedHorizonCheckpointAudit(obj,request)
%CONFIGURENATIVEREDUCEDHORIZONCHECKPOINTAUDIT Bind a saved Case-B owner.
    cfg = struct('ctrl',struct( ...
        'nativeReducedHorizonRtiAudit',request));
    obj.configureNativeReducedHorizonRtiAudit(cfg);
end
%----------------------------------------------------------------------
function configureHighLeverageRuntimeCheckpointAudit(obj,request)
%CONFIGUREHIGHLEVERAGERUNTIMECHECKPOINTAUDIT Bind an approved saved E3 owner.
    cfg = struct('ctrl',struct('highLeverageRuntimeAudit',request));
    obj.configureHighLeverageRuntimeAudit(cfg);
end
%----------------------------------------------------------------------
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
        obj.latentHistoryCondensationActive && ...
        obj.terminalStmCondensationActive && ...
        obj.realtimeRtiSolutionOwnerEnabled && ...
        obj.realtimeRtiIterationCount==1, ...
        'nMHEv2:PreparedHorizonDataReuseCheckpointRequest', ...
        'Exact data reuse requires the qualified saved V73 estimator.');
    obj.preparedHorizonDataReuseRequested = true;
    obj.preparedHorizonDataReuseActive = true;
end
%----------------------------------------------------------------------
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
        obj.latentHistoryCondensationActive && ...
        obj.terminalStmCondensationActive && ...
        obj.realtimeRtiSolutionOwnerEnabled && ...
        obj.realtimeRtiIterationCount==1, ...
        'nMHEv2:AcceptedReplayReuseCheckpointRequest', ...
        'Exact replay reuse requires the qualified saved V73 estimator.');
    obj.acceptedReplayReuseRequested = true;
    obj.acceptedReplayReuseActive = true;
end
%----------------------------------------------------------------------
end  % public methods

    %======================================================================
    methods (Access = private)
    %----------------------------------------------------------------------
    function configureRealtimeRtiAudit(obj,cfg)
    %CONFIGUREREALTIMERTIAUDIT Install the approved default-inactive seed.
        legacyField = 'realtimeNmheCondensedRtiAudit';
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
            'nMHEv2:RealtimeRtiRequestConflict', ...
            'Legacy RTI seeding and RTI solution ownership are exclusive.');
        if ownerPresent
            request = cfg.ctrl.(ownerField);
        else
            request = cfg.ctrl.(legacyField);
        end
        assert(isstruct(request) && isscalar(request) && ...
            isfield(request,'enabled'), 'nMHEv2:RealtimeRtiRequest', ...
            'The nMHE RTI request requires enabled.');
        enabled = request.enabled;
        assert(isscalar(enabled) && (islogical(enabled) || ...
            (isnumeric(enabled) && isfinite(enabled) && ...
            ismember(enabled,[0 1]))), 'nMHEv2:RealtimeRtiRequest', ...
            'The nMHE RTI enabled selector must be logical.');
        if ~logical(enabled)
            return
        end
        if ownerPresent
            required = {'auditOnly','changeId','condensingMode', ...
                'solutionOwner','nmheIterationCount','fallbackPolicy'};
            assert(all(isfield(request,required)) && ...
                ~logical(request.auditOnly) && ...
                string(request.changeId)== ...
                    "phase18c-v17a-condensed-rti-compiled-runtime-owner-promotion-v1" && ...
                string(request.condensingMode)=="full_state" && ...
                string(request.solutionOwner)=="qualified_rti" && ...
                string(request.fallbackPolicy)=="hold_previous", ...
                'nMHEv2:RealtimeRtiRequest', ...
                'The enabled nMHE RTI solution-owner request is unauthorized.');
            iterationCount = double(request.nmheIterationCount);
            obj.realtimeRtiSolutionOwnerEnabled = true;
        else
            required = {'auditOnly','changeId','condensingMode', ...
                'iterationCount','correctionSolver'};
            assert(all(isfield(request,required)) && ...
                logical(request.auditOnly) && ...
                string(request.changeId)== ...
                    "phase18c-v17a-realtime-nmhe-full-condensing-rti-audit-v1" && ...
                string(request.condensingMode)=="full_state" && ...
                string(request.correctionSolver)=="configured_full_nlp", ...
                'nMHEv2:RealtimeRtiRequest', ...
                'The enabled nMHE RTI audit request is not authorized.');
            iterationCount = double(request.iterationCount);
        end
        assert(isscalar(iterationCount) && isfinite(iterationCount) && ...
            iterationCount==fix(iterationCount) && ...
            iterationCount>=1 && iterationCount<=3, ...
            'nMHEv2:RealtimeRtiIterationCount', ...
            'The nMHE RTI audit permits one to three iterations.');
        assert(obj.solverName=="fmincon" && obj.method=="multiple", ...
            'nMHEv2:RealtimeRtiSolver', ...
            ['The condensed nMHE RTI path requires the fmincon-configured ', ...
             'multiple-shooting formulation as its offline reference.']);

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
                'nMHEv2:RealtimeRtiQpMode', ...
                'The condensed nMHE RTI QP mode is unsupported.');
            rtiOptions.CondensedRtiQpMode = qpMode;
        end
        obj.realtimeRtiSolver = AeroFlex.optim.SQPSolver(rtiOptions);
        obj.realtimeRtiEnabled = true;
        obj.realtimeRtiChangeId = string(request.changeId);
        obj.realtimeRtiIterationCount = iterationCount;
        if isfield(cfg.ctrl,'nmheCondensedRtiNonlinearBacktrackingAudit') && ...
                ~isempty(cfg.ctrl.nmheCondensedRtiNonlinearBacktrackingAudit)
            guard = cfg.ctrl.nmheCondensedRtiNonlinearBacktrackingAudit;
            requiredGuard = {'enabled','auditOnly','changeId', ...
                'maxTrials','beta'};
            assert(isstruct(guard) && isscalar(guard) && ...
                all(isfield(guard,requiredGuard)) && logical(guard.enabled) && ...
                logical(guard.auditOnly) && string(guard.changeId)== ...
                "phase18c-v17a-caseb-nmhe-rti-nonlinear-backtracking-audit-v1", ...
                'nMHEv2:RealtimeRtiBacktrackingRequest', ...
                'The nMHE nonlinear-backtracking audit request is invalid.');
            assert(isscalar(guard.maxTrials) && guard.maxTrials==fix(guard.maxTrials) && ...
                guard.maxTrials>=1 && guard.maxTrials<=20 && ...
                isscalar(guard.beta) && isfinite(guard.beta) && ...
                guard.beta>0 && guard.beta<1, ...
                'nMHEv2:RealtimeRtiBacktrackingRequest', ...
                'The nMHE nonlinear-backtracking parameters are invalid.');
            obj.realtimeRtiNonlinearBacktracking = struct( ...
                'enabled',true,'maxTrials',double(guard.maxTrials), ...
                'beta',double(guard.beta));
        end
        if isfield(cfg.ctrl,'nmheCondensedRtiFullCorrectionAudit') && ...
                ~isempty(cfg.ctrl.nmheCondensedRtiFullCorrectionAudit)
            fallback = cfg.ctrl.nmheCondensedRtiFullCorrectionAudit;
            assert(isstruct(fallback) && isscalar(fallback) && ...
                isfield(fallback,'enabled') && isfield(fallback,'auditOnly') && ...
                isfield(fallback,'changeId') && logical(fallback.enabled) && ...
                logical(fallback.auditOnly) && string(fallback.changeId)== ...
                "phase18c-v17a-caseb-nmhe-rti-full-correction-fallback-audit-v1", ...
                'nMHEv2:RealtimeRtiFullCorrectionRequest', ...
                'The nMHE full-correction fallback request is invalid.');
            obj.realtimeRtiFullCorrectionFallback = struct('enabled',true);
        end
    end

    %----------------------------------------------------------------------
    function configureNativeReducedHorizonRtiAudit(obj,cfg)
    %CONFIGURENATIVEREDUCEDHORIZONRTIAUDIT Bind scheduled-only H2 dispatch.
        fieldName = 'nativeReducedHorizonRtiAudit';
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl,fieldName) || ...
                isempty(cfg.ctrl.(fieldName))
            return
        end
        request = cfg.ctrl.(fieldName);
        required = {'enabled','auditOnly','changeId','estimatorEnabled', ...
            'estimatorBinarySha256'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)) && logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1", ...
            'nMHEv2:NativeReducedHorizonRequest', ...
            'The native reduced-horizon estimator request is unauthorized.');
        obj.nativeReducedHorizonRtiRequested = ...
            logical(request.enabled) && logical(request.estimatorEnabled);
        if ~obj.nativeReducedHorizonRtiRequested
            return
        end
        if ~obj.realtimeRtiEnabled || obj.Ne~=8 || obj.nx~=74 || obj.nw~=1 || ...
                ~obj.reciprocalCondensedEnabled || ...
                ~isa(obj.reciprocalProvider, ...
                    'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            obj.nativeReducedHorizonRtiLastFallback = ...
                "unsupported estimator architecture; exact RTI retained";
            return
        end
        functionName = ...
            'AeroFlex_ctrl_scheduledReciprocalEstimatorReducedTangentHorizonAudit_mex';
        kernelPath = string(which(functionName));
        expectedHash = lower(string(request.estimatorBinarySha256));
        if kernelPath=="" || strlength(expectedHash)~=64 || ...
                obj.localFileHash(kernelPath)~=expectedHash
            obj.nativeReducedHorizonRtiLastFallback = ...
                "estimator horizon MEX unavailable or hash mismatch";
            return
        end
        obj.nativeReducedHorizonRtiKernel = str2func(functionName);
        obj.nativeReducedHorizonRtiKernelIdentity = functionName+"|"+ ...
            expectedHash+"|"+string(computer('arch'))+"|"+string(mexext);
        obj.nativeReducedHorizonRtiActive = true;
    end

    %----------------------------------------------------------------------
    function configureHighLeverageRuntimeAudit(obj,cfg)
    %CONFIGUREHIGHLEVERAGERUNTIMEAUDIT Bind the scheduled-only E3 value MEX.
        fieldName = 'highLeverageRuntimeAudit';
        if ~isfield(cfg,'ctrl') || ~isfield(cfg.ctrl,fieldName) || ...
                isempty(cfg.ctrl.(fieldName))
            return
        end
        request = cfg.ctrl.(fieldName);
        required = {'enabled','auditOnly','changeId','e3Enabled', ...
            'estimatorValueBinarySha256'};
        assert(isstruct(request) && isscalar(request) && ...
            all(isfield(request,required)) && logical(request.auditOnly) && ...
            string(request.changeId)== ...
                "phase18c-v17a-casebc-high-leverage-runtime-audit-v1", ...
            'nMHEv2:HighLeverageRuntimeRequest', ...
            'The high-leverage estimator request is unauthorized.');
        obj.terminalStmCondensationRequested = false;
        obj.terminalStmCondensationActive = false;
        obj.latentHistoryCondensationRequested = false;
        obj.latentHistoryCondensationActive = false;
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
        if ~obj.nativeReducedHorizonRtiActive || obj.Ne~=8 || ...
                obj.nx~=74 || obj.nw~=1 || ...
                ~obj.reciprocalCondensedEnabled || ...
                ~isa(obj.reciprocalProvider, ...
                    'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            obj.nativeValueHorizonLastFallback = ...
                "unsupported estimator architecture; H2 replay retained";
            return
        end
        functionName = ...
            'AeroFlex_ctrl_scheduledReciprocalEstimatorValueHorizonAudit_mex';
        kernelPath = string(which(functionName));
        expectedHash = lower(string(request.estimatorValueBinarySha256));
        if kernelPath=="" || strlength(expectedHash)~=64 || ...
                obj.localFileHash(kernelPath)~=expectedHash
            obj.nativeValueHorizonLastFallback = ...
                "estimator value-horizon MEX unavailable or hash mismatch";
            return
        end
        obj.nativeValueHorizonKernel = str2func(functionName);
        obj.nativeValueHorizonKernelIdentity = functionName+"|"+ ...
            expectedHash+"|"+string(computer('arch'))+"|"+string(mexext);
        obj.nativeValueHorizonActive = true;
        if isfield(request,'estimatorCausalRolloutEnabled') && ...
                logical(request.estimatorCausalRolloutEnabled)
            assert(isfield(request,'estimatorCausalRolloutBinarySha256'), ...
                'nMHEv2:NativeCausalRolloutRequest', ...
                'The estimator causal-rollout binary hash is required.');
            obj.nativeCausalRolloutRequested = true;
            causalFunctionName = ...
                'AeroFlex_ctrl_scheduledReciprocalEstimatorCausalRolloutAudit_mex';
            causalKernelPath = string(which(causalFunctionName));
            causalExpectedHash = lower(string( ...
                request.estimatorCausalRolloutBinarySha256));
            if causalKernelPath=="" || strlength(causalExpectedHash)~=64 || ...
                    obj.localFileHash(causalKernelPath)~=causalExpectedHash
                obj.nativeCausalRolloutLastFallback = ...
                    "estimator causal-rollout MEX unavailable or hash mismatch";
            else
                obj.nativeCausalRolloutKernel = str2func(causalFunctionName);
                obj.nativeCausalRolloutKernelIdentity = ...
                    causalFunctionName+"|"+causalExpectedHash+"|"+ ...
                    string(computer('arch'))+"|"+string(mexext);
                obj.nativeCausalRolloutActive = true;
            end
        end
        if isfield(request,'latentHistoryCondensationEnabled') && ...
                logical(request.latentHistoryCondensationEnabled)
            assert(obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled, ...
                'nMHEv2:LatentHistoryCondensationRequest', ...
                ['Exact latent-history condensation requires the qualified ', ...
                 'online condensed-RTI solution owner.']);
            obj.latentHistoryCondensationRequested = true;
            obj.latentHistoryCondensationActive = true;
        end
        if isfield(request,'terminalStmCondensationEnabled') && ...
                logical(request.terminalStmCondensationEnabled)
            assert(obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled, ...
                'nMHEv2:TerminalStmCondensationRequest', ...
                ['Exact terminal-STM condensation requires the qualified ', ...
                 'online condensed-RTI solution owner.']);
            obj.terminalStmCondensationRequested = true;
            obj.terminalStmCondensationActive = true;
        end
        if isfield(request,'preparedHorizonDataReuseEnabled') && ...
                logical(request.preparedHorizonDataReuseEnabled)
            assert(isfield(request,'preparedHorizonDataReuseScope') && ...
                string(request.preparedHorizonDataReuseScope)== ...
                    "saved_state_discriminator" && ...
                obj.nativeCausalRolloutActive && ...
                obj.latentHistoryCondensationActive && ...
                obj.terminalStmCondensationActive && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                isfield(request,'r1Enabled') && ...
                logical(request.r1Enabled) && ...
                isfield(request,'nmheIterationCount') && ...
                double(request.nmheIterationCount)==1, ...
                'nMHEv2:PreparedHorizonDataReuseRequest', ...
                ['Exact prepared-horizon data reuse requires the ', ...
                 'qualified V73 estimator owner.']);
            obj.preparedHorizonDataReuseRequested = true;
            obj.preparedHorizonDataReuseActive = true;
        end
        if isfield(request,'acceptedReplayReuseEnabled') && ...
                logical(request.acceptedReplayReuseEnabled)
            assert(isfield(request,'acceptedReplayReuseScope') && ...
                string(request.acceptedReplayReuseScope)== ...
                    "saved_state_discriminator" && ...
                obj.latentHistoryCondensationActive && ...
                obj.terminalStmCondensationActive && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                isfield(request,'r1Enabled') && ...
                logical(request.r1Enabled) && ...
                isfield(request,'nmheIterationCount') && ...
                double(request.nmheIterationCount)==1, ...
                'nMHEv2:AcceptedReplayReuseRequest', ...
                ['Exact accepted-replay reuse requires the qualified ', ...
                 'V73 estimator owner.']);
            obj.acceptedReplayReuseRequested = true;
            obj.acceptedReplayReuseActive = true;
        end
        if isfield(request,'r1Enabled') && logical(request.r1Enabled)
            assert(isfield(request,'nmheIterationCount') && ...
                isfield(request,'r1Scope') && ...
                string(request.r1Scope)=="saved_state_discriminator" && ...
                obj.realtimeRtiEnabled && ...
                obj.realtimeRtiSolutionOwnerEnabled && ...
                obj.realtimeRtiIterationCount==3 && ...
                isscalar(request.nmheIterationCount) && ...
                double(request.nmheIterationCount)==1, ...
                'nMHEv2:HighLeverageR1Request', ...
                ['The one-correction estimator discriminator requires ', ...
                 'the qualified three-correction H2/E3 owner.']);
            obj.realtimeRtiIterationCount = 1;
        end
    end

    %----------------------------------------------------------------------
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

    %----------------------------------------------------------------------
    function lambda = emptySolverMultipliers(~)
        lambda = struct('eqnonlin',zeros(0,1), ...
            'ineqnonlin',zeros(0,1),'lower',zeros(0,1), ...
            'upper',zeros(0,1));
    end

    %----------------------------------------------------------------------
    function [candidate,info] = buildRealtimeRtiSeed(obj,nlp,zInitial)
    %BUILDREALTIMERTISEED Construct an exact condensed RTI full-NLP seed.
        assert(obj.realtimeRtiEnabled && ~isempty(obj.realtimeRtiSolver), ...
            'nMHEv2:RealtimeRtiDisabled', ...
            'The condensed nMHE RTI seed requires its approved selector.');
        if obj.nativeReducedHorizonRtiActive
            [candidate,info] = obj.buildNativeReducedHorizonRtiSeed( ...
                nlp,zInitial);
            return
        end
        stateCount = obj.nx*(obj.Ne+1);
        eliminatedColumns = obj.nx+(1:obj.nx*obj.Ne);
        request = struct('auditOnly',true, ...
            'transformation',speye(numel(zInitial)), ...
            'eliminatedColumns',eliminatedColumns, ...
            'eliminatedRows',1:(obj.nx*obj.Ne), ...
            'iterationCount',obj.realtimeRtiIterationCount, ...
            'constraintTolerance',1e-8, ...
            'inequalityTolerance',1e-10, ...
            'nonlinearBacktracking',obj.realtimeRtiNonlinearBacktracking);
        assert(stateCount+obj.nw*obj.Ne==numel(zInitial) && ...
            numel(eliminatedColumns)==obj.nx*obj.Ne, ...
            'nMHEv2:RealtimeRtiDimension', ...
            'The condensed nMHE decision ordering is incompatible.');
        [candidate,solverInfo] = obj.realtimeRtiSolver.solveCondensedRti( ...
            nlp.cost,zInitial,nlp.lb,nlp.ub,nlp.nonl,nlp.nonl,request);
        info = solverInfo;
        info.enabled = true;
        info.attempted = true;
        info.fallbackToFullInitial = false;
        info.changeId = obj.realtimeRtiChangeId;
        info.condensingMode = "full_state";
        info.correctionSolver = obj.solverName;
        info.fullDecisionCount = numel(zInitial);
        info.eliminatedStateCount = numel(eliminatedColumns);
        info.retainedDecisionCount = ...
            numel(zInitial)-numel(eliminatedColumns);
        info.nativeReducedHorizon = false;
        info.nativeKernelIdentity = "disabled";
        info.horizonPreparationSeconds = 0;
        info.causalWarmStartCondensationApplied = ...
            obj.nativeCausalRolloutActive;
    end

    %----------------------------------------------------------------------
    function [candidate,info] = buildNativeReducedHorizonRtiSeed( ...
            obj,nlp,zInitial)
    %BUILDNATIVEREDUCEDHORIZONRTISEED Apply exact causal Ne=8 condensing.
        totalTimer = tic;
        preparationTimer = tic;
        if obj.preparedHorizonDataReuseActive
            assert(isfield(nlp,'valueHorizonData'), ...
                'nMHEv2:PreparedHorizonDataReuseMissing', ...
                'The shared estimator horizon data is unavailable.');
            data = nlp.valueHorizonData;
        else
            data = obj.buildNativeEstimatorValueHorizonData();
        end
        horizonPreparationSeconds = toc(preparationTimer);
        packetFunction = @(z)obj.nativeReducedEstimatorPacket(z,data);
        request = struct('auditOnly',true, ...
            'changeId', ...
                "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1", ...
            'iterationCount',obj.realtimeRtiIterationCount, ...
            'constraintTolerance',1e-8, ...
            'inequalityTolerance',1e-10, ...
            'nonlinearBacktracking',obj.realtimeRtiNonlinearBacktracking);
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
        info.condensingMode = "native_reduced_horizon";
        info.correctionSolver = obj.solverName;
        info.fullDecisionCount = numel(zInitial);
        info.eliminatedStateCount = obj.nx*obj.Ne;
        info.retainedDecisionCount = obj.nx+obj.nw*obj.Ne;
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

    %----------------------------------------------------------------------
    function packet = nativeReducedEstimatorPacket(obj,z,data)
    %NATIVEREDUCEDESTIMATORPACKET Assemble the exact retained-coordinate map.
        idx = obj.buildIndexMaps();
        nativeNodes = reshape(z(1:obj.nx*(obj.Ne+1)), ...
            obj.nx,obj.Ne+1);
        disturbances = reshape(z(obj.nx*(obj.Ne+1)+1:end),1,obj.Ne);
        [~,~,stateExpansion,stateOffset] = ...
            obj.nativeReducedHorizonRtiKernel( ...
                nativeNodes,data.latentInitial,disturbances, ...
                data.wingTotal,data.wingIncrement,data.elevator, ...
                data.thrust,data.rigid,data.packets);
        reducedCount = obj.nx+obj.nw*obj.Ne;
        expansion = spalloc(numel(z),reducedCount, ...
            obj.nx*reducedCount*(obj.Ne+1)+obj.Ne);
        offset = zeros(numel(z),1);
        for node = 1:obj.Ne+1
            expansion(idx.x{node},:) = stateExpansion(:,:,node);
            offset(idx.x{node}) = stateOffset(:,node);
        end
        disturbanceRows = [idx.w{:}];
        expansion(disturbanceRows,obj.nx+(1:obj.Ne)) = speye(obj.Ne);
        if obj.knownChiContextEnabled
            equality = nativeNodes(obj.knownChiIndex,1)- ...
                obj.knownChiHistory(:,1);
            Aeq = expansion(obj.knownChiIndex,:);
            beq = -equality-offset(obj.knownChiIndex);
        else
            Aeq = sparse(0,reducedCount);
            beq = zeros(0,1);
        end
        packet = struct('expansion',expansion,'offset',offset, ...
            'A',sparse(0,reducedCount),'b',zeros(0,1), ...
            'Aeq',Aeq,'beq',beq);
    end

    %----------------------------------------------------------------------
    function z = realtimeRtiAuditInitial(obj)
    %REALTIMERTIAUDITINITIAL Reproduce the estimator's pre-solve seed owner.
        z = obj.z0;
        index = obj.buildIndexMaps();
        if obj.firstFullSolve
            assert(obj.coldStartStrategy=="trim", ...
                'nMHEv2:RealtimeRtiAuditColdStart', ...
                'The saved-checkpoint RTI audit requires the trim cold start.');
            z([index.w{:}]) = 0;
        else
            z = obj.applyDisturbanceWarmStart(z,index);
        end
        if obj.reciprocalProviderEnabled
            z = obj.buildReciprocalFeasibleWarmStart(z);
        end
        z = obj.bindKnownChiNodes(z);
    end

    %----------------------------------------------------------------------
    function info = emptyRealtimeRtiInfo(obj)
        info = struct('enabled',obj.realtimeRtiEnabled, ...
            'attempted',false,'qualified',false, ...
            'fallbackToFullInitial',false,'message',"not attempted", ...
            'identifier',"",'changeId',obj.realtimeRtiChangeId, ...
            'condensingMode',"disabled", ...
            'correctionSolver',obj.solverName,'iterationCount',0, ...
            'completedIterations',0,'iterations',struct([]), ...
            'nonlinearEqualityInf',nan,'nonlinearInequalityMax',nan, ...
            'boundViolationMax',nan,'constraintViolationInf',nan, ...
            'objective',nan,'valueReplaySeconds',0,'totalSeconds',0, ...
            'fullDecisionCount',numel(obj.z0), ...
            'eliminatedStateCount',0, ...
            'retainedDecisionCount',numel(obj.z0), ...
            'nativeReducedHorizon',false, ...
            'nativeKernelIdentity',"disabled", ...
            'horizonPreparationSeconds',0);
        info.causalWarmStartCondensationApplied = ...
            obj.nativeCausalRolloutActive;
        info.preparedHorizonDataReuseApplied = ...
            obj.preparedHorizonDataReuseActive;
        info.horizonDataBuildsPerSample = ...
            double(obj.preparedHorizonDataReuseActive);
    end

    %----------------------------------------------------------------------
    function result = solveFullNlpAudit(obj,nlp,initial)
    %SOLVEFULLNLPAUDIT Compare a seed without changing estimator history.
        assert(obj.solverName=="fmincon" && all(isfinite(initial)), ...
            'nMHEv2:RealtimeRtiAuditCorrection', ...
            'The read-only correction requires a finite fmincon seed.');
        timer = tic;
        try
            [decision,cost,exitflag,output,lambda] = fmincon( ...
                nlp.cost,initial,[],[],[],[],nlp.lb,nlp.ub, ...
                nlp.nonl,obj.solverOpts);
            failureIdentifier = "";
            failureMessage = "";
        catch solveException
            decision = zeros(0,1);
            cost = inf;
            exitflag = -99;
            output = struct('message',string(solveException.message));
            lambda = struct();
            failureIdentifier = string(solveException.identifier);
            failureMessage = string(solveException.message);
        end
        seconds = toc(timer);
        if isempty(decision) || numel(decision)~=numel(initial) || ...
                any(~isfinite(decision)) || ~isfinite(cost)
            violation = inf;
        else
            [inequality,equality] = nlp.nonl(decision);
            violation = obj.constraintViolation( ...
                inequality,equality,decision,nlp.lb,nlp.ub);
        end
        result = struct('decision',decision,'cost',cost, ...
            'exitflag',exitflag,'output',output,'lambda',lambda, ...
            'seconds',seconds,'constraintViolationInf',violation, ...
            'failureIdentifier',failureIdentifier, ...
            'failureMessage',failureMessage, ...
            'accepted',exitflag>0 && isfinite(violation) && ...
                violation<=obj.constraintTolerance());
    end

    %----------------------------------------------------------------------
    function xTrim = configureReciprocalProvider(obj,cfg,nativeTrimState)
        xTrim = nativeTrimState(:);
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
            'nMHEv2:ReciprocalProviderConflict', ...
            ['Select either the exact-source or scheduled reciprocal ', ...
             'provider, not both.']);
        if hasScheduledRequest
            request = cfg.ctrl. ...
                scheduledReciprocalControllerModelProviderAudit;
            required = {'enabled','auditOnly','changeId','initialPacket', ...
                'initialContext','formulation'};
            assert(isstruct(request) && isscalar(request), ...
                'nMHEv2:ScheduledReciprocalProviderRequest', ...
                'The scheduled reciprocal provider request must be scalar.');
            missing = required(~isfield(request,required));
            assert(isempty(missing), ...
                'nMHEv2:ScheduledReciprocalProviderRequest', ...
                'The scheduled reciprocal provider request lacks: %s.', ...
                strjoin(missing,', '));
            assert(logical(request.auditOnly) && string(request.changeId)== ...
                "phase18c-v17a-casebc-scheduled-reciprocal-runtime-integration-v1", ...
                'nMHEv2:ScheduledReciprocalProviderApproval', ...
                'The scheduled reciprocal provider approval binding changed.');
            obj.reciprocalProviderEnabled = logical(request.enabled);
            if ~obj.reciprocalProviderEnabled
                return
            end
            formulation = lower(string(request.formulation));
            assert(isscalar(formulation) && ...
                formulation=="condensed_internal_memory" && ...
                obj.scheduledPackageHistoryEnabled && ...
                ~obj.knownChiContextEnabled && ...
                obj.coldStartStrategy=="trim", ...
                'nMHEv2:ScheduledReciprocalProviderScope', ...
                ['The scheduled provider requires condensed memory, causal ', ...
                 'scheduled-package history, native chi, and trim startup.']);
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
                if isfield(compiledRequest,'estimatorEnabled')
                    compiledRequest.enabled = logical( ...
                        compiledRequest.enabled) && logical( ...
                        compiledRequest.estimatorEnabled);
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
            obj.reciprocalCondensedEnabled = true;
            obj.reciprocalLatentIndex = obj.nativeStateCount+ ...
                (1:obj.reciprocalProvider.hiddenStateCount);
            augmentedTrim = obj.reciprocalProvider.initialize(nativeTrimState);
            obj.reciprocalLatentTrim = ...
                augmentedTrim(obj.reciprocalLatentIndex);
            obj.reciprocalLatentHistory = repmat( ...
                obj.reciprocalLatentTrim,1,obj.Ne+1);
            obj.reciprocalContextHistory = ...
                repmat({request.initialContext},1,obj.Ne);
            obj.reciprocalPendingContext = request.initialContext;
            obj.reciprocalPendingContextValid = true;
            if isfield(request,'unboundedAerodynamicLagStates') && ...
                    ~isempty(request.unboundedAerodynamicLagStates)
                selector = request.unboundedAerodynamicLagStates;
                assert(isscalar(selector) && (islogical(selector) || ...
                    (isnumeric(selector) && isfinite(selector) && ...
                    ismember(selector,[0 1]))), ...
                    'nMHEv2:ReciprocalAerodynamicLagBounds', ...
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
            'nMHEv2:ReciprocalProviderRequest', ...
            'The reciprocal provider request is incomplete or unauthorized.');
        obj.reciprocalProviderEnabled = logical(request.enabled);
        if ~obj.reciprocalProviderEnabled
            return
        end
        assert(~obj.scheduledPackageHistoryEnabled && ...
            ~obj.knownChiContextEnabled && ...
            obj.coldStartStrategy=="trim", ...
            'nMHEv2:ReciprocalProviderScope', ...
            ['The exact Case-A reciprocal audit cannot coexist with ', ...
             'scheduled package history, known-chi replacement, or a ', ...
             'linear-batch cold start.']);
        candidatePath = char(string(request.candidatePath));
        expectedHash = lower(string(request.candidateSha256));
        assert(isfile(candidatePath) && ...
            obj.localFileHash(candidatePath)==expectedHash, ...
            'nMHEv2:ReciprocalProviderCandidateHash', ...
            'The reciprocal provider candidate is unavailable or changed.');
        loaded = load(candidatePath,'candidate');
        assert(isfield(loaded,'candidate'), ...
            'nMHEv2:ReciprocalProviderCandidate', ...
            'The reciprocal provider artifact does not contain candidate.');
        obj.reciprocalProvider = ...
            AeroFlex.ctrl.ReciprocalControllerModelProvider( ...
                obj.model,loaded.candidate,obj.uModelTrim);
        if isfield(cfg.ctrl,'compiledReciprocalIntervalProvider') && ...
                ~isempty(cfg.ctrl.compiledReciprocalIntervalProvider)
            obj.reciprocalProvider.configureFixedIntervalKernelAudit( ...
                cfg.ctrl.compiledReciprocalIntervalProvider);
        end
        obj.reciprocalLatentIndex = ...
            obj.nativeStateCount+(1:obj.reciprocalProvider.hiddenStateCount);
        if isfield(request,'formulation') && ~isempty(request.formulation)
            formulation = lower(string(request.formulation));
        else
            formulation = "augmented_reference";
        end
        assert(isscalar(formulation) && ismember(formulation, ...
            ["augmented_reference","condensed_internal_memory"]), ...
            'nMHEv2:ReciprocalProviderFormulation', ...
            'The reciprocal provider formulation is unsupported.');
        obj.reciprocalProviderFormulation = formulation;
        obj.reciprocalCondensedEnabled = ...
            formulation=="condensed_internal_memory";
        if isfield(request,'disturbanceWarmStart') && ...
                ~isempty(request.disturbanceWarmStart)
            disturbanceWarmStart = lower(string( ...
                request.disturbanceWarmStart));
        else
            disturbanceWarmStart = "hold_last";
        end
        assert(isscalar(disturbanceWarmStart) && ismember( ...
            disturbanceWarmStart,["hold_last","bounded_secant"]), ...
            'nMHEv2:ReciprocalDisturbanceWarmStart', ...
            'The reciprocal disturbance warm-start policy is unsupported.');
        obj.reciprocalDisturbanceWarmStart = disturbanceWarmStart;
        if isfield(request,'unboundedAerodynamicLagStates') && ...
                ~isempty(request.unboundedAerodynamicLagStates)
            unboundedLagStates = ...
                request.unboundedAerodynamicLagStates;
            assert(isscalar(unboundedLagStates) && ...
                (islogical(unboundedLagStates) || ...
                (isnumeric(unboundedLagStates) && ...
                isfinite(unboundedLagStates) && ...
                ismember(unboundedLagStates,[0 1]))), ...
                'nMHEv2:ReciprocalAerodynamicLagBounds', ...
                ['The reciprocal aerodynamic-lag bound selector must ', ...
                 'be a logical scalar.']);
            obj.reciprocalUnboundedAerodynamicLagStates = ...
                logical(unboundedLagStates);
        end
        augmentedTrim = obj.reciprocalProvider.initialize(nativeTrimState);
        obj.reciprocalLatentTrim = ...
            augmentedTrim(obj.reciprocalLatentIndex);
        if obj.reciprocalCondensedEnabled
            % The source-owned latent coordinates are deterministic model
            % memory, not independent estimation variables.
            obj.nx = obj.nativeStateCount;
            xTrim = nativeTrimState(:);
            obj.reciprocalLatentHistory = repmat( ...
                obj.reciprocalLatentTrim,1,obj.Ne+1);
        else
            obj.nx = obj.reciprocalProvider.stateCount;
            obj.measurementContract = obj.augmentMeasurementContract();
            xTrim = augmentedTrim;
            obj.reciprocalLatentHistory = zeros(0,obj.Ne+1);
        end
        obj.reciprocalContextHistory = ...
            repmat({request.initialContext},1,obj.Ne);
        obj.reciprocalPendingContext = request.initialContext;
        obj.reciprocalPendingContextValid = true;
    end

    function context = prepareScheduledReciprocalContext(obj,context)
    %PREPARESCHEDULEDRECIPROCALCONTEXT Retain one interval's native model.
        if ~isa(obj.reciprocalProvider, ...
                'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
            return
        end
        required = {'scheduledPacket','scheduledPackage'};
        assert(isstruct(context) && isscalar(context) && ...
            all(isfield(context,required)), ...
            'nMHEv2:ScheduledReciprocalContext', ...
            ['The completed scheduled interval requires its packet and ', ...
             'native package.']);
        package = obj.validateScheduledIntervalPackage( ...
            context.scheduledPackage);
        integrator = obj.cfg.modelHandle( ...
            obj.cfg,package.beam,package.aero,package.base);
        context.scheduledModel = AeroFlex.sched.applyToROMIntegrator( ...
            integrator,package,obj.cfg);
    end

    %----------------------------------------------------------------------
    function z = applyDisturbanceWarmStart(obj,z,index)
    %APPLYDISTURBANCEWARMSTART Seed only the newest causal disturbance.
        if obj.Ne <= 1
            return
        end
        latest = z(index.w{end-1});
        switch obj.reciprocalDisturbanceWarmStart
            case "hold_last"
                terminal = latest;
            case "bounded_secant"
                if obj.Ne > 2
                    predecessor = z(index.w{end-2});
                    terminal = latest+(latest-predecessor);
                else
                    terminal = latest;
                end
                terminal = min(max(terminal,obj.wL),obj.wU);
            otherwise
                error('nMHEv2:ReciprocalDisturbanceWarmStart', ...
                    'Unhandled reciprocal disturbance warm-start policy.');
        end
        z(index.w{end}) = terminal;
    end

    %----------------------------------------------------------------------
    function contract = augmentMeasurementContract(obj)
        contract = obj.nativeMeasurementContract;
        contract.measure = @(state) ...
            obj.nativeMeasurementContract.measure( ...
                state(1:obj.nativeStateCount));
        contract.jacobian = @(state) [ ...
            obj.nativeMeasurementContract.jacobian( ...
                state(1:obj.nativeStateCount)), ...
            zeros(obj.nativeMeasurementContract.ny, ...
                numel(obj.reciprocalLatentIndex))];
        contract.linearizationAtTrim.jacobian = [ ...
            obj.nativeMeasurementContract.linearizationAtTrim.jacobian, ...
            zeros(obj.nativeMeasurementContract.ny, ...
                numel(obj.reciprocalLatentIndex))];
        if isfield(contract,'measureInterval')
            contract.measureInterval = @(previous,current,disturbance,input,Ts) ...
                obj.nativeMeasurementContract.measureInterval( ...
                    previous(1:obj.nativeStateCount), ...
                    current(1:obj.nativeStateCount),disturbance,input,Ts);
            contract.jacobianInterval = ...
                @(previous,current,disturbance,input,Ts) ...
                obj.augmentedIntervalJacobian(previous,current, ...
                    disturbance,input,Ts);
        end
        if isfield(contract,'measureIntervalWithContext')
            contract.measureIntervalWithContext = ...
                @(previous,current,disturbance,input,Ts,context) ...
                obj.nativeMeasurementContract.measureIntervalWithContext( ...
                    previous(1:obj.nativeStateCount), ...
                    current(1:obj.nativeStateCount),disturbance,input,Ts,context);
            contract.jacobianIntervalWithContext = ...
                @(previous,current,disturbance,input,Ts,context) ...
                obj.augmentedIntervalJacobianWithContext(previous,current, ...
                    disturbance,input,Ts,context);
        end
    end

    %----------------------------------------------------------------------
    function z = buildReciprocalFeasibleWarmStart(obj,z,valueHorizonData)
    %BUILDRECIPROCALFEASIBLEWARMSTART Condense recorded interval dynamics.
        assert(obj.reciprocalProviderEnabled && ...
            numel(obj.reciprocalContextHistory)==obj.Ne, ...
            'nMHEv2:ReciprocalFeasibleWarmStart', ...
            'The reciprocal feasible rollout requires a complete history.');
        index = obj.buildIndexMaps();
        z = z(:);
        if obj.reciprocalCondensedEnabled
            assert(isequal(size(obj.reciprocalLatentHistory), ...
                [obj.reciprocalProvider.hiddenStateCount,obj.Ne+1]) && ...
                all(isfinite(obj.reciprocalLatentHistory),'all'), ...
                'nMHEv2:ReciprocalLatentHistory', ...
                'The condensed reciprocal latent history is invalid.');
            latent = obj.reciprocalLatentHistory(:,1);
            state = [z(index.x{1});latent];
        else
            z(index.x{1}(obj.reciprocalLatentIndex)) = ...
                obj.Xhist(obj.reciprocalLatentIndex,1);
            state = z(index.x{1});
        end
        if obj.nativeCausalRolloutActive
            for interval = 1:obj.Ne
                context = obj.reciprocalContextForInterval(interval);
                expectedTotal = obj.uModelTrim+obj.Uhist(:,interval);
                assert(norm(context.wingTotal-expectedTotal,inf) <= ...
                    100*eps(max(1,norm(context.wingTotal,inf))), ...
                    'nMHEv2:ReciprocalFeasibleWarmStartInput', ...
                    ['The completed reciprocal context and retained applied ', ...
                     'input history do not describe the same wing command.']);
            end
            if nargin < 3 || isempty(fieldnames(valueHorizonData))
                data = obj.buildNativeEstimatorValueHorizonData();
            else
                data = valueHorizonData;
            end
            disturbances = reshape(z([index.w{:}]),obj.nw,obj.Ne);
            endpoints = obj.evaluateNativeEstimatorCausalRollout( ...
                state(1:obj.nativeStateCount),disturbances,data);
            for interval = 1:obj.Ne
                z(index.x{interval+1}) = ...
                    endpoints(1:obj.nativeStateCount,interval);
            end
            return
        end
        for interval = 1:obj.Ne
            context = obj.reciprocalContextForInterval(interval);
            expectedTotal = obj.uModelTrim+obj.Uhist(:,interval);
            assert(norm(context.wingTotal-expectedTotal,inf) <= ...
                100*eps(max(1,norm(context.wingTotal,inf))), ...
                'nMHEv2:ReciprocalFeasibleWarmStartInput', ...
                ['The completed reciprocal context and retained applied ', ...
                 'input history do not describe the same wing command.']);
            state = obj.reciprocalProvider.propagateInterval( ...
                state,z(index.w{interval}),context,false);
            if obj.reciprocalCondensedEnabled
                z(index.x{interval+1}) = ...
                    state(1:obj.nativeStateCount);
            else
                z(index.x{interval+1}) = state;
            end
        end
    end

    %----------------------------------------------------------------------
    function latentHistory = rebuildReciprocalLatentHistory( ...
            obj,X,W,valueHorizonData,acceptedReplayEndpoints)
    %REBUILDRECIPROCALLATENTHISTORY Persist exact accepted internal memory.
        assert(obj.reciprocalCondensedEnabled && ...
            isequal(size(X),[obj.nativeStateCount,obj.Ne+1]) && ...
            isequal(size(W),[obj.nw,obj.Ne]), ...
            'nMHEv2:ReciprocalLatentRebuild', ...
            'The accepted condensed reciprocal horizon is incompatible.');
        latentHistory = zeros( ...
            obj.reciprocalProvider.hiddenStateCount,obj.Ne+1);
        latentHistory(:,1) = obj.reciprocalLatentHistory(:,1);
        if obj.latentHistoryCondensationActive
            if nargin >= 5 && ~isempty(acceptedReplayEndpoints)
                endpoints = acceptedReplayEndpoints;
            else
                if nargin < 4 || isempty(fieldnames(valueHorizonData))
                    data = obj.buildNativeEstimatorValueHorizonData();
                else
                    data = valueHorizonData;
                end
                endpoints = obj.evaluateNativeEstimatorValueHorizon(X,W,data);
            end
            assert(isequal(size(endpoints), ...
                [obj.reciprocalProvider.stateCount,obj.Ne]) && ...
                all(isfinite(endpoints),'all'), ...
                'nMHEv2:AcceptedReplayReuseOutput', ...
                'The accepted estimator replay endpoints are invalid.');
            latentHistory(:,2:obj.Ne+1) = ...
                endpoints(obj.reciprocalLatentIndex,:);
            assert(all(isfinite(latentHistory),'all'), ...
                'nMHEv2:ReciprocalLatentRebuild', ...
                'The accepted reciprocal latent history is nonfinite.');
            return
        end
        for interval = 1:obj.Ne
            augmented = obj.reciprocalProvider.propagateInterval( ...
                [X(:,interval);latentHistory(:,interval)], ...
                W(:,interval), ...
                obj.reciprocalContextForInterval(interval),false);
            latentHistory(:,interval+1) = ...
                augmented(obj.reciprocalLatentIndex);
        end
        assert(all(isfinite(latentHistory),'all'), ...
            'nMHEv2:ReciprocalLatentRebuild', ...
            'The accepted reciprocal latent history is nonfinite.');
    end

    %----------------------------------------------------------------------
    function terminalStm = refreshReciprocalTerminalStmExact( ...
            obj,X,W,latentHistory)
    %REFRESHRECIPROCALTERMINALSTMEXACT Reproduce the retained local STM.
        assert(obj.terminalStmCondensationActive && ...
            obj.reciprocalCondensedEnabled && ...
            isequal(size(X),[obj.nativeStateCount,obj.Ne+1]) && ...
            isequal(size(W),[obj.nw,obj.Ne]) && ...
            isequal(size(latentHistory), ...
                [obj.reciprocalProvider.hiddenStateCount,obj.Ne+1]), ...
            'nMHEv2:TerminalStmCondensationState', ...
            'The exact terminal-STM refresh input is incompatible.');
        terminalStart = [X(:,obj.Ne);latentHistory(:,obj.Ne)];
        [~,localSensitivity] = obj.reciprocalProvider.propagateInterval( ...
            terminalStart,W(:,obj.Ne), ...
            obj.reciprocalContextForInterval(obj.Ne),true);
        disturbanceColumns = obj.reciprocalProvider.stateCount+(1:obj.nw);
        terminalStm = [ ...
            localSensitivity(1:obj.nx,1:obj.nx), ...
            localSensitivity(1:obj.nx,disturbanceColumns), ...
            zeros(obj.nx,obj.nu)];
        assert(isequal(size(terminalStm), ...
            [obj.nx,obj.nx+obj.nw+obj.nu]) && ...
            all(isfinite(terminalStm),'all'), ...
            'nMHEv2:TerminalStmCondensationOutput', ...
            'The exact terminal-STM refresh output is invalid.');
    end

    %----------------------------------------------------------------------
    function state = reciprocalValidationState(obj)
    %RECIPROCALVALIDATIONSTATE Return the provider-compatible trim state.
        if obj.reciprocalCondensedEnabled
            state = [obj.stateModelTrim;obj.reciprocalLatentTrim];
        else
            state = obj.stateModelTrim;
        end
    end

    %----------------------------------------------------------------------
    function [previousJacobian,currentJacobian,disturbanceJacobian] = ...
            augmentedIntervalJacobian(obj,previous,current,disturbance,input,Ts)
        [previousNative,currentNative,disturbanceJacobian] = ...
            obj.nativeMeasurementContract.jacobianInterval( ...
                previous(1:obj.nativeStateCount), ...
                current(1:obj.nativeStateCount),disturbance,input,Ts);
        zeroColumns = zeros(size(previousNative,1), ...
            numel(obj.reciprocalLatentIndex));
        previousJacobian = [previousNative,zeroColumns];
        currentJacobian = [currentNative,zeroColumns];
    end

    %----------------------------------------------------------------------
    function [previousJacobian,currentJacobian,disturbanceJacobian] = ...
            augmentedIntervalJacobianWithContext(obj,previous,current, ...
            disturbance,input,Ts,context)
        [previousNative,currentNative,disturbanceJacobian] = ...
            obj.nativeMeasurementContract.jacobianIntervalWithContext( ...
                previous(1:obj.nativeStateCount), ...
                current(1:obj.nativeStateCount),disturbance,input,Ts,context);
        zeroColumns = zeros(size(previousNative,1), ...
            numel(obj.reciprocalLatentIndex));
        previousJacobian = [previousNative,zeroColumns];
        currentJacobian = [currentNative,zeroColumns];
    end

    %----------------------------------------------------------------------
    function xhat = outputEstimate(obj)
        if obj.reciprocalProviderEnabled
            xhat = obj.xhat(1:obj.nativeStateCount);
        else
            xhat = obj.xhat;
        end
    end

    %----------------------------------------------------------------------
    function context = reciprocalContextForInterval(obj,interval)
        if ~obj.reciprocalProviderEnabled
            context = struct();
            return
        end
        assert(isscalar(interval) && interval>=1 && interval<=obj.Ne && ...
            numel(obj.reciprocalContextHistory)==obj.Ne && ...
            isstruct(obj.reciprocalContextHistory{interval}), ...
            'nMHEv2:ReciprocalContextHistory', ...
            'The reciprocal completed-interval history is incomplete.');
        context = obj.reciprocalContextHistory{interval};
    end

    %----------------------------------------------------------------------
    function digest = localFileHash(~,path)
        engine = javaMethod( ...
            'getInstance','java.security.MessageDigest','SHA-256');
        fileId = fopen(path,'rb');
        assert(fileId>=0,'nMHEv2:ReciprocalProviderHashOpen', ...
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
    %----------------------------------------------------------------------
    function [prepared,info] = consumePreparedFeedback(obj,u_prev,zSolve0)
        prepared = struct();
        info = struct('enabled',obj.preparationFeedbackAudit, ...
            'available',false,'used',false,'discarded',false, ...
            'preparationSeconds',0,'reason',"disabled");
        if ~obj.preparationFeedbackAudit
            return
        end
        cache = obj.preparedFeedback;
        obj.preparedFeedback = struct('valid',false);
        if ~isstruct(cache) || ~isfield(cache,'valid') || ~cache.valid
            info.reason = "not_prepared";
            return
        end
        info.available = true;
        if ~isfield(cache,'expectedInput') || ...
                ~isequaln(cache.expectedInput(:),u_prev(:))
            info.discarded = true;
            info.reason = "applied_input_mismatch";
            return
        end
        if ~isfield(cache,'z') || ~isequaln(cache.z(:),zSolve0(:))
            info.discarded = true;
            info.reason = "warm_start_mismatch";
            return
        end
        if ~isfield(cache,'modelL') || ~isequaln(cache.modelL,obj.model.L)
            info.discarded = true;
            info.reason = "prediction_model_mismatch";
            return
        end
        if ~isfield(cache,'modelParConst') || ...
                ~isequaln(cache.modelParConst,obj.model.parConst)
            info.discarded = true;
            info.reason = "prediction_package_mismatch";
            return
        end
        if ~isfield(cache,'prepared') || ~isstruct(cache.prepared)
            info.discarded = true;
            info.reason = "malformed_cache";
            return
        end
        prepared = cache.prepared;
        info.used = true;
        info.reason = "prepared_feedback";
        if isfield(cache,'preparationSeconds')
            info.preparationSeconds = cache.preparationSeconds;
        end
    end
    %----------------------------------------------------------------------
    function contract = buildMeasurementContract(obj,cfg,trimState)
        if isfield(cfg,'mheMeasurement') && ~isempty(cfg.mheMeasurement)
            contract = cfg.mheMeasurement;
            required = {'measure','jacobian','linearizationAtTrim'};
            assert(isstruct(contract) && all(isfield(contract,required)), ...
                'nMHEv2:MeasurementContract', ...
                'mheMeasurement requires measure, jacobian, and linearizationAtTrim.');
            assert(isfield(contract.linearizationAtTrim,'measurement') && ...
                isfield(contract.linearizationAtTrim,'jacobian'), ...
                'nMHEv2:MeasurementContract', ...
                'The measurement trim linearization is incomplete.');
            if isfield(contract,'intervalRows') && ~isempty(contract.intervalRows)
                intervalRequired = {'measureInterval','jacobianInterval'};
                assert(all(isfield(contract,intervalRequired)), ...
                    'nMHEv2:IntervalMeasurementContract', ...
                    ['An interval measurement requires measureInterval and ', ...
                     'jacobianInterval handles.']);
            end
        else
            C = zeros(size(obj.sensor.PhiY,1),obj.nx);
            C(:,obj.model.idx.q1) = obj.sensor.PhiY;
            contract = struct();
            contract.measure = @(state) obj.sensor.measure(state(obj.model.idx.q1));
            contract.jacobian = @(~) C;
            contract.linearizationAtTrim = struct( ...
                'measurement',contract.measure(trimState), ...
                'jacobian',C);
        end
        if isfield(contract,'contextWidth') && ~isempty(contract.contextWidth)
            assert(isscalar(contract.contextWidth) && ...
                isfinite(contract.contextWidth) && ...
                contract.contextWidth == floor(contract.contextWidth) && ...
                contract.contextWidth >= 0, ...
                'nMHEv2:MeasurementContextContract', ...
                'contextWidth must be a finite nonnegative integer scalar.');
        else
            contract.contextWidth = 0;
        end
        if contract.contextWidth > 0
            requiredContextFields = {'measureIntervalWithContext', ...
                'jacobianIntervalWithContext','validateContext'};
            assert(all(isfield(contract,requiredContextFields)) && ...
                isfield(contract.linearizationAtTrim,'context'), ...
                'nMHEv2:MeasurementContextContract', ...
                ['A context-aware measurement requires interval measurement, ', ...
                 'Jacobian, validator, and trim-context fields.']);
            trimContext = contract.linearizationAtTrim.context(:);
            assert(numel(trimContext) == contract.contextWidth && ...
                all(isfinite(trimContext)), ...
                'nMHEv2:MeasurementContextContract', ...
                'The trim measurement context has invalid dimensions or values.');
            contract.validateContext(trimContext);
            contract.linearizationAtTrim.context = trimContext;
        else
            contract.linearizationAtTrim.context = zeros(0,1);
        end
        contract.ny = numel(contract.linearizationAtTrim.measurement);
        assert(contract.ny > 0 && ...
            isequal(size(contract.linearizationAtTrim.jacobian),[contract.ny,obj.nx]), ...
            'nMHEv2:MeasurementContract','Measurement linearization dimensions are invalid.');
        assert(all(isfinite(contract.linearizationAtTrim.measurement)) && ...
            all(isfinite(contract.linearizationAtTrim.jacobian),'all'), ...
            'nMHEv2:MeasurementContract','Measurement trim linearization must be finite.');
        if isfield(contract,'intervalRows') && ~isempty(contract.intervalRows)
            contract.intervalRows = unique(contract.intervalRows(:));
            assert(all(isfinite(contract.intervalRows)) && ...
                all(contract.intervalRows == floor(contract.intervalRows)) && ...
                all(contract.intervalRows >= 1) && ...
                all(contract.intervalRows <= contract.ny), ...
                'nMHEv2:IntervalMeasurementContract', ...
                'intervalRows must contain valid unique measurement rows.');
            nodeRows = setdiff((1:contract.ny).',contract.intervalRows);
            assert(isequal(size(cfg.Qe),[contract.ny,contract.ny]), ...
                'nMHEv2:IntervalMeasurementWeight', ...
                'cfg.Qe dimensions must match the interval measurement contract.');
            crossWeight = cfg.Qe(nodeRows,contract.intervalRows);
            assert(norm(crossWeight,'fro') <= 100*eps* ...
                max(1,norm(cfg.Qe,'fro')), ...
                'nMHEv2:IntervalMeasurementWeight', ...
                'Qe must not couple node and interval measurement rows.');
        else
            contract.intervalRows = zeros(0,1);
        end
    end
    %----------------------------------------------------------------------
    function y = measureState(obj,state)
        y = obj.measurementContract.measure(state(:));
        y = y(:);
        assert(numel(y) == obj.ny && all(isfinite(y)), ...
            'nMHEv2:MeasurementContract','Measurement output is invalid.');
    end
    %----------------------------------------------------------------------
    function jacobian = measurementJacobian(obj,state)
        jacobian = obj.measurementContract.jacobian(state(:));
        assert(isequal(size(jacobian),[obj.ny,obj.nx]) && ...
            all(isfinite(jacobian),'all'), ...
            'nMHEv2:MeasurementContract','Measurement Jacobian is invalid.');
    end
    %----------------------------------------------------------------------
    function active = hasIntervalMeasurement(obj)
        active = ~isempty(obj.measurementContract.intervalRows);
    end
    %----------------------------------------------------------------------
function [measurementContext,knownChi,intervalPackage] = parseEstimateContexts(obj,varargin)
    intervalPackage = [];
    expectedPackageCount = double(obj.scheduledPackageHistoryEnabled);
    if ~obj.knownChiContextEnabled
        if ~obj.scheduledPackageHistoryEnabled
            % Preserve the established context-free/measurement-only API.
            measurementContext = obj.parseMeasurementContext(varargin{:});
            knownChi = zeros(0,1);
            return
        end
        measurementWidth = obj.measurementContract.contextWidth;
        expectedCount = double(measurementWidth > 0)+expectedPackageCount;
        assert(numel(varargin) == expectedCount, ...
            'nMHEv2:EstimateContextCount', ...
            'The active estimator contract requires %d ordered context arguments.', ...
            expectedCount);
        cursor = 0;
        if measurementWidth > 0
            cursor = cursor+1;
            measurementContext = obj.parseMeasurementContext(varargin{cursor});
        else
            measurementContext = zeros(0,1);
        end
        if obj.scheduledPackageHistoryEnabled
            cursor = cursor+1;
            intervalPackage = obj.validateScheduledIntervalPackage(varargin{cursor});
        end
        knownChi = zeros(0,1);
        return
    end
    measurementWidth = obj.measurementContract.contextWidth;
    expectedCount = double(measurementWidth > 0)+ ...
            double(obj.knownChiContextEnabled)+expectedPackageCount;
        assert(numel(varargin) == expectedCount, ...
            'nMHEv2:EstimateContextCount', ...
            ['The active estimator contract requires %d ordered context ', ...
             'arguments (measurement context first, known chi second).'], ...
            expectedCount);
        cursor = 0;
        if measurementWidth > 0
            cursor = cursor+1;
            measurementContext = obj.parseMeasurementContext(varargin{cursor});
        else
            measurementContext = zeros(0,1);
        end
    if obj.knownChiContextEnabled
            cursor = cursor+1;
            knownChi = varargin{cursor};
            assert(isnumeric(knownChi) && isreal(knownChi), ...
                'nMHEv2:KnownChiContext', ...
                'Known chi context must be a real numeric vector.');
            knownChi = knownChi(:);
            assert(numel(knownChi) == 3 && all(isfinite(knownChi)), ...
                'nMHEv2:KnownChiContext', ...
                'Known chi context must contain three finite values.');
            assert(norm(knownChi([1 3]),inf) <= 1e-12, ...
                'nMHEv2:KnownChiLateralScope', ...
                ['The approved known-chi estimator contract is restricted ', ...
                 'to symmetric-longitudinal operation.']);
    else
        knownChi = zeros(0,1);
    end
    if obj.scheduledPackageHistoryEnabled
        cursor = cursor+1;
        intervalPackage = obj.validateScheduledIntervalPackage(varargin{cursor});
    end
end
    %----------------------------------------------------------------------
    function context = parseMeasurementContext(obj,varargin)
        width = obj.measurementContract.contextWidth;
        if width == 0
            assert(isempty(varargin), ...
                'nMHEv2:UnexpectedMeasurementContext', ...
                'This measurement contract does not accept exogenous context.');
            context = zeros(0,1);
            return
        end
        assert(isscalar(varargin), ...
            'nMHEv2:MissingMeasurementContext', ...
            'This measurement contract requires one exogenous context vector.');
        context = varargin{1};
        assert(isnumeric(context) && isreal(context), ...
            'nMHEv2:MeasurementContext', ...
            'Measurement context must be a real numeric vector.');
        context = context(:);
        assert(numel(context) == width && all(isfinite(context)), ...
            'nMHEv2:MeasurementContext', ...
            'Measurement context must contain %d finite values.',width);
        obj.measurementContract.validateContext(context);
    end
    %----------------------------------------------------------------------
    function package = validateScheduledIntervalPackage(obj,package)
    %VALIDATESCHEDULEDINTERVALPACKAGE Check one completed physical interval.
        required = {'L','idx','parConst','beam','aero','base','mu'};
        assert(isstruct(package) && isscalar(package) && ...
            all(isfield(package,required)), ...
            'nMHEv2:ScheduledPackageContract', ...
            'A completed interval must provide the complete source package.');
        assert(isequal(size(package.L),[obj.nx,obj.nx]) && ...
            all(isfinite(package.L),'all') && isequal(package.idx,obj.model.idx), ...
            'nMHEv2:ScheduledPackageStateContract', ...
            'The historical package state order or dynamics are incompatible.');
        assert(isnumeric(package.mu) && isreal(package.mu) && ...
            ~isempty(package.mu) && all(isfinite(package.mu(:))), ...
            'nMHEv2:ScheduledPackageCoordinate', ...
            'The historical package coordinate must be finite.');
        hasChart = isfield(package,'scheduleStateCoordinate') && ...
            isstruct(package.scheduleStateCoordinate) && ...
            isfield(package.scheduleStateCoordinate,'q1ToPhysical') && ...
            isfield(package.scheduleStateCoordinate,'q2ToPhysical') && ...
            isfield(package.scheduleStateCoordinate,'qxiToPhysical') && ...
            isfield(package.scheduleStateCoordinate,'qGammaToFull') && ...
            isfield(package.scheduleStateCoordinate,'qGammaFromFull');
        hasExactFixedBasisMarker = isfield(package,'coordinateProvenance') && ...
            string(package.coordinateProvenance) == "exact_source_fixed_basis";
        assert(hasChart || hasExactFixedBasisMarker, ...
            'nMHEv2:ScheduledPackageChart', ...
            ['The historical package has no qualified state-chart contract ', ...
             'or exact fixed-basis provenance.']);
    end
    %----------------------------------------------------------------------
    function entry = makeScheduledPackageHistoryEntry(obj,package)
    %MAKESCHEDULEDPACKAGEHISTORYENTRY Bind a source package to active chart.
        package = obj.validateScheduledIntervalPackage(package);
        stateToActive = obj.incomingPackageStateToActive;
        assert(isequal(size(stateToActive),[obj.nx,obj.nx]) && ...
            all(isfinite(stateToActive),'all'), ...
            'nMHEv2:IncomingScheduledPackageMap', ...
            'The completed-interval package has no valid active-chart map.');
        entry = struct('package',package,'stateToActive',stateToActive, ...
            'stateFromActive',stateToActive\eye(obj.nx));
    end
    %----------------------------------------------------------------------
    function models = buildIntervalModels(obj)
    %BUILDINTERVALMODELS Return the exact model that owned each interval.
        if ~obj.scheduledPackageHistoryEnabled
            models = repmat({obj.model},1,obj.Ne);
            return
        end
        assert(obj.scheduledPackageHistoryValid && ...
            numel(obj.scheduledPackageHistory) == obj.Ne, ...
            'nMHEv2:ScheduledPackageHistoryIncomplete', ...
            'Scheduled-horizon propagation requires one package per interval.');
        models = cell(1,obj.Ne);
        for interval = 1:obj.Ne
            entry = obj.scheduledPackageHistory{interval};
            assert(isstruct(entry) && isfield(entry,'package') && ...
                isfield(entry,'stateToActive') && ...
                isfield(entry,'stateFromActive'), ...
                'nMHEv2:ScheduledPackageHistoryEntry', ...
                'A retained scheduled package entry is incomplete.');
            package = obj.validateScheduledIntervalPackage(entry.package);
            assert(isequal(size(entry.stateToActive),[obj.nx,obj.nx]) && ...
                isequal(size(entry.stateFromActive),[obj.nx,obj.nx]) && ...
                all(isfinite(entry.stateToActive),'all') && ...
                all(isfinite(entry.stateFromActive),'all') && ...
                norm(entry.stateToActive*entry.stateFromActive-eye(obj.nx),inf) ...
                    <= 1e-9, ...
                'nMHEv2:ScheduledPackageHistoryMap', ...
                'A retained package cannot map between source and active charts.');
            integrator = obj.cfg.modelHandle( ...
                obj.cfg,package.beam,package.aero,package.base);
            integrator = AeroFlex.sched.applyToROMIntegrator( ...
                integrator,package,obj.cfg);
            assert(isequal(integrator.idx,obj.model.idx) && ...
                isequal(size(integrator.L),[obj.nx,obj.nx]) && ...
                all(isfinite(integrator.L),'all'), ...
                'nMHEv2:ScheduledPackageModel', ...
                'A historical package did not construct a compatible model.');
            models{interval} = struct('integrator',integrator, ...
                'stateToActive',entry.stateToActive, ...
                'stateFromActive',entry.stateFromActive);
        end
    end
    %----------------------------------------------------------------------
    function z = bindKnownChiNodes(obj,z)
        if ~obj.knownChiContextEnabled
            return
        end
        index = obj.buildIndexMaps();
        assert(isequal(size(obj.knownChiHistory),[3,obj.Ne+1]) && ...
            all(isfinite(obj.knownChiHistory),'all'), ...
            'nMHEv2:KnownChiHistory', ...
            'The known chi history must be finite and 3-by-(Ne+1).');
        for node = 1:obj.Ne+1
            state = z(index.x{node});
            state(obj.knownChiIndex) = obj.knownChiHistory(:,node);
            z(index.x{node}) = state;
        end
    end
    %----------------------------------------------------------------------
function [state,sensitivity] = propagateKnownChiInterval(obj,state, ...
            input,disturbance,sensitivity,needJacobian,chiStart,chiEnd,model, ...
            reciprocalContext)
        if nargin < 9 || isempty(model)
            model = obj.model;
        end
        if obj.reciprocalProviderEnabled
            assert(nargin>=10 && isstruct(reciprocalContext) && ...
                ~isstruct(model), ...
                'nMHEv2:ReciprocalIntervalContext', ...
                ['The exact Case-A reciprocal provider requires one ', ...
                 'completed-interval context and cannot use a scheduled model.']);
            input = input(:);
            assert(numel(input)==obj.nu && ...
                norm(reciprocalContext.wingTotal-input,inf) <= ...
                    100*eps(max(1,norm(input,inf))), ...
                'nMHEv2:ReciprocalIntervalInputOwner', ...
                ['The retained nMHE input and completed substep commands ', ...
                 'do not describe the same realized interval.']);
            [state,localSensitivity] = ...
                obj.reciprocalProvider.propagateInterval( ...
                    state,disturbance,reciprocalContext,needJacobian);
            if needJacobian
                localState = localSensitivity(:,1:obj.nx);
                localDisturbance = localSensitivity(:,obj.nx+1:end);
                if isempty(sensitivity)
                    sensitivity = [localState,localDisturbance, ...
                        zeros(obj.nx,obj.nu)];
                else
                    sensitivity = localState*sensitivity;
                    sensitivity(:,obj.nx+(1:obj.nw)) = ...
                        sensitivity(:,obj.nx+(1:obj.nw))+ ...
                        localDisturbance;
                end
            else
                sensitivity = [];
            end
            return
        end
        if isstruct(model)
            assert(isfield(model,'integrator') && ...
                isfield(model,'stateToActive') && ...
                isfield(model,'stateFromActive'), ...
                'nMHEv2:ScheduledPackageModel', ...
                'The scheduled interval model is incomplete.');
            integrator = model.integrator;
            stateToActive = model.stateToActive;
            stateFromActive = model.stateFromActive;
        else
            integrator = model;
            stateToActive = eye(obj.nx);
            stateFromActive = eye(obj.nx);
        end
        assert(isprop(integrator,'dt') && isfinite(integrator.dt) && ...
            integrator.dt > 0 && ...
            isequal(size(stateToActive),[obj.nx,obj.nx]) && ...
            isequal(size(stateFromActive),[obj.nx,obj.nx]), ...
            'nMHEv2:ScheduledPackageModel', ...
            'Every interval model must provide a positive finite step.');
        nSubsteps = round(obj.Ts/integrator.dt);
        sampleTolerance = 100*eps(max([1,obj.Ts,integrator.dt]));
        assert(nSubsteps >= 1 && ...
            abs(nSubsteps*integrator.dt-obj.Ts) <= sampleTolerance, ...
            'nMHEv2:ScheduledPackageSampleTime', ...
            'A historical scheduled package does not cover one estimator sample.');
        if ~obj.knownChiContextEnabled
            for substep = 1:nSubsteps
                [state,sensitivity] = obj.stepInActiveChart( ...
                    integrator,stateToActive,stateFromActive,state,input, ...
                    disturbance,sensitivity,needJacobian);
            end
            return
        end
        chiStart = chiStart(:);
        chiEnd = chiEnd(:);
        assert(numel(chiStart) == 3 && numel(chiEnd) == 3 && ...
            all(isfinite(chiStart)) && all(isfinite(chiEnd)), ...
            'nMHEv2:KnownChiInterval', ...
            'Known chi interval endpoints must be finite three-vectors.');
        for substep = 1:nSubsteps
            startFraction = (substep-1)/nSubsteps;
            endFraction = substep/nSubsteps;
            state(obj.knownChiIndex) = ...
                (1-startFraction)*chiStart+startFraction*chiEnd;
            if needJacobian && ~isempty(sensitivity)
                sensitivity(obj.knownChiIndex,:) = 0;
            end
            [state,sensitivity] = obj.stepInActiveChart( ...
                integrator,stateToActive,stateFromActive,state,input, ...
                disturbance,sensitivity,needJacobian);
            state(obj.knownChiIndex) = ...
                (1-endFraction)*chiStart+endFraction*chiEnd;
            if needJacobian && ~isempty(sensitivity)
                sensitivity(obj.knownChiIndex,:) = 0;
            end
        end
    end
    %----------------------------------------------------------------------
    function [state,sensitivity] = stepInActiveChart(~,integrator, ...
            stateToActive,stateFromActive,state,input,disturbance, ...
            sensitivity,needJacobian)
    %STEPINACTIVECHART Propagate a source package in the active estimator chart.
        stateSource = stateFromActive*state;
        if needJacobian && ~isempty(sensitivity)
            sensitivitySource = stateFromActive*sensitivity;
        else
            sensitivitySource = sensitivity;
        end
        [stateSource,sensitivitySource] = integrator.step( ...
            stateSource,input,disturbance,sensitivitySource,needJacobian);
        state = stateToActive*stateSource;
        if needJacobian && ~isempty(sensitivitySource)
            sensitivity = stateToActive*sensitivitySource;
        else
            sensitivity = sensitivitySource;
        end
    end
    %----------------------------------------------------------------------
    function [y,jacobianPrevious,jacobianCurrent,jacobianDisturbance] = ...
            measureInterval(obj,statePrevious,stateCurrent,disturbance,input,context)
        if obj.measurementContract.contextWidth > 0
            context = context(:);
            obj.measurementContract.validateContext(context);
            y = obj.measurementContract.measureIntervalWithContext( ...
                statePrevious(:),stateCurrent(:),disturbance(:), ...
                input(:),obj.Ts,context);
            [jacobianPrevious,jacobianCurrent,jacobianDisturbance] = ...
                obj.measurementContract.jacobianIntervalWithContext( ...
                statePrevious(:),stateCurrent(:),disturbance(:), ...
                input(:),obj.Ts,context);
        else
            y = obj.measurementContract.measureInterval( ...
                statePrevious(:),stateCurrent(:),disturbance(:),input(:),obj.Ts);
            [jacobianPrevious,jacobianCurrent,jacobianDisturbance] = ...
                obj.measurementContract.jacobianInterval( ...
                statePrevious(:),stateCurrent(:),disturbance(:),input(:),obj.Ts);
        end
        y = y(:);
        assert(numel(y) == obj.ny && all(isfinite(y)), ...
            'nMHEv2:IntervalMeasurementContract', ...
            'Interval measurement output is invalid.');
        assert(isequal(size(jacobianPrevious),[obj.ny,obj.nx]) && ...
            isequal(size(jacobianCurrent),[obj.ny,obj.nx]) && ...
            isequal(size(jacobianDisturbance),[obj.ny,obj.nw]) && ...
            all(isfinite(jacobianPrevious),'all') && ...
            all(isfinite(jacobianCurrent),'all') && ...
            all(isfinite(jacobianDisturbance),'all'), ...
            'nMHEv2:IntervalMeasurementContract', ...
            'Interval measurement Jacobian dimensions or values are invalid.');
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function [z,details] = buildLinearBatchColdStart(obj)
    % Build a model-derived cold horizon, then impose nonlinear continuity.
        details = obj.coldStartInfo("linear_batch");
        details.attempted = true;
        z = obj.z0;

        stateCount = obj.nx;
        disturbanceCount = obj.nw;
        inputCount = obj.nu;
        horizon = obj.Ne;
        nUnknown = stateCount+horizon*disturbanceCount;
        try
            sensitivity = [eye(stateCount), ...
                zeros(stateCount,disturbanceCount+inputCount)];
            state = obj.stateModelTrim;
            [state,sensitivity] = obj.propagateKnownChiInterval( ...
                state,obj.uModelTrim,zeros(disturbanceCount,1), ...
                sensitivity,true,obj.knownChiHistory(:,1), ...
                obj.knownChiHistory(:,2),[], ...
                obj.reciprocalContextForInterval(1));
            Ad = sensitivity(:,1:stateCount);
            Bw = sensitivity(:,stateCount+(1:disturbanceCount));
            Bu = sensitivity(:,stateCount+disturbanceCount+ ...
                (1:inputCount));

            C = obj.measurementContract.linearizationAtTrim.jacobian;
            stateMaps = cell(horizon+1,1);
            knownStates = cell(horizon+1,1);
            stateMaps{1} = [eye(stateCount), ...
                zeros(stateCount,horizon*disturbanceCount)];
            knownStates{1} = zeros(stateCount,1);
            for interval = 1:horizon
                stateMaps{interval+1} = Ad*stateMaps{interval};
                columns = stateCount+(interval-1)*disturbanceCount+ ...
                    (1:disturbanceCount);
                stateMaps{interval+1}(:,columns) = ...
                    stateMaps{interval+1}(:,columns)+Bw;
                knownStates{interval+1} = Ad*knownStates{interval}+ ...
                    Bu*obj.Uhist(:,interval);
            end

            if ~obj.hasIntervalMeasurement()
                measurementMap = zeros((horizon+1)*obj.ny,nUnknown);
                measurementDeviation = zeros((horizon+1)*obj.ny,1);
                yTrim = ...
                    obj.measurementContract.linearizationAtTrim.measurement;
                for node = 1:horizon+1
                    rows = (node-1)*obj.ny+(1:obj.ny);
                    measurementMap(rows,:) = C*stateMaps{node};
                    measurementDeviation(rows) = obj.Yhist(:,node)- ...
                        yTrim-C*knownStates{node};
                end
                measurementWeight = kron( ...
                    speye(horizon+1),sparse(obj.Qe));
            else
                intervalRows = obj.measurementContract.intervalRows;
                nodeRows = setdiff((1:obj.ny).',intervalRows);
                nodeRowCount = numel(nodeRows);
                intervalRowCount = numel(intervalRows);
                totalRows = (horizon+1)*nodeRowCount+ ...
                    horizon*intervalRowCount;
                measurementMap = zeros(totalRows,nUnknown);
                measurementDeviation = zeros(totalRows,1);
                yTrim = ...
                    obj.measurementContract.linearizationAtTrim.measurement;
                cursor = 0;
                for node = 1:horizon+1
                    rows = cursor+(1:nodeRowCount);
                    measurementMap(rows,:) = ...
                        C(nodeRows,:)*stateMaps{node};
                    measurementDeviation(rows) = ...
                        obj.Yhist(nodeRows,node)-yTrim(nodeRows)- ...
                        C(nodeRows,:)*knownStates{node};
                    cursor = cursor+nodeRowCount;
                end
                for interval = 1:horizon
                    input = obj.uModelTrim+obj.Uhist(:,interval);
                    [yInterval,Cprevious,Ccurrent,Cdisturbance] = ...
                        obj.measureInterval(obj.stateModelTrim, ...
                            obj.stateModelTrim, ...
                            zeros(disturbanceCount,1),input, ...
                            obj.measurementContextHistory(:,interval+1));
                    rows = cursor+(1:intervalRowCount);
                    intervalMap = ...
                        Cprevious(intervalRows,:)*stateMaps{interval}+ ...
                        Ccurrent(intervalRows,:)*stateMaps{interval+1};
                    columns = stateCount+(interval-1)* ...
                        disturbanceCount+(1:disturbanceCount);
                    intervalMap(:,columns) = ...
                        intervalMap(:,columns)+ ...
                        Cdisturbance(intervalRows,:);
                    measurementMap(rows,:) = intervalMap;
                    measurementDeviation(rows) = ...
                        obj.Yhist(intervalRows,interval+1)- ...
                        yInterval(intervalRows)- ...
                        Cprevious(intervalRows,:)* ...
                            knownStates{interval}- ...
                        Ccurrent(intervalRows,:)* ...
                            knownStates{interval+1};
                    cursor = cursor+intervalRowCount;
                end
                measurementWeight = blkdiag( ...
                    kron(speye(horizon+1), ...
                        sparse(obj.Qe(nodeRows,nodeRows))), ...
                    kron(speye(horizon), ...
                        sparse(obj.Qe(intervalRows,intervalRows))));
            end

            if all(obj.Uhist == 0,"all") && ...
                    all(measurementDeviation == 0)
                % Retain the accepted locked-trim horizon when the complete
                % measurement batch has exactly zero modeled deviation.
                details.exitflag = 1;
                details.measurementResidualTwo = 0;
                details.stateBoundViolationInf = max([0, ...
                    max(obj.xL-obj.stateModelTrim,[],"all"), ...
                    max(obj.stateModelTrim-obj.xU,[],"all")]);
                zeroDisturbance = zeros(disturbanceCount,1);
                details.disturbanceBoundViolationInf = max([0, ...
                    max(obj.wL-zeroDisturbance,[],"all"), ...
                    max(zeroDisturbance-obj.wU,[],"all")]);
                details.disturbanceInfinityNorm = 0;
                details.accepted = ...
                    details.stateBoundViolationInf == 0 && ...
                    details.disturbanceBoundViolationInf == 0;
                return
            end

            regularization = sparse(nUnknown,nUnknown);
            regularization(1:stateCount,1:stateCount) = sparse(obj.Pe);
            firstW = stateCount+(1:disturbanceCount);
            regularization(firstW,firstW) = ...
                regularization(firstW,firstW)+sparse(obj.Re);

            if isfield(obj.cfg,'RdWe') && ~isempty(obj.cfg.RdWe)
                for interval = 2:horizon
                    previousW = stateCount+ ...
                        (interval-2)*disturbanceCount+ ...
                        (1:disturbanceCount);
                    currentW = stateCount+ ...
                        (interval-1)*disturbanceCount+ ...
                        (1:disturbanceCount);
                    rows = (1:disturbanceCount).';
                    differenceMap = sparse([rows;rows], ...
                        [previousW(:);currentW(:)], ...
                        [-ones(disturbanceCount,1); ...
                         ones(disturbanceCount,1)], ...
                        disturbanceCount,nUnknown);
                    regularization = regularization + ...
                        differenceMap.'*sparse(obj.cfg.RdWe)* ...
                        differenceMap;
                end
            end
            if isfield(obj.cfg,'mhe') && ...
                    isfield(obj.cfg.mhe,'RwTerminal') && ...
                    ~isempty(obj.cfg.mhe.RwTerminal)
                terminalW = stateCount+ ...
                    (horizon-1)*disturbanceCount+ ...
                    (1:disturbanceCount);
                regularization(terminalW,terminalW) = ...
                    regularization(terminalW,terminalW)+ ...
                    sparse(obj.cfg.mhe.RwTerminal);
            end

            hessian = measurementMap.'*measurementWeight* ...
                measurementMap+regularization;
            hessian = 0.5*(hessian+hessian.');
            gradient = -measurementMap.'*measurementWeight* ...
                measurementDeviation;
            lowerBound = [obj.xL-obj.stateModelTrim; ...
                repmat(obj.wL,horizon,1)];
            upperBound = [obj.xU-obj.stateModelTrim; ...
                repmat(obj.wU,horizon,1)];

            if horizon > 1
                holdLast = sparse(disturbanceCount,nUnknown);
                previousW = stateCount+ ...
                    (horizon-2)*disturbanceCount+ ...
                    (1:disturbanceCount);
                terminalW = stateCount+ ...
                    (horizon-1)*disturbanceCount+ ...
                    (1:disturbanceCount);
                holdLast(:,previousW) = -speye(disturbanceCount);
                holdLast(:,terminalW) = speye(disturbanceCount);
                holdRightHandSide = zeros(disturbanceCount,1);
            else
                holdLast = [];
                holdRightHandSide = [];
            end

            % The initializer QP is small relative to the nonlinear NLP;
            % tight tolerances prevent its numerical error from consuming
            % the nonlinear feasibility allowance.
            options = optimoptions('quadprog', ...
                'Algorithm','interior-point-convex', ...
                'Display','off', ...
                'OptimalityTolerance',1e-9, ...
                'ConstraintTolerance',1e-9);
            [decision,~,exitflag] = quadprog( ...
                hessian,gradient,[],[],holdLast,holdRightHandSide, ...
                lowerBound,upperBound,[],options);
            details.exitflag = exitflag;
            if isempty(decision) || exitflag <= 0 || ...
                    ~all(isfinite(decision))
                details.message = ...
                    "The initializer QP did not return a finite solution.";
                return
            end

            initialState = obj.stateModelTrim+decision(1:stateCount);
            disturbance = reshape( ...
                decision(stateCount+1:end),disturbanceCount,horizon);
            candidate = zeros( ...
                (horizon+1)*stateCount+horizon*disturbanceCount,1);
            for node = 1:horizon+1
                stateRows = (node-1)*stateCount+(1:stateCount);
                candidate(stateRows) = initialState;
                if node <= horizon
                    state = initialState;
                    input = obj.uModelTrim+obj.Uhist(:,node);
                    [state,~] = obj.propagateKnownChiInterval( ...
                        state,input,disturbance(:,node),[],false, ...
                        obj.knownChiHistory(:,node), ...
                        obj.knownChiHistory(:,node+1),[], ...
                        obj.reciprocalContextForInterval(node));
                    initialState = state;
                end
            end
            candidate((horizon+1)*stateCount+1:end) = disturbance(:);
            stateMatrix = reshape( ...
                candidate(1:(horizon+1)*stateCount), ...
                stateCount,horizon+1);
            details.measurementResidualTwo = norm( ...
                measurementMap*decision-measurementDeviation,2);
            details.stateBoundViolationInf = max([0, ...
                max(obj.xL-stateMatrix,[],"all"), ...
                max(stateMatrix-obj.xU,[],"all")]);
            details.disturbanceBoundViolationInf = max([0, ...
                max(obj.wL-disturbance,[],"all"), ...
                max(disturbance-obj.wU,[],"all")]);
            details.disturbanceInfinityNorm = norm(disturbance,inf);
            details.accepted = all(isfinite(candidate)) && ...
                isfinite(details.measurementResidualTwo) && ...
                details.stateBoundViolationInf == 0 && ...
                details.disturbanceBoundViolationInf == 0;
            if details.accepted
                z = candidate;
            else
                details.message = ...
                    "The nonlinear initializer rollout violated a bound.";
            end
        catch exception
            details.exitflag = -1;
            details.message = string(exception.message);
        end
    end

    %----------------------------------------------------------------------
    function info = coldStartInfo(~,strategy)
    % Return a stable, compact initializer information contract.
        info = struct( ...
            'strategy',string(strategy), ...
            'attempted',false, ...
            'accepted',false, ...
            'fallbackToTrim',false, ...
            'exitflag',NaN, ...
            'seconds',0, ...
            'measurementResidualTwo',NaN, ...
            'stateBoundViolationInf',NaN, ...
            'disturbanceBoundViolationInf',NaN, ...
            'disturbanceInfinityNorm',NaN, ...
            'message',"");
    end

    %----------------------------------------------------------------------
    function idx = buildIndexMaps(obj)
        nx = obj.nx;   nw = obj.nw;   N = obj.Ne;
        idx.x = cell(N+1,1);
        for j = 1:N+1,   idx.x{j} = (j-1)*nx + (1:nx);  end
        startW = (N+1)*nx;
        idx.w  = cell(N,1);
        for j = 1:N,     idx.w{j} = startW + (j-1)*nw + (1:nw);  end
    end
    %----------------------------------------------------------------------
    % function zShift = shiftGuess(obj,p)
    %     idx = obj.buildIndexMaps();  N = obj.Ne;
    %     zShift = zeros(size(p));
    %     for j = 1:N
    %         zShift(idx.x{j}) = p(idx.x{j+1});   % x_{k+1} becomes new x_k
    %         zShift(idx.w{j}) = p(idx.w{j+1});   % w_{k+1} → w_k
    %     end
    %     zShift(idx.x{N+1}) = p(idx.x{N+1});     % keep last state
    % end

    function zShift = shiftGuess(obj,p)
        idx = obj.buildIndexMaps(); nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
        zShift=zeros(size(p));

        % Shift states.
        for j = 1:Ne
            zShift(idx.x{j}) = p(idx.x{j+1});
        end
    
        % Shift disturbances.
        for j = 1:Ne-1
            zShift(idx.w{j}) = p(idx.w{j+1});
        end
        zShift(idx.w{Ne}) = p(idx.w{Ne});   % hold-last
    
        % Predict terminal state using last shifted state, last known control,
        % and held disturbance. 
        % Trying this out to see if it could get a better warm-start
        x = zShift(idx.x{Ne});
        w = zShift(idx.w{Ne});
        u = obj.uModelTrim + obj.Uhist(:,end);
    
        Sdummy = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
    
        if obj.knownChiContextEnabled
            terminalChi = obj.knownChiHistory(:,end);
        else
            terminalChi = zeros(0,1);
        end
        if obj.reciprocalCondensedEnabled
            augmented = obj.reciprocalProvider.propagateInterval( ...
                [x;obj.reciprocalLatentHistory(:,end)],w, ...
                obj.reciprocalContextForInterval(obj.Ne),false);
            x = augmented(1:obj.nativeStateCount);
        else
            [x,Sdummy] = obj.propagateKnownChiInterval( ...
                x,u,w,Sdummy,false,terminalChi,terminalChi,[], ...
                obj.reciprocalContextForInterval(obj.Ne));
        end
    
        zShift(idx.x{Ne+1}) = x;
    end
    %----------------------------------------------------------------------
    function data = buildNativeEstimatorValueHorizonData(obj)
    %BUILDNATIVEESTIMATORVALUEHORIZONDATA Pack invariant E3 interval data.
        contexts = reshape(obj.reciprocalContextHistory,1,obj.Ne);
        obj.reciprocalProvider.configureHorizonSensitivityStencilAudit( ...
            contexts);
        packets = AeroFlex.ctrl. ...
            buildScheduledReciprocalHorizonPacketAudit(contexts);
        wingTotal = zeros(4,14,obj.Ne);
        wingIncrement = zeros(4,14,obj.Ne);
        elevator = zeros(1,14,obj.Ne);
        thrust = zeros(1,14,obj.Ne);
        rigid = zeros(9,14,obj.Ne);
        for interval = 1:obj.Ne
            context = contexts{interval};
            wingTotal(:,:,interval) = context.wingTotal;
            wingIncrement(:,:,interval) = context.wingIncrement;
            elevator(:,:,interval) = context.elevatorIncrement;
            thrust(:,:,interval) = context.thrustIncrement;
            rigid(:,:,interval) = context.rigidState;
        end
        data = struct('packets',packets,'wingTotal',wingTotal, ...
            'wingIncrement',wingIncrement,'elevator',elevator, ...
            'thrust',thrust,'rigid',rigid, ...
            'latentInitial',obj.reciprocalLatentHistory(:,1));
    end

    %----------------------------------------------------------------------
    function endpoints = evaluateNativeEstimatorValueHorizon(obj,X,W,data)
    %EVALUATENATIVEESTIMATORVALUEHORIZON Execute exact E3 value replay.
        assert(obj.nativeValueHorizonActive && ...
            isequal(size(X),[obj.nx,obj.Ne+1]) && ...
            isequal(size(W),[obj.nw,obj.Ne]), ...
            'nMHEv2:NativeValueHorizonInput', ...
            'The estimator E3 value-horizon input is incompatible.');
        endpoints = obj.nativeValueHorizonKernel( ...
            X(:,1:obj.Ne),data.latentInitial,W, ...
            data.wingTotal,data.wingIncrement,data.elevator, ...
            data.thrust,data.rigid,data.packets);
        assert(isequal(size(endpoints),[obj.reciprocalProvider.stateCount, ...
            obj.Ne]) && all(isfinite(endpoints),'all'), ...
            'nMHEv2:NativeValueHorizonOutput', ...
            'The estimator E3 value-horizon output is invalid.');
    end

    %----------------------------------------------------------------------
    function endpoints = evaluateNativeEstimatorCausalRollout(obj,x0,W,data)
    %EVALUATENATIVEESTIMATORCAUSALROLLOUT Execute exact feasible rollout.
        assert(obj.nativeCausalRolloutActive && ...
            numel(x0)==obj.nx && isequal(size(W),[obj.nw,obj.Ne]), ...
            'nMHEv2:NativeCausalRolloutInput', ...
            'The estimator causal-rollout input is incompatible.');
        endpoints = obj.nativeCausalRolloutKernel( ...
            x0,data.latentInitial,W,data.wingTotal,data.wingIncrement, ...
            data.elevator,data.thrust,data.rigid,data.packets);
        assert(isequal(size(endpoints),[obj.reciprocalProvider.stateCount, ...
            obj.Ne]) && all(isfinite(endpoints),'all'), ...
            'nMHEv2:NativeCausalRolloutOutput', ...
            'The estimator causal-rollout output is invalid.');
    end

    %----------------------------------------------------------------------
    function nlp = assembleWindow( ...
            obj,uHistoryOverride,valueHorizonDataOverride)
    % Build cost, constraints, bounds for fmincon
    %----------------------------------------------------------------------
        nx = obj.nx;  nw = obj.nw;  N = obj.Ne;
        idx = obj.buildIndexMaps();
        S_last = obj.Sprev;                          % STM at arrival
        if nargin < 2 || isempty(uHistoryOverride)
            Ulocal = obj.Uhist(:,1:obj.Ne);
        else
            Ulocal = uHistoryOverride;
            assert(isequal(size(Ulocal),[obj.nu,obj.Ne]) && ...
                all(isfinite(Ulocal),'all'), ...
                'nMHEv2:PreparationHistory', ...
                'Prepared input history must be finite and nu-by-Ne.');
        end
        valueHorizonData = struct();
        if obj.nativeValueHorizonActive
            if nargin >= 3 && ~isempty(fieldnames(valueHorizonDataOverride))
                valueHorizonData = valueHorizonDataOverride;
            else
                valueHorizonData = ...
                    obj.buildNativeEstimatorValueHorizonData();
            end
            intervalModels = cell(1,obj.Ne);
        else
            intervalModels = obj.buildIntervalModels();
        end
        valueReplayCacheDecision = zeros(0,1);
        valueReplayCacheInequality = zeros(0,1);
        valueReplayCacheEquality = zeros(0,1);
        valueReplayCacheEndpoints = zeros(0,0);
        valueReplayCacheHits = 0;
        valueReplayCacheMisses = 0;

        % -------------------------- COST ---------------------------------
        function [J,gradJ,hessianJ] = costFun(p)
            J = 0;
            gradJ = zeros(numel(p),1);
            buildGaussNewtonHessian = nargout >= 3;
            if buildGaussNewtonHessian
                hessianJ = sparse(numel(p),numel(p));
            end

            % arrival penalty --------------------------------------------
            arrivalX   = obj.Xhist(:,  end-obj.Ne);   % <— oldest kept sample
            % arrivalW   = obj.Whist(:,  end-obj.Ne+1);   % <— oldest kept sample
            arrivalW   = obj.Whist(:,1);                % oldest retained interval
            % dx   = p(idx.x{1}) - obj.Xhist(:,1);
            dx   = p(idx.x{1}) - arrivalX;
            J    = J + 0.5*dx'*obj.Pe*dx;
            gradJ(idx.x{1}) = obj.Pe*dx;
            if buildGaussNewtonHessian
                hessianJ(idx.x{1},idx.x{1}) = ...
                    hessianJ(idx.x{1},idx.x{1})+sparse(obj.Pe);
            end

            % first-slot disturbance penalty
            % dw   = p(idx.w{1}) - obj.Whist(:,1);
            % dw   = p(idx.w{1}) - arrivalW;
            % J    = J + 0.5*dw'*obj.Re*dw;
            % gradJ(idx.w{1}) = obj.Re*dw;

            % innovation terms ------------------------------------------
            Xmat  = reshape(p([idx.x{:}]),nx,[]);
            if ~obj.hasIntervalMeasurement()
                for j = 1:N+1
                    y = obj.measureState(Xmat(:,j));
                    Cf = obj.measurementJacobian(Xmat(:,j));
                    e = obj.Yhist(:,j) - y;
                    J = J + 0.5*e'*obj.Qe*e;
                    gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cf'*(obj.Qe*e);
                    if buildGaussNewtonHessian
                        hessianJ(idx.x{j},idx.x{j}) = ...
                            hessianJ(idx.x{j},idx.x{j})+ ...
                            sparse(Cf'*obj.Qe*Cf);
                    end
                end
            else
                intervalRows = obj.measurementContract.intervalRows;
                nodeRows = setdiff((1:obj.ny).',intervalRows);
                if ~isempty(nodeRows)
                    nodeWeight = obj.Qe(nodeRows,nodeRows);
                    for j = 1:N+1
                        y = obj.measureState(Xmat(:,j));
                        Cf = obj.measurementJacobian(Xmat(:,j));
                        e = obj.Yhist(nodeRows,j)-y(nodeRows);
                        J = J + 0.5*e'*nodeWeight*e;
                        gradJ(idx.x{j}) = gradJ(idx.x{j}) - ...
                            Cf(nodeRows,:)'*(nodeWeight*e);
                        if buildGaussNewtonHessian
                            hessianJ(idx.x{j},idx.x{j}) = ...
                                hessianJ(idx.x{j},idx.x{j})+ ...
                                sparse(Cf(nodeRows,:)'*nodeWeight*Cf(nodeRows,:));
                        end
                    end
                end
                intervalWeight = obj.Qe(intervalRows,intervalRows);
                Wmat = reshape(p([idx.w{:}]),obj.nw,[]);
                for j = 1:N
                    input = obj.uModelTrim+Ulocal(:,j);
                    [y,Cprevious,Ccurrent,Cdisturbance] = obj.measureInterval( ...
                        Xmat(:,j),Xmat(:,j+1),Wmat(:,j),input, ...
                        obj.measurementContextHistory(:,j+1));
                    e = obj.Yhist(intervalRows,j+1)-y(intervalRows);
                    J = J + 0.5*e'*intervalWeight*e;
                    gradJ(idx.x{j}) = gradJ(idx.x{j}) - ...
                        Cprevious(intervalRows,:)'*(intervalWeight*e);
                    gradJ(idx.x{j+1}) = gradJ(idx.x{j+1}) - ...
                        Ccurrent(intervalRows,:)'*(intervalWeight*e);
                    gradJ(idx.w{j}) = gradJ(idx.w{j}) - ...
                        Cdisturbance(intervalRows,:)'*(intervalWeight*e);
                    if buildGaussNewtonHessian
                        residualJacobian = [ ...
                            Cprevious(intervalRows,:), ...
                            Ccurrent(intervalRows,:), ...
                            Cdisturbance(intervalRows,:)];
                        localIndex = [idx.x{j},idx.x{j+1},idx.w{j}];
                        hessianJ(localIndex,localIndex) = ...
                            hessianJ(localIndex,localIndex)+ ...
                            sparse(residualJacobian'*intervalWeight*residualJacobian);
                    end
                end
            end
            
            if isfield(obj.cfg,'mhe') && isfield(obj.cfg.mhe,'penalizeAllW') && obj.cfg.mhe.penalizeAllW
                for j = 1:N
                    wj = p(idx.w{j});
                    J = J + 0.5*wj.'*obj.Re*wj;
                    gradJ(idx.w{j}) = gradJ(idx.w{j}) + obj.Re*wj;
                    if buildGaussNewtonHessian
                        hessianJ(idx.w{j},idx.w{j}) = ...
                            hessianJ(idx.w{j},idx.w{j})+sparse(obj.Re);
                    end
                end
            else
                dw = p(idx.w{1}) - arrivalW;
                J = J + 0.5*dw.'*obj.Re*dw;
                gradJ(idx.w{1}) = gradJ(idx.w{1}) + obj.Re*dw;
                if buildGaussNewtonHessian
                    hessianJ(idx.w{1},idx.w{1}) = ...
                        hessianJ(idx.w{1},idx.w{1})+sparse(obj.Re);
                end
            end

            % Smoothness Term
            if isfield(obj.cfg,'RdWe')
                RdW = obj.cfg.RdWe;
            else
                RdW = [];
            end
            
            if ~isempty(RdW)
                Wmat = reshape(p([idx.w{:}]),obj.nw,[]);
            
                for j = 2:N
                    dwSmooth = Wmat(:,j) - Wmat(:,j-1);
            
                    J = J + 0.5*dwSmooth.'*RdW*dwSmooth;
            
                    gradJ(idx.w{j})   = gradJ(idx.w{j})   + RdW*dwSmooth;
                    gradJ(idx.w{j-1}) = gradJ(idx.w{j-1}) - RdW*dwSmooth;
                    if buildGaussNewtonHessian
                        localIndex = [idx.w{j-1},idx.w{j}];
                        differenceJacobian = [-eye(obj.nw),eye(obj.nw)];
                        hessianJ(localIndex,localIndex) = ...
                            hessianJ(localIndex,localIndex)+ ...
                            sparse(differenceJacobian'*RdW*differenceJacobian);
                    end
                end
            end

            % --------------------------------------------------------------
            % Terminal disturbance penalty:
            % discourages nonzero estimated gust at the newest/current slot.
            % --------------------------------------------------------------
            if isfield(obj.cfg,'mhe') && isfield(obj.cfg.mhe,'RwTerminal') && ...
                    ~isempty(obj.cfg.mhe.RwTerminal)
            
                RwT = obj.cfg.mhe.RwTerminal;
            
                wN = p(idx.w{N});
            
                J = J + 0.5*wN.'*RwT*wN;
                gradJ(idx.w{N}) = gradJ(idx.w{N}) + RwT*wN;
                if buildGaussNewtonHessian
                    hessianJ(idx.w{N},idx.w{N}) = ...
                        hessianJ(idx.w{N},idx.w{N})+sparse(RwT);
                end
            end
            if buildGaussNewtonHessian
                hessianJ = sparse(0.5*(hessianJ+hessianJ.'));
            end
        end

        % ---------------------- CONTINUITY ------------------------------
        function [c,ceq,gradc,gradceq] = nonlFun(z)
            % if obj.method == "single"
            %     c=[];ceq=[];gradc=[];gradceq=[];return
            % end


            needJacobian = nargout >= 3;
            currentReplayEndpoints = zeros(0,0);
            if obj.acceptedReplayReuseActive && ~needJacobian
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

            knownInitialCount = 3*double(obj.knownChiContextEnabled);
            latentInitialCount = numel(obj.reciprocalLatentIndex)* ...
                double(obj.reciprocalProviderEnabled && ...
                    ~obj.reciprocalCondensedEnabled);
            constraintCount = nx*N+knownInitialCount+latentInitialCount;
            ceq = zeros(constraintCount,1);
            c = [];
            gradc = [];
            gradceq = sparse(constraintCount,numel(z));
            sensitivityTemplate = [];
            if needJacobian
                sensitivityTemplate = ...
                    [eye(nx),zeros(nx,nw+obj.nu)];
                Sx0 = sensitivityTemplate(:,1:nx);
                Sw0 = sensitivityTemplate(:,nx+1:nx+nw);
                if obj.reciprocalCondensedEnabled
                    % Condensing latent memory couples later defects to
                    % earlier native states and gust decisions.
                    nnzPerInterval = nx*numel(z);
                else
                    nnzPerInterval = nx + nnz(Sx0) + nnz(Sw0);
                end
                gradceq = spalloc(constraintCount,numel(z), ...
                    N*nnzPerInterval+knownInitialCount+latentInitialCount);
            end

            X = reshape(z(1:(N+1)*nx),nx,N+1);
            W = reshape(z((N+1)*nx+1:end),nw,N);
            % S = [eye(nx) , zeros(nx,nw+obj.nu)];
            if needJacobian
                S_last = sensitivityTemplate;
            end
            % S = obj.Sprev;           % <-- STM that ended the *previous* RT call  = I at k=0 only
            % ceq(1:nx,:) = obj.xhat - X(:, 1);
            % gradceq(1:nx, 1:nx) = -S(1:nx, 1:nx);
            % ceqArrival       = obj.xlast - X(:,1) ;         % (nx×1)
            % ceqArrival       = X(:,1) - obj.xlast;         % (nx×1)
            % ceq(1:nx) = ceqArrival;
            % gradceq(1:nx, 1:nx) = eye(nx,nx);

            % subplot(211), stairs(obj.Uhist(1,:)), title('Uhist after shift')
            % subplot(212), stairs(W), title('W decision order')
            row0 = 0;
            % row0 = row0 + nx;

            Nend = N-1;
            if obj.reciprocalCondensedEnabled
                if isa(obj.reciprocalProvider, ...
                        'AeroFlex.ctrl.ScheduledReciprocalControllerModelProvider')
                    obj.reciprocalProvider. ...
                        configureHorizonSensitivityStencilAudit( ...
                            obj.reciprocalContextHistory);
                end
                latent = obj.reciprocalLatentHistory(:,1);
                if needJacobian
                    latentDerivative = sparse( ...
                        obj.reciprocalProvider.hiddenStateCount,numel(z));
                end
            end
            % Nend = N-2;
            % Nend = N-3;
            useNativeValueHorizon = obj.nativeValueHorizonActive && ...
                obj.reciprocalCondensedEnabled && ~needJacobian;
            if useNativeValueHorizon
                endpoints = obj.evaluateNativeEstimatorValueHorizon( ...
                    X,W,valueHorizonData);
                currentReplayEndpoints = endpoints;
                ceq(1:nx*N) = reshape( ...
                    X(:,2:N+1)-endpoints(1:nx,:),[],1);
                row0 = nx*N;
            else
            for k = 0:Nend
                % disp(k+1)
                x = X(:,k+1);

                % x_n = x;
                wk = W(:,k+1);
                % uk = obj.Uhist(:,k+1);
                uk = obj.uModelTrim + Ulocal(:,k+1);
                if obj.reciprocalCondensedEnabled
                    context = obj.reciprocalContextForInterval(k+1);
                    assert(norm(context.wingTotal-uk,inf) <= ...
                        100*eps(max(1,norm(context.wingTotal,inf))), ...
                        'nMHEv2:ReciprocalCondensedInputOwner', ...
                        ['The condensed reciprocal context and retained ', ...
                         'input history describe different wing commands.']);
                    [augmentedEnd,localSensitivity] = ...
                        obj.reciprocalProvider.propagateInterval( ...
                            [x;latent],wk,context,needJacobian);
                    x_end = augmentedEnd(1:obj.nativeStateCount);

                    rows = row0+(1:nx);
                    ceq(rows) = X(:,k+2)-x_end;
                    if needJacobian
                        augmentedDerivative = sparse( ...
                            obj.reciprocalProvider.stateCount,numel(z));
                        idx_xk = k*nx+(1:nx);
                        idx_xkp1 = (k+1)*nx+(1:nx);
                        idx_wk = (N+1)*nx+k*nw+(1:nw);
                        augmentedDerivative(1:nx,idx_xk) = speye(nx);
                        augmentedDerivative( ...
                            obj.reciprocalLatentIndex,:) = latentDerivative;
                        nextDerivative = ...
                            sparse(localSensitivity(:,1: ...
                                obj.reciprocalProvider.stateCount))* ...
                                augmentedDerivative;
                        nextDerivative(:,idx_wk) = ...
                            nextDerivative(:,idx_wk)+sparse( ...
                                localSensitivity(:, ...
                                    obj.reciprocalProvider.stateCount+ ...
                                    (1:obj.nw)));
                        gradceq(rows,:) = -nextDerivative(1:nx,:);
                        gradceq(rows,idx_xkp1) = ...
                            gradceq(rows,idx_xkp1)+speye(nx);
                        latentDerivative = nextDerivative( ...
                            obj.reciprocalLatentIndex,:);
                        S_last = [ ...
                            localSensitivity(1:nx,1:nx), ...
                            localSensitivity(1:nx, ...
                                obj.reciprocalProvider.stateCount+ ...
                                (1:obj.nw)), ...
                            zeros(nx,obj.nu)];
                    end
                    latent = augmentedEnd(obj.reciprocalLatentIndex);
                    row0 = row0+nx;
                    continue
                end
                % S = S_last;
                S = sensitivityTemplate;
                [x,S] = obj.propagateKnownChiInterval( ...
                    x,uk,wk,S,needJacobian, ...
                    obj.knownChiHistory(:,k+1), ...
                    obj.knownChiHistory(:,k+2), ...
                    intervalModels{k+1}, ...
                    obj.reciprocalContextForInterval(k+1));
                x_end = x;

                rows = row0 + (1:nx);
                % disp(rows(1));
                % ceq(rows) = X(:,k+1) - x;
                ceq(rows) = X(:,k+2) - x_end;
                % ceq(rows) = x_end - X(:,k+2) ;
                % ceq(rows) = X(:,k+1) - x_end;
                % ceq(rows) = X(:,k+2) - S*[x_n ; wk; uk];
                % ceq(rows) = X(:,k+2) - S_end*[x_end; wk; uk]-x_end;
                % ceq(rows) = X(:,k+2) - S_end*[x_n; wk; uk]-x_end;

                if needJacobian
                    idx_xk = k*nx + (1:nx);
                    idx_xkp1 = (k+1)*nx + (1:nx);
                    idx_wk = (N+1)*nx + k*nw + (1:nw);
                    Sx = S(:,1:nx);
                    Sw = S(:,nx+1:nx+nw);
                    gradceq(rows,idx_xkp1) = speye(nx);
                    gradceq(rows,idx_xk) = -sparse(Sx);
                    gradceq(rows,idx_wk) = -sparse(Sw);
                    S_last = S;
                end
               
                row0 = row0 + nx;
            end
            end
            if obj.knownChiContextEnabled
                rows = nx*N+(1:3);
                ceq(rows) = X(obj.knownChiIndex,1)- ...
                    obj.knownChiHistory(:,1);
                if needJacobian
                    gradceq(rows,obj.knownChiIndex) = speye(3);
                end
            end
            if obj.reciprocalProviderEnabled && ...
                    ~obj.reciprocalCondensedEnabled
                rows = nx*N+knownInitialCount+(1:latentInitialCount);
                ceq(rows) = X(obj.reciprocalLatentIndex,1)- ...
                    obj.Xhist(obj.reciprocalLatentIndex,1);
                if needJacobian
                    gradceq(rows,obj.reciprocalLatentIndex) = ...
                        speye(latentInitialCount);
                end
            end
            % ceq = [ceqArrival ; ceq];
            % gArrive = [eye(nx), zeros(nx, numel(z)-nx)];
            % gradceq = [gArrive ; gradceq];
            if needJacobian
                gradceq = gradceq.';
            elseif obj.acceptedReplayReuseActive
                valueReplayCacheDecision = z;
                valueReplayCacheInequality = c;
                valueReplayCacheEquality = ceq;
                valueReplayCacheEndpoints = currentReplayEndpoints;
            end


            % % --- inequality:  |Δw_k| ≤ dWmax ------------------
            % dW  = diff([obj.Whist(:,1) , W],1,2);           % w₀-w_{-1} … w_N-w_{N-1}
            % c   = abs(dW)' - obj.wU;                         % vectorised (N×1)
            % gradc = zeros(nw*N , length(z));               % build once
            % 
            % for k = 1:N
            %     signk = sign(dW(k));            % ±1
            %     gradc(k, idx.w{k}) =  signk;    % ∂/∂w_k
            % end
        end

        % return STM for next warm-start
        function S_out = getSTM(),  S_out = S_last;  end

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

        % -------- bounds -------------------------------------------------
        xLrep = repmat(obj.xL,N+1,1);
        xUrep = repmat(obj.xU,N+1,1);
        wLrep = repmat(obj.wL,N,1);
        wUrep = repmat(obj.wU,N,1);

        nlp.cost   = @costFun;
        nlp.nonl   = @nonlFun;
        nlp.getAcceptedReplay = @getAcceptedReplay;
        if obj.preparedHorizonDataReuseActive
            nlp.valueHorizonData = valueHorizonData;
        end
        nlp.getSTM = @getSTM;
        nlp.lb     = [xLrep; wLrep];
        nlp.ub     = [xUrep; wUrep];
    end
    %----------------------------------------------------------------------
    function uModelTrim = buildModelTrim(obj,trim)
    % Resolve the total model trim input in radians and radians per second.
        nSurf = obj.cfg.ctrl.n_surf;
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
            error('nMHEv2:TrimControl', ...
                'A total wing trim control is required.');
        end

        assert(all(isfinite(uModelTrim)), 'nMHEv2:TrimControl', ...
            'The total model trim control must be finite.');
        assert(norm(uModelTrim(nSurf+1:end),inf) <= 100*eps, ...
            'nMHEv2:TrimRate', ...
            'The locked trim must have zero surface-rate channels.');
    end

    %----------------------------------------------------------------------
    function v = expandToLength(~,v,n,name)
    % Convert scalar/vector input to a validated n-by-1 column.
        v = v(:);
        if isempty(v)
            error('nMHEv2:Config','%s is empty.',name);
        end
        if isscalar(v)
            v = repmat(v,n,1);
            return
        end
        if numel(v) ~= n
            error('nMHEv2:Dimension', ...
                '%s must be scalar or length %d. Got length %d.', ...
                name,n,numel(v));
        end
    end

    %----------------------------------------------------------------------
    function tolerance = constraintTolerance(obj)
    % Return the active solver's declared feasibility tolerance.
        switch obj.solverName
            case "fmincon"
                tolerance = obj.solverOpts.ConstraintTolerance;
            case "custom_sqp"
                tolerance = obj.sqpSolver.options.ConstraintTolerance;
            otherwise
                error('nMHEv2:Solver', ...
                    'Unhandled solverName = "%s".',obj.solverName);
        end
        assert(isscalar(tolerance) && isfinite(tolerance) && tolerance > 0, ...
            'nMHEv2:ConstraintTolerance', ...
            'The solver constraint tolerance must be a positive finite scalar.');
    end

    %----------------------------------------------------------------------
    function violation = constraintViolation(~,inequality,equality, ...
            decision,lowerBound,upperBound)
    % Return the infinity-norm violation used by candidate acceptance.
        if isempty(inequality)
            inequalityViolation = 0;
        else
            inequalityViolation = max(max(inequality),0);
        end
        if isempty(equality)
            equalityViolation = 0;
        else
            equalityViolation = norm(equality,inf);
        end
        violation = max(equalityViolation,inequalityViolation);
        if nargin >= 6
            decision = decision(:);
            lowerBound = lowerBound(:);
            upperBound = upperBound(:);
            assert(numel(decision)==numel(lowerBound) && ...
                numel(decision)==numel(upperBound), ...
                'nMHEv2:ConstraintViolationBounds', ...
                'The decision and bound dimensions are incompatible.');
            boundViolation = max([0;lowerBound-decision; ...
                decision-upperBound]);
            violation = max(violation,boundViolation);
        end
    end

    %----------------------------------------------------------------------
    function message = outputMessage(~,output)
    % Extract a solver message without assuming a solver-specific structure.
        message = '';
        if isstruct(output) && isfield(output,'message') && ...
                ~isempty(output.message)
            message = char(string(output.message));
        end
    end

    %----------------------------------------------------------------------
    function debugPlots(obj,t_k,info)
        wH = info.wHorizon;
         % ---- initialise persistent store on first call --------------------
        if ~isfield(obj.dbg,'t')
            obj.dbg.t = [];          % Ns×1   sample times
            obj.dbg.W = [];          % Ns×Ne  horizon matrix (rows = samples)
        end
        if isempty(obj.dbg)
            obj.dbg.t=t_k; 
            obj.dbg.W=wH(:)';
            obj.dbg.cont=norm(info.wHorizon);
        else
            obj.dbg.t(end+1)=t_k; obj.dbg.W(end+1,1:obj.Ne)=wH(:)';
            obj.dbg.cont(end+1)=info.continuity;
    
        end
        figure(2001),clf,hold on
        for s=1:obj.Ne, stairs(obj.dbg.t,obj.dbg.W(:,s)); end
        hold off,grid on,xlabel('t [s]'),ylabel('w'),title('horizon slots'), legend; %drawnow
        
        figure(2004),clf, hold on,plot(obj.dbg.t,obj.dbg.cont(2:end)); hold off;grid on
        xlabel('t [s]'); ylabel('|continuity|'); title('MHE feasibility'), legend
    end

    function debugPlots2(obj,t_k,info)
        %DEBUGPLOTS Time-aligned MHE disturbance horizon plot.
        
            wH = info.wHorizon;   % expected Ne x nw or Ne x 1 for nw=1
        
            if size(wH,2) > size(wH,1) && obj.nw == 1
                wH = wH(:);
            end
        
            if ~isfield(obj.dbg,'tSolve')
                obj.dbg.tSolve = [];
                obj.dbg.tSlot  = [];
                obj.dbg.Wslot  = [];
                obj.dbg.cont   = [];
            end
        
            obj.dbg.tSolve(end+1,1) = t_k;
            obj.dbg.cont(end+1,1)   = info.continuity;
        
            % Slot physical times: oldest to newest.
            tSlots = zeros(obj.Ne,1);
            for s = 1:obj.Ne
                tSlots(s) = t_k - (obj.Ne - s)*obj.Ts;
            end
        
            obj.dbg.tSlot = [obj.dbg.tSlot; tSlots(:)];
            obj.dbg.Wslot = [obj.dbg.Wslot; wH(:)];
        
            figure(3001); clf;
            plot(obj.dbg.tSlot,obj.dbg.Wslot,'.-','LineWidth',1.2);
            grid on;
            xlabel('physical time associated with MHE slot [s]');
            ylabel('\hat{w}');
            title('Time-aligned MHE disturbance estimates');
        
            figure(3004); clf;
            semilogy(obj.dbg.tSolve,max(obj.dbg.cont,eps),'LineWidth',1.2);
            grid on;
            xlabel('solve time [s]');
            ylabel('||c_{eq}||');
            title('MHE feasibility');
    end

    function localCheckMHEEqualityGradientBlocks(obj,nlp,z,lb,ub)
    %LOCALCHECKMHEEQUALITYGRADIENTBLOCKS Check MHE dynamic equality Jacobian
    % by perturbing only X variables and only W variables separately.
    
        idx = obj.buildIndexMaps();
    
        nx = obj.nx;
        nw = obj.nw;
        Ne = obj.Ne;
    
        z  = z(:);
        lb = lb(:);
        ub = ub(:);
        n  = numel(z);
    
        z = min(max(z,lb),ub);
    
        [~,ceq0,~,Gceq] = nlp.nonl(z);
    
        if size(Gceq,1) ~= n
            Gceq = Gceq.';
        end
    
        % If ceq includes X0-arrival equality, dynamic rows start after nx.
        if numel(ceq0) >= nx*(Ne+1)
            rowsDyn = nx+1:numel(ceq0);
        else
            rowsDyn = 1:numel(ceq0);
        end
    
        h = 1e-6;
    
        % X-only direction.
        dX = zeros(n,1);
        rng(301);
        for k = 1:Ne+1
            dX(idx.x{k}) = randn(nx,1);
        end
        dX = dX / max(norm(dX),eps);
    
        % W-only direction.
        dW = zeros(n,1);
        rng(302);
        for k = 1:Ne
            dW(idx.w{k}) = randn(nw,1);
        end
        dW = dW / max(norm(dW),eps);
    
        errX = obj.localDirectionalEqualityError(nlp,z,dX,h,lb,ub,Gceq,rowsDyn);
        errW = obj.localDirectionalEqualityError(nlp,z,dW,h,lb,ub,Gceq,rowsDyn);
    
        fprintf('\n');
        fprintf('======================================================================\n');
        fprintf(' NMHE DYNAMIC-EQUALITY BLOCK GRADIENT CHECK\n');
        fprintf('======================================================================\n');
        fprintf('  dynamic defect rel error, X-only direction : %.3e\n', errX);
        fprintf('  dynamic defect rel error, W-only direction : %.3e\n', errW);
        fprintf('======================================================================\n\n');
    end

    function err = localDirectionalEqualityError(obj,nlp,z,d,h,lb,ub,Gceq,rowsDyn)
    %LOCALDIRECTIONALEQUALITYERROR Directional FD check for selected equality rows.
    
    
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

    function [stop,options,optchanged] = fminconOutputFunction(obj,~,values,state)
    %FMINCONOUTPUTFUNCTION Capture audit-only fmincon progress diagnostics.
        stop = false;
        options = [];
        optchanged = false;
        phase = string(state);
        if phase == "init"
            obj.fminconIterationTrace = struct('iteration',{},'objective',{}, ...
                'constraintViolation',{},'firstOrderOptimality',{}, ...
                'stepSize',{},'cgIterations',{});
            obj.fminconIterationPlot = obj.createFminconIterationPlot();
        end
        if ~(phase == "init" || phase == "iter" || phase == "done")
            return
        end
        entry = struct('iteration',obj.fminconDiagnosticField(values, ...
            'iteration',0),'objective',obj.fminconDiagnosticField(values, ...
            'fval',nan),'constraintViolation',obj.fminconDiagnosticField( ...
            values,'constrviolation',nan),'firstOrderOptimality', ...
            obj.fminconDiagnosticField(values,'firstorderopt',nan), ...
            'stepSize',obj.fminconDiagnosticField(values,'stepsize',nan), ...
            'cgIterations',obj.fminconDiagnosticField(values,'cgiterations',nan));
        obj.fminconIterationTrace(end+1) = entry;
        obj.updateFminconIterationPlot(entry);
    end

    function value = fminconDiagnosticField(~,values,name,defaultValue)
        value = defaultValue;
        if isstruct(values) && isfield(values,name) && ...
                isscalar(values.(name)) && isfinite(values.(name))
            value = values.(name);
        end
    end

    function plotData = createFminconIterationPlot(~)
        figureHandle = figure('Name','fmincon iteration diagnostics', ...
            'NumberTitle','off','Color','w');
        layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        titles = ["Objective" "Constraint violation" ...
            "First-order optimality" "Step size"];
        lines = gobjects(4,1);
        for index = 1:4
            nexttile(layout,index); hold on; grid on;
            title(titles(index)); xlabel('fmincon iteration');
            set(gca,'YScale','log'); ylim([1e-10 inf]);
            lines(index) = animatedline('LineWidth',1.2,'Marker','o', ...
                'MarkerSize',4);
        end
        title(layout,'Live fmincon convergence (audit diagnostic only)');
        plotData = struct('enabled',true,'figure',figureHandle,'lines',lines);
    end

    function updateFminconIterationPlot(obj,entry)
        plotData = obj.fminconIterationPlot;
        if ~isstruct(plotData) || ~isfield(plotData,'enabled') || ...
                ~plotData.enabled || ~isgraphics(plotData.figure)
            return
        end
        plotFloor = 1e-10;
        values = [entry.objective entry.constraintViolation ...
            entry.firstOrderOptimality entry.stepSize];
        for index = 1:numel(values)
            addpoints(plotData.lines(index),entry.iteration, ...
                max(abs(values(index)),plotFloor));
        end
        drawnow limitrate
    end
    %----------------------------------------------------------------------
    end   % private methods
end   % classdef

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
    'nMHEv2:ScheduledAggregateSelector', ...
    'The scheduled aggregate selector must be logical.');
enabled = logical(selector);
end
