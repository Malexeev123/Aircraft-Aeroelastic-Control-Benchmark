function cfg = nominalConfig
    %NOMINALCONFIG  Return a structure with default simulation parameters.
    %   Call this at the start of any example or test:
    %       cfg = nominalConfig;
    %   The fields are grouped logically (structural sizes, simulation
    %   settings, flight conditions, paths) so downstream classes
    %   (BeamModel, AeroROM, SimRunner, Controllers) can rely on a
    %   single source‑of‑truth.
    cfg.NetworkPath = 0;
    
    % ---------- STRUCTURAL & AERODYNAMIC SIZES ------------------------------
    % These are defaults, they are overwritten by the
    % simulation_parameters.mat
    cfg.Nm = 20;   % Number of structural modes  (also Γ length)
    cfg.Na = 80;   % Combined aerodynamic Γ + χ Roger states
    % ---------- SIMULATION CLOCK -------------------------------------------

    %%% THESE ARE DEFAULTS TO INITIALIZE NOW. DO NOT REMOVE FOR ROBUSTNESS.
        cfg.sim.dt    = 6.25e-4;  % [s] fixed IMEX step (0.000625)
        % cfg.sim.t_end = 1.0;      % [s] total horizon
        cfg.sim.t_end = 5.0;      % [s] total horizon
        cfg.plotOpts.debugPlots = true;
        cfg.plotOpts.wingTip = true;
        cfg.case = 'openLoop'; 
    %%%

    % ------ Custom Sim Overwrite Logic. Comment out as needed and set.
    % To overwrite SHARPy input deck. Good for coupledfull sims
    % future versions will have more robust sim inputs.
    cfg.overwrite.doOverwrite = false; % false | true
    % cfg.overwrite.t_end = 5.0 ; % [s] to overwrite t_end time
    cfg.overwrite.t_end = 0.165  ; % [s] to overwrite t_end time
    cfg.overwrite.gustOldLength = false;

    % ---------- ACTUATOR LOGIC --------------------------------------------
    % To turn on actuator dynamics in Sim
    cfg.actuator.enable = true; % true | false

    % ---------- FLIGHT / TRIM CONDITIONS -----------------------------------
    cfg.struct.retain_rigid = false;
    cfg.struct.rigid_body_modes = false;
    % cfg.struct.damping = false; % This is for Structural damping
    cfg.struct.damping = true;
    cfg.struct.discrete = false;   % need to explore implications of this
    cfg.struct.L = 1.1; % [m] % Wing Span
    
    % ---------- FLIGHT / TRIM CONDITIONS -----------------------------------
    cfg.flight.U_inf = 40;    % [m/s] free‑stream velocity
    cfg.flight.rho   = 1.225; % [kg/m³] air density (ISA sea‑level)
    cfg.flight.b_ref = 0.05;  % [m] reference semi‑chord (0.1 / 2) like for artola
    % cfg.flight.b_ref = 0.1;  % [m] reference chord 
    cfg.flight.aoa_deg = 0; % in degrees
    cfg.flight.S_ref = cfg.struct.L *cfg.flight.b_ref; % [m^2] Since it is rectangular b*c

    q_inf = 0.5*cfg.flight.rho*cfg.flight.U_inf^2;
    cfg.flight.a = q_inf*cfg.flight.b_ref^2;
    cfg.flight.t_inf = cfg.flight.b_ref/cfg.flight.U_inf;

    cfg.force_map.scale =  cfg.flight.a;

    cfg.flight.a = 1;
    cfg.flight.t_inf = 1;


    cfg.trim.do      = true; % run pre‑trim solver? (true/false)
    cfg.trim.tolSteady  = 1e-5;      % e.g. 1e‑7
    % cfg.trim.tolSteady  = 1e-6;      % e.g. 1e‑7
    % cfg.trim.tolNewton  = 1e-4;      % e.g. 1e‑6
    cfg.trim.tolNewton  = 5e-4;      % e.g. 1e‑6
    cfg.trim.alphaDeg  = 0;
    cfg.trim.deltaDeg = 0;
    cfg.trim.thrust   = zeros(cfg.Nm,1);         % modal thrust N
    cfg.trim.itersInt   = 10000;         % N
    cfg.trim.dtQS = 1e-4;
    cfg.aero.discrete = false;   % false | true → use obj.ROM_dsc
    % cfg.library.enable = true;

    cfg.library.noExtrapolate = true;
    cfg.library.updateMode = 'perPlantStep'; % perPlantStep | frozenTrim

    % cfg.trim.alphaBoundsDeg = [-2 10];
    cfg.trim.alphaBoundsDeg = [-4 10];
    cfg.trim.wingDeltaBoundsDeg = [-20 20];
    cfg.trim.elevatorBoundsDeg = [-25 25];
    
    % Use a physically positive thrust range unless reverse thrust is modeled.
    rParams = RigidBody.methods.paramsRigid_PazyUAV();
    cfg.tail = rParams.tail.addition;
    cfg.trim.thrustBoundsN = [0, 2*rParams.mass*9.807];
    
    cfg.trim.alphaGuessDeg = cfg.flight.aoa_deg;
    % cfg.trim.tolSteady = 2.5e-6;
    % cfg.trim.tolSteady = 2.5e-6;
    % cfg.trim.tolSteady = 5e-6;
    cfg.trim.tolSteady = 2.5e-5;
    % cfg.trim.tolSteady = 5e-5;
    % cfg.trim.tolNewtonScaled = 1e-5;
    cfg.trim.tolNewtonAccept = 1e-4;
    cfg.trim.maxNewtonIter = 60;
    cfg.trim.itersInt = 20000;
    cfg.trim.lambda0 = 1e-2;
    
    cfg.trim.thrustActsOnWing = false;
    % The calibrated dual-force solve already returns the full-span wrench.
    cfg.trim.mirrorWingClamp = false;
    cfg.trim.shiftWingMomentToCG = false;
    cfg.debug.trimDerivatives = true;
    % ----------- GUST SETTINGS ---------------------------------------------
    cfg.gust.on      = true;  % include Dryden gust input
    cfg.gust.gust_length= 10.0;
    cfg.gust.gust_intensity = 0.02;
    cfg.gust.gust_offset= 0;
     cfg.gust.gust_offset_ratio = .15; %start at 15% of sim
    cfg.gust.scaleMag = 0.0005; % scaling the gust by. Should be 1.
    % cfg.gust.scaleMag = 0.005; % scaling the gust by. Should be 1.
    % cfg.gust.scaleMag = 1; % scaling the gust by. Should be 1.
    % cfg.gust.scaleMag = .01; % scaling the gust by. Should be 1.

    cfg.gust.a = q_inf*cfg.flight.b_ref^2;
    cfg.gust.t_inf = cfg.flight.b_ref/cfg.flight.U_inf;
    cfg.rigid.a = cfg.gust.a;
    cfg.rigid.t_inf = cfg.gust.t_inf;

    % cfg.gust.a = 1;
    % cfg.gust.t_inf = 1;

    % ---------- CONTROL-SURFACE SETTINGS -----------------------------------
    cfg.ctrl.n_surf  = 2;     % e.g., two aileron flaps
    cfg.ctrl.var_per = 2;     % δ and δ̇ per surface
    cfg.ctrl.Nc       = 6;                     
    cfg.ctrl.Ne       = 6; 
    % cfg.control     = struct();        % unused
    % ---------- DEBUG / VERBOSITY ------------------------------------------
    cfg.debug.level  =1;     % 0 = silent, 3 = verbose + plots
    % cfg.debug.plt_scale = .55/2;
    cfg.debug.plt_scale = 1.1/2;
    cfg.debug.trim = 1;
   

    
    % Optional symmetry handling if using one semi-wing to represent both wings.
    cfg.sim.symmetricWingPair = false;
    cfg.sim.wingMultiplicity  = 1;
    
    % Conservative default: do not overwrite chi unless you explicitly want
    % rigid-body alpha injected into the ROM attitude states.
    cfg.sim.coupleRigidAttitudeIntoWingChi = false;
    
    cfg.library.holdLastOnExtrapolate = false;

    cfg.sim.commandsAreTrimIncrements = true;
    
    % Match trim and runtime convention.
    cfg.sim.mirrorWingClamp = cfg.trim.mirrorWingClamp;
    cfg.sim.shiftWingMomentToCG = cfg.trim.shiftWingMomentToCG;

    cfg.debug.schedule = true;
    cfg.library.debug = true;
    cfg.library.beamMapTol = 1e-8;
    cfg.library.allowIncompatibleBeamMaps = false;
    cfg.library.interpMode = "linear"; % linear | nearest
    % cfg.library.allowIncompatibleBeamMaps = true;   % debug only

    
   % ---------- COUPLED TRIM SOLVER -----------------------------------------
    % Nested trim:
    %   inner  : projected/clamped flexible steady state, ||Pz*q1dot||
    %   outer  : rigid residual [Fx; Fz; My]
    cfg.trim.thrustActsOnWing  = false;
    
    cfg.trim.tolSteady = 1e-6;
    cfg.trim.innerAcceptTol = 5e-6;
    cfg.trim.fdInnerAcceptTol = 2e-5;
    
    cfg.trim.tolPhysical = [ ...
        1.0e-4;    % Fx [N]
        1.0e-4;    % Fz [N]
        1.0e-4];   % My [N*m]
    
    cfg.trim.maxNewtonIter = 80;
    cfg.trim.itersInt = 20000;
    
    % Trust-region LM controls.
    cfg.trim.lambda0 = 1e-2;
    cfg.trim.lambdaMin = 1e-4;
    cfg.trim.lambdaMax = 1e6;
    cfg.trim.lambdaDown = 2;
    cfg.trim.lambdaUp = 5;
    cfg.trim.lambdaUpSoft = 2;
    cfg.trim.lineSearchMax = 10;
    cfg.trim.maxStall = 10;
    
    % Smooth transition from force-dominated to full moment-sensitive trim.
    cfg.trim.forceBlendHi = 0.30;
    cfg.trim.forceBlendLo = 0.03;
    cfg.trim.stepBlendHi = 0.30;
    cfg.trim.stepBlendLo = 0.03;
    cfg.trim.momentArmScale = 0.10;
    
    % Variable metric.  Alpha is intentionally more expensive than control
    % surfaces, because the non-unique trim should not use AoA as the cheap knob.
    cfg.trim.varScale = [ ...
        deg2rad(0.25);
        deg2rad(2.00);
        deg2rad(3.00);
        0.25];
    
    cfg.trim.regWeight = [ ...
        1e-3;
        1e-6;
        1e-7;
        0];
    
    % Residual-dependent trust region.  Coarse values are used far from trim;
    % fine values are blended in near convergence.
    cfg.trim.maxStepCoarse = [ ...
        deg2rad(0.08);
        deg2rad(0.35);
        deg2rad(0.75);
        0.15];
    
    cfg.trim.maxStepFine = [ ...
        deg2rad(0.015);
        deg2rad(0.08);
        deg2rad(0.20);
        0.04];
    
    % Finite-difference steps.  These are intentionally larger than the final
    % tolerance because each residual call includes an inner projected solve.
    cfg.trim.fdStep = [ ...
        deg2rad(0.04);
        deg2rad(0.04);
        deg2rad(0.03);
        0.02];
    
    cfg.trim.fdMaxHalves = 3;
    
    cfg.trim.thrustPolish = true;
    cfg.trim.thrustPolishTol = [1e-4; 1e-4];
    cfg.trim.thrustPolishMyTol = 2e-3;
    
    cfg.trim.endgamePolish = true;
    cfg.trim.endgamePolishNorm = 1.0e-2;
    cfg.trim.endgamePolishVars = [2 3];     % wing, elevator
    cfg.trim.endgamePolishMaxStep = [ ...
        0;
        deg2rad(0.03);
        deg2rad(0.08);
        0];
    cfg.trim.endgamePolishLineMax = 5;
    cfg.trim.endgameAcceptFrac = 0.75;
    cfg.trim.endgameFxRelax = 5e-2;
    
    cfg.library.debug = false;
    cfg.debug.trimSchedule = false;


    % ---------- RUNTIME SCHEDULING SAFETY -----------------------------------
    % perPlantStep is a hybrid model switch.  The candidate ROM is accepted only
    % if it is locally admissible at the current x,u,w.  For no-gust open-loop,
    % these settings should make perPlantStep indistinguishable from frozenTrim
    % unless the aircraft genuinely moves through the scheduling grid.
    cfg.library.updateDeadbandU        = 0.05;   % [m/s]
    % cfg.library.updateDeadbandAlphaDeg = 0.005;   % [deg]
    cfg.library.updateDeadbandAlphaDeg = 0.02;   % [deg]
    cfg.library.minUpdateIntervalSteps = 25;     % avoids repeated LU rebuilds
    % cfg.library.minUpdateIntervalSteps = 5;     % avoids repeated LU rebuilds
    

    
    % Load-jump gate for the current longitudinal coupled model:
    %   [dFx; dFz; dMy]
    cfg.library.loadJumpTol = [ ...
        1.0e-3;    % [N]
        1.0e-3;    % [N]
        5.0e-4];   % [N*m]
    
    % Diagnostics.  Rejections are printed initially, then rate-limited.
    cfg.library.affineStateDiagnostic = false;
    cfg.library.scheduleAcceptedPrintEvery = 25;
    cfg.library.scheduleRejectPrintMax = 12;

    cfg.library.printDeadbandHolds = false;  % set true only for diagnostics
    cfg.library.deadbandPrintEvery = 100;

    cfg.library.inputOffsetDiagnostic = false;
    

    % Runtime scheduler:
    %   reference : use commanded/trajectory alpha when provided; otherwise trim
    %   filtered  : use slowly filtered body alpha, mainly for diagnostic manoeuvres
    %   trimHold  : gust-only validation around one trim point
    %   instant   : debug only
    cfg.library.scheduleAlphaMode = 'reference'; 
    cfg.library.scheduleTauU = 0.0;
    cfg.library.scheduleTauAlpha = 0.10;
    
    % forces_0 is an affine preload/trim offset.  Hold it during runtime
    % schedule switches unless a separate bumpless equilibrium update is added.
    cfg.library.forceOffsetMode = 'holdLast';     % holdLast | blend | scheduled
    cfg.library.forceOffsetBlend = 0.02;
    
    % Candidate schedule admissibility.
    % cfg.library.admissiblePzAbsTol = 5e-5;
    % cfg.library.admissiblePzGrowth = 5.0;
    cfg.library.admissiblePzAbsTol = 5e-3;
    cfg.library.admissiblePzGrowth = 10.0;
    
    % Absolute floor plus relative dynamic allowance.
    % cfg.library.loadJumpTol = [1e-3; 1e-3; 5e-4];
    cfg.library.loadJumpTol = [1e-3; 1e-3; 1e-3];
    % cfg.library.loadJumpRelTol = [0.02; 0.02; 0.02];
    cfg.library.loadJumpRelTol = [0.02; 0.02; 0.04];
    cfg.library.loadJumpRefFloor = [1.0; 1.0; 0.10];
    
    % Keep diagnostics off for normal runs.
    cfg.library.forceOffsetRawDiagnostic = false;
    cfg.library.blockJumpDiagnostic = false;
    cfg.library.blockJumpPrintMax = 6;

    cfg.library.rejectBackoffSteps = max(1, round(0.01/cfg.sim.dt));
    cfg.library.rejectRetryFactor = 1.0;

end
