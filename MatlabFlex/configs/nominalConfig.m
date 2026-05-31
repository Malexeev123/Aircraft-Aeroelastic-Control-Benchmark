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
    cfg.sim.dt    = 6.25e-4;  % [s] fixed IMEX step (0.000625)
    cfg.sim.t_end = 1.0;      % [s] total horizon
    cfg.plotOpts.debugPlots = true;
    cfg.plotOpts.wingTip = true;
    cfg.case = 'openLoop'; 


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
    cfg.trim.tolSteady  = 1e-6;      % e.g. 1e‑7
    % cfg.trim.tolNewton  = 1e-4;      % e.g. 1e‑6
    cfg.trim.tolNewton  = 2.5e-4;      % e.g. 1e‑6
    cfg.trim.alphaDeg  = 0;
    cfg.trim.deltaDeg = 0;
    cfg.trim.thrust   = zeros(cfg.Nm,1);         % modal thrust N
    cfg.trim.itersInt   = 10000;         % N
    cfg.trim.dtQS = 1e-4;
    cfg.aero.discrete = false;   % false | true → use obj.ROM_dsc
    % ----------- GUST SETTINGS ---------------------------------------------
    cfg.gust.on      = true;  % include Dryden gust input
    cfg.gust.gust_length= 10.0;
    cfg.gust.gust_intensity = 0.02;
    cfg.gust.gust_offset= 2;
     cfg.gust.gust_offset_ratio = .15; %start at 15% of sim


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
    

end
