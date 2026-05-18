% function cfg = GLAConfig(cfg)
% cfg          = nominalConfig;          % inherit old fields
cfg.ctrl.case     = 'GLA';                  % activates closed‑loop
% cfg.ctrl.Ts      = 0.0025;
% cfg.ctrl.Ts      = 0.005;
cfg.ctrl.Ts      = 0.01;
% cfg.ctrl.Ts      = 0.025;


cfg.forceRealGust = false;             % forces gust input to be true gust (perfect est)
cfg.forceRealSense = true;            % forces sensor data to be true data (does not load)                             
% horizon sizes ----------------------------------------------------------
% cfg.ctrl.Nc       = 24;                     % MPC steps  (0.06 s)
cfg.ctrl.Ne       = 8;                     % MHE steps  (0.06 s)
cfg.ctrl.Nc       = 10;                     
% cfg.ctrl.Nc       = 12;                     
% cfg.ctrl.Nc       = 4;                     
% cfg.ctrl.Ne       = 4;                     

 cfg.ctrl.mhe_method = 'multiple';   % 'multiple' or 'single'
% beam length needed by sensor ------------------------------------------
cfg.struct.L = 0.55;                   % [m] Pazy wing span

% tuning weights (Artola Table V‑1)  -------------------------------------
Nm  = cfg.Nm;
% cfg.Qc = blkdiag(eye(Nm),zeros(2*Nm+1+3));   % velocity states only
% cfg.Pc = zeros(2*Nm);
% cfg.Rc = 1e-3;                         % actuator penalty

% cfg.Qe = 0.1*eye(5*6);
% cfg.Qe = 0.5*eye(5*6);
% cfg.Qe = 1*eye(5*6);
% cfg.Qe(4:6,:) = 0;
% cfg.Qe(10:12 ,:) = 0;
% cfg.Qe(16:18 ,:) = 0;
% cfg.Qe(22:24 ,:) = 0;
% cfg.Qe(28:30 ,:) = 0;

% cfg.Qe = 1e-1*eye(4*6);
% cfg.Qe = 1e-2*eye(4*6);
% cfg.Qe = 2*eye(4*6);
% cfg.Qe = 5*eye(4*6);
% cfg.Qe = .5*eye(4*6);
% cfg.Qe = .85*eye(4*6);
cfg.Qe = eye(4*6);
cfg.Qe(4:6,:) = 0;
cfg.Qe(10:12 ,:) = 0;
cfg.Qe(16:18 ,:) = 0;
cfg.Qe(22:24 ,:) = 0;

x_len = 3*cfg.Nm+1+cfg.Na+3;

% cfg.Pe =  blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe =  2*blkdiag(zeros(Nm),eye(x_len-Nm));
cfg.Pe = 1* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe = .1* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe = 1e-2* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe = 5e-2* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe = .2* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Pe = .05* blkdiag(zeros(Nm),eye(x_len-Nm));
% cfg.Re = 5e-3;
% cfg.Re = 1e-3;
cfg.Re = 1e-2; % This is the best I believe
% cfg.Re = 1e-1;
% cfg.Re = 1;

% cfg.RdWe = 1e-2*eye(1);     % start
cfg.RdWe = 5e-3*eye(1);     % start
% cfg.RdWe = 1e-2*eye(cfg.nw);     % start
% cfg.RdWe = 1e-1*eye(cfg.nw);   % smoother

% cfg.mhe.penalizeAllW = true;
cfg.mhe.penalizeAllW = false;
cfg.mhe.RwTerminal = 1e-3;
%% MPC 

cfg.Qc = .5*blkdiag(eye(2*cfg.Nm),zeros(cfg.Nm+1+cfg.Na+3));   % velocity states & Force states only      Current this
% cfg.Qc = 2*blkdiag(eye(2*cfg.Nm),zeros(cfg.Nm+1+cfg.Na+3));   % velocity states & Force states only
% cfg.Qc = .5*blkdiag(eye(2*cfg.Nm),zeros(cfg.Nm+1+cfg.Na+3));   % velocity states & Force states only


% cfg.Qc = 5*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Qc = 2*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Qc = .1*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Qc = blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Qc = 0.85*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Qc = .5*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only
% cfg.Pc = cfg.Qc;
cfg.Pc = 2*cfg.Qc;
% cfg.Pc = 2.5*cfg.Qc; % Current this
% cfg.Pc = .5*cfg.Qc;
% cfg.Pc = 5*cfg.Qc;
% cfg.Pc = 0*cfg.Qc;
% cfg.Pc = 2.5*cfg.Qc;
% cfg.Pc = .1*cfg.Qc;

% cfg.Pc = 5*blkdiag(eye(cfg.Nm),zeros(2*cfg.Nm+1+cfg.Na+3));   % velocity states only

% cfg.Pc = 2.5*cfg.Qc;
% cfg.Rc = 1e-2;                          % actuator penalty
% cfg.Rc = 1e-3;        % Current this                 % actuator penalty
% cfg.Rc = diag([5e-4; 5e-4;5e-4;5e-4]);                         % actuator penalty
cfg.Rc = diag([5e-3; 5e-3 ;5e-3;5e-3]);                         % actuator penalty
% cfg.Rc = 1e-4;                         % actuator penalty
% cfg.Rc = 5e-5;                         % actuator penalty
% cfg.Rc = 5e-4;                         % actuator penalty
% % cfg.Rc = 1e-1;                         % actuator penalty
% cfg.Rc = 1;                         % actuator penalty
% cfg.Rc = 5e-3;                         % actuator penalty
% cfg.ctrl.aT = 1e-1;    % Current Tthis
cfg.ctrl.aT = 1e-2;
% cfg.ctrl.aT = 1e-6;
% % cfg.ctrl.aT = 5e-4;
% cfg.ctrl.aT = 5e-5;
% cfg.ctrl.aT = 1e-3;
% cfg.ctrl.aT = 1e-2;
% cfg.ctrl.aT = 1e-1;


cfg.ctrl.actuatorDeflectionAlpha = 0.5; 
% control limits ---------------------------------------------------------
cfg.uL = deg2rad(-15);
cfg.uU = deg2rad( 15);
cfg.urateL = deg2rad(-100);
cfg.urateU = deg2rad( 100);
% cfg.urateL = deg2rad(-60);  % Current this
% cfg.urateU = deg2rad( 60);  % Current this
% cfg.urateL = deg2rad(-50);
% cfg.urateU = deg2rad( 50);
% cfg.wL = -1;  % Need to edit these vals
cfg.wL = -1e-6;  % Need to edit these vals
cfg.wU = 10 ;
% cfg.wL = -100;  % Need to edit these vals
% cfg.wU = 100 ;
cfg.ctrl.n_surf  = 2;     % e.g., two aileron flaps
cfg.ctrl.var_per = 2;     % δ and δ̇ per surface
% solver options ---------------------------------------------------------
cfg.sqpOpts = struct('maxIter',20,'tol',1e-6);
cfg.sim.storeSens = true;

% the following three lines are ONLY handles; they receive (beam,cfg) later
cfg.sensorHandle    = @(beam,cfg) AeroFlex.sensor.WingVelSensor(beam,cfg);
cfg.estimatorHandle = @(cfg)      AeroFlex.ctrl.nMHE(cfg);
cfg.controllerHandle= @(cfg)      AeroFlex.ctrl.nMPC(cfg);
% end


% Debug: this should be kept as false

cfg.ctrl.skipMHE = false;
cfg.ctrl.skipMPC = false;
