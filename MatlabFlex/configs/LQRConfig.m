%======================================================================
% OUTER LQR CONFIGURATION
%======================================================================

% Enable/disable outer-loop.
cfg.outerLQR.enable = true;

% Outer loop should be slower than the nMHE/nMPC inner loop.
cfg.ctrl.outerEvery = 4;                     % update LQR every 4 inner samples
cfg.ctrl.outerTs    = cfg.ctrl.outerEvery * cfg.ctrl.Ts;

% Use discrete-time LQR at the outer-loop update rate.
% cfg.outerLQR.useDiscrete = false;
cfg.outerLQR.useDiscrete = true;

% Rigid-body state used by OuterLQR:
%   x_rb = [u; v; w; phi; theta; psi; p; q; r]
cfg.outerLQR.stateNames = {'u','v','w','phi','theta','psi','p','q','r'};

% Outer LQR command vector:
%   u_rb = [delta_e; delta_a; delta_r; thrust]
%
% delta_e : elevator / symmetric pitch surface
% delta_a : differential aileron / roll command
% delta_r : rudder / yaw command
% thrust  : axial thrust command
cfg.outerLQR.inputNames = {'delta_e','delta_a','delta_r','thrust'};

% Trim rigid-body input.
cfg.outerLQR.uTrim = [ ...
    0;                              % delta_e [rad]
    trim.deltaDeg(1);                              % delta_a [rad]
    0;                              % delta_r [rad]
    getfield_safe_trim_thrust(cfg, trim)];  % thrust [N]

% If the linear model is not supplied, OuterLQR will finite-difference its
% internal rigid-body model.
% Optional:
% cfg.outerLQR.A = A_rb;
% cfg.outerLQR.B = B_rb;

% Aileron roll-moment gain used by fallback finite-difference model.
% Should tune/identify this from coupled model or SHARPy.
% Units: [N*m/rad]
cfg.outerLQR.aileronRollMomentGain = 0.02;

% Optional aerodynamic yaw/roll damping guesses for fallback model.
cfg.outerLQR.Lp = -0.01;    % roll damping moment per p, [N*m/(rad/s)]
cfg.outerLQR.Nr = -0.01;    % yaw damping moment per r,  [N*m/(rad/s)]
cfg.outerLQR.Mq = -0.02;    % pitch damping moment per q,[N*m/(rad/s)]

% Bryson-style maximum acceptable rigid-body errors.
% These determine Q = diag(1/xMax.^2).
cfg.outerLQR.xMax = [ ...
    2.0;                    % u error [m/s]
    0.5;                    % v error [m/s]
    1.0;                    % w error [m/s]
    deg2rad(10);             % roll error [rad]
    deg2rad(8);              % pitch error [rad]
    deg2rad(15);             % yaw error [rad]
    deg2rad(20);             % p error [rad/s]
    deg2rad(20);             % q error [rad/s]
    deg2rad(20)];            % r error [rad/s]

% Maximum allowed rigid-body control effort.
% These determine R = diag(1/uMax.^2).
cfg.outerLQR.uMax = [ ...
    deg2rad(10);             % elevator [rad]
    deg2rad(8);              % aileron differential command [rad]
    deg2rad(10);             % rudder [rad]
    0.5];                    % thrust perturbation [N]

% Optional manual Q/R override:
% cfg.outerLQR.Q = diag([...]);
% cfg.outerLQR.R = diag([...]);

% Shared wing-surface allocation for outer differential aileron command.
% For two wing surfaces:
%   delta_a > 0 maps to [left; right] = [+delta_a; -delta_a]
cfg.outerLQR.aileronSurfaceMap = [1; -1];

% Optional symmetric elevator-like wing command if using shared symmetric surfaces.
cfg.outerLQR.elevatorSurfaceMap = [1; 1];

% Use these maps only if those commands should enter the wing ROM surface vector.
cfg.outerLQR.useAileronOnWingSurfaces  = true;
cfg.outerLQR.useElevatorOnWingSurfaces = false;

%======================================================================
% COMMAND FUSION
%======================================================================
cfg.fusion.enable = true;

% Crossover below the first flexible mode and above rigid-body bandwidth.
% Outer bandwidth stays below wing modes.
% Start around 2-3 rad/s if first flexible mode is ~30 rad/s.
cfg.fusion.omega_c = 2.5;     % [rad/s]

% If nMPC returns a total command around trim, keep false.
% If nMPC returns only correction delta_u, set true.
cfg.fusion.innerIsDeviation = false;

%======================================================================
% ACTUATOR CHAIN
%======================================================================
cfg.actuator.enable = true;
cfg.actuator.omega = 40;                 % first-order servo bandwidth [rad/s]
cfg.actuator.deltaL = cfg.uL;
cfg.actuator.deltaU = cfg.uU;
cfg.actuator.rateL  = cfg.urateL;
cfg.actuator.rateU  = cfg.urateU;


% cfg.outerLQR.designStateIdx = [4 5 7 8];   % phi, p
% cfg.outerLQR.designInputIdx = [2 1];     % delta_a

cfg.outerLQR.designStateIdx = [4 7];   % phi, p
cfg.outerLQR.designInputIdx = [2];     % delta_a

cfg.outerLQR.aileronRollMomentGain = 0.05;
% cfg.outerLQR.elevatorPitchMomentGain = -0.05;
cfg.outerLQR.Lp = -0.02;
% cfg.outerLQR.Mp = -0.04;

function T = getfield_safe_trim_thrust(cfg, trim)
    if exist('trim','var') && isstruct(trim) && isfield(trim,'thrust')
        T = trim.thrust;
    elseif isfield(cfg,'trim') && isfield(cfg.trim,'thrust')
        T = cfg.trim.thrust;
    else
        T = 0;
    end
end