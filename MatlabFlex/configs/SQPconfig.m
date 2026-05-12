%======================================================================
% Custom SQP NMPC solver
%======================================================================
cfg.ctrl.mpcSolver = "custom_sqp";   % "fmincon" | "custom_sqp"

cfg.ctrl.sqp = struct();

cfg.ctrl.sqp.MaxIterations       = 15;
cfg.ctrl.sqp.OptimalityTolerance = 1e-4;
cfg.ctrl.sqp.ConstraintTolerance = 1e-6;
cfg.ctrl.sqp.StepTolerance       = 1e-7;

cfg.ctrl.sqp.Display = "iter";       % use "none" for speed runs

% Keep elastic mode on while debugging. Turn off later if the QPs are
% consistently feasible.
cfg.ctrl.sqp.ElasticMode        = true;
cfg.ctrl.sqp.ElasticPenaltyEq   = 1e5;
cfg.ctrl.sqp.ElasticPenaltyIneq = 1e5;

% BFGS matches the Artola-style SQP description.
cfg.ctrl.sqp.HessianMode = "bfgs";

% Use identity if BFGS dense Hessian becomes too slow:
% cfg.ctrl.sqp.HessianMode = "identity";

cfg.ctrl.sqp.CheckGradientsOnce = true;
cfg.ctrl.sqp.GradientCheckDirections = 12;
cfg.ctrl.sqp.GradientCheckStep = 1e-6;

cfg.ctrl.sqp.QPOptions = optimoptions('quadprog', ...
    'Algorithm','interior-point-convex', ...
    'Display','off', ...
    'OptimalityTolerance',1e-8, ...
    'ConstraintTolerance',1e-8, ...
    'StepTolerance',1e-10, ...
    'MaxIterations',200);