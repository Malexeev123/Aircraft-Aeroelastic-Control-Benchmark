%======================================================================
% Custom SQP NMPC solver
%======================================================================
% cfg.ctrl.mpcSolver = "custom_sqp";   % "fmincon" | "custom_sqp"
cfg.ctrl.mpcSolver = "fmincon";   % "fmincon" | "custom_sqp"

cfg.ctrl.sqp = struct();

cfg.ctrl.sqp.MaxIterations       = 5;
cfg.ctrl.sqp.OptimalityTolerance = 1e-4;
cfg.ctrl.sqp.ConstraintTolerance = 1e-6;
cfg.ctrl.sqp.StepTolerance       = 1e-7;

cfg.ctrl.sqp.Display = "none";       % use "none" for speed runs

% Keep elastic mode on while debugging. Turn off later if the QPs are
% consistently feasible.
cfg.ctrl.sqp.ElasticMode        = true;
cfg.ctrl.sqp.ElasticPenaltyEq   = 1e5;
cfg.ctrl.sqp.ElasticPenaltyIneq = 1e5;

% BFGS matches the Artola-style SQP description.
cfg.ctrl.sqp.HessianMode = "bfgs";
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

%%
%======================================================================
% nMHE solver selection
%======================================================================
cfg.ctrl.mheSolver = "custom_sqp";   % "fmincon" | "custom_sqp"
% cfg.ctrl.mheSolver = "fmincon";   % "fmincon" | "custom_sqp"

cfg.ctrl.mheSqp = struct();

% RTI-like settings?
% cfg.ctrl.mheSqp.MaxIterations       = 1;
cfg.ctrl.mheSqp.MaxIterations       = 10;
% cfg.ctrl.mheSqp.OptimalityTolerance = 1e-3;
cfg.ctrl.mheSqp.OptimalityTolerance = 5e-3;
cfg.ctrl.mheSqp.ConstraintTolerance = 1e-6;
cfg.ctrl.mheSqp.StepTolerance       = 1e-7;

cfg.ctrl.mheSqp.Display = "none";

% Since your MHE continuity is already ~1e-9, elastic slacks are not needed.
cfg.ctrl.mheSqp.ElasticMode = false;

% Start cheap. BFGS can be enabled after the estimator is behaving.
cfg.ctrl.mheSqp.HessianMode = "bfgs"; % identity | bfgs
% cfg.ctrl.mheSqp.HessianMode = "identity"; % identity | bfgs
% cfg.ctrl.mheSqp.InitialHessianScale = 1.0;
% cfg.ctrl.mheSqp.HessianRegularization = 1e-8;


cfg.ctrl.mheSqp.CheckGradientsOnce = true;
cfg.ctrl.mheSqp.GradientCheckDirections = 12;
cfg.ctrl.mheSqp.GradientCheckStep = 1e-6;

% % Keep line-search permissive for RTI-like behavior.
cfg.ctrl.mheSqp.MeritPenalty = 1e3;
cfg.ctrl.mheSqp.MaxLineSearch = 20;
cfg.ctrl.mheSqp.LineSearchBeta = 0.5;
cfg.ctrl.mheSqp.LineSearchC1 = 1e-4;
cfg.ctrl.mheSqp.MeritAcceptTolerance = 1e-6;

cfg.ctrl.mheSqp.QPOptions = optimoptions('quadprog', ...
    'Algorithm','interior-point-convex', ...
    'Display','off', ...
    'OptimalityTolerance',1e-7, ...
    'ConstraintTolerance',1e-7, ...
    'StepTolerance',1e-9, ...
    'MaxIterations',100);