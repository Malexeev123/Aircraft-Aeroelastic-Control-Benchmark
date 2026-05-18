%% test_custom_sqp_qualification.m
%==========================================================================
% CUSTOM SQP QUALIFICATION TESTS
%==========================================================================
% Purpose:
%   Qualify AeroFlex.optim.SQPSolver independently from the full aeroelastic
%   nMPC problem.
%
% Tests:
%   1. Linear first-order MPC multiple-shooting problem.
%      - Convex QP with linear constraints.
%      - Compare custom SQP to direct quadprog solution.
%      - With exact Hessian inserted, SQP should converge in one/few steps.
%
%   2. Nonlinear second-order oscillator MPC multiple-shooting problem.
%      - Nonlinear dynamics with analytic Jacobians.
%      - Compare custom SQP to fmincon if available.
%      - This test is structurally close to nMPC:
%            z = [X_0; ...; X_N; U_0; ...; U_{N-1}]
%            X_0 = x_meas
%            X_{k+1} - Phi(X_k,U_k) = 0
%
% Required:
%   +AeroFlex/+optim/SQPSolver.m
%   Optimization Toolbox for quadprog; fmincon comparison is optional.
%==========================================================================
clear classes;
rehash;
clear; clc; close all;

addpath(genpath('MatlabFlex\')); addpath(genpath('TestBenchPazy\'));
addpath(genpath('plots\'));

fprintf('\n');
fprintf('======================================================================\n');
fprintf(' CUSTOM SQP QUALIFICATION SUITE\n');
fprintf('======================================================================\n');

assert(exist('AeroFlex.optim.SQPSolver','class') == 8, ...
    ['Could not find AeroFlex.optim.SQPSolver. Make sure ', ...
     '+AeroFlex/+optim/SQPSolver.m is on the MATLAB path.']);

assert(exist('quadprog','file') == 2, ...
    'quadprog is required because SQPSolver uses quadprog for QP subproblems.');

hasFmincon = exist('fmincon','file') == 2;

if hasFmincon
    fprintf('[setup] fmincon found. Nonlinear test will compare against fmincon.\n');
else
    fprintf('[setup] fmincon not found. Nonlinear test will run custom SQP only.\n');
end

%% ========================================================================
% TEST 1: LINEAR FIRST-ORDER MPC/QP
%==========================================================================
fprintf('\n');
fprintf('======================================================================\n');
fprintf(' TEST 1: LINEAR FIRST-ORDER MPC / EXACT QP\n');
fprintf('======================================================================\n');

linCfg = struct();
linCfg.N  = 20;
linCfg.dt = 0.1;

% Continuous first-order system:
%   xdot = -a*x + b*u
%
% Exact ZOH discretization:
%   x_{k+1} = Ad*x_k + Bd*u_k
linCfg.a  = 1.2;
linCfg.b  = 1.0;
linCfg.Ad = exp(-linCfg.a*linCfg.dt);
linCfg.Bd = (1 - exp(-linCfg.a*linCfg.dt)) * linCfg.b / linCfg.a;

linCfg.x0   = 1.0;
linCfg.xRef = 0.0;
linCfg.uRef = 0.0;

linCfg.Q = 5.0;
linCfg.R = 0.05;
linCfg.P = 20.0;

linCfg.uMax = 2.0;

nlpLin = assembleLinearFirstOrderMPC(linCfg);

% Feasible warm start from zero input rollout.
z0Lin = buildLinearFeasibleGuess(nlpLin,linCfg);

% One-time gradient check.
sqpOpts = AeroFlex.optim.SQPSolver.defaultOptions();
sqpOpts.Display = "iter";
sqpOpts.MaxIterations = 10;
sqpOpts.OptimalityTolerance = 1e-8;
sqpOpts.ConstraintTolerance = 1e-9;
sqpOpts.StepTolerance = 1e-10;
sqpOpts.ElasticMode = false;
sqpOpts.HessianMode = "identity";
sqpOpts.CheckGradientsOnce = false;
sqpOpts.GradientCheckDirections = 10;
sqpOpts.GradientCheckStep = 1e-6;

%% Testing heere
solverLin = AeroFlex.optim.SQPSolver(sqpOpts);

fprintf('\n[Test 1] Directional gradient check:\n');
gradReportLin = solverLin.checkGradients(nlpLin.cost,nlpLin.nonl,z0Lin,nlpLin.lb,nlpLin.ub);

% Inject exact Hessian. This makes the SQP QP subproblem the exact original
% QP for this linear-quadratic test.
solverLin.H = nlpLin.H;

[zSQP1,fSQP1,efSQP1,outSQP1,lambdaSQP1] = solverLin.solve( ...
    nlpLin.cost,z0Lin,nlpLin.lb,nlpLin.ub,nlpLin.nonl); 

% Direct quadprog reference solution.
qpOpts = optimoptions('quadprog', ...
    'Algorithm','interior-point-convex', ...
    'Display','off', ...
    'OptimalityTolerance',1e-10, ...
    'ConstraintTolerance',1e-10, ...
    'StepTolerance',1e-12);

[zQP1,~,efQP1,outQP1] = quadprog( ...
    nlpLin.H,nlpLin.h, ...
    [],[], ...
    nlpLin.Aeq,nlpLin.beq, ...
    nlpLin.lb,nlpLin.ub, ...
    [],qpOpts); 

% [fQP1,~] = nlpLin.cost(zQP2);
[fQP1,~] = nlpLin.cost(zQP1);

errZ1 = norm(zSQP1-zQP1,inf);
errF1 = abs(fSQP1-fQP1);
[~,ceqSQP1] = nlpLin.nonl(zSQP1);

fprintf('\n[Test 1 Results]\n');
fprintf('  SQP exitflag             : %+d\n', efSQP1);
fprintf('  SQP message              : %s\n', outSQP1.message);
fprintf('  QP reference exitflag    : %+d\n', efQP1);
fprintf('  ||z_SQP-z_QP||_inf       : %.3e\n', errZ1);
fprintf('  |J_SQP-J_QP|             : %.3e\n', errF1);
fprintf('  ||ceq_SQP||_inf          : %.3e\n', norm(ceqSQP1,inf));
fprintf('  SQP constr violation     : %.3e\n', outSQP1.constrviolation);
fprintf('  SQP first-order opt      : %.3e\n', outSQP1.firstorderopt);
fprintf('  SQP iterations           : %d\n', outSQP1.iterations);

pass1 = efSQP1 > 0 && errZ1 < 1e-6 && norm(ceqSQP1,inf) < 1e-8;

if pass1
    fprintf('  TEST 1 STATUS            : PASS\n');
else
    fprintf('  TEST 1 STATUS            : FAIL / INVESTIGATE\n');
end

plotLinearTest(nlpLin,linCfg,zSQP1,zQP1);

%% ========================================================================
% TEST 2: NONLINEAR SECOND-ORDER OSCILLATOR MPC
%==========================================================================
fprintf('\n');
fprintf('======================================================================\n');
fprintf(' TEST 2: NONLINEAR SECOND-ORDER OSCILLATOR MPC\n');
fprintf('======================================================================\n');

nlCfg = struct();
nlCfg.N  = 25;
nlCfg.dt = 0.05;

% Nonlinear oscillator:
%
%   p_dot = v
%   v_dot = -wn^2*p - 2*zeta*wn*v - cdrag*v*abs(v) + u
%
% Forward Euler discrete map:
%
%   p_{k+1} = p_k + dt*v_k
%   v_{k+1} = v_k + dt*(-wn^2*p_k - 2*zeta*wn*v_k
%                       - cdrag*v_k*abs(v_k) + u_k)
%
nlCfg.wn    = 2.0;
nlCfg.zeta  = 0.15;
nlCfg.cdrag = 0.25;

nlCfg.x0   = [1.0; 0.0];
nlCfg.xRef = [0.0; 0.0];
nlCfg.uRef = 0.0;

nlCfg.Q = diag([10.0, 0.5]);
nlCfg.R = 0.05;
nlCfg.P = diag([25.0, 2.0]);

nlCfg.uMax = 4.0;

nlpNL = assembleNonlinearSecondOrderMPC(nlCfg);

% Feasible initial guess from zero-input rollout.
z0NL = buildNonlinearFeasibleGuess(nlpNL,nlCfg);

sqpOpts2 = AeroFlex.optim.SQPSolver.defaultOptions();
sqpOpts2.Display = "iter";
% sqpOpts2.MaxIterations = 25;
sqpOpts2.MaxIterations = 100;
sqpOpts2.OptimalityTolerance = 1e-5;
sqpOpts2.ConstraintTolerance = 1e-8;
sqpOpts2.StepTolerance = 1e-9;
sqpOpts2.ElasticMode = false;
sqpOpts2.HessianMode = "bfgs";
sqpOpts2.InitialHessianScale = 1.0;
sqpOpts2.CheckGradientsOnce = false;
sqpOpts2.GradientCheckDirections = 20;
sqpOpts2.GradientCheckStep = 1e-6;
sqpOpts2.MeritPenalty = 1e3;

solverNL = AeroFlex.optim.SQPSolver(sqpOpts2);

fprintf('\n[Test 2] Directional gradient check:\n');
gradReportNL = solverNL.checkGradients(nlpNL.cost,nlpNL.nonl,z0NL,nlpNL.lb,nlpNL.ub);

[zSQP2,fSQP2,efSQP2,outSQP2,lambdaSQP2] = solverNL.solve( ...
    nlpNL.cost,z0NL,nlpNL.lb,nlpNL.ub,nlpNL.nonl);

[~,ceqSQP2] = nlpNL.nonl(zSQP2);

fprintf('\n[Test 2 SQP Results]\n');
fprintf('  SQP exitflag             : %+d\n', efSQP2);
fprintf('  SQP message              : %s\n', outSQP2.message);
fprintf('  SQP objective            : %.12e\n', fSQP2);
fprintf('  ||ceq_SQP||_inf          : %.3e\n', norm(ceqSQP2,inf));
fprintf('  SQP constr violation     : %.3e\n', outSQP2.constrviolation);
fprintf('  SQP first-order opt      : %.3e\n', outSQP2.firstorderopt);
fprintf('  SQP iterations           : %d\n', outSQP2.iterations);

zFMIN2 = [];
fFMIN2 = nan;
efFMIN2 = nan;
outFMIN2 = struct();

if hasFmincon
    fminOpts = optimoptions('fmincon', ...
        'Algorithm','sqp', ...
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true, ...
        'Display','none', ...
        'OptimalityTolerance',1e-8, ...
        'ConstraintTolerance',1e-9, ...
        'StepTolerance',1e-10, ...
        'MaxIterations',200, ...
        'MaxFunctionEvaluations',50000);

    [zFMIN2,fFMIN2,efFMIN2,outFMIN2] = fmincon( ...
        nlpNL.cost,z0NL, ...
        [],[],[],[], ...
        nlpNL.lb,nlpNL.ub, ...
        nlpNL.nonl,fminOpts);

    [~,ceqFMIN2] = nlpNL.nonl(zFMIN2);

    errZ2 = norm(zSQP2-zFMIN2,inf);
    errF2 = abs(fSQP2-fFMIN2);

    fprintf('\n[Test 2 fmincon Comparison]\n');
    fprintf('  fmincon exitflag         : %+d\n', efFMIN2);
    fprintf('  fmincon objective        : %.12e\n', fFMIN2);
    fprintf('  fmincon iterations       : %d\n', outFMIN2.iterations);
    fprintf('  ||ceq_fmincon||_inf      : %.3e\n', norm(ceqFMIN2,inf));
    fprintf('  ||z_SQP-z_fmincon||_inf  : %.3e\n', errZ2);
    fprintf('  |J_SQP-J_fmincon|        : %.3e\n', errF2);
end

pass2Grad = gradReportNL.objectiveRelErrMax < 1e-6 && ...
            gradReportNL.eqRelErrMax < 1e-5;

pass2Feas = efSQP2 > 0 && norm(ceqSQP2,inf) < 1e-7;

if hasFmincon
    pass2Compare = abs(fSQP2-fFMIN2) < 1e-4;
else
    pass2Compare = true;
end

pass2Exit = efSQP2 > 0 || ...
            (outSQP2.iterations >= sqpOpts2.MaxIterations && pass2Feas && pass2Compare);

pass2 = pass2Grad && pass2Feas && pass2Compare && pass2Exit;

% pass2 = pass2Grad && pass2Feas && pass2Compare;

fprintf('\n[Test 2 Pass Criteria]\n');
fprintf('  gradient pass            : %d\n', pass2Grad);
fprintf('  feasibility pass         : %d\n', pass2Feas);
fprintf('  fmincon comparison pass  : %d\n', pass2Compare);
if pass2
    fprintf('  TEST 2 STATUS            : PASS\n');
else
    fprintf('  TEST 2 STATUS            : FAIL / INVESTIGATE\n');
end

plotNonlinearTest(nlpNL,nlCfg,zSQP2,zFMIN2);

%% ========================================================================
% FINAL SUMMARY
%==========================================================================
fprintf('\n');
fprintf('======================================================================\n');
fprintf(' CUSTOM SQP QUALIFICATION SUMMARY\n');
fprintf('======================================================================\n');
fprintf('  Test 1 linear exact-QP status      : %s\n', passFailString(pass1));
fprintf('  Test 2 nonlinear MPC status        : %s\n', passFailString(pass2));
fprintf('----------------------------------------------------------------------\n');
fprintf('  Test 1 objective grad error        : %.3e\n', gradReportLin.objectiveRelErrMax);
fprintf('  Test 1 equality grad error         : %.3e\n', gradReportLin.eqRelErrMax);
fprintf('  Test 2 objective grad error        : %.3e\n', gradReportNL.objectiveRelErrMax);
fprintf('  Test 2 equality grad error         : %.3e\n', gradReportNL.eqRelErrMax);
fprintf('======================================================================\n');

if pass1 && pass2
    fprintf('\nQualification result: SQP core is behaving correctly on simple MPC problems.\n');
    fprintf('If the aeroelastic nMPC still has high equality-gradient error, the issue is likely in ROM sensitivity/Jacobian assembly, not the SQP wrapper.\n\n');
else
    fprintf('\nQualification result: SQP core needs debugging before returning to aeroelastic nMPC.\n\n');
end

%% ========================================================================
% LOCAL FUNCTIONS
%==========================================================================

function nlp = assembleLinearFirstOrderMPC(cfg)
%ASSEMBLELINEARFIRSTORDERMPC Linear first-order MPC problem.
%
% Decision vector:
%   z = [x_0; ...; x_N; u_0; ...; u_{N-1}]
%
% Equality constraints:
%   x_0 = x_meas
%   x_{k+1} - Ad*x_k - Bd*u_k = 0
%
% Cost:
%   sum_{k=0}^{N-1} 1/2 Q (x_k-xRef)^2
%                 + 1/2 R (u_k-uRef)^2
%   + 1/2 P (x_N-xRef)^2

    N = cfg.N;
    nX = 1;
    nU = 1;

    nVar = (N+1)*nX + N*nU;

    idx.x = cell(N+1,1);
    for k = 1:N+1
        idx.x{k} = k;
    end

    idx.u = cell(N,1);
    startU = N+1;
    for k = 1:N
        idx.u{k} = startU + k;
    end

    H = sparse(nVar,nVar);
    h = zeros(nVar,1);

    for k = 1:N
        ix = idx.x{k};
        iu = idx.u{k};

        H(ix,ix) = H(ix,ix) + cfg.Q;
        h(ix) = h(ix) - cfg.Q*cfg.xRef;

        H(iu,iu) = H(iu,iu) + cfg.R;
        h(iu) = h(iu) - cfg.R*cfg.uRef;
    end

    ixN = idx.x{N+1};
    H(ixN,ixN) = H(ixN,ixN) + cfg.P;
    h(ixN) = h(ixN) - cfg.P*cfg.xRef;

    H = 0.5*(H+H.');

    % Direct linear equality matrices for quadprog reference.
    nEq = N + 1;
    Aeq = sparse(nEq,nVar);
    beq = zeros(nEq,1);

    Aeq(1,idx.x{1}) = 1;
    beq(1) = cfg.x0;

    row = 2;
    for k = 1:N
        Aeq(row,idx.x{k+1}) = 1;
        Aeq(row,idx.x{k})   = -cfg.Ad;
        Aeq(row,idx.u{k})   = -cfg.Bd;
        row = row + 1;
    end

    lb = -inf(nVar,1);
    ub =  inf(nVar,1);

    for k = 1:N
        lb(idx.u{k}) = -cfg.uMax;
        ub(idx.u{k}) =  cfg.uMax;
    end

    function [J,g] = costFun(z)
        z = z(:);

        J = 0.5*z.'*H*z + h.'*z;
        g = H*z + h;
    end

    function [c,ceq,gradc,gradceq] = nonlFun(z)
        z = z(:);

        c = [];
        gradc = [];

        ceq = Aeq*z - beq;
        gradceq = Aeq.';
    end

    nlp = struct();
    nlp.cost = @costFun;
    nlp.nonl = @nonlFun;
    nlp.lb = lb;
    nlp.ub = ub;
    nlp.H = H;
    nlp.h = h;
    nlp.Aeq = Aeq;
    nlp.beq = beq;
    nlp.idx = idx;
end

function z0 = buildLinearFeasibleGuess(nlp,cfg)
%BUILDLINEARFEASIBLEGUESS Zero-input feasible rollout.

    N = cfg.N;

    z0 = zeros(numel(nlp.lb),1);

    x = cfg.x0;
    z0(nlp.idx.x{1}) = x;

    for k = 1:N
        u = 0;
        z0(nlp.idx.u{k}) = u;

        x = cfg.Ad*x + cfg.Bd*u;
        z0(nlp.idx.x{k+1}) = x;
    end
end

function nlp = assembleNonlinearSecondOrderMPC(cfg)
%ASSEMBLENONLINEARSECONDORDERMPC Nonlinear second-order MPC.
%
% Dynamics:
%   p_{k+1} = p_k + dt*v_k
%   v_{k+1} = v_k + dt*(-wn^2*p_k - 2*zeta*wn*v_k
%                       - cdrag*v_k*abs(v_k) + u_k)
%
% Decision vector:
%   z = [X_0; ...; X_N; U_0; ...; U_{N-1}]
%
% where X_k = [p_k; v_k].

    N = cfg.N;
    nx = 2;
    nu = 1;

    nVar = (N+1)*nx + N*nu;

    idx.x = cell(N+1,1);
    for k = 1:N+1
        idx.x{k} = (k-1)*nx + (1:nx);
    end

    startU = (N+1)*nx;
    idx.u = cell(N,1);
    for k = 1:N
        idx.u{k} = startU + (k-1)*nu + 1;
    end

    lb = -inf(nVar,1);
    ub =  inf(nVar,1);

    for k = 1:N
        lb(idx.u{k}) = -cfg.uMax;
        ub(idx.u{k}) =  cfg.uMax;
    end

    Q = 0.5*(cfg.Q+cfg.Q.');
    P = 0.5*(cfg.P+cfg.P.');
    R = cfg.R;

    xRef = cfg.xRef(:);
    uRef = cfg.uRef;

    function [J,g] = costFun(z)
        z = z(:);

        J = 0;
        g = zeros(nVar,1);

        for k = 1:N
            xk = z(idx.x{k});
            uk = z(idx.u{k});

            dx = xk - xRef;
            du = uk - uRef;

            J = J + 0.5*dx.'*Q*dx + 0.5*R*du^2;

            g(idx.x{k}) = g(idx.x{k}) + Q*dx;
            g(idx.u{k}) = g(idx.u{k}) + R*du;
        end

        xN = z(idx.x{N+1});
        dxN = xN - xRef;

        J = J + 0.5*dxN.'*P*dxN;
        g(idx.x{N+1}) = g(idx.x{N+1}) + P*dxN;
    end

    function [c,ceq,gradc,gradceq] = nonlFun(z)
        z = z(:);

        c = [];
        gradc = [];

        nEq = nx*(N+1);

        ceq = zeros(nEq,1);
        Jceq = spalloc(nEq,nVar,nx + N*(nx + nx*nx + nx*nu));

        % Initial condition:
        %   X0 - x_meas = 0
        rows = 1:nx;
        ceq(rows) = z(idx.x{1}) - cfg.x0(:);
        Jceq(rows,idx.x{1}) = speye(nx);

        row0 = nx;

        for k = 1:N
            xk = z(idx.x{k});
            uk = z(idx.u{k});

            [xNext,A,B] = nonlinearStepAndJacobians(xk,uk,cfg);

            rows = row0 + (1:nx);

            ceq(rows) = z(idx.x{k+1}) - xNext;

            Jceq(rows,idx.x{k+1}) = speye(nx);
            Jceq(rows,idx.x{k})   = -A;
            Jceq(rows,idx.u{k})   = -B;

            row0 = row0 + nx;
        end

        gradceq = Jceq.';
    end

    nlp = struct();
    nlp.cost = @costFun;
    nlp.nonl = @nonlFun;
    nlp.lb = lb;
    nlp.ub = ub;
    nlp.idx = idx;
end

function z0 = buildNonlinearFeasibleGuess(nlp,cfg)
%BUILDNONLINEARFEASIBLEGUESS Feasible zero-input rollout for nonlinear test.

    N = cfg.N;

    z0 = zeros(numel(nlp.lb),1);

    x = cfg.x0(:);
    z0(nlp.idx.x{1}) = x;

    for k = 1:N
        u = 0;
        z0(nlp.idx.u{k}) = u;

        [x,~,~] = nonlinearStepAndJacobians(x,u,cfg);
        z0(nlp.idx.x{k+1}) = x;
    end
end

function [xNext,A,B] = nonlinearStepAndJacobians(x,u,cfg)
%NONLINEARSTEPANDJACOBIANS Forward Euler nonlinear oscillator map.

    p = x(1);
    v = x(2);

    dt = cfg.dt;
    wn = cfg.wn;
    zeta = cfg.zeta;
    cdrag = cfg.cdrag;

    acc = -wn^2*p - 2*zeta*wn*v - cdrag*v*abs(v) + u;

    xNext = [ ...
        p + dt*v;
        v + dt*acc];

    % d/dv of v*abs(v) = 2*abs(v), except at v=0 where the derivative is
    % not classical. The central finite-difference derivative at v=0 is 0,
    % and 2*abs(0)=0, so this expression is consistent enough for testing.
    dDrag_dv = 2*abs(v);

    A = [ ...
        1,              dt;
       -dt*wn^2,        1 + dt*(-2*zeta*wn - cdrag*dDrag_dv)];

    B = [0; dt];
end

function plotLinearTest(nlp,cfg,zSQP,zQP)
%PLOTLINEARTEST Plot linear MPC qualification result.

    N = cfg.N;
    tX = (0:N)*cfg.dt;
    tU = (0:N-1)*cfg.dt;

    xSQP = zeros(1,N+1);
    uSQP = zeros(1,N);

    xQP = zeros(1,N+1);
    uQP = zeros(1,N);

    for k = 1:N+1
        xSQP(k) = zSQP(nlp.idx.x{k});
        xQP(k) = zQP(nlp.idx.x{k});
    end

    for k = 1:N
        uSQP(k) = zSQP(nlp.idx.u{k});
        uQP(k) = zQP(nlp.idx.u{k});
    end

    figure('Color','w','Name','Test 1 Linear MPC');
    tiledlayout(2,1,'TileSpacing','tight');

    nexttile;
    plot(tX,xSQP,'LineWidth',2); hold on;
    plot(tX,xQP,'--','LineWidth',2);
    yline(cfg.xRef,':','LineWidth',1.5);
    grid on;
    ylabel('x');
    legend('custom SQP','quadprog','reference','Location','best');
    title('Test 1: Linear first-order MPC');

    nexttile;
    stairs(tU,uSQP,'LineWidth',2); hold on;
    stairs(tU,uQP,'--','LineWidth',2);
    yline(cfg.uMax,':','LineWidth',1.5);
    yline(-cfg.uMax,':','LineWidth',1.5);
    grid on;
    xlabel('t [s]');
    ylabel('u');
    legend('custom SQP','quadprog','bounds','Location','best');
end

function plotNonlinearTest(nlp,cfg,zSQP,zFMIN)
%PLOTNONLINEARTEST Plot nonlinear oscillator MPC result.

    N = cfg.N;
    tX = (0:N)*cfg.dt;
    tU = (0:N-1)*cfg.dt;

    xSQP = zeros(2,N+1);
    uSQP = zeros(1,N);

    for k = 1:N+1
        xSQP(:,k) = zSQP(nlp.idx.x{k});
    end

    for k = 1:N
        uSQP(k) = zSQP(nlp.idx.u{k});
    end

    hasFmin = ~isempty(zFMIN);

    if hasFmin
        xFM = zeros(2,N+1);
        uFM = zeros(1,N);

        for k = 1:N+1
            xFM(:,k) = zFMIN(nlp.idx.x{k});
        end

        for k = 1:N
            uFM(k) = zFMIN(nlp.idx.u{k});
        end
    end

    figure('Color','w','Name','Test 2 Nonlinear MPC');
    tiledlayout(3,1,'TileSpacing','tight');

    nexttile;
    plot(tX,xSQP(1,:),'LineWidth',2); hold on;
    if hasFmin
        plot(tX,xFM(1,:),'--','LineWidth',2);
    end
    yline(cfg.xRef(1),':','LineWidth',1.5);
    grid on;
    ylabel('p');
    title('Test 2: Nonlinear second-order oscillator MPC');
    if hasFmin
        legend('custom SQP','fmincon','reference','Location','best');
    else
        legend('custom SQP','reference','Location','best');
    end

    nexttile;
    plot(tX,xSQP(2,:),'LineWidth',2); hold on;
    if hasFmin
        plot(tX,xFM(2,:),'--','LineWidth',2);
    end
    yline(cfg.xRef(2),':','LineWidth',1.5);
    grid on;
    ylabel('v');

    nexttile;
    stairs(tU,uSQP,'LineWidth',2); hold on;
    if hasFmin
        stairs(tU,uFM,'--','LineWidth',2);
    end
    yline(cfg.uMax,':','LineWidth',1.5);
    yline(-cfg.uMax,':','LineWidth',1.5);
    grid on;
    xlabel('t [s]');
    ylabel('u');
end

function s = passFailString(flag)
    if flag
        s = 'PASS';
    else
        s = 'FAIL';
    end
end