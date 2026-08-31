classdef SQPSolver < handle
%======================================================================
% SQPSOLVER
%======================================================================
% Lightweight Sequential Quadratic Programming solver for NMPC/NMHE
% multiple-shooting NLPs of the form:
%
%   min_z     f(z)
%   s.t.      c(z)   <= 0
%             ceq(z) == 0
%             lb <= z <= ub
%
% Required function signatures:
%
%   [f,g] = costFun(z)
%
%   [c,ceq,gradc,gradceq] = nonlFun(z)
%
% where gradients follow MATLAB fmincon convention:
%
%   gradc    : nVar x nIneq
%   gradceq  : nVar x nEq
%
% Each SQP iteration solves:
%
%   min_d  g' d + 1/2 d' H d
%   s.t.   ceq + J_eq d = 0
%          c   + J_in d <= 0
%          lb-z <= d <= ub-z
%
% Optional elastic slacks are included to avoid QP infeasibility during
% early iterations.
%
% Exit flags:
%    1 : first-order and feasibility tolerances satisfied
%    2 : step tolerance satisfied and feasible
%    0 : maximum iterations reached
%   -1 : QP subproblem failed
%   -2 : line search failed
%   -3 : invalid function/gradient evaluation
%
% This is not intended to be a universal replacement for fmincon. It is
% intended for repeatedly solved warm-started NMPC/NMHE multiple-shooting
% problems with user-supplied gradients.
%======================================================================

    properties
        options struct
        H
        lastOutput struct
        lastLambda struct
        lastGradCheck struct
        lastAuditCapture struct
        lastPreparedInitial struct
        lbfgsSHistory cell
        lbfgsYHistory cell
    end

    methods(Static)
        function opts = defaultOptions()
            opts = struct();

            % SQP major iterations.
            opts.MaxIterations       = 25;
            opts.OptimalityTolerance = 1e-4;
            opts.ConstraintTolerance = 1e-6;
            opts.StepTolerance       = 1e-7;

            % Merit line-search.
            opts.MeritPenalty        = 1e3;
            opts.LineSearchBeta      = 0.5;
            opts.LineSearchC1        = 1e-4;
            opts.MaxLineSearch       = 15;
            opts.MinAlpha            = 1e-6;
            opts.MeritAcceptTolerance = 1e-8;
            % "fixed_beta" is the established path.  The audit-only
            % quadratic option uses the rejected merit value to choose a
            % safeguarded variable contraction instead of always halving.
            opts.LineSearchUpdate = "fixed_beta";

            % Default-disabled audit-only second-order correction (SOC).
            % A curved equality manifold can reject a useful tangential QP
            % step through the Maratos effect.  SOC restores nonlinear
            % feasibility at a trial point before applying the existing
            % acceptance test; it never relaxes the final NLP constraints.
            opts.SecondOrderCorrectionEnabled = false;
            opts.SecondOrderCorrectionMaxRelativeNorm = 1.0;

            % QP elastic mode. This makes the QP subproblem feasible even
            % when the linearized constraints are inconsistent.
            opts.ElasticMode         = true;
            opts.ElasticPenaltyEq    = 1e5;
            opts.ElasticPenaltyIneq  = 1e5;
            opts.ElasticSlackQuad    = 1e-8;
            opts.ElasticSlackTol     = 1e-6;

            % Default-disabled audit path for an infeasible hard QP.  It
            % relaxes only the linearized equalities for one restoration
            % step; final nonlinear acceptance remains entirely hard.
            opts.FeasibilityRestorationEnabled = false;
            opts.FeasibilityRestorationPenaltyL1 = 1e5;
            opts.FeasibilityRestorationPenaltyL2 = 1e-8;
            opts.FeasibilityRestorationPostIterations = 1;

            % Hessian approximation.
            %   "bfgs"         : damped BFGS Hessian of Lagrangian
            %   "lbfgs"        : limited-memory damped BFGS Hessian
            %   "identity"     : fixed identity/scaled identity Hessian
            %   "gauss_newton" : analytic residual-Jacobian curvature
            % opts.HessianMode         = "bfgs";
            opts.HessianMode         = "bfgs";
            opts.InitialHessianScale = 1.0;
            opts.HessianRegularization = 1e-9;
            opts.GaussNewtonDamping = 0;
            opts.ResetHessianOnFailure = true;
            opts.LimitedMemoryBfgsPairs = 5;

            % Gradient check.
            opts.CheckGradientsOnce  = true;
            opts.GradientCheckType   = "directional";
            opts.GradientCheckDirections = 12;
            opts.GradientCheckStep   = 1e-6;
            opts.GradientCheckVerbose = true;

            % Display.
            opts.Display = "none";     % "none" | "iter" | "final"

            % Optional solver-phase timing. Disabled execution does not
            % emit timing fields or change the numerical algorithm.
            opts.CollectTiming = false;

            % Disabled-by-default capture for bounded solver audits.
            opts.AuditCapture = false;
            opts.LiveIterationPlot = false;

            % Audit-only merit evaluation acceleration. When enabled, the
            % current-point merit reuses already evaluated NLP values, and
            % compatible constraint callbacks may skip trial Jacobians.
            opts.UseValueOnlyMeritEvaluation = false;

            % Default-inactive condensed-RTI QP selection. The direct path
            % is used only when its unconstrained reduced Newton step is
            % finite and satisfies every transformed constraint; otherwise
            % it falls back to a warm-started active-set QP and finally the
            % configured QP algorithm.
            opts.CondensedRtiQpMode = "configured_quadprog";

            % Audit-only line-search selection. The production default is
            % the established L1-merit path; "filter" enables the
            % Fletcher-Leyffer feasibility/objective acceptance audit.
            opts.LineSearchStrategy = "l1_merit";
            opts.FilterGammaTheta = 1e-5;
            opts.FilterGammaObjective = 1e-5;
            opts.FilterDelta = 1.0;
            opts.FilterSTheta = 1.1;
            opts.FilterSObjective = 2.3;
            opts.FilterEtaObjective = 1e-4;
            opts.FilterThetaMaximumFactor = 1e4;
            opts.FilterThetaMinimumFactor = 1e-4;

            % quadprog options.
            opts.QPOptions = optimoptions('quadprog', ...
                'Algorithm','interior-point-convex', ...
                'Display','off', ...
                'OptimalityTolerance',1e-8, ...
                'ConstraintTolerance',1e-8, ...
                'StepTolerance',1e-10, ...
                'MaxIterations',200);
        end
    end

    methods
        %--------------------------------------------------------------
        function obj = SQPSolver(opts)
            if nargin < 1 || isempty(opts)
                opts = AeroFlex.optim.SQPSolver.defaultOptions();
            else
                opts = obj.mergeOptions(AeroFlex.optim.SQPSolver.defaultOptions(),opts);
            end

            if ~(isscalar(opts.CollectTiming) && ...
                    (islogical(opts.CollectTiming) || ...
                     (isnumeric(opts.CollectTiming) && ...
                      isfinite(opts.CollectTiming) && ...
                      ismember(opts.CollectTiming,[0,1]))))
                error('SQPSolver:CollectTiming', ...
                    'CollectTiming must be a scalar logical or binary value.');
            end
            opts.CollectTiming = logical(opts.CollectTiming);

            if ~(isscalar(opts.LiveIterationPlot) && ...
                    (islogical(opts.LiveIterationPlot) || ...
                     (isnumeric(opts.LiveIterationPlot) && ...
                      isfinite(opts.LiveIterationPlot) && ...
                      ismember(opts.LiveIterationPlot,[0,1]))))
                error('SQPSolver:LiveIterationPlot', ...
                    'LiveIterationPlot must be a scalar logical or binary value.');
            end
            opts.LiveIterationPlot = logical(opts.LiveIterationPlot);

            if ~(isscalar(opts.UseValueOnlyMeritEvaluation) && ...
                    (islogical(opts.UseValueOnlyMeritEvaluation) || ...
                     (isnumeric(opts.UseValueOnlyMeritEvaluation) && ...
                      isfinite(opts.UseValueOnlyMeritEvaluation) && ...
                      ismember(opts.UseValueOnlyMeritEvaluation,[0,1]))))
                error('SQPSolver:UseValueOnlyMeritEvaluation', ...
                    ['UseValueOnlyMeritEvaluation must be a scalar ', ...
                     'logical or binary value.']);
            end
            opts.UseValueOnlyMeritEvaluation = ...
                logical(opts.UseValueOnlyMeritEvaluation);

            if ~(isscalar(opts.SecondOrderCorrectionEnabled) && ...
                    (islogical(opts.SecondOrderCorrectionEnabled) || ...
                     (isnumeric(opts.SecondOrderCorrectionEnabled) && ...
                      isfinite(opts.SecondOrderCorrectionEnabled) && ...
                      ismember(opts.SecondOrderCorrectionEnabled,[0,1]))))
                error('SQPSolver:SecondOrderCorrectionEnabled', ...
                    ['SecondOrderCorrectionEnabled must be a scalar ', ...
                     'logical or binary value.']);
            end
            opts.SecondOrderCorrectionEnabled = ...
                logical(opts.SecondOrderCorrectionEnabled);
            obj.validatePositiveOption(opts.SecondOrderCorrectionMaxRelativeNorm, ...
                'SecondOrderCorrectionMaxRelativeNorm');

            if ~(isscalar(opts.FeasibilityRestorationEnabled) && ...
                    (islogical(opts.FeasibilityRestorationEnabled) || ...
                     (isnumeric(opts.FeasibilityRestorationEnabled) && ...
                      isfinite(opts.FeasibilityRestorationEnabled) && ...
                      ismember(opts.FeasibilityRestorationEnabled,[0,1]))))
                error('SQPSolver:FeasibilityRestorationEnabled', ...
                    ['FeasibilityRestorationEnabled must be a scalar ', ...
                     'logical or binary value.']);
            end
            opts.FeasibilityRestorationEnabled = ...
                logical(opts.FeasibilityRestorationEnabled);
            obj.validatePositiveOption(opts.FeasibilityRestorationPenaltyL1, ...
                'FeasibilityRestorationPenaltyL1');
            obj.validatePositiveOption(opts.FeasibilityRestorationPenaltyL2, ...
                'FeasibilityRestorationPenaltyL2');
            if ~(isnumeric(opts.FeasibilityRestorationPostIterations) && ...
                    isscalar(opts.FeasibilityRestorationPostIterations) && ...
                    isfinite(opts.FeasibilityRestorationPostIterations) && ...
                    opts.FeasibilityRestorationPostIterations >= 0 && ...
                    opts.FeasibilityRestorationPostIterations == ...
                    floor(opts.FeasibilityRestorationPostIterations))
                error('SQPSolver:FeasibilityRestorationPostIterations', ...
                    ['FeasibilityRestorationPostIterations must be a ', ...
                     'nonnegative integer.']);
            end

            if ~(isstring(opts.LineSearchStrategy) || ...
                    ischar(opts.LineSearchStrategy)) || ...
                    ~isscalar(string(opts.LineSearchStrategy))
                error('SQPSolver:LineSearchStrategy', ...
                    'LineSearchStrategy must be a scalar string or character vector.');
            end
            opts.LineSearchStrategy = lower(string(opts.LineSearchStrategy));
            if ~ismember(opts.LineSearchStrategy, ...
                    ["l1_merit","filter","safeguarded_hybrid"])
                error('SQPSolver:LineSearchStrategy', ...
                ['LineSearchStrategy must be "l1_merit", "filter", ', ...
                     'or "safeguarded_hybrid".']);
            end

            obj.validateOpenUnitOption(opts.FilterGammaTheta, ...
                'FilterGammaTheta');
            obj.validateOpenUnitOption(opts.FilterGammaObjective, ...
                'FilterGammaObjective');
            obj.validatePositiveOption(opts.FilterDelta,'FilterDelta');
            obj.validatePositiveOption(opts.FilterSTheta,'FilterSTheta');
            obj.validatePositiveOption(opts.FilterSObjective, ...
                'FilterSObjective');
            obj.validateOpenUnitOption(opts.FilterEtaObjective, ...
                'FilterEtaObjective');
            obj.validatePositiveOption(opts.FilterThetaMaximumFactor, ...
                'FilterThetaMaximumFactor');
            obj.validatePositiveOption(opts.FilterThetaMinimumFactor, ...
                'FilterThetaMinimumFactor');
            if opts.FilterSTheta <= 1 || ...
                    opts.FilterSObjective <= 2*opts.FilterSTheta
                error('SQPSolver:FilterExponents', ...
                    ['FilterSTheta must exceed 1 and FilterSObjective ', ...
                     'must exceed 2*FilterSTheta.']);
            end
            if opts.FilterEtaObjective >= 0.5 || ...
                    opts.FilterThetaMaximumFactor <= 1 || ...
                    opts.FilterThetaMinimumFactor >= 1
                error('SQPSolver:FilterOptionRange', ...
                    ['FilterEtaObjective must be below 0.5, the maximum ', ...
                     'factor must exceed 1, and the minimum factor must ', ...
                     'be below 1.']);
            end
            if opts.LineSearchStrategy == "filter" && ...
                    opts.UseValueOnlyMeritEvaluation
                error('SQPSolver:FilterValueOnlyCombination', ...
                    ['The filter audit requires full-Jacobian trials. ', ...
                     'Set UseValueOnlyMeritEvaluation to false.']);
            end
            if opts.LineSearchStrategy == "safeguarded_hybrid" && ...
                    ~opts.UseValueOnlyMeritEvaluation
                error('SQPSolver:HybridValueOnlyRequirement', ...
                    ['The safeguarded hybrid audit requires the qualified ', ...
                     'value-only merit evaluation path.']);
            end
            if ~(isstring(opts.LineSearchUpdate) || ...
                    ischar(opts.LineSearchUpdate)) || ...
                    ~isscalar(string(opts.LineSearchUpdate))
                error('SQPSolver:LineSearchUpdate', ...
                    'LineSearchUpdate must be a scalar string or character vector.');
            end
            opts.LineSearchUpdate = lower(string(opts.LineSearchUpdate));
            if ~ismember(opts.LineSearchUpdate,["fixed_beta" "quadratic_safeguarded"])
                error('SQPSolver:LineSearchUpdate', ...
                    ['LineSearchUpdate must be "fixed_beta" or ', ...
                     '"quadratic_safeguarded".']);
            end

            if ~(isstring(opts.HessianMode) || ischar(opts.HessianMode)) || ...
                    ~isscalar(string(opts.HessianMode))
                error('SQPSolver:HessianMode', ...
                    'HessianMode must be a scalar string or character vector.');
            end
            opts.HessianMode = lower(string(opts.HessianMode));
            if ~ismember(opts.HessianMode, ...
                    ["bfgs","lbfgs","identity","gauss_newton"])
                error('SQPSolver:HessianMode', ...
                    ['HessianMode must be "bfgs", "lbfgs", "identity", ', ...
                     'or "gauss_newton".']);
            end
            if ~(isnumeric(opts.LimitedMemoryBfgsPairs) && ...
                    isscalar(opts.LimitedMemoryBfgsPairs) && ...
                    isfinite(opts.LimitedMemoryBfgsPairs) && ...
                    opts.LimitedMemoryBfgsPairs >= 1 && ...
                    opts.LimitedMemoryBfgsPairs == ...
                    floor(opts.LimitedMemoryBfgsPairs))
                error('SQPSolver:LimitedMemoryBfgsPairs', ...
                    'LimitedMemoryBfgsPairs must be a positive integer.');
            end
            if ~(isnumeric(opts.GaussNewtonDamping) && ...
                    isscalar(opts.GaussNewtonDamping) && ...
                    isfinite(opts.GaussNewtonDamping) && ...
                    opts.GaussNewtonDamping >= 0)
                error('SQPSolver:GaussNewtonDamping', ...
                    'GaussNewtonDamping must be a finite nonnegative scalar.');
            end

            obj.options = opts;
            obj.H = [];
            obj.lastOutput = struct();
            obj.lastLambda = struct();
            obj.lastGradCheck = struct();
            obj.lastPreparedInitial = struct();
            obj.lbfgsSHistory = cell(0,1);
            obj.lbfgsYHistory = cell(0,1);
        end

        %--------------------------------------------------------------
        function prepared = prepareInitialConstraints(obj,z,lb,ub,nonlFun)
        %PREPAREINITIALCONSTRAINTS Cache one exact initial constraint evaluation.
        % The cache is valid only for the identical decision and simple bounds.
            z = z(:);
            n = numel(z);
            if isempty(lb), lb = -inf(n,1); else, lb = lb(:); end
            if isempty(ub), ub = inf(n,1); else, ub = ub(:); end
            assert(numel(lb) == n && numel(ub) == n, ...
                'SQPSolver:PreparedBounds', ...
                'Prepared constraint bounds must match the decision length.');
            timer = tic;
            [c,ceq,gradc,gradceq] = nonlFun(z);
            seconds = toc(timer);
            c = c(:);
            ceq = ceq(:);
            if isempty(gradc), gradc = zeros(n,0); end
            if isempty(gradceq), gradceq = zeros(n,0); end
            if size(gradc,1) ~= n && size(gradc,2) == n, gradc = gradc.'; end
            if size(gradceq,1) ~= n && size(gradceq,2) == n
                gradceq = gradceq.';
            end

            assert(size(gradc,1) == n && size(gradc,2) == numel(c), ...
                'SQPSolver:PreparedGradient', ...
                'Prepared inequality Jacobian has inconsistent dimensions.');
            assert(size(gradceq,1) == n && size(gradceq,2) == numel(ceq), ...
                'SQPSolver:PreparedGradient', ...
                'Prepared equality Jacobian has inconsistent dimensions.');
            assert(all(isfinite(c)) && all(isfinite(ceq)) && ...
                all(isfinite(gradc),'all') && all(isfinite(gradceq),'all'), ...
                'SQPSolver:PreparedFinite', ...
                'Prepared constraints and Jacobians must be finite.');
            prepared = struct('valid',true,'z',z,'lb',lb,'ub',ub, ...
                'c',c,'ceq',ceq,'gradc',gradc,'gradceq',gradceq, ...
                'seconds',seconds);
            obj.lastPreparedInitial = prepared;
        end

        %--------------------------------------------------------------
        function [z,output] = solveCondensedRti( ...
                obj,costFun,z0,lb,ub,nonlFun,valueNonlFun,request)
        %SOLVECONDENSEDRTI Apply bounded relinearized condensed QP steps.
        % This audit-only method eliminates selected equality-constrained
        % coordinates exactly at each linearization. It does not perform a
        % line search or relax the final nonlinear feasibility checks.
            required = {'auditOnly','transformation','eliminatedColumns', ...
                'eliminatedRows','iterationCount','constraintTolerance', ...
                'inequalityTolerance'};
            assert(isstruct(request) && isscalar(request) && ...
                all(isfield(request,required)) && logical(request.auditOnly), ...
                'SQPSolver:CondensedRtiRequest', ...
                'The condensed RTI request is incomplete or not audit-only.');

            z = z0(:);
            lb = lb(:);
            ub = ub(:);
            n = numel(z);
            assert(numel(lb)==n && numel(ub)==n && all(lb<=ub), ...
                'SQPSolver:CondensedRtiBounds', ...
                'The condensed RTI bounds are invalid.');
            transformation = sparse(request.transformation);
            assert(size(transformation,1)==n && ...
                all(isfinite(nonzeros(transformation))), ...
                'SQPSolver:CondensedRtiTransformation', ...
                'The condensed RTI transformation is invalid.');
            eliminatedColumns = double(request.eliminatedColumns(:).');
            eliminatedRows = double(request.eliminatedRows(:).');
            assert(all(isfinite(eliminatedColumns)) && ...
                all(eliminatedColumns==fix(eliminatedColumns)) && ...
                all(eliminatedColumns>=1) && ...
                all(eliminatedColumns<=size(transformation,2)) && ...
                numel(unique(eliminatedColumns))==numel(eliminatedColumns), ...
                'SQPSolver:CondensedRtiElimination', ...
                'The eliminated RTI coordinate indices are invalid.');
            assert(isscalar(request.iterationCount) && ...
                isfinite(request.iterationCount) && ...
                request.iterationCount==fix(request.iterationCount) && ...
                request.iterationCount>=1 && request.iterationCount<=3, ...
                'SQPSolver:CondensedRtiIterationCount', ...
                'The condensed RTI iteration count must be an integer in [1,3].');
            assert(isscalar(request.constraintTolerance) && ...
                isfinite(request.constraintTolerance) && ...
                request.constraintTolerance>0 && ...
                isscalar(request.inequalityTolerance) && ...
                isfinite(request.inequalityTolerance) && ...
                request.inequalityTolerance>=0, ...
                'SQPSolver:CondensedRtiTolerance', ...
                'The condensed RTI tolerances are invalid.');
            backtracking = struct('enabled',false,'maxTrials',0,'beta',nan);
            if isfield(request,'nonlinearBacktracking')
                backtracking = request.nonlinearBacktracking;
                assert(isstruct(backtracking) && isscalar(backtracking) && ...
                    isfield(backtracking,'enabled') && ...
                    isscalar(backtracking.enabled), ...
                    'SQPSolver:CondensedRtiBacktracking', ...
                    'The nonlinear-backtracking request is invalid.');
                if logical(backtracking.enabled)
                    assert(isfield(backtracking,'maxTrials') && ...
                        isfield(backtracking,'beta') && ...
                        isscalar(backtracking.maxTrials) && ...
                        backtracking.maxTrials==fix(backtracking.maxTrials) && ...
                        backtracking.maxTrials>=1 && backtracking.maxTrials<=20 && ...
                        isscalar(backtracking.beta) && ...
                        isfinite(backtracking.beta) && backtracking.beta>0 && ...
                        backtracking.beta<1, ...
                        'SQPSolver:CondensedRtiBacktracking', ...
                        'The nonlinear-backtracking parameters are invalid.');
                end
            end

            iterationTemplate = struct('iteration',0, ...
                'derivativeSeconds',0,'condensingSeconds',0, ...
                'objectiveDerivativeSeconds',0, ...
                'constraintDerivativeSeconds',0, ...
                'qpSeconds',0,'qpExitflag',nan,'qpIterations',nan, ...
                'qpFirstOrderOptimality',nan,'reducedKktInf',nan, ...
                'qpComplementarity',nan,'qpMinimumLinearSlack',nan, ...
                'qpMaximumInactiveMultiplier',nan, ...
                'qpMethod',"", ...
                'decisionCount',0,'equalityCount',0, ...
                'inequalityCount',0,'linearEqualityResidualInf',inf, ...
                'linearInequalityResidualMax',inf,'stepInf',inf);
            iterations = repmat(iterationTemplate,request.iterationCount,1);
            totalTimer = tic;
            allQpAccepted = true;
            message = "Condensed RTI completed.";
            for iteration = 1:request.iterationCount
                derivativeTimer = tic;
                objectiveDerivativeTimer = tic;
                [~,gradient,hessian] = costFun(z);
                objectiveDerivativeSeconds = toc(objectiveDerivativeTimer);
                constraintDerivativeTimer = tic;
                [inequality,equality,gradientInequality,gradientEquality] = ...
                    nonlFun(z);
                constraintDerivativeSeconds = toc(constraintDerivativeTimer);
                inequality = inequality(:);
                equality = equality(:);
                if isempty(inequality)
                    gradientInequality = sparse(n,0);
                elseif size(gradientInequality,1)~=n && ...
                        size(gradientInequality,2)==n
                    gradientInequality = gradientInequality.';
                end
                if isempty(equality)
                    gradientEquality = sparse(n,0);
                elseif size(gradientEquality,1)~=n && ...
                        size(gradientEquality,2)==n
                    gradientEquality = gradientEquality.';
                end
                derivativeSeconds = toc(derivativeTimer);
                audit = struct('gradient',gradient(:),'hessian',hessian, ...
                    'inequality',inequality,'equality',equality, ...
                    'gradientInequality',gradientInequality, ...
                    'gradientEquality',gradientEquality, ...
                    'lowerBound',lb,'upperBound',ub);
                [step,stepInfo] = obj.solveSelectedCondensedQp( ...
                    audit,z,transformation,eliminatedColumns, ...
                    eliminatedRows);
                iterations(iteration) = stepInfo;
                iterations(iteration).iteration = iteration;
                iterations(iteration).derivativeSeconds = derivativeSeconds;
                iterations(iteration).objectiveDerivativeSeconds = ...
                    objectiveDerivativeSeconds;
                iterations(iteration).constraintDerivativeSeconds = ...
                    constraintDerivativeSeconds;
                if stepInfo.qpExitflag<=0 || isempty(step) || ...
                        any(~isfinite(step))
                    allQpAccepted = false;
                    message = "A condensed RTI QP failed; the seed was rejected.";
                    break
                end
                if logical(backtracking.enabled)
                    acceptedTrial = false;
                    alpha = 1;
                    for trialIndex = 1:double(backtracking.maxTrials)
                        trial = min(max(z+alpha*step,lb),ub);
                        [trialInequality,trialEquality] = valueNonlFun(trial);
                        trialViolation = max([norm(trialEquality(:),inf); ...
                            max([trialInequality(:);0]); ...
                            max([lb-trial;trial-ub;0])]);
                        if isfinite(trialViolation) && ...
                                trialViolation<=request.constraintTolerance
                            z = trial;
                            acceptedTrial = true;
                            break
                        end
                        alpha = alpha*backtracking.beta;
                    end
                    if ~acceptedTrial
                        allQpAccepted = false;
                        message = "Condensed RTI nonlinear backtracking rejected the step.";
                        break
                    end
                else
                    z = z+step;
                end
            end

            valueTimer = tic;
            [finalInequality,finalEquality] = valueNonlFun(z);
            valueReplaySeconds = toc(valueTimer);
            boundViolation = max([lb-z;z-ub;0]);
            equalityInf = norm(finalEquality(:),inf);
            inequalityMax = max([finalInequality(:);0]);
            finalViolation = max([equalityInf;inequalityMax;boundViolation]);
            qualified = allQpAccepted && all(isfinite(z)) && ...
                equalityInf<=request.constraintTolerance && ...
                inequalityMax<=request.inequalityTolerance && ...
                boundViolation<=request.inequalityTolerance;
            if ~qualified && allQpAccepted
                message = "Condensed RTI nonlinear replay exceeded tolerance.";
            end
            [finalObjective,~] = costFun(z);
            output = struct('qualified',qualified,'message',message, ...
                'iterationCount',request.iterationCount, ...
                'completedIterations',sum([iterations.qpExitflag]>0), ...
                'iterations',iterations, ...
                'nonlinearEqualityInf',equalityInf, ...
                'nonlinearInequalityMax',inequalityMax, ...
                'boundViolationMax',boundViolation, ...
                'constraintViolationInf',finalViolation, ...
                'objective',finalObjective, ...
                'valueReplaySeconds',valueReplaySeconds, ...
                'totalSeconds',toc(totalTimer));
        end

        %--------------------------------------------------------------
        function [z,output] = solvePreparedReducedRti(obj,costFun,z0, ...
                lb,ub,valueNonlFun,reducedPacketFun,request)
        %SOLVEPREPAREDREDUCEDRTI Solve an exactly prepared reduced RTI QP.
        %   The caller supplies the affine full-step expansion produced by
        %   causal horizon tangent recursion. The existing QP backend and
        %   exact nonlinear replay remain the acceptance owners.
            required = {'auditOnly','changeId','iterationCount', ...
                'constraintTolerance','inequalityTolerance'};
            assert(isstruct(request) && isscalar(request) && ...
                all(isfield(request,required)) && logical(request.auditOnly) && ...
                string(request.changeId)== ...
                    "phase18c-v17a-casebc-native-reduced-horizon-rti-audit-v1", ...
                'SQPSolver:PreparedReducedRtiRequest', ...
                'The prepared reduced RTI request is unauthorized.');
            assert(isscalar(request.iterationCount) && ...
                request.iterationCount==fix(request.iterationCount) && ...
                request.iterationCount>=1 && request.iterationCount<=3 && ...
                isscalar(request.constraintTolerance) && ...
                isfinite(request.constraintTolerance) && ...
                request.constraintTolerance>0 && ...
                isscalar(request.inequalityTolerance) && ...
                isfinite(request.inequalityTolerance) && ...
                request.inequalityTolerance>=0, ...
                'SQPSolver:PreparedReducedRtiOptions', ...
                'The prepared reduced RTI options are invalid.');

            z = z0(:);
            lb = lb(:);
            ub = ub(:);
            n = numel(z);
            assert(numel(lb)==n && numel(ub)==n && all(lb<=ub), ...
                'SQPSolver:PreparedReducedRtiBounds', ...
                'The prepared reduced RTI bounds are invalid.');
            backtracking = struct('enabled',false,'maxTrials',0,'beta',nan);
            if isfield(request,'nonlinearBacktracking')
                backtracking = request.nonlinearBacktracking;
            end
            if logical(backtracking.enabled)
                assert(isfield(backtracking,'maxTrials') && ...
                    isfield(backtracking,'beta') && ...
                    isscalar(backtracking.maxTrials) && ...
                    backtracking.maxTrials==fix(backtracking.maxTrials) && ...
                    backtracking.maxTrials>=1 && ...
                    backtracking.maxTrials<=20 && ...
                    isscalar(backtracking.beta) && ...
                    isfinite(backtracking.beta) && ...
                    backtracking.beta>0 && backtracking.beta<1, ...
                    'SQPSolver:PreparedReducedRtiBacktracking', ...
                    'The prepared reduced RTI backtracking request is invalid.');
            end

            iterationTemplate = struct('iteration',0, ...
                'derivativeSeconds',0,'condensingSeconds',0, ...
                'objectiveDerivativeSeconds',0, ...
                'constraintDerivativeSeconds',0, ...
                'qpSeconds',0,'qpExitflag',nan,'qpIterations',nan, ...
                'qpFirstOrderOptimality',nan,'reducedKktInf',nan, ...
                'qpComplementarity',nan,'qpMinimumLinearSlack',nan, ...
                'qpMaximumInactiveMultiplier',nan,'qpMethod',"", ...
                'decisionCount',0,'equalityCount',0, ...
                'inequalityCount',0,'linearEqualityResidualInf',inf, ...
                'linearInequalityResidualMax',inf,'stepInf',inf);
            iterations = repmat(iterationTemplate,request.iterationCount,1);
            totalTimer = tic;
            allQpAccepted = true;
            message = "Prepared reduced RTI completed.";
            for iteration = 1:request.iterationCount
                derivativeTimer = tic;
                objectiveTimer = tic;
                [~,gradient,hessian] = costFun(z);
                objectiveSeconds = toc(objectiveTimer);
                packetTimer = tic;
                packet = reducedPacketFun(z);
                packetSeconds = toc(packetTimer);
                requiredPacket = {'expansion','offset','A','b','Aeq','beq'};
                assert(isstruct(packet) && isscalar(packet) && ...
                    all(isfield(packet,requiredPacket)), ...
                    'SQPSolver:PreparedReducedRtiPacket', ...
                    'The prepared reduced RTI packet is incomplete.');
                expansion = sparse(packet.expansion);
                offset = packet.offset(:);
                reducedCount = size(expansion,2);
                assert(size(expansion,1)==n && numel(offset)==n && ...
                    all(isfinite(nonzeros(expansion))) && ...
                    all(isfinite(offset)) && ...
                    isequal(size(hessian),[n n]) && ...
                    numel(gradient)==n && all(isfinite(gradient)) && ...
                    all(isfinite(hessian),'all'), ...
                    'SQPSolver:PreparedReducedRtiPacketFinite', ...
                    'The prepared reduced RTI packet is invalid or nonfinite.');
                preparationTimer = tic;
                reducedHessian = sparse(expansion.'*hessian*expansion);
                reducedHessian = 0.5*(reducedHessian+reducedHessian.');
                reducedGradient = expansion.'*(hessian*offset+gradient(:));
                A = sparse(packet.A);
                b = packet.b(:);
                Aeq = sparse(packet.Aeq);
                beq = packet.beq(:);
                assert(size(A,2)==reducedCount && numel(b)==size(A,1) && ...
                    size(Aeq,2)==reducedCount && ...
                    numel(beq)==size(Aeq,1) && ...
                    all(isfinite(nonzeros(A))) && all(isfinite(b)) && ...
                    all(isfinite(nonzeros(Aeq))) && all(isfinite(beq)), ...
                    'SQPSolver:PreparedReducedRtiLinearConstraints', ...
                    'The prepared reduced constraints are invalid.');
                upperStep = ub-z;
                lowerStep = z-lb;
                finiteUpper = isfinite(upperStep);
                finiteLower = isfinite(lowerStep);
                packetInequalityCount = size(A,1);
                upperCount = nnz(finiteUpper);
                lowerCount = nnz(finiteLower);
                boundedA = spalloc( ...
                    packetInequalityCount+upperCount+lowerCount, ...
                    reducedCount,nnz(A)+nnz(expansion(finiteUpper,:))+ ...
                    nnz(expansion(finiteLower,:)));
                boundedB = zeros( ...
                    packetInequalityCount+upperCount+lowerCount,1);
                packetRows = 1:packetInequalityCount;
                upperRows = packetInequalityCount+(1:upperCount);
                lowerRows = packetInequalityCount+upperCount+(1:lowerCount);
                boundedA(packetRows,:) = A;
                boundedB(packetRows) = b;
                boundedA(upperRows,:) = expansion(finiteUpper,:);
                boundedB(upperRows) = ...
                    upperStep(finiteUpper)-offset(finiteUpper);
                boundedA(lowerRows,:) = -expansion(finiteLower,:);
                boundedB(lowerRows) = ...
                    lowerStep(finiteLower)+offset(finiteLower);
                A = boundedA;
                b = boundedB;
                condensingSeconds = toc(preparationTimer);
                derivativeSeconds = toc(derivativeTimer);

                qpTimer = tic;
                [reducedStep,qpExitflag,qpOutput,qpLambda,qpMethod] = ...
                    obj.solveCondensedReducedQp( ...
                        reducedHessian,reducedGradient,A,b,Aeq,beq);
                qpSeconds = toc(qpTimer);
                stepInfo = iterationTemplate;
                stepInfo.iteration = iteration;
                stepInfo.derivativeSeconds = derivativeSeconds;
                stepInfo.objectiveDerivativeSeconds = objectiveSeconds;
                stepInfo.constraintDerivativeSeconds = packetSeconds;
                stepInfo.condensingSeconds = condensingSeconds;
                stepInfo.qpSeconds = qpSeconds;
                stepInfo.qpExitflag = qpExitflag;
                stepInfo.qpMethod = qpMethod;
                stepInfo.decisionCount = reducedCount;
                stepInfo.equalityCount = size(Aeq,1);
                stepInfo.inequalityCount = size(A,1);
                if isstruct(qpOutput) && isfield(qpOutput,'iterations') && ...
                        isscalar(qpOutput.iterations)
                    stepInfo.qpIterations = double(qpOutput.iterations);
                end
                if isstruct(qpOutput) && ...
                        isfield(qpOutput,'firstorderopt') && ...
                        isscalar(qpOutput.firstorderopt)
                    stepInfo.qpFirstOrderOptimality = ...
                        double(qpOutput.firstorderopt);
                end
                if isempty(reducedStep) || any(~isfinite(reducedStep)) || ...
                        qpExitflag<=0
                    iterations(iteration) = stepInfo;
                    allQpAccepted = false;
                    message = "A prepared reduced RTI QP failed.";
                    break
                end
                fullStep = offset+expansion*reducedStep;
                stepInfo.linearEqualityResidualInf = ...
                    norm(Aeq*reducedStep-beq,inf);
                stepInfo.linearInequalityResidualMax = ...
                    max([A*reducedStep-b;0]);
                stepInfo.stepInf = norm(fullStep,inf);
                reducedKkt = reducedHessian*reducedStep+reducedGradient;
                if isstruct(qpLambda)
                    if isfield(qpLambda,'ineqlin') && ...
                            numel(qpLambda.ineqlin)==size(A,1)
                        reducedKkt = reducedKkt+A.'*qpLambda.ineqlin(:);
                        slack = b-A*reducedStep;
                        multiplier = qpLambda.ineqlin(:);
                        stepInfo.qpComplementarity = ...
                            max(abs(slack.*multiplier));
                        stepInfo.qpMinimumLinearSlack = min(slack);
                        inactive = slack>1e-7;
                        if any(inactive)
                            stepInfo.qpMaximumInactiveMultiplier = ...
                                max(abs(multiplier(inactive)));
                        else
                            stepInfo.qpMaximumInactiveMultiplier = 0;
                        end
                    end
                    if isfield(qpLambda,'eqlin') && ...
                            numel(qpLambda.eqlin)==size(Aeq,1)
                        reducedKkt = reducedKkt+Aeq.'*qpLambda.eqlin(:);
                    end
                end
                stepInfo.reducedKktInf = norm(reducedKkt,inf);
                iterations(iteration) = stepInfo;

                if logical(backtracking.enabled)
                    acceptedTrial = false;
                    alpha = 1;
                    for trialIndex = 1:double(backtracking.maxTrials)
                        trial = min(max(z+alpha*fullStep,lb),ub);
                        [trialInequality,trialEquality] = valueNonlFun(trial);
                        trialViolation = max([norm(trialEquality(:),inf); ...
                            max([trialInequality(:);0]); ...
                            max([lb-trial;trial-ub;0])]);
                        if isfinite(trialViolation) && ...
                                trialViolation<=request.constraintTolerance
                            z = trial;
                            acceptedTrial = true;
                            break
                        end
                        alpha = alpha*backtracking.beta;
                    end
                    if ~acceptedTrial
                        allQpAccepted = false;
                        message = "Prepared reduced RTI backtracking rejected the step.";
                        break
                    end
                else
                    z = z+fullStep;
                end
            end

            valueTimer = tic;
            [finalInequality,finalEquality] = valueNonlFun(z);
            valueReplaySeconds = toc(valueTimer);
            boundViolation = max([lb-z;z-ub;0]);
            equalityInf = norm(finalEquality(:),inf);
            inequalityMax = max([finalInequality(:);0]);
            finalViolation = max([equalityInf;inequalityMax;boundViolation]);
            qualified = allQpAccepted && all(isfinite(z)) && ...
                equalityInf<=request.constraintTolerance && ...
                inequalityMax<=request.inequalityTolerance && ...
                boundViolation<=request.inequalityTolerance;
            if ~qualified && allQpAccepted
                message = "Prepared reduced RTI nonlinear replay exceeded tolerance.";
            end
            [finalObjective,~] = costFun(z);
            output = struct('qualified',qualified,'message',message, ...
                'iterationCount',request.iterationCount, ...
                'completedIterations',sum([iterations.qpExitflag]>0), ...
                'iterations',iterations, ...
                'nonlinearEqualityInf',equalityInf, ...
                'nonlinearInequalityMax',inequalityMax, ...
                'boundViolationMax',boundViolation, ...
                'constraintViolationInf',finalViolation, ...
                'objective',finalObjective, ...
                'valueReplaySeconds',valueReplaySeconds, ...
                'totalSeconds',toc(totalTimer));
        end

        %--------------------------------------------------------------
        function resetHessian(obj,nVar)
            obj.H = obj.options.InitialHessianScale * speye(nVar);
            obj.lbfgsSHistory = cell(0,1);
            obj.lbfgsYHistory = cell(0,1);
        end

        %--------------------------------------------------------------
        function [z,fval,exitflag,output,lambda] = solve(obj,costFun,z0,lb,ub,nonlFun,varargin)
        %SOLVE Run SQP.
        %
        % Inputs:
        %   costFun : function handle [f,g] = costFun(z)
        %   z0      : initial decision vector
        %   lb, ub  : simple bounds
        %   nonlFun : function handle [c,ceq,gradc,gradceq] = nonlFun(z)
        %
        % Outputs match fmincon-style enough for nMPC diagnostics.

            opts = obj.options;
            obj.lastAuditCapture = struct();
            livePlot = obj.createLiveIterationPlot(opts.LiveIterationPlot);
            livePlotCleanup = onCleanup(@() obj.closeLiveIterationPlot(livePlot));

            preparedInitial = struct();
            if ~isempty(varargin)
            assert(isscalar(varargin) && isstruct(varargin{1}), ...
                    'SQPSolver:PreparedInput', ...
                    'The optional prepared input must be one scalar struct.');
                preparedInitial = varargin{1};
            end

            z = z0(:);
            n = numel(z);

            lb = lb(:);
            ub = ub(:);

            if isempty(lb), lb = -inf(n,1); end
            if isempty(ub), ub =  inf(n,1); end

            assert(numel(lb)==n,'SQPSolver: lb length mismatch.');
            assert(numel(ub)==n,'SQPSolver: ub length mismatch.');

            % Project initial guess into simple bounds.
            z = min(max(z,lb),ub);

            if isempty(obj.H) || size(obj.H,1) ~= n || size(obj.H,2) ~= n
                obj.resetHessian(n);
            end

            funcCount = 0;
            constrCount = 0;

            exitflag = 0;
            message = "Maximum SQP iterations reached.";

            lambda = struct('eqnonlin',[],'ineqnonlin',[],'lower',[],'upper',[]);
            qpExitflag = nan;
            qpOutput = struct();

            alpha = 1;
            stepInf = inf;
            slackEqInf = 0;
            slackIneqInf = 0;

            timing = struct();
            totalTimer = [];
            if opts.CollectTiming
                timing = obj.emptyTiming(opts.MaxIterations + ...
                    opts.FeasibilityRestorationPostIterations);
                totalTimer = tic;
            end

            preparedInitialUsed = obj.preparedInitialMatches( ...
                preparedInitial,z,lb,ub);
            try
                if preparedInitialUsed
                    objectiveTimer = [];
                    if opts.CollectTiming, objectiveTimer = tic; end
                    [f,g] = costFun(z);
                    evalTiming = struct( ...
                        'objectiveEvaluationSeconds',obj.elapsed(objectiveTimer,opts.CollectTiming), ...
                        'constraintEvaluationSeconds',0);
                    g = g(:);
                    c = preparedInitial.c;
                    ceq = preparedInitial.ceq;
                    gradc = preparedInitial.gradc;
                    gradceq = preparedInitial.gradceq;
                else
                    [f,g,c,ceq,gradc,gradceq,evalTiming] = ...
                        obj.evalAll(costFun,nonlFun,z);
                end
                timing = obj.addEvaluationTiming(timing,evalTiming);
                funcCount = funcCount + 1;
                constrCount = constrCount + 1;
            catch ME
                fval = nan;
                exitflag = -3;
                output = obj.makeOutput(0,funcCount,constrCount,nan,nan,nan,nan, ...
                    "Initial function evaluation failed: " + string(ME.message), ...
                    nan,struct(),nan,nan);
                output = obj.attachTiming(output,timing,totalTimer,0);
                output.preparedInitialUsed = preparedInitialUsed;
                obj.lastOutput = output;
                return
            end

            feasInf = obj.constraintViolation(c,ceq,z,lb,ub);

            filterEnabled = opts.LineSearchStrategy == "filter";
            hybridEnabled = ...
                opts.LineSearchStrategy == "safeguarded_hybrid";
            filterTheta = nan(opts.MaxIterations,1);
            filterObjective = nan(opts.MaxIterations,1);
            filterCount = 0;
            thetaCurrent = obj.l1InfeasibilityFromValues(c,ceq,z,lb,ub);
            thetaScale = max(1,thetaCurrent);
            thetaMaximum = opts.FilterThetaMaximumFactor*thetaScale;
            thetaMinimum = opts.FilterThetaMinimumFactor*thetaScale;

            % Use a simple first-order check initially.
            % For equality-constrained problems, this is conservative because it ignores
            % the possibility that Lagrange multipliers balance the objective gradient.
            kktInf = norm(g,inf);
            
            % --------------------------------------------------------------
            % Initial-point exit check.
            % This is essential for warm-started real-time NMPC.
            % If the repaired warm start is already feasible and nearly stationary,
            % do not waste time solving a QP.
            % --------------------------------------------------------------
            if feasInf <= opts.ConstraintTolerance && ...
               kktInf <= opts.OptimalityTolerance
            
                fval = f;
                exitflag = 1;
                message = "SQP converged at initial point: feasibility and optimality satisfied.";
            
                output = obj.makeOutput(0,funcCount,constrCount,fval,feasInf,kktInf, ...
                    0,message,nan,struct(),1,0,0);
                output = obj.attachTiming(output,timing,totalTimer,0);
                output.preparedInitialUsed = preparedInitialUsed;
            
                lambda = struct('eqnonlin',zeros(numel(ceq),1), ...
                                'ineqnonlin',zeros(numel(c),1), ...
                                'lower',[], ...
                                'upper',[]);
            
                obj.lastOutput = output;
                obj.lastLambda = lambda;
            
                if opts.Display == "final" || opts.Display == "iter"
                    fprintf('[SQP] exitflag=%d | iter=%d | f=%.6e | feas=%.3e | opt=%.3e | step=%.3e | %s\n', ...
                        exitflag,0,fval,feasInf,kktInf,0,message);
                end
            
                return
            end

            if opts.Display == "iter"
                obj.printHeader();
                obj.printIter(0,f,feasInf,kktInf,0,1,0,0,nan);
            end
            obj.updateLiveIterationPlot(livePlot,0,f,feasInf,kktInf,0,1);

            zPrev = z;
            gPrev = g;
            gradcPrev = gradc;
            gradceqPrev = gradceq;
            lambdaEqPrev = zeros(numel(ceq),1);
            lambdaIneqPrev = zeros(numel(c),1);

            restorationUsed = false;
            for iter = 1:(opts.MaxIterations + ...
                    opts.FeasibilityRestorationPostIterations)
                if iter > opts.MaxIterations && ~restorationUsed
                    break
                end

                if opts.CollectTiming
                    timing.iteration(iter).iteration = iter;
                    timing.iteration(iter).objective = f;
                    timing.iteration(iter).constraintViolationInf = feasInf;
                    timing.iteration(iter).firstOrderOptimalityInf = kktInf;
                end

                % -----------------------------------------------------
                % 1. Build and solve QP subproblem
                % -----------------------------------------------------
                hessianTimer = [];
                if opts.CollectTiming
                    hessianTimer = tic;
                end
                if opts.HessianMode == "gauss_newton"
                    Hraw = obj.objectiveGaussNewtonHessian(costFun,z,n);
                    Hq = obj.makePositiveDefinite( ...
                        obj.dampGaussNewtonHessian(Hraw));
                elseif opts.HessianMode == "lbfgs"
                    Hraw = [];
                    Hq = obj.makePositiveDefinite( ...
                        obj.limitedMemoryBfgsHessian(n));
                else
                    Hraw = [];
                    Hq = obj.makePositiveDefinite(obj.H);
                end
                hessianSeconds = obj.elapsed(hessianTimer,opts.CollectTiming);
                timing = obj.addPhaseTiming( ...
                    timing,iter,"hessianWorkSeconds",hessianSeconds);

                qpAssemblyTimer = [];
                if opts.CollectTiming
                    qpAssemblyTimer = tic;
                end
                [Haug,faug,Aineq,bineq,Aeq,beq,lbq,ubq,slackInfo] = ...
                    obj.buildQP(Hq,g,c,ceq,gradc,gradceq,z,lb,ub);
                if opts.AuditCapture
                    if opts.HessianMode == "gauss_newton"
                        gaussNewtonAudit = obj.checkGaussNewtonHessian( ...
                            costFun,z,lb,ub,Hraw);
                    else
                        gaussNewtonAudit = struct();
                    end
                    obj.lastAuditCapture = struct( ...
                        'iteration',iter,'seed',z,'objective',f, ...
                        'gradientInfinityNorm',norm(g,inf), ...
                        'nonlinearInequalityInfinityNorm',norm(c,inf), ...
                        'nonlinearEqualityInfinityNorm',norm(ceq,inf), ...
                        'hessianMinimumEigenvalue',min(eig(full(Hq))), ...
                        'hessianMaximumEigenvalue',max(eig(full(Hq))), ...
                        'gaussNewtonDamping',opts.GaussNewtonDamping, ...
                        'gaussNewtonHessianAudit',gaussNewtonAudit, ...
                        'qpVariableCount',numel(faug), ...
                        'qpInequalityCount',size(Aineq,1), ...
                        'qpEqualityCount',size(Aeq,1), ...
                        'qpBoundFeasible',all(lbq <= ubq), ...
                        'qpAineq',Aineq,'qpBineq',bineq, ...
                        'qpAeq',Aeq,'qpBeq',beq, ...
                        'qpLb',lbq,'qpUb',ubq,'qpHessian',Haug, ...
                        'qpGradient',faug);
                end
                qpAssemblySeconds = obj.elapsed( ...
                    qpAssemblyTimer,opts.CollectTiming);
                timing = obj.addPhaseTiming( ...
                    timing,iter,"qpAssemblySeconds",qpAssemblySeconds);

                qpSolveTimer = [];
                if opts.CollectTiming
                    qpSolveTimer = tic;
                end
                try
                    [dAug,~,qpExitflag,qpOutput,qpLambda] = quadprog( ...
                        Haug,faug,Aineq,bineq,Aeq,beq,lbq,ubq,[],opts.QPOptions);
                catch ME
                    qpSolveSeconds = obj.elapsed( ...
                        qpSolveTimer,opts.CollectTiming);
                    timing = obj.addPhaseTiming( ...
                        timing,iter,"qpSolveSeconds",qpSolveSeconds);
                    exitflag = -1;
                    message = "quadprog call failed: " + string(ME.message);
                    break
                end
                qpSolveSeconds = obj.elapsed(qpSolveTimer,opts.CollectTiming);
                timing = obj.addPhaseTiming( ...
                    timing,iter,"qpSolveSeconds",qpSolveSeconds);
                if opts.CollectTiming
                    timing.iteration(iter).qpExitflag = qpExitflag;
                end

                restoration = struct('attempted',false,'accepted',false, ...
                    'exitflag',nan,'output',struct(),'slackEqInf',nan);
                if isempty(dAug) || qpExitflag <= 0
                    if opts.AuditCapture
                        obj.lastAuditCapture.qpExitflag = qpExitflag;
                        obj.lastAuditCapture.qpOutput = qpOutput;
                    end
                    if opts.FeasibilityRestorationEnabled
                        [dAug,qpExitflag,qpOutput,slackInfo,restoration] = ...
                            obj.solveFeasibilityRestoration( ...
                            c,ceq,gradc,gradceq,z,lb,ub,opts.QPOptions);
                        if opts.AuditCapture
                            obj.lastAuditCapture.restoration = restoration;
                        end
                    end
                    if isempty(dAug) || qpExitflag <= 0
                        exitflag = -1;
                        message = "QP subproblem failed. quadprog exitflag = " + string(qpExitflag);

                        if opts.ResetHessianOnFailure
                            obj.resetHessian(n);
                        end

                        break
                    end
                end

                d = dAug(1:n);

                [slackEqInf,slackIneqInf] = obj.extractSlackNorms(dAug,slackInfo);
                if restoration.attempted
                    slackEqInf = restoration.slackEqInf;
                    restorationUsed = restoration.accepted;
                end

                % QP multipliers corresponding to original nonlinear constraints.
                % Initialize them before any successful early-return branch.
                lambdaEq = zeros(numel(ceq),1);
                lambdaIneq = zeros(numel(c),1);
                lambdaLower = zeros(n,1);
                lambdaUpper = zeros(n,1);

                if isfield(qpLambda,'eqlin') && ~isempty(qpLambda.eqlin)
                    lambdaEq = qpLambda.eqlin(1:numel(ceq));
                end

                if isfield(qpLambda,'ineqlin') && ...
                        ~isempty(qpLambda.ineqlin) && numel(c) > 0
                    lambdaIneq = qpLambda.ineqlin(1:numel(c));
                end

                if isfield(qpLambda,'lower') && ~isempty(qpLambda.lower)
                    lambdaLower = qpLambda.lower(1:n);
                end

                if isfield(qpLambda,'upper') && ~isempty(qpLambda.upper)
                    lambdaUpper = qpLambda.upper(1:n);
                end
                
                stepInf = norm(d,inf);
                % --------------------------------------------------------------
                % Tiny-step exit before merit line search.
                % If the current NLP point is feasible and the QP says to move by an
                % insignificant amount, accept the current point as a usable RT-NMPC
                % solution instead of failing the line search.
                % --------------------------------------------------------------
                if stepInf <= opts.StepTolerance && feasInf <= opts.ConstraintTolerance
                
                    exitflag = 2;
                    message = "SQP stopped: QP step below StepTolerance and current point is feasible.";
                
                    fval = f;
                
                    output = obj.makeOutput(iter,funcCount,constrCount,fval,feasInf,kktInf, ...
                        stepInf,message,qpExitflag,qpOutput,1,slackEqInf,slackIneqInf);
                    if opts.CollectTiming
                        timing.iteration(iter).stepInf = stepInf;
                    end
                    output = obj.attachTiming(output,timing,totalTimer,iter);
                
                    lambda.eqnonlin = lambdaEq;
                    lambda.ineqnonlin = lambdaIneq;
                
                    obj.lastOutput = output;
                    obj.lastLambda = lambda;
                
                    if opts.Display == "final" || opts.Display == "iter"
                        fprintf('[SQP] exitflag=%d | iter=%d | f=%.6e | feas=%.3e | opt=%.3e | step=%.3e | %s\n', ...
                            exitflag,iter,fval,feasInf,kktInf,stepInf,message);
                    end
                
                    return
                end
                % -----------------------------------------------------
                % 2. Merit line search
                % -----------------------------------------------------
                lineSearchTimer = [];
                if opts.CollectTiming
                    lineSearchTimer = tic;
                end
                if opts.UseValueOnlyMeritEvaluation
                    phi0 = obj.meritFromValues( ...
                        f,c,ceq,z,lb,ub,opts.MeritPenalty);
                else
                    phi0 = obj.merit( ...
                        costFun,nonlFun,z,lb,ub,opts.MeritPenalty);
                end

                directionalObj = g.'*d;
                directionalTheta = obj.l1InfeasibilityDirectional( ...
                    c,ceq,z,lb,ub,gradc,gradceq,d);
                directionalMerit = directionalObj + ...
                    opts.MeritPenalty*directionalTheta;
                predictedDecrease = max(1e-12,-directionalObj);

                accepted = false;
                alpha = 1.0;
                phiTrial = nan;
                lineSearchTrialCount = 0;
                thetaTrial = nan;
                filterDecision = obj.emptyFilterDecision();
                lineSearchOwner = opts.LineSearchStrategy;
                if hybridEnabled
                    if directionalMerit < 0
                        lineSearchOwner = "l1_merit";
                    else
                        lineSearchOwner = "filter";
                    end
                end

                for ls = 1:opts.MaxLineSearch
                    lineSearchTrialCount = ls;
                    zTrial = z + alpha*d;
                    zTrial = min(max(zTrial,lb),ub);
                    [zTrial,socInfo] = obj.secondOrderCorrectTrial( ...
                        nonlFun,z,zTrial,d,alpha,lb,ub);

                    if filterEnabled || hybridEnabled
                        [fTrial,cTrial,ceqTrial] = obj.evalValues( ...
                            costFun,nonlFun,zTrial);
                        thetaTrial = obj.l1InfeasibilityFromValues( ...
                            cTrial,ceqTrial,zTrial,lb,ub);
                        phiTrial = fTrial + opts.MeritPenalty*thetaTrial;
                        if lineSearchOwner == "filter"
                            filterDecision = obj.filterTrialDecision( ...
                                f,thetaCurrent,directionalObj,alpha, ...
                                fTrial,thetaTrial,thetaMinimum,thetaMaximum, ...
                                filterTheta(1:filterCount), ...
                                filterObjective(1:filterCount));
                            accepted = filterDecision.accepted;
                        else
                            meritTol = opts.MeritAcceptTolerance * ...
                                max(1,abs(phi0));
                            accepted = phiTrial <= phi0 + ...
                                opts.LineSearchC1*alpha*directionalMerit + ...
                                meritTol;
                        end
                    else
                        phiTrial = obj.merit( ...
                            costFun,nonlFun,zTrial,lb,ub,opts.MeritPenalty);

                        meritTol = opts.MeritAcceptTolerance * max(1,abs(phi0));

                        accepted = ...
                            phiTrial <= phi0 - ...
                                opts.LineSearchC1*alpha*predictedDecrease + meritTol || ...
                            phiTrial <= phi0 + meritTol;
                    end

                    if accepted
                        break
                    end

                    alpha = obj.nextBacktrackingAlpha( ...
                        alpha,phi0,phiTrial,directionalMerit,opts);

                    if alpha < opts.MinAlpha
                        break
                    end
                end

                lineSearchSeconds = obj.elapsed( ...
                    lineSearchTimer,opts.CollectTiming);
                timing = obj.addPhaseTiming( ...
                    timing,iter,"lineSearchSeconds",lineSearchSeconds);
                if opts.CollectTiming
                    currentMeritEvaluationCount = ...
                        double(~opts.UseValueOnlyMeritEvaluation);
                    timing.lineSearchTrialCount = ...
                        timing.lineSearchTrialCount + lineSearchTrialCount;
                    timing.meritEvaluationCount = ...
                        timing.meritEvaluationCount + ...
                        currentMeritEvaluationCount + lineSearchTrialCount;
                    timing.currentMeritReuseCount = ...
                        timing.currentMeritReuseCount + ...
                        double(opts.UseValueOnlyMeritEvaluation);
                    timing.valueOnlyMeritTrialCount = ...
                        timing.valueOnlyMeritTrialCount + ...
                        double(opts.UseValueOnlyMeritEvaluation)* ...
                        lineSearchTrialCount;
                    timing.iteration(iter).lineSearchTrialCount = ...
                        lineSearchTrialCount;
                    timing.iteration(iter).meritEvaluationCount = ...
                        currentMeritEvaluationCount + lineSearchTrialCount;
                    timing.iteration(iter).currentMeritReused = ...
                        opts.UseValueOnlyMeritEvaluation;
                    timing.iteration(iter).valueOnlyMeritTrialCount = ...
                        double(opts.UseValueOnlyMeritEvaluation)* ...
                        lineSearchTrialCount;
                    timing.iteration(iter).initialMerit = phi0;
                    timing.iteration(iter).finalMerit = phiTrial;
                    timing.iteration(iter).acceptedAlpha = alpha;
                    timing.iteration(iter).lineSearchAccepted = accepted;
                    timing.iteration(iter).lineSearchStrategy = ...
                        opts.LineSearchStrategy;
                    timing.iteration(iter).lineSearchUpdate = ...
                        opts.LineSearchUpdate;
                    timing.iteration(iter).lineSearchOwner = ...
                        lineSearchOwner;
                    timing.iteration(iter).initialInfeasibilityL1 = ...
                        thetaCurrent;
                    timing.iteration(iter).finalInfeasibilityL1 = ...
                        thetaTrial;
                    timing.iteration(iter).objectiveDirectionalDerivative = ...
                        directionalObj;
                    timing.iteration(iter).infeasibilityDirectionalDerivative = ...
                        directionalTheta;
                    timing.iteration(iter).meritDirectionalDerivative = ...
                        directionalMerit;
                    timing.iteration(iter).meritDirectionIsDescent = ...
                        directionalMerit < 0;
                    timing.iteration(iter).filterAcceptable = ...
                        filterDecision.filterAcceptable;
                    timing.iteration(iter).filterSwitching = ...
                        filterDecision.switching;
                    timing.iteration(iter).filterArmijoObjective = ...
                        filterDecision.armijoObjective;
                    timing.iteration(iter).filterSufficientFeasibility = ...
                        filterDecision.sufficientFeasibility;
                    timing.iteration(iter).filterSufficientObjective = ...
                        filterDecision.sufficientObjective;
                    timing.iteration(iter).filterSize = filterCount;
                    timing.iteration(iter).filterRestorationRequired = ...
                        lineSearchOwner == "filter" && ~accepted;
                    timing.iteration(iter).secondOrderCorrectionApplied = ...
                        socInfo.applied;
                    timing.iteration(iter).secondOrderCorrectionNormInf = ...
                        socInfo.correctionNormInf;
                    timing.iteration(iter).rawTrialInfeasibilityL1 = ...
                        socInfo.rawInfeasibilityL1;
                    timing.filterTrialEvaluationCount = ...
                        timing.filterTrialEvaluationCount + ...
                        double(lineSearchOwner == "filter")* ...
                        lineSearchTrialCount;
                    timing.filterAcceptedTrialCount = ...
                        timing.filterAcceptedTrialCount + ...
                        double(lineSearchOwner == "filter" && accepted);
                    timing.filterRestorationRequiredCount = ...
                        timing.filterRestorationRequiredCount + ...
                        double(lineSearchOwner == "filter" && ~accepted);
                    timing.hybridMeritOwnedIterationCount = ...
                        timing.hybridMeritOwnedIterationCount + ...
                        double(hybridEnabled && lineSearchOwner == "l1_merit");
                    timing.hybridFilterOwnedIterationCount = ...
                        timing.hybridFilterOwnedIterationCount + ...
                        double(hybridEnabled && lineSearchOwner == "filter");
                end

                if ~accepted
                
                    % If the line search failed, but the current point is already feasible,
                    % this is not necessarily a catastrophic MPC failure. It usually means
                    % the QP direction did not provide a useful merit decrease from an
                    % already feasible warm-started trajectory.
                    if feasInf <= opts.ConstraintTolerance
                
                        exitflag = 2;
                        message = "SQP stopped: line search failed, but current point is feasible. Returning current point.";
                
                        fval = f;
                
                        output = obj.makeOutput(iter,funcCount,constrCount,fval,feasInf,kktInf, ...
                            stepInf,message,qpExitflag,qpOutput,alpha,slackEqInf,slackIneqInf);
                        if opts.CollectTiming
                            timing.iteration(iter).stepInf = stepInf;
                        end
                        output = obj.attachTiming(output,timing,totalTimer,iter);
                
                        lambda.eqnonlin = lambdaEq;
                        lambda.ineqnonlin = lambdaIneq;
                
                        obj.lastOutput = output;
                        obj.lastLambda = lambda;
                
                        if opts.Display == "final" || opts.Display == "iter"
                            fprintf('[SQP] exitflag=%d | iter=%d | f=%.6e | feas=%.3e | opt=%.3e | step=%.3e | %s\n', ...
                                exitflag,iter,fval,feasInf,kktInf,stepInf,message);
                        end
                
                        return
                    end
                
                    % If infeasible, then this really is a line-search failure.
                    exitflag = -2;
                    if lineSearchOwner == "filter"
                        message = ["Filter line search found no acceptable step ", ...
                            "and current point is infeasible; restoration is not ", ...
                            "enabled in this audit scope."];
                    else
                        message = "Merit line search failed and current point is infeasible.";
                    end
                
                    if opts.ResetHessianOnFailure
                        obj.resetHessian(n);
                    end
                
                    break
                end

                if lineSearchOwner == "filter" && ...
                        ~(filterDecision.switching && ...
                          filterDecision.armijoObjective)
                    filterCount = filterCount + 1;
                    filterTheta(filterCount) = ...
                        (1-opts.FilterGammaTheta)*thetaCurrent;
                    filterObjective(filterCount) = ...
                        f-opts.FilterGammaObjective*thetaCurrent;
                end
                % -----------------------------------------------------
                % 3. Accept step and evaluate new NLP data
                % -----------------------------------------------------
                % Preserve the exact trial accepted by the globalization
                % test.  In SOC mode zTrial includes the nonlinear
                % feasibility correction; rebuilding z+alpha*d here would
                % silently discard that accepted correction.
                zNew = zTrial;

                try
                    [fNew,gNew,cNew,ceqNew,gradcNew,gradceqNew,evalTiming] = ...
                        obj.evalAll(costFun,nonlFun,zNew);
                    timing = obj.addEvaluationTiming(timing,evalTiming,iter);
                    funcCount = funcCount + 1;
                    constrCount = constrCount + 1;
                catch ME
                    exitflag = -3;
                    message = "Function evaluation after step failed: " + string(ME.message);
                    break
                end

                feasInfNew = obj.constraintViolation(cNew,ceqNew,zNew,lb,ub);
                thetaCurrent = obj.l1InfeasibilityFromValues( ...
                    cNew,ceqNew,zNew,lb,ub);

                % -----------------------------------------------------
                % 4. Damped BFGS Hessian update
                % -----------------------------------------------------
                hessianTimer = [];
                if opts.CollectTiming
                    hessianTimer = tic;
                end
                if opts.HessianMode == "bfgs"
                    s = zNew - z;

                    % gradLagOld = obj.lagrangianGradient(g,gradc,gradceq,lambdaIneq,lambdaEq);
                    gradLagOld = obj.lagrangianGradient(g,gradc,gradceq,lambdaIneq,lambdaEq,lambdaLower,lambdaUpper);
                    % gradLagNew = obj.lagrangianGradient(gNew,gradcNew,gradceqNew,lambdaIneq,lambdaEq);
                    gradLagNew = obj.lagrangianGradient(gNew,gradcNew,gradceqNew,lambdaIneq,lambdaEq,lambdaLower,lambdaUpper);

                    y = gradLagNew - gradLagOld;

                    obj.updateBFGS(s,y);
                elseif opts.HessianMode == "lbfgs"
                    s = zNew-z;
                    gradLagOld = obj.lagrangianGradient( ...
                        g,gradc,gradceq,lambdaIneq,lambdaEq, ...
                        lambdaLower,lambdaUpper);
                    gradLagNew = obj.lagrangianGradient( ...
                        gNew,gradcNew,gradceqNew,lambdaIneq,lambdaEq, ...
                        lambdaLower,lambdaUpper);
                    obj.updateLimitedMemoryBfgs(s,gradLagNew-gradLagOld,Hq);
                end
                hessianSeconds = obj.elapsed(hessianTimer,opts.CollectTiming);
                timing = obj.addPhaseTiming( ...
                    timing,iter,"hessianWorkSeconds",hessianSeconds);

                % -----------------------------------------------------
                % 5. Update iterate
                % -----------------------------------------------------
                zPrev = z;
                gPrev = g;
                gradcPrev = gradc;
                gradceqPrev = gradceq;
                lambdaEqPrev = lambdaEq;
                lambdaIneqPrev = lambdaIneq;

                z = zNew;
                f = fNew;
                g = gNew;
                c = cNew;
                ceq = ceqNew;
                gradc = gradcNew;
                gradceq = gradceqNew;

                lambda.eqnonlin = lambdaEq;
                lambda.ineqnonlin = lambdaIneq;
                lambda.lower = lambdaLower;
                lambda.upper = lambdaUpper;

                gradLag = obj.lagrangianGradient( ...
                    g,gradc,gradceq,lambdaIneq,lambdaEq,lambdaLower,lambdaUpper);
                kktInf = norm(gradLag,inf);
                feasInf = feasInfNew;

                if opts.CollectTiming
                    timing.iteration(iter).objective = f;
                    timing.iteration(iter).constraintViolationInf = feasInf;
                    timing.iteration(iter).firstOrderOptimalityInf = kktInf;
                    timing.iteration(iter).stepInf = stepInf;
                end

                if opts.Display == "iter"
                    obj.printIter(iter,f,feasInf,kktInf,stepInf,alpha, ...
                        slackEqInf,slackIneqInf,qpExitflag);
                end
                obj.updateLiveIterationPlot(livePlot,iter,f,feasInf,kktInf, ...
                    stepInf,alpha);

                % -----------------------------------------------------
                % 6. Exit checks
                % -----------------------------------------------------
                if feasInf <= opts.ConstraintTolerance && ...
                   kktInf <= opts.OptimalityTolerance

                    exitflag = 1;
                    message = "SQP converged: first-order optimality and feasibility satisfied.";
                    break
                end

                if feasInf <= opts.ConstraintTolerance && ...
                   stepInf <= opts.StepTolerance

                    exitflag = 2;
                    message = "SQP converged: step tolerance satisfied and feasible.";
                    break
                end

                if stepInf <= opts.StepTolerance && feasInf > opts.ConstraintTolerance
                    exitflag = -2;
                    message = "SQP stopped: step tolerance reached but constraints not satisfied.";
                    break
                end
            end

            fval = f;
            reportedIterations = min(iter,opts.MaxIterations + ...
                double(restorationUsed)* ...
                opts.FeasibilityRestorationPostIterations);

            if opts.Display == "final" || opts.Display == "iter"
                fprintf('[SQP] exitflag=%d | iter=%d | f=%.6e | feas=%.3e | opt=%.3e | step=%.3e | %s\n', ...
                    exitflag,reportedIterations,fval,feasInf,kktInf, ...
                    stepInf,message);
            end

            output = obj.makeOutput(reportedIterations,funcCount,constrCount, ...
                fval,feasInf,kktInf, ...
                stepInf,message,qpExitflag,qpOutput,alpha,slackEqInf,slackIneqInf);
            if opts.AuditCapture
                output.auditCapture = obj.lastAuditCapture;
            end
            output = obj.attachTiming( ...
                output,timing,totalTimer,reportedIterations);

            obj.lastOutput = output;
            obj.lastLambda = lambda;
        end

        %--------------------------------------------------------------
        function report = checkGradients(obj,costFun,nonlFun,z,lb,ub)
        %CHECKGRADIENTS Directional finite-difference check.
        %
        % This is intentionally directional instead of full coordinate
        % finite-difference because NMPC decision vectors can be very large.

            opts = obj.options;
            z = z(:);
            n = numel(z);

            lb = lb(:);
            ub = ub(:);

            if isempty(lb), lb = -inf(n,1); end
            if isempty(ub), ub =  inf(n,1); end

            z = min(max(z,lb),ub);

            [f,g,c,ceq,gradc,gradceq] = obj.evalAll(costFun,nonlFun,z);

            nDir = opts.GradientCheckDirections;
            h0 = opts.GradientCheckStep;

            objErr = nan(nDir,1);
            cErr   = nan(nDir,1);
            ceqErr = nan(nDir,1);

            rng(1);

            for j = 1:nDir
                d = randn(n,1);
                d = d / max(norm(d),eps);

                h = obj.boundSafeStep(z,d,h0,lb,ub);

                if h <= eps
                    continue
                end

                zp = z + h*d;
                zm = z - h*d;

                fp = costFun(zp);
                fm = costFun(zm);

                objFD = (fp - fm)/(2*h);
                objAN = g.'*d;

                denObj = max([1, abs(objFD), abs(objAN)]);
                objErr(j) = abs(objFD - objAN) / denObj;
                
                [cp,ceqp] = nonlFun(zp);
                [cm,ceqm] = nonlFun(zm);

                if ~isempty(c)
                    cFD = (cp(:) - cm(:))/(2*h);
                    cAN = gradc.'*d;

                    denC = max([1, norm(cFD,inf), norm(cAN,inf)]);
                    cErr(j) = norm(cFD - cAN,inf) / denC;
                end

                if ~isempty(ceq)
                    ceqFD = (ceqp(:) - ceqm(:))/(2*h);
                    ceqAN = gradceq.'*d;

                    denCeq = max([1, norm(ceqFD,inf), norm(ceqAN,inf)]);
                    ceqErr(j) = norm(ceqFD - ceqAN,inf) / denCeq;
                end
            end

            report = struct();
            report.objectiveRelErrMax = max(objErr,[],'omitnan');
            report.ineqRelErrMax      = max(cErr,[],'omitnan');
            report.eqRelErrMax        = max(ceqErr,[],'omitnan');
            report.objectiveRelErr    = objErr;
            report.ineqRelErr         = cErr;
            report.eqRelErr           = ceqErr;

            obj.lastGradCheck = report;

            if opts.GradientCheckVerbose
                fprintf('\n');
                fprintf('======================================================================\n');
                fprintf(' CUSTOM SQP DIRECTIONAL GRADIENT CHECK\n');
                fprintf('======================================================================\n');
                fprintf('  directions checked       : %d\n', nDir);
                fprintf('  max objective rel error  : %.3e\n', report.objectiveRelErrMax);
                fprintf('  max inequality rel error : %.3e\n', report.ineqRelErrMax);
                fprintf('  max equality rel error   : %.3e\n', report.eqRelErrMax);
                fprintf('======================================================================\n\n');
            end
        end
    end

    methods(Access=private)
        %--------------------------------------------------------------
        function tf = preparedInitialMatches(~,prepared,z,lb,ub)
            tf = isstruct(prepared) && isscalar(prepared) && ...
                isfield(prepared,'valid') && logical(prepared.valid) && ...
                all(isfield(prepared,{'z','lb','ub','c','ceq','gradc','gradceq'}));
            if ~tf
                return
            end
            tf = isequaln(prepared.z(:),z(:)) && ...
                isequaln(prepared.lb(:),lb(:)) && ...
                isequaln(prepared.ub(:),ub(:)) && ...
                all(isfinite(prepared.c(:))) && all(isfinite(prepared.ceq(:))) && ...
                all(isfinite(prepared.gradc),'all') && ...
                all(isfinite(prepared.gradceq),'all') && ...
                size(prepared.gradc,1) == numel(z) && ...
                size(prepared.gradc,2) == numel(prepared.c) && ...
                size(prepared.gradceq,1) == numel(z) && ...
                size(prepared.gradceq,2) == numel(prepared.ceq);
        end

        %--------------------------------------------------------------
        function opts = mergeOptions(~,defaults,user)
            opts = defaults;

            if isempty(user)
                return
            end

            f = fieldnames(user);

            for i = 1:numel(f)
                opts.(f{i}) = user.(f{i});
            end

            if ~isfield(opts,'QPOptions') || isempty(opts.QPOptions)
                opts.QPOptions = defaults.QPOptions;
            end
        end

        %--------------------------------------------------------------
        function [f,g,c,ceq,gradc,gradceq,timing] = ...
                evalAll(obj,costFun,nonlFun,z)
            timing = struct( ...
                'objectiveEvaluationSeconds',0, ...
                'constraintEvaluationSeconds',0);

            objectiveTimer = [];
            if obj.options.CollectTiming
                objectiveTimer = tic;
            end
            [f,g] = costFun(z);
            timing.objectiveEvaluationSeconds = obj.elapsed( ...
                objectiveTimer,obj.options.CollectTiming);

            g = g(:);

            constraintTimer = [];
            if obj.options.CollectTiming
                constraintTimer = tic;
            end
            [c,ceq,gradc,gradceq] = nonlFun(z);
            timing.constraintEvaluationSeconds = obj.elapsed( ...
                constraintTimer,obj.options.CollectTiming);

            c = c(:);
            ceq = ceq(:);

            n = numel(z);

            if isempty(gradc)
                gradc = zeros(n,0);
            end

            if isempty(gradceq)
                gradceq = zeros(n,0);
            end

            if size(gradc,1) ~= n && size(gradc,2) == n
                gradc = gradc.';
            end

            if size(gradceq,1) ~= n && size(gradceq,2) == n
                gradceq = gradceq.';
            end

            if size(gradc,1) ~= n
                error('SQPSolver:GradientDimension', ...
                    'gradc must be nVar x nIneq.');
            end

            if size(gradceq,1) ~= n
                error('SQPSolver:GradientDimension', ...
                    'gradceq must be nVar x nEq.');
            end

            if size(gradc,2) ~= numel(c)
                error('SQPSolver:GradientDimension', ...
                    'gradc has %d columns but c has %d entries.', ...
                    size(gradc,2),numel(c));
            end

            if size(gradceq,2) ~= numel(ceq)
                error('SQPSolver:GradientDimension', ...
                    'gradceq has %d columns but ceq has %d entries.', ...
                    size(gradceq,2),numel(ceq));
            end
        end

        %--------------------------------------------------------------
        function timing = emptyTiming(obj,maxIterations)
            iterationTemplate = struct( ...
                'iteration',nan, ...
                'objectiveEvaluationSeconds',0, ...
                'constraintEvaluationSeconds',0, ...
                'hessianWorkSeconds',0, ...
                'qpAssemblySeconds',0, ...
                'qpSolveSeconds',0, ...
                'lineSearchSeconds',0, ...
                'lineSearchTrialCount',0, ...
                'meritEvaluationCount',0, ...
                'currentMeritReused',false, ...
                'valueOnlyMeritTrialCount',0, ...
                'initialMerit',nan, ...
                'finalMerit',nan, ...
                'acceptedAlpha',nan, ...
                'lineSearchAccepted',false, ...
                'lineSearchStrategy',"l1_merit", ...
                'lineSearchUpdate',"fixed_beta", ...
                'lineSearchOwner',"l1_merit", ...
                'initialInfeasibilityL1',nan, ...
                'finalInfeasibilityL1',nan, ...
                'objectiveDirectionalDerivative',nan, ...
                'infeasibilityDirectionalDerivative',nan, ...
                'meritDirectionalDerivative',nan, ...
                'meritDirectionIsDescent',false, ...
                'filterAcceptable',false, ...
                'filterSwitching',false, ...
                'filterArmijoObjective',false, ...
                'filterSufficientFeasibility',false, ...
                'filterSufficientObjective',false, ...
                'filterSize',0, ...
                'filterRestorationRequired',false, ...
                'secondOrderCorrectionApplied',false, ...
                'secondOrderCorrectionNormInf',0, ...
                'rawTrialInfeasibilityL1',nan, ...
                'objective',nan, ...
                'constraintViolationInf',nan, ...
                'firstOrderOptimalityInf',nan, ...
                'stepInf',nan, ...
                'qpExitflag',nan);
            timing = struct( ...
                'objectiveEvaluationSeconds',0, ...
                'constraintEvaluationSeconds',0, ...
                'hessianWorkSeconds',0, ...
                'qpAssemblySeconds',0, ...
                'qpSolveSeconds',0, ...
                'lineSearchSeconds',0, ...
                'lineSearchTrialCount',0, ...
                'meritEvaluationCount',0, ...
                'currentMeritReuseCount',0, ...
                'valueOnlyMeritTrialCount',0, ...
                'filterTrialEvaluationCount',0, ...
                'filterAcceptedTrialCount',0, ...
                'filterRestorationRequiredCount',0, ...
                'hybridMeritOwnedIterationCount',0, ...
                'hybridFilterOwnedIterationCount',0, ...
                'measuredPhaseSeconds',0, ...
                'unattributedSeconds',0, ...
                'totalSolveSeconds',0, ...
                'iteration',repmat(iterationTemplate,maxIterations,1));
            if ~obj.options.CollectTiming
                timing = struct();
            end
        end

        %--------------------------------------------------------------
        function timing = addEvaluationTiming(obj,timing,increment,iter)
            if ~obj.options.CollectTiming
                return
            end
            timing.objectiveEvaluationSeconds = ...
                timing.objectiveEvaluationSeconds + ...
                increment.objectiveEvaluationSeconds;
            timing.constraintEvaluationSeconds = ...
                timing.constraintEvaluationSeconds + ...
                increment.constraintEvaluationSeconds;
            if nargin >= 4 && iter >= 1
                timing.iteration(iter).objectiveEvaluationSeconds = ...
                    timing.iteration(iter).objectiveEvaluationSeconds + ...
                    increment.objectiveEvaluationSeconds;
                timing.iteration(iter).constraintEvaluationSeconds = ...
                    timing.iteration(iter).constraintEvaluationSeconds + ...
                    increment.constraintEvaluationSeconds;
            end
        end

        %--------------------------------------------------------------
        function timing = addPhaseTiming(obj,timing,iter,field,seconds)
            if ~obj.options.CollectTiming
                return
            end
            timing.(field) = timing.(field) + seconds;
            timing.iteration(iter).(field) = ...
                timing.iteration(iter).(field) + seconds;
        end

        %--------------------------------------------------------------
        function seconds = elapsed(~,timerValue,enabled)
            if enabled
                seconds = toc(timerValue);
            else
                seconds = 0;
            end
        end

        %--------------------------------------------------------------
        function output = attachTiming(obj,output,timing,totalTimer,iter)
            if ~obj.options.CollectTiming
                return
            end

            timing.totalSolveSeconds = toc(totalTimer);
            phaseFields = ["objectiveEvaluationSeconds", ...
                "constraintEvaluationSeconds","hessianWorkSeconds", ...
                "qpAssemblySeconds","qpSolveSeconds", ...
                "lineSearchSeconds"];
            timing.measuredPhaseSeconds = 0;
            for ii = 1:numel(phaseFields)
                timing.measuredPhaseSeconds = timing.measuredPhaseSeconds + ...
                    timing.(phaseFields(ii));
            end
            timing.unattributedSeconds = max(0, ...
                timing.totalSolveSeconds-timing.measuredPhaseSeconds);

            if iter > 0
                timing.iteration = timing.iteration(1:iter);
            else
                timing.iteration = timing.iteration([]);
            end
            output.timing = timing;
        end

        %--------------------------------------------------------------
        function [fullStep,info] = solveSelectedCondensedQp( ...
                obj,audit,linearizationPoint,transformation, ...
                eliminatedColumns,eliminatedRows)
        %SOLVESELECTEDCONDENSEDQP Eliminate a square equality-coordinate block.
            preparationTimer = tic;
            n = numel(linearizationPoint);
            assert(isequal(size(audit.hessian),[n n]) && ...
                numel(audit.gradient)==n && ...
                size(audit.gradientEquality,1)==n && ...
                size(audit.gradientInequality,1)==n, ...
                'SQPSolver:CondensedRtiPacket', ...
                'The condensed RTI derivative packet has invalid dimensions.');
            assert(all(isfinite(audit.hessian),'all') && ...
                all(isfinite(audit.gradient)) && ...
                all(isfinite(audit.equality)) && ...
                all(isfinite(audit.gradientEquality),'all') && ...
                all(isfinite(audit.gradientInequality),'all'), ...
                'SQPSolver:CondensedRtiPacketFinite', ...
                'The condensed RTI derivative packet contains nonfinite data.');
            Hchart = sparse(transformation.'*audit.hessian*transformation);
            Hchart = 0.5*(Hchart+Hchart.');
            gChart = transformation.'*audit.gradient(:);
            AeqChart = sparse(audit.gradientEquality.'*transformation);
            beqChart = -audit.equality(:);
            assert(all(eliminatedRows>=1) && ...
                all(eliminatedRows<=size(AeqChart,1)) && ...
                numel(unique(eliminatedRows))==numel(eliminatedRows), ...
                'SQPSolver:CondensedRtiElimination', ...
                'The eliminated RTI equality-row indices are invalid.');
            inequality = audit.inequality(:);
            inactiveInequality = isinf(inequality) & inequality<0;
            assert(~any(isnan(inequality)) && ...
                ~any(isinf(inequality) & inequality>0), ...
                'SQPSolver:CondensedRtiInequalityFinite', ...
                ['Only negative-infinite, inactive nonlinear inequality ', ...
                 'rows are permitted.']);
            activeInequality = ~inactiveInequality;
            if ~any(activeInequality)
                Achart = sparse(0,size(transformation,2));
                bChart = zeros(0,1);
            else
                Achart = sparse(audit.gradientInequality(:, ...
                    activeInequality).'*transformation);
                bChart = -inequality(activeInequality);
            end
            upperStep = audit.upperBound(:)-linearizationPoint(:);
            lowerStep = linearizationPoint(:)-audit.lowerBound(:);
            finiteUpper = isfinite(upperStep);
            finiteLower = isfinite(lowerStep);
            Achart = [Achart;transformation(finiteUpper,:); ...
                -transformation(finiteLower,:)];
            bChart = [bChart;upperStep(finiteUpper);lowerStep(finiteLower)];

            keptColumns = setdiff(1:size(transformation,2), ...
                eliminatedColumns,'stable');
            keptRows = setdiff(1:size(AeqChart,1), ...
                eliminatedRows,'stable');
            Aee = AeqChart(eliminatedRows,eliminatedColumns);
            Aek = AeqChart(eliminatedRows,keptColumns);
            assert(size(Aee,1)==size(Aee,2), ...
                'SQPSolver:CondensedRtiSquareBlock', ...
                'The selected condensed RTI elimination block is not square.');
            eliminatedOffset = Aee\beqChart(eliminatedRows);
            eliminatedMap = -(Aee\Aek);
            assert(all(isfinite(eliminatedOffset)) && ...
                all(isfinite(nonzeros(eliminatedMap))), ...
                'SQPSolver:CondensedRtiEliminationFinite', ...
                'The condensed RTI elimination map is nonfinite.');
            offset = zeros(size(transformation,2),1);
            offset(eliminatedColumns) = eliminatedOffset;
            expansion = spalloc(size(transformation,2),numel(keptColumns), ...
                numel(keptColumns)+nnz(eliminatedMap));
            expansion(keptColumns,:) = speye(numel(keptColumns));
            expansion(eliminatedColumns,:) = eliminatedMap;

            reducedHessian = sparse(expansion.'*Hchart*expansion);
            reducedHessian = 0.5*(reducedHessian+reducedHessian.');
            g = expansion.'*(Hchart*offset+gChart);
            A = Achart*expansion;
            b = bChart-Achart*offset;
            Aeq = AeqChart(keptRows,:)*expansion;
            beq = beqChart(keptRows)-AeqChart(keptRows,:)*offset;
            [Aeq,beq] = obj.dropNumericallyZeroEqualities(Aeq,beq);
            condensingSeconds = toc(preparationTimer);

            qpTimer = tic;
            [reducedStep,qpExitflag,qpOutput,qpLambda,qpMethod] = ...
                obj.solveCondensedReducedQp( ...
                    reducedHessian,g,A,b,Aeq,beq);
            qpSeconds = toc(qpTimer);
            if isempty(reducedStep) || any(~isfinite(reducedStep))
                fullStep = [];
                linearEqualityResidual = inf;
                linearInequalityResidual = inf;
            else
                chartStep = offset+expansion*reducedStep;
                fullStep = transformation*chartStep;
                fullEqualityResidual = ...
                    audit.gradientEquality.'*fullStep+audit.equality(:);
                if isempty(audit.inequality)
                    fullInequalityResidual = zeros(0,1);
                else
                    fullInequalityResidual = ...
                        audit.gradientInequality.'*fullStep+ ...
                        audit.inequality(:);
                end
                linearEqualityResidual = norm(fullEqualityResidual,inf);
                linearInequalityResidual = max([fullInequalityResidual; ...
                    fullStep-upperStep;-fullStep-lowerStep;0]);
            end
            qpIterations = nan;
            qpFirstOrderOptimality = nan;
            reducedKktInf = nan;
            qpComplementarity = nan;
            qpMinimumLinearSlack = nan;
            qpMaximumInactiveMultiplier = nan;
            if isstruct(qpOutput) && isfield(qpOutput,'iterations') && ...
                    isscalar(qpOutput.iterations)
                qpIterations = double(qpOutput.iterations);
            end
            if isstruct(qpOutput) && isfield(qpOutput,'firstorderopt') && ...
                    isscalar(qpOutput.firstorderopt)
                qpFirstOrderOptimality = double(qpOutput.firstorderopt);
            end
            if ~isempty(reducedStep) && all(isfinite(reducedStep)) && ...
                    isstruct(qpLambda)
                reducedGradient = reducedHessian*reducedStep+g;
                if isfield(qpLambda,'ineqlin') && ~isempty(qpLambda.ineqlin)
                    reducedGradient = reducedGradient+ ...
                        A.'*qpLambda.ineqlin(:);
                end
                if isfield(qpLambda,'eqlin') && ~isempty(qpLambda.eqlin)
                    reducedGradient = reducedGradient+ ...
                        Aeq.'*qpLambda.eqlin(:);
                end
                if isfield(qpLambda,'lower') && ~isempty(qpLambda.lower)
                    reducedGradient = reducedGradient-qpLambda.lower(:);
                end
                if isfield(qpLambda,'upper') && ~isempty(qpLambda.upper)
                    reducedGradient = reducedGradient+qpLambda.upper(:);
                end
                reducedKktInf = norm(reducedGradient,inf);
                if ~isempty(A) && isfield(qpLambda,'ineqlin') && ...
                        numel(qpLambda.ineqlin)==size(A,1)
                    linearSlack = b-A*reducedStep;
                    linearMultiplier = qpLambda.ineqlin(:);
                    qpComplementarity = max(abs( ...
                        linearSlack.*linearMultiplier));
                    qpMinimumLinearSlack = min(linearSlack);
                    inactive = linearSlack>1e-7;
                    if any(inactive)
                        qpMaximumInactiveMultiplier = ...
                            max(abs(linearMultiplier(inactive)));
                    else
                        qpMaximumInactiveMultiplier = 0;
                    end
                end
            end
            info = struct('iteration',0,'derivativeSeconds',0, ...
                'objectiveDerivativeSeconds',0, ...
                'constraintDerivativeSeconds',0, ...
                'condensingSeconds',condensingSeconds, ...
                'qpSeconds',qpSeconds,'qpExitflag',qpExitflag, ...
                'qpIterations',qpIterations, ...
                'qpFirstOrderOptimality',qpFirstOrderOptimality, ...
                'reducedKktInf',reducedKktInf, ...
                'qpComplementarity',qpComplementarity, ...
                'qpMinimumLinearSlack',qpMinimumLinearSlack, ...
                'qpMaximumInactiveMultiplier', ...
                    qpMaximumInactiveMultiplier, ...
                'qpMethod',qpMethod, ...
                'decisionCount',numel(keptColumns), ...
                'equalityCount',size(Aeq,1), ...
                'inequalityCount',size(A,1), ...
                'linearEqualityResidualInf',linearEqualityResidual, ...
                'linearInequalityResidualMax',linearInequalityResidual, ...
                'stepInf',norm(fullStep,inf));
        end

        %--------------------------------------------------------------
        function [step,exitflag,output,lambda,method] = ...
                solveCondensedReducedQp(obj,H,g,A,b,Aeq,beq)
        %SOLVECONDENSEDREDUCEDQP Select the default-inactive RTI QP path.
            mode = string(obj.options.CondensedRtiQpMode);
            assert(isscalar(mode) && ismember(mode,[ ...
                "configured_quadprog","direct_then_active_set"]), ...
                'SQPSolver:CondensedRtiQpMode', ...
                'The condensed RTI QP mode is invalid.');
            if mode=="direct_then_active_set" && isempty(Aeq)
                directTimer = tic;
                direct = -lsqminnorm(full(H),g,1e-12);
                residual = norm(H*direct+g,inf);
                violation = max([A*direct-b;0]);
                if all(isfinite(direct)) && residual<=1e-8 && ...
                        violation<=1e-10
                    step = direct;
                    exitflag = 1;
                    output = struct('iterations',0, ...
                        'firstorderopt',residual, ...
                        'algorithm','direct reduced Newton', ...
                        'message','The direct reduced Newton step is feasible.', ...
                        'directSolveSeconds',toc(directTimer));
                    lambda = struct('ineqlin',zeros(size(A,1),1), ...
                        'eqlin',zeros(0,1),'lower',zeros(numel(g),1), ...
                        'upper',zeros(numel(g),1));
                    method = "direct_newton";
                    return
                end
            end

            zero = zeros(numel(g),1);
            zeroFeasible = max([A*zero-b;0])<=1e-10 && ...
                norm(Aeq*zero-beq,inf)<=1e-10;
            if mode=="direct_then_active_set" && zeroFeasible
                activeOptions = optimoptions('quadprog', ...
                    'Algorithm','active-set','Display','off', ...
                    'ConstraintTolerance',1e-10, ...
                    'StepTolerance',1e-12,'MaxIterations',500);
                [step,~,exitflag,output,lambda] = quadprog( ...
                    H,g,A,b,Aeq,beq,[],[],zero,activeOptions);
                if ~isempty(step) && exitflag>0 && all(isfinite(step))
                    method = "active_set";
                    return
                end
            end

            [step,~,exitflag,output,lambda] = quadprog( ...
                H,g,A,b,Aeq,beq,[],[],[],obj.options.QPOptions);
            method = "configured_quadprog";
        end

        %--------------------------------------------------------------
        function [Aeq,beq] = dropNumericallyZeroEqualities(~,Aeq,beq)
            rowScale = full(max(abs(Aeq),[],2));
            nonzeroRightHandSide = rowScale<=1e-13 & abs(beq)>1e-12;
            assert(~any(nonzeroRightHandSide), ...
                'SQPSolver:CondensedRtiInconsistentEquality', ...
                ['A zero condensed RTI equality row has a nonzero right-', ...
                 'hand side. Repair the warm start before condensing.']);
            keep = rowScale>1e-13 | abs(beq)>1e-12;
            Aeq = Aeq(keep,:);
            beq = beq(keep);
        end

        %--------------------------------------------------------------
        function [Haug,faug,Aineq,bineq,Aeq,beq,lbq,ubq,slackInfo] = ...
                buildQP(obj,H,g,c,ceq,gradc,gradceq,z,lb,ub)

            opts = obj.options;

            n = numel(z);
            meq = numel(ceq);
            mineq = numel(c);

            Aeq0 = gradceq.';
            beq0 = -ceq;

            Aineq0 = gradc.';
            bineq0 = -c;

            dlb = lb - z;
            dub = ub - z;

            if ~opts.ElasticMode
                Haug = H;
                faug = g;

                Aeq = Aeq0;
                beq = beq0;

                Aineq = Aineq0;
                bineq = bineq0;

                lbq = dlb;
                ubq = dub;

                slackInfo = struct('n',n,'meq',meq,'mineq',mineq, ...
                                   'elastic',false,'idxSp',[], ...
                                   'idxSm',[],'idxSi',[]);
                return
            end

            nSlack = 2*meq + mineq;
            nAug = n + nSlack;

            idxSp = n + (1:meq);
            idxSm = n + meq + (1:meq);
            idxSi = n + 2*meq + (1:mineq);

            Haug = blkdiag(sparse(H), opts.ElasticSlackQuad*speye(nSlack));

            faug = zeros(nAug,1);
            faug(1:n) = g;

            if meq > 0
                faug(idxSp) = opts.ElasticPenaltyEq;
                faug(idxSm) = opts.ElasticPenaltyEq;
            end

            if mineq > 0
                faug(idxSi) = opts.ElasticPenaltyIneq;
            end

            if meq > 0
                Aeq = spalloc(meq,nAug,nnz(Aeq0)+2*meq);
                Aeq(:,1:n) = sparse(Aeq0);
                Aeq(:,idxSp) = speye(meq);
                Aeq(:,idxSm) = -speye(meq);
                beq = beq0;
            else
                Aeq = [];
                beq = [];
            end

            if mineq > 0
                Aineq = spalloc(mineq,nAug,nnz(Aineq0)+mineq);
                Aineq(:,1:n) = sparse(Aineq0);
                Aineq(:,idxSi) = -speye(mineq);
                bineq = bineq0;
            else
                Aineq = [];
                bineq = [];
            end

            lbq = [-inf(nAug,1)];
            ubq = [ inf(nAug,1)];

            lbq(1:n) = dlb;
            ubq(1:n) = dub;

            if nSlack > 0
                lbq(n+1:end) = 0;
            end

            slackInfo = struct('n',n,'meq',meq,'mineq',mineq, ...
                               'elastic',true,'idxSp',idxSp, ...
                               'idxSm',idxSm,'idxSi',idxSi);
        end

        %--------------------------------------------------------------
        function [slackEqInf,slackIneqInf] = extractSlackNorms(~,dAug,slackInfo)
            slackEqInf = 0;
            slackIneqInf = 0;

            if ~slackInfo.elastic
                return
            end

            if slackInfo.meq > 0
                sp = dAug(slackInfo.idxSp);
                sm = dAug(slackInfo.idxSm);
                slackEqInf = max(abs(sp-sm),[],'omitnan');
            end

            if slackInfo.mineq > 0
                si = dAug(slackInfo.idxSi);
                slackIneqInf = max(abs(si),[],'omitnan');
            end
        end

        %--------------------------------------------------------------
        function [dAug,exitflag,qpOutput,slackInfo,restoration] = ...
                solveFeasibilityRestoration(obj,c,ceq,gradc,gradceq,z,lb,ub,qpOptions)
            %SOLVEFEASIBILITYRESTORATION One bounded Phase-I equality step.
            % The nonlinear candidate is still checked by the unchanged hard
            % SQP acceptance path after this step.

            n = numel(z);
            meq = numel(ceq);
            mineq = numel(c);
            Aeq0 = gradceq.';
            beq0 = -ceq;
            Aineq0 = gradc.';
            bineq0 = -c;
            nAug = n + 2*meq;
            idxSp = n + (1:meq);
            idxSm = n + meq + (1:meq);

            Haug = blkdiag(speye(n), ...
                obj.options.FeasibilityRestorationPenaltyL2*speye(2*meq));
            faug = zeros(nAug,1);
            faug([idxSp idxSm]) = obj.options.FeasibilityRestorationPenaltyL1;

            if meq > 0
                Aeq = spalloc(meq,nAug,nnz(Aeq0)+2*meq);
                Aeq(:,1:n) = sparse(Aeq0);
                Aeq(:,idxSp) = -speye(meq);
                Aeq(:,idxSm) = speye(meq);
            else
                Aeq = [];
                beq0 = [];
            end
            if mineq > 0
                Aineq = spalloc(mineq,nAug,nnz(Aineq0));
                Aineq(:,1:n) = sparse(Aineq0);
            else
                Aineq = [];
                bineq0 = [];
            end

            lbq = [lb-z; zeros(2*meq,1)];
            ubq = [ub-z; inf(2*meq,1)];
            slackInfo = struct('n',n,'meq',0,'mineq',0, ...
                'elastic',false,'idxSp',[],'idxSm',[],'idxSi',[]);
            restoration = struct('attempted',true,'accepted',false, ...
                'exitflag',nan,'output',struct(),'slackEqInf',nan, ...
                'qpVariableCount',nAug,'qpInequalityCount',size(Aineq,1), ...
                'qpEqualityCount',size(Aeq,1));
            try
                [dAug,~,exitflag,qpOutput] = quadprog(Haug,faug,Aineq, ...
                    bineq0,Aeq,beq0,lbq,ubq,[],qpOptions);
            catch ME
                dAug = [];
                exitflag = -1;
                qpOutput = struct('message',ME.message);
            end
            restoration.exitflag = exitflag;
            restoration.output = qpOutput;
            restoration.accepted = ~isempty(dAug) && exitflag > 0;
            if restoration.accepted && meq > 0
                restoration.slackEqInf = norm(dAug(idxSp)-dAug(idxSm),inf);
            end
        end

        %--------------------------------------------------------------
        function Hpd = makePositiveDefinite(obj,H)
            opts = obj.options;

            Hpd = 0.5*(H+H.');

            if isempty(Hpd)
                return
            end

            n = size(Hpd,1);

            reg = opts.HessianRegularization;

            for attempt = 1:8
                [~,p] = chol(Hpd + reg*speye(n));

                if p == 0
                    Hpd = Hpd + reg*speye(n);
                    return
                end

                reg = max(10*reg,1e-8);
            end

            Hpd = opts.InitialHessianScale*speye(n);
        end

        %--------------------------------------------------------------
        function hessian = objectiveGaussNewtonHessian(~,costFun,z,n)
        %OBJECTIVEGAUSSNEWTONHESSIAN Obtain the analytic least-squares curvature.
        % The cost callback owns this derivative because only it knows the
        % residual contract.  It must be rebuilt at every SQP major iterate.
            try
                [~,~,hessian] = costFun(z);
            catch ME
                error('SQPSolver:GaussNewtonHessian', ...
                    ['HessianMode "gauss_newton" requires costFun to return ', ...
                     '[cost,gradient,hessian]. %s'],ME.message);
            end
            assert(isequal(size(hessian),[n,n]) && ...
                all(isfinite(hessian),'all'), ...
                'SQPSolver:GaussNewtonHessian', ...
                'The analytic Gauss-Newton Hessian must be finite and n-by-n.');
            hessian = sparse(0.5*(hessian+hessian.'));
        end

        %--------------------------------------------------------------
        function hessian = dampGaussNewtonHessian(obj,hessian)
        %DAMPGAUSSNEWTONHESSIAN Apply optional Marquardt diagonal damping.
            damping = obj.options.GaussNewtonDamping;
            if damping == 0
                return
            end
            diagonal = max(abs(diag(hessian)),obj.options.InitialHessianScale);
            hessian = hessian+damping*spdiags(diagonal,0, ...
                numel(diagonal),numel(diagonal));
        end

        %--------------------------------------------------------------
        function report = checkGaussNewtonHessian(obj,costFun,z,lb,ub,hessian)
        %CHECKGAUSSNEWTOnHESSIAN Audit local curvature against grad differences.
        % A nonzero discrepancy is expected when nonlinear residual second
        % derivatives are material; this check records that approximation
        % rather than treating it as an analytic-gradient failure.
            directionCount = min(6,obj.options.GradientCheckDirections);
            step = obj.options.GradientCheckStep;
            generator = RandStream('mt19937ar','Seed',1801);
            relativeError = nan(directionCount,1);
            for directionIndex = 1:directionCount
                direction = randn(generator,numel(z),1);
                direction = direction/max(norm(direction),eps);
                h = obj.boundSafeStep(z,direction,step,lb,ub);
                if h <= 0
                    continue
                end
                [~,gradientPlus] = costFun(z+h*direction);
                [~,gradientMinus] = costFun(z-h*direction);
                finiteDifference = (gradientPlus(:)-gradientMinus(:))/(2*h);
                approximation = hessian*direction;
                relativeError(directionIndex) = norm( ...
                    finiteDifference-approximation,inf)/ ...
                    max([1,norm(finiteDifference,inf),norm(approximation,inf)]);
            end
            eigenvalues = eig(full(hessian));
            report = struct( ...
                'directions',directionCount, ...
                'step',step, ...
                'gradientDifferenceRelativeError',relativeError, ...
                'maximumGradientDifferenceRelativeError', ...
                    max(relativeError,[],"omitnan"), ...
                'minimumEigenvalue',min(eigenvalues), ...
                'maximumEigenvalue',max(eigenvalues), ...
                'conditionEstimate',max(abs(eigenvalues))/ ...
                    max(min(abs(eigenvalues)),eps));
        end

        %--------------------------------------------------------------
        function updateBFGS(obj,s,y)
            opts = obj.options;

            s = s(:);
            y = y(:);

            if norm(s) <= eps || norm(y) <= eps
                return
            end

            H = obj.H;
            H = 0.5*(H+H.');

            Hs = H*s;

            sHs = s.'*Hs;
            sy  = s.'*y;

            if sHs <= eps
                obj.resetHessian(numel(s));
                return
            end

            % Powell damping to preserve positive definiteness.
            if sy < 0.2*sHs
                theta = (0.8*sHs)/(sHs - sy);
                y = theta*y + (1-theta)*Hs;
                sy = s.'*y;
            end

            if sy <= eps
                return
            end

            Hnew = H - (Hs*Hs.')/sHs + (y*y.')/sy;
            Hnew = 0.5*(Hnew+Hnew.');

            if any(~isfinite(Hnew(:)))
                obj.resetHessian(numel(s));
            else
                obj.H = Hnew;
            end
        end

        %--------------------------------------------------------------
        function H = limitedMemoryBfgsHessian(obj,n)
        %LIMITEDMEMORYBFGSHESSIAN Reconstruct H from retained curvature pairs.
            H = obj.options.InitialHessianScale*speye(n);
            for pairIndex = 1:numel(obj.lbfgsSHistory)
                H = obj.dampedBfgsUpdate(H,obj.lbfgsSHistory{pairIndex}, ...
                    obj.lbfgsYHistory{pairIndex});
            end
        end

        %--------------------------------------------------------------
        function updateLimitedMemoryBfgs(obj,s,y,Hcurrent)
        %UPDATELIMITEDMEMORYBFGS Store one damped curvature pair.
            Hnext = obj.dampedBfgsUpdate(Hcurrent,s,y);
            if isempty(Hnext)
                return
            end
            s = s(:);
            Hs = Hcurrent*s;
            sHs = s.'*Hs;
            sy = s.'*y(:);
            if sHs <= eps || norm(s) <= eps || norm(y) <= eps
                return
            end
            if sy < 0.2*sHs
                theta = (0.8*sHs)/(sHs-sy);
                y = theta*y(:)+(1-theta)*Hs;
            else
                y = y(:);
            end
            if s.'*y <= eps
                return
            end
            obj.lbfgsSHistory{end+1,1} = s;
            obj.lbfgsYHistory{end+1,1} = y;
            maximumPairs = obj.options.LimitedMemoryBfgsPairs;
            if numel(obj.lbfgsSHistory) > maximumPairs
                obj.lbfgsSHistory = obj.lbfgsSHistory(end-maximumPairs+1:end);
                obj.lbfgsYHistory = obj.lbfgsYHistory(end-maximumPairs+1:end);
            end
            obj.H = Hnext;
        end

        %--------------------------------------------------------------
        function Hnew = dampedBfgsUpdate(~,H,s,y)
        %DAMPEDBFGSUPDATE Return one positive-curvature Powell-damped update.
            s = s(:); y = y(:);
            if norm(s) <= eps || norm(y) <= eps
                Hnew = [];
                return
            end
            H = 0.5*(H+H.');
            Hs = H*s;
            sHs = s.'*Hs;
            sy = s.'*y;
            if sHs <= eps
                Hnew = [];
                return
            end
            if sy < 0.2*sHs
                theta = (0.8*sHs)/(sHs-sy);
                y = theta*y+(1-theta)*Hs;
                sy = s.'*y;
            end
            if sy <= eps
                Hnew = [];
                return
            end
            Hnew = 0.5*(H-(Hs*Hs.')/sHs+(y*y.')/sy);
            if any(~isfinite(Hnew),'all')
                Hnew = [];
            end
        end

        %--------------------------------------------------------------
        % function gradLag = lagrangianGradient(~,g,gradc,gradceq,lambdaIneq,lambdaEq)
        %     gradLag = g(:);
        % 
        %     if ~isempty(gradceq) && ~isempty(lambdaEq)
        %         gradLag = gradLag + gradceq*lambdaEq(:);
        %     end
        % 
        %     if ~isempty(gradc) && ~isempty(lambdaIneq)
        %         gradLag = gradLag + gradc*lambdaIneq(:);
        %     end
        % end

        function gradLag = lagrangianGradient(~,g,gradc,gradceq,lambdaIneq,lambdaEq,lambdaLower,lambdaUpper)
        %LAGRANGIANGRADIENT Full KKT stationarity residual including bounds.
        %
        % Constraints:
        %   c(z) <= 0
        %   ceq(z) = 0
        %   lb <= z <= ub
        %
        % Bound convention:
        %   lower: lb - z <= 0  -> derivative -I
        %   upper: z - ub <= 0  -> derivative +I
        %
        % Therefore:
        %   grad L = grad f + gradc*lambdaIneq + gradceq*lambdaEq
        %            - lambdaLower + lambdaUpper
        
            gradLag = g(:);
        
            if nargin < 7 || isempty(lambdaLower)
                lambdaLower = zeros(size(gradLag));
            end
        
            if nargin < 8 || isempty(lambdaUpper)
                lambdaUpper = zeros(size(gradLag));
            end
        
            if ~isempty(gradceq) && ~isempty(lambdaEq)
                gradLag = gradLag + gradceq*lambdaEq(:);
            end
        
            if ~isempty(gradc) && ~isempty(lambdaIneq)
                gradLag = gradLag + gradc*lambdaIneq(:);
            end
        
            gradLag = gradLag - lambdaLower(:) + lambdaUpper(:);
        end

        %--------------------------------------------------------------
        function feas = constraintViolation(~,c,ceq,z,lb,ub)
            vals = 0;

            if ~isempty(ceq)
                vals = max(vals,norm(ceq,inf));
            end

            if ~isempty(c)
                vals = max(vals,max([0;c(:)]));
            end

            vals = max(vals,max([0;lb(:)-z(:)]));
            vals = max(vals,max([0;z(:)-ub(:)]));

            feas = vals;
        end

        %--------------------------------------------------------------
        function phi = merit(obj,costFun,nonlFun,z,lb,ub,rho)
            f = costFun(z);

            [c,ceq] = nonlFun(z);

            phi = obj.meritFromValues(f,c,ceq,z,lb,ub,rho);
        end

        %--------------------------------------------------------------
        function phi = meritFromValues(~,f,c,ceq,z,lb,ub,rho)

            v = 0;

            if ~isempty(ceq)
                v = v + norm(ceq(:),1);
            end

            if ~isempty(c)
                v = v + sum(max(0,c(:)));
            end

            v = v + sum(max(0,lb(:)-z(:)));
            v = v + sum(max(0,z(:)-ub(:)));

            phi = f + rho*v;
        end

        %--------------------------------------------------------------
        function [f,c,ceq] = evalValues(~,costFun,nonlFun,z)
            f = costFun(z);
            [c,ceq] = nonlFun(z);
            c = c(:);
            ceq = ceq(:);
        end

        %--------------------------------------------------------------
        function theta = l1InfeasibilityFromValues(~,c,ceq,z,lb,ub)
            theta = 0;
            if ~isempty(ceq)
                theta = theta + norm(ceq(:),1);
            end
            if ~isempty(c)
                theta = theta + sum(max(0,c(:)));
            end
            theta = theta + sum(max(0,lb(:)-z(:)));
            theta = theta + sum(max(0,z(:)-ub(:)));
        end

        %--------------------------------------------------------------
        function dTheta = l1InfeasibilityDirectional( ...
                ~,c,ceq,z,lb,ub,gradc,gradceq,d)
            dTheta = 0;

            if ~isempty(ceq)
                dceq = gradceq.'*d;
                positive = ceq > 0;
                negative = ceq < 0;
                zero = ~(positive | negative);
                dTheta = dTheta + sum(dceq(positive)) - ...
                    sum(dceq(negative)) + sum(abs(dceq(zero)));
            end

            if ~isempty(c)
                dc = gradc.'*d;
                positive = c > 0;
                zero = c == 0;
                dTheta = dTheta + sum(dc(positive)) + ...
                    sum(max(0,dc(zero)));
            end

            lowerResidual = lb-z;
            lowerDerivative = -d;
            positive = lowerResidual > 0;
            zero = lowerResidual == 0;
            dTheta = dTheta + sum(lowerDerivative(positive)) + ...
                sum(max(0,lowerDerivative(zero)));

            upperResidual = z-ub;
            upperDerivative = d;
            positive = upperResidual > 0;
            zero = upperResidual == 0;
            dTheta = dTheta + sum(upperDerivative(positive)) + ...
                sum(max(0,upperDerivative(zero)));
        end

        %--------------------------------------------------------------
        function [zCorrected,info] = secondOrderCorrectTrial(obj, ...
                nonlFun,zCurrent,zTrial,~,~,lb,ub)
        %SECONDORDERCORRECTTRIAL Recover nonlinear feasibility at a trial.
        % The correction is a bounded minimum-norm solve of the trial-point
        % constraint linearization.  It is applied only in the explicit
        % audit mode and only when it reduces L1 infeasibility.

            zCorrected = zTrial;
            info = struct('applied',false,'correctionNormInf',0, ...
                'rawInfeasibilityL1',nan);
            if ~obj.options.SecondOrderCorrectionEnabled
                return
            end

            try
                [cRaw,ceqRaw,gradcRaw,gradceqRaw] = nonlFun(zTrial);
                cRaw = cRaw(:);
                ceqRaw = ceqRaw(:);
                n = numel(zTrial);
                if isempty(gradcRaw), gradcRaw = zeros(n,0); end
                if isempty(gradceqRaw), gradceqRaw = zeros(n,0); end
                if size(gradcRaw,1) ~= n && size(gradcRaw,2) == n
                    gradcRaw = gradcRaw.';
                end
                if size(gradceqRaw,1) ~= n && size(gradceqRaw,2) == n
                    gradceqRaw = gradceqRaw.';
                end
            catch
                return
            end

            rawTheta = obj.l1InfeasibilityFromValues( ...
                cRaw,ceqRaw,zTrial,lb,ub);
            info.rawInfeasibilityL1 = rawTheta;
            if rawTheta <= obj.options.ConstraintTolerance
                return
            end

            activeIneq = cRaw > 0;
            if isempty(ceqRaw) && ~any(activeIneq)
                return
            end
            if size(gradceqRaw,1) ~= n || size(gradcRaw,1) ~= n
                return
            end

            if any(activeIneq)
                Aineq = gradcRaw(:,activeIneq).';
                bineq = -cRaw(activeIneq);
            else
                Aineq = [];
                bineq = [];
            end
            if isempty(ceqRaw)
                Aeq = [];
                beq = [];
            else
                Aeq = gradceqRaw.';
                beq = -ceqRaw;
            end

            lower = lb-zTrial;
            upper = ub-zTrial;
            try
                [correction,~,socExitflag] = quadprog(speye(n),zeros(n,1), ...
                    Aineq,bineq,Aeq,beq,lower,upper,[], ...
                    obj.options.QPOptions);
            catch
                return
            end
            if socExitflag <= 0 || isempty(correction) || ...
                    any(~isfinite(correction))
                return
            end

            correctionNorm = norm(correction,inf);
            primaryNorm = max(norm(zTrial-zCurrent,inf),eps);
            if correctionNorm > ...
                    obj.options.SecondOrderCorrectionMaxRelativeNorm*primaryNorm
                return
            end
            candidate = min(max(zTrial+correction,lb),ub);
            try
                [cCandidate,ceqCandidate] = nonlFun(candidate);
            catch
                return
            end
            candidateTheta = obj.l1InfeasibilityFromValues( ...
                cCandidate,ceqCandidate,candidate,lb,ub);
            if candidateTheta > rawTheta + ...
                    10*eps(max(1,rawTheta))
                return
            end

            zCorrected = candidate;
            info.applied = true;
            info.correctionNormInf = correctionNorm;
        end

        %--------------------------------------------------------------
        function alphaNext = nextBacktrackingAlpha(~,alpha,phi0,phiTrial, ...
                directionalMerit,opts)
        %NEXTBACKTRACKINGALPHA Safeguarded merit interpolation update.
        % The quadratic model is used only for a finite descent direction;
        % otherwise the established fixed-beta contraction is retained.

            alphaNext = opts.LineSearchBeta*alpha;
            if opts.LineSearchUpdate ~= "quadratic_safeguarded" || ...
                    ~isfinite(phi0) || ~isfinite(phiTrial) || ...
                    ~isfinite(directionalMerit) || directionalMerit >= 0
                return
            end
            denominator = 2*(phiTrial-phi0-directionalMerit*alpha);
            if ~isfinite(denominator) || denominator <= 0
                return
            end
            quadraticTrial = -directionalMerit*alpha^2/denominator;
            lowerSafeguard = 0.1*alpha;
            upperSafeguard = 0.8*alpha;
            alphaNext = min(max(quadraticTrial,lowerSafeguard),upperSafeguard);
        end

        %--------------------------------------------------------------
        function decision = filterTrialDecision(obj,f,theta,directionalObj, ...
                alpha,fTrial,thetaTrial,thetaMinimum,thetaMaximum, ...
                filterTheta,filterObjective)
            opts = obj.options;
            decision = obj.emptyFilterDecision();

            decision.filterAcceptable = thetaTrial < thetaMaximum && ...
                all(thetaTrial < filterTheta | fTrial < filterObjective);
            if ~decision.filterAcceptable
                return
            end

            decision.switching = theta <= thetaMinimum && ...
                directionalObj < 0 && ...
                alpha*(-directionalObj)^opts.FilterSObjective > ...
                    opts.FilterDelta*theta^opts.FilterSTheta;
            decision.armijoObjective = fTrial <= ...
                f + opts.FilterEtaObjective*alpha*directionalObj;
            decision.sufficientFeasibility = thetaTrial <= ...
                (1-opts.FilterGammaTheta)*theta;
            decision.sufficientObjective = fTrial <= ...
                f-opts.FilterGammaObjective*theta;

            if decision.switching
                decision.accepted = decision.armijoObjective;
            else
                decision.accepted = decision.sufficientFeasibility || ...
                    decision.sufficientObjective;
            end
        end

        %--------------------------------------------------------------
        function decision = emptyFilterDecision(~)
            decision = struct( ...
                'accepted',false, ...
                'filterAcceptable',false, ...
                'switching',false, ...
                'armijoObjective',false, ...
                'sufficientFeasibility',false, ...
                'sufficientObjective',false);
        end

        %--------------------------------------------------------------
        function validateOpenUnitOption(~,value,name)
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
                    value > 0 && value < 1)
                error('SQPSolver:OptionRange', ...
                    '%s must be a finite scalar strictly between 0 and 1.',name);
            end
        end

        %--------------------------------------------------------------
        function validatePositiveOption(~,value,name)
            if ~(isnumeric(value) && isscalar(value) && isfinite(value) && ...
                    value > 0)
                error('SQPSolver:OptionRange', ...
                    '%s must be a finite positive scalar.',name);
            end
        end

        %--------------------------------------------------------------
        function h = boundSafeStep(~,z,d,h0,lb,ub)
            h = h0;

            active = abs(d) > 0;

            if any(active)
                distUpper = inf(size(z));
                distLower = inf(size(z));

                idx = active & isfinite(ub);
                distUpper(idx) = (ub(idx)-z(idx))./abs(d(idx));

                idx = active & isfinite(lb);
                distLower(idx) = (z(idx)-lb(idx))./abs(d(idx));

                hBound = min([distUpper(active);distLower(active)]);

                if isfinite(hBound)
                    h = min(h,0.49*hBound);
                end
            end

            if ~isfinite(h) || h <= 0
                h = 0;
            end
        end

        %--------------------------------------------------------------
        function output = makeOutput(~,iter,funcCount,constrCount,fval,feas,kkt,step,message, ...
                                    qpExitflag,qpOutput,alpha,slackEqInf,slackIneqInf)
            output = struct();

            output.iterations = iter;
            output.funcCount = funcCount;
            output.constrCount = constrCount;
            output.fval = fval;
            output.constrviolation = feas;
            output.firstorderopt = kkt;
            output.stepsize = step;
            output.message = char(message);
            output.algorithm = 'custom-sqp-bfgs-quadprog';
            output.qpExitflag = qpExitflag;
            output.qpOutput = qpOutput;
            output.alpha = alpha;
            output.slackEqInf = slackEqInf;
            output.slackIneqInf = slackIneqInf;
        end

        %--------------------------------------------------------------
        function printHeader(~)
            fprintf('\n');
            fprintf('  iter |      fval      |  feas_inf  |   opt_inf  |  step_inf  | alpha | slackEq | slackIneq | qp\n');
            fprintf('-------------------------------------------------------------------------------------------------\n');
        end

        %--------------------------------------------------------------
        function plotData = createLiveIterationPlot(~,enabled)
        %CREATELIVEITERATIONPLOT Default-disabled, audit-only SQP dashboard.
            plotData = struct('enabled',false);
            if ~enabled
                return
            end
            figureHandle = figure('Name','SQP iteration diagnostics', ...
                'NumberTitle','off','Color','w');
            layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
            titles = ["Objective" "Feasibility infinity norm" ...
                "First-order residual infinity norm" "Step size and alpha"];
            lines = gobjects(5,1);
            for index = 1:4
                nexttile(layout,index); hold on; grid on;
                title(titles(index)); xlabel('SQP iteration');
                set(gca,'YScale','log');
                ylim([1e-10 inf]);
                if index < 4
                    lines(index) = animatedline('LineWidth',1.2, ...
                        'Marker','o','MarkerSize',4);
                else
                    lines(4) = animatedline('LineWidth',1.2, ...
                        'Marker','o','MarkerSize',4,'DisplayName','step');
                    lines(5) = animatedline('LineWidth',1.2, ...
                        'Marker','x','MarkerSize',4,'DisplayName','alpha');
                    legend('Location','best');
                end
            end
            title(layout,'Live SQP convergence (audit diagnostic only)');
            plotData = struct('enabled',true,'figure',figureHandle, ...
                'lines',lines);
        end

        %--------------------------------------------------------------
        function updateLiveIterationPlot(~,plotData,iteration,objective, ...
                feasibility,optimality,step,alpha)
            if ~plotData.enabled || ~isgraphics(plotData.figure)
                return
            end
            plotFloor = 1e-10;
            addpoints(plotData.lines(1),iteration,max(abs(objective),plotFloor));
            addpoints(plotData.lines(2),iteration,max(abs(feasibility),plotFloor));
            addpoints(plotData.lines(3),iteration,max(abs(optimality),plotFloor));
            addpoints(plotData.lines(4),iteration,max(abs(step),plotFloor));
            addpoints(plotData.lines(5),iteration,max(abs(alpha),plotFloor));
            drawnow limitrate
        end

        %--------------------------------------------------------------
        function closeLiveIterationPlot(~,plotData)
            if isstruct(plotData) && isfield(plotData,'enabled') && ...
                    plotData.enabled && isgraphics(plotData.figure)
                close(plotData.figure);
            end
        end

        %--------------------------------------------------------------
        function printIter(~,iter,f,feas,kkt,step,alpha,slackEq,slackIneq,qpflag)
            fprintf('  %4d | %14.6e | %9.2e | %9.2e | %9.2e | %5.2f | %7.1e | %9.1e | %3g\n', ...
                iter,f,feas,kkt,step,alpha,slackEq,slackIneq,qpflag);
        end
    end
end
