classdef SQPSolver < handle
%SQPSolver  Lightweight SQP with BFGS Hessian & back‑tracking
%
%   sol = SQPSolver(opts).solve(problem,xInit);
%   • problem must expose  objFun(p,x1) and nonlin(p,x1)
%   • Supports analytic gradients if problem.gradFun / jacFun set
%
%   Options (struct):
%     maxIter, tol,  stepRule ('armijo'|'filter'),  qpBackend (handle)
%   opts = struct('maxIter',20,'tol',1e-6,'verbose',true,'qpBackend',@quadprog)


    properties, opts,   H,  qp,  end
    methods
        function obj = SQPSolver(opts)
            obj.opts = opts;
            if ~isfield(opts,'qpBackend') 
                obj.qp = @quadprog; 
            else 
                obj.qp = opts.qpBackend; 
            end
        end

        function [p,info] = solve(obj,prob,x1,p0)
            n = numel(p0);   obj.H = eye(n);
            p = p0;  
            info = struct;
            info.iter=0;
            info.conv=false;

            [lb,ub] = prob.bounds;

            for k = 1:obj.opts.maxIter
                info.iter = k;
                [F,g] = prob.objFun(p,x1);
                [~,ceq,~,Aeq] = prob.nonlin(p,x1);
                % Aeq = prob.jacFun(p,x1) if ~isempty(prob.jacFun); %#ok<*NOPRT> 
                % --- QP sub‑problem ----------------------------------
                d = obj.qp(obj.H, g, [],[], Aeq, -ceq, lb-p, ub-p, [], ...
                           optimoptions('quadprog'));
                           % optimoptions('quadprog','Display','off'));

                % --- line‑search / filter ----------------------------
                alpha = obj.lineSearch(prob,p,x1,d,F,g);

                pNew = p + alpha*d;
                % if norm(d)<obj.opts.tol, p=pNew; break; end

                % --- BFGS update -------------------------------------
                s = pNew - p;
                [~,gNew] = prob.objFun(pNew,x1);
                y = gNew - g;
                if s'*y > 1e-8
                    rho = 1/(y'*s);
                    obj.H = (eye(n)-rho*s*y')*obj.H*(eye(n)-rho*y*s') + rho*(s*s');
                end

                if norm(d) < obj.opts.tol
                    info.conv=true; info.iter=k; break
                end
                p = pNew;
            end
            info.finalCost = prob.objFun(p,x1);
        end

        function alpha = lineSearch(obj,prob,p,x1,d,F,g)
            alpha = 1;
            while alpha>1e-4
                pTry = p + alpha*d;
                [Ftry,~] = prob.objFun(pTry,x1);
                if Ftry < F + 1e-4*alpha*g'*d, break, end
                alpha = 0.5*alpha;
            end
        end
    end
end
