%%
% control/@nMPC/nMPC.m
classdef nMPCold < ControllerBase
% NMPC with multiple shooting & disturbance ramp‑down, following
% Artola et al. AIAA‑J 2021, Eq.(37)  :contentReference[oaicite:7]{index=7}

    properties (Access = private)
        f
        dt, Nc
        Qc, Pc, Rc
        uL, uU, xL, xU
        wRamp                    % predicted w profile
        sensFcn
        opts
        lastSolve
    end

    methods
        function obj = nMPCold(cfg)
            obj@ControllerBase(cfg);
            obj.f  = cfg.f;  obj.dt = cfg.dt;  obj.Nc = cfg.Nc;
            obj.Qc = cfg.Qc; obj.Pc = cfg.Pc;  obj.Rc = cfg.Rc;
            obj.uL = cfg.uL; obj.uU = cfg.uU;
            obj.xL = cfg.xL; obj.xU = cfg.xU;
            obj.wRamp = @(w0) [linspace(w0,0,obj.Nc/2), zeros(1,obj.Nc/2)];
            obj.sensFcn = @(x0,u) aeroelasticAdjointSens(obj.f,x0,u,0,obj.dt);
            obj.opts = optimoptions('fmincon','Algorithm','sqp', ...
                       'Display','none','MaxIterations',cfg.maxIter);
        end

        function u = computeControl(obj,xhat,t)
            
            % % 1) build the prediction NLP --------------------------------
            % nlp  = MultipleShootingNLP( buildMPCcfg(obj,xhat) );
            %    %  • M  = obj.cfg.Nc
            %    %  • f  = same modal RHS handle
            %    %  • L  = stage cost using Qc,Rc
            %    %  • Jterm = terminal cost with Pc
            %    %  • bounds on x,u
            % 
            % % 2) solve ----------------------------------------------------
            % [p,info] = obj.sqp.solve(nlp,xhat,obj.p0);
            % 
            % % 3) first control in optimal sequence -----------------------
            % m   = obj.cfg.nu;
            % u   = p( (obj.nx*(obj.M-1))+ (1:m) );
            % obj.p0           = p;    % warm start
            % obj.lastSolveInfo= info; % for log/debug
            % obj.prevU        = u;

            % decision vector p = [u0 … uNc‑1]
            m  = numel(obj.prevU);
            p0 = repmat(obj.prevU,m*obj.Nc,1);

            lb = repmat(obj.uL,obj.Nc,1);
            ub = repmat(obj.uU,obj.Nc,1);

            cost = @(p) obj.costFun(p,xhat);
            nonl = @(p) obj.nonlinCon(p,xhat);
            
            % % inside computeControl() of nMPC_Artola
            % % -------------------------------------
            % nlp = MultipleShootingNLP(buildMPCcfg(obj,xhat));   % helper returns struct
            % [p,info] = obj.sqp.solve(nlp,xhat,p0);              %#ok<NASGU>
            % u = obj.saturate( p(1:obj.cfg.nu) );
            % obj.sqp = SQPSolver(obj.cfg.sqpOpts);



            p  = fmincon(cost,p0,[],[],[],[],lb,ub,nonl,obj.opts);
            obj.lastSolve = p;
            u = p(1:m);
            u = obj.saturate(u);
            obj.prevU = u;
        end
    %---------------------------------------------------------------
        function J = costFun(obj,p,x0)
            m = numel(obj.uL);
            uSeq = reshape(p,m,[]);
            wSeq = obj.wRamp(0);           % vanishing disturbance assumption
            xk = x0;  J = 0;
            for k = 1:obj.Nc
                uk = uSeq(:,k);    wk = wSeq(k);
                ek = xk - obj.cfg.xRef;
                J  = J + 0.5*ek.'*obj.Qc*ek + 0.5*uk.'*obj.Rc*uk;
                xk = obj.intStep(xk,uk,wk);
            end
            J = J + 0.5*(xk-obj.cfg.xRef).'*obj.Pc*(xk-obj.cfg.xRef);
        end
    %---------------------------------------------------------------
        function [c,ceq] = nonlinCon(obj,p,x0)
            m = numel(obj.uL);
            uSeq = reshape(p,m,[]);
            wSeq = obj.wRamp(0);
            xk = x0;    ceq = [];
            for k = 1:obj.Nc
                xNext = obj.intStep(xk,uSeq(:,k),wSeq(k));
                ceq   = [ceq ; xNext - obj.guard(k)]; %#ok<AGROW>
                xk = xNext;
            end
            c = [];   % box handled
        end
    %---------------------------------------------------------------
        function xn = intStep(obj,x,u,w)
            xn = x + obj.f(x,u,w,0)*obj.dt;
        end
        function xg = guard(~,~), xg = zeros(0,1); end
    end
end






% classdef nMPC < ControllerBase
% %NMPC  Nonlinear model‑predictive controller (multiple shooting SQP).
% %
% %   Requires an external optimisation backend passed via cfg.solver.
% %   Only the call signature is enforced here; heavy lifting is delegated
% %   to the user‑supplied solver function for real‑time guarantees.
% 
%     properties (Access = private)
%         solver   % function handle   u = solver(xhat,t,prevU,cfg)
%     end
% 
%     methods
%         function obj = nMPC(cfg)
%             obj@ControllerBase(cfg);
%             obj.solver = cfg.solver;
%         end
% 
%         function u = computeControl(obj,xhat,t)
%             u = obj.solver(xhat,t,obj.prevU,obj.cfg);
%             u = obj.saturate(u);
%             obj.prevU = u;
%         end
%     end
% end
