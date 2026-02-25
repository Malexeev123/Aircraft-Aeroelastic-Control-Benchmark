


%%
% control/@nMHE/nMHE.m
classdef nMHE2 < AeroFlex.ctrl.EstimatorBase
    
    
% Non‑linear MHE as in Artola 2021, Eq.(38) – multiple shooting + SQP
% References: Artola et al., AIAA J‑59‑7, pp.2697‑2713, 2021  

    properties % (Access = private)
        f  % dynamics  handle   xdot = f(x,u,w,t)
        h  % output    handle   y    = h(x)
        dt
        Ne              % window length (# samples)
        Qe, Pe, Re        % weight matrices
        wL, wU, xL, xU     % bounds
        solverOpts
        sqp
        cfg
        % --- circular buffers ---
        Y , U,  T
        idx
        % --- sensitivities cache ---
        sensFcn
    end

    methods
        function obj = nMHE2(cfg)
            obj@AeroFlex.ctrl.EstimatorBase(cfg);
            % obj.sqp = SQPSolver(cfg.sqpOpts);
            % % obj.nx  = numel(obj.xhat);
            % % obj.p0  = zeros(obj.nx + cfg.nu*cfg.Ne + cfg.Ne,1); % warm start

            obj.cfg.trim = cfg.trim;
            obj.f  = cfg.f;   obj.h = cfg.h;
            obj.dt = cfg.dt;  obj.Ne = cfg.Ne;
            obj.Qe = cfg.Qe;  obj.Pe = cfg.Pe; obj.Re = cfg.Re;
            obj.wL = cfg.wL;  obj.wU = cfg.wU;
            obj.xL = cfg.xL;  obj.xU = cfg.xU;
            obj.solverOpts = optimoptions('fmincon','Algorithm','sqp', ...
                              'MaxIterations',cfg.maxIter,'Display','none');
            obj.Y = nan(cfg.ny,obj.Ne);  obj.U = nan(cfg.nu,obj.Ne);
            obj.T = nan(1,obj.Ne);       obj.idx = 0;

            % sensitivity function (analytical adjoint) – Eq.(53) :contentReference[oaicite:3]{index=3}
            obj.sensFcn = @(x0,u,w)        ...
                 aeroelasticAdjointSens(obj.f,x0,u,w,obj.dt);
        end

        function xhat = estimate(obj,z,u,t)
            % % 1) build Multiple‑Shooting NLP ------------------------------
            nlp = MultipleShootingNLP( buildMHEcfg(obj,zk,uk) );
            %    %  • buildMHEcfg fills: f, L (stage cost), Jterm, bounds
            %    %  • M  = obj.cfg.Ne       (window length)
            %    %  • dt = obj.cfg.dt
            % 
            % % 2) solve with generic SQP engine ----------------------------
            [p,info]    = obj.sqp.solve(nlp, obj.xhat, obj.p0);
            % 
            % % 3) unpack solution -----------------------------------------
            % obj.xhat    = p(1:obj.nx);         % initial state of window
            % obj.yhat    = zk - info.innov;     % last predicted output
            % obj.p0      = p;                   % warm‑start next call
            % xhat        = obj.xhat;

            obj.idx = obj.idx + 1;
            k = mod(obj.idx-1,obj.Ne)+1;
            obj.Y(:,k) = z;  obj.U(:,k) = u;  obj.T(k)=t;

            if obj.idx < obj.Ne                                       %#ok<*BDSCI>
                obj.yhat = obj.h(obj.xhat);
                xhat = obj.xhat;  return
            end

            % ----- decision vector p = [xE ; wHist] ------------------
            n  = numel(obj.xhat);           m = numel(u);
            p0 = [obj.xhat ; zeros(m*obj.Ne,1)];

            lb = [obj.xL ; repmat(obj.wL,obj.Ne,1)];
            ub = [obj.xU ; repmat(obj.wU,obj.Ne,1)];

            cost = @(p) obj.costFun(p);
            nonl = @(p) obj.nonlinCon(p);

            p  = fmincon(cost,p0,[],[],[],[],lb,ub,nonl,obj.solverOpts);

            obj.xhat = p(1:n);
            obj.yhat = obj.h(obj.xhat);
            xhat     = obj.xhat;
        end
    %---------------------------------------------------------------
        function J = costFun(obj,p)
            n = numel(obj.xhat);   m = numel(obj.wL);
            xE = p(1:n);
            wHist = reshape(p(n+1:end),m,[]);

            % first two quadratic penalties – Eq.(38a) :contentReference[oaicite:4]{index=4}
            J  = 0.5*(xE-obj.cfg.ref.states).' * obj.Pe * (xE-obj.cfg.ref.states);
            J  = J + 0.5*(wHist(:,1)-obj.cfg.ref.w_r).' ...
                        * obj.Re * (wHist(:,1)-obj.cfg.ref.w_r);

            % output mismatch Σ_{k=-Ne}^0 0.5 * e^T Qe e
            xk = xE;   wk = wHist(:,1);
            for j = obj.Ne:-1:1
                uk = obj.U(:,j);
                zk = obj.Y(:,j);
                yk = obj.h(xk);
                J  = J + 0.5*(zk-yk).'*obj.Qe*(zk-yk);

                % roll window backwards through dynamics
                xk = obj.intStepBack(xk,uk,wk);
                if j>1, wk = wHist(:,j); end
            end
        end
    %---------------------------------------------------------------
        function [c,ceq] = nonlinCon(obj,p)
            % Multiple shooting continuity constraints – Eq.(38c) :contentReference[oaicite:5]{index=5}
            n = numel(obj.xhat);   m = numel(obj.wL);
            xE = p(1:n);
            wHist = reshape(p(n+1:end),m,[]);
            ceq = [];

            xk = xE;   wk = wHist(:,1);
            for j = obj.Ne:-1:1
                uk = obj.U(:,j);
                xkPrev = obj.intStepBack(xk,uk,wk);
                ceq = [ceq ; xkPrev - obj.Y(1:n,j)]; %#ok<AGROW>
                if j>1, wk = wHist(:,j); end
                xk = xkPrev;
            end
            c = [];  % box handled by fmincon bounds
        end
    %---------------------------------------------------------------
        function xprev = intStepBack(obj,x,u,w)
            % reverse‑time integration over dt (Euler backward)
            xprev = x - obj.f(x,u,w,0)*obj.dt;
        end
    end
end
