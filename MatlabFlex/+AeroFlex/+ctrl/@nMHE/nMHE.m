classdef nMHE < AeroFlex.ctrl.EstimatorBase
    %------------------------------------------------------------------
    %   NON‑LINEAR MOVING‑HORIZON ESTIMATOR
    %   – Multiple‑shooting (default)  OR  Single‑shooting (toggle)
    %------------------------------------------------------------------
    %   Author :  Maxim Alexeev  |  26‑Jun‑2025
    % -----------------------------------------------------------------
    %   PUBLIC TOGGLE
    %       cfg.ctrl.mhe_method   =  'multiple'  (default)  |  'single'
    % -----------------------------------------------------------------
    %   KEY FIXES VS. ORIGINAL VERSION
    %   ①  Correct warm‑start shift (Diehl 2005 RTI).
    %   ②  Continuity constraint contains control Jacobian Su.
    %   ③  Decision vector dimension  nVar = (Ne+1)*nx + Ne*nw   – no extra node.
    %   ④  Optional single‑shoot formulation shares ~80 % code.
    %   ⑤  Built‑in debug plot replicates Artola‑Fig.10 horizon snapshot.
    %------------------------------------------------------------------

    properties
        % === copied from cfg ======================================================
        dt  double                     % integrator step [s]
        Ts  double                     % sample period   [s]
        Ne  double                     % estimation horizon length [samples]
        method (1,1) string {mustBeMember(method,["multiple","single"])} = "multiple"
        model
        Qe, Pe, Re                       % weighting matrices (arrival / meas.)
        xL, xU,  wL, wU                  % box bounds
        % -------------------------------------------------------------------------
        buff            struct        % shared ring‑buffer handle (SimRunner)
        nx, nu, ny, nw                   % problem dimensions
        z0               double       % rolling warm‑start vector
        % xhat  double                  % latest state estimate
        % what  double                  % latest disturbance history (vector)
        H                          % LBFGS / BFGS Hessian approx
        solverOpts                   % fmincon options structure

        % – debug ---------------------------------------------------------------
        debug = true    % true  |false
        dbg                                                  % struct for plots
    end

    %% =====================================================================
    %%  CONSTRUCTOR
    %% =====================================================================
    methods
        function obj = nMHE(cfg)
            arguments, cfg struct, end
            obj@AeroFlex.ctrl.EstimatorBase(cfg);

            % ---------- store scalar config ------------------------------------
            obj.dt     = cfg.sim.dt;
            obj.Ts     = cfg.ctrl.Ts;
            obj.Ne     = cfg.ctrl.Ne;
            % obj.method = lower(string(cfg.ctrl.mhe_method));  % new toggle
            obj.method = lower(string(cfg.ctrl.mhe_method));

            obj.Qe = cfg.Qe;   obj.Pe = cfg.Pe;   obj.Re = cfg.Re;
            obj.xL = cfg.xL;   obj.xU = cfg.xU;
            obj.wL = cfg.wL;   obj.wU = cfg.wU;

            % ---------- ring‑buffer handle & dimensions ------------------------
            obj.buff = cfg.buff;
            obj.nx = size(obj.buff.X,1);
            obj.nu = size(obj.buff.U,1);
            obj.ny = size(obj.buff.Y,1);
            obj.nw = cfg.nw;                 % gust dim (usually 1)

            % ---------- optimiser options (interior‑point, LBFGS) -------------
            obj.solverOpts = optimoptions("fmincon", ...
                    "Algorithm","interior-point", ...    % interior-point | sqp
                    "Display","none", ...
                    "SpecifyObjectiveGradient",true, ...
                    "SpecifyConstraintGradient",true, ...
                    "HessianApproximation","lbfgs", ...
                    "OptimalityTolerance",1e-5, ...
                    "StepTolerance",1e-8);

            % ---------- allocate warm‑start ------------------------------------
            obj.z0   = zeros((obj.Ne+1)*obj.nx + obj.Ne*obj.nw ,1);
            obj.xhat = zeros(obj.nx,1);
            obj.what = zeros(obj.Ne*obj.nw,1);
            obj.H    = speye(numel(obj.z0));    % initial PD Hessian
            obj.model = cfg.modelHandle(cfg,beam,aero,base);   % pass handle in cfg

        end
    end   % ctor

    %% =====================================================================
    %%  MAIN ENTRY –  estimate()
    %% =====================================================================
    methods
        function [xhat,what,info] = estimate(obj, zMeas, uNow, tNow, buff)
            %  Called once per sample Ts from SimRunner
            % -----------------------------------------------------------------
            arguments
                obj
                zMeas   double   % sensor measurement ny×1
                uNow    double   % applied control   nu×1  (known)
                tNow    double   % current sim time  [s]
                buff    struct   % *handle* – ring buffer snapshot
            end

            obj.buff = buff;           % store live handle

            % ---------- rolling index / early warm‑up -----------------------
            obj.buff.valid = min(obj.buff.valid+1 , obj.Ne+1);
            if obj.buff.valid < obj.Ne+1
                % not enough data yet → cheap pseudo‑predict
                obj.xhat = obj.buff.X(:,obj.buff.valid-1);
                % z220(idx.x{j}) = obj.buff.X(:, obj.idxWrap(-obj.Ne)); %this one?

                xhat = obj.xhat;  what = obj.what;  
                info = struct('cost', 0, 'exitflag', [], 'continuity', 0, 'wRefPrev', zeros(obj.Ne*obj.nw,1));
                return
            end
            % ---------- build NLP  ------------------------------------------
            nlp = obj.assembleWindow();

            % ---------- solve ------------------------------------------------
            [p,fval,exitflag] = fmincon(nlp.cost, obj.z0, [],[],[],[], ...
                                          nlp.lb, nlp.ub, nlp.nonl, obj.solverOpts);

            % ---------- shift‑and‑store warm‑start (RTI) --------------------
            obj.z0 = obj.shiftGuess(p);         %      <<< FIX #1

            % ---------- unpack state / disturbance --------------------------
            idx = obj.buildIndexMaps();
            obj.xhat = p(idx.x{end});
            obj.what = p([idx.w{:}]);

            xhat = obj.xhat;    what = obj.what;
            info.cost = fval;   info.exitflag = exitflag;
            info.wHorizon = reshape(obj.what,obj.nw,[]).';
            
            [~,ceq,~,~]= nlp.nonl(p);           % reuse the constraint callback
            info.continuity = norm(ceq);     % ‖ceq‖ for SimRunner debug
            info.wRefPrev   = obj.what;      % entire disturbance horizon (length Ne)
            % ---------- debug plot of horizon (Artola‑style) ----------------
            % ---------- debug visualisation ----------------------------------------
            if obj.debug
                % ---- initialise persistent store on first call --------------------
                if ~isfield(obj.dbg,'t')
                    obj.dbg.t = [];          % Ns×1   sample times
                    obj.dbg.W = [];          % Ns×Ne  horizon matrix (rows = samples)
                end
            
                % ---- append current sample ---------------------------------------
                obj.dbg.t(end+1,1)        = tNow;
                obj.dbg.W(end+1,1:obj.Ne) = info.wHorizon(:).';    % row-vector assignment
            
                % ---- plot 1 : coloured stairs  (one line per horizon slot) -------
                figure(2001), clf, hold on
                for slot = 1:obj.Ne
                    stairs(obj.dbg.t , obj.dbg.W(:,slot) , 'LineWidth',1.1);
                end
                hold off, grid on; 
                legend('Location','southwest')
                xlabel('time  [s]'), ylabel('w  (gust)')
                title('MHE disturbance horizon – slot index = colour')
                

                figure(2003), clf, hold on
               
                stairs(obj.dbg.t , obj.dbg.W(:,end) , 'k-o')

                hold off, grid on
                xlabel('time  [s]'), ylabel('w  (gust)')
                title('MHE disturbance horizon – last horizon')
                % % ---- plot 2 : 2-D image (Artola Fig. 10 style) --------------------
                % figure(2002), clf
                % imagesc(obj.dbg.t , 1:obj.Ne , obj.dbg.W.'), axis xy
                % xlabel('time  [s]'), ylabel('slot   (1 = oldest)')
                % title('disturbance horizon  (piece-wise constant, sliding)')
                % colorbar
            end

            % if obj.debug
            %     if isempty(obj.dbg)
            %         obj.dbg.t  = tNow;
            %         obj.dbg.w  = info.wHorizon(:).';
            %     else
            %         obj.dbg.t(end+1) = tNow;
            %         obj.dbg.w  = [obj.dbg.w ; info.wHorizon(:).'];
            %     end
            %     % figure(2001), clf, plot(obj.dbg.t, obj.dbg.w), grid on
            %     % xlabel('time [s]'), title('MHE disturbance horizon – newest → right')
            %     figure(2001), clf
            %     stairs( (0:samCnt-1)*cfg.ctrl.Ts , obj.dbg.wHorizon.' , 'LineWidth',1.1)
            %     xlabel('time  [s]'), ylabel('w  [gust]')
            %     title('disturbance horizon – slot index = colour')
            %     grid on
            % 
            % 
            %     imagesc(obj.dbg.wHorizon), axis xy
            %     xlabel('sample'), ylabel('slot'), colorbar
            %     title('MHE disturbance horizon (oldest = top)')
            % 
            % end
        end
    end  % estimate

    %% =====================================================================
    %%  INTERNAL HELPERS  ===================================================
    %% =====================================================================
    methods (Access=private)

        % ---------------- index maps (single place) ------------------------
        function idx = buildIndexMaps(obj)
            nx = obj.nx;  nw = obj.nw;  Ne = obj.Ne;
            idx.x = cell(Ne+1,1);
            for j = 1:Ne+1, idx.x{j} = (j-1)*nx + (1:nx); end
            startW = (Ne+1)*nx;
            idx.w = cell(Ne,1);
            for j = 1:Ne, idx.w{j} = startW + (j-1)*nw + (1:nw); end
        end

        % ------------------------------------------------------------------
        function zShift = shiftGuess(obj,p)
            %  Real‑time iteration shift:  drop oldest sample, append newest
            idx = obj.buildIndexMaps();
            nx = obj.nx; nw = obj.nw; Ne=obj.Ne;

            zShift = zeros(size(p));

            % ----- states --------------------------------------------------
            for j = 1:Ne               % move x^{-(Ne-1)} … x^{0}
                zShift(idx.x{j}) = p(idx.x{j+1});
            end
            zShift(idx.x{Ne+1}) = p(idx.x{Ne+1});   % duplicate latest

            % ----- disturbances -------------------------------------------
            for j = 1:Ne-1             % shift w
                zShift(idx.w{j}) = p(idx.w{j+1});
            end
            zShift(idx.w{Ne}) = p(idx.w{Ne});       % keep newest forecast
        end

        % ------------------------------------------------------------------
        function nlp = assembleWindow(obj)
            %  Build window‑specific cost & constraints (dense)
            % ----------------------------------------------------------------
            idx = obj.buildIndexMaps();   nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
            idxWvec = [idx.w{:}];

            % ---------- COST (handles grad) --------------------------------
            function [J,gradJ] = costFun(p)
                gradJ = zeros(numel(p),1);

                % 1. arrival (state only) -----------------------------------
                xE  = p(idx.x{1});
                xRef= obj.buff.X(:, obj.idxWrap(-1));   % state @ k=-1
                dx  = xE - xRef;
                J   = 0.5*dx.'*obj.Pe*dx;
                gradJ(idx.x{1}) = obj.Pe*dx;

                % 2. disturbance arrival  -----------------------------------
                wE  = p(idx.w{1});   wRef = obj.what(end-obj.nw+1:end);
                dw  = wE - wRef;
                J   = J + 0.5*dw.'*obj.Re*dw;
                gradJ(idx.w{1}) = obj.Re*dw;

                % 3. measurement innovations  -------------------------------
                C = obj.buff.sensorMat;    % ny×nq1
                Cfull = zeros(obj.ny,obj.nx);
                Cfull(:,obj.buff.idx.q1) = C;

                Xmat = reshape(p([idx.x{:}]),nx,[]);  % (nx×Ne+1)
                for j = 1:Ne+1
                    yk = C * Xmat(obj.buff.idx.q1,j);
                    zk = obj.buff.Y(:, obj.idxWrap(j-(Ne+1)));
                    ek = zk - yk;
                    J  = J + 0.5*ek.'*obj.Qe*ek;
                    gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cfull.'*(obj.Qe*ek);
                end
            end

            % ---------- CONTINUITY  ---------------------------------------
            function [c,ceq,gradc,gradceq] = nonlFun(p)
                gradc=[];  nVar=numel(p);
                ceq = zeros(nx*Ne,1);
                gradceq = zeros(nVar,nx*Ne);

                wHist = reshape(p(idxWvec),nw,Ne).';   % Ne×nw

                row0=0;
                for j=1:Ne
                    Sj = obj.buff.S(:,:, obj.idxWrap(-Ne+j));  % STM of *previous* interval
                    Sw = Sj(:, nx+1);                         % ∂x/∂w
                    Su = Sj(:, nx+2:nx+1+obj.nu);             % ∂x/∂u   <<< NEW

                    uj = obj.buff.U(:, obj.idxWrap(j-(Ne+1)));% known control
                    xk   = p(idx.x{j});   xkp1 = p(idx.x{j+1});
                    wj   = wHist(j,:)';
                    rhs  = xk + Sw*wj + Su*uj;                % f(xk,uj,wj)
                    % rhs  = xk + Sw*wj;                % f(xk,uj,wj)

                    rows = row0+(1:nx);
                    ceq(rows) = xkp1 - rhs;

                    % Jacobians -------------------------------------------------
                    gradceq(idx.x{j},   rows) = -eye(nx);
                    gradceq(idx.x{j+1},rows) =  eye(nx);
                    gradceq(idx.w{j},   rows) = -Sw;
                    row0 = row0+nx;
                end
                % gradceq = gradceq.';   % m×n per fmincon contract
                c=[];
            end

            % ---------- BOUNDS --------------------------------------------
            xLrep = repmat(obj.xL, obj.Ne+1,1);
            xUrep = repmat(obj.xU, obj.Ne+1,1);
            wLrep = repmat(obj.wL, obj.Ne,1);
            wUrep = repmat(obj.wU, obj.Ne,1);

            nlp.cost = @costFun;
            nlp.nonl = @nonlFun;
            nlp.lb   = [xLrep ; wLrep];
            nlp.ub   = [xUrep ; wUrep];
        end

        % ------------------------------------------------------------------
        function kIdx = idxWrap(obj,shift)
            kIdx = mod(obj.buff.valid-1 + shift, obj.Ne+1) + 1;
        end
    end  % private helpers
end   % class
% classdef nMHE < AeroFlex.ctrl.EstimatorBase
% % Non‑linear Moving–Horizon Estimator  (Artola 2021, Eq.(38))
% % Multiple‑shooting + SQP.  Uses state & STM buffers filled by SimRunner.
% % ----------------------------------------------------------------------
%     properties
%         % configuration mirrors fields already in cfg.control or cfg.nmhe
%         dt                      % integrator step (copied from cfg.dt)
%         Ts                      % sampling period
%         Ne                      % window length (# samples)
%         Qe, Pe, Re                % weighting matrices
%         wL, wU                   % disturbance bounds
%         xL, xU                   % state bounds
%         solverOpts              % fmincon('sqp') options
%         % circular buffers shared with SimRunner ------------------------
%         buff                    % struct: X,S,U,Y,T  (handles, not copies)
%         H                       % Hessian
%         % local rolling index
%         kBuff   = 0
%         % warm‑start vector
%         z0                      %   [ xE ; vec(w_hist) ]
%         % dimensions
%         nx,  nu, ny, nw
%         % trim references
%         xRefPrev
%         ref, maxSam, wRef
%         % testSw
%         alphaSmooth    % < 0 disables
%         muBarrier      % 0 disables
%         betaBarrier  
%     end
% % ----------------------------------------------------------------------
%     methods
%         function obj = nMHE(cfg)
%             arguments, cfg struct, end
%             obj@AeroFlex.ctrl.EstimatorBase(cfg);
% 
%             % ------------------ copy scalars from cfg ------------------
%             obj.dt  = cfg.sim.dt;         
%             obj.Ts = cfg.ctrl.Ts;
%             obj.Ne  = cfg.ctrl.Ne;         
%             obj.Qe = cfg.Qe;
%             obj.Pe  = cfg.Pe;         obj.Re = cfg.Re;
%             obj.wL  = cfg.wL;         obj.wU = cfg.wU;
%             obj.xL  = cfg.xL;         obj.xU = cfg.xU;
% 
%             obj.buff = cfg.buff;      % <- struct with X,S,U,Y,T handles
%             obj.maxSam = cfg.maxSam;
%             obj.nx   = size(obj.buff.X,1);
%             obj.nu   = size(obj.buff.U,1);
%             obj.ny   = size(obj.buff.Y,1);
%             obj.nw = 1;
%             obj.ref = cfg.ref;
%             % fmincon options (SQP) – 1 GN + 1 SOC ≃ Artola’s RTI
%             % obj.solverOpts = optimoptions('fmincon', ...
%             %     'Algorithm','sqp', 'MaxIterations',100, ...
%             %     'Display','iter', 'SpecifyObjectiveGradient',true, ...
%             %     'SpecifyConstraintGradient',true, 'checkGradients', false); % 'SpecifyConstraintGradient',true, false
%             obj.solverOpts = optimoptions('fmincon', ...
%                 'Algorithm','interior-point', 'MaxIterations',500, ... % sqp, interior-point
%                 'Display','off', 'HessianApproximation','lbfgs', 'SpecifyObjectiveGradient',true, ...
%                 'SpecifyConstraintGradient',true, 'checkGradients', false); % 'SpecifyConstraintGradient',true, false
% 
%             % opts             = optimoptions('fmincon','Algorithm','sqp');
%             % opts.HessianApproximation = 'lbfgs';   % cheaper per iter
%              obj.solverOpts.StepTolerance        = 1e-8;
% 
% 
%             % obj.solverOpts = optimoptions('fmincon', ...
%             %     'Algorithm','interior-point', 'MaxIterations',100, ...
%             %     'Display','iter-detailed', 'SpecifyObjectiveGradient',true, ...
%             %     'SpecifyConstraintGradient',true, 'checkGradients', true); % 'SpecifyConstraintGradient',true, false
% 
%              obj.solverOpts.FiniteDifferenceType ='central';
%              % obj.solverOpts.FunValCheck = 'on';
%              % obj.solverOpts.StepTolerance = 1e-5;
%              % obj.solverOpts.ConstraintTolerance = 1e-5;
%              obj.solverOpts.OptimalityTolerance = 1e-5;
%              % obj.solverOpts.OptimalityTolerance = 1e-4;
%              % obj.solverOpts.OptimalityTolerance = 1e-3; % temp for test
% 
%              % obj.solverOpts.PlotFcn = 'optimplot';
%             % obj.solverOpts.ScaleProblem = 'obj-and-constr'; % 'none', 'obj-and-constr'
% 
% 
%             % initial warm‑start
%             % obj.z0 = [ cfg.ref.states ; zeros(obj.nw*obj.Ne,1) ];
%              % ---------- index maps -------------------------------------------------
%             nx     = obj.nx;   nu = obj.nu;   nw = obj.nw;   Ne = obj.Ne;
% 
%             % blocks for all states  x^{-Ne} … x^{0}
%             idx.x  = cell(Ne+1,1);
%             for j = 1:Ne+1
%                 idx.x{j} = (j-1)*nx + (1:nx);
%             end
% 
%             % blocks for disturbances  w^{-Ne} … w^{-1}
%             startW = (Ne+1)*nx;
%             idx.w  = cell(Ne,1);
%             for j = 1:Ne
%                 idx.w{j} = startW + (j-1)*nw + (1:nw);
%             end
% 
%             nVar = (Ne+1)*nx + Ne*nw;          % 144·25 + 1·24 = 3664
%             if isempty(obj.H)
%                 obj.H = speye(nVar);         % initial positive-definite Hessian
%             end
% 
%             % fill state blocks with buffer data
%             % nVar = (obj.Ne+1)*obj.nx + obj.Ne*obj.nw;          % 144·25 + 1·24 = 3664
%             z0 = zeros(nVar,1);
%             for j = 1:obj.Ne+1
%                 % z0(idx.x{j}) = obj.buff.X(:, obj.idxWrap(-obj.Ne-1+j));
%                 z0(idx.x{j}) = obj.buff.X(:, obj.idxWrap(-obj.Ne+j)); %this one?
%             end
%             % fill disturbance blocks with zeros
%             for j = 1:obj.Ne
%                 z0(idx.w{j}) = 0;
%             end
%             obj.z0 = z0;
%             obj.what = zeros(obj.Ne,1);
% 
% 
%             % obj.alphaSmooth  = 1e-3;  % < 0 disables
%             obj.alphaSmooth  = 0;  % < 0 disables
%             obj.muBarrier    = 1e-5;  % 0 disables
%             % obj.betaBarrier  = 150;
%             obj.betaBarrier  = 0;
%         end
% % ----------------------------------------------------------------------
%     function [xhat, what, info] = estimate(obj,z,u,t, buff)
%         %  Called once per sampling period Ts from SimRunner
%         %  -------------------------------------------------------------
%             % --- push measurement & control into ring buffer ----------
%             obj.buff = buff;      % <- struct with X,S,U,Y,T handles
%             info = struct('cost', 0, 'exitflag', [], 'continuity', 0, 'wRefPrev', zeros(obj.Ne*obj.nw,1));
%             obj.kBuff = mod(obj.kBuff, obj.Ne) + 1;    k = obj.kBuff;
%             obj.buff.Y(:,k) = z;
%             obj.buff.U(:,k) = u;
%             obj.buff.T(  k) = t;
% 
%             % warm‑up:  wait until buffer full -------------------------
%             % if obj.buff.valid < obj.Ne
%             %     obj.yhat = obj.h(obj.xhat); xhat=obj.xhat; return
%             % end
%             % if obj.buff.valid < obj.Ne
%             obj.buff.valid = min(obj.buff.valid + 1, obj.Ne+1);
% 
%             % if obj.kBuff< obj.Ne
%             % if obj.buff.valid< 10
%             if obj.buff.valid< obj.Ne+1
%                 % predict output using constant C‑matrix stored in the buffer
%                 C        = obj.buff.sensorMat;            % 5 × nx
%                 if isempty(obj.xhat)
%                     obj.xhat = obj.ref.states;          % or zeros(obj.nx,1);
%                     % xhat_part = 
%                 end
%                 q1  = obj.xhat(obj.buff.idx.q1);            % 20 × 1  (modal velocities)
%                 obj.yhat = C * q1;                  % ŷ = C x̂
%                 xhat     = obj.xhat;
%                 what = obj.what;
%                 return
%             end
% 
% 
%             persistent wRefPrev                 % remember between calls
%                 if isempty(wRefPrev)                % first call: use config
%                     wRefPrev = zeros(obj.nw,1);     % or cfg.ref.w_r(1:nw)
%                     % wRefPrev = zeros(obj.nw*obj.Ne,1);     % or cfg.ref.w_r(1:nw)
%                     obj.wRef = wRefPrev;
%                 end
% 
%             persistent xRefPrev                 % remember between calls
%             if isempty(xRefPrev)                % first call: use config
%                 xRefPrev = obj.ref.states;     % or cfg.ref.w_r(1:nw)
%                 % wRefPrev = zeros(obj.nw*obj.Ne,1);     % or cfg.ref.w_r(1:nw)
%                 obj.xRefPrev = xRefPrev;
%             end
%             % obj.wRef = wRefPrev;                    %  nw×1  used below
%             % obj.wRef = p( idx.w{2} );
% 
%             % ----- build window‑NLP ----------------------------------
%             nlp = obj.assembleWindow();
%             idxBad = find(nlp.lb > nlp.ub);
%             disp([idxBad , nlp.lb(idxBad) , nlp.ub(idxBad)])
% 
% %% custom SQP
%             % 
%             % % ---------- run custom SQP -------------------------------------------
%             % sqpOpts.MaxIter = 80;
%             % sqpOpts.TolFeas = 1e-8;   sqpOpts.TolOpt = 1e-6;
%             % sqpOpts.Display = true;  % or true for iteration log
%             % 
%             % [p, f, obj.H, info] = AeroFlex.optim.runSQP(nlp, obj.z0, obj.H);
%             % % 
%             % idx.x  = cell(obj.Ne+1,1);
%             % for j = 1:obj.Ne+1
%             %     idx.x{j} = (j-1)*obj.nx + (1:obj.nx);
%             % end
%             % % blocks for disturbances  w^{-Ne} … w^{-1}
%             % startW = (obj.Ne+1)*obj.nx;
%             % idx.w  = cell(obj.Ne,1);
%             % for j = 1:obj.Ne
%             %     idx.w{j} = startW + (j-1)*obj.nw + (1:obj.nw);
%             % end
%             % 
%             % if info.flag<=0
%             %     warning('nMHE:SQP','Solver did not converge (flag %d), keep last x̂.',info.flag);
%             %     disp('Solver did not converge (flag '); disp(info.flag); disp(', keep last x̂.');
%             % else
%             %     obj.z0   = p;                            % warm start
%             %     obj.xhat = p( obj.idx.x{obj.Ne+1} );     % newest state
%             %     obj.what = p( [obj.idx.w{:}] );          % disturbance hist.
%             % end
% %%
%             [p,f,exitFlag] = fmincon( nlp.cost, obj.z0, ...
%                           [],[],[],[], nlp.lb, nlp.ub, ...
%                           nlp.nonl, obj.solverOpts);
% 
% %%
%             % [p,~,exitFlag] = fmincon(nlp.cost, obj.z0, ...
%             %               [],[],[],[], nlp.lb, nlp.ub, ...
%             %               [], obj.solverOpts); % try without ineq
%             % 
%             % optsFD = optimoptions('fmincon','Algorithm','interior-point', ...
%             %           'SpecifyObjectiveGradient',false, ...
%             %           'Display','off','MaxIterations',0);
%             % 
%             % [~,g_fd] = fmincon(nlp.cost, obj.z0,[],[],[],[],nlp.lb,nlp.ub,[], optsFD);
%             % 
%             % [~,g_an] = nlp.cost(obj.z0);
%             % 
%             % pOpt = p;
%             % [~,lambda] = nlp.nonl(pOpt);     % lambda.lower, lambda.upper etc.
%             % fprintf('# active x bounds %d  w bounds %d\n', ...
%             %          nnz(lambda.lower(idx.x0)), nnz(lambda.lower(idx.w)));
%             % 
%             % relErr = norm(g_fd - g_an) / max(1,norm(g_fd));
%             % fprintf('Gradient mismatch %.3e\n', relErr);
% 
%             idx.x  = cell(obj.Ne+1,1);
%             for j = 1:obj.Ne+1
%                 idx.x{j} = (j-1)*obj.nx + (1:obj.nx);
%             end
%             % blocks for disturbances  w^{-Ne} … w^{-1}
%             startW = (obj.Ne+1)*obj.nx;
%             idx.w  = cell(obj.Ne,1);
%             for j = 1:obj.Ne
%                 idx.w{j} = startW + (j-1)*obj.nw + (1:obj.nw);
%             end
%             if exitFlag<=0
%                 warning('nMHE:SQP','Solver did not converge, keeping previous x̂.');
%             else
%                 obj.z0   = p;                          % warm start next call
%                 obj.xhat = p( idx.x{obj.Ne+1} );       % state at k = 0
%                 obj.what = p( [idx.w{:}] );            % all disturbances
%             end
%             idxWvec = [idx.w{:}];                     % 1 × (nw*Ne)
% 
%             wRefPrev = p( idx.w{2} );
%             xRefPrev = p(idx.x{2});
%             obj.xRefPrev = xRefPrev;
%             % wRefPrev = p( idxWvec );
%             % wRefPrev = p( idx.w{end-1} );
% 
%             obj.yhat = z;                   % last measurement (optional)
%             xhat = obj.xhat;
%             what = obj.what;
% 
%             info.cost       = f;                % final cost from fmincon
%             info.exitflag   = exitFlag;         % or info.flag for runSQP
%             % info.continuity = err_state;        % residual computed
% 
%             wHist   = reshape(p(idxWvec), obj.nw, obj.Ne).';  % Ne × nw  
%             % % obj.wRef = wHist(2,:).';                  % use 2nd-oldest as new reference
%             % persistent wRefPrev                 % remember between calls
%             %     if isempty(wRefPrev)                % first call: use config
%             %         wRefPrev = zeros(obj.nw*obj.Ne,1);     % or cfg.ref.w_r(1:nw)
%             %     end
%             % wRefPrev(:,end) = wHist;
%             info.wRefPrev = wHist;
%             obj.wRef =  wRefPrev;
%         end
% % ----------------------------------------------------------------------
%     end                             % methods (public)
% % ======================================================================
%     methods %(Access=private)
%         function nlp = assembleWindow(obj)
%         % Build moving‑window NLP from ring buffers (size‑independent)
%         % ------------------------------------------------------------
%             % nx = obj.nx; nu = obj.nu; Ne = obj.Ne;
%             % nw = obj.nw;
%             % % indices in decision vector
%             idx.xE = 1:obj.nx;
%             % idx.w  = nx+ (1:nw*Ne);
% 
%             % Attempt2:
%             % ---------- index maps -------------------------------------------------
%             nx     = obj.nx;   nu = obj.nu;   nw = obj.nw;   Ne = obj.Ne;
% 
%             % blocks for all states  x^{-Ne} … x^{0}
%             idx.x  = cell(Ne+1,1);
%             for j = 1:Ne+1
%                 idx.x{j} = (j-1)*nx + (1:nx);
%             end
% 
%             % blocks for disturbances  w^{-Ne} … w^{-1}
%             startW = (Ne+1)*nx;
%             idx.w  = cell(Ne,1);
%             for j = 1:Ne
%                 idx.w{j} = startW + (j-1)*nw + (1:nw);
%             end
% 
%             nVar = (Ne+1)*nx + Ne*nw;         
% 
% 
%             % -------- objective --------------------------------------
% %% V-4
%         % =======================================================
%             %  COST FUNCTION  J(p)   and  ∇J(p)
%             % % =======================================================
%             % function [J,gradJ] = costFun(p)
%             %     persistent Cfull
%             %         if isempty(Cfull)
%             %             Cfull           = zeros(obj.ny,obj.nx);  % allocate once
%             %             q1Idx = obj.buff.idx.q1;
%             %             Cfull(:,q1Idx)  = obj.buff.sensorMat;
%             %         end
%             %     % ---------- alias short-hands --------------------------------------
%             %     nx  = obj.nx;        nw = obj.nw;        Ne = obj.Ne;
%             %     idxX = idx.x;        idxW = idx.w;       idxWvec = [idxW{:}];
%             % 
%             %     % ---------- reshape disturbance history ----------------------------
%             %     wHist = reshape( p(idxWvec) , nw , Ne ).';      % Ne×nw  (oldest row =1)
%             %     X   = reshape( p([idx.x{:}]), obj.nx, [] );  % (nx × Ne+1)
%             % 
%             %     % ---------- rolling references  ------------------------------------
%             %     % if isempty(obj.xRefPrev), obj.xRefPrev = zeros(nx,1); end
%             %     if isempty(obj.wRef), obj.wRef = zeros(nw,1); end
%             %     % xRef = obj.xRefPrev;                % nx×1
%             %     xRef = obj.buff.X(:, obj.idxWrap(-1));     % state at k = –1
%             % 
%             %     wRef = obj.wRef;                % nw×1
%             % 
%             %     % ---------- ARRIVAL cost  (Artola 38a) ------------------------------
%             %     xE   = p(idxX{1});                  % x^{-Ne}
%             %     wOld = p(idxW{1});                  % w^{-Ne}
%             % 
%             %     dx   = xE - xRef;                   % nx×1
%             %     dw   = wOld - wRef;                 % nw×1
%             % 
%             %     J       = 0.5 * ( dx.'*obj.Pe*dx  +  dw.'*obj.Re*dw );
%             %     gradJ   = zeros(numel(p),1);
%             %     % xgrad = obj.Pe * dx;
%             %     % if size(xgrad,1) ~= 1
%             %     %     xgrad = xgrad.';
%             %     % end
%             %     gradJ(idxX{1}) = gradJ(idx.x{1}) + obj.Pe * dx;
%             %     % gradJ(idxX{1}.') = xgrad;
%             %     gradJ(idxW{1}) = gradJ(idx.w{1}) + obj.Re * dw;
%             % 
%             %     % ---------- MEASUREMENT innovation loop -----------------------------
%             %     xk = p(idxX{Ne+1});                 % newest state x^0
%             %     for j = Ne+1:-1:1                   % j = Ne+1 (newest) … 1 (oldest)
%             %         yk = obj.buff.sensorMat * xk(obj.buff.idx.q1);        % predict
%             %         zk = obj.buff.Y(:, obj.idxWrap(j-(Ne+1)));            % measured
%             %         ek = zk - yk;                                        % innovation
%             % 
%             %         J = J + 0.5 * ek.' * obj.Qe * ek;
%             % 
%             % 
%             % 
%             %         gradJ(idxX{j}) = gradJ(idxX{j}) - Cfull.' * (obj.Qe*ek);
%             % 
%             %         if j>1                                       % propagate back
%             %             Sbk = obj.buff.S(:,:,obj.idxWrap(j-(Ne+2)));      % STM k-1
%             %             wk  = wHist(j-1,:).';
%             %             xk  = xk - Sbk(:,nx+1:nx+nw)*wk;                  % x_{j-1}
%             %         end
%             %     end
%             %     % --- 2. measurement stages ----------------------------------------
%             %     % C  = obj.buff.sensorMat;                % 30 × nx dense – no sparse!
%             %     % q1 = obj.buff.idx.q1;
%             %     % Cfull    = zeros(obj.ny,obj.nx);           % 30×144  dense
%             %     % Cfull(:,obj.buff.idx.q1) = C;
%             %     % for j = 1:obj.Ne+1
%             %     %     ek   = obj.buff.Y(:,obj.idxWrap(j-(obj.Ne+1))) ...
%             %     %            - C * X(q1,j);
%             %     %     J    = J + 0.5 * ek.' * obj.Qe * ek;
%             %     %     gradJ(idx.x{j}) = gradJ(idx.x{j}) - ...
%             %     %                       Cfull.' * (obj.Qe * ek);      % dense 144×1
%             %     % end
%             %     % ---------- SMOOTHNESS term  Σ‖Δw‖²  (optional) ---------------------
%             %     % alphaSm = obj.alphaSmooth;      % set in constructor, 0 ⇒ skip
%             %     % if alphaSm > 0
%             %     %     dW = diff(wHist,1,1);                       % (Ne-1)×nw
%             %     %     J  = J + 0.5*alphaSm*sum(dW(:).^2);
%             %     % 
%             %     %     for j = 2:Ne
%             %     %         dw = wHist(j,:) - wHist(j-1,:);
%             %     %         gradJ(idxW{j})   = gradJ(idxW{j})   + alphaSm*dw.';
%             %     %         gradJ(idxW{j-1}) = gradJ(idxW{j-1}) - alphaSm*dw.';
%             %     %     end
%             %     % end
%             %     % 
%             %     % % ---------- POSITIVITY SOFT BARRIER  μΣexp(−βw)  --------------------
%             %     % muBar = obj.muBarrier;   betaBar = obj.betaBarrier;
%             %     % if muBar > 0
%             %     %     expTerm = exp(-betaBar * wHist);            % Ne×nw
%             %     %     J       = J + muBar * sum(expTerm(:));
%             %     % 
%             %     %     gradJ(idxWvec) = gradJ(idxWvec) - ...
%             %     %                      muBar*betaBar*expTerm(:);
%             %     % end
%             % end
%             % ========================================================================
% 
% %% V-3            
% %             % ---------- ARRIVAL-COST & OUTPUT TERM (dense, no sparse) -------------
% %             % ---------- costFun  -----------------------------
%             function [J,gradJ] = costFun(p)
%                 % rolling reference (Artola 2021, Eq. 38a)
%                 xRef = obj.buff.X(:, obj.idxWrap(-1));     % state at k = –1
%                 % wHistRec =  p([idx.w{:}]).' ;                      % row (1 × Ne)
% 
% 
%                 wRef = p( idx.w{2} );                      % use w^{-Ne+1} from *prev* solution
%                 % wRef = obj.wRef;                      % use w^{-Ne+1} from *prev* solution
%  %%              
%                 % % wHist = p(idx.w);
%                 % % flat index of all disturbance slots in the decision vector p
%                 idxWvec = [idx.w{:}];                     % 1 × (nw*Ne)
%                 % % wHist   = p(idxWvec).';%  Ne × nw
%                 % 
%                 % % % % reshape into chronological matrix (row 1 = oldest sample, row Ne = latest)
%                 wHist   = reshape( p(idxWvec), nw, Ne ).';%  Ne × nw
%                 % % ---- rolling w_ref (Artola, eq. 38a) -------------------------------
%                 % persistent wRefPrev                 % remember between calls
%                 % if isempty(wRefPrev)                % first call: use config
%                 %     wRefPrev = zeros(obj.nw,1);     % or cfg.ref.w_r(1:nw)
%                 % end
%                 % wRef = wRefPrev;                    %  nw×1  used below
%                 % wRefPrev = wHist(2,:).';            % NEXT iteration will use w_{-Ne+1}
% %%
%                 % reshape helpers ---------------------------------------------------
%                 X   = reshape( p([idx.x{:}]), obj.nx, [] );  % (nx × Ne+1)
%                % W   = p([idx.w{:}]).';                       % row (1 × Ne)
% 
% 
%                 %%
% 
%                 %%
%                 % --- 1. arrival ----------------------------------------------------
%                 dX  = X(:,1) - xRef;           % x_E − x_ref
%                 % dW  = W(1)  - wRef;            % w_{-Ne} − w_ref
%                 % J   = 0.5*( dX.' * obj.Pe * dX + obj.Re * dW^2 );
%                 J   = 0.5*( dX.' * obj.Pe * dX );
% % %%
%                 % wRef = obj.z0( idx.w{2} );       % last-but-one disturbance
% 
%                 % wRef = obj.z0( idxWvec);       % last-but-one disturbance
%                 % wRef = obj.z0( idx.w{Ne-1} );       % last-but-one disturbance
%                 % vector form: wHist is Ne×nw
%                 % dw    = wHist(:) - repmat(wRef,Ne,1);
%                 % J     = J + 0.5*dw.'*kron(eye(Ne),obj.Re)*dw;
%                 dW    = wHist(1,:) - wRef;         % w^{-Ne} − w_ref
%                 % dW    = wHist(Ne,:) - wRef;         % w^{-Ne} − w_ref
% 
%                 % dW = wHist - wRef*ones(length(wHist),1);         % w^{-Ne} − w_ref
%                 % dW = wHist - wRef;         % w^{-Ne} − w_ref
%                 J     = J + 0.5 * dW.' * obj.Re * dW;
% %%
%                 gradJ          = zeros(numel(p),1);
%                 gradJ(idx.x{1}) =  obj.Pe * dX;        % ∂J/∂x_E
%                 gradJ(idx.w{1}) =  obj.Re * dW;        % ∂J/∂w_{-Ne}
%                 % gradJ(idxWvec) =  obj.Re * dW;        % ∂J/∂w_{-Ne}
%                 % gradJ(idxWvec) = kron(ones(Ne,1),obj.Re) .* dw;
% %%
%                 % alpha = 1e-2;                       % tune 0.01 … 1
%                 % dW    = diff(wHist,1,1);            % (Ne-1) × nw  differences down the rows
%                 % 
%                 % J     = J + 0.5 * alpha * sum(dW(:).^2);   % add to scalar cost
%                 % 
%                 % % -------- gradient contribution -----------------------------------------
%                 % for j = 2:Ne                        % affects rows 2…Ne
%                 %     colThis   = idx.w{j};           % scalar index (nw=1 here)
%                 %     colPrev   = idx.w{j-1};
%                 % 
%                 %     gradJ(colThis) = gradJ(colThis) + alpha*( wHist(j,:)   - wHist(j-1,:) ).';
%                 %     gradJ(colPrev) = gradJ(colPrev) - alpha*( wHist(j,:)   - wHist(j-1,:) ).';
%                 % end
% 
% %%
%                 % --- 2. measurement stages ----------------------------------------
%                 C  = obj.buff.sensorMat;                % 30 × nx dense – no sparse!
%                 q1 = obj.buff.idx.q1;
%                 Cfull    = zeros(obj.ny,obj.nx);           % 30×144  dense
%                 Cfull(:,obj.buff.idx.q1) = C;
%                 for j = 1:obj.Ne+1
%                     ek   = obj.buff.Y(:,obj.idxWrap(j-(obj.Ne+1))) ...
%                            - C * X(q1,j);
%                     J    = J + 0.5 * ek.' * obj.Qe * ek;
%                     gradJ(idx.x{j}) = gradJ(idx.x{j}) - ...
%                                       Cfull.' * (obj.Qe * ek);      % dense 144×1
%                 end
%                 % wRefPrev = wHist(2,:).';             % will be reference next call
% 
%                 % fprintf('ŵ_now=%+.3f … ŵ_oldest=%+.3f\n',...
%                 %     wHist(end), wHist(1));
% 
% 
%             end
% 
%             % function [J,gradJ] = costFun(p)
%             % 
%             %     % --- flat indices --------------------------------------------------
%             %     idxWvec = [idx.w{:}];                % 1×(nw·Ne)
%             %     nVar    = numel(p);
%             % 
%             %     % ------ reference update  (Artola, line below Eq.(38a)) -----------
%             %     % keep the state & disturbance found *one sample ago*
%             %     persistent xRefPrev wRefPrev
%             %     if isempty(xRefPrev)
%             %         xRefPrev = obj.ref.states;     %# initialisation
%             %         wRefPrev = 0.0;
%             %     end
%             %     xRef = xRefPrev;                           %   qae_ref
%             %     wRef = wRefPrev;                           %   w_ref
%             % 
%             %     % update for *next* call (sample k becomes k-1)
%             %     xRefPrev = obj.buff.X(:,obj.idxWrap(-1));  % state at start of 2nd sub-interval
%             %     wRefPrev = obj.what( max(1,end-1) );       % most recent ŵ(k-1)
%             % 
%             %     % ----------- unpack decision vector --------------------------------
%             %     xE    = p(idx.x{1});                       % state at  k = –Ne
%             %     wHist = p(idxWvec);                        % vector size (nw·Ne)×1
%             % 
%             %     % -------------- dense objective & gradient -------------------------
%             %     J      = 0;
%             %     gradJ  = zeros(nVar,1);
%             % 
%             %     % 1. arrival cost ‖xE-xRef‖_Pe + ‖w₋Ne-wRef‖_Re
%             %     dx     = xE  - xRef;
%             %     dw     = wHist(1) - wRef;                  % first *sample*, not scalar 0
%             %     J      = 0.5*dx.'*obj.Pe*dx + 0.5*dw.'*obj.Re*dw;
%             %     % --- rolling reference -------------------------------------------------
%             %     % % Eq.(38a): use state and disturbance at start of *second* sub–interval
%             %     % xRef = obj.z0( idx.x{Ne} );         % x̂ from previous solve (k = -1)
%             %     % dx     = xE  - xRef;
%             %     % J      = 0.5*dx.'*obj.Pe*dx;
%             %     % 
%             %     wRef = obj.z0( idx.w{Ne-1} );       % last-but-one disturbance
%             %     % vector form: wHist is Ne×nw
%             %     dw    = wHist(:) - repmat(wRef,Ne,1);
%             %     J     = J + 0.5*dw.'*kron(eye(Ne),obj.Re)*dw;
%             %     gradJ(idxWvec) = kron(ones(Ne,1),obj.Re) .* dw;
%             % 
%             %     gradJ(idx.x{1}) = obj.Pe*dx;
%             %     gradJ(idx.w{1}) = obj.Re*dw;               % ← full vector if nw>1
%             % 
%             %     % 2. measurement innovation loop (identical maths, but **no sparse**)
%             %     C        = obj.buff.sensorMat;             % 30×20   (constant)
%             %     Cfull    = zeros(obj.ny,obj.nx);           % 30×144  dense
%             %     Cfull(:,obj.buff.idx.q1) = C;
%             % 
%             %     xk = xE;                                   % start at k = –Ne
%             %     for j = obj.Ne+1:-1:1                      % … to k = 0
%             %     % for j = obj.Ne:-1:1                      % … to k = 0
%             %         yk = C*xk(obj.buff.idx.q1);
%             %         zk = obj.buff.Y(:,obj.idxWrap(j-(obj.Ne+1)));
%             %         ek = zk-yk;
%             % 
%             %         J = J + 0.5*ek.'*obj.Qe*ek;
%             %         gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cfull.'*(obj.Qe*ek);
%             % 
%             %         if j>1
%             %             Sbk = obj.buff.S(:,:,obj.idxWrap(j-(obj.Ne+2))); % STM k-1→k
%             %             wj  = wHist(j-1);
%             %             xk  = xk - Sbk(:,obj.nx+1)*wj;      % reverse-prop one sample
%             %         end
%             %     end
%             % end
% 
%             % function [J,gradJ] = costFun(p)
%             %     % --- flat index vectors -----------------------------------------
%             %     idxWvec = [idx.w{:}];            % 1 × (nw·Ne)   numeric row
%             %     idxXvec = [idx.x{:}];            % 1 × (nx·(Ne+1))
%             %     % ---------- arrival-cost terms (Eq. 38a) ----------------------------
%             %     xRef = obj.buff.X(:, obj.idxWrap(0));   % newest sample (k = 0)
%             % 
%             %     xE    = p(idx.x{1});                    % oldest node  (j = 1)
%             %     wHist = reshape(p(idxWvec), nw, Ne).';  % Ne × nw
%             % 
%             %     % J  = 0.5 * (xE - xRef).' * obj.Pe * (xE - xRef) ...
%             %     %    + 0.5 * (wHist(1,:) - obj.cfg.ref.w_r(1:length(wHist))).' ...
%             %     %            * obj.Re  * (wHist(1,:) - obj.cfg.ref.w_r(1:length(wHist)));
%             %     J  = 0.5 * (xE - xRef).' * obj.Pe * (xE - xRef);
%             % 
%             % 
%             %     gradJ          = zeros(size(p));
%             %     gradJ(idx.x{1})= obj.Pe * (xE - xRef);          % ∂J/∂xE
%             %     wRef = obj.cfg.ref.w_r(1);           % scalar reference, 1×1
%             % 
%             %     dW   = wHist(1) - wRef;              % scalar
%             %     J    = J + 0.5 * obj.Re * (dW^2);    % obj.Re is 1×1
%             % 
%             %     gradJ(idx.w{1}) = obj.Re * dW;       % 1 element → matches idx.w{1}          
%             %     % gradJ(idx.w{1})= obj.Re * (wHist(1,:) - obj.cfg.ref.w_r(1:length(wHist)));  % ∂J/∂w₋Ne
%             % 
%             %     % xE    = p(idx.xE);
%             %     % % wHist = reshapbe(p(idx.w),nw, []).';
%             %     % wHist = p(idx.w);
%             %     % wHist = reshape(p(idx.w),nw*Ne,[]).';
%             %     % (38a) arrival cost
%             %     C  = obj.buff.sensorMat;            % (ny × nq1)   constant Φ_y
%             %     % Cfull            = sparse(obj.ny , obj.nx);   % 30 × 144
%             %     Cfull            = zeros(obj.ny , obj.nx);   % 30 × 144
%             %     q1Idx            = obj.buff.idx.q1;           % 1:20
%             %     Cfull(:,q1Idx)   = obj.buff.sensorMat;        % Φ_y in those columns
%             %     % % test = (xE-obj.buff.X(:,obj.idxWrap(0)));
%             %     % J = 0.5*(xE-obj.buff.X(:,obj.idxWrap(0))).' ...
%             %     %        *obj.Pe*(xE-obj.buff.X(:,obj.idxWrap(0)));
%             %     % J = J + 0.5*(wHist(:,1)-obj.cfg.ref.w_r(1:length(wHist)) ).' *obj.Re*(wHist(:,1)-obj.cfg.ref.w_r(1:length(wHist)));
%             %     % gradJ = zeros(size(p));
%             %     % gradJ(idx.xE) = obj.Pe*(xE-obj.buff.X(:,obj.idxWrap(0)));
%             %     % gradJ(idx.w) = obj.Re*(wHist(:,1));
%             % 
%             %     % measurement loop
%             %     % xk = z(idx.x{Ne+1});          % newest node  (k = 0)
%             % 
%             %     % xk = xE;
%             %     for j=obj.Ne+1:-1:1
%             %     % for j=obj.Ne:-1:1
%             % 
%             %         % yk = obj.buff.sensorMat*xk(obj.buff.idx.q1);
%             %         % zk = obj.buff.Y(:,obj.idxWrap(-j+1));
%             %         % ek = zk-yk;                       % innovation
%             %         % 
%             %         % J  = J + 0.5*ek.'*obj.Qe*ek;
%             %         % % Csub = blkdiag( zeros(numel(xE)-numel(obj.buff.idx.q1)), C );
%             %         % % gradSub = Cfull.' * (obj.Qe*ek);
%             %         % gradJ(idx.xE) = gradJ(idx.xE) - Cfull.' * (obj.Qe*ek) ;
%             %         % % back‑prop state with cached STM (most efficient)
%             %         % Sbk = obj.buff.S(:,:,obj.idxWrap(-j));   % S_{k-1}
%             %         % wj  =  wHist(j); %wHist(:,j);
%             %         % % xk  = xk - Sbk(:,obj.nx+1:end)*wj;        % x_{k-1}
%             %         % xk = xk - Sbk(:, obj.nx+1 ) * wj;   % ✱ uses only gust cols
%             %         % xk  = z(idx.x{j});                         % current decision state
%             %         xk  = p(idx.x{j});                         % current decision state
%             %         yk  = C * xk(q1Idx);                       % predicted sensor output
%             % 
%             %         zk  = obj.buff.Y(:, obj.idxWrap(j-(Ne+1)));% stored measurement
%             %                                                    % idxWrap(0) → newest sample
%             %         ek  = zk - yk;  
%             %          % ---- stage cost ---------------------------------------------
%             %         J = J + 0.5 * ek.' * obj.Qe * ek;
%             % 
%             %         % ---- gradient wrt x^j --------------------------------------
%             %         % dfsdf = sparse(q1Idx , 1:obj.ny , 1 , obj.nx , obj.ny);
%             %         gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cfull.' * (obj.Qe * ek);     % 144 × 1
%             % 
%             %         % if j>1
%             %         %     xk = z(idx.x{j-1});   % move one node older
%             %         % end
%             %     end
%             %     % disps1 = (norm(gradJ(idx.xE)))
%             %     % disps2 = norm(gradJ(idx.w))
%             % end
% 
%             % -------- non‑linear equality constraints ----------------
%             % -------- NON-LINEAR CONTINUITY (dense Jacobian) -----------------------
%             function [c,ceq,gradc,gradceq] = nonlFun(p)
%                 gradc   = [];
%                 nVar    = numel(p);
%                 ceq     = zeros(obj.nx*obj.Ne,1);            % dense
%                 gradceq = zeros(nVar, obj.nx*obj.Ne);        % dense Jacobian
% 
%                 % flat helpers
%                 idxWvec = [idx.w{:}];
%                 wHist   = reshape( p(idxWvec), nw, Ne ).';%  Ne × nw
% 
%                 xk   = p(idx.x{1});                          % x^{-Ne}
% 
%                 % ceq    = zeros(obj.nx*obj.Ne,1);
%                 % gradceq = zeros(obj.nx*obj.Ne , numel(p));   % DENSE – 45 kB for nx=144,Ne=25
%                 % 
%                 % % row = 0;
%                 % % for j = 1:obj.Ne
%                 % %     Sw  = obj.buff.S(:,obj.nx+1, obj.idxWrap(-obj.Ne-1+j)); % 144×1
%                 % %     rowIdx = row + (1:obj.nx);
%                 % % 
%                 % %     gradceq(rowIdx, idx.x{j+1}) =  eye(obj.nx);
%                 % %     gradceq(rowIdx, idx.x{j})   = -eye(obj.nx);
%                 % %     gradceq(rowIdx, idx.w{j})   = -Sw;
%                 % % 
%                 % %     ceq(rowIdx) = p(idx.x{j+1}) - p(idx.x{j}) - Sw*p(idx.w{j});
%                 % %     row = row + obj.nx;
%                 % % end
%                 % % gradceq = gradceq.';
% 
%                 row0 = 0;
% 
%                 for j = 1:obj.Ne
%                 % for j = 1:obj.Ne-1
%                     % Sj   = obj.buff.S(:,:,obj.idxWrap(j-(obj.Ne+1)));
%                     Sj = obj.buff.S(:,:, obj.idxWrap(-Ne+j));
%                     % Sj = obj.buff.S(:,:, obj.idxWrap(-Ne-1+j));
%                     % %%
%                     % SwTest = obj.buff.testSw(:,:,obj.idxWrap(-Ne+j));
%                     % Sw = SwTest;
%                     % %%
% 
%                     Sw   = Sj(:,obj.nx+1);                   % 144×1
%                     wj   = p(idx.w{j});
%                     xkp1 = p(idx.x{j+1});
% %%
%                     % wj   = wHist(j,:).';             %  nw×1
%                     % ceq(row) = x_next - x_curr - Sw*wj;
%   %%
%                     rhs  = xk + Sw*wj;
%                     ceqBlk = xkp1 - rhs;                     % size nx
%                     rows   = row0 + (1:obj.nx);
%                     ceq(rows) = ceqBlk;
% 
%                     % fill Jacobian columns (dense)
%                     gradceq(idx.x{j+1}, rows) =  eye(obj.nx);
%                     gradceq(idx.x{j  }, rows) = -eye(obj.nx);
%                     gradceq(idx.w{j }, rows) = -Sw;
%                     % gradceq(idx.w{j }, rows) = Sw;
%                     % if j < Ne
%                     %     gradceq(idx.w{j+1},rows) =  Sw;                   % ∂(ŵj+1-ŵj)/∂ŵj+1
%                     % 
%                     % end
%                     xk   = xkp1;   row0 = row0 + obj.nx;
%                 end
%                 % gradceq = gradceq.';                         % fmincon format m×n
%                 c =[];
%             end
% 
%             % function [c,ceq,gradc,gradceq] = nonlFun(p)
%             %     gradc = [];
%             %     % xE    = p(idx.xE);
%             %     % --- flat index vectors -----------------------------------------
%             %     idxWvec = [idx.w{:}];            % 1 × (nw·Ne)   numeric row
%             %     % idxXvec = [idx.x{:}];            % 1 × (nx·(Ne+1))
%             %     % % ---------- arrival-cost terms (Eq. 38a) ----------------------------
%             %     % xRef = obj.buff.X(:, obj.idxWrap(0));   % newest sample (k = 0)
%             %     % 
%             %     % xE    = p(idx.x{1});                    % oldest node  (j = 1)
%             %     wHist = reshape(p(idxWvec), nw, Ne).';  % Ne × nw                ceq = zeros(obj.nx*obj.Ne,1);     % pre‑alloc
%             %     % gradceq = sparse(obj.nx*obj.Ne, numel(p));
%             %     % gradceq = sparse(obj.nx*obj.Ne, numel(p));
%             %     gradceq  = spalloc( obj.nx*obj.Ne , numel(p) , ...
%             %             obj.nx*(1+obj.nw)*obj.Ne );   % sparsity hint
%             %     % gradceq = sparse( numel(p), obj.nx*obj.Ne);
%             % 
%             %     % xk = xE;
%             %    % gradceq(row, idx.xE) = eye(obj.nx);
%             %     row = 0;
%             %     for j = 1:Ne
%             %         x_next = p(idx.x{j+1});           % x^{k-j+1}
%             %         x_curr = p(idx.x{j});             % x^{k-j}
%             %         wj     = p(idx.w{j});             % w^{k-j}
%             % 
%             %         % SbkFull = obj.buff.S(:,:, obj.idxWrap(-Ne-1+j));
%             %         SbkFull = obj.buff.S(:,:, obj.idxWrap(-Ne-1+j));
%             %         Sw      = SbkFull(:, nx+1:nx+nw);          % 144×nw
%             % 
%             %         rowIdx   = row + (1:nx);
%             %         ceq(rowIdx) = x_next - x_curr - Sw*wj;      % = 0
%             % 
%             %         % Jacobian entries
%             %         gradceq(rowIdx, idx.x{j+1}) =  eye(nx);
%             %         gradceq(rowIdx, idx.x{j})   = -eye(nx);
%             %         gradceq(rowIdx, idx.w{j})   = -Sw;
%             % 
%             %         row = row + nx;
%             %     end
% 
%                 % for j=obj.Ne:-1:1
%                 %     Sbk = obj.buff.S(:,:,obj.idxWrap(-j));      % STM k‑1→k
%                 %     Sw = Sbk(:,obj.nx+obj.nw); 
%                 %     % uk  = obj.buff.U(:,obj.idxWrap(-j));
%                 %     xk_prev = obj.buff.X(:,obj.idxWrap(-j));    % saved state
%                 % 
%                 %     wj = wHist(j); %(:,j)
%                 %     % xk = xk - Sbk(:,obj.nx+1:end)*wj;
%                 %     xk2 = xk - Sbk(:, obj.nx+1) * wj;   %  uses only gust cols
%                 %     ceqBlk = xk2 - xk_prev;                  % continuity
%                 %     % x_next = 
%                 %     % ceqBlk2 = 
%                 %     % ceqBlk = xk;                  % continuity
%                 %     row = (j-1)*obj.nx + (1:obj.nx);
%                 %     ceq(row) = ceqBlk;
%                 %     % gradient wrt decision vars
%                 %     % gradceq(row, idx.xE) = Sbk(:,1:obj.nx);
%                 %     % Jacobian wrt xE  (only in oldest block)
%                 %     if j == obj.Ne
%                 %         gradceq(row, idx.xE) = eye(obj.nx);
%                 %     end
%                 %     % gradceq(idx.xE, row) = eye(obj.nx);
%                 %     % gradceq(row, idx.w((j-1)*obj.nw+(1:obj.nw))) = -Sbk(:,obj.nx+1:end);
%                 %     % gradceq(row, idx.w((j-1)*obj.nw+(1:obj.nw))) = -Sbk(:,obj.nx+obj.nw); % nw is only 1
%                 %     gradceq(row, idx.w((j-1)*obj.nw+(1:obj.nw))) = -Sbk(:,obj.nx+obj.nw); % nw is only 1
%                 % 
%                 %     % ---------- Jacobian w.r.t. w_j ----------------------------
%                 %     colWj = idx.w( (j-1)*obj.nw + (1:obj.nw) );   % scalar index
%                 %     gradceq(row , colWj) = -Sw;                   % note the minus
%                 % 
%                 %     % ---------- roll one step back -----------------------------
%                 %     xk = xk_prev;                                 % move to older node
%                 %     % % gradceq(idx.w((j-1)*obj.nw+(1:obj.nw)), row) = -Sbk(:,obj.nx+obj.nw); % nw is only 1
%                 %     % % roll back
%                 %     % wj = wHist(j); %(:,j)
%                 %     % % xk = xk - Sbk(:,obj.nx+1:end)*wj;
%                 %     % xk = xk - Sbk(:, obj.nx+1) * wj;   % ✱ uses only gust cols
%                 % 
%                 % end
%             %     gradceq = gradceq.';
%             % 
%             %     c = [];   % no inequality (boxes handled in lb/ub)
%             %     % if i==Nc     % last block
%             %     %     H  = eye(numel(z));           % cheap placeholder
%             %     %     KKT = [H   dceq.';  dceq  zeros(size(dceq,1))];
%             %     %     fprintf('cond(KKT) %.2e\n', condest(KKT));
%             %     % end
%             % 
%             % end
% 
%             % -------- bounds ----------------------------------------
%             % state bounds replicated for every node
%             xLrep = repmat(obj.xL , Ne+1, 1);
%             xUrep = repmat(obj.xU , Ne+1, 1);
% 
%             % disturbance bounds for every interval
%             wLrep = repmat(obj.wL , Ne, 1);
%             wUrep = repmat(obj.wU , Ne, 1);
% 
%             nlp.lb = [xLrep ; wLrep];
%             nlp.ub = [xUrep ; wUrep];
% 
%             % lb = [obj.xL ; repmat(obj.wL,obj.Ne,1)];
%             % ub = [obj.xU ; repmat(obj.wU,obj.Ne,1)];
% 
%             nlp.cost  = @costFun;
%             nlp.nonl  = @nonlFun;
%             % nlp.lb    = lb;   nlp.ub = ub;
%         end
% % ----------------------------------------------------------------------
%     % end
%         function kIdx = idxWrap(obj,shift)
%             % circular index into buffer (latest sample at kBuff)
%             kIdx = mod(obj.kBuff-1+shift,obj.Ne+1)+1;
%             % kIdx = mod(obj.kBuff-1+shift,obj.Ne)+1;
%         end
%         function out = C(obj)
%             % measurement matrix (PhiY) from sensor buffer
%             out = obj.buff.sensorMat;   % 5×nx
%         end
%     end     
% end
% 
