classdef nMPC < AeroFlex.ctrl.ControllerBase
%========================================================
%  NON-LINEAR  MODEL-PREDICTIVE  CONTROLLER (multiple-shoot)
%  – matches Artola 2021, Eq.(37)
%  – internal ROMIntegrator used for continuity Jacobian
%  – vanishing-disturbance assumption (first Nc/2 slots)
%========================================================

    properties
        % ---------- timing & horizon ----------------------------------
        dt double
        Ts double
        Nc double                                % # samples in horizon
    
        % ---------- dimensions ----------------------------------------
        nx double
        nu double
        nw double = 1                            % single gust DOF
    
        % ---------- weights / boxes -----------------------------------
        Qc double; Rc double; Pc double
        xL double; xU double
        uL double; uU double
        urL double; urU double
        wL double; wU double                     % bounds on first Nc/2 slots
    
        % ---------- internal model ------------------------------------
        model          % same ROMIntegrator class you use in MHE
        Sprev          % STM at end of last call (warm-start)
    
        % ---------- rolling histories ---------------------------------
        Xhist          % (nx × Nc+1)
        Uhist          % (nu × Nc)
        Whist          % (nw × Nc)   (zeros after Nc/2)
        uPrev
        k  double = 0  % sample counter
    
        % ---------- optimiser state -----------------------------------
        z0  double     % warm start
        H   double     % LBFGS Hessian
        solverOpts
        % ---------- debug ---------------------------------------------
        debug logical = true
        dbg struct
        wHat      double         % last estimated gust at t-k  (scalar nw = 1)
        Nd           double         % # slots with linear decay  (=⌈Nc/2⌉)
        aT           double         % terminal-set radius in (37f)
        uTrim
        xTrim
    end
    %==================================================================
    methods
    %------------------------------------------------------------------
    function obj = nMPC(cfg,beam,aero,base)
        obj@AeroFlex.ctrl.ControllerBase(cfg);
    
        % basic numbers ------------------------------------------------
        obj.dt  = cfg.sim.dt;
        obj.Ts  = cfg.ctrl.Ts;
        obj.Nc  = cfg.ctrl.Nc;
    
        % model identical to estimator's --------------------------------
        obj.model = cfg.modelHandle(cfg,beam,aero,base);
        obj.nx    = size(obj.model.L,1);
        obj.nu    = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    
        % weights / boxes ----------------------------------------------
        obj.Qc = cfg.Qc; obj.Rc = cfg.Rc; obj.Pc = cfg.Pc;
        obj.xL = cfg.xL; obj.xU = cfg.xU;
        obj.uL = cfg.uL; obj.uU = cfg.uU;
        obj.wL = cfg.wL; obj.wU = cfg.wU;
        obj.urL = cfg.urateL; obj.urU = cfg.urateU;
    
        % histories ----------------------------------------------------
        obj.Xhist = zeros(obj.nx,obj.Nc+1);
        obj.Uhist = zeros(obj.nu,obj.Nc);
        obj.Whist = zeros(obj.nw,obj.Nc);        % will stay ≡0 after Nc/2
    
        % initial state / guess ----------------------------------------
        xTrim        = cfg.trim.states;
        obj.xTrim = xTrim;
        obj.uTrim        = cfg.trim.Uinpt;
        obj.Xhist(:) = repmat(xTrim,1,obj.Nc+1);
        % obj.z0       = [repmat(xTrim,obj.Nc+1,1); ...
        %                 zeros(obj.nw*ceil(obj.Nc/2),1); ...
        %                 zeros(obj.nu*obj.Nc,1)];
        obj.z0       = [repmat(xTrim,obj.Nc+1,1); ...
                        zeros(obj.nu*obj.Nc,1)];
        obj.H        = speye(numel(obj.z0));
        obj.Nd   = ceil(obj.Nc/2);
        obj.aT   = cfg.ctrl.aT;          % small positive number (≈10-4 … 10-3)
        obj.wHat = 0;                 % will be overwritten by first MHE call
    
        % STM identity -------------------------------------------------
        obj.Sprev = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
    
        % solver options ----------------------------------------------
        obj.solverOpts = optimoptions('fmincon', ...
            'Algorithm','interior-point', ... % interior-point | sqp
            'SpecifyObjectiveGradient',true, ...
            'SpecifyConstraintGradient',true, ...
            'HessianApproximation','lbfgs', ...
            'Display','iter', ...
            'OptimalityTolerance',1e-5, ...
            'StepTolerance',1e-8);
        obj.solverOpts.CheckGradients = true;
    
        obj.uPrev = obj.uTrim;
        obj.dbg = struct('t', [], 'U', [], 'cont', 0);
    
    end
    %------------------------------------------------------------------
    function [uk, Uinfo] = computeControl(obj,xhat,whatEst,uPrev,t_k)
    
        % --------------------------------------------------------------
        %  xMeas   : state estimate from the MHE   (size nx)
        %  whatEst : gust estimate from the MHE    (size nw = 1)
        %  uPrev   : actuator value actually sent to the plant at t-k
        % --------------------------------------------------------------
        obj.wHat = whatEst;                   % cache latest gust
    
        % ---------- roll local buffers ----------------------------
        if obj.k >= obj.Nc
            obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
            obj.Uhist(:,1:end-1) = obj.Uhist(:,2:end);
            obj.Whist(:,1:end-1) = obj.Whist(:,2:end);
        else
            obj.k = obj.k + 1;
        end
        obj.Xhist(:,end) = xhat;                % newest state sample
        obj.Whist(:,end) = whatEst;              % newest gust sample
        % Uhist(:,end) will be filled **after** optimisation --------
        % ---- 2. until horizon is full, just output copy of last estimate
        if obj.k < obj.Nc
            % obj.Xhist(:,end) = obj.xhat;
            % xhat = obj.xhat;
            % what = obj.what;
            uk = obj.uPrev;
            obj.Uhist(:,end) = obj.uPrev;
            Uinfo = struct('cost',0,'continuity',0,'flag',0,...
                          'uHorizon',zeros(obj.Nc*obj.nu,1));
            % obj.z0(idx.x{obj.k}) = xhat;
            return
        end
    %     % ---------- shift buffers ------------------------------------
    %     if obj.k >= obj.Nc
    %         obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
    %         obj.Uhist(:,1:end-1) = obj.Uhist(:,2:end);
    %         obj.Whist(:,1:end-1) = obj.Whist(:,2:end);
    %     else
    %         obj.k = obj.k + 1;
    %     end
    %     obj.Xhist(:,end) = xhat;          % newest measurement
    %     obj.Uhist(:,end) = uPrev;          % for information only
    
        % ---------- assemble & solve NLP -----------------------------
        nlp = obj.assembleWindow();
        [zOpt,fval,flag] = fmincon(nlp.cost,obj.z0,[],[],[],[], ...
                                nlp.lb,nlp.ub,nlp.nonl,obj.solverOpts);
    
        if flag<=0
            warning('nMPC:solve','fmincon failed – hold last input.');
            uk = obj.uPrev;
            Uinfo.cost       = fval;
            Uinfo.exitflag   = flag;
            % Uinfo.uHorizon   = Uopt';

            % Uinfo.uHorizon   = Uopt;
            Uinfo.uHorizon   = obj.uPrev;
            [~,ceq]         = nlp.nonl(zOpt);
            Uinfo.continuity = norm(ceq);
            return
        end
    
        % ---------- extract & shift warm-start -----------------------
        idx = obj.buildIndexMaps();
        % idxU = idx.u{:};
        uk  = zOpt(idx.u{1});              % control to plant
        obj.z0 = obj.shiftGuess(zOpt);
        obj.Uhist(:,end) = uk;             % store actually sent control
        obj.uPrev = uk;                         % store for ZOH & cost


        Uopt = reshape(zOpt((obj.Nc+1)*obj.nx+1:end),obj.nu,obj.Nc);

        Uinfo.cost       = fval;
        Uinfo.exitflag   = flag;
        % Uinfo.uHorizon   = Uopt';
        Uinfo.uHorizon   = Uopt;
        [~,ceq]         = nlp.nonl(zOpt);
        Uinfo.continuity = norm(ceq);
    
        if obj.debug,  obj.debugPlots(t_k,Uinfo);  end
    
    end
    %------------------------------------------------------------------
    end % public
    %==================================================================
    methods(Access=private)
    function wPred = buildPredictedGust(obj)
        % returns 1×Nc vector   w^0 … w^{Nc-1}
        if obj.wHat == 0
            wPred = zeros(1,obj.Nc);
            return
        end
        w0 = obj.wHat;
        wPred = [ linspace(w0,0,obj.Nd+1) , zeros(1,obj.Nc-obj.Nd-1) ];
    end
    
    %------------------------------------------------------------------
    function idx = buildIndexMaps(obj)
        nx=obj.nx; nu=obj.nu; nw=obj.nw; N=obj.Nc;
        
        idx.x = cell(N+1,1);
        for k=1:N+1, idx.x{k} = (k-1)*nx + (1:nx); end
        % disturbance slots only for k=0…Nc/2-1
        % Nd = ceil(N/2);
        % startW=(N+1)*nx;
        % idx.w = cell(Nd,1);
        % for k=1:Nd, idx.w{k}=startW+(k-1)*nw+(1:nw); end
        % startU = startW + Nd*nw;
        startU = (N+1)*nx;
        idx.u = cell(N,1);
        for k=1:N, idx.u{k}=startU+(k-1)*nu+(1:nu); end
    end
    %------------------------------------------------------------------
    function zShift = shiftGuess(obj,z)
        idx=obj.buildIndexMaps(); 
        nx=obj.nx; nu=obj.nu; nw=obj.nw;
        N=obj.Nc; Nc = N; 
        % Nd=numel(idx.w);
        % idx.x = arrayfun(@(k) k*nx+(1:nx),       0:Nc,   'uni',0);
        % idx.u = arrayfun(@(k) (Nc+1)*nx+k*nu+(1:nu), 0:Nc-1, 'uni',0);
        zShift = zeros(size(z));
        % states -------------------------------------------------------
        for k=1:N, zShift(idx.x{k}) = z(idx.x{k+1}); end
        zShift(idx.x{N+1}) = z(idx.x{N+1});
        % % disturbances (vanish after horizon) -------------------------
        % for k=1:Nd-1, zShift(idx.w{k}) = z(idx.w{k+1}); end
        % zShift(idx.w{Nd}) = 0;             % new right-most slot = 0
        % controls ----------------------------------------------------
        for k=1:N-1, zShift(idx.u{k}) = z(idx.u{k+1}); end
        zShift(idx.u{N}) = z(idx.u{N});    % repeat last control
    end
    %------------------------------------------------------------------
    function nlp = assembleWindow(obj)
        nx=obj.nx; nu=obj.nu; nw=obj.nw; Nc=obj.Nc; Nd=ceil(Nc/2);
        % idx = obj.buildIndexMaps();
        Nint = round(obj.Ts/obj.dt);
        idx.x = arrayfun(@(k) k*nx+(1:nx),       0:Nc,   'uni',0);
        idx.u = arrayfun(@(k) (Nc+1)*nx+k*nu+(1:nu), 0:Nc-1, 'uni',0);
        % nVar  = (Nc+1)*nx + Nc*nu;
    
        % ---------- cached gust profile for this horizon ----------
        wHorz = obj.buildPredictedGust();        % 1×Nc
        nVar  = (Nc+1)*nx + Nc*nu;
    
        % --------------------- COST ----------------------------------
        % function [J,g] = costFun(p)
        %     g = zeros(numel(p),1); J=0;
        %     % loop over stages ----------------------------------------
        %     for k=0:N-1
        %         xk = p(idx.x{k+1});
        %         uk = p(idx.u{k+1});
        %         J = J + 0.5*xk.'*obj.Qc*xk + 0.5*uk.'*obj.Rc*uk;
        %         g(idx.x{k+1}) = obj.Qc*xk;
        %         g(idx.u{k+1}) = obj.Rc*uk;
        %     end
        %     % terminal ------------------------------------------------
        %     xN = p(idx.x{N+1});
        %     J  = J + 0.5*xN.'*obj.Pc*xN;
        %     g(idx.x{N+1}) = g(idx.x{N+1}) + obj.Pc*xN;
        % end
        % ========== COST ==========================================
        function [J,g] = costFun(z)
            g = zeros(nVar,1);  
            J = 0;
    
            % for k = 0:Nc-2
            for k = 0:Nc
            % for k = 0:Nc-1
            % for k = 1:Nc-1
            % for k = 0:Nc
                xk = z(idx.x{k+1});
              
                xRef = obj.xTrim;                % use your reference here
                dX   = xk - xRef;
                % dU   = uk - obj.prevU;
    
                
                % J    = J + 0.5*dX.'*obj.Qc*dX + 0.5*dU.'*obj.Rc*dU;
                J    = J + 0.5*dX.'*obj.Qc*dX;
                if k~=Nc
                    uk = z(idx.u{k+1});
                    dU = uk - obj.uTrim;   % deviation from equilibrium setting
                    g(idx.u{k+1}) = g(idx.u{k+1}) + obj.Rc*dU;
                    J = J+ 0.5*dU.'*obj.Rc*dU;
                end
                g(idx.x{k+1}) = g(idx.x{k+1}) + obj.Qc*dX;

            end
            % ---- terminal cost -----------------------------------
            xN   = z(idx.x{Nc+1});
            % xN   = z(idx.x{Nc});
            dXT  = xN - obj.xTrim;
            J    = J + 0.5*dXT.'*obj.Pc*dXT;
            g(idx.x{Nc+1}) = g(idx.x{Nc+1}) + obj.Pc*dXT;
        end
    
        % ---------------- CONTINUITY --------------------------------
        function [c,ceq,dc,dcEq] = nonlFun(p)
            c=[]; dc=[];
            dc    = zeros(nVar,1);             % gradient column
    
            ceq = zeros(nx*Nc,1);
            dcEq = zeros(nx*Nc,numel(p));
            % dcEq  = spalloc(nVar , nx*Nc , nx*(nx+nu)*Nc);
            % dcEq  = spalloc(nVar , nx*Nc , nx*(nx+nu)*Nc).';
    
            % dcEq = zeros(numel(p),nx*Nc);
    
            row0=0;
            for k=0:Nc-1
                idx_xk   = k*nx     + (1:nx);
                idx_xkp1 = (k+1)*nx + (1:nx);
                idx_uk   = (Nc+1)*nx + k*nu + (1:nu);
                % xk   = p(idx.x{k+1});
                % uk   = p(idx.u{k+1});
                xk   = p(idx_xk);
    
                uk   = p(idx_uk);
                % disturbance mapping --------------------------------
    
                wk   = wHorz(k+1);              % ← exogenous gust
    
                % integrate Ts --------------------------------------
                x = xk;  
                S = [eye(nx), zeros(nx,nw+nu)];
                for m=1:Nint
                    [x,S] = obj.model.step(x,uk,wk,S,true);
                end
                Sx = S(:,1:nx);  
                % Sw = S(:,nx+1);  
                Su = S(:,nx+nw+1:end);
                % %%
                % %Temp
                % % ctr1 = ctrbf(Sx, Su, eye(obj.nx));
                % % [Abar,Bbar,Cbar,T,k] = ctrbf(Sx, S(:,obj.nx+1:end), eye(144));
                % [Abar,Bbar,Cbar,T,k] = ctrbf(Sx, Su, eye(144));
                % 
                % [T, Ut] = schur(Abar);
                % Anew = T.'*Abar*T;
                % Bnew = T.'*Bbar;
                % % TB = schur(Bbar);
                %%
                xPred = x;
                rowIdx = row0+(1:nx);
                rows = rowIdx;
                ceq(rowIdx) = p(idx.x{k+2}) - xPred;
    
              
                dcEq(rowIdx,idx_xkp1) =  eye(nx);
                dcEq(rowIdx,idx_xk) = -Sx;
                dcEq(rowIdx,idx_uk) = -Su;
    
    
                % if k<Nd
                %     dcEq(idx.w{k+1},rowIdx) = -Sw;
                % end
                row0 = row0 + nx;
            end
            dcEq = dcEq.';        % fmincon expects m×n Jacobian
    
         
            % % ---- inequality: terminal set (37f) ---------------------------------
            xN = p( idx.x{Nc+1} );
            % xN = p( idx.x{Nc} );
            xNerr = xN - obj.xTrim;
            term = 0.5 * xNerr.' * obj.Pc * xNerr;
            c     = term - obj.aT;             % c ≤ 0
            % dc    = zeros(nVar,1);             % gradient column
            dc( idx.x{Nc+1} ) = obj.Pc * xNerr;   % ∂c/∂x_N
    
        end
    
        % ========== BOUNDS ========================================
        lbu_sub = [obj.uL; obj.uL; obj.urL; obj.urL];
        ubu_sub = [obj.uU; obj.uU; obj.urU; obj.urU];
        lb = [ repmat(obj.xL ,Nc+1,1);
           repmat(lbu_sub ,Nc  ,1)];
        ub = [ repmat(obj.xU ,Nc+1,1);
           repmat(ubu_sub ,Nc  ,1)];
        
        nlp = struct('cost',@costFun, ...
                 'nonl',@nonlFun, ...
                 'lb',lb,'ub',ub);
        % xLrep = repmat(obj.xL,N+1,1);
        % xUrep = repmat(obj.xU,N+1,1);
        % wLrep = repmat(obj.wL,Nd,1);
        % wUrep = repmat(obj.wU,Nd,1);
        % uLrep = repmat(obj.uL,N,1);
        % uUrep = repmat(obj.uU,N,1);
        % 
        % lb = [xLrep; wLrep; uLrep];
        % ub = [xUrep; wUrep; uUrep];
        % 
        % nlp.cost = @costFun;
        % nlp.nonl = @nonlFun;
        % nlp.lb   = lb;
        % nlp.ub   = ub;
    end
    %----------------------------------------------------------------------
        % function debugPlots(obj,t_k,info)
        %     uH = info.uHorizon;
        %      % ---- initialise persistent store on first call --------------------
        %     if ~isfield(obj.dbg,'t')
        %         obj.dbg.t = [];          % Ns×1   sample times
        %         obj.dbg.U = [];          % Ns×Nc  horizon matrix (rows = samples)
        %     end
        %     if isempty(obj.dbg)
        %         obj.dbg.t=t_k; 
        %         obj.dbg.U=uH(:)';
        %         obj.dbg.cont=norm(info.uHorizon);
        %     else
        %         obj.dbg.t(end+1)=t_k; 
        %         % obj.dbg.U(end+1,1:obj.Nc)=uH(:)';
        %         obj.dbg.U(end+1,1:obj.Nc)=uH.';
        %         obj.dbg.cont(end+1)=info.continuity;
        % 
        %     end
        %     % figure(2001),clf,hold on
        %     % for s=1:obj.Nc, stairs(obj.dbg.t,obj.dbg.U(:,s)); end
        %     % hold off,grid on,xlabel('t [s]'),ylabel('u'),title('control slots'), legend
        %     figure(2001),clf,hold on
        %     for s=1:obj.Nc, stairs(obj.dbg.t,obj.dbg.U(:,s)); end
        %     hold off,grid on,xlabel('t [s]'),ylabel('u'),title('control slots'), legend
        % 
        %     figure(2004),clf, hold on,plot(obj.dbg.t,obj.dbg.cont(2:end)); hold off;grid on
        %     xlabel('t [s]'); ylabel('|continuity|'); title('MHE feasibility'), legend
        % end
        %------------------------------------------------------------------
        function debugPlots(obj,t_k,info)
        % Shows how the Nc control moves (uHorizon) evolve over time.
        % Each row in dbg.U is the **flattened** horizon at one sample.
        %------------------------------------------------------------------
            uH = info.uHorizon;               % size = (nu × Nc)
        
            % ---------- initialise the persistent store ------------------
            if ~isfield(obj.dbg,'t')
                obj.dbg.t   = [];                         % 1 × Ns
                obj.dbg.U   = [];                         % Ns × (nu*Nc)
                obj.dbg.cont= [];                         % 1 × Ns-1
            end
        
            % ---------- append this sample --------------------------------
            obj.dbg.t(end+1)    = t_k;
            obj.dbg.U(end+1,:)  = uH(:).';                % flatten → 1 × (nu*Nc)
            obj.dbg.cont(end+1) = info.continuity;
        
            % ---------- plotting ------------------------------------------
            figure(3001); clf
            for ii = 1:obj.nu          % one subplot per actuator
                subplot(obj.nu,1,ii); hold on
                cols = ii:obj.nu:size(obj.dbg.U,2);       % pick horizon slots for actuator-ii
                for s = 1:obj.Nc
                    stairs(obj.dbg.t , obj.dbg.U(:,cols(s)) , 'LineWidth',1.0);
                end
                hold off, grid on,legend
                ylabel(sprintf('u_%d [-]',ii));
                if ii==1, title('Control horizons (oldest → left)'); end
            end
            xlabel('time  [s]'); 
            % drawnow;
        
            % ---------- continuity trace ----------------------------------
            figure(3004); clf
            % plot( obj.dbg.t(2:end) , obj.dbg.cont(2:end) , 'LineWidth',1.0 );
            plot( obj.dbg.t , obj.dbg.cont(2:end) , 'LineWidth',1.0 );
            grid on; xlabel('time [s]'); ylabel('|continuity|');
            title('MPC feasibility');
        end
        %------------------------------------------------------------------
    %------------------------------------------------------------------
    end % private
end % class

% %% Old BUFF nMPC method
% % control/@nMPC/nMPC.m
% classdef nMPC < AeroFlex.ctrl.ControllerBase
% % Non‑linear MPC with multiple shooting (Artola 2021, Eq.(37))
% %  • prediction horizon  T_c = Nc·Ts
% %  • cost: stage (Qc,Rc) + terminal (Pc)
% %  • Actuator and state boxes from cfg.uL/uU, xL/xU
% 
%     properties
%         dt, Ts, Nc                              % timing
%         Qc, Rc, Pc                              % weights
%         xL, xU, uL, uU, urL, urU                % boxes
%         nx, nw, nu                                  % sizes
%         buff                                    % ↯ shared with SimRunner
%         kBuff = 0                               % circular pointer
% 
%         solverOpts                              % fmincon options
%         H                                       % Hessian (for SQP)
% 
%         z0                                      % warm start [x0 ; vec(u)]
%         % prevU                                   % last control for ZOH
%         ref                                     % trim references
%         xTrim
%     end
% % ---------------------------------------------------------------------
%     methods
%         function obj = nMPC(cfg)
%             obj@AeroFlex.ctrl.ControllerBase(cfg);
%             % ---------- copy scalars from cfg ---------------------------------
%             obj.dt  = cfg.sim.dt;
%             obj.Ts  = cfg.ctrl.Ts;
%             obj.Nc  = cfg.ctrl.Nc;
% 
%             obj.Qc  = cfg.Qc;   obj.Rc = cfg.Rc;   obj.Pc = cfg.Pc;
%             obj.xL  = cfg.xL;   obj.xU = cfg.xU;
%             obj.uL  = cfg.uL;   obj.uU = cfg.uU;
%             obj.urL = cfg.urateL; obj.urU = cfg.urateU;
%             obj.buff = cfg.buff;                        % ring buffers
%             obj.nx   = size(obj.buff.X,1);
%             obj.nu   = size(obj.buff.U,1);
%             obj.nw = 1;
%             % obj.ref  = cfg.ref;
%             obj.xTrim = cfg.ref.states;
%             obj.ref  = cfg.ref.deltaDeg*ones(4,1);
%             cfg.trim.u = cfg.ref.deltaDeg*ones(4,1);
% 
% 
%             % obj.prevU = zeros(obj.nu,1);                % start with zero actuation
% 
%             % ---------- fmincon (SQP) options ----------------------------------
%             % obj.solverOpts = optimoptions('fmincon',...
%             %     'Algorithm','sqp', ...
%             %     'SpecifyObjectiveGradient',true,...
%             %     'SpecifyConstraintGradient',true,...
%             %     'MaxIterations',50,...
%             %     'Display','iter-detailed', ...
%             %     'checkGradients', false);
%             obj.solverOpts = optimoptions('fmincon',...
%                 'Algorithm','interior-point', ...
%                 'SpecifyObjectiveGradient',true,...
%                 'SpecifyConstraintGradient',true,...
%                 'MaxIterations',200,...
%                 'HessianApproximation','lbfgs', ...
%                 'Display','off', ...  % iter-detailed
%                 'checkGradients', false);
% 
%              % obj.solverOpts.StepTolerance        = 1e-8;
% 
%              % obj.solverOpts.FiniteDifferenceType ='central';
%              % obj.solverOpts.FunValCheck = 'on';
%              % obj.solverOpts.StepTolerance = 1e-5;
%              % obj.solverOpts.ConstraintTolerance = 1e-5;
%              obj.solverOpts.OptimalityTolerance = 1e-5;
%              % obj.solverOpts.OptimalityTolerance = 1e-4;
% 
% 
%            % ---------- warm-start vector  [x0 ; u(-Nc)…u(-1)] ---------------
%             % obj.z0 = [cfg.ref.states ; repmat(obj.prevU , obj.Nc,1)];
%             % initial Hessian for custom SQP
%             nv   = numel(obj.z0);
%             obj.H = speye(nv);
%         end
% % ---------------------------------------------------------------------
%         function uk = computeControl(obj, xhat, t, buff)
%         % called once every sampling period by SimRunner
%         % -------------------------------------------------------------
%             obj.buff = buff;
%             obj.z0 = [repmat(xhat, obj.Nc+1, 1) ; repmat(obj.prevU , obj.Nc,1)];
% 
%             % obj.kBuff = obj.kBuff + 1;              % circular pointer
%             obj.kBuff = mod(obj.kBuff , obj.Nc) + 1;
%             obj.buff.valid = min(obj.buff.valid+1 , obj.Nc+1);
% 
%             % --------------- warm-up (fill buffer first) ------------------
%             if obj.buff.valid < obj.Nc+1
%                 uk = obj.prevU;    return
%             end
% 
%             % build & solve window NLP
%             % nlp = obj.assembleWindow(xhat);                 
%             nlp = obj.assembleWindow_full(xhat);   % Currently runs?
%             %%
%             [zOpt,~,exitFlag] = fmincon(nlp.cost, obj.z0, ...
%                                  [],[],[],[], nlp.lb, nlp.ub, ...
%                                  nlp.nonl, obj.solverOpts);
% 
%             % ---- custom SQP, comment-in when stable ----------
%             % sqpOpts.MaxIter = 80;  sqpOpts.Display = true;
%             % [zOpt,obj.H,info] = AeroFlex.optim.runSQP(nlp,obj.z0,obj.H,sqpOpts);
%             % exitFlag = info.flag;
% 
%             if exitFlag<=0
%                 warning('nMPC','fmincon failed – holding previous control.')
%                 uk = obj.prevU;
%             else
%                 obj.z0 = zOpt;                                  % warm start
%                 uk      = zOpt(obj.nx + (1:obj.nu));            % first move
%                 obj.prevU = uk;
%             end
%         end
%     end
% % =====================================================================
%     methods%(Access=private)
%         function nlp = assembleWindow(obj,x0)
%         % build cost / constraints from shared buffers & current x0
%         % -----------------------------------------------------------------
%         % decision vector  z = [ x0 ; u⁰ ; u¹ ; … ; u^{Nc-1} ]
%         % -----------------------------------------------------------------
%         nx=obj.nx; nu=obj.nu; Nc=obj.Nc; buff=obj.buff;
% 
%         % if strcmp(obj.mode,'full')        % ---------- FULL MS ----------
%             % 1. index maps ------------------------------------------------
%             idx.x = cell(Nc+1,1);
%             for k = 1:Nc+1
%                 idx.x{k} = (k-1)*nx + (1:nx);
%             end
%             idx.u = cell(Nc,1);
%             for k = 1:Nc
%                 idx.u{k} = (Nc+1)*nx + (k-1)*nu + (1:nu);
%             end
%             nVar   = (Nc+1)*nx + Nc*nu;
% 
%             % 2. cost ------------------------------------------------------
%             function [J,g] = costFun(p)
%                 J = 0;                     % accumulate
%                 g = zeros(nVar,1);
% 
%                 xk = p(idx.x{1});          % first state
%                 for k = 1:Nc                   % stage loop
%                     uk   = p(idx.u{k});
%                     % xRef = buff.X(:,obj.idxWrap(-k+1));      % “preview” ref.
%                     xRef = obj.xTrim;
%                     J    = J + 0.5*(xk-xRef).'*obj.Qc*(xk-xRef) ...
%                                + 0.5*(uk-obj.prevU).'*obj.Rc*(uk-obj.prevU);
% 
%                     g(idx.x{k}) = g(idx.x{k}) + obj.Qc*(xk-xRef);
%                     g(idx.u{k}) = obj.Rc*(uk-obj.prevU);
% 
%                     xk = p(idx.x{k+1});                    % next free state
%                 end
%                 % terminal weight
%                 xRef = buff.X(:,obj.idxWrap(0));
%                 J    = J + 0.5*(xk-xRef).'*obj.Pc*(xk-xRef);
%                 g(idx.x{Nc+1}) = g(idx.x{Nc+1}) + obj.Pc*(xk-xRef);
%             end
% 
%             % 3. continuity constraints  x^{k+1}-f(x^k,u^k)=0 -------------
%             function [c,ceq,dc,dcEq] = nonlFun(p)
%                 c   = [];         dc  = [];      % no inequalities
%                 ceq = zeros(nx*Nc,1);
%                 dcEq = spalloc(nx*Nc,nVar,nx*(nx+nu)*Nc);
% 
%                 row = 0;
%                 for k = 1:Nc
%                     xk   = p(idx.x{k});
%                     xkp1 = p(idx.x{k+1});
%                     uk   = p(idx.u{k});
% 
%                     % --- use pre-integrated STM (buff.S) -----------------
%                     Sk = buff.S(:,:, obj.idxWrap(-(k-1)));   % k-th STM
%                     dx = xkp1 - ( xk + Sk(:,nx+1+1:nx+1+nu)*(uk-obj.prevU) );
% 
%                     rIdx      = row + (1:nx);
%                     ceq(rIdx) = dx;
% 
%                     % Jacobians
%                     dcEq(rIdx,idx.x{k+1}) =  eye(nx);
%                     dcEq(rIdx,idx.x{k})   = -eye(nx);
%                     dcEq(rIdx,idx.u{k})   = -Sk(:,nx+1+1:nx+1+nu);
% 
%                     row = row + nx;
%                 end
% 
%             dcEq = dcEq.';     % fmincon expects Jacobian transposed
%             end
%              % 4. bounds ----------------------------------------------------
%             lbu_sub = [obj.uL; obj.uL; obj.urL; obj.urL];
%             ubu_sub = [obj.uU; obj.uU; obj.urU; obj.urU];
%             % ---- bounds ---------------------------------------------
% 
%             lb = [ repmat(obj.xL ,Nc+1,1) ;
%                    repmat(lbu_sub ,Nc  ,1) ];
%             ub = [ repmat(obj.xU ,Nc+1,1) ;
%                    repmat(ubu_sub ,Nc  ,1) ];
% 
%             nlp.cost  = @costFun;
%             nlp.nonl  = @nonlFun;
%             nlp.lb    = lb;
%             nlp.ub    = ub;
% 
%             % --------- warm-start shift ----------------------------------
%             %   [x^0 … x^{Nc}  u^0 … u^{Nc-2}]  ← shift left
%             xBlock = reshape(obj.z0( [idx.x{:}] ), nx,[]);   % (nx × Nc+1)
%             uBlock = reshape(obj.z0( [idx.u{:}] ), nu,[]);   % (nu × Nc)
%             obj.z0 = [ xBlock(:,2:end);        % drop oldest
%                        xBlock(:,end);          % repeat last state
%                        uBlock(:,2:end);        % shift controls
%                        uBlock(:,end) ];    
% 
%         % else
% 
%         %     nx=obj.nx; nu=obj.nu; Nc=obj.Nc; buff=obj.buff;
%         % 
%         %     % decision vector z = [x0 ; uSeq(:)]
%         %     idx.x0 = 1:nx;
%         %     % idx.u  = nx + (1:nu*Nc);
%         %     idx.u  = cell(Nc,1);
%         %     for k = 1:Nc
%         %         idx.u{k} = nx + (k-1)*nu + (1:nu);
%         %     end
%         %     idx.uVec = [idx.u{:}];
%         % 
%         %     % ---- cost ------------------------------------------------
%         %     function [J,g] = costFun(z)
%         %         xk   = z(idx.x0);                          % state at node 0
%         %         uSeq = reshape(z(idx.uVec),nu,[]);         % nu × Nc
%         %         g    = zeros(size(z));   J = 0;
%         % 
%         %         for i=1:obj.Nc
%         %             Qi = obj.Qc;  Ri = obj.Rc;
%         %             refX = buff.X(:, obj.idxWrap(i-1));    % use latest ref
%         %             dX   = xk - refX;
%         %             dU   = uSeq(:,i) - obj.prevU;
%         % 
%         %             J = J + 0.5*dX.'*Qi*dX + 0.5*dU.'*Ri*dU;
%         % 
%         %             % gradients
%         %             g(idx.x0)              = g(idx.x0) + Qi*dX;
%         %             g(idx.u{i})            =           Ri*dU;
%         % 
%         %             % ---- propagate with STM columns for **u** ---------
%         %             ScolU = buff.S(:,:, obj.idxWrap(i-1));
%         %             S_ctrl = ScolU(:, nx+1+obj.nw : end);   % 144×4
%         %             xk = xk + S_ctrl*dU;   
%         %         end
%         %         % ---- terminal cost -----------------------------------
%         %         refX = buff.X(:, obj.idxWrap(obj.Nc));         % ref at end
%         %         dXT  = xk - refX;
%         %         J    = J + 0.5*dXT.'*obj.Pc*dXT;
%         %         g(idx.x0) = g(idx.x0) + obj.Pc*dXT;        % ∂terminal/∂x0
%         %     end
%         % 
%         %     % ---- continuity constraints (multiple shooting) ----------
%         %     function [c,ceq,gradC,gradceq] = nonlFun(z)
%         %         xk   = z(idx.x0);
%         %         uSeq = reshape(z(idx.uVec),nu,[]);
%         %         ceq  = zeros(nx*Nc,1);
%         %         gradceq = spalloc(nx*Nc , numel(z) , nx*(1+nu)*Nc);
%         %         gradC =[];
%         %         for i = 1:Nc
%         %             ScolU   = buff.S(:,:, obj.idxWrap(i-1));
%         %             S_ctrl  = ScolU(:, nx+1+obj.nw : end);   % 144×4
%         %             dU      = uSeq(:,i) - obj.prevU;
%         % 
%         %             xPred   = xk + S_ctrl*dU;               % node i
%         %             row     = (i-1)*nx + (1:nx);
%         %             ceq(row)= xPred - buff.X(:,obj.idxWrap(i));
%         % 
%         %             % Jacobians
%         %             gradceq(row, idx.x0)   = eye(nx);
%         %             gradceq(row, idx.u{i}) = S_ctrl;
%         % 
%         %             xk = xPred;                              % shift
%         %         end
%         %         c = [];
%         %         gradceq = gradceq.';                         % to (nVar×m)
%         %     end
%         %     lbu_sub = [obj.uL; obj.uL; obj.urL; obj.urL];
%         %     ubu_sub = [obj.uU; obj.uU; obj.urU; obj.urU];
%         %     % ---- bounds ---------------------------------------------
%         %     lb = [obj.xL ; repmat(lbu_sub,obj.Nc,1)];
%         %     ub = [obj.xU ; repmat(ubu_sub, obj.Nc,1)];
%         % 
%         %     nlp.cost  = @costFun;
%         %     nlp.nonl  = @nonlFun;
%         %     nlp.lb    = lb;
%         %     nlp.ub    = ub;
%         % % end
%         % end
%         end
% 
%         function nlp = assembleWindow_full(obj,x0)
%             % variables: z = [ x0; x1; … ; xNc ; u0; … ; uNc-1 ]
%             nx=obj.nx; nu=obj.nu; Nc=obj.Nc;
% 
%             % block indices ----------------------------------------------------
%             idx.x  = arrayfun(@(k) (k*nx)+(1:nx),       0:Nc, 'uni',0);
%             idx.u  = arrayfun(@(k) ((Nc+1)*nx)+k*nu+(1:nu), 0:Nc-1,'uni',0);
% 
%             % cost -------------------------------------------------------------
%             function [J,g] = costFun(z)
%                 g = zeros(numel(z),1);   J = 0;
%                 for k=0:Nc-1
%                     xk = z(idx.x{k+1});            % state k
%                     uk = z(idx.u{k+1});            % control k
%                     xRef = obj.buff.X(:,obj.idxWrap(-Nc+k));     % measured/target
%                     J  = J + 0.5*(xk-xRef).'*obj.Qc*(xk-xRef) ...
%                             + 0.5*(uk-obj.prevU).'*obj.Rc*(uk-obj.prevU);
%                     g(idx.x{k+1}) = g(idx.x{k+1}) + obj.Qc*(xk-xRef);
%                     g(idx.u{k+1}) = g(idx.u{k+1}) + obj.Rc*(uk-obj.prevU);
%                 end
%                 % terminal
%                 xN  = z(idx.x{Nc+1});
%                 xRef= obj.buff.X(:,obj.idxWrap(0));
%                 J   = J + 0.5*(xN-xRef).'*obj.Pc*(xN-xRef);
%                 g(idx.x{Nc+1}) = g(idx.x{Nc+1}) + obj.Pc*(xN-xRef);
%             end
% 
%             % continuity -------------------------------------------------------
%             function [c,ceq,gradC,gradc] = nonlFun(z)
%                 gradC = [];
%                 ceq  = zeros(nx*Nc,1);
%                 gradc= spalloc(nx*Nc,numel(z),nx*(2+nu)*Nc);
%                 for k=0:Nc-1
%                     xk   = z(idx.x{k+1});
%                     uk   = z(idx.u{k+1});
%                     STM  = obj.buff.S(:,:,obj.idxWrap(-Nc+k)); % contains Φx,Φu
%                     xkp1 = z(idx.x{k+2});
%                     model = xk + STM(:,nx+1+1:nx+1+nu)*(uk-obj.prevU); % Φx≈I by construction
%                     row = k*nx + (1:nx);
%                     ceq(row) = xkp1 - model;
%                     gradc(row,idx.x{k+2}) = eye(nx);
%                     gradc(row,idx.x{k+1}) = -eye(nx);
%                     gradc(row,idx.u{k+1}) = -STM(:,nx+1+1:nx+1+nu);
%                 end
%                 c=[]; gradc = gradc.';
%             end
% 
%             % bounds -----------------------------------------------------------
%             lbu_sub = [obj.uL; obj.uL; obj.urL; obj.urL];
%             ubu_sub = [obj.uU; obj.uU; obj.urU; obj.urU];
%             % ---- bounds ---------------------------------------------
% 
%             lb = [ repmat(obj.xL ,Nc+1,1) ;
%                    repmat(lbu_sub ,Nc  ,1) ];
%             ub = [ repmat(obj.xU ,Nc+1,1) ;
%                    repmat(ubu_sub ,Nc  ,1) ];
% 
%             % lb = [repmat(obj.xL,Nc+1,1); repmat(obj.uL,Nc,1)];
%             % ub = [repmat(obj.xU,Nc+1,1); repmat(obj.uU,Nc,1)];
% 
%             nlp = struct('cost',@costFun,'nonl',@nonlFun,'lb',lb,'ub',ub);
%             end
% 
%         % ---------------------------------------------------------------------
%         function kIdx = idxWrap(obj,shift)
%             % circular pointer into ring buffer      (latest = 0)
%             kIdx = mod(obj.kBuff-1+shift , obj.Nc) + 1;
%         end
%     end
% end


