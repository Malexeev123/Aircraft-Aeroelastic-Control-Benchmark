classdef nMHEv2 < AeroFlex.ctrl.EstimatorBase
%======================================================================
%  NON‑LINEAR MOVING–HORIZON ESTIMATOR (buffer–free, RT‑integrator)
%  -------------------------------------------------------------------
%  • Multiple‑shooting   *or*   single‑shooting  (cfg.ctrl.mhe_method)
%  • No external ring‑buffer:   past X,U,Y are stored locally.
%  • Uses its own ROMIntegrator instance for continuity Jacobians.
%  • Warm‑start shift (RTI) implemented exactly as Diehl 2005.
%======================================================================

    properties %(SetAccess = private)
        % ------------ configuration scalars --------------------------
        dt   double            % integrator step [s]
        Ts   double            % sample period   [s]
        Ne   double            % horizon length  [samples]
        method (1,1) string {mustBeMember(method,["multiple","single"])} = "multiple" % multiple, single
        
        % ------------ dimensions -------------------------------------
        nx  double
        nu  double
        ny  double
        nw  double

        % ------------ weights / bounds -------------------------------
        Qe  double
        Pe  double
        Re  double
        xL  double;  xU double
        wL  double;  wU double

        % ------------ model & data -----------------------------------
        model            % handle → ROMIntegrator (private copy)
        sensor
        Sprev
        %  local rolling buffers (chronological)
        Xhist            % (nx×Ne+1)
        Uhist            % (nu×Ne  )
        Yhist            % (ny×Ne+1)
        Whist            % (nw×Ne  )
        % index helpers
        k                 % current sample counter (starts at 0)

        % ------------ optimiser state --------------------------------
        z0   double      % warm‑start vector
        H    double      % Hessian approximation (LBFGS)
        solverOpts       % fmincon options

        % ------------ outputs ----------------------------------------
        % xhat double
        % what double

        % ------------ debug store ------------------------------------
        debug logical = true
        dbg   struct
    end

%======================================================================
methods
%--------------------------------------------------------------------
function obj = nMHEv2(cfg,beam,aero,base)
    arguments
        cfg struct
        beam 
        aero 
        base 
    end
    obj@AeroFlex.ctrl.EstimatorBase(cfg);

    % ---------- scalars ---------------------------------------------
    obj.dt  = cfg.sim.dt;           obj.Ts = cfg.ctrl.Ts;
    obj.Ne  = cfg.ctrl.Ne;          obj.method = lower(string(cfg.ctrl.mhe_method));

    % ---------- dimensions ------------------------------------------
    obj.model = cfg.modelHandle(cfg,beam,aero,base);      % ROM copy
    obj.sensor = cfg.sensorHandle(beam, cfg); 
    obj.nx = size(obj.model.L,1);
    obj.nu = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    obj.ny = cfg.ny;          obj.nw = cfg.nw;

    % ---------- weights / bounds ------------------------------------
    obj.Qe = cfg.Qe;   obj.Pe = cfg.Pe;   obj.Re = cfg.Re;
    obj.xL = cfg.xL;   obj.xU = cfg.xU;   obj.wL = cfg.wL; obj.wU = cfg.wU;

    % ---------- optimiser opts --------------------------------------
    % interior-point, true
    obj.solverOpts = optimoptions('fmincon', 'Algorithm','interior-point ', ...   
        'SpecifyObjectiveGradient',true, ...
        'SpecifyConstraintGradient',true, ...
        'HessianApproximation','lbfgs','Display','iter', ...
        'OptimalityTolerance',1e-5, 'StepTolerance',1e-8);

    % obj.solverOpts.FiniteDifferenceType ='central';
    obj.solverOpts.FiniteDifferenceStepSize = 1e-8;
    % obj.solverOpts.EnableFeasibilityMode = true;
    obj.solverOpts.CheckGradients = true;

    % ---------- local buffers ---------------------------------------
    obj.Xhist = zeros(obj.nx,obj.Ne+1);
    obj.Uhist = zeros(obj.nu,obj.Ne);
    obj.Yhist = zeros(obj.ny,obj.Ne+1);
    obj.Whist = zeros(obj.nw,obj.Ne);
    obj.k = 0;

    % ---------- warm‑start ------------------------------------------
    idx = obj.buildIndexMaps();

    obj.z0 = zeros((obj.Ne+1)*obj.nx + obj.Ne*obj.nw,1);
    obj.H  = speye(numel(obj.z0));

    % obj.xhat = zeros(obj.nx,1);
    obj.what = zeros(obj.nw*obj.Ne,1);
    xTrim = cfg.trim.states;              % vector nx×1
    obj.xhat = xTrim;                     % <-- initial state estimate
    obj.Xhist(:,1:end) = repmat(xTrim,1,obj.Ne+1);

    % obj.z0(idx.x{1}) = xTrim;
    % obj.z0(1:obj.nx*(obj.Ne+1)) = repmat(xTrim,obj.Ne+1,1);
    % 
    % set a neutral disturbance history (usually zero)
    obj.Whist(:, :) = 0;                  % or cfg.trim.w_r if you have one

    S0 = [eye(numel( obj.xhat)) , zeros(numel(obj.xhat),1+obj.nu)];
    obj.Sprev  = S0;
    obj.dbg = struct('t', [], 'W', [], 'cont', 0);
end
%--------------------------------------------------------------------
function [xhat,what,info] = estimate(obj,z_k,u_prev,t_k)
%   Called exactly once per Ts from the real‑time loop.
%--------------------------------------------------------------------
    % ---------- roll buffers ----------------------------------------
    if obj.k>=obj.Ne          % buffer full ⇒ shift left one slot 
        obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
        obj.Yhist(:,1:end-1) = obj.Yhist(:,2:end);
        obj.Whist(:,1:end-1)= obj.Whist(:,2:end);
        obj.Uhist(:,1:end-1)= obj.Uhist(:,2:end);
    else
        obj.k = obj.k+1;
    end

    % store placeholder for newest slot (x_k unknown yet) ------------
    obj.Yhist(:,end) = z_k;
    obj.Uhist(:,end) = u_prev;

    % ---------- early exit until buffer full ------------------------
    if obj.k < obj.Ne
        % g_k = cfg.gust_profile(t_k);          % true gust at this instant
        % [obj.xhat,~] = obj.model.step(obj.xhat,u_prev,g_k,[],[]);
    
        obj.Xhist(:,end) = obj.xhat;          % update newest slot
        obj.Whist(:,end) = obj.what(end-obj.nw+1:end);   % newest slot only
        xhat = obj.xhat; 
        what=obj.what; 
        info=struct('cost',0,'continuity',0,'exitflag',0,'wHorizon',zeros(obj.Ne,1));
        return
    end
    idx = obj.buildIndexMaps();

    % 1. arrival state x₀ := current best estimate
    obj.z0(idx.x{1}) = obj.xhat;          
    
    % 2. forward-predict states x₁ … xᴺ with zero gust
    x_tmp = obj.xhat;
    S_tmp = obj.Sprev;                    % identity STM
    for j = 1:obj.Ne
        uj = obj.Uhist(:,j);              % control during slot j (currently zeros)
        wj = 0;                           % initial gust guess
        for m = 1:round(obj.Ts/obj.dt)
            [x_tmp,S_tmp] = obj.model.step(x_tmp, uj, wj, S_tmp, false);
        end
        obj.z0(idx.x{j+1}) = x_tmp;
    end
    % ---------- build & solve NLP -----------------------------------
    nlp = obj.assembleWindow();
    [p,fval,exitflag] = fmincon(nlp.cost,obj.z0,[],[],[],[],nlp.lb,nlp.ub,nlp.nonl,obj.solverOpts);

    % ---------- shift warm‑start ------------------------------------
    obj.z0 = obj.shiftGuess(p);

    % ---------- unpack solution -------------------------------------
    idx = obj.buildIndexMaps();
    obj.xhat = p(idx.x{end});
    obj.what = p([idx.w{:}]);
    xhat = obj.xhat;  what=obj.what;

    % ---------- update latest history slot with estimate ------------
    obj.Xhist(:,end) = xhat;
    obj.Whist(:,end) = what(end-obj.nw+1:end);
    obj.Sprev = nlp.getSTM();

    % ---------- diagnostics ----------------------------------------
    info.cost=fval; info.exitflag=exitflag; info.wHorizon=reshape(obj.what,obj.nw,[])';
    [~,ceq] = nlp.nonl(p); info.continuity=norm(ceq);

    if obj.debug, obj.debugPlots(t_k,info); end
    % obj.Sprev = S_n;
end
%--------------------------------------------------------------------
end  % public methods

%======================================================================
methods (Access = private)
%--------------------------------------------------------------------
function idx = buildIndexMaps(obj)
    nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
    idx.x=cell(Ne+1,1); for j=1:Ne+1, idx.x{j}=(j-1)*nx+(1:nx); end
    startW=(Ne+1)*nx; idx.w=cell(Ne,1);
    for j=1:Ne, idx.w{j}=startW+(j-1)*nw+(1:nw); end
end
%--------------------------------------------------------------------
function zShift = shiftGuess(obj,p)
    idx = obj.buildIndexMaps(); nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
    zShift=zeros(size(p));
    for j=1:Ne, zShift(idx.x{j}) = p(idx.x{j+1}); end
    zShift(idx.x{Ne+1}) = p(idx.x{Ne+1});
    for j=1:Ne-1, zShift(idx.w{j}) = p(idx.w{j+1}); end
    zShift(idx.w{Ne})   = p(idx.w{Ne});
end
%--------------------------------------------------------------------
function nlp = assembleWindow(obj)
    nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
    idx = obj.buildIndexMaps(); idxWvec=[idx.w{:}];
    Nint=round(obj.Ts/obj.dt);
    S_last = obj.Sprev; 
    % ------------------------------------------------ COST/GRAD ------
    function [J,gradJ] = costFun(p)
        gradJ = zeros(numel(p),1); 
        J=0;
        % arrival -----------------------------------------------------
        xE = p(idx.x{1}); 
        dx = xE-obj.Xhist(:,1);
        J = 0.5*dx'*obj.Pe*dx; 
        gradJ(idx.x{1}) = obj.Pe*dx;
        
        wE = p(idx.w{1}); 
        dw = wE-obj.Whist(:,1);
        J = J+0.5*dw'*obj.Re*dw; 
        gradJ(idx.w{1}) = obj.Re*dw;
        % measurement innovations ------------------------------------
        C=obj.sensor.PhiY; 
        Cfull=zeros(obj.ny,obj.nx);
        Cfull(:,obj.model.idx.q1)=C;
        Xmat=reshape(p([idx.x{:}]),nx,[]);
        for j=1:Ne+1
            e=obj.Yhist(:,j)-C*Xmat(obj.model.idx.q1,j);
            J=J+0.5*e'*obj.Qe*e;
            gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cfull'*(obj.Qe*e);
        end
    end
%% V1
    % ------------------------------------------------ CONTINUITY -----
    % function [c,ceq,gradc,gradceq] = nonlFun(p)
    %     if obj.method=="single"
    %         c=[]; ceq=[]; gradc=[]; gradceq=[]; return
    %     end
    %     gradc=[];
    %     ceq=zeros(nx*Ne,1); gradceq=zeros(numel(p),nx*Ne);
    %     wHist=reshape(p(idxWvec),nw,Ne).'; row0=0;
    %     x_prev=p(idx.x{1});
    %     Nint  = round(obj.Ts / obj.dt);          % <— one-time constant
    % 
    %     S_n   = S_last;                       % STM at start of interval
    %     x_prv = p(idx.x{1});
    %     Dw = zeros(nx, nw*Ne);          % derivative block  [initialised once]
    % 
    %     for j = 1:Ne
    %     % for j = 1:Ne-1
    %         wj   = wHist(j);
    %         uj   = obj.Uhist(:,j);
    %         x_k = x_prv;
    %         % integrate across the full Ts interval ---------------------------
    %         for m = 1:Nint
    %             [x_prv,S_n] = obj.model.step(x_prv , uj , wj , S_n , true);
    % 
    %             % Dw = S_n(:,1:nx) * Dw;          % Φ_x * existing
    %             % col_j = (j-1)*nw + (1:nw);
    %             % Dw(:,col_j) = Dw(:,col_j) + S_n(:,nx+1);   % add current Φ_w
    %         end
    % 
    %         Sw     = S_n(:, obj.nx+1            );   % gust col
    %         Su     = S_n(:, obj.nx+2 : obj.nx+1+obj.nu);
    % 
    %         rhs  = x_k + Sw*wj + Su*uj;                % f(xk,uj,wj)
    % 
    %         rows = row0+(1:nx);
    %         x_next = x_prv;                      % now at k-j+1
    % 
    %         ceq(rows) = p(idx.x{j+1}) - x_next;
    %         % ceq(rows) = p(idx.x{j+1}) - rhs;
    % 
    %         gradceq(idx.x{j},  rows) = -eye(nx);
    %         gradceq(idx.x{j+1},rows) =  eye(nx); 
    %         if j <Ne 
    %             gradceq(idx.w{j},  rows) = -Sw;
    %             gradceq(idx.w{j+1},  rows) = Sw;
    %         end
    %         % rows = row0+(1:nx);
    %         % gradceq(idx.w{j},rows) = -Dw(:, (j-1)*nw+(1:nw));   % direct
    %         % gradceq(idx.w{1:j-1},rows) = -Dw(:,1:(j-1)*nw);     % indirect
    %         row0  = row0+nx;
    %     end
  %% V-2
    % function [c,ceq,gradc,gradceq] = nonlFun(p)
    %     if obj.method=="single"
    %         c=[]; ceq=[]; gradc=[]; gradceq=[]; return;
    %     end
    %     gradc=[];
    %     ceq=zeros(nx*Ne,1);
    %     gradceq=zeros(numel(p),nx*Ne);
    %     % reshape history
    %     wHist=reshape(p(idxWvec),nw,Ne).';
    %     % initial state and STM
    %     x_prev=p(idx.x{1}); S_n=S_last;
    %     % initialise Dw for indirect gradients
    %     Dw=zeros(nx,nw*Ne);
    %     row0=0;
    %     for j=1:Ne
    %         uj=obj.Uhist(:,j); wj=wHist(j);
    %         x_k=x_prev;
    %         % propagate through Nint sub-steps
    %         for m=1:Nint
    %             [x_prev,S_n]=obj.model.step(x_prev,uj,wj,S_n,true);
    %             % chain-rule propagation of Dw
    %             Dw = S_n(:,1:nx)*Dw;
    %             col_idx=(j-1)*nw+(1:nw);
    %             Dw(:,col_idx) = Dw(:,col_idx) + S_n(:,nx+1:nx+nw);
    %         end
    %         x_next=x_prev;
    %         rows=row0+(1:nx);
    %         ceq(rows)=p(idx.x{j+1})-x_next;
    %         % analytical Jacobian blocks
    %         gradceq(idx.x{j},rows)   = -eye(nx);
    %         gradceq(idx.x{j+1},rows) =  eye(nx);
    %         % direct+indirect disturbance gradient
    %         % direct+indirect disturbance gradient
    %         idxW_j = [idx.w{1:j}];
    %         gradceq(idxW_j, rows) = -Dw(:,1:(j*nw))';
    %         % debug: compare direct vs indirect gains
    % 
    %         % gradceq(idx.w{1:j},rows) = -Dw(:,1:j*nw);
    %         % % debug probes
    %         if obj.debug
    %             fprintf('interval %d: |Sw|=%.3e, |Dw(:,end)|=%.3e\n', j, norm(S_n(:,nx+1)), norm(Dw(:,j*nw)));
    %             % evaluate FD for first component of p(idx.w{j})
    %             % (optional) insert here a small finite-diff check
    %         end
    %            row0=row0+nx;
  %         end

   %% V3
    % function [c,ceq,gradc,gradceq] = nonlFun(p)
    %     if obj.method=="single"
    %         c=[]; ceq=[]; gradc=[]; gradceq=[];   % single-shoot → no constraints
    %         return
    %     end
    % 
    %     % ----- convenience --------------------------------------------------
    %     nVar   = numel(p);
    %     ceq    = zeros(nx*Ne,1);
    %     gradceq= zeros(nVar , nx*Ne);          % (var × constraint-row)
    %     idxW   = [idx.w{:}];                   % flat indices of all w–slots
    %     Nint   = round(obj.Ts/obj.dt);
    % 
    %     % ----- reshape decision vector --------------------------------------
    %     wHist  = reshape(p(idxW),nw,Ne).';     % (Ne × nw)  oldest → newest
    %     x_prev = p(idx.x{1});                  % x^{-Ne}
    % 
    %     % ----- sensitivity accumulator  Dw  ---------------------------------
    %     Dw = zeros(nx,nw*Ne);                  % ∂x_prev / ∂w^{-Ne:-1}
    %     S_n = S_last;                          % STM from previous solve
    %     Sw     = S_n(:, obj.nx+1);   % gust col
    % %         Su     = S_n(:, obj.nx+2 : obj.nx+1+obj.nu);
    %     row0 = 0;
    %     for j = 1:Ne
    %         uj = obj.Uhist(:,j);               % known control
    %         wj = wHist(j);                     % scalar (nw = 1 here)
    % 
    %         % -- integrate through the whole sub-interval --------------------
    %         for m = 1:Nint
    %             [x_prev,S_n] = obj.model.step(x_prev,uj,wj,S_n,true);
    %             % chain rule  Dw ← Φ_x*Dw  +  Φ_w
    %             % Dw = S_n(:,1:nx) * Dw;                    % Φ_x * accumulated
    %             col_j = (j-1)*nw + (1:nw);               % current w column(s)
    %             % Dw(:,col_j) = Dw(:,col_j) + S_n(:,nx+1:nx+nw);   % add Φ_w
    %             Dw = S_n(:,1:nx) * Dw + S_n(:,nx+1:nx+nw);
    % 
    %         end
    % 
    %         % -------------- build constraint row ----------------------------
    %         rows = row0 + (1:nx);
    %         ceq(rows) = p(idx.x{j+1}) - x_prev;          % x^{j+1} − Φ(...)
    % 
    %         % Jacobian wrt states -------------------------------------------
    %         gradceq(idx.x{j}  ,rows) = -eye(nx);
    %         gradceq(idx.x{j+1},rows) =  eye(nx);
    % 
    %         % Jacobian wrt ALL disturbances seen so far ---------------------
    %         idxW_block           = [idx.w{1:j}];         % 1 … j   (vector)
    %         % gradceq(idxW_block , rows) = -Dw(:,1:j*nw)'; % note the transpose
    %         % gradceq(idx.w{j} , rows) = -Sw.'; % note the transpose
    %         gradceq(idx.w{j} , rows) = -Sw.'; % note the transpose
    % 
    %         % -------------- on-line FD probe (first state) ------------------
    %         if obj.debug && j==1                         % probe only once
    %             % epsFD = 1e-6;
    %             epsFD = 1e-6;
    %             pFD   = p;  
    %             pFD(idx.w{1}) = pFD(idx.w{1}) + epsFD;
    %             % integrate one FD step
    %             xFD   = p(idx.x{1});  Sfd = S_last;
    %             for m = 1:Nint
    %                 [xFD,Sfd] = obj.model.step(xFD,uj,wj+epsFD,Sfd,false);
    %             end
    %             dNum = (xFD(1)-x_prev(1))/epsFD;         % FD slope (1st state)
    %             dAn  = -Dw(1,1);                         % analytic in same row
    %             fprintf('[FD-check] row-1 col-1  analytic = %+9.6f\n\r ', dAn);
    %             fprintf('FD = %+9.6f  rel.err = %.2e\n\r', ...
    %                     dNum,abs(dAn-dNum)/max(1,abs(dNum)));
    %         end
    % 
    %         row0 = row0 + nx;            % advance to next continuity block
    % 
    % 
    %     end   
    %     % gradceq=gradceq'; 
    %     gradc = [];
    %     c=[];         % m×n
    %     S_last = S_n;           % capture STM for outer accessor
    % end
   %% V-4
    % % function [c,ceq,gradc,gradceq] = nonlFun(p)
    % %     if obj.method=="single"
    % %         c=[]; ceq=[]; gradc=[]; gradceq=[];  return
    % %     end
    % %     gradc=[];
    % %     % sizes ----------------------------------------------------------------
    % %     nVar   = numel(p);
    % %     ceq    = zeros(nx*Ne,1);
    % %     gradceq= zeros(nVar , nx*Ne);          % (var × row)
    % %     idxW   = [idx.w{:}];                   % flat disturbance indices
    % %     Nint   = round(obj.Ts / obj.dt);       % steps inside one shooting node
    % % 
    % %     % unpack decision vector ------------------------------------------------
    % %     wHist  = reshape(p(idxW),nw,Ne).';     % Ne × nw
    % %     x_prev = p(idx.x{1});                  % x^{-Ne}
    % %     S_n    = S_last;                       % STM at start of horizon
    % %     % x_n = x_prev;
    % %     Dw = zeros(nx,nw*Ne);                  % ∂x_prev / ∂w^{-Ne:−1}
    % %     % Sw     = S_n(:, obj.nx+1);   % gust col
    % %     Sw     = S_last(:, obj.nx+1);   % gust col
    % %     row0 = 0;
    % %     % for j = 1:Ne
    % %     %     % uj = obj.Uhist(:,j);       wj = wHist(j);
    % %     %     % % --- integrate whole sub-interval ---------------------------------
    % %     %     % % for m = 1:Nint  
    % %     %     % %     [x_prev,S_n] = obj.model.step(x_prev, uj, wj, S_n, true);
    % %     %     % % 
    % %     %     % %     % ----- Dw update (chain-rule) ---------------------------------
    % %     %     % %     Dw = S_n(:,1:nx) * Dw + S_n(:,nx+1:nx+nw);   % Φ_x*Dw  + Φ_w
    % %     %     % % 
    % %     %     % %     if obj.debug && j==1 && m==1
    % %     %     % %         fprintf('[dbg] first Δt: col_j = %d (should be 1)\r\n', ...
    % %     %     % %                  (j-1)*nw+1);
    % %     %     % %     end
    % %     %     % % end
    % %     %     % Dw_interval = zeros(nx,nw);           % clear for this interval
    % %     %     % for m = 1:Nint
    % %     %     %     [x_prev,S_n] = obj.model.step(x_prev, uj, wj, S_n, true);
    % %     %     %     Dw_interval = Dw_interval + S_n(:,nx+1:nx+nw);   % accumulate Φ_w
    % %     %     % end
    % %     %     % Dw = S_n(:,1:nx)*Dw + Dw_interval;    % ONE Φ_x per interval
    % %     %     uj = obj.Uhist(:,j);
    % %     %     wj = wHist(j);
    % %     % 
    % %     %     % ---- integrate whole sub-interval & collect DIRECT Φ_w -------------
    % %     %     Sw_sum = zeros(nx,nw);          % Σ Φ_x^{k} Φ_w  over Δt steps
    % %     %     for m = 1:Nint
    % %     %         [x_prev,S_n] = obj.model.step(x_prev, uj, wj, S_n, true);
    % %     %         Sw_sum       = S_n(:,nx+1:nx+nw) + S_n(:,1:nx)*Sw_sum;
    % %     %     end
    % %     %     if j==1
    % %     %         fprintf('[check] after j=1  ‖x_prev‖ = %.3e\n',norm(x_prev));
    % %     %         % numerical derivative over this interval only
    % %     %         eps = 1e-6;
    % %     %         xFD1 = x_prev;             % state already integrated
    % %     %         % re-integrate with +eps gust for interval-1 only
    % %     %         [xTmp, ~] = obj.model.step( p(idx.x{1}), uj, wj+eps, S_last, false);
    % %     %         for m2 = 2:Nint
    % %     %             [xTmp, ~] = obj.model.step( xTmp, uj, wj   , [],      false);
    % %     %         end
    % %     %         dNum1 = (xTmp(1) - x_prev(1))/eps;     % first state component
    % %     %         fprintf('[check] Sw_sum(1) = %.3e   FD interval-1 = %.3e\n',...
    % %     %                 Sw_sum(1), dNum1);
    % %     %     end
    % %     %     % Sw_sum now = Σ Φ_x^{k} Φ_w  (k=0…Nint-1)  for this interval
    % %     % 
    % %     %     % ---- propagate OLD sensitivities one level -------------------------
    % %     %     Dw = S_n(:,1:nx) * Dw;          % Φ_x^(interval) * previous Dw
    % %     % 
    % %     %     % ---- write DIRECT column for current w_j ---------------------------
    % %     %     col_j = (j-1)*nw + (1:nw);      % column(s) for w^{-Ne+j-1}
    % %     %     Dw(:,col_j) = Dw(:,col_j) + Sw_sum;
    % %     %     % ---------------- continuity equation -----------------------------
    % %     %     rows                 = row0 + (1:nx);
    % %     %     ceq(rows)            = p(idx.x{j+1}) - x_prev;
    % %     %     % ceq(rows)            = p(idx.x{j+1}) - x_n;
    % %     %     gradceq(idx.x{j},rows)   = -eye(nx);
    % %     %     gradceq(idx.x{j+1},rows) =  eye(nx);
    % %     % 
    % %     %     % ---------- Jacobian wrt ALL disturbances up to j -----------------
    % %     %     idxW_block               = [idx.w{1:j}];      % 1 … j
    % %     %     idxWa               = [idx.w{:}];      % 1 … j
    % %     %     % gradceq(idxW_block,rows) = -Dw(:,1:j*nw)';    % (len(idxW_block) × nx)
    % %     %     gradceq(idx.w{j},rows) = -Sw';    % (len(idxW_block) × nx)
    % %     %     % gradceq(idxWa((j-1)*obj.nw+(1:obj.nw)), rows) = -S_n(:,obj.nx+obj.nw); % nw is only 1
    % %     % 
    % %     %     % ---------- light debug dump --------------------------------------
    % %     %     if obj.debug && j==1
    % %     %         fprintf('Sw(1:5)  = %s\r\n', mat2str(S_n(1:5,nx+1)',5));
    % %     %         fprintf('Dw(1:5)  = %s\r\n', mat2str(Dw(1:5,1)',5));
    % %     %     end
    % %     for j = 1:Ne
    % %         uj = obj.Uhist(:,j);  wj = wHist(j);
    % %         Sw_sum = zeros(nx,nw);
    % % 
    % %         for m = 1:Nint
    % %             [x_prev,S_n] = obj.model.step(x_prev, uj, wj, S_n, true);
    % %             Sw_sum       = S_n(:,nx+1:nx+nw) + S_n(:,1:nx)*Sw_sum;
    % %         end
    % %         Sw     = S_n(:, obj.nx+1);   % gust col
    % % 
    % %         % propagate older columns once
    % %         Dw = S_n(:,1:nx) * Dw;
    % %         % store newest column
    % %         col_j = (j-1)*nw + (1:nw);
    % %         Dw(:,col_j) = Sw_sum;
    % % 
    % %         % continuity row j -----------------------------------------------
    % %         rows = row0 + (1:nx);
    % %         ceq(rows) = p(idx.x{j+1}) - x_prev;
    % %         gradceq(idx.x{j},rows)   = -eye(nx);
    % %         gradceq(idx.x{j+1},rows) =  eye(nx);
    % %         gradceq(idx.w{j},rows)   = -Sw;     % only NEW column
    % % 
    % %         % quick sanity check
    % %         if obj.debug && j==1
    % %             epsFD = 1e-6;
    % %             x_fd  = p(idx.x{1});
    % %             for m2=1:Nint, [x_fd,~]=obj.model.step(x_fd,uj,wj+epsFD,[],false); end
    % %             fprintf('[sens] Sw_sum(1)=%.3e  FD=%.3e\n', Sw_sum(1), (x_fd(1)-x_prev(1))/epsFD);
    % %         end
    % % 
    % %         row0 = row0 + nx;         % advance row pointer
    % %     end
    % % 
    % %     % ---------- exact FD probe for row-13 / col-13 ------------------------
    % %     if obj.debug
    % %         rowProbe = nx + 1;                       % 1st state in block-2
    % %         colProbe = (Ne+1)*nx + 1;                % first disturbance
    % %         dAn = gradceq(colProbe,rowProbe);        % analytic entry
    % % 
    % %         epsFD = 1e-6;
    % %         pFD   = p;  pFD(colProbe) = pFD(colProbe) + epsFD;
    % %         ceqFD = ceqOnly(pFD, S_last);                    % helper below
    % %         dNum  = (ceqFD(rowProbe)-ceq(rowProbe))/epsFD;
    % % 
    % %         fprintf('[FD-exact] row %d col %d  an = %+9.6e  FD = %+9.6e  rel = %.2e\r\n', ...
    % %                  rowProbe,colProbe,dAn,dNum,abs(dAn-dNum)/max(1,abs(dNum)));
    % %     end
    % % 
    % %     % gradceq = gradceq.';   % (m × n) orientation
    % %     c       = [];
    % %     S_last  = S_n;         % capture STM for next outer call
    % % 
    % %     nz = find(abs(gradceq(:))>1e-12);
    % %     fprintf('[map] first 20 non-zeros (row,col):\n');
    % %     for kd = 1:min(20,numel(nz))
    % %         [r,c] = ind2sub(size(gradceq), nz(kd));
    % %         fprintf('  (%d,%d) = %.2e\n', r,c, gradceq(r,c));
    % %     end
    % % end
    % function [c,ceq,gradc,gradceq] = nonlFun(p)
    % 
    %     if obj.method=="single"
    %         c=[]; ceq=[]; gradc=[]; gradceq=[]; return
    %     end
    %     nVar   = numel(p);
    %     % wHist  = reshape(p(idxW),nw,Ne).';     % Ne × nw
    % 
    %     gradc   = [];                         % no inequality constraints
    %     ceq     = zeros(nx*Ne,1);
    %     gradceq = zeros(nVar , nx*Ne);        % nVar-by-constraints  ✅
    % 
    % %     nVar   = numel(p);
    % %     ceq    = zeros(nx*Ne,1);
    % %     gradceq= zeros(nVar , nx*Ne);          % (var × row)
    %     idxW   = [idx.w{:}];                   % flat disturbance indices
    %     Nint   = round(obj.Ts / obj.dt);       % steps inside one shooting node
    % % 
    % %     % unpack decision vector ------------------------------------------------
    %     wHist  = reshape(p(idxW),nw,Ne).';     % Ne × nw
    % %     x_prev = p(idx.x{1});                  % x^{-Ne}
    % %     S_n    = S_last;   
    % 
    % 
    %     Dw      = zeros(nx, nw*Ne);           % all disturbance columns
    %     x_prev  = p(idx.x{1});
    %     S_n     = S_last;                     % STM at start of horizon
    %     row0    = 0;
    %     x_n = x_prev;
    %     % for j = 1:Ne
    %     for j = Ne:-1:1
    %         uj  = obj.Uhist(:,j);
    %         wj  = wHist(j);
    %         Sw_sum = zeros(nx,nw);
    %         Sw     = S_n(:, obj.nx+1);   % gust col
    %         Sx     = S_n(:,1: obj.nx);   % state col
    %         for m = 1:Nint
    %             [x_next,S_n] = obj.model.step(x_prev, uj, wj, S_n, true);
    %             x_prev = x_next;
    %             Sw_sum       = S_n(:,nx+1:nx+nw) + S_n(:,1:nx)*Sw_sum;
    %         end
    %         % Sw     = S_n(:, obj.nx+1);   % gust col
    %         % Sx     = S_n(:,1: obj.nx);   % state col
    %         % debug
    %         if j == 1
    %             % ----- 2a: verify Φx, Φw dimensions & rough magnitude ------------
    %             fprintf('[Φx] diag(1)=%.4f  should NOT be 1.0\n', Sx(1,1));
    %             fprintf('[Φw] (1)=%.4e     must match FD below\n', Sw(1));
    % 
    %             % ----- 2b: forward finite-difference on w ------------------------
    %             epsFD = 1e-4;
    %             % xFD   = x_prev;
    %             xFD   = x_n;
    %             for m = 1:Nint
    %                 [xFD,~] = obj.model.step(xFD , uj , wj+epsFD , [] , false);
    %             end
    %             dFD = (xFD-x_next)/epsFD;         % nx×1
    %             fprintf('[FD-check] Φw(1)=%.4e  FD(1)=%.4e  ratio=%.2f\n',...
    %                     -Sw(1), -dFD(1), (-Sw(1))/(-dFD(1)));
    %         end
    %         % x_n = x_prev;
    %         % propagate old columns, add new one
    %         Dw = S_n(:,1:nx)*Dw;
    %         col_j = (j-1)*nw + (1:nw);
    %         Dw(:,col_j) = Sw_sum;
    % 
    %         rows        = row0 + (1:nx);
    %         % ce = ceqOnly( p, S_last);
    %         % ceq = ce+ceq;
    %         ceq(rows)   = p(idx.x{j+1}) - x_prev;
    %         % ceq(rows)   = p(idx.x{j+ 1}) - x_n;
    %         % ceq(rows)   = x_n-p(idx.x{j+1}) ;
    %         x_n = x_prev;
    % 
    %         % ceq(rows)   = p(idx.x{j+1}) + x_prev;
    %         % ceq(rows)   = x_prev-p(idx.x{j+1});
    % 
    %         % gradceq(idx.x{j},rows)   = -eye(nx)+gradceq(idx.x{j},rows);
    %         % gradceq(idx.x{j},rows)   = -eye(nx);
    %         % gradceq(idx.x{j+1},rows) =  eye(nx)+gradceq(idx.x{j+1},rows);           
    %         gradceq(idx.x{j+1},rows) =  eye(nx);           
    %         % gradceq(idx.x{j+1},rows)   = -Sx;
    %         gradceq(idx.x{j},rows)   = -Sx;
    %         % gradceq(idx.x{j},rows)   = -Sx+gradceq(idx.x{j},rows);
    %         % gradceq(idx.x{j+1},rows)   = Sx;
    %         % gradceq(idx.x{j+1},rows) =  Sx+gradceq(idx.x{j+1},rows);
    % 
    %         % gradceq(idx.w{j},rows)   = -Sw_sum;    % newest column only
    %         gradceq(idx.w{j},rows)   = -Sw;    % newest column only
    %         % gradceq(idx.w{j},rows)   = Sw;    % newest column only
    % 
    %         fprintf('[map]  row %3d  uses x-block %2d  (col %d)\n',...
    %         rows(1), j, idx.x{j}(1));
    % 
    % 
    %         % --- pick the suspicious row/column reported by CheckGradients ----------
    %         r = 29;            % the row index in ceq   (e.g. 50)
    %         c = idx.x{1}(9);  % the matching column in p (x{1}, 50-th state)
    % 
    %         % --- analytic derivative you have just written into gradceq -------------
    %         dAn = gradceq(c, r);        % note: nVar × nCeq orientation
    % 
    %         % --- finite difference of that single element ----------------------------
    %         epsFD = 1e-7;
    %         pFD   = p;           pFD(c) = pFD(c) + epsFD;
    % 
    %         % recompute *only* ceq(r) at pFD
    %         % 1. unpack the variables needed for interval j=1
    %         x_j   = pFD(idx.x{1});          % perturbed x₀
    %         w_j   = pFD(idx.w{1});          % perturbed w₀
    %         u_j   = obj.Uhist(:,1);
    % 
    %         x_tmp = x_j;
    %         S_tmp = [];                     % no STM needed here
    %         for m = 1:Nint                  % integrate Ts
    %             [x_tmp,~] = obj.model.step(x_tmp, u_j, w_j, S_tmp, false);
    %         end
    %         ceqPert = pFD(idx.x{2}(r)) - x_tmp(r);   % *same residual formula*
    %         ceqBase =  p (idx.x{2}(r)) - x_next(r);
    % 
    %         dFD = (ceqPert - ceqBase)/epsFD;
    % 
    %         fprintf('analytic %.6f   FD %.6f   rel.err %.2e\n', dAn, dFD, ...
    %                 abs(dAn-dFD)/max(1,abs(dFD)));
    % 
    %         row0 = row0 + nx;
    %     end
    %     obj.Sprev = S_n;          
    %     % S_last = S_n;          
    %     c = [];                     % no inequalities
    %     % -------- NO transpose here --------
    % end
%% V5
    % function [c, ceq, gradc, gradceq] = nonlFun( z)
    % % NONLINEARCONSTRAINTS Compute continuity constraints and Jacobians for MHE (multiple shooting).
    % % Inputs:
    % %   obj : nMHEv2 object containing system dynamics and horizon configuration.
    % %   z   : Decision vector [x0, x1, ..., x_N, w0, ..., w_{N-1}] (concatenated states and disturbances).
    % % Outputs:
    % %   c      : Inequality constraints (empty for this MHE formulation).
    % %   ceq    : Equality constraint residuals (continuity errors between shooting nodes).
    % %   gradc  : Gradient of inequality constraints (empty).
    % %   gradceq: Gradient of equality constraints (Jacobian of continuity constraints).
    % 
    %     % Unpack problem dimensions
    %     Nx   = obj.nx;             % State dimension
    %     Nw   = obj.nw;             % Disturbance input dimension per interval (0 if none)
    %     Nseg = obj.Ne;           % Number of shooting segments (horizon intervals)
    %     Nint   = round(obj.Ts / obj.dt);
    %     % Split decision vector into state nodes X and disturbance inputs W
    %     % X is size (Nx x (Nseg+1)), where columns 1..(Nseg+1) are x0...x_N.
    %     X = reshape(z(1 : Nx*(Nseg+1)), Nx, Nseg+1);
    %     if Nw > 0
    %         % W is size (Nw x Nseg), columns 1..Nseg are w0...w_{N-1}
    %         W = reshape(z(Nx*(Nseg+1) + 1 : Nx*(Nseg+1) + Nw*Nseg), Nw, Nseg);
    %     else
    %         W = zeros(0, Nseg);
    %     end
    % 
    %     % Initialize outputs
    %     c    = [];  % no inequality constraints
    %     ceq  = zeros(Nx * Nseg, 1);
    %     % Preallocate Jacobian (as dense matrix). We fill only the nonzero structure.
    %     gradc   = [];
    %     gradceq = zeros(length(z), Nx * Nseg);
    % 
    %     % Time step and integration settings
    %     dt    = obj.dt;           % total duration of each segment
    %     Nsub  = Nint;         % number of integration sub-steps per segment (for accuracy)
    %     h     = dt / Nsub;        % sub-step size
    % 
    %     % Loop over each segment in the horizon
    %     for i = 1:Nseg
    %     % for i = Nseg:-1:1
    %         % Indices for this segment's constraint block in the output vector
    %         idxC_start = (i-1)*Nx + 1;             % starting index in ceq for segment i
    %         idxC_end   = i * Nx;                  % ending index in ceq for segment i
    % 
    %         % Extract states and disturbance for this segment
    %         % Initial state x_i and target state x_{i+1}:
    %         x0 = X(:, i);                         % state at beginning of segment i
    %         x_next = X(:, i+1);                   % state at end of segment (the optimization variable to match integration)
    %         % Disturbance (or unknown input) over this segment:
    %         if Nw > 0
    %             w = W(:, i);
    %         else
    %             w = [];                           % no disturbance estimated
    %         end
    %         % Known control input over this segment (e.g., measured actuator commands):
    %         if isprop(obj, 'U') && ~isempty(obj.U)
    %             u = obj.U(:, i);
    %         else
    %             u = [];                           % no control or control not provided
    %         end
    % 
    %         % Initialize state and sensitivity for integration
    %         x = x0; 
    %         % S_n = obj.Sprev;                          % STM: identity at start of segment (∂x/∂x0)
    %         if Nw > 0
    %             R = zeros(Nx, Nw);                % Input sensitivity (∂x/∂w); zero at start (no initial dependence on w)
    %         end
    % 
    %         % Integrate dynamics over the segment [t_i, t_{i+1}] with RK4, propagating state, S, and R
    %         uj = obj.Uhist(:,i);       % known control
    %         % wj = wHist(i);             % decision-var disturbance
    %         S_n = [eye(numel( obj.xhat)) , zeros(numel(obj.xhat),1+obj.nu)];         % initial STM (nx×(nx+nw))
    % 
    %         % for m = 1:Nint
    %             % [x,S_n] = obj.model.step(x, uj, w, S_n, true);
    %         % end
    % 
    %         % sensitivities ------------------------------------------------
    %         Sx = S_n(:,1:nx);          % Φ_x
    %         Sw = S_n(:,nx+1);          % Φ_w   (nw == 1)
    % 
    %         % After integration, 'x' is the model-predicted state at t_{i+1}
    %         x_end = x;
    %         % Compute continuity residual: difference between optimization variable and integrated state
    %         ceq_seg = x_next - x_end; 
    %         ceq(idxC_start:idxC_end) = ceq_seg;
    % 
    %         % Compute Jacobian blocks for this segment's constraints:
    %         % Indices of variable blocks in the decision vector
    %         x_i_idx    = (i-1)*Nx + 1 : i*Nx;             % indices of x_i in z
    %         x_next_idx = (i)*Nx + 1     : (i+1)*Nx;       % indices of x_{i+1} in z
    %         if Nw > 0
    %             % Starting index of the w block in z = Nx*(Nseg+1)+1
    %             w_block_offset = Nx * (Nseg+1);
    %             w_i_idx = w_block_offset + (i-1)*Nw + 1 : w_block_offset + i*Nw;  % indices of w_i
    %         end
    % 
    %         % Fill Jacobian for continuity constraint of segment i:
    %         % ∂ceq/∂x_i = - (∂x_end/∂x_i) = -S_x (state transition matrix)
    %         % gradceq(x_i_idx,    idxC_start:idxC_end) = - Sx.';   % use S' so that gradceq(:,j) gives ∂(ceq_j)/∂x correctly
    %         gradceq(x_i_idx,    idxC_start:idxC_end) = - Sx.';   % use S' so that gradceq(:,j) gives ∂(ceq_j)/∂x correctly
    %         % ∂ceq/∂x_{i+1} = I (current segment's next state appears with +1)
    %         gradceq(x_next_idx, idxC_start:idxC_end) =  eye(Nx);
    %         % ∂ceq/∂w_i = - (∂x_end/∂w) = -R (if disturbances are estimated)
    %         if Nw > 0
    %             % gradceq(w_i_idx, idxC_start:idxC_end) = - R';
    %             gradceq(w_i_idx, idxC_start:idxC_end) = - Sw;
    %         end
    % 
    %         % DEBUG: verbose output for residuals and Jacobian magnitudes
    %         if isprop(obj, 'debug') && obj.debug
    %             fprintf('Segment %d continuity residual norm = %.3e (max comp = %.3e)\n', ...
    %                     i, norm(ceq_seg), max(abs(ceq_seg)));
    %             % Check if any large residual remains (should be near zero at solution)
    %             % Print sample of sensitivity matrices
    %             fprintf('  STM (∂x_end/∂x_i) diag range: [%.3f, %.3f]\n', min(diag(S_n)), max(diag(S_n)));
    %             if Nw > 0
    %                 fprintf('  Input sens (∂x_end/∂w_i) range: [%.3e, %.3e]\n', min(R(:)), max(R(:)));
    %             end
    %             % Optional finite-difference check for first state component:
    %             % if obj.debug > 1
    %             %     % Perturb initial state slightly to numerically verify gradient
    %             %     epsVal = 1e-6;
    %             %     e = zeros(Nx,1); e(1) = epsVal;
    %             %     % Re-integrate with perturbed x_i:
    %             %     x_pert = X(:, i) + e;
    %             %     x_test = x_pert; S_test = eye(Nx);
    %             %     if Nw > 0, R_test = zeros(Nx, Nw); end
    %             %     for jj = 1:Nsub
    %             %         [xd1, A1, Bw1] = obj.modelFcn(x_test, u, w);
    %             %         k1St = A1 * S_test;
    %             %         if Nw > 0, k1Rt = A1 * R_test + Bw1; end
    %             %         x_midt = x_test + 0.5*h*xd1;
    %             %         S_midt = S_test + 0.5*h*k1St;
    %             %         if Nw > 0, R_midt = R_test + 0.5*h*k1Rt; end
    %             %         [xd2, A2, Bw2] = obj.modelFcn(x_midt, u, w);
    %             %         k2St = A2 * S_midt;
    %             %         if Nw > 0, k2Rt = A2 * R_midt + Bw2; end
    %             %         x_midt2 = x_test + 0.5*h*xd2;
    %             %         S_midt2 = S_test + 0.5*h*k2St;
    %             %         if Nw > 0, R_midt2 = R_test + 0.5*h*k2Rt; end
    %             %         [xd3, A3, Bw3] = obj.modelFcn(x_midt2, u, w);
    %             %         k3St = A3 * S_midt2;
    %             %         if Nw > 0, k3Rt = A3 * R_midt2 + Bw3; end
    %             %         x_endt = x_test + h*xd3;
    %             %         S_endt = S_test + h*k3St;
    %             %         if Nw > 0, R_endt = R_test + h*k3Rt; end
    %             %         [xd4, A4, Bw4] = obj.modelFcn(x_endt, u, w);
    %             %         k4St = A4 * S_endt;
    %             %         if Nw > 0, k4Rt = A4 * R_endt + Bw4; end
    %             %         x_test = x_test + (h/6)*(xd1 + 2*xd2 + 2*xd3 + xd4);
    %             %         S_test = S_test + (h/6)*(k1St + 2*k2St + 2*k3St + k4St);
    %             %         if Nw > 0
    %             %             R_test = R_test + (h/6)*(k1Rt + 2*k2Rt + 2*k3Rt + k4Rt);
    %             %         end
    %             %     end
    %             %     % Compute perturbed final state and residual:
    %             %     x_end_pert = x_test;
    %             %     res_pert = x_next - x_end_pert;
    %             %     % Finite-difference approximation of first column of gradient:
    %             %     grad_fd_col1 = (res_pert - ceq_seg) / epsVal;   % numeric ∂ceq/∂(x_i(1))
    %             %     fprintf('  Finite-diff vs Analytic for ∂ceq/∂x_i(1): [');
    %             %     fprintf(' %.3e', grad_fd_col1);
    %             %     fprintf(' ] vs [');
    %             %     fprintf(' %.3e', gradceq(x_i_idx, idxC_start));
    %             %     fprintf(' ]\n');
    %             % end
    %         end
    %     end  % end for each segment
    % 
    %     % (gradc remains empty as there are no inequality constraints)
    % end
%% V6
        % function [c, ceq, gradc, gradceq] = nonlFun(z)
        %         % Dimensions
        %         nx = obj.nx;               % state dimension
        %         nw = obj.nw;               % disturbance input dimension
        %         N  = obj.Ne;                % number of shooting intervals
        %         gradc=[];
        %         % Initialize outputs
        %         ceq    = zeros(nx * N, 1);
        %         gradceq = zeros(length(z), nx * N);
        %         % Loop over each shooting interval
        %         S0 = [eye(numel( obj.xhat)) , zeros(numel(obj.xhat),1+obj.nu)];         % initial STM (nx×(nx+nw))
        %         for k = 0:(N-1)
        %         % for k = (N-1):-1:0
        %             % Indices for variables in z
        %             idx_xk   = (k * nx + 1) : ((k+1) * nx);         % x_k indices
        %             idx_xkp1 = ((k+1) * nx + 1) : ((k+2) * nx);     % x_{k+1} indices
        %             idx_wk   = (N+1) * nx + (k * nw + 1 : (k+1) * nw);  % w_k indices (assuming states x_0...x_N come first)
        %             % Extract current values
        %             xk   = z(idx_xk);
        %             xkp1 = z(idx_xkp1);
        %             wk   = z(idx_wk);
        %             uk   = obj.Uhist(:, k+1);  % known control input for interval k
        %             % Propagate dynamics for interval k (with sensitivity)
        %             S0 = [eye(numel( obj.xhat)) , zeros(numel(obj.xhat),1+obj.nu)];         % initial STM (nx×(nx+nw))
        %             % for m = 1:Nint
        %             % 
        %             %     [x_end, S_end] = obj.model.step(xk, uk, wk, S0,true);
        %             %     if m ~=Nint
        %             %         S0 = S_end;
        %             %         xk = x_end;
        %             %     end
        %             % end
        %             % x_end = xk;
        %             % S_end = S0;
        %             [x_end, S_end] = obj.model.step(xk, uk, wk, S0, true);
        %             % Compute residual for this interval: ceq_k = x_{k+1} - x_end
        %             ceq_k = xkp1 - x_end;
        %             ceq((k*nx+1):(k+1)*nx) = ceq_k;
        %             % Fill gradient blocks:
        %             % ∂ceq/∂x_{k+1}: identity (nx×nx)
        %             gradceq(idx_xkp1, (k*nx+1):(k+1)*nx) = eye(nx);
        %             % ∂ceq/∂x_k: - (∂x_end/∂x_k)  => use -S_end(:,1:nx)^T
        %             Sx = S_end(:, 1:nx);               % size nx×nx
        %             gradceq(idx_xk, (k*nx+1):(k+1)*nx) = - Sx.';   % transpose to nx×nx, then assign
        %             % gradceq(idx_xk, (k*nx+1):(k+1)*nx) = - Sx';   % transpose to nx×nx, then assign
        %             % ∂ceq/∂w_k: - (∂x_end/∂w_k)  => use -S_end(:,nx+1:nx+nw)^T
        %             Sw = S_end(:, nx+1:nx+nw);         % size nx×nw
        %             gradceq(idx_wk, (k*nx+1):(k+1)*nx) = - Sw.';   % transpose to nw×nx
        %             % gradceq(idx_wk, (k*nx+1):(k+1)*nx) = - Sw';   % transpose to nw×nx
        %             % gradceq(idx_wk, (k*nx+1):(k+1)*nx) = - Sw;   % transpose to nw×nx
        %         end
        %         c = [];   % no inequality constraints
        %     end
        
%% V7
                % ---------- Continuity constraints & exact Jacobian --------------
        function [c,ceq,gradc,gradceq] = nonlFun(z)

            if obj.method == "single"
                c = [];  ceq = [];  gradc = [];  gradceq = [];
                return                                               % ← early exit
            end

            % ---- 1. convenient sizes ------------------------------------
            nx = obj.nx;   nw = obj.nw;   N = obj.Ne;
            nVar = numel(z);

            ceq     = zeros(nx*N,1);
            gradceq = zeros(nVar , nx*N);        % (row = variable , col = constraint)

            % ---- 2. reshape decision vector -----------------------------
            X = reshape(z(1:(N+1)*nx)         , nx , N+1);          % x₀ … x_N
            W = reshape(z((N+1)*nx+1:end)     , nw , N      );      % w₀ … w_{N-1}

            % ---- 3. loop over shooting intervals ------------------------
            row0   = 0;                      % running row pointer in ceq/grad
            % S_last = obj.Sprev;              % STM that arrived from previous solve
            S = S_last;
            Nint   = round(obj.Ts/obj.dt);   % # internal integrator steps
            S = [eye(nx) , zeros(nx,nw+obj.nu)];     %  <-- **reset here**
            for k = 0:N-1
            % for k = 0:N-2
                % indices ----------------------------------------------------------
                idx_xk   = k*nx      + (1:nx);
                idx_xkp1 = (k+1)*nx  + (1:nx);
                idx_wk   = (N+1)*nx  + k*nw + (1:nw);
            
                % variables --------------------------------------------------------
                xk   = X(:,k+1);
                wk   = W(:,k+1);
                % xk   = obj.Xhist(:,k+2);
                % wk   = obj.Whist(:,k+2);
                uk   = obj.Uhist(:,k+1);
            
                % -------- START OF INTERVAL: fresh sensitivities ------------------
                S = [eye(nx) , zeros(nx,nw+obj.nu)];     %  <-- **reset here**
                x = xk;
            
                % integrate Ts -----------------------------------------------------
                % for m = 1:Nint+1
                for m = 1:Nint
                    [x,S] = obj.model.step(x,uk,wk,S,true);
                end
            
                x_end = x;  S_end = S;                    % @ t_{k+1}
            
                % residual & Jacobian ---------------------------------------------
                rows = row0 + (1:nx);
                ceq(rows)        = X(:,k+2) - x_end;
                % ceq(rows)        = obj.Xhist(:,k+3) - x_end;
            
                Sx = S_end(:,1:nx);           % size nx×nx
                Sw = S_end(:,nx+1:nx+nw);     % size nx×nw
            
                % gradceq(idx_xkp1,rows) =  Sx.';
                gradceq(idx_xkp1,rows) =  eye(nx);
                gradceq(idx_xk  ,rows) = -Sx.';          % (nx × nx)
                % gradceq(idx_xk  ,rows) = -eye(nx);          % (nx × nx)
                % gradceq(idx_xk  ,rows) = -Sx;          % (nx × nx)
                gradceq(idx_wk  ,rows) = -Sw.';          % (nw × nx)
                % gradceq(idx_wk  ,rows) = -Sw;          % (nw × nx)
            
                % if obj.debug %&& k==0
                %     fprintf('[dbg]  k=%d  ‖res‖=%.2e  max|Sx|=%.2e  max|Sw|=%.2e\n', ...
                %              k,norm(ceq(rows)),max(abs(Sx(:))),max(abs(Sw(:))));
                % end
                S_last = S_end;
                row0 = row0 + nx;             % next block
            end


            % ---- 4. housekeeping ----------------------------------------
            c        = [];    gradc = [];
            % obj.Sprev  = S_last;              % warm-start for *next* call
        end


        function S_out = getSTM(), S_out = S_last; end
    % ========== helper: ceq without Jacobian (for FD probe) ================
    function ce = ceqOnly( pLocal, S_last)
        idx = obj.buildIndexMaps();
        ce = zeros(obj.nx*obj.Ne,1);
        x_tmp = pLocal(idx.x{1});
        S_tmp = S_last;
        out   = 0;
        for jj = 1:obj.Ne
            uj = obj.Uhist(:,jj);
            wj = pLocal(idx.w{jj});
            for mm = 1:Nint
                [x_tmp,S_tmp] = obj.model.step(x_tmp, uj, wj, S_tmp, false);
            end
            ce(out+(1:obj.nx)) = pLocal(idx.x{jj+1}) - x_tmp;
            out = out + obj.nx;
        end
    end
    % ------------------------------------------------ BOUNDS ---------
    xLrep=repmat(obj.xL,Ne+1,1); xUrep=repmat(obj.xU,Ne+1,1);
    wLrep=repmat(obj.wL,Ne,1);  wUrep=repmat(obj.wU,Ne,1);

    nlp.cost=@costFun; 
    nlp.nonl=@nonlFun;
    nlp.getSTM = @getSTM;
    nlp.lb=[xLrep;wLrep]; nlp.ub=[xUrep;wUrep];


end
%--------------------------------------------------------------------
function debugPlots(obj,t_k,info)
    wH = info.wHorizon;
     % ---- initialise persistent store on first call --------------------
    if ~isfield(obj.dbg,'t')
        obj.dbg.t = [];          % Ns×1   sample times
        obj.dbg.W = [];          % Ns×Ne  horizon matrix (rows = samples)
    end
    if isempty(obj.dbg)
        obj.dbg.t=t_k; 
        obj.dbg.W=wH(:)';
        obj.dbg.cont=norm(info.wHorizon);
    else
        obj.dbg.t(end+1)=t_k; obj.dbg.W(end+1,1:obj.Ne)=wH(:)';
        obj.dbg.cont(end+1)=info.continuity;

    end
    figure(2001),clf,hold on
    for s=1:obj.Ne, stairs(obj.dbg.t,obj.dbg.W(:,s)); end
    hold off,grid on,xlabel('t [s]'),ylabel('w'),title('horizon slots'), legend

    figure(2004),clf, hold on,plot(obj.dbg.t,obj.dbg.cont(2:end)); hold off;grid on
    xlabel('t [s]'); ylabel('|continuity|'); title('MHE feasibility'), legend
end
%--------------------------------------------------------------------

end
end
