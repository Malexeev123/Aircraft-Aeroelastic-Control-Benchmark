classdef nMHEv2 < AeroFlex.ctrl.EstimatorBase
%========================================================
%  NON-LINEAR  MOVING–HORIZON ESTIMATOR  (multiple-shoot)
%========================================================
%  * buffer-free – stores its own X,U,Y,W histories
%  * Diehl-style warm start (shift & predict)
%  * Exact Jacobian of continuity constraints
%========================================================

    properties
        % ---- scalar config
        dt  double
        Ts  double
        Ne  double
        method (1,1) string {mustBeMember(method,["multiple","single"])} = "multiple"
        % ---- dimensions
        nx  double
        nu  double
        ny  double
        nw  double
        % ---- weights / bounds
        Qe  double;  Pe double;  Re double
        xL  double;  xU double
        wL  double;  wU double
        % ---- model & buffers
        model    % ROMIntegrator
        sensor
        Sprev
        Xhist
        Uhist
        Yhist
        Whist
        k   double                         % sample counter
        % ---- optimiser state
        z0  double
        H   double
        solverOpts
        firstFullSolve
        % ---- outputs
        % xhat double
        % what double
        % ---- debug
        debug logical = true
        dbg   struct
        xlast

        solverName string = "fmincon"
        sqpSolver
        sqpCheckDone logical = false
    end
%======================================================================
methods
%----------------------------------------------------------------------
function obj = nMHEv2(cfg,beam,aero,base, trim)
% Constructor
%----------------------------------------------------------------------
    obj@AeroFlex.ctrl.EstimatorBase(cfg, trim);

    % -------- 1. basic numbers ---------------------------------------
    obj.dt   = cfg.sim.dt;
    obj.Ts   = cfg.ctrl.Ts;
    obj.Ne   = cfg.ctrl.Ne;
    obj.method = lower(string(cfg.ctrl.mhe_method));

    obj.model  = cfg.modelHandle(cfg,beam,aero,base);
    obj.sensor = cfg.sensorHandle(beam,cfg);

    obj.nx = size(obj.model.L,1);
    obj.nu = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    obj.ny = cfg.ny;
    obj.nw = cfg.nw;

    % -------- 2. weights / bounds ------------------------------------
    obj.Qe = cfg.Qe;     obj.Pe = cfg.Pe;     obj.Re = cfg.Re;
    xTrim  = trim.states;                 %  n × 1
    % margin = 1.60;                            % 30 % either side
    margin = 10;                            % 30 % either side
    % cfg.xL = (1 - margin) .* xTrim;
    % cfg.xU = (1 + margin) .* xTrim;

    % absFloor = 1e-1;  
    absFloor = 1;  

    delta   = max(abs(margin .* xTrim) , absFloor);   % n×1 non‑negative
    cfg.xL      = xTrim - delta;
    cfg.xU      = xTrim + delta;
    obj.xL = cfg.xL;     obj.xU = cfg.xU;
    obj.wL = cfg.wL;     obj.wU = cfg.wU;

    % -------- 3. optimiser options -----------------------------------
    obj.solverOpts = optimoptions('fmincon',...
        'Algorithm','sqp',...     % interior-point | sqp
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true,...
        'HessianApproximation','lbfgs',...
        'Display','none',...
        'OptimalityTolerance',1e-4, ...
        'StepTolerance',1e-7, ...
        'ConstraintTolerance',1e-6);

    % obj.solverOpts.FiniteDifferenceType ='central';
    % obj.solverOpts.CheckGradients = true;
            % 'FiniteDifferenceStepSize',1e-8,...

    %======================================================================
    % nMHE solver selection
    %======================================================================
    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'mheSolver')
        obj.solverName = lower(string(cfg.ctrl.mheSolver));
    else
        obj.solverName = "fmincon";
    end
    
    switch obj.solverName
    
        case "fmincon"
            obj.sqpSolver = [];
    
            % Make sure gradient checking is not accidentally active online.
            if isprop(obj.solverOpts,'CheckGradients')
                obj.solverOpts.CheckGradients = false;
            end
    
        case "custom_sqp"
            sqpOpts = AeroFlex.optim.SQPSolver.defaultOptions();
    
            if isfield(cfg.ctrl,'mheSqp') && ~isempty(cfg.ctrl.mheSqp)
                f = fieldnames(cfg.ctrl.mheSqp);
                for ii = 1:numel(f)
                    sqpOpts.(f{ii}) = cfg.ctrl.mheSqp.(f{ii});
                end
            end
    
            obj.sqpSolver = AeroFlex.optim.SQPSolver(sqpOpts);
    
        otherwise
            error('nMHEv2:Solver', ...
                  'Unknown cfg.ctrl.mheSolver = "%s". Use "fmincon" or "custom_sqp".', ...
                  obj.solverName);
    end
    
    obj.sqpCheckDone = false;

    % -------- 4. rolling buffers  ------------------------------------
    obj.Xhist = zeros(obj.nx,obj.Ne+1);
    obj.Uhist = zeros(obj.nu,obj.Ne);
    obj.Yhist = zeros(obj.ny,obj.Ne+1);
    obj.Whist = zeros(obj.nw,obj.Ne);
    obj.k     = 0;

    % -------- 5. initial guess ---------------------------------------
    xTrim           = trim.states;        % equilibrium state
    obj.xhat        = xTrim;
    obj.Xhist(:,:)  = repmat(xTrim,1,obj.Ne+1);
    obj.Whist(:,:)  = 0;                      % zero-gust history
    obj.what = zeros(obj.Ne*obj.nw,1);

    idx   = obj.buildIndexMaps();
    obj.z0           = zeros((obj.Ne+1)*obj.nx + obj.Ne*obj.nw,1);
    for j = 1:obj.Ne+1
        obj.z0(idx.x{j}) = xTrim;
    end
    % obj.z0(idx.x{:}) = xTrim;                 % arrival state
    % obj.z0(idx.x{end}) = xTrim;                 % arrival state
    obj.H            = speye(numel(obj.z0));  % LBFGS initial Hessian

    % initial STM (∂x/∂x0  and  ∂x/∂w) – identity / zeros
    obj.Sprev = [eye(obj.nx) , zeros(obj.nx,obj.nw+obj.nu)];
    obj.dbg = struct('t', [], 'W', [], 'cont', 0);
    obj.firstFullSolve = true;
end
%----------------------------------------------------------------------
function [xhat,what,info] = estimate(obj,z_k,u_prev,t_k)
% Real-time entry point (called every Ts)
%----------------------------------------------------------------------
    if ~isempty(obj.xlast) | obj.k < obj.Ne
        obj.xlast = obj.xhat;
    % else
        % obj.xlast = obj.Xhist(:,1);

    end
    % ---- 1. roll/append histories every call ----------------------------
    obj.Xhist(:,1:end-1) = obj.Xhist(:,2:end);
    obj.Yhist(:,1:end-1) = obj.Yhist(:,2:end);
    obj.Whist(:,1:end-1) = obj.Whist(:,2:end);
    obj.Uhist(:,1:end-1) = obj.Uhist(:,2:end);
    
    obj.Yhist(:,end) = z_k;
    obj.Uhist(:,end) = u_prev;
    
    % obj.k = min(obj.k + 1, obj.Ne);
    obj.k = min(obj.k + 1, obj.Ne+1);
    
    % ---- 2. until horizon is sufficiently populated ---------------------
    % if obj.k < obj.Ne
    if obj.k < obj.Ne+1
        obj.Xhist(:,end) = obj.xhat;
    
        xhat = obj.xhat;
        what = obj.what;
        idx = obj.buildIndexMaps();

        obj.z0(idx.x{obj.k}) = xhat;
        info = struct( ...
            'cost',0, ...
            'continuity',0, ...
            'exitflag',0, ...
            'wHorizon',zeros(obj.Ne,obj.nw));
    
        return
    end

    idx = obj.buildIndexMaps();
    % 
    % I think this stays commented out
    % So do not overwrite an internal shooting node before the solve.
    % The shifted solution is already stored in obj.z0 from the previous call.
    % obj.z0(idx.x{obj.k}) = obj.xhat;

    % % current arrival state -------------------------
    % obj.z0(idx.x{1}) = obj.xhat;
    % Nint=round(obj.Ts/obj.dt);
    % 
    % % forward-predict the horizon with zero gust ----
    % x_tmp = obj.xhat;  S_tmp = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
    %     for j = 1:obj.Ne
    %         uj = obj.Uhist(:,j);  wj = obj.z0(idx.w{j});
    %         for m = 1:Nint
    %             [x_tmp,S_tmp] = obj.model.step(x_tmp, uj, wj, S_tmp, true);
    %         end
    %         obj.z0(idx.x{j+1}) = x_tmp;
    %     end
    % obj.Sprev = S_tmp;      

    % initial disturbance guess ---------------------
    % obj.z0([idx.w{:}]) = 0;            % => keep the same dimension
    if obj.firstFullSolve
        obj.z0([idx.w{:}]) = 0;
        obj.firstFullSolve = false;
    else
        % Keep shifted disturbance guess.
        % Optionally:
        obj.z0(idx.w{end}) = obj.z0(idx.w{end-1});  % hold-last
        % or
        % obj.z0(idx.w{end}) = zeros(obj.nw,1);     % only new slot zero
    end
    % % ---- 4. assemble NLP & solve ------------------------------------
    % nlp = obj.assembleWindow();
    % [p,fval,exitflag] = fmincon(nlp.cost,obj.z0,...
    %                             [],[],[],[],nlp.lb,nlp.ub,...
    %                             nlp.nonl,obj.solverOpts);
    %======================================================================
    % Assemble and solve nMHE multiple-shooting NLP
    %======================================================================
    nlp = obj.assembleWindow();
    
    switch obj.solverName
    
        case "fmincon"
    
            [p,fval,exitflag,output,lambda] = fmincon( ...
                nlp.cost,obj.z0, ...
                [],[],[],[], ...
                nlp.lb,nlp.ub, ...
                nlp.nonl,obj.solverOpts);
    
        case "custom_sqp"
    
            if ~obj.sqpCheckDone && obj.sqpSolver.options.CheckGradientsOnce
    
                % Global directional check.
                obj.sqpSolver.checkGradients(nlp.cost,nlp.nonl,obj.z0,nlp.lb,nlp.ub);
    
                % MHE-specific block check: separates X and W sensitivity errors.
                obj.localCheckMHEEqualityGradientBlocks(nlp,obj.z0,nlp.lb,nlp.ub);
    
                obj.sqpCheckDone = true;
            end
    
            [p,fval,exitflag,output,lambda] = obj.sqpSolver.solve( ...
                nlp.cost,obj.z0,nlp.lb,nlp.ub,nlp.nonl);
    
        otherwise
            error('nMHEv2:Solver','Unhandled solverName = "%s".', obj.solverName);
    end

    % ---- 5. post-process solution -----------------------------------
    % obj.z0  = obj.shiftGuess(p);          % warm start for next call
    % % obj.xhat = p(idx.x{end});
    % obj.xhat = p(idx.x{end});
    % obj.what = p([idx.w{:}]);
    % % obj.xlast = obj.xhat;
    % obj.xlast = p(idx.x{1});


    % obj.Xhist(:,end) = obj.xhat;
    % obj.Whist(:,end) = obj.what(end-obj.nw+1:end);
    Xsol = reshape(p([idx.x{:}]), obj.nx, obj.Ne+1);
    Wsol = reshape(p([idx.w{:}]), obj.nw, obj.Ne);
    
    

    obj.Xhist = Xsol;
    obj.Whist = Wsol;
    
    obj.xhat = Xsol(:,end);
    obj.what = Wsol(:);
    obj.xlast = p(idx.x{1});

    obj.z0 = obj.shiftGuess(p);


    obj.Sprev        = nlp.getSTM();

    xhat = obj.xhat;   
    what = obj.what;

    info.cost       = fval;
    info.exitflag   = exitflag;
    info.output     = output;
    info.lambda     = lambda;
    info.wHorizon   = reshape(obj.what,obj.nw,[])';
    
    [~,ceq]         = nlp.nonl(p);
    info.continuity = norm(ceq);
    
    if isfield(output,'constrviolation')
        info.constrviolation = output.constrviolation;
    end
    
    if isfield(output,'firstorderopt')
        info.firstorderopt = output.firstorderopt;
    end
    
    if isfield(output,'stepsize')
        info.stepsize = output.stepsize;
    end
    
    if isfield(output,'qpExitflag')
        info.qpExitflag = output.qpExitflag;
    end
    if obj.debug  
        obj.debugPlots(t_k,info);  
        % obj.debugPlots2(t_k,info);  
    
    end
    % 
    % wH = info.wHorizon(:);
    % 
    % fprintf('[MHE w check] t = %.4f | min(w)=%.4e | max(w)=%.4e | wL=%.4e | wU=%.4e | hitL=%d | hitU=%d\n', ...
    % t_k, min(wH), max(wH), obj.wL(1), obj.wU(1), ...
    % any(wH <= obj.wL(1) + 1e-6), ...
    % any(wH >= obj.wU(1) - 1e-6));
end
%----------------------------------------------------------------------
end  % public methods

    %======================================================================
    methods (Access = private)
    %----------------------------------------------------------------------
    function idx = buildIndexMaps(obj)
        nx = obj.nx;   nw = obj.nw;   N = obj.Ne;
        idx.x = cell(N+1,1);
        for j = 1:N+1,   idx.x{j} = (j-1)*nx + (1:nx);  end
        startW = (N+1)*nx;
        idx.w  = cell(N,1);
        for j = 1:N,     idx.w{j} = startW + (j-1)*nw + (1:nw);  end
    end
    %----------------------------------------------------------------------
    % function zShift = shiftGuess(obj,p)
    %     idx = obj.buildIndexMaps();  N = obj.Ne;
    %     zShift = zeros(size(p));
    %     for j = 1:N
    %         zShift(idx.x{j}) = p(idx.x{j+1});   % x_{k+1} becomes new x_k
    %         zShift(idx.w{j}) = p(idx.w{j+1});   % w_{k+1} → w_k
    %     end
    %     zShift(idx.x{N+1}) = p(idx.x{N+1});     % keep last state
    % end

    function zShift = shiftGuess(obj,p)
        idx = obj.buildIndexMaps(); nx=obj.nx; nw=obj.nw; Ne=obj.Ne;
        zShift=zeros(size(p));

        % Shift states.
        for j = 1:Ne
            zShift(idx.x{j}) = p(idx.x{j+1});
        end
    
        % Shift disturbances.
        for j = 1:Ne-1
            zShift(idx.w{j}) = p(idx.w{j+1});
        end
        zShift(idx.w{Ne}) = p(idx.w{Ne});   % hold-last
    
        % Predict terminal state using last shifted state, last known control,
        % and held disturbance. 
        % Trying this out to see if it could get a better warm-start
        x = zShift(idx.x{Ne});
        w = zShift(idx.w{Ne});
        u = obj.Uhist(:,end);
    
        Nsub = round(obj.Ts/obj.dt);
        Sdummy = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
    
        for m = 1:Nsub
            [x,Sdummy] = obj.model.step(x,u,w,Sdummy,false);
        end
    
        zShift(idx.x{Ne+1}) = x;
    end
    %----------------------------------------------------------------------
    function nlp = assembleWindow(obj)
    % Build cost, constraints, bounds for fmincon
    %----------------------------------------------------------------------
        nx = obj.nx;  nw = obj.nw;  N = obj.Ne;
        idx = obj.buildIndexMaps();
        Nsub = round(obj.Ts/obj.dt);                 % integrator steps/slot
        S_last = obj.Sprev;                          % STM at arrival

        % -------------------------- COST ---------------------------------
        function [J,gradJ] = costFun(p)
            J = 0;  
            gradJ = zeros(numel(p),1);

            % arrival penalty --------------------------------------------
            arrivalX   = obj.Xhist(:,  end-obj.Ne);   % <— oldest kept sample
            % arrivalW   = obj.Whist(:,  end-obj.Ne+1);   % <— oldest kept sample
            arrivalW   = obj.Whist(:,  end-obj.Ne+1);   % <— oldest kept sample
            % dx   = p(idx.x{1}) - obj.Xhist(:,1);
            dx   = p(idx.x{1}) - arrivalX;
            J    = J + 0.5*dx'*obj.Pe*dx;
            gradJ(idx.x{1}) = obj.Pe*dx;

            % first-slot disturbance penalty
            % dw   = p(idx.w{1}) - obj.Whist(:,1);
            % dw   = p(idx.w{1}) - arrivalW;
            % J    = J + 0.5*dw'*obj.Re*dw;
            % gradJ(idx.w{1}) = obj.Re*dw;

            % innovation terms ------------------------------------------
            C     = obj.sensor.PhiY;
            Cf    = zeros(obj.ny,obj.nx);
            Cf(:,obj.model.idx.q1) = C;

            Xmat  = reshape(p([idx.x{:}]),nx,[]);
            for j = 1:N+1
                e = obj.Yhist(:,j) - C*Xmat(obj.model.idx.q1,j);
                J = J + 0.5*e'*obj.Qe*e;
                gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cf'*(obj.Qe*e);
            end
            
            if isfield(obj.cfg,'mhe') && isfield(obj.cfg.mhe,'penalizeAllW') && obj.cfg.mhe.penalizeAllW
                for j = 1:N
                    wj = p(idx.w{j});
                    J = J + 0.5*wj.'*obj.Re*wj;
                    gradJ(idx.w{j}) = gradJ(idx.w{j}) + obj.Re*wj;
                end
            else
                dw = p(idx.w{1}) - arrivalW;
                J = J + 0.5*dw.'*obj.Re*dw;
                gradJ(idx.w{1}) = gradJ(idx.w{1}) + obj.Re*dw;
            end

            % Smoothness Term
            if isfield(obj.cfg,'RdWe')
                RdW = obj.cfg.RdWe;
            else
                RdW = [];
            end
            
            if ~isempty(RdW)
                Wmat = reshape(p([idx.w{:}]),obj.nw,[]);
            
                for j = 2:N
                    dwSmooth = Wmat(:,j) - Wmat(:,j-1);
            
                    J = J + 0.5*dwSmooth.'*RdW*dwSmooth;
            
                    gradJ(idx.w{j})   = gradJ(idx.w{j})   + RdW*dwSmooth;
                    gradJ(idx.w{j-1}) = gradJ(idx.w{j-1}) - RdW*dwSmooth;
                end
            end

            % --------------------------------------------------------------
            % Terminal disturbance penalty:
            % discourages nonzero estimated gust at the newest/current slot.
            % --------------------------------------------------------------
            if isfield(obj.cfg,'mhe') && isfield(obj.cfg.mhe,'RwTerminal') && ...
                    ~isempty(obj.cfg.mhe.RwTerminal)
            
                RwT = obj.cfg.mhe.RwTerminal;
            
                wN = p(idx.w{N});
            
                J = J + 0.5*wN.'*RwT*wN;
                gradJ(idx.w{N}) = gradJ(idx.w{N}) + RwT*wN;
            end
        end

        % ---------------------- CONTINUITY ------------------------------
        function [c,ceq,gradc,gradceq] = nonlFun(z)
            % if obj.method == "single"
            %     c=[];ceq=[];gradc=[];gradceq=[];return
            % end


            ceq     = zeros(nx*(N),1);
            % gradceq = zeros(numel(z),nx*N);
            % gradceq = zeros(nx*(N), numel(z));
            S = [eye(nx) , zeros(nx,nw+obj.nu)];
            Sx0 = S(:,1:nx);
            Sw0 = S(:,nx+1:nx+nw);

            nnzPerInterval = nx + nnz(Sx0) + nnz(Sw0);
            gradceq = spalloc(nx*N, numel(z), N*nnzPerInterval);

            X = reshape(z(1:(N+1)*nx),nx,N+1);
            W = reshape(z((N+1)*nx+1:end),nw,N);
            % S = [eye(nx) , zeros(nx,nw+obj.nu)];
            S_last = [eye(nx) , zeros(nx,nw+obj.nu)];
            % S = obj.Sprev;           % <-- STM that ended the *previous* RT call  = I at k=0 only
            % ceq(1:nx,:) = obj.xhat - X(:, 1);
            % gradceq(1:nx, 1:nx) = -S(1:nx, 1:nx);
            % ceqArrival       = obj.xlast - X(:,1) ;         % (nx×1)
            % ceqArrival       = X(:,1) - obj.xlast;         % (nx×1)
            % ceq(1:nx) = ceqArrival;
            % gradceq(1:nx, 1:nx) = eye(nx,nx);

            Ulocal = obj.Uhist(:,1:obj.Ne);      % exactly Ne columns
            % Ulocal = obj.Uhist(:,1:obj.Ne+1);      % exactly Ne columns
            % subplot(211), stairs(obj.Uhist(1,:)), title('Uhist after shift')
            % subplot(212), stairs(W), title('W decision order')
            row0 = 0;
            % row0 = row0 + nx;

            Nend = N-1;
            % Nend = N-2;
            % Nend = N-3;
            S = [eye(nx) , zeros(nx,nw+obj.nu)];

            for k = 0:Nend
                idx_xk   = k*nx     + (1:nx);
                idx_xkp1 = (k+1)*nx + (1:nx);
                idx_wk   = (N+1)*nx + k*nw + (1:nw);
                idx_wkp1   = (N+1)*nx + (k+1)*nw + (1:nw);
                % disp(k+1)
                x = X(:,k+1);

                % x_n = x;
                wk = W(:,k+1);
                % uk = obj.Uhist(:,k+1);
                uk = Ulocal(:,k+1);
                % S = S_last;
                S = [eye(nx) , zeros(nx,nw+obj.nu)];
                for m = 1:Nsub
                    [x,S] = obj.model.step(x,uk,wk,S,true);
                end
                x_end = x;   S_end = S;

                Sx = S_end(:,1:nx);
                Sw = S_end(:,nx+1:nx+nw);

                rows = row0 + (1:nx);
                % disp(rows(1));
                % ceq(rows) = X(:,k+1) - x;
                ceq(rows) = X(:,k+2) - x_end;
                % ceq(rows) = x_end - X(:,k+2) ;
                % ceq(rows) = X(:,k+1) - x_end;
                % ceq(rows) = X(:,k+2) - S*[x_n ; wk; uk];
                % ceq(rows) = X(:,k+2) - S_end*[x_end; wk; uk]-x_end;
                % ceq(rows) = X(:,k+2) - S_end*[x_n; wk; uk]-x_end;


                % gradceq(rows, idx_xkp1) =  eye(nx);     % d/dx_{k+1}
                gradceq(rows, idx_xkp1) =  speye(nx);     % d/dx_{k+1}

                % gradceq(rows, idx_xk) = -Sx; % d/dx_k
                % gradceq(rows, idx_wk) = -Sw;
                gradceq(rows, idx_xk) = -sparse(Sx); % d/dx_k
                gradceq(rows, idx_wk) = -sparse(Sw);
               
                row0 = row0 + nx;
                S_last = S_end;                 
                % S_last = S;                    
            end
            % ceq = [ceqArrival ; ceq];
            % gArrive = [eye(nx), zeros(nx, numel(z)-nx)];
            % gradceq = [gArrive ; gradceq];
            c=[]; gradc=[];
            gradceq = gradceq.';


            % % --- inequality:  |Δw_k| ≤ dWmax ------------------
            % dW  = diff([obj.Whist(:,1) , W],1,2);           % w₀-w_{-1} … w_N-w_{N-1}
            % c   = abs(dW)' - obj.wU;                         % vectorised (N×1)
            % gradc = zeros(nw*N , length(z));               % build once
            % 
            % for k = 1:N
            %     signk = sign(dW(k));            % ±1
            %     gradc(k, idx.w{k}) =  signk;    % ∂/∂w_k
            % end
        end

        % return STM for next warm-start
        function S_out = getSTM(),  S_out = S_last;  end

        % -------- bounds -------------------------------------------------
        xLrep = repmat(obj.xL,N+1,1);
        xUrep = repmat(obj.xU,N+1,1);
        wLrep = repmat(obj.wL,N,1);
        wUrep = repmat(obj.wU,N,1);

        nlp.cost   = @costFun;
        nlp.nonl   = @nonlFun;
        nlp.getSTM = @getSTM;
        nlp.lb     = [xLrep; wLrep];
        nlp.ub     = [xUrep; wUrep];
    end
    %----------------------------------------------------------------------
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
        hold off,grid on,xlabel('t [s]'),ylabel('w'),title('horizon slots'), legend; %drawnow
        
        figure(2004),clf, hold on,plot(obj.dbg.t,obj.dbg.cont(2:end)); hold off;grid on
        xlabel('t [s]'); ylabel('|continuity|'); title('MHE feasibility'), legend
    end

    function debugPlots2(obj,t_k,info)
        %DEBUGPLOTS Time-aligned MHE disturbance horizon plot.
        
            wH = info.wHorizon;   % expected Ne x nw or Ne x 1 for nw=1
        
            if size(wH,2) > size(wH,1) && obj.nw == 1
                wH = wH(:);
            end
        
            if ~isfield(obj.dbg,'tSolve')
                obj.dbg.tSolve = [];
                obj.dbg.tSlot  = [];
                obj.dbg.Wslot  = [];
                obj.dbg.cont   = [];
            end
        
            obj.dbg.tSolve(end+1,1) = t_k;
            obj.dbg.cont(end+1,1)   = info.continuity;
        
            % Slot physical times: oldest to newest.
            tSlots = zeros(obj.Ne,1);
            for s = 1:obj.Ne
                tSlots(s) = t_k - (obj.Ne - s)*obj.Ts;
            end
        
            obj.dbg.tSlot = [obj.dbg.tSlot; tSlots(:)];
            obj.dbg.Wslot = [obj.dbg.Wslot; wH(:)];
        
            figure(3001); clf;
            plot(obj.dbg.tSlot,obj.dbg.Wslot,'.-','LineWidth',1.2);
            grid on;
            xlabel('physical time associated with MHE slot [s]');
            ylabel('\hat{w}');
            title('Time-aligned MHE disturbance estimates');
        
            figure(3004); clf;
            semilogy(obj.dbg.tSolve,max(obj.dbg.cont,eps),'LineWidth',1.2);
            grid on;
            xlabel('solve time [s]');
            ylabel('||c_{eq}||');
            title('MHE feasibility');
    end

    function localCheckMHEEqualityGradientBlocks(obj,nlp,z,lb,ub)
    %LOCALCHECKMHEEQUALITYGRADIENTBLOCKS Check MHE dynamic equality Jacobian
    % by perturbing only X variables and only W variables separately.
    
        idx = obj.buildIndexMaps();
    
        nx = obj.nx;
        nw = obj.nw;
        Ne = obj.Ne;
    
        z  = z(:);
        lb = lb(:);
        ub = ub(:);
        n  = numel(z);
    
        z = min(max(z,lb),ub);
    
        [~,ceq0,~,Gceq] = nlp.nonl(z);
    
        if size(Gceq,1) ~= n
            Gceq = Gceq.';
        end
    
        % If ceq includes X0-arrival equality, dynamic rows start after nx.
        if numel(ceq0) >= nx*(Ne+1)
            rowsDyn = nx+1:numel(ceq0);
        else
            rowsDyn = 1:numel(ceq0);
        end
    
        h = 1e-6;
    
        % X-only direction.
        dX = zeros(n,1);
        rng(301);
        for k = 1:Ne+1
            dX(idx.x{k}) = randn(nx,1);
        end
        dX = dX / max(norm(dX),eps);
    
        % W-only direction.
        dW = zeros(n,1);
        rng(302);
        for k = 1:Ne
            dW(idx.w{k}) = randn(nw,1);
        end
        dW = dW / max(norm(dW),eps);
    
        errX = obj.localDirectionalEqualityError(nlp,z,dX,h,lb,ub,Gceq,rowsDyn);
        errW = obj.localDirectionalEqualityError(nlp,z,dW,h,lb,ub,Gceq,rowsDyn);
    
        fprintf('\n');
        fprintf('======================================================================\n');
        fprintf(' NMHE DYNAMIC-EQUALITY BLOCK GRADIENT CHECK\n');
        fprintf('======================================================================\n');
        fprintf('  dynamic defect rel error, X-only direction : %.3e\n', errX);
        fprintf('  dynamic defect rel error, W-only direction : %.3e\n', errW);
        fprintf('======================================================================\n\n');
    end

    function err = localDirectionalEqualityError(obj,nlp,z,d,h,lb,ub,Gceq,rowsDyn)
    %LOCALDIRECTIONALEQUALITYERROR Directional FD check for selected equality rows.
    
    
        z  = z(:);
        d  = d(:);
        lb = lb(:);
        ub = ub(:);
    
        active = abs(d) > 0;
        hBound = inf;
    
        idxU = active & isfinite(ub);
        if any(idxU)
            hBound = min(hBound,min((ub(idxU)-z(idxU))./abs(d(idxU))));
        end
    
        idxL = active & isfinite(lb);
        if any(idxL)
            hBound = min(hBound,min((z(idxL)-lb(idxL))./abs(d(idxL))));
        end
    
        if isfinite(hBound)
            h = min(h,0.49*hBound);
        end
    
        if h <= eps
            err = nan;
            return
        end
    
        [~,ceqp] = nlp.nonl(z + h*d);
        [~,ceqm] = nlp.nonl(z - h*d);
    
        ceqFD = (ceqp(:)-ceqm(:))/(2*h);
        ceqAN = Gceq.'*d;
    
        rowsDyn = rowsDyn(:);
    
        den = max([1,norm(ceqFD(rowsDyn),inf),norm(ceqAN(rowsDyn),inf)]);
    
        err = norm(ceqFD(rowsDyn)-ceqAN(rowsDyn),inf)/den;
    end
    %----------------------------------------------------------------------
    end   % private methods
end   % classdef