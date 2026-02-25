classdef SimRunner
    properties
        cfg, beam, aero, base
        L                       % linear operator
        Lfac, Ufac, piv                 % LU factors of (I-γΔt L)   (implicit solve)
        sens                    % sensitivity matrix S (optional)
        controller              % (optional)
        estimator
        sensor

        idx
        parConst

        log, buff
        % time
    end

    % properties (Constant, Hidden)
    %     gamma = (2 - sqrt(2))/2;
    %     delta = -2*sqrt(2)/3;
    % end
    % 
    methods
        function obj = SimRunner(cfg, beam, aero, base, trim, controller)

            if  exist('trim', 'var')
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
                % absMargin = [1e-3;   1e-2; …];   % vector of per‑state minima
                % cfg.xL = xTrim - max( abs(margin*xTrim) , absMargin );
                % cfg.xU = xTrim + max( abs(margin*xTrim) , absMargin );
            end

            if nargin<6 && strcmp(cfg.case, 'openLoop')
                controller = []; estimator  =[];
            else
                % estimator  = cfg.estimatorHandle(cfg);
                % controller = cfg.controllerHandle(cfg);            
            end
            obj.cfg  = cfg;  obj.beam = beam;  obj.aero = aero;  obj.base = base;
            % obj.controller = controller; obj.estimator = estimator;
            [obj.L,~] = AeroFlex.core.assemble_L_matrix(cfg,beam,aero,base);

            % --- IMEX DIRK(2,3,2) coefficients (Kennedy–Carpenter 2003) ----

            gamma = (2 - sqrt(2))/2;
            % delta = -2*sqrt(2)/3;

            I   = speye(size(obj.L));
            A  = I - gamma*cfg.sim.dt*obj.L;
            [obj.Lfac, obj.Ufac, obj.piv] = lu(A,'vector');
            
            obj.idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
            
            % forces0 = beam.eta_e + aero.DataMatrix.forces_aero_beam_dof;
            forces0 = beam.eta_e + aero.DataMatrix.forces_aero_beam_dof;
            % forces0 = beam.eta_e;
             %── build constant parameter sub‑struct -----------------------
            obj.parConst = struct( ...
                'Gamma1',   beam.Gamma1, ...
                'Gamma2',   beam.Gamma2, ...
                'Gamma_g',  base.Gamma_g, ...
                'Gamma_xi', base.Gamma_xi, ...
                'forces_0', forces0, ...   % steady aero load Nm×1
                'scaleAero', cfg.flight.a*cfg.flight.t_inf, ...
                't_inf', cfg.flight.t_inf, ...
                'scaleA', cfg.flight.a, ...
                'gustSet', cfg.gust, ...
                'Bw', aero.forceMap.Bw, ...
                'Dw', aero.forceMap.Dw, ...
                'Bdel',  aero.forceMap.B_delta, ...
                'Ddel',  aero.forceMap.D_delta, ...
                'Bddel', aero.forceMap.B_ddelta, ...
                'Dddel', aero.forceMap.D_ddelta, ...
                'gust_input', aero.gust_input, ...
                'dt', cfg.sim.dt, ...
                'Na', aero.Na);

            ctrlSurfDefs = size(obj.parConst.Bdel,2) + size(obj.parConst.Bddel,2);
            obj.sensor = AeroFlex.sensor.WingVelSensor(beam, cfg);
            obj.cfg.ny = size(obj.sensor.PhiY,1);
            obj.cfg.nu = ctrlSurfDefs;
            obj.cfg.nw = size(obj.parConst.Bw,2);

            if  exist('trim', 'var')
                obj.log.x(:,1) = trim.states;
                obj.log.y(:,1) = trim.sensor.y(:, end);
                trim.w_r = zeros(size(aero.gust_input));
                obj.cfg.ref = trim;


            else
                obj.log.x(:,1) = zeros(size(I,1), 1);
                obj.log.y(:,1) = zeros(size(obj.sensor.PhiY,1), 1);
                obj.cfg.ref = [];

            end
            

        end
        % ------------------------------------------------------------------
        function [t,X, log, S] = run(obj,x0,tEnd,storeSens)
            arguments
                obj
                x0           %     (:,1)  double
                tEnd        %                 double
                storeSens   %(1,1) logical = false
            end

            if nargin<4, storeSens=false; end
            gamma = (2 - sqrt(2))/2;
            delta = -2*sqrt(2)/3;
            % delta =1 - 1/(2*gamma);
            

            % ---------- aliases -----------------------------------------------------
            cfg    = obj.cfg;             dt = obj.cfg.sim.dt;  
            Nx     = size(obj.L,1);       %Nm   = cfg.control.Nc;
            if ~strcmp(obj.cfg.case, 'openLoop')
                Ts     = cfg.ctrl.Ts;  
                Nint = round(Ts/dt);

            end

            t  = 0:dt:tEnd;  
            Nt = numel(t);
            X  = zeros(numel(x0), Nt);  
            X(:,1)=x0;
            % ---------- allocate ----------------------------------------------------
            % Nt     = floor((tEnd-0)/dt) + 1;
            % t      = zeros(1,Nt);         t(1) = 0;
            % X      = zeros(Nx,Nt);        X(:,1) = x0;
            % S      = eye(Nx, Nx+cfg.nu);                % STM  t0
            obj.log.U = zeros(cfg.nu, Nt);
            obj.log.Y = zeros(cfg.ny, Nt);              % measurement buffer
            obj.log.t = t;              obj.log.Scond = zeros(1,Nt);
            
                        

            parBase = obj.parConst;          % copy‑by‑value once
            parBase.u_ctrl = zeros(4,1);     % default open‑loop
            parBase.gust   = 0;
            gust_input = parBase.gust_input;
            uPrev = parBase.u_ctrl;
            if storeSens  % Need to implement
                % S = eye(numel(x0)); obj.sens = zeros(size(S,1),size(S,2),Nt); obj.sens(:,:,1)=S;
                S = [eye(numel(x0)), zeros(numel(x0), 1+size(parBase.u_ctrl,1))]; 
                obj.sens = zeros(size(S,1),size(S,2),Nt); obj.sens(:,:,1)=S;
                cfg.buff.S(:,:,1) = S;

            else
                S = [];
                cfg.buff.S(:,:,1) = S;

            end

            % ---------- shared buffer ----------------------------------------------
            cfg.buff = struct( ...
               'X'     , X , ...
               'S'     , obj.sens , ...
               'U'     , obj.log.U , ...
               'Y'     , obj.log.Y , ...
               't'     , obj.log.t , ...
               'sensorMat', obj.sensor.PhiY , ...
               'valid', 0 );
            cfg.buff.idx = obj.idx;         % struct with q1,q2,qxi,…

            if ~strcmp(cfg.case, 'openLoop')
                estimator  = cfg.estimatorHandle(cfg);
                controller = cfg.controllerHandle(cfg); 
                obj.estimator  = estimator;
                obj.controller = controller;
            end

            kSam = 0;   uPrev = zeros(cfg.nu,1);

            % Main Time Loop
            for k=1:Nt-1

                x_n = X(:,k);
                

                time_k  = t(k);
                disp('time'); disp(time_k);

                yk = obj.sensor.measure(x_n(obj.idx.q1));           % v_y at 5 positions

                % if ~strcmp(obj.cfg.case, 'openLoop')
                %     % ----- synthetic measurement -----
                %     % yk = obj.sensor.measure(x_n(obj.idx.q1));           % v_y at 5 positions
                %     % obj.log.y(:,k) = yk;
                % 
                %     % ----- estimator -------------------------------------------------
                %     xhat = obj.estimator.estimate(obj.cfg, yk,uPrev,time_k);
                %     obj.log.est(:,k) = xhat;
                %     % ----- controller ------------------------------------------------
                %     uk   = obj.controller.computeControl(xhat,time_k);
                %     parBase.u_ctrl = uk;
                %     % ----- plant step -----------------------------------------------
                %     % x    = obj.integrate(x,uk,obj.dt);      
                %     uPrev = uk;
                % end
                

                %── runtime update of gust & control ---------------------
                parNL = parBase;             % cheap struct copy
                % parNL.gust  = gustSignal(k, gust_input);  % plug if available
                parNL.gust  = gust_input(k);  % plug if available
                obj.log.wTrue(:,k) = gust_input(k);    

                % parNL.u_ctrl = controlSignal(k);   % plug if available
                if storeSens
                    Sn = S(:,:,k);                        % STM at current node
                end

                % ---------- stage 1  (explicit) ------------------------------
                k1N = obj.Nonlinear(x_n, obj.idx, parNL)*dt;
                if storeSens
                [Nq1,Nu1] = AeroFlex.sim.nonlinearJacobian(x_n, obj.idx, parNL);
                % [Nq1,Nu1] = AeroFlex.sim.buildNonlinearJacobian(x_n, obj.idx, parNL);

                    hN1 = Nq1*Sn*dt;
                    j1N = (Nq1*Sn+Nu1)*dt;
                    rhsS1 = obj.L *(Sn + gamma*j1N)*dt;
                    j1L  = obj.solveImplicit(rhsS1);          % (I-γΔtL)⁻¹ rhs
                    S1 = Sn + gamma*( j1L + j1N );
                end
                % ---------- stage 1  (implicit) ------------------------------
                rhs = obj.L*( x_n + gamma*k1N)*dt;

                k1L   =  obj.solveImplicit(rhs);          % y1 = x + γΔt L y1
                y1 = x_n + gamma*(k1L + k1N);

                % ---------- stage 2  (explicit) ------------------------------
                k2N = obj.Nonlinear( y1, obj.idx, parNL)*dt;

                if storeSens
                    [Nq2,Nu2] = AeroFlex.sim.nonlinearJacobian(y1, obj.idx, parNL);
                    % [Nq2,Nu2] = AeroFlex.sim.buildNonlinearJacobian(y1, obj.idx, parNL);

                    j2N = dt * ( Nq2*S1 + Nu2 );
                    s_for_j2L = Sn + (1-gamma)*j1L+delta*j1N + (1-delta)*j2N;
                    rhsS2  = obj.L * s_for_j2L*dt;
                    j2L  = obj.solveImplicit(rhsS2);
                    S2   = Sn + (1-gamma)*( j1L ) + gamma*( j2L )+delta*j1N + (1-delta)*j2N;
                end

                % ---------- stage 2  (implicit) ------------------------------
                x_for_k2L = x_n + (1-gamma)*k1L + delta*k1N + (1-delta)*k2N;
                rhs2        = obj.L*x_for_k2L*dt;
                k2L   = obj.solveImplicit(rhs2);
        
                % ------- Stage 3 (explicit) -----------------------------
                y2 = x_n + (1-gamma)*k1L + gamma*k2L + delta*k1N + (1 - delta)*k2N;
                k3N = obj.Nonlinear(y2, obj.idx, parNL)*dt;
                if storeSens
                    [Nq3,Nu3] = AeroFlex.sim.nonlinearJacobian(y2, obj.idx, parNL);
                    % [Nq3,Nu3] = AeroFlex.sim.buildNonlinearJacobian(y2, obj.idx, parNL);
                    
                    j3N = dt * ( Nq3*S2 + Nu3 );          % dt attached
                    S_np1  = Sn + (1-gamma)*(j1L+j2N)+(gamma)*(j2L+j3N);
                end
                % ------- Combine ---------------------------------------  
                x_np1 = x_n + ((1-gamma)*(k1L + k2N) + gamma*(k2L + k3N));                
               
                if storeSens
                    S(:,:,k+1) = S_np1;
                    cfg.buff.S(:,:,k+1) = S_np1;
                    obj.log.Scond(k+1) = cond(S_np1);
                    if k == 1      % print only once
                        fprintf('‖∂x/∂x0‖   %g\n', norm(S(:,       1:144), 'fro'));
                        fprintf('‖∂x/∂u ‖   %g\n', norm(S(:, 144+1), 'fro'));
                        fprintf('‖∂x/∂w ‖   %g\n', norm(S(:, 144+2:end), 'fro'));
                    end
                    currTest = S_np1 - Sn;

                end
                    
                if any(isnan(x_np1)) || any(isinf(x_np1))
                    error('Nonlinear terms exploded at iteration %d', k);
                end  
                disp(norm(x_np1));
                X(:,k+1) = x_np1;
                
                % ----- every Ts seconds call high‑level layer ---------------------
                if ~strcmp(obj.cfg.case, 'openLoop')
                % cfg.buff.X(:,kSam) = x_n;
                % cfg.buff.X(:,k) = x_np1;

                if mod(k, Nint) == 0
                    kSam          = kSam + 1;
                    cfg.buff.X(:,kSam) = x_n;

                    % before writing the new sample ------------------------
                    kBuff = kSam;
                    prevIdx = mod(kBuff-2 , cfg.ctrl.Ne+1) + 1;
                    if kBuff == 1          % wrapped – copy previous state into slot 1
                        cfg.buff.X(:,kBuff)   = cfg.buff.X(:,prevIdx);
                        cfg.buff.S(:,:,kBuff) = cfg.buff.S(:,:,prevIdx);
                        cfg.buff.U(:,kBuff)   = cfg.buff.U(:,prevIdx);
                        cfg.buff.t(kBuff)     = cfg.buff.t(prevIdx) + cfg.sim.dt;
                    end
                    % now write Y,U,T as usual

                    % cfg.buff.Xn1(:,kSam) = x_np1;
                    cfg.buff.U(:,kSam) = uPrev;
                    % zk = obj.sensor.measure(x_np1(obj.idx.q1));           % v_y at 5 positions
                    zk = obj.sensor.measure(x_np1(obj.idx.q1));           % v_y at 5 positions
                    cfg.buff.Y(:,kSam) = zk;
                    cfg.buff.t(kSam)   = t(k+1);
                    cfg.buff.valid     = min(kSam, cfg.ctrl.Ne);
                    % cfg.buff.valid = min(cfg.buff.valid+1 , cfg.ctrl.Ne+1);   % N = max(Ne,Nc)+1

            
                    % ---------- estimation & control ------------------------------
                    
                    [xhat, what] = estimator.estimate(zk, uPrev, t(k+1), cfg.buff);
                    % ---- quick probe ------------------------------------------------------
                    k  = estimator.idxWrap(-1);          % newest sample  (k = 0)
                    km = estimator.idxWrap(-2);          % previous sample (k = -1)
                    
                    x_now = estimator.buff.X(:,k);       % state at t(k)
                    x_prev= estimator.buff.X(:,km);      % state at t(k-1)
                    Sw     = estimator.buff.S(:,estimator.nx+1 , km);  % 144×1 disturbance column
                    wk     = estimator.what(end);        % last (oldest) ŵ  (= w_{k-1})
                    
                    err_state = norm( x_now - x_prev - Sw*wk );
                    disp(['continuity residual  ||x_k - f(x_{k-1},w_{k-1})|| = ',num2str(err_state)])
                    idxWin = estimator.idxWrap(-(0:cfg.ctrl.Ne));   % e.g. [4 3 2 1 7 6 5]
                    idxWin = idxWin(end:-1:1);                      % sort to chronological
                    if kSam ==1
                    figure
                    end
                    plot(cfg.buff.t,cfg.buff.X(1,:),'k.-'), hold on
                    plot(cfg.buff.t(idxWin), cfg.buff.X(1,idxWin),'ro')
                    title('Ring-buffer: black = full, red = window Ne')

                    % uCmd = controller.computeControl(xhat, t(k+1), cfg.buff)
                    % temp 
                    uCmd = zeros(4,1);
                    obj.log.U(:,kSam) = uCmd;     % store control at sample rate
                    obj.log.Y(:,kSam) = zk;
                    obj.log.wEst(:,kSam) = what( end ) ;   % keep only w_k (newest)
                    obj.log.xhat(:,kSam) = xhat;
                    uPrev = uCmd;                 % zero‑order hold for next window
                end
                end
                % ---------- diagnostics ------------------------------------
                                   
                % --- logging ----------------------------------------------------
                obj.log.x(:,k+1)  = x_np1;
                obj.log.y(:,k+1)  = yk;
                % if ~strcmp(obj.cfg.case, 'openLoop')
                %     obj.log.u(:,k)  = uk;
                %     obj.log.sqp(k)  = obj.controller.lastSolveInfo; % cost, #it, etc.
                % end
                
            end
            log = obj.log;
        end

        % ------------------------------------------------------------------
        function y = solveImplicit(obj,rhs)
            %SOLVEIMPLICIT  Solve (I-γΔt L) y = rhs using LU factors.
            %   For LU with 'vector' option we have P*A = L*U, where P is
            %   the permutation matrix defined by piv.  To solve Ay=rhs:
            %       z = L \ (P*rhs) ;  y = U \ z.
            z = rhs(obj.piv,:);           % apply permutation P*rhs
            z = obj.Lfac \ z;            % forward substitution (L)
            y = obj.Ufac \ z;            % backward substitution (U)
        end

        function f = Nonlinear(obj ,x, idx, parNL)
            % Wrapper around fast vectorised intrinsic‑beam nonlinear terms.
            % Nm  = obj.beam.Nm;     Na2 = obj.cfg.Na;
            f   = AeroFlex.sim.nonlinear_terms(x, parNL, idx);
            % f = f/obj.cfg.sim.dt;   % convert to derivative form
        end

        function J = NJacobian(obj,x, idx, parNL)
            % h  = 1e-6*max(1,norm(x)); n = numel(x); fx = obj.Nonlinear(x); J = zeros(n);
            % for i = 1:n
            %     xp = x; xp(i) = xp(i)+h; J(:,i) = (obj.Nonlinear(xp)-fx)/h;
            % end
            %NJACOBIAN Finite‑difference Jacobian of Nonlinear() –
            % expensive but only used when sensitivities are requested.
            h  = 1e-6; n = numel(x);
            fx = obj.Nonlinear(x, idx, parNL);
            J  = zeros(n);
            for i = 1:n
                xp     = x; xp(i) = xp(i)+h;
                J(:,i) = (obj.Nonlinear(xp, idx, parNL) - fx)/h;
            end
        end

        % function gust_k = gustSignal(k, gust)
        %     gust_k = gust(k);
        % end
        
    end
end
