classdef nMHEv3 < AeroFlex.ctrl.EstimatorBase
    properties
       % ---- scalar config
        dt  double
        Ts  double
        Ne  double
        Nm
        idx
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
        % ---- outputs
        % xhat double
        % what double
        % ---- debug
        debug logical = true
        dbg   struct
        
        % Buffers for past data and previous solution (for warm start)
        yBuffer          % Measured outputs buffer (cell or array of length N)
        uBuffer          % Input buffer (length N)
        % Storage for last optimized trajectory (for warm start)
        lastStateSeq     % [nx x (N+1)] state trajectory from last solution
        lastDistSeq      % [nw x N] disturbance sequence from last solution
        xRef        % Reference prior state for current horizon start (for equality constraint)
        wRef        % Reference disturbance for horizon initial step (for cost penalization)
    end
    methods
        function obj = nMHEv3(cfg,beam,aero,base)
            % Constructor: initialize model, horizon, weights, bounds, etc.
        %----------------------------------------------------------------------

            obj@AeroFlex.ctrl.EstimatorBase(cfg);
        
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
            obj.Nm = cfg.Nm;
            % -------- 2. weights / bounds ------------------------------------
            obj.Qe = cfg.Qe;     obj.Pe = cfg.Pe;     obj.Re = cfg.Re;
            obj.xL = cfg.xL;     obj.xU = cfg.xU;
            obj.wL = cfg.wL;     obj.wU = cfg.wU;
            
            % -------- 3. optimiser options -----------------------------------
            obj.solverOpts = optimoptions('fmincon',...
                'Algorithm','interior-point',... % interior-point
                'SpecifyObjectiveGradient',true,...
                'SpecifyConstraintGradient',true,...
                'HessianApproximation','lbfgs',...
                'Display','iter-detailed',...
                'FiniteDifferenceStepSize',1e-8,...
                'CheckGradients',true,...
                'OptimalityTolerance',1e-5,...
                'StepTolerance',1e-6);
            obj.solverOpts.FiniteDifferenceType = 'central';
            % -------- 4. rolling buffers  ------------------------------------
            obj.Xhist = zeros(obj.nx,obj.Ne+1);
            obj.Uhist = zeros(obj.nu,obj.Ne);
            obj.Yhist = zeros(obj.ny,obj.Ne+1);
            obj.Whist = zeros(obj.nw,obj.Ne);
            obj.k     = 0;
        
            % -------- 5. initial guess ---------------------------------------
            xTrim           = cfg.trim.states;        % equilibrium state
            obj.xhat        = xTrim;
            obj.Xhist(:,:)  = repmat(xTrim,1,obj.Ne+1);
            obj.Whist(:,:)  = 0;                      % zero-gust history
            obj.what = zeros(obj.Ne*obj.nw,1);
        
            obj.idx   = obj.buildIndexMaps();
            obj.z0           = zeros((obj.Ne+1)*obj.nx + obj.Ne*obj.nw,1);
            % obj.z0(idx.x{1}) = xTrim;                 % arrival state
            obj.H            = speye(numel(obj.z0));  % LBFGS initial Hessian
        
            % initial STM (∂x/∂x0  and  ∂x/∂w) – identity / zeros
            obj.Sprev = [eye(obj.nx) , zeros(obj.nx,obj.nw+obj.nu)];
            obj.dbg = struct('t', [], 'W', [], 'cont', 0);

            % obj.xBounds = xBounds;   % expect an [nx x 2] matrix [lower, upper]
            % obj.wBounds = wBounds;   % [1 x 2] for single disturbance
            % Initialize buffers
            obj.yBuffer = {}; 
            obj.uBuffer = {};
            obj.lastStateSeq = []; 
            obj.lastDistSeq = [];
            % % Initialize prior references (to be set when first data available)
            obj.xRef = xTrim;        % initial state prior (assumed known initial condition)
            obj.wRef = 0;              % initial disturbance reference (assume zero initially)
        end

        function [xhat, what, info] = estimate(obj, y_meas, u_prev, t_k) 
            % ESTIMATE: perform one MHE solve given new measurement and input.
            % Append the latest measurement and input to buffers
            obj.yBuffer{end+1} = y_meas;
            obj.uBuffer{end+1} = u_prev;
            if length(obj.yBuffer) > obj.Ne
                % If buffer exceeds horizon length, drop the oldest entry (rolling window)
                obj.yBuffer(1) = [];
                obj.uBuffer(1) = [];
            end
            idx = obj.idx;
            % current arrival state -------------------------
            Nint=round(obj.Ts/obj.dt);
            % If not enough data to fill one horizon window, return initial state (no estimation yet)
            if length(obj.yBuffer) < obj.Ne

            % if obj.buff.valid < obj.Ne + 1
                % Use linear output model if available (C * q1) for quick prediction
                C     = obj.sensor.PhiY;
                % Predict output from current state estimate (assuming output depends on q1 subset of state)
                q1_idx = idx.q1;
                obj.yhat = C * obj.xhat(q1_idx);
                xhat     = obj.xhat;
                what     = zeros(obj.Ne * obj.nw, 1);  % no disturbance estimate yet
                info = struct('cost', [], 'exitflag', [], 'continuity', [], 'grad_err', []);
                return;
            end
        %%
            % If we have a previous solution, use it to shift the guess (warm start) for the new window
            if ~isempty(obj.z0)
                obj.shiftGuess();   % shift obj.z0, and update reference state/disturbance (obj.xRef, obj.wRef)
            else
                % First time solving: initialize decision vector guess [x_E; w_history]
                obj.z0 = [ repmat(obj.xhat, obj.Ne, 1); zeros(obj.nw * obj.Ne, 1) ];d
                obj.xRef = obj.xhat;                 % initial reference state
                obj.wRef = zeros(obj.nw, 1);               % initial reference disturbance
            end
        %%
            % obj.z0 = obj.getInitialGuess();
        %%
            % Assemble the multiple-shooting NLP (cost and constraints with gradients)
            nlp = obj.assembleWindow();

            % Solve the nonlinear program using SQP (fmincon)
            % [p, fval, exitFlag] = fmincon(nlp.cost, obj.z0, [], [], [], [], nlp.lb, nlp.ub, nlp.nonl, obj.solverOpts);
            %%
            % 
            % 
            % % Setup bounds for decision vector: [states(1..N), disturbances(0..N-1)]
            % % There are N unknown state vectors (excluding the fixed initial state) and N disturbance scalars.
            % nx = obj.nx; Ne = obj.Ne;
            % % State bounds (for states 1 to N):
            % x_lower = repmat(obj.xBounds(:,1), Ne, 1);
            % x_upper = repmat(obj.xBounds(:,2), Ne, 1);
            % % Disturbance bounds (for w0 to w_{N-1}):
            % w_lower = repmat(obj.wBounds(1), Ne, 1);
            % w_upper = repmat(obj.wBounds(2), Ne, 1);
            % % Combine into vector form [state1..stateN; w0..w_{N-1}]
            % lb = [x_lower; w_lower];
            % ub = [x_upper; w_upper];
            % 
            % % Define objective and constraint function handles for fmincon
            % funObj = @(z) obj.costFunction(z);
            % funCon = @(z) obj.constraintFunction(z);
            % 

            try
                [p, fval, exitFlag] = fmincon(nlp.cost, obj.z0, [], [], [], [], nlp.lb, nlp.ub, nlp.nonl, obj.solverOpts);            
            catch solveError
                error('nMHEv2: NLP solver error: %s', solveError.message);
            end
            if exitFlag <= 0
                warning('nMHEv2: fmincon did not converge (exit flag %d). Results may be invalid.', exitFlag);
                p = obj.z0;

            end
            
            % Unpack solution into state trajectory and disturbances
            [x_seq, w_seq] = obj.unpackSolution(p);
            % The estimated current state is the last state in the trajectory (time k)
            xhat = x_seq(:, end);
            % Also return the estimated disturbance sequence (if needed)
            what = w_seq;
            
            % Update prior references for next horizon:
            % Use the state and disturbance at the *start of the second interval* of this solution as new references:contentReference[oaicite:19]{index=19}.
            % That corresponds to state at index 2 (i.e., time k-N_e+1) and disturbance at interval 1.
            % obj.xRef = x_seq(:, 2);   % state at horizon index 1 (start of second sub-interval)
            % if obj.nw > 0
            %     obj.wRef = w_seq(1);  % disturbance of second sub-interval (index 1, since w_seq indexed from 0)
            % end
            
            % Store the solution for warm start in subsequent calls
            obj.lastStateSeq = x_seq;
            obj.lastDistSeq = w_seq;
            
            % Debug logging: residual norm, sensitivity magnitudes, continuity residuals
            if obj.debug
            %     % Compute measurement residuals (y_predicted - y_measured) over the horizon
            %     [y_pred_seq] = obj.predictOutputs(x_seq);  % predicted outputs from estimated states
            %     % Form residual vector for all measured outputs in the window
            %     res_all = [];
            %     for i = 1:length(obj.yBuffer)
            %         y_pred_i = y_pred_seq{i};
            %         y_meas_i = obj.yBuffer{i};
            %         res_all = [res_all; (y_pred_i - y_meas_i)]; %#ok<AGROW>
            %     end
            %     res_norm = norm(res_all);  % 2-norm of all residuals
            %     % Compute sensitivity matrices and continuity residuals for each interval
            %     maxSens = zeros(1, Ne);        % max |element| of sensitivity matrix per interval
            %     contResiduals = zeros(1, Ne);  % norm of continuity constraint residual per interval
            %     % Note: state sequence x_seq includes initial state as x_seq(:,1).
            %     % We enforce x_seq(:,1) = xRef, so continuity from that to x_seq(:,2) uses w_seq(1).
            %     for k = 1:Ne  % k=1 corresponds to interval 0 (initial->1), k=N corresponds to last interval
            %         if k == 1
            %             x_start = obj.xRef;       % known initial state
            %             x_end   = x_seq(:, k+1);       % x_seq(:,2)
            %             w_k     = w_seq(k);            % w_seq(1) corresponds to interval 0
            %             u_k     = obj.uBuffer{1};      % first control in buffer
            %         else
            %             x_start = x_seq(:, k);         % state at interval start (k index in seq corresponds to start of interval k, since seq(1)=initial)
            %             x_end   = x_seq(:, k+1);
            %             w_k     = w_seq(k);
            %             u_k     = obj.uBuffer{k};      % corresponding control
            %         end
            %         % Continuity residual (should be ~0 if constraint satisfied)
            %         x_pred = obj.propagateState(x_start, u_k, w_k);
            %         cont_res_vec = x_end - x_pred;
            %         contResiduals(k) = norm(cont_res_vec, 2);  % norm of difference
            %         % Sensitivity: compute Jacobian of f wrt state at this interval
            %         J = obj.computeStateSensitivity(x_start, u_k, w_k);
            %         maxSens(k) = max(max(abs(J)));  % maximum absolute entry in Jacobian
            %     end
                % % Log the results
                % fprintf('MHE Debug: Residual norm = %.3e, Max |S| per step = [', res_norm);
                % fprintf('%.2e ', maxSens);
                % fprintf('], Continuity residuals = [');
                % fprintf('%.1e ', contResiduals);
                % fprintf(']\n');
                %%
                % % Check solver result and update estimates
                % if exitFlag <= 0
                %     warning('nMHE:estimate', 'Solver did not converge (flag %d), keeping previous estimate.', exitFlag);
                %     % If solver fails, do not update state/disturbance estimates (use last)
                %     p = obj.z0;
                % end
                % Extract estimated current state (latest node) and disturbance history from solution
                idx = obj.idx;  % index mapping for decision vector
                obj.xhat = p( idx.x{obj.Ne} );            % estimated state at time t (k=0)
                % obj.yhat = z;                               % use current measurement as "predicted" output (for output)
                obj.what = p( [idx.w{:}] );                % estimated disturbance history over window (concatenated)
            
                % Store solution for warm start next time
                obj.z0 = p;
            
                % Update reference for next horizon (per Artola et al., update to state/disturbance at second interval start)
                obj.xRef = p( idx.x{2} );                  % reference state = state at start of 2nd sub-interval (previous solution)
                if obj.nw > 0
                    obj.wRef = p( idx.w{2} );              % reference disturbance = disturbance in second-oldest interval
                end
            
                % Prepare info output
                info = struct();
                info.cost       = fval; 
                info.exitflag   = exitFlag;
                % Compute continuity residuals for each interval (norm of equality constraint violation)
                ceq = obj.ceqOnly(p);
                info.continuity = zeros(obj.Ne, 1);
                for j = 1:obj.Ne
                    % Norm of state continuity residual on interval j
                    rows = (j-1)*obj.nx + (1:obj.nx);
                    info.continuity(j) = norm( ceq(rows) );
                    if obj.debug
                        fprintf('Interval %d continuity residual norm = %.3e\n', j, info.continuity(j));
                    end
                end
            
                % % Sensitivity matrix diagnostics: check condition numbers of state STM for each interval (from buffers)
                % if obj.cfg.debug
                %     for j = 1:obj.Ne
                %         Sj = obj.buff.S(:, 1:obj.nx, obj.idxWrap(-obj.Ne + j - 1));  % STM for interval j (∂x_{j+1}/∂x_j)
                %         condSj = cond(Sj);
                %         fprintf('Interval %d STM condition = %.3e\n', j, condSj);
                %     end
                % end
            
                % Finite-difference vs analytical gradient check for one variable (debug)
                if obj.debug
                    % Choose a representative decision variable (e.g., first element of state vector) to perturb
                    kVar = 1;
                    delta = 1e-6 * (1 + abs(p(kVar)));
                    dp    = zeros(size(p)); dp(kVar) = delta;
                    ceq_plus  = obj.ceqOnly(p + dp);
                    ceq_minus = obj.ceqOnly(p - dp);
                    fd_col = (ceq_plus - ceq_minus) / (2*delta);           % finite-difference approximation of gradient column
                    % Get analytical constraint gradients at solution
                    [~, ~, ~, gradceq] = nlp.nonl(p);
                    analytic_col = gradceq(kVar,:);                       % column corresponding to variable kVar
                    grad_err = norm(fd_col - analytic_col) / max(1, norm(fd_col));
                    fprintf('Gradient check for variable %d: rel. error = %.3e\n', kVar, grad_err);
                    info.grad_err = grad_err;
                end
            
                % Return estimates
                xhat = obj.xhat;
                what = obj.what;
            
                %%
            end

        end
            function nlp = assembleWindow(obj)
                % Assemble the MHE optimization window (objective and constraints with gradients)
            
                % Pre-compute index mappings for decision vector components
                idx = obj.idx;  % struct with idx.x{1...Ne+1} and idx.w{1...Ne}
                nx = obj.nx;  nw = obj.nw;  Ne = obj.Ne;
                ny = obj.ny;

                % Precompute full output matrix (for computing h(x) quickly)
                persistent Cfull
                if isempty(Cfull)
                    Cfull = zeros(ny, nx);
                    Cfull(:, obj.idx.q1) = obj.sensor.PhiY;  % embed sensor matrix into full state dimension
                end

                % % ** Nested cost function ** (objective J and its gradient)
                function [J, gradJ] = costFun(p)
                    % Extract state and disturbance variables from p
                    % X_mat is (nx x (Ne+1)) matrix of all state nodes, w_hist is (Ne x nw) matrix of all disturbances
                    % X_mat = reshape( p([idx.x{:}]), nx, Ne+1 );
                    X_mat = reshape( p([idx.x{:}]), nx, Ne );
                    if nw > 0
                        w_hist = reshape( p([idx.w{:}]), nw, Ne ).';   % size (Ne x nw), each row j is w_j (oldest j=1)
                    else
                        w_hist = zeros(Ne, 0);
                    end

                    % Arrival cost (penalize deviation of oldest state and oldest disturbance from references)
                    % x_E = x_{0|0} (oldest state), w_old = w_{-Ne} (oldest disturbance)
                    % References: obj.xRef (q_{0|0,A}) and obj.wRef (w_A) updated from previous solution
                    J = 0;
                    gradJ = zeros(numel(p), 1);
                    % State arrival term: 0.5 * (x_E - xRef)^T Pe * (x_E - xRef)

                    xE = X_mat(:, 1);
                    dx = xE - obj.xRef;
                    J = J + 0.5 * (dx.' * obj.Pe * dx);
                    gradJ(idx.x{1}) = obj.Pe * dx;  % gradient w.rt x_E

                    % Disturbance arrival term: 0.5 * (w_old - wRef)^T Re * (w_old - wRef)
                    if nw > 0 
                        w_old = p(idx.w{1});
                        dw = w_old - obj.wRef;
                        J = J + 0.5 * (dw.' * obj.Re * dw);
                        gradJ(idx.w{1}) = obj.Re * dw;  % gradient w.rt oldest w
                    end

                    % Output measurement error terms (for each k = -Ne,...,0)
                    % Loop backward through time to use available S in buffer (optional; here we compute directly)
                    % for j = 1:Ne+1
                    for j = 1:Ne
                        % Predicted output at node j and measurement error
                        y_pred = Cfull * X_mat(:, j);
                        % Measurement corresponding to state j (note: idxWrap maps relative index to buffer index)
                        % z_meas = obj.buff.Y(:, obj.idxWrap(j - (Ne+1)));
                        z_meas = obj.yBuffer{j}.';
                        e = z_meas.' - y_pred;                             % innovation e_k = y_k(actual) - y_k(predicted)
                        % Add contribution to cost
                        J = J + 0.5 * e.' * obj.Qe * e;
                        % Gradient w.rt state x_j: - (∂y/∂x)^T Qe * e  (since y_pred = Cfull * x)
                        gradJ(idx.x{j}) = gradJ(idx.x{j}) - Cfull.' * (obj.Qe * e);
                    end

                    % (Optional) Disturbance smoothness regularization: minimize ∑‖Δw‖^2 
                    % if obj.alphaSmooth > 0
                    %     dW = diff(w_hist, 1, 1);                   % (Ne-1 x nw) differences between successive w's
                    %     J = J + 0.5 * obj.alphaSmooth * sum(dW(:).^2);
                    %     % Gradient of smoothness term: contributes ±alphaSmooth * (w_j - w_{j-1}) for each interior w
                    %     for j = 2:Ne
                    %         dw = w_hist(j,:) - w_hist(j-1,:);
                    %         gradJ(idx.w{j})   = gradJ(idx.w{j})   + obj.alphaSmooth * dw.';
                    %         gradJ(idx.w{j-1}) = gradJ(idx.w{j-1}) - obj.alphaSmooth * dw.';
                    %     end
                    % end

                    % (Optional) Barrier terms to keep disturbances positive (if needed)
                    % muBar = obj.muBarrier;  betaBar = obj.betaBarrier;
                    % if muBar > 0
                    %     expTerm = exp( - betaBar * w_hist );         % elementwise
                    %     J = J + muBar * sum(expTerm(:));
                    %     gradJ([idx.w{:}]) = gradJ([idx.w{:}]) - muBar * betaBar * expTerm(:);
                    % end
                    %% Other
                    % [x_vars, w_vars] = obj.unpackSolution(p);
                    % % Note: x_vars includes states 1..N (initial state is fixed and not in z).
                    % % Compute output prediction errors and process noise penalties
                    % % Predicted outputs for each measurement in the window:
                    % y_pred_seq = obj.predictOutputs([obj.xRef, x_vars]);  % include the fixed initial state at start
                    % % Assemble cost: measurement residuals + process noise deviations
                    % % Measurement residual cost: 0.5 * sum (y_meas - y_pred)' Qe (y_meas - y_pred)
                    % J = 0;
                    % for i = 1:obj.Ne  % there are N measurements in buffer (assuming full horizon of data)
                    %     y_pred = y_pred_seq{i};     % predicted output at interval i (buffer index i)
                    %     y_meas = obj.yBuffer{i}.';    % actual measurement
                    %     e = (y_meas - y_pred);
                    %     % Cost contribution (treat Qe as weighting matrix)
                    %     J = J + 0.5 * (e' * obj.Qe * e);
                    % end
                    % % Prior state term: if initial state is not exactly fixed, penalize deviation from reference
                    % % (In this implementation initial state is equality-constrained, so this term is optional)
                    % % e_x0 = obj.xPriorRef - obj.xPriorRef (always zero due to hard constraint)
                    % % J = J + 0.5 * (e_x0' * P * e_x0);   % P can be taken as inverse covariance of prior estimate
                    % 
                    % % Process noise penalty: 0.5 * sum (w_k - w_ref_k)' Rw (w_k - w_ref_k)
                    % % We include a reference for initial disturbance (from previous window) and assume 0 reference for subsequent.
                    % for k = 1:obj.Ne
                    %     w_k = w_vars(k);
                    %     if k == 1
                    %         w_ref = obj.wPriorRef;   % reference for initial interval disturbance
                    %     else
                    %         w_ref = 0;               % assume no expected disturbance in later intervals
                    %     end
                    %     dw = (w_k - w_ref);
                    %     J = J + 0.5 * (dw' * obj.Rw * dw);
                    % end
                    %%
                end  % costFun
            
                % ** Nested nonlinear constraint function ** (continuity constraints and their Jacobian)
                function [c, ceq, gradc, gradceq] = nonlFun(p)
                    % No inequality constraints in MHE (all bounds handled via lb/ub)
                    c = [];
                    % Allocate equality constraint vector and Jacobian
                    ceq = zeros(obj.nx * obj.Ne, 1);
                    % gradceq = zeros(numel(p), obj.nx * obj.Ne);
                    gradceq = zeros(obj.nx * obj.Ne, numel(p));
                    S_interval = [ eye(obj.nx), zeros(obj.nx, obj.nw+obj.nu) ];
                    S0 = S_interval;
                    % Loop through each interval j = 1,...,Ne to enforce x_{j+1} = F(x_j, u_j, w_j)
                    % (multiple shooting continuity constraints)
                    x_curr = p( idx.x{1} );  % start with oldest state (x_1 in decision vector, corresponds to t=-Ne)
                    % x_curr = p( idx.x{2} );  % start with oldest state (x_1 in decision vector, corresponds to t=-Ne)
                    % for j = 1:Ne
                    % for j = 1:Ne-1
                    % for j = 2:Ne-1
                    % for j = Ne-1:-1:1
                    for j = 1:Ne-2
                        % Identify control applied over interval j (assumed known and stored in buffer)
                        % u_j corresponds to interval from time index (horizon start + j-1) to (horizon start + j)
                        % disp(j)
                        % u_j = obj.uBuffer(:, obj.idxWrap(-Ne + j));  
                        u_j = cell2mat(obj.uBuffer(:, j));  
                        % Disturbance for this interval (decision variable)
                        w_j = (nw > 0) * p( idx.w{j} );   % if nw==0, this yields 0; if nw>0, the actual scalar/vector

                        % Integrate dynamics over Ts (divided into Nint smaller steps) to get predicted next state and STM
                        x_pred = x_curr;
                        % Initialize STM for this interval (nx x (nx+nw)): identity for state part, zero for disturbance part
                        S_interval = [ eye(obj.nx), zeros(obj.nx, obj.nw+obj.nu) ];
                        % Integrate in Nint steps using ROMIntegrator (with STM propagation)
                        Nint = round(obj.Ts / obj.dt);
                        for k = 1:Nint
                            [x_pred, S_interval] = obj.model.step(x_pred, u_j, w_j, S_interval, true);
                        end

                        Sx = S_interval(:, 1:obj.nx);
                        Sw = S_interval(:, obj.nx+1 : obj.nx+nw);
                        % Continuity constraint: x_next (decision var) minus predicted x_next must be zero
                        x_next_var = p( idx.x{j+1} );           % the (j+1)-th state in decision vector
                        % x_next_var = p( idx.x{j+2} );           % the (j+1)-th state in decision vector
                        
                        ceq_block = x_next_var - x_pred;

                        % rows = (j-1)*obj.nx + (1:obj.nx);
                        rows = (j)*obj.nx + (1:obj.nx);
                        ceq(rows) = ceq_block;
                        % disp(rows(1));
                        % Compute Jacobian for this interval's constraints:
                        % ∂ceq_block/∂x_{j+1} = I, ∂ceq_block/∂x_j = -∂x_pred/∂x_j, ∂ceq_block/∂w_j = -∂x_pred/∂w_j
                        gradceq(rows, idx.x{j+1}) =  eye(obj.nx);
                        % gradceq(rows, idx.x{j+1}) =  -eye(obj.nx);
                        % gradceq(rows, idx.x{j+1}) =  Sx;
                        gradceq(rows, idx.x{j}) = -Sx;
                        % gradceq(rows, idx.x{j}) = -eye(obj.nx);
                        if nw > 0
                            gradceq(rows, idx.w{j}) = -Sw;
                        end

                        % gradceq(idx.x{j+1}, rows) =  eye(obj.nx);
                        % gradceq(idx.x{j},   rows) = -Sx.';
                        % if nw > 0
                        %     gradceq(idx.w{j}, rows) = -Sw.';
                        % end

                        % Prepare for next interval
                        x_curr = x_next_var;
                        %%
                        % [x_vars, w_vars] = obj.unpackSolution(z);
                        % % Note: x_vars is [state1,...,stateN] (each state is nx-length), and 
                        % % w_vars = [w0,...,w_{N-1}] (each w scalar here).
                        % % We build ceq as a concatenation of all continuity constraint residuals.
                        % N = obj.Ne; nx = obj.nx;
                        % ceq = zeros(nx * N, 1);
                        % % Loop over each interval k = 0...N-1
                        % for k = 0:(N-1)
                        %     if k == 0
                        %         % Interval 0: from known initial state xPriorRef to x_vars(:,1) using w_vars(1)
                        %         x_start = obj.xPriorRef;
                        %         x_target = x_vars(:, 1);
                        %         w_k = w_vars(1);
                        %         u_k = obj.uBuffer{1};
                        %     else
                        %         % Interval k: from x_vars(:,k) to x_vars(:,k+1) using w_vars(k+1)
                        %         x_start = x_vars(:, k);       % state at beginning of interval k
                        %         x_target = x_vars(:, k+1);    % state at end of interval k
                        %         w_k = w_vars(k+1);
                        %         u_k = obj.uBuffer{k+1};
                        %     end
                        %     % Propagate state from x_start over one sample period Ts
                        %     x_pred = obj.propagateState(x_start, u_k, w_k);
                        %     % Continuity residual = x_target - x_pred (should equal zero)
                        %     idx = (k * nx + 1):(k * nx + nx);
                        %     ceq(idx) = x_target - x_pred;
                        % end
                        % % Additionally, enforce explicitly the initial equality constraint x0 = xPriorRef.
                        % % (Since we do not include x0 in decision variables, this is inherently satisfied.)
                        % % If x0 were part of z, we would add: ceq_init = x0_var - obj.xPriorRef.
                        %%
                    end
            
                    % Transpose gradient for fmincon requirements (constraints×variables)
                    gradceq = gradceq.';
                    gradc   = [];  % no inequality constraints
                end  % nonlFun
            
                % Define variable bounds (linear inequality constraints not used since bounds suffice)
                % State bounds apply to each state node, disturbance bounds to each w in history
                % xLrep = repmat(obj.xL, obj.Ne+1, 1);
                % xUrep = repmat(obj.xU, obj.Ne+1, 1);
                xLrep = repmat(obj.xL, obj.Ne, 1);
                xUrep = repmat(obj.xU, obj.Ne, 1);
                wLrep = repmat(obj.wL, obj.Ne, 1);
                wUrep = repmat(obj.wU, obj.Ne, 1);
            
                % Return NLP definition
                nlp.cost = @(p) costFun(p);
                nlp.nonl = @(p) nonlFun(p);
                nlp.lb   = [xLrep; wLrep];
                nlp.ub   = [xUrep; wUrep];
            end
            function idx = buildIndexMaps(obj)
                % Build index maps for state and disturbance components in the decision vector
                % idx.x{j} gives indices of x_j in the vector, for j=1...Ne+1
                % idx.w{j} gives indices of w_j in the vector, for j=1...Ne
                nx = obj.nx;  nw = obj.nw;  Ne = obj.Ne;
                % idx.x = cell(Ne+1, 1);
                idx.x = cell(Ne, 1);
                % for j = 1:Ne+1
                for j = 1:Ne
                    idx.x{j} = (j-1)*nx + (1:nx);
                end
                idx.w = cell(Ne, 1);
                % startW = (Ne+1)*nx;
                startW = (Ne)*nx;
                for j = 1:Ne
                    idx.w{j} = startW + (j-1)*nw + (1:nw);
                end
                idx.q1 = 1:obj.Nm; 
                idx.q2 = obj.Nm+1:obj.Nm*2;
                idx.qxi = obj.Nm*2+1 : obj.Nm*3+1;
                % idx.qg = obj.Nm*3+1+1:
            end
            
            function shiftGuess(obj)
                % Shift the previous solution (obj.z0) to warm-start the next estimation window
                % Drops the oldest state/w from the window and adds an initial guess for the new interval at the end.
                if isempty(obj.z0)
                    return;
                end
                p_old = obj.z0;
                idx = obj.idx;
                Ne = obj.Ne; nw = obj.nw;
                % Construct new decision vector guess
                new_p = zeros(size(p_old));
                % Shift states: x^1_new (oldest) = old x^2, ..., x^Ne_new = old x^{Ne+1}
                for j = 1:Ne-1
                    new_p( idx.x{j} ) = p_old( idx.x{j+1} );
                end
                % For the new last state x^{Ne+1}, use model prediction from old last state (or repeat last state as guess)
                x_last_old = p_old( idx.x{Ne} );
                u_last = obj.uBuffer(:, end);    % control applied in last interval (from t=-1 to t=0)
                if nw > 0
                    w_last = p_old( idx.w{Ne} );          % last interval disturbance from previous solution
                else
                    w_last = [];
                end
                % Predict state one step ahead (from t=0 to t=Ts) for initial guess of new state at t=Ts
                x_guess = x_last_old;
                S = [eye(obj.nx), zeros(obj.nx, obj.nw+obj.nu)];
                Nint = round(obj.Ts / obj.dt);
                for k = 1:Nint
                    [x_guess, S] = obj.model.step(x_guess, u_last, zeros(size(w_last)), S, false);
                end
                % new_p( idx.x{Ne+1} ) = x_guess;
                new_p( idx.x{Ne} ) = x_guess;
                % Shift disturbances: w^1_new = old w^2, ..., w^{Ne-1}_new = old w^Ne
                for j = 1:Ne-1
                    new_p( idx.w{j} ) = p_old( idx.w{j+1} );
                end
                % New last disturbance w^Ne_new: assume decays to zero (vanishing disturbance beyond horizon)
                if nw > 0
                    new_p( idx.w{Ne} ) = 0;
                end
                % Update stored guess and references
                obj.z0 = new_p;
                % (obj.xRef and obj.wRef will be updated after solving, in estimate, using previous solution's second interval)
            end
            
            function ceq = ceqOnly(obj, p)
                % Compute only the equality constraint residual (continuity) for a given decision vector p
                % This is used for diagnostics (norms, gradient checks) outside of the optimizer.
                Ne = obj.Ne; nx = obj.nx; nw = obj.nw;
                idx = obj.idx;
                ceq = zeros(nx * Ne, 1);
                x_curr = p( idx.x{1} );
                for j = 1:Ne-1
                    % Use the same integration as in nonlFun to get x_pred
                    % u_j = obj.buff.U(:, obj.idxWrap(-obj.Ne + j));
                    u_j = obj.uBuffer(:, j);
                    if nw > 0
                        w_j = p( idx.w{j} );
                    else
                        w_j = [];
                    end
                    x_pred = x_curr;
                    Nint = round(obj.Ts / obj.dt);
                    for k = 1:Nint
                        [x_pred, ~] = obj.model.step(x_pred, u_j, w_j, [], false);
                    end
                    % Residual = x_next_var - x_pred
                    x_next_var = p( idx.x{j+1} );
                    ceq_block = x_next_var - x_pred;
                    ceq((j-1)*nx + (1:nx)) = ceq_block;
                    x_curr = x_next_var;
                end
            end
            function z0 = getInitialGuess(obj)
            % Generate a warm-start initial guess for the decision variables (states 1..N and w0..w_{N-1}).
            N = obj.Ne; nx = obj.nx;
            z0 = zeros(nx * N + obj.nw * N, 1);
            if ~isempty(obj.lastStateSeq)
                % Use previous solution shifted as initial guess
                % lastStateSeq is [x0, x1, ..., xN] from previous horizon
                % We will use [x1,...,xN] from last solution as guess for [x1,...,x_{N-1}] in new horizon.
                % And use last state's copy or one-step prediction for new x_N guess.
                x_last = obj.lastStateSeq;
                w_last = obj.lastDistSeq;
                % Shift state guess: use previous x(2)...x(N) as guess for new state1...state(N-1)
                for k = 1:(N-1)
                    x_guess = x_last(:, k+1);  % previous horizon state at k+1
                    % Place into z0 (states block)
                    idx = (k-1)*nx + 1;
                    z0(idx:idx+nx-1) = x_guess;
                end
                % Guess for new last state x_N: extrapolate from last known state
                % We propagate the last state using the latest input and assuming zero disturbance
                new_u = obj.uBuffer{end};      % newest input
                last_x_est = x_last(:, end);   % previous horizon final state
                xN_guess = obj.propagateState(last_x_est, new_u, 0);
                idx = (N-1)*nx + 1;
                z0(idx:idx+nx-1) = xN_guess;
                % Shift disturbance guess: use previous w(2)...w(N-1) as guess for new w0...w_{N-2}
                for k = 2:(N-1)
                    w_guess = w_last(k);  % previous horizon disturbance at interval k (0-indexed in w_last)
                    z0(nx*N + k - 1) = w_guess;
                end
                % Guess new last disturbance w_{N-1} as zero (assuming no new disturbance)
                z0(nx*N + N) = 0;
                % For initial w0 (new horizon's first interval): we can guess it as previous horizon's second interval disturbance
                if length(w_last) >= 2
                    z0(nx*N + 1) = w_last(2);
                else
                    z0(nx*N + 1) = w_last(1);
                end
            else
                % For states, use prior state as baseline
                Nint=round(obj.Ts/obj.dt);

                % forward-predict the horizon with zero gust ----
                x_tmp = obj.xhat;  S_tmp = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nu)];
                obj.z0(1:nx) = x_tmp;

                    for j = 1:obj.Ne-1
                        uj = obj.Uhist(:,j);  wj = obj.what(j);
                        for m = 1:Nint
                            [x_tmp,~] = obj.model.step(x_tmp, uj, wj, S_tmp, false);
                        end
                        % obj.z0(idx.x{j+1}) = x_tmp;
                        % idx = (j-1)*nx + 1;
                        idx = (j)*nx + 1;

                        obj.z0(idx:idx+nx-1) = x_tmp;
                    end
                % obj.Sprev = S_tmp;      
                
                % for k = 1:N
                %     idx = (k-1)*nx + 1;
                %     z0(idx:idx+nx-1) = obj.xRef;  % start with prior state as guess for all states
                % end
                % For disturbances, start with zero
                z0(nx*N + 1:end) = 0;
            end
        end
        
        function [x_seq, w_seq] = unpackSolution(obj, z)
            % UNPACKSOLUTION: split decision vector z into state sequence and disturbance sequence.
            % State sequence will include states 1..N (not including initial known state).
            % Disturbance sequence includes w0..w_{N-1}.
            N = obj.Ne; nx = obj.nx;
            % Extract state variables (concatenated)
            x_seq_vars = reshape(z(1:nx*N), nx, N);
            % Prepend the fixed initial state to reconstruct full trajectory (optional)
            % Here we do include the initial known state at the beginning for completeness:
            x_seq = [obj.xRef, x_seq_vars];
            % Extract disturbances
            w_seq = z(nx*N+1 : nx*N + N).';
            % w_seq will be a 1xN row (we transpose to row for convenience; it's scalar per interval here)
        end
        
        function x_next = propagateState(obj, x_current, u_current, w_current)
            % PROPAGATESTATE: propagate state x_current over one sampling period Ts using model.step.
            % This function ensures the model integrates over the full Ts (possibly using smaller dt steps).
            % If model.step integrates one internal time-step, loop until Ts is covered.
            x_next = x_current;
            % Determine internal time step
            if isprop(obj.model, 'dt')
                dt = obj.model.dt;
            else
                dt = obj.Ts;  % assume model.step covers full Ts if dt not provided
            end
            if abs(dt - obj.Ts) < 1e-9
                % Single step covers the whole period
                x_next = obj.model.step(x_current, u_current, w_current);
            else
                % Integrate in sub-steps to span Ts
                steps = round(obj.Ts / dt);
                if steps < 1
                    steps = 1;
                end
                for i = 1:steps
                    % Assume w_current is piecewise constant over Ts (use same w for each sub-step)
                    x_next = obj.model.step(x_next, u_current, w_current);
                end
            end
        end
        
        function J = computeStateSensitivity(obj, x, u, w)
            % COMPUTESTATESENSITIVITY: Numerically compute the Jacobian d(f)/d(x) for one interval.
            % (This is used for debug/analysis; in a real implementation, model's analytic sensitivity could be used.)
            Nx = length(x);
            fx = obj.propagateState(x, u, w);
            J = zeros(Nx, Nx);
            epsVal = 1e-6;
            for j = 1:Nx
                x_pert = x;
                x_pert(j) = x_pert(j) + epsVal;
                f_pert = obj.propagateState(x_pert, u, w);
                J(:, j) = (f_pert - fx) / epsVal;
            end
        end
        
        function y_pred_seq = predictOutputs(obj, x_seq)
            % PREDICTOUTPUTS: Compute predicted outputs for each state in a sequence.
            % If model has an output function H or C matrix, use it; otherwise assume full state is measured.
            y_pred_seq = cell(1, length(obj.yBuffer));
            for i = 1:length(obj.yBuffer)
                % We assume measurement at time corresponding to state x_seq(:, i+1) 
                % because yBuffer(1) corresponds to the output after first interval.
                % If yBuffer(1) is output at time k-N_e+1, that corresponds to state x_seq(:,2).
                % So index offset: yBuffer{i} corresponds to x_seq(:, i+1).
                state_index = i + 1;
                if state_index > size(x_seq, 2)
                    state_index = size(x_seq, 2); % guard (should not happen if buffers aligned)
                end
                x_i = x_seq(:, state_index);
                % Compute output from state (assuming linear output: y = C * x, or use model output function if defined)
                if isprop(obj.model, 'C')
                    y_pred = obj.model.C * x_i;
                elseif ismethod(obj.model, 'output')
                    y_pred = obj.model.output(x_i);
                else
                    y_pred = x_i;  % assume full state observable if no output mapping given
                end
                y_pred_seq{i} = y_pred;
            end
        end
    end % methods

end   % classdef