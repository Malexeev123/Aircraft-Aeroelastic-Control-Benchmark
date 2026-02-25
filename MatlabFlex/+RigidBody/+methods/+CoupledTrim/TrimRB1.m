function trim = trimFlexibleAircraft(cfg, beam, aero, base)
%TRIMFLEXIBLEAIRCRAFT  Solve for the coupled trim equilibrium of a flexible aircraft.
%  trim = trimFlexibleAircraft(cfg, beam, aero, base) finds the steady-state 
%  equilibrium (trim) of the aircraft with flexible wing. It returns a struct 
%  with trim results including the state vector for simulation.
%
% Inputs:
%   cfg   - configuration structure (contains flight conditions, etc.)
%   beam  - flexible wing structural model (BeamModel object with modal data)
%   aero  - aerodynamic model (AeroROM object or aerodynamic data for wing)
%   base  - coupling data structure (from build_xi_coupling, with matrices etc.)
%
% Outputs:
%   trim  - structure with trimmed state:
%           .q        = [Nm x 1] modal deflections at trim
%           .alpha    = [1 x 1] trimmed angle of attack (rad)
%           .delta_e  = [1 x 1] trimmed elevator deflection (rad)
%           .states   = full state vector for simulation initial condition
%           (Additional fields like .forces or .iterations can be included as needed)
%
% The function solves for: 
%   - q (modal elastic coordinates), 
%   - alpha (angle of attack), 
%   - delta_e (elevator deflection), 
%   - [optional] thrust (if included as unknown for force balance),
% that satisfy both rigid-body force/moment balance and flexible-wing modal force balance.
%
% The system of equations is solved using fsolve (Newton method) for accuracy.
% Convergence and the final wing shape are plotted for diagnostics.

    Nm = beam.Nm;          % number of flexible modes in structural model
    % Set initial guesses for unknowns:
    alpha0  = deg2rad(cfg.trim.alphaDeg);    % initial guess for angle of attack (rad)
    delta0  = deg2rad(cfg.trim.deltaDeg);    % initial guess for elevator deflection (rad)
    q0      = zeros(Nm,1);                   % initial guess for modal deflections (start with no deflection)
    X0 = [q0; alpha0; delta0];
    % If including thrust as unknown, append initial thrust guess:
    if isfield(cfg.trim, 'thrust') && ~isempty(cfg.trim.thrust)
        T0 = cfg.trim.thrust;  % initial thrust (possibly 0 for glider)
        X0 = [X0; T0];
        solveForThrust = true;
    else
        solveForThrust = false;
    end

    % Set up fsolve options for convergence tolerance and output display
    options = optimoptions('fsolve', 'Display', 'iter', ...
                           'FunctionTolerance', cfg.trim.tolNewton, ...
                           'StepTolerance',    1e-9, ...           % tighter step tol for accuracy
                           'OptimalityTolerance', cfg.trim.tolNewton, ...
                           'MaxIterations', 200, 'MaxFunctionEvaluations', 1000);
    % Use an OutputFcn to record residual norm at each iteration for plotting
    iterCount = 0;
    resHist = [];   % will record norm of residual at each iteration
    options.OutputFcn = @outFcn;
    
    % Nested output function to capture fsolve iteration data
    function stop = outFcn(x, optimValues, state)
        if strcmp(state, 'iter')
            iterCount = iterCount + 1;
            % Record the norm of the residual (or sum of squares)
            resnorm = norm(optimValues.fval); 
            resHist(iterCount,1) = resnorm;
        end
        stop = false;  % never stop early
    end

    % Define the residual function for trim (combines rigid-body and flexible dynamics)
    function R = trimResidual(X)
        % X contains [q(1..Nm); alpha; delta; (optional thrust)]
        q     = X(1:Nm);
        alpha = X(Nm+1);
        delta = X(Nm+2);
        if solveForThrust
            T = X(Nm+3);
        else
            T = cfg.trim.thrust;  % use given thrust from config if not solving for it
        end
        % Compute flexible wing forces: modal and clamping (root) forces
        [Fstruct, Faero, F_clamp] = evalFlexibleROM(q, alpha, aero, beam, base);
        % Fstruct = [Nm x 1] structural internal forces for each mode
        % Faero   = [Nm x 1] aerodynamic generalized forces on each mode
        % F_clamp = struct with fields F and M for wing root forces/moments on fuselage
        %
        % Compute rigid-body force/moment residuals (including gravity, thrust, tail, fin, wing)
        R_rb = evalRigidBodyEOM(F_clamp, alpha, delta, T, cfg);
        % R_rb is a 3x1 residual vector: [X-force; Z-force; pitch-moment] (all should be zero at trim)
        %
        % Compute flexible modal equilibrium residuals: internal - external = 0 for each mode
        R_flex = Fstruct - Faero;
        % (Each element i represents the net force in mode i; at trim this should be zero)
        %
        % Combine residuals into one vector for fsolve
        R = [R_flex; R_rb];
    end

    % Solve the nonlinear system for trim
    fprintf('[trimFlexibleAircraft] Solving for trim equilibrium...\n');
    X_sol = fsolve(@trimResidual, X0, options);
    % fsolve returns when residual is within tolerance or max iterations reached.

    % Extract solution
    q_sol     = X_sol(1:Nm);
    alpha_sol = X_sol(Nm+1);
    delta_sol = X_sol(Nm+2);
    if solveForThrust
        T_sol = X_sol(Nm+3);
    else
        T_sol = cfg.trim.thrust;
    end

    % Populate output structure
    trim.q       = q_sol;
    trim.alpha   = alpha_sol;
    trim.delta_e = delta_sol;
    trim.thrust  = T_sol;
    % Assemble full state vector for simulation initial condition
    % Note: state layout must match AeroFlex SimRunner expectations.
    % Typically: [flexible positions (Nm); flexible velocities (Nm); (possibly placeholder); aero states (Na); Euler angles (3)]
    Na = aero.Na;
    state = zeros(3*Nm + 1 + Na + 3, 1);  % pre-allocate state vector
    % Fill flexible modal states:
    state(1:Nm)         = q_sol;             % modal displacements
    state(Nm+1:2*Nm)    = 0;                % modal velocities (zero at trim)
    state(2*Nm+1:3*Nm)  = 0;                % (if a third block of modal states exists, e.g. in intrinsic formulation – set to zero)
    % Fill aerodynamic states (assumed zero initial perturbation from steady solution):
    state(3*Nm+2 : 3*Nm+1+Na) = zeros(Na,1);
    % Fill rigid-body states:
    % If we include one positional state (e.g., altitude or forward position) as indicated by '+1', set it:
    state(3*Nm+1) = 0;   % e.g., initial altitude (0) or placeholder if not used
    % Body velocity: Assuming body axes aligned with flight path initially
    % (these may be embedded in the aerodynamic states or not explicitly stored; here we assume constant velocity as no acceleration)
    % For completeness, one could include [u; v; w; p; q; r] if needed. In AeroFlex, these might be derived, so we leave them implicitly handled.
    % Euler angles:
    euler = [0; alpha_sol; 0];   % [roll; pitch; yaw], assuming small pitch = angle of attack, level wings, no sideslip
    state(end-2:end) = euler;
    trim.states = state;

    % --- Plot convergence diagnostics ---
    if ~isempty(resHist)
        figure('Name','Trim Convergence'); clf;
        semilogy(resHist, '-o'); grid on;
        xlabel('Iteration'); ylabel('Residual Norm');
        title('Trim Solution Convergence');
    end

    % --- Plot trimmed wing deformation ---
    if isfield(cfg.plotOpts, 'wingTip') && cfg.plotOpts.wingTip
        % Reconstruct approximate wing shape along the span using mode shapes.
        % Here we assume mode shapes are available or can be approximated.
        span = cfg.struct.L;  % half-span length [m] (from cfg.struct, e.g., L = wing half-span)
        y = linspace(0, span, 50)';    % spanwise coordinate from root (y=0) to tip (y=L)
        deflection = zeros(size(y));
        % If mode shapes (e.g., base.Phi_y) are available, use them. Otherwise, approximate:
        if isfield(base, 'Phi') 
            % base.Phi could be a matrix giving displacement shape of each mode along span (rows: modes, cols: positions)
            % (Assume base.Phi_z(:,j) might give mode shape in z-direction for node j; this is a placeholder.)
        end
        % Simple approximation: use first bending mode as sine shape, others ignored or add if known.
        if Nm >= 1
            % Assume first mode is a bending mode, approximate shape as a sine curve
            phi1 = sin(pi/2 * (y/span));  % clamped-free first mode shape (approximate)
            deflection = deflection + q_sol(1) * phi1;
        end
        if Nm >= 2
            % If second mode is a torsional mode, we could approximate twist (not plotted here).
            % (Could compute twist distribution and convert to tip deflection etc., omitted for simplicity.)
        end
        % Plot wing deflected shape (vertical deflection vs span)
        figure('Name','Trimmed Wing Shape'); clf;
        plot(y, deflection, 'b-', 'LineWidth', 2); hold on;
        plot([0 span], [0 0], 'k--');
        xlabel('Spanwise position (m)'); ylabel('Vertical deflection (m)');
        title('Trimmed Wing Bending Deflection');
        legend('Wing deflection','Undeformed');
        grid on;
    end

    fprintf('[trimFlexibleAircraft] Trim solution found: alpha = %.3f deg, delta_e = %.3f deg\n', ...
             rad2deg(alpha_sol), rad2deg(delta_sol));
    if solveForThrust
        fprintf('    (Required thrust: %.2f N)\n', T_sol);
    end
end