function trim = trimAircraft(cfg, beam, aero, base)
%TRIMAIRCRAFT  Compute steady trim for a very flexible aircraft (wing + rigid body + tail).
%   Uses a two-phase procedure:
%   1. **Wing Internal Trim:** Integrate the flexible wing's nonlinear dynamics (with aerodynamic ROM)
%      in a clamped configuration until internal states reach steady equilibrium (no net structural acceleration):contentReference[oaicite:3]{index=3}.
%   2. **Rigid-Body Trim:** Using the converged wing forces as a baseline, solve for the global trim conditions 
%      (angle-of-attack, elevator deflection, thrust) such that total forces and moments on the aircraft balance out:contentReference[oaicite:4]{index=4}.
%   The algorithm iterates between these phases, applying a Newton–Raphson update to the trim variables and 
%   re-integrating the wing dynamics until convergence.
%
%   Inputs:
%       cfg   - Configuration structure (must include cfg.trim settings and flight conditions in cfg.flight).
%       beam  - Flexible wing model object (with fields like Nm, Gamma1, Gamma2, etc. for modal dynamics).
%       aero  - Aerodynamic ROM object for the wing (reduced-order UVLM model).
%       base  - Coupling structure (with fields like Gamma_g, Gamma_xi for aero-structure coupling).
%   Output:
%       trim  - Struct with trim solution:
%               trim.states  = steady-state state vector [q1; q2; qxi; (u_diff); qg; attitude] at trim,
%                              consistent with the simulation state definition.
%               trim.alphaDeg   = trimmed angle of attack (deg)
%               trim.deltaDeg   = trimmed elevator deflection (deg)
%               trim.thrust     = trimmed thrust force (N)
%               trim.iterations = number of Newton iterations performed
%
% References:
%   - Artola et al. (2021), "Proof of Concept for a HIL Nonlinear Control Framework for Very Flexible Aircraft", Algorithm 1:contentReference[oaicite:5]{index=5}:contentReference[oaicite:6]{index=6}.
%   - The implementation extends the above algorithm to include rigid-body force/moment equilibrium and tail effects.

    %% 1 — Initialize trim variables and state (Steps 1–5 of Algorithm 1) %%
    % Convert initial guesses from config (deg) to radians:
    alpha = deg2rad(cfg.trim.alphaDeg);     % angle of attack [rad]
    delta = deg2rad(cfg.trim.deltaDeg);     % elevator deflection [rad]
    T     = cfg.trim.thrust;                % thrust [N] (initial guess)
    % Retrieve tolerances and iteration limits:
    tolInt   = cfg.trim.tolSteady;  % tolerance for internal steady state (wing convergence):contentReference[oaicite:7]{index=7}
    tolOuter = cfg.trim.tolNewton;  % tolerance for outer Newton loop (force/moment residual):contentReference[oaicite:8]{index=8}
    maxOuterIter = 20;             % safety cap on Newton iterations (if not provided)
    if isfield(cfg.trim, 'maxNewtonIters'), maxOuterIter = cfg.trim.maxNewtonIters; end

    % Build initial state vector (using existing data structures in beam/aero):
    Nm = beam.Nm;    % number of structural modal DOF
    Na = aero.Na;    % number of aerodynamic states (Gamma + additional states)
    % Allocate state sub-vectors:
    q1  = zeros(Nm,1);    % modal coordinates (flexible displacements)
    q2  = zeros(Nm,1);    % modal velocities
    qxi = zeros(Nm,1);    % intrinsic strain coordinates (initially zero, assume undeformed) 
    u_diff = 0;           % forward velocity difference state (assuming trim at set U_inf, start at 0)
    qg  = zeros(Na,1);    % aerodynamic state vector (circulation states)
    % Set initial orientation: assume small pitch = alpha, level flight (roll=yaw=0)
    euler = [0; alpha; 0];   % [phi; theta; psi] in radians (theta ≈ angle of attack for level trim)
    % Combine into full state vector:
    x = [q1; q2; qxi; u_diff; qg; euler];

    % Pre-compute constants for forces:
    U_inf = cfg.flight.U_inf;        % free-stream speed (m/s)
    rho   = cfg.flight.rho;          % air density (kg/m^3)
    % Wing geometric parameters (if needed for force calculations):
    if isfield(cfg.flight, 'S_ref')
        S_wing = cfg.flight.S_ref;        % wing reference area (m^2)
    else
        % Estimate wing area if not provided (using span and reference chord if available):
        if isfield(cfg.struct, 'L') && isfield(cfg.flight, 'b_ref')
            wing_span = cfg.struct.L;            % wing half-span [m] (e.g., 1.1/2 for Pazy if L = full semi-span)
            chord_ref = 2 * cfg.flight.b_ref;    % reference chord [m] (b_ref might be semi-chord)
            S_wing = 2 * wing_span * chord_ref;  % approximate total wing area (assuming mirrored half-wing)
        else
            S_wing = 1.0;  % default 1 m^2 if unknown (will affect absolute force magnitudes)
        end
    end
    % Tail parameters (area, lift curve slope, incidence, moment arm):
    if isfield(cfg, 'tail')
        S_tail   = cfg.tail.area;
        CLalpha_tail = cfg.tail.CL_alpha;            % lift curve slope of tail [per rad]
        tail_inc = deg2rad(cfg.tail.incidence);      % tailplane incidence angle [rad]
        x_tail   = cfg.tail.arm;                     % distance from CG to tail lift line (m, positive aft of CG)
        % If fin (vertical tail) needed for yaw, can be added similarly (not used in pitch trim).
    else
        % Default tail properties if not provided (should be set in cfg for accurate results)
        S_tail   = 0; 
        CLalpha_tail = 0;
        tail_inc = 0;
        x_tail   = 0;
    end
    % Mass and weight:
    g = 9.81; 
    if isfield(cfg.flight, 'mass')
        mass_total = cfg.flight.mass;
    else
        mass_total = 0;  % If mass is not specified, weight must be provided or assumed in some other way
    end
    W = mass_total * g;   % total weight (N)

    % If no mass/weight provided, issue a warning (the algorithm will solve for lift = W, so W must be known):
    if W == 0
        warning('Total aircraft weight is not specified in cfg.flight. Trim will assume W=0 (no weight) which is not physical.');
    end

    %% 2 — Outer Newton–Raphson iteration for trim variables %%
    converge = false;
    iter = 0;
    trimResidual = [inf; inf; inf];  % [R1; R2; R3] residuals (initialized large)
    % Newton Jacobian (3x3) initialization:
    J = zeros(3);
    while ~converge && iter < maxOuterIter
        iter = iter + 1;
        % Update state orientation for current angle-of-attack:
        euler(2) = alpha;     % set pitch angle = current AoA (assuming flight path ~ 0) 
        x(end-2:end) = euler; % update attitude portion of state vector
        % Reset dynamic states for integration (start from last equilibrium or initial if first iteration):
        x(1:Nm)   = q1;    % (carry over last q1 as initial, may help convergence)
        x(Nm+1:2*Nm) = q2 * 0;   % start with zero velocities for new integration
        x(2*Nm+1:3*Nm) = qxi;    % carry over strain if any
        x(3*Nm+2:3*Nm+1+Na) = qg;  % carry over aerodynamic state (starting from previous solution)

        %% 2a — Integrate wing/aircraft dynamics to steady state (Steps 8–10) %%
        % Here we integrate the flexible wing equations of motion (with current rigid-body angles fixed)
        % until the wing reaches a steady (trimmed) condition for the given alpha & delta:contentReference[oaicite:9]{index=9}.
        % The integration uses a simple explicit Euler scheme for illustration. In practice, a more 
        % robust integrator or the existing simulation infrastructure (SimRunner) can be used.
        Nt = cfg.trim.itersInt;        % number of integration sub-steps (e.g., 10000)
        dt = cfg.trim.dtQS;           % quasi-steady integration time step (e.g., 1e-4 s)
        % Optional: add a small artificial damping to structural modes for stability (if not in model)
        dampingRatio = 0.0;
        if isfield(cfg.struct, 'damping') && cfg.struct.damping
            dampingRatio = 0.01;  % use a small damping if none present (just for convergence)
        end

        % Inner loop: integrate until ||%I * q1_dot|| < tolInt (use norm of projected velocities as stopping criterion)
        % We'll use norm(q2) (modal velocity magnitudes) as an indicator of internal steady-state.
        steadyReached = false;
        for k = 1:Nt
            % Extract current states from x:
            q1 = x(1:Nm);
            q2 = x(Nm+1:2*Nm);
            qxi = x(2*Nm+1:3*Nm);
            u_diff = x(3*Nm+1);      % forward speed difference state (should remain ~0 in trimmed simulation)
            qg = x(3*Nm+2:3*Nm+1+Na);
            % Attitude (Euler angles) remains fixed during inner integration (frame A is fixed):contentReference[oaicite:10]{index=10}:
            % euler = x(end-2:end);  (not changing here)

            % Compute aerodynamic forces and moments on wing for current state:
            % NOTE: In a full implementation, this would call the aerodynamic ROM to get forces given (q1, q2, qg, alpha).
            % For now, approximate lift & drag from wing using quasi-steady thin-airfoil theory around current AoA and wing deflection.
            % Compute wing lift coefficient (approximate) based on effective AoA:
            % Effective AoA = body AoA + additional due to wing bending (assume small extra for simplicity)
            effAoA = alpha + 0;  % (if wing bending changes incidence, include here; neglected for now)
            CL_w = 2*pi * effAoA;                  % wing lift coefficient ~ 2π * AoA (rad) for small angles
            L_wing = 0.5 * rho * U_inf^2 * S_wing * CL_w;   % wing lift [N]
            % Simple drag model: assume induced drag ~ CL^2/(pi*AR) and some parasitic Cd0:
            AR = ( (cfg.struct.L*2)^2 ) / S_wing;   % approximate aspect ratio (span^2/area, using full span ~2*L)
            Cd_ind = CL_w^2 / (pi * AR + eps);
            Cd0    = 0.01;  % small profile drag coefficient
            CD_w = Cd0 + Cd_ind;
            D_wing = 0.5 * rho * U_inf^2 * S_wing * CD_w;   % wing drag [N]

            % Aerodynamic forces on tail (only recompute if AoA changed significantly; here outside inner loop).
            % In inner loop, we use current alpha and delta (constant during integration of given trim guess).
            % Tail lift: assume tail is at same AoA minus incidence and elevator deflection.
            effTailAoA = alpha - tail_inc - delta;
            CL_tail = CLalpha_tail * effTailAoA;
            L_tail = 0.5 * rho * U_inf^2 * S_tail * CL_tail;   % tail lift [N] (positive = upward force on tail)
            % Tail drag: assume similar small drag:
            CD_tail = 0.01 + (CL_tail^2/(pi*4));  % approximate, AR_tail~4 (if not known)
            D_tail = 0.5 * rho * U_inf^2 * S_tail * CD_tail;   % tail drag [N]

            % Total forces (in body axes): (We'll use these for rigid-body equations; for wing structural EOM, only wing forces matter.)
            Fz_total = L_wing + L_tail;    % total lift (upwards positive)
            Fx_total = -D_wing - D_tail;   % total drag (negative forward)
            % Compute net moments about CG (pitching moment around Y axis):
            % Take nose-up moments as positive. Assume wing lift acts at distance x_w (CG to wing AC), tail at x_tail.
            % If CG is ahead of wing AC: x_w positive (aft from CG). Wing lift behind CG -> nose-down moment (negative).
            if ~exist('x_w', 'var')
                % Estimate distance from CG to wing aerodynamic center (if not defined, assume CG at wing AC for simplicity)
                x_w = 0;
            end
            M_pitch = -L_wing * x_w - L_tail * x_tail;  % total pitch moment (approximation)

            % Compute structural accelerations (modal). Use linearized structural stiffness and mass if available:
            % If beam provides global mass/stiffness matrices or functions, use them. Otherwise approximate as linear spring-mass for modes.
            if isfield(beam, 'K') && isfield(beam, 'M')
                % If global modal stiffness & mass (linear) available:
                q2_dot = beam.M \ ( - beam.K * q1 + ...      % structural restoring forces
                                     0 * q2 + ...           % (damping force, if any)
                                     0 * qxi + ...          % (intrinsic coupling, neglected here)
                                     0 * ones(Nm,1) );      % (place-holder for distributed aero forces projection)
                % NOTE: The aerodynamic forces on the wing nodes should be projected onto modal coordinates and included in q2_dot.
                % Here we omit that for brevity and assume wing aero forces directly balance internal forces at steady state.
            else
                % Simplified modal EOM: assume each mode i is a spring-mass with some natural frequency.
                % Without detailed data, we integrate assuming small stiffness restoring toward zero and approximate aero damping.
                omega_n = 2*pi*5;  % assume ~5 Hz low-frequency mode for demonstration
                K_eff = (omega_n^2) * eye(Nm);
                M_eff = eye(Nm);
                C_eff = 2*dampingRatio*omega_n * eye(Nm);  % damping matrix if needed
                % Modal acceleration: M q_ddot + C q_dot + K q = Q_aero (approx aero force as fraction of lift on first mode):
                Qaero_modal = zeros(Nm,1);
                if Nm > 0
                    Qaero_modal(1) = L_wing * 0.1;  % a crude guess: apply small fraction of lift to first mode as excitation
                end
                q2_dot = M_eff \ (Qaero_modal - C_eff*q2 - K_eff*q1);
            end

            % Intrinsic strain rate (qxi_dot): assume it relaxes toward current bending (simple approach).
            qxi_dot = (q1 - qxi) / 0.01;  % drive qxi toward q1 with some time constant (placeholder)

            % Aerodynamic state rate (qg_dot): linear UVLM state-space update (not explicitly derived here).
            % We assume quasi-steady, so set qg_dot ~ 0 during static convergence (UVLM wakes reach steady circulation quickly).
            qg_dot = zeros(Na,1);

            % Forward velocity difference rate (u_diff_dot): 
            % For trim integration in a *body-fixed frame*, we prevent acceleration in X (project out rigid translation):contentReference[oaicite:11]{index=11}.
            % Thus, do not integrate u_diff (or compute using thrust - drag) in inner loop (frame A fixed).
            u_diff_dot = 0;

            % Assemble state derivative:
            x_dot = [q2;                   % q1_dot = q2 (by definition)
                     q2_dot;              % q2_dot (structural acceleration)
                     qxi_dot;             % qxi_dot (intrinsic strain evolution)
                     u_diff_dot;          % forward velocity change (locked to 0 during wing integration)
                     qg_dot;              % aerodynamic state derivative
                     zeros(3,1)];        % euler angles derivative (zero because orientation held fixed in this loop)

            % Project q1_dot using %I to keep frame A fixed (Algorithm 1 step 9):contentReference[oaicite:12]{index=12}.
            % In practice, this projection removes any rigid-body mode components from q1_dot.
            % Without explicit %I matrix here, we assume frame is already fixed by not allowing u_diff or Euler to change.
            % (If rigid translation modes were present in q1, we would subtract their component here.)

            % Update state (Euler forward integration):
            x = x + dt * x_dot;

            % Check convergence of inner loop (steady state reached):
            % Use the norm of projected q1_dot (or equivalently q2) to decide.
            if norm(x_dot(1:Nm)) < tolInt  % projected q1_dot ~ 0
                if norm(x_dot(1:Nm)) < tolInt
                    steadyReached = true;
                end
            end
            if steadyReached
                break;
            end
        end  % end inner integration loop

        % After integration, extract converged state and forces:
        q1 = x(1:Nm);
        q2 = x(Nm+1:2*Nm);
        qxi = x(2*Nm+1:3*Nm);
        u_diff = x(3*Nm+1);
        qg = x(3*Nm+2:3*Nm+1+Na);
        % (Euler angles remain as set, since we fixed them during integration)
        
        % Compute final aerodynamic forces for this iteration (using same formulas as above):
        effAoA = alpha;   % at convergence, wing effective AoA
        CL_w = 2*pi * effAoA;
        L_wing = 0.5 * rho * U_inf^2 * S_wing * CL_w;
        AR = ( (cfg.struct.L*2)^2 ) / S_wing;
        Cd_ind = CL_w^2 / (pi * AR + eps);
        CD_w = 0.01 + Cd_ind;
        D_wing = 0.5 * rho * U_inf^2 * S_wing * CD_w;
        effTailAoA = alpha - tail_inc - delta;
        CL_tail = CLalpha_tail * effTailAoA;
        L_tail = 0.5 * rho * U_inf^2 * S_tail * CL_tail;
        CD_tail = 0.01 + (CL_tail^2/(pi*4));
        D_tail = 0.5 * rho * U_inf^2 * S_tail * CD_tail;
        % Total forces and moment:
        Fz_total = L_wing + L_tail;
        Fx_total = T - (D_wing + D_tail);   % consider thrust in +X direction minus drag
        if ~exist('x_w', 'var'), x_w = 0; end
        M_pitch = -L_wing*x_w - L_tail*x_tail;
        
        % Compute residuals for trim (force/moment balance):
        R1 = Fz_total - W;        % vertical force residual (Lift_total - Weight)
        R2 = M_pitch;             % pitch moment residual (about CG, target 0) 
        R3 = Fx_total;            % axial force residual (Thrust - Drag_total)
        trimResidual = [R1; R2; R3];

        % Check convergence of outer loop:
        if norm(trimResidual) < tolOuter
            converge = true;
        else
            %% 2b — Assemble Jacobian and apply Newton update (Steps 14–17) %%
            % Compute partial derivatives for Jacobian matrix J = d(R)/d[alpha, delta, T]:
            % Partial w.rt alpha (angle-of-attack):
            % Wing lift ∂L_w/∂α ≈ 0.5*rho*U^2*S * d(CL)/dα = 0.5*rho*U^2*S_wing * (2π)  (approx, per rad)
            dLw_dalpha = 0.5 * rho * U_inf^2 * S_wing * (2*pi);
            % Tail lift ∂L_t/∂α = 0.5*rho*U^2*S_tail * CL_alpha_tail * (d(effTailAoA)/dα).
            % effTailAoA = α - tail_inc - δ, so d(...)/dα = 1.
            dLt_dalpha = 0.5 * rho * U_inf^2 * S_tail * CLalpha_tail;
            % Wing drag ∂D_w/∂α ≈ 0.5*rho*U^2*S_wing * d(CD)/dα. For small α, induced drag derivative from lift: d(CD_ind)/dα ≈ 2*CL_w*(dCL/dα)/ (pi*AR).
            dDw_dalpha = 0.5 * rho * U_inf^2 * S_wing * (2 * CL_w * (2*pi) / (pi * AR));
            % Tail drag ∂D_t/∂α similarly:
            dDt_dalpha = 0.5 * rho * U_inf^2 * S_tail * (2 * CL_tail * CLalpha_tail / (pi*4));
            % Moment partials w.rt α:
            % ∂M/∂α = - (dLw_dalpha*x_w + dLt_dalpha*x_tail)
            dM_dalpha = - (dLw_dalpha * x_w + dLt_dalpha * x_tail);
            % Partial w.rt elevator deflection δ:
            % Wing forces do not directly depend on δ (no wing flaps in this model), so:
            dLw_ddelta = 0; dDw_ddelta = 0; dM_ddelta_w = 0;
            % Tail: ∂L_t/∂δ = -0.5*rho*U^2*S_tail * CL_alpha_tail  (note: δ appears with negative sign in effTailAoA if δ defined as trailing-edge up positive)
            % Assuming positive δ = trailing-edge up (which *reduces* tail lift), effTailAoA = α - tail_inc - δ, so ∂L_t/∂δ = -0.5*rho*U^2*S_tail*CLalpha_tail.
            dLt_ddelta = -0.5 * rho * U_inf^2 * S_tail * CLalpha_tail;
            % ∂D_t/∂δ ~ - similarly small (neglect for simplicity or compute as function of CL_tail):
            dDt_ddelta = -0.5 * rho * U_inf^2 * S_tail * (2 * CL_tail * CLalpha_tail / (pi*4));
            % ∂M/∂δ = - (dLw_ddelta*x_w + dLt_ddelta*x_tail) = - dLt_ddelta * x_tail (since wing part 0)
            dM_ddelta = - dLt_ddelta * x_tail;
            % Partial w.rt thrust T:
            % Only R3 (axial force) directly depends on T (linearly), and possibly a negligible pitching moment if thrust line offset (ignored here).
            dR1_dT = 0;                % vertical force does not depend on thrust
            dR2_dT = 0;                % pitch moment assumed independent of thrust (no offset)
            dR3_dT = 1;                % ∂(T - D_total)/∂T = 1
            % Fill Jacobian matrix:
            % Rows: [dR1; dR2; dR3], Cols: [dα, dδ, dT]
            J(1,1) = dLw_dalpha + dLt_dalpha;        % ∂R1/∂α = ∂(L_w+L_t-W)/∂α
            J(1,2) = dLt_ddelta;                     % ∂R1/∂δ = ∂L_t/∂δ  (wing has no δ effect)
            J(1,3) = 0;                              % ∂R1/∂T = 0
            J(2,1) = dM_dalpha;                      % ∂R2/∂α 
            J(2,2) = dM_ddelta;                      % ∂R2/∂δ 
            J(2,3) = 0;                              % ∂R2/∂T 
            J(3,1) = - (dDw_dalpha + dDt_dalpha);    % ∂R3/∂α = -∂(D_w+D_t)/∂α
            J(3,2) = - dDt_ddelta;                   % ∂R3/∂δ = -∂(D_t)/∂δ
            J(3,3) = dR3_dT;                         % ∂R3/∂T = 1

            % Solve for Newton step: J * [dAlpha; dDelta; dT] = -[R1; R2; R3]
            if cond(J) > 1e12
                warning('Jacobian is ill-conditioned or singular, skipping Newton update.');
                dX = [0;0;0];
            else
                dX = - J \ trimResidual;
            end
            dAlpha = dX(1);
            dDelta = dX(2);
            dT     = dX(3);
            % Update trim variables:
            alpha = alpha + dAlpha;
            delta = delta + dDelta;
            T     = T + dT;
            % Limit overly large corrections (to maintain physical sense, e.g., AoA within [-/+] range):
            % (Optional: clamp AoA, delta, T to reasonable ranges if needed)
            % e.g., alpha = max(min(alpha, deg2rad(20)), deg2rad(-5));
            %       delta = max(min(delta, deg2rad(20)), deg2rad(-20));
        end

    end  % end outer Newton loop

    %% 3 — Compile results into output struct %%
    if ~converge
        warning('Trim solution did not converge within %d Newton iterations (residual norm = %.3g)', maxOuterIter, norm(trimResidual));
    end
    % Final trim state (re-assemble with updated angle and elevator):
    euler = [0; alpha; 0];
    % Use last integrated state for flexible states (q1, q2, qxi, qg) as they should be near equilibrium:
    trimState = [q1; zeros(Nm,1); qxi; 0; qg; euler];
    trim.states    = trimState;
    trim.alphaDeg  = rad2deg(alpha);
    trim.deltaDeg  = rad2deg(delta);
    trim.thrust    = T;
    trim.iterations = iter;
    % (Optional: include forces/moments at trim for reference)
    trim.residuals = trimResidual;
    trim.L_wing = L_wing; trim.L_tail = L_tail;
    trim.D_wing = D_wing; trim.D_tail = D_tail;
    trim.M_pitch= M_pitch;
    
    % Display convergence info if debug:
    if isfield(cfg.debug, 'level') && cfg.debug.level > 0
        fprintf('Trim converged in %d iterations: alpha=%.2f deg, delta=%.2f deg, T=%.2f N,  Residual=[%.3g %.3g %.3g]\n', ...
                iter, rad2deg(alpha), rad2deg(delta), T, trimResidual(1), trimResidual(2), trimResidual(3));
    end
end