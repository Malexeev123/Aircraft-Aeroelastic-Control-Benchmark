function R_rb = EvalRigidBodyEOM(F_clamp, alpha, delta_e, T, cfg)
%EVALRIGIDBODYEOM  Compute rigid-body force and moment residuals for trim.
%  R_rb = evalRigidBodyEOM(F_clamp, alpha, delta_e, T, cfg) returns a 3x1 residual 
%  vector [R_X; R_Z; R_M] representing the net forces (X-forward, Z-down) and 
%  pitching moment about CG for the rigid-body dynamics. These are zero at trim.
%
% Inputs:
%   F_clamp - struct with wing root reaction forces and moments, typically:
%             F_clamp.F = [F_x; F_y; F_z] at wing root (body axes)
%             F_clamp.M = [M_x; M_y; M_z] moment at wing root (about root, body axes)
%   alpha   - angle of attack (rad)
%   delta_e - elevator deflection (rad) (positive = trailing-edge up)
%   T       - thrust force (N) (positive forward along body X-axis)
%   cfg     - configuration (contains mass, inertia, geom parameters, etc.)
%
% Outputs:
%   R_rb    - [3x1] residual: [X-force; Z-force; pitch-moment]
%             (should be [0;0;0] at equilibrium)

    % Unpack some config parameters
    rParams = RigidBody.methods.paramsRigid_PazyUAV();
    m   = rParams.mass;           % mass of aircraft (kg) - user should set cfg.mass
    % m   = rParams.mNoWing/1000;           % mass of aircraft (kg) - user should set cfg.mass
    g   = 9.81;               % gravitational acceleration (m/s^2)
    W   = m * g;              % weight (N)
    % If not provided in cfg, user must define aircraft geometry for moments:
    if isfield(cfg, 'geom')
        x_wing = cfg.geom.x_wing;   % [m] CG to wing aerodynamic center (positive if behind CG)
        x_tail = cfg.geom.x_tail;   % [m] CG to tail lift application (positive if behind CG)
    else
        % Placeholder: assume some geometry (user should replace with actual values)
        % x_wing = -0.1;  % wing AC 0.1 m *ahead* of CG (negative sign means behind CG is negative here)
        x_wing = 0.1;  % wing AC 0.1 m *ahead* of CG (negative sign means behind CG is negative here)
        % x_tail = 0.5;   % tail located 0.5 m behind CG (positive nose to tail)
        x_tail = -0.5;   % tail located 0.5 m behind CG (positive nose to tail)
    end

    % Compute tail and fin aerodynamic forces (body axes) given alpha & delta
    % aeroForces = aircraftAeroModel(alpha, delta_e, cfg);
    aeroForces = RigidBody.methods.rigidAeroDyn(alpha, 0, cfg);
    % aeroForces is a struct, e.g.:
    % aeroForces.F_tail = [Fxt; Fyt; Fzt] (N) tail forces
    % aeroForces.F_fin  = [Fxf; Fyf; Fzf] (N) fin forces
    % aeroForces.M_tail = [Mxt; Myt; Mzt] (N*m) tail moments about CG
    % aeroForces.M_fin  = [Mxf; Myf; Mzf] (N*m) fin moments about CG
    %
    % For symmetric trim, we expect Y forces ~0, and tail/fin moments mostly about Y (pitch) or Z (yaw).
    % Here we'll use the X,Z forces and Y (pitch) moment from these.

    % Wing aerodynamic forces: can take from wing root reaction (clamping) force.
    F_wing = F_clamp(1:3);   % [F_x; F_y; F_z] at wing root (body axes)
    M_wing_root = F_clamp(1:3);  % [M_x; M_y; M_z] moment about wing root

    % Sum forces in X direction (forward)
    % Positive X-force residual means excess thrust or insufficient drag.
    % We include wing, tail, fin drag (forward forces negative).
    % Drag components: negative of X-components of aero forces (since X forward).
    % Note: F_wing(1) includes wing drag (likely negative or zero if lift only in Z), 
    % tail/fin X-forces likely negative (drag).
    Drag_total = -F_wing(1) ...              % wing drag (wing F_x is typically negative if drag)
                 - aeroForces.F_tail(1) ...  % tail drag
                 - aeroForces.F_fin(1);      % fin drag
    R_X = T - Drag_total;
    % Sum forces in Z direction (downwards)
    % Positive Z is downward. Weight = +W (down). Lift (up) appears as negative force in +Z axis.
    % Wing + tail lift: wing F_z + tail F_z (both likely negative if producing upward lift, 
    % tail F_z could be positive if tail is pushing down).
    Lift_total = -F_wing(3) - aeroForces.F_tail(3);
    % (We take upward lift as positive quantity: wing upward lift = -F_wing(3) since F_z pos is down)
    % Tail force: if elevator up, tail produces downforce -> F_tail(3) positive, then -F_tail(3) is negative (negative "lift").
    % Lift_total effectively = L_wing + L_tail (taking L_tail negative if it's a downforce).
    % R_Z = W - Lift_total;
    R_Z = W + Lift_total;
    % R_Z = -(W - Lift_total);
    % Now moments about CG (pitching moment, Y-axis).
    % We need wing moment about CG: combine wing root moment and the moment from wing net lift about CG.
    % Convert wing root forces into moments about CG:
    % If wing root is offset from CG (e.g., wing AC at x_wing behind CG):
    % wing lift (F_wing(3) upward negative) at x_wing produces a moment about CG: M_wing_CG = -x_wing * F_wing(3).
    % (Using right-hand rule: if x_wing is negative (AC ahead of CG), an upward wing lift (negative F_z) yields nose-down moment -> negative M_y.)
    M_wing_CG = - x_wing * F_wing(3);
    % Include wing root pitching moment (M_y component). Note: wing root moment M_y is nose-up positive as given.
    M_wing_total = M_wing_root(2) + M_wing_CG;
    % Tail pitching moment about CG: tail force (mostly vertical) at distance x_tail:
    % If tail produces lift up (negative F_tail(3)), at behind CG (x_tail > 0) -> nose-up or nose-down?
    % Tail downforce (F_tail(3) positive) * behind CG gives nose-up moment.
    M_tail =  x_tail * (-aeroForces.F_tail(3));
    % (-F_tail(3) is upward lift on tail; if tail force is downward (positive F_tail(3)), this yields -F_tail(3) negative, times +x_tail = negative moment (nose-down), 
    % which is counter-intuitive because a tail downforce actually gives nose-up. 
    % Let's double-check: If tail downforce = positive F_tail(3), 
    % torque = r x F, r (back) positive, F (down) positive -> nose-up (positive M_y).
    % Our sign above: M_tail = x_tail * (-F_tail(3)) yields for F_tail(3)>0 (downforce): -F_tail(3) <0 times +x_tail = negative -> nose-down, which is opposite.
    % We need to correct the sign: instead, define positive delta_e yields negative lift (downforce) and we want that to give positive nose-up moment.
    % So let's do: M_tail = - x_tail * F_tail(3).)
    M_tail = - x_tail * F_wing(3);  % *** This line is incorrect logically; using aeroForces below instead.
    M_tail =  x_tail * aeroForces.F_tail(3);  % Revised: if F_tail(3) is positive (downwards), moment = positive * positive = positive nose-up.
    % Actually, correction: For clarity, compute tail moment from forces using cross product concept:
    % Nose-up moment (about Y) from tail = (tail lift * lever arm) with appropriate sign.
    % We expect: if tail force is downward (positive F_tail(3)), at x_tail behind CG, this produces nose-up -> positive M_y.
    % Thus: M_tail = + x_tail * F_tail(3). (Because F_tail(3) positive yields positive moment).
    M_tail =  x_tail * F_clamp;  % << placeholder, see above calculation >>
    %   ^ Correction needed: the above is incomplete in code context. We have to implement the logic described:
    if aeroForces.F_tail(3) >= 0 
        % tail downforce
        M_tail =  x_tail * aeroForces.F_tail(3);   % positive nose-up
    else
        % tail upward lift
        M_tail =  x_tail * aeroForces.F_tail(3);   % this will be negative (nose-down)
    end
    % Sum wing and tail contributions for total pitching moment (M_y):
    M_total = M_wing_total + M_tail + aeroForces.M_tail(2) + aeroForces.M_fin(2);
    M_total = M_wing_total  + aeroForces.M_tail(2) + aeroForces.M_fin(2);
    % (We also add any pitching moment from tailplane aerodynamic moment and fin (though fin likely contributes little to pitch). 
    % aeroForces.M_tail(2) could include tail's own aerodynamic pitching moment about its quarter-chord, etc. Fin M_y might be small.)
    R_M = M_total;
    % The residuals are defined such that they should be zero at equilibrium:
    R_rb = [R_X; R_Z; R_M];
end