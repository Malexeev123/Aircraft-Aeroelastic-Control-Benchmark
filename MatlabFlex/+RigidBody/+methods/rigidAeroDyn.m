function aeroForces = rigidAeroDyn(alpha, delta_e, cfg)
%AIRCRAFTAEROMODEL  Compute steady aerodynamic forces for tail and fin.
%  aeroForces = aircraftAeroModel(alpha, delta_e, cfg) returns a struct with tail 
%  and fin forces/moments in body axes for given angle of attack and elevator deflection.
%
% Inputs:
%   alpha   - angle of attack (rad)
%   delta_e - elevator deflection (rad) (positive = trailing-edge up)
%   cfg     - config structure containing aerodynamic coefficients and geometry:
%             cfg.flight.rho, cfg.flight.U_inf (air density, speed)
%             cfg.tail.S, cfg.tail.CL_alpha, cfg.tail.CL_delta, cfg.tail.CD0 (tail area, lift slopes, drag)
%             cfg.tail.x (longitudinal distance from CG, m, positive aft)
%             cfg.fin.S, cfg.fin.CD0 (fin area, drag coeff)
%             cfg.fin.x, cfg.fin.z (position of fin for moment calc if needed)
%
% Output:
%   aeroForces - struct with fields:
%                F_tail = [Fxt; Fyt; Fzt] tail force in body axes (N)
%                M_tail = [Mxt; Myt; Mzt] tail moment about CG (N*m)
%                F_fin  = [Fxf; Fyf; Fzf] fin force in body axes (N)
%                M_fin  = [Mxf; Myf; Mzf] fin moment about CG (N*m)
%
% Notes: Positive Z axis is downwards. Tail lift upward will appear as negative Fzt.

    % Dynamic pressure
    rho = cfg.flight.rho;
    V   = cfg.flight.U_inf;
    q_inf = 0.5 * rho * V^2;
    % Tail parameters (with some defaults if not provided)
    if ~isfield(cfg, 'tail')
        % Default tail parameters (user should specify these in cfg for accuracy)
        % cfg.tail.S         = 0.2 * cfg.flight.S_ref;   % tail area (assume 20% of wing area)
        cfg.tail.S         = 0.45 * cfg.flight.S_ref;   % tail area (assume 20% of wing area)
        cfg.tail.CL_alpha  = 4*pi;                    % lift curve slope (per rad)
        % cfg.tail.CL_alpha  = 2*pi;                    % lift curve slope (per rad)
        cfg.tail.CL_delta  = 2*pi * 0.8;              % elevator effectiveness (0.8 factor)
        cfg.tail.CD0       = 0.01;                    % zero-lift drag
        cfg.tail.e         = 0.8;                     % Oswald efficiency for tail
        % cfg.tail.x         =  +0.5;                   % tail 0.5 m behind CG
        cfg.tail.x         =  -0.5;                   % tail 0.5 m behind CG
        cfg.tail.Z_offset  =  0.0;                    % tail on body axis (no vertical offset)
    end
    if ~isfield(cfg, 'fin')
        % Default fin parameters 
        cfg.fin.S   = 0.1 * cfg.flight.S_ref;  % fin area (10% of wing area)
        cfg.fin.CD0 = 0.02;                   % fin drag coeff
        % cfg.fin.x   =  0.5;                   % fin behind CG
        cfg.fin.x   =  -0.5;                   % fin behind CG
        cfg.fin.Z_offset = 0.2;              % fin height above CG (for yaw moment if needed)
    end

    % Effective angle at tail (considering downwash – here we assume none or include a simple factor):
    alpha_tail = alpha * 0.9;  % assume tail sees 90% of fuselage AoA due to wing downwash (simple estimate)
    % Compute tail lift coefficient
    CL_tail = cfg.tail.CL_alpha * alpha_tail + cfg.tail.CL_delta * ( - delta_e );
    % (Note: minus sign for delta_e because positive delta_e (TE up) reduces lift.)
    % Tail lift force (positive upward). In body axes (X forward, Z down), upward lift = -F_z.
    % L_tail = (cfg.flight.b_ref/2)^2*q_inf * cfg.tail.S * CL_tail;
    L_tail = cfg.tail.S * CL_tail;
    % L_tail = (cfg.flight.b_ref)^2*q_inf * cfg.tail.S * CL_tail;
    % L_tail = q_inf * cfg.tail.S * CL_tail;
    % Tail drag coefficient (assume quadratic or use zero-lift + induced)
    CD_tail = cfg.tail.CD0 + (CL_tail^2)/(pi * cfg.tail.e * 4);  % simple induced drag, AR~4 (placeholder)
    % D_tail = (cfg.flight.b_ref/2)^2*q_inf * cfg.tail.S * CD_tail;
    % D_tail = (cfg.flight.b_ref/2)^2*q_inf * cfg.tail.S * CD_tail;
    % D_tail = (cfg.flight.b_ref)^2*q_inf * cfg.tail.S * CD_tail;
    D_tail = q_inf * cfg.tail.S * CD_tail;
    % D_tail =   cfg.tail.S * CD_tail;
    % Tail forces in body axes:
    aeroForces.F_tail = [ -D_tail;  0;  -L_tail ];
    % aeroForces.F_tail = [ -D_tail;  0;  L_tail ];
    % Tail moment about CG:
    % Tail located cfg.tail.x behind CG (positive x means behind). Tail also could be at some vertical offset (if not on same level as CG).
    % Pitching moment from tail lift: M_y = x_tail * L_tail (nose-up positive if L_tail is upward and tail behind CG).
    M_tail_y = cfg.tail.x * L_tail; 
    if L_tail < 0
        % If L_tail is negative (downforce), plug in formula carefully: downforce (negative L_tail) at tail behind CG yields nose-down (-) moment.
        M_tail_y = cfg.tail.x * L_tail;
    end
    aeroForces.M_tail = [0; M_tail_y; 0];
    % (We ignore any vertical offset of tail for pitch moment; if tail is above/below CG, it would also create a pitching moment due to lift being off CG height.)
    % No roll or yaw moments from tail in symmetric flight.

    % Fin forces:
    % At zero sideslip (beta=0), fin has no lift. Only drag.
    CD_fin = cfg.fin.CD0;
    % D_fin = (cfg.flight.b_ref)^2*q_inf * cfg.fin.S * CD_fin;
    D_fin =  cfg.fin.S * CD_fin;
    % D_fin = (cfg.flight.b_ref/2)^2*q_inf * cfg.fin.S * CD_fin;
    % D_fin = q_inf * cfg.fin.S * CD_fin;
    aeroForces.F_fin = [ -D_fin; 0; 0 ];
    % Fin moment about CG:
    % Fin drag acting above CG could create a pitching moment (nose-down if fin above CG). We'll include if Z_offset:
    M_fin_y = 0;
    if isfield(cfg.fin, 'Z_offset')
        % Drag force (backwards) at some vertical offset -> creates a pitch moment: M_y = D_fin * Z_offset (nose-down if fin above CG, which is positive M_y by RH rule? Actually drag backward above CG tends to pitch nose-up? Let's reason:
        % Drag at fin above CG tries to pull the top backward, pitching nose up (think of force above center pulling back -> nose goes up).
        % If Z_offset positive (fin above CG), drag (X-) -> moment around Y = -Z_offset * Drag (negative means nose-down? Need cross product: r (0,0,Z_offset), F(-D,0,0), r x F about Y = r_x*F_z - r_z*F_x = 0 - Z_offset*(-D) = + Z_offset*D -> positive M_y means nose-up.)
        M_fin_y = cfg.fin.Z_offset * D_fin;
    end
    aeroForces.M_fin = [0; M_fin_y; 0];


    % % 2b — Assemble Jacobian and apply Newton update (Steps 14–17) %%
    %         Compute residuals for trim (force/moment balance):
    %     R2 = Fz_total - W;        % vertical force residual (Lift_total - Weight)
    %     R3 = M_pitch;             % pitch moment residual (about CG, target 0) 
    %     R1 = Fx_total;            % axial force residual (Thrust - Drag_total)
    % Compute partial derivatives for Jacobian matrix J = d(R)/d[alpha, delta, T]:
    % Partial w.rt alpha (angle-of-attack):
    % Wing lift ∂L_w/∂α ≈ 0.5*rho*U^2*S * d(CL)/dα = 0.5*rho*U^2*S_wing * (2π)  (approx, per rad)
    % dLw_dalpha = 0.5 * rho * U_inf^2 * S_wing * (2*pi);
    % Tail lift ∂L_t/∂α = 0.5*rho*U^2*S_tail * CL_alpha_tail * (d(effTailAoA)/dα).
    % effTailAoA = α - tail_inc - δ, so d(...)/dα = 1.
    dLt_dalpha = 0.5 * rho * V^2 * cfg.tail.S * alpha_tail;
    % dLf_dalpha = 0.5 * rho * V^2 * cfg.fin.S * alpha_fin;
    % Wing drag ∂D_w/∂α ≈ 0.5*rho*U^2*S_wing * d(CD)/dα. For small α, induced drag derivative from lift: d(CD_ind)/dα ≈ 2*CL_w*(dCL/dα)/ (pi*AR).
    % dDw_dalpha = 0.5 * rho * V^2 * S_wing * (2 * CL_w * (2*pi) / (pi * AR));
    % Tail drag ∂D_t/∂α similarly:
    dDt_dalpha = 0.5 * rho * V^2 * cfg.tail.S * (2 * CL_tail * alpha_tail / (pi*4));
    dDf_dalpha = 0.5 * rho * V^2 * cfg.fin.S * (2 * CD_fin  / (pi*4));
    % Moment partials w.rt α:
    % ∂M/∂α = - (dLw_dalpha*x_w + dLt_dalpha*x_tail)
    dM_dalpha = - ( dLt_dalpha * cfg.tail.x) - cfg.fin.Z_offset*dDf_dalpha;
    % Partial w.rt elevator deflection δ:
    % Wing forces do not directly depend on δ (no wing flaps in this model), so:
    dLw_ddelta = 0; dDw_ddelta = 0; dM_ddelta_w = 0;
    % Tail: ∂L_t/∂δ = -0.5*rho*U^2*S_tail * CL_alpha_tail  (note: δ appears with negative sign in effTailAoA if δ defined as trailing-edge up positive)
    % Assuming positive δ = trailing-edge up (which *reduces* tail lift), effTailAoA = α - tail_inc - δ, so ∂L_t/∂δ = -0.5*rho*U^2*S_tail*CLalpha_tail.
    dLt_ddelta = -0.5 * rho * V^2 * cfg.tail.S * alpha_tail;
    % ∂D_t/∂δ ~ - similarly small (neglect for simplicity or compute as function of CL_tail):
    dDt_ddelta = -0.5 * rho * V^2 * cfg.tail.S * (2 * CL_tail * alpha_tail / (pi*4));
    dDf_ddelta = -0.5 * rho * V^2 * cfg.fin.S * (2 * CD_fin  / (pi*4));
    % ∂M/∂δ = - (dLw_ddelta*x_w + dLt_ddelta*x_tail) = - dLt_ddelta * x_tail (since wing part 0)
    dM_ddelta = - dLt_ddelta * cfg.tail.x;
    % Partial w.rt thrust T:
    % Only R3 (axial force) directly depends on T (linearly), and possibly a negligible pitching moment if thrust line offset (ignored here).
    dR3_dT = 1;                % ∂(T - D_total)/∂T = 1
    dR1_dT = 0;                % vertical force does not depend on thrust
    dR2_dT = 0;                % pitch moment assumed independent of thrust (no offset)
    % Fill Jacobian matrix:
    % Rows: [dR1; dR2; dR3], Cols: [dα, dδ, dT]
    J(2,1) = dLt_dalpha;        % ∂R1/∂α = ∂(L_w+L_t-W)/∂α
    J(2,2) = dLt_ddelta;                     % ∂R1/∂δ = ∂L_t/∂δ  (wing has no δ effect)
    J(2,3) = 0;                              % ∂R1/∂T = 0
    J(3,1) = dM_dalpha;                      % ∂R2/∂α 
    % J(3,2) = dM_ddelta;                      % ∂R2/∂δ 
    J(3,2) = 0;                      % ∂R2/∂δ 
    J(3,3) = 0;                              % ∂R2/∂T 
    J(1,1) = - (dDf_dalpha +dDt_dalpha);    % ∂R3/∂α = -∂(D_w+D_t)/∂α
    % J(1,2) = - (dDf_ddelta+dDt_ddelta);                   % ∂R3/∂δ = -∂(D_t)/∂δ
    J(1,2) = 0;                   % ∂R3/∂δ = -∂(D_t)/∂δ
    J(1,3) = dR3_dT;                         % ∂R3/∂T = 1
    aeroForces.J = J;
            % dX = - J \ trimResidual;
    end