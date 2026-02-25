function [q1,q2,qxi,qg,chi] = oneEulerStep( ...
        q1,q2,qxi,qg,chi, deltaE,Tset,alpha, beam, aero, cfg, base)
% oneEulerStep  –  forward‑Euler update of intrinsic modal states
%
% This lightweight step is used *only* inside the trim loop,
% where Δt is chosen ≪ aero t‑inf so that the solution
% drifts quasi‑statically toward equilibrium.
%
% INPUT
%   q1,q2,qxi,qg,chi : current modal states              (column vectors)
%   deltaE           : elevator deflection    [rad]
%   Tset             : thrust scalar          [N]
%   alpha            : angle‑of‑attack        [rad]
%   beam             : BeamModel  (contains Omega, Γ‑tensors,…)
%   aero             : AeroROM   (contains ROM_cts, Bw, Bδ,…)
%   cfg              : master configuration struct
%
% OUTPUT (updated states after one Δt)
%   q1,q2,qxi,qg,chi
%
% NOTES
%   • Uses the same equations (Artola (28)) already
%     coded in your IMEX integrator, but with Δt_trim
%     and explicit Euler (no LU solve needed here).
%   • Stiff aero γ‑states are integrated implicitly INSIDE
%     the ROM (AΓ diagonal) → stable even with large eigenvalues.
%

% ---------- short symbols ---------------------------------------------
Nm  = beam.Nm;
Na2 = cfg.Na;

dt  = cfg.trim.dtQS;          % e.g. 2e‑4  (set in nominalConfig)

Om  = beam.Omega;             % diag(Nm)
G1  = beam.Gamma1;            % Nm×Nm×Nm
G2  = beam.Gamma2;
Gg  = base.Gamma_g;
Gxi = base.Gamma_xi;

% aerodynamic scaling
qInf = 0.5*cfg.flight.rho*cfg.flight.U_inf^2;
tInf = cfg.flight.b_ref / cfg.flight.U_inf;
aScl = qInf * cfg.flight.b_ref^2;

% ROM matrices already scaled to continuous time
Agamma   = aero.forceMap.A_Gamma;               % Na2×Na2
B1   = aero.forceMap.B1;
B0   = aero.forceMap.B0;
Cgamma   = aero.forceMap.C_Gamma;
D0   = aero.forceMap.D0;
D1   = aero.forceMap.D1;
Bw   = aero.forceMap.Bw; Dw = aero.forceMap.Dw;
Bdelta   = aero.forceMap.B_delta; Ddelta = aero.forceMap.D_delta;
Bddelta  = aero.forceMap.B_ddelta; Dddelta = aero.forceMap.D_ddelta;

% ---------- nonlinear & external terms --------------------------------
% (i) control / disturbance vectors
u_delta  = deltaE * ones(length(aero.aeroMesh.surface_m),1);         % 2×1
u_dDelta = 0;                                         % rate = 0 here
gust     = 0;

% (ii) modal thrust -> η term  (simple lookup)
eta = beam.phi1.' * Tset*ones(beam.nFlex,1)+ aScl*aero.DataMatrix.forces_aero_beam_dof;

% ---------- structural RHS --------------------------------------------
% eq.(28a): q1̇
% begin with Ω q2
q1dot = Om * q2;

% − q1_l Γ1^l q1   − q2_l Γ2^l q2   + qxi_l Γg^l qxi
for l = 1:Nm
    q1dot = q1dot ...
        - q1(l)* (G1(:,:,l) * q1) ...
        - q2(l)* (G2(:,:,l) * q2) ...
        + qxi(l)* (Gg(:,:,l) * qxi);
end

% + aerodynamic modal forces (linear part)
q1dot = q1dot + aScl*( Cgamma*qg ...
                    + D1*q1 ...
                    - D0*(Om\q2) ...
                    + Ddelta*u_delta ...
                    + tInf*Dddelta*u_dDelta ...
                    + Dw*gust);

% + thrust modal load
q1dot = q1dot + eta;

% eq.(28b): q2̇
q2dot = - Om*q1;
for l = 1:Nm
    q2dot = q2dot + q2(l)*(G2(:,:,l).')*q1;
end

% eq.(28c): qxi̇   (length Nm+1)
qxidot = zeros(Nm+1,1);
for l = 1:Nm+1
    qxidot = qxidot + qxi(l)*(Gxi(:,:,l) * q1);
end

% eq.(28d): qΓ̇  (integrated implicitly via (I‑dt/tInf AΓ)⁻¹)
rhs_gamma = B1*q1 ...
     - (B0*(Om\q2))/tInf ...
     + (Bdelta*u_delta)/tInf + Bddelta*u_dDelta ...
     + (base.FM.Bchi/tInf)*chi;
qg   = (eye(Na2)-dt/tInf*Agamma)\(qg + dt*rhs_gamma);

% eq.(28e): χ̇  = X q1
chi   = chi + dt*(aero.forceMap.X * q1);

% ---------- explicit Euler updates for remaining states ---------------
q1  = q1  + dt*q1dot;
q2  = q2  + dt*q2dot;
qxi = qxi + dt*qxidot;
q1 = q1(:, end);
q2 = q2(:, end);
qxi = qxi(:, end);
end
