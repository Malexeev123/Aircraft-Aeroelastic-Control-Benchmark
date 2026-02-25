%NONLINEAR_TERMS  Vectorised evaluation of intrinsic beam nonlinearities.
%   N = nonlinear_terms(idx, S, par)
%
%   Faster drop‑in replacement for the earlier loop‑based version. The
%   routine is almost allocation‑free and avoids explicit for‑loops over
%   the modal index *l* by exploiting tensor expansion (R2023b pagemtimes)
%   and compile‑time persistent indices.
%
%   INPUT (grouped in two structs for brevity)
%   ----------------------------------------------------------
%   idx    : struct with integer index ranges into state vector
%            generated once in SimRunner ctor via buildIndexStruct()
%
%            .q1    [1 Nm]
%            .q2    [1 Nm]
%            .qxi   [1 Nm+1]
%            .qGam  [1 Na2]
%            .chi   [1 3]
%
%   S      : state vector at current time step  (column)
%
%   par    : struct with all constant tensors
%            .Gamma1      (Nm×Nm×Nm)
%            .Gamma2      (Nm×Nm×Nm)
%            .Gamma_g     (Nm×Nm×Nm+1)
%            .Gamma_xi    (Nm+1×Nm×Nm+1)
%            .forces_0    (Nm×1)   aerodynamic linear load coeff.
%            .Bw, Dw      gust maps (Na2×1)
%            .Bdel, Ddel, Bddel, Dddel  control surf maps
%            .scaleAero   scalar a t_∞ from cfg for quick switching
%            .dt          timestep (scalar)
%
%   OUTPUT
%   ----------------------------------------------------------
%   N_terms : (Nx×1) nonlinear increment *already multiplied by dt*
%%  -------------------------------------------------------------------------
function N_terms = nonlinear_terms(S, par, idx)

% ---------- extract modal sub‑vectors (no find/ismember in loop) ---------
q1      = S(idx.q1 );            % Nm
q2      = S(idx.q2 );            % Nm
q_xi    = S(idx.qxi);            % Nm+1
% q_Gamma = S(idx.qGam);           %#ok<NASGU>  % Na2 (unused in open loop)
% chi_hat = S(idx.chi );           %#ok<NASGU>  % 3   (unused in open loop)

Nm = numel(idx.q1);
Nx = numel(S);
N_terms = zeros(Nx,1);
% ---------- Γ₁ & Γ₂  contractions  (vectorised) ---------------------------
%    −Σ_l  q1_l Γ₁(:,:,l) q1    = - Σ_l  diag(q1) *
%                                  vec(Gamma1)^T * kron(q1,I)
% Use pagemtimes (R2023b) :
G1q1 = squeeze(pagemtimes(par.Gamma1, 'none', q1, 'none'));  % Nm×Nm×1  <- Σ over 3rd dim
% G1q1 = squeeze(pagemtimes(permute(par.Gamma1,[1 3 2]), 'none', q1, 'none'));  % Nm×Nm×1  <- Σ over 3rd dim
N1_G1 = -(G1q1(:,:,1) * q1);
% N(idx.q1) = N(idx.q1) + N1_G1;

G2q2 = squeeze(pagemtimes(par.Gamma2, 'none', q2, 'none'));
% G2q2 = squeeze(pagemtimes(permute(par.Gamma2,[1 3 2]), 'none', q2, 'none'));
N1_G2 = -(G2q2(:,:,1) * q2);
% N(idx.q2) = N(idx.q2) + N1_G2;

% ---------- Γ_g  term -----------------------------------------------------
Ggqx = squeeze(pagemtimes(par.Gamma_g, 'none', q_xi, 'none'));

% Ggqx = squeeze(pagemtimes(permute(par.Gamma_g, [1 3 2]), 'none', q_xi, 'none'));
N1_Gg =  (Ggqx(:,:,1) * q_xi);
% N(idx.q2) = N(idx.q2) + N1_Gg;


% Gxiq1 = squeeze(pagemtimes(permute(par.Gamma_xi, [1 3 2]), 'none', q_xi, 'none'));
% N_qxi =  (Gxiq1(:,:,1) * q1);

% ---------- Γ_xi  term  (Nm+1×1) -----------------------------------------
% Gxiq1 = squeeze(pagemtimes(par.Gamma_xi, 'none', q1, 'none'));
% N_qxi =  (Gxiq1(:,:,1) * q_xi);
%
% disp(size(par.Gamma_xi))
Gxiq1 = squeeze(pagemtimes(par.Gamma_xi, 'none', q_xi, 'none'));
% disp(size(Gxiq1))

N_qxi =  (Gxiq1(:,:,1) * q1);

% ---------- Γ2^T q1  (Eq 9b) ----------------------------------------------
GT = squeeze(pagemtimes(par.Gamma2, 'transpose', q1, 'none'));
N_GT =  GT(:,:,1)*q2;     % note: transpose on first arg

% ---------- aerodynamic linear forces  -----------------------------------
% N_force = par.gustSet.a*par.gustSet.t_inf  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = par.scaleAero  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = par.scaleA  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = 0.5*par.scaleA  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = .5*par.scaleAero  * par.forces_0; % Current Correct
% N_force = .45*par.scaleAero  * par.forces_0; % Current Correct
% N_force = .465*par.scaleAero  * par.forces_0; % Current Correct
% N_force = .475*par.scaleAero  * par.forces_0; % Current Correct
% N_force = par.scaleAero  * par.forces_0; % Current Correct
% N_force = .5*par.scaleA  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = 1/(par.scaleA)  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = .435*par.scaleAero  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force = .4082*par.scaleAero  * par.forces_0;  % Nm × 1  (time‑invariant in open loop)
% N_force =  par.forces_0;  % Nm × 1  (time‑invariant in open loop)

N_force = par.Fscale *par.forces_0; % Current Correct

% ---------- assemble -----------------------------------------------------
N_terms(idx.q1 ) = N1_G1 + N1_G2 + N1_Gg + N_force+par.N_Thrust;  % Eq 9a
% N_terms(idx.q1 ) = N1_G1 + N1_G2  + N_force+par.N_Thrust;  % Eq 9a
% N_terms(idx.q1 ) = N1_G1 + N1_G2 + N1_Gg +par.N_Thrust;  % Eq 9a
N_terms(idx.q2 ) =  N_GT; % Eq 9b
N_terms(idx.qxi) = N_qxi;
% N_terms(idx.qxi) = zeros(size(N_qxi));

% ---------- gust contribution (if provided) --------------------------------
if isfield(par,'gust') && ~isempty(par.gust)
    g   = par.gust;
    % g   = 0.5*par.gust;
    % g   = 0.525*par.gust;
    % g   = 0.5125*par.gust;
    NgB = par.Bw * g;                 % circulatory part
    % NgD = par.scaleAero * par.Dw * g; % non‑circulatory part
    % NgD = par.gustSet.a*par.gustSet.t_inf * par.Dw * g; % non‑circulatory part
    % NgD = par.scaleAero * par.Dw * g; % non‑circulatory part
    NgD = par.Fscale * par.Dw * g; % non‑circulatory part
    
    % NgD =  par.Dw * g; % non‑circulatory part
    N_terms(idx.qGam) = N_terms(idx.qGam) + NgB;  % Γ̇ input
    N_terms(idx.q1 ) = N_terms(idx.q1 ) + NgD;    % modal forces
end

% ---------- control‑surface contribution (if provided) ---------------------
if isfield(par,'u_ctrl') && ~isempty(par.u_ctrl)
    %% Temp 
    % t_inf = .0013;
    % % t_inf = 1;
    % q_inf = 980;
    % aeroScale = 2.45;

    u   = par.u_ctrl;    % [δ₁ δ₂ δ̇₁ δ̇₂]ᵀ
    NcB = (1/par.t_inf)*par.Bdel*u(1:2) + par.Bddel*u(3:4);
    NcD = ( par.scaleA*par.Ddel*u(1:2) + par.scaleAero*par.Dddel*u(3:4) );
    % NcB = (1/t_inf)*par.Bdel*u(1:2) + par.Bddel*u(3:4);
    % NcD = ( aeroScale*par.Ddel*u(1:2) + aeroScale*t_inf*par.Dddel*u(3:4) );

    N_terms(idx.qGam) = N_terms(idx.qGam) + NcB;
    % test rate projection for trim
    % if isfield(par, 'RateProject')
    %     if par.RateProject.projSet
    %         NcD = par.RateProject.Pz*NcD;
    %     end
    % 
    % end
    N_terms(idx.q1 ) = N_terms(idx.q1 ) + NcD;
end
                                          

% N_terms = par.dt * N_terms;   % IMEX increment
end





% function N_terms = nonlinear_terms(num_modes, state, Gamma1, Gamma2, Gamma_g, ...
%     Gamma_xi, dt, num_aero_states2, forces_aero_beam_dof, input_vars, gust_input_i)
%     % This function computes the "nonlinear" part of the IMEX scheme
%     % for the full state vector, with the indexing matching Artola's setup:
%     %
%     %   state = [ q1; q2; q_xi; q_Gamma; hat_chi ]
%     %
%     % where:
%     %   q1 has num_modes entries
%     %   q2 has num_modes entries
%     %   q_xi has (num_modes + 1) entries
%     %   q_Gamma has num_aero_states2 entries
%     %   hat_chi has 3 entries (Euler angles, etc.)
%     %
%     % Equations of motion from eq. (9a)-(9c) in the reference, for example.
% 
%     % --- Extract the relevant pieces from 'state' ---
%     q1_start      = 1;
%     q1_end        = num_modes;
%     q2_start      = q1_end + 1;
%     q2_end        = q1_end + num_modes;           % 2*num_modes
%     q_xi_start    = q2_end + 1;                   % 2*num_modes+1
%     q_xi_end      = q2_end + (num_modes+1);       % 3*num_modes+1
%     q_gamma_start = q_xi_end + 1;                 % 3*num_modes+2
%     q_gamma_end   = q_xi_end + num_aero_states2;  % 3*num_modes+1 + num_aero_states2
%     hat_chi_start = q_gamma_end + 1;              % 3*num_modes+1 + num_aero_states2 +1
%     hat_chi_end   = hat_chi_start + 2;            % up to +2 -> last 3
% 
%     % Now parse them out:
%     q1      = state(q1_start : q1_end);
%     q2      = state(q2_start : q2_end);
%     q_xi    = state(q_xi_start : q_xi_end);
%     q_gamma = state(q_gamma_start : q_gamma_end);
%     hat_chi = state(hat_chi_start : hat_chi_end);
% 
%     % Control input parse:
%     icontrol_input = input_vars.control_inputs;
%     if size(icontrol_input,1) == 1 % column
%         icontrol_input = icontrol_input';
%     end
%     u_delta_i = icontrol_input(1:2);
%     u_ddelta_i = icontrol_input(3:4);
% 
% 
%     eta_thrust= input_vars.eta_thrust;
% 
%     N_len = 2*num_modes + (num_modes+1) + num_aero_states2 + 3;
%     N1    = zeros(N_len,1);
%     %===== Γ₁ Term: -q1_l * Γ₁^l * q1 (Equation 9a) =====
%     N1_Gamma1 = zeros(num_modes, 1);
%     for l = 1:num_modes
%         Gamma1_l = Gamma1(:, :, l);  % Slice of Γ₁ for mode l
%         N1_Gamma1 = N1_Gamma1 - q1(l) * (Gamma1_l * q1);
%     end
% 
%     %===== Γ₂ Term: -q2_l * Γ₂^l * q2 (Equation 9a) =====
%     N1_Gamma2 = zeros(num_modes, 1);
%     for l = 1:num_modes
%         Gamma2_l = Gamma2(:, :, l);  % Slice of Γ₂ for mode l
%         N1_Gamma2 = N1_Gamma2 - q2(l) * (Gamma2_l * q2);
%     end
% 
%     %===== Γ_g Term: +qξ_l * Γ_g^l * qξ (Equation 9a) =====
%     N1_GammaG = zeros(num_modes, 1);
%     for l = 1:num_modes %+ 1  % Loop over num_modes + 1 for q_xi
%         GammaG_l = Gamma_g(:, :, l);
%         N1_GammaG = N1_GammaG + q_xi(l) * (GammaG_l * q_xi);
%     end
% 
%     %===== Γ_xi Term: +qξ_l * Γ_xi^l * q1 (Equation 9c) =====
%     N1_GammaXi = zeros(num_modes + 1, 1);  % Correct size for q_xi
%     for l = 1:num_modes+ 1
%         GammaXi_l = Gamma_xi(:, :, l);  % Slice of Γ_xi for mode l
%         N1_GammaXi = N1_GammaXi + q_xi(l) * (GammaXi_l * q1);
%     end
% 
%     %===== Γ₂^⊤ Term: q2_l*(Γ₂^l)^⊤ * q1 (Equation 9b) =====
%     N2_Gamma2T = zeros(num_modes, 1);
%     for l = 1:num_modes
%         Gamma2_l = Gamma2(:, :, l);
%         N2_Gamma2T = N2_Gamma2T + q2(l) * (Gamma2_l' * q1);
%     end
% 
%     %===== Gust to N1 =====
%     Bw = input_vars.Bw;
%     Dw = input_vars.Dw;
%     aero_scale = input_vars.aero_scale;
%     t_inf = input_vars.t_inf;
%     if input_vars.scale_gust_choice ==1
% 
%         N_gustB = Bw* gust_input_i;
%         N_gustD = aero_scale*t_inf*Dw * gust_input_i;
%     else
%         N_gustB = Bw* gust_input_i;
%         N_gustD = Dw * gust_input_i;
%         aero_scale =1;
%     end
% 
%     %===== Controls to N1 =====
%     B_delta = input_vars.B_delta;
%     D_delta= input_vars.D_delta;
%     B_ddelta = input_vars.B_ddelta;
%     D_ddelta= input_vars.D_ddelta;
% 
%     if input_vars.scale_gust_choice ==1
% 
%         control_ddeltaB = B_ddelta*u_ddelta_i;
%         control_deltaB = (1/t_inf)* B_delta*u_delta_i;
% 
%         control_ddeltaD = aero_scale*t_inf*D_ddelta*u_ddelta_i;
%         control_deltaD = aero_scale*D_delta*u_delta_i;
%     else
%         control_ddeltaB = B_ddelta*u_ddelta_i;
%         control_deltaB = B_delta*u_delta_i;
% 
%         control_ddeltaD = D_ddelta*u_ddelta_i;
%         control_deltaD = D_delta*u_delta_i;
%         aero_scale =1;
%     end
%     N_controlB = control_deltaB+ control_ddeltaB;
%     N_controlD = control_ddeltaD+control_deltaD;
%     %===== Assign to N1 =====
%     N1(q1_start:q1_end) = N1_Gamma1 + N1_Gamma2 + N1_GammaG + eta_thrust + aero_scale*forces_aero_beam_dof + N_gustD+N_controlD;  % Equation 9a
%     % N1(q1_start:q1_end) = N1_Gamma1 + N1_Gamma2 + N1_GammaG;  % Equation 9a
%     N1(q2_start:q2_end) = N2_Gamma2T;             % Equation 9b
%     N1(q_xi_start: q_xi_end) = N1_GammaXi;  % Equation 9c
%     N1(q_gamma_start:q_gamma_end) = N_gustB + N_controlB;
%     % disp(norm(N1))
%     N_terms = N1 * dt;
%     % % N_terms = N1 ;
% end
