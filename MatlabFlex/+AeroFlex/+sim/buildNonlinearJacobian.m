% function fJac = buildNonlinearJacobian(G1,G2,Gg,Gxi,idx)

function [jacobian_nl, Nu] = buildNonlinearJacobian(x_n,idx, par)
    % from solve_projected_system()
    % --- Extract the relevant pieces from 'state' ---
    % q1_start      = 1;
    % q1_end        = num_modes;
    % q2_start      = q1_end + 1;
    % q2_end        = q1_end + num_modes;           % 2*num_modes
    % q_xi_start    = q2_end + 1;                   % 2*num_modes+1
    % q_xi_end      = q2_end + (num_modes+1);       % 3*num_modes+1
    % q_gamma_start = q_xi_end + 1;                 % 3*num_modes+2
    % q_gamma_end   = q_xi_end + num_aero_states2;  % 3*num_modes+1 + num_aero_states2
    % hat_chi_start = q_gamma_end + 1;              % 3*num_modes+1 + num_aero_states2 +1
    % hat_chi_end   = hat_chi_start + 2;            % up to +2 -> last 3
    q1      = x_n(idx.q1 );            % Nm
    q2      = x_n(idx.q2 );            % Nm
    qxi    = x_n(idx.qxi);
    G1 = par.Gamma1;
    G2 = par.Gamma2;
    Gg = par.Gamma_g;
    Gxi = par.Gamma_xi;
    
    Naero = par.Na;
    num_aero_states2 = Naero;
    Nm   = numel(q1);
    Nxi  = Nm+1;
    nx   = 3*Nm + 1 + Naero+ 3;              % total state dimension
    % m    = 2;                             % # inputs  (δ̇ , w_g)
    ctrlIdx1 = size(par.Bdel, 2);
    ctrlIdx2 = size(par.Bddel,2);
    m    = size(par.Bw,2) + ctrlIdx1 + ctrlIdx2;                             % # inputs  (δ̇ , w_g)
    % % Now parse them out:
    % q1      = state(q1_start : q1_end);
    % q2      = state(q2_start : q2_end);
    % q_xi    = state(q_xi_start : q_xi_end);
    % q_gamma = state(q_gamma_start : q_gamma_end);
    % hat_chi = state(hat_chi_start : hat_chi_end);
    % 
    % % Control input parse:
    % icontrol_input = input_vars.control_inputs;
    % u_delta_i = icontrol_input(1:2);
    % u_ddelta_i = icontrol_input(3:4);

    num_modes = Nm;

    
% %%  
% 
%     q1   = x_n(1:Nm);
%     q2   = x_n(Nm+1:2*Nm);
%     qxi  = x_n(2*Nm+1 : 3*Nm+1);   
    % q_Gamma, hat_chi are purely linear → zero in ∂N/∂x.
    
    % -------------------------------------------------------------------------
    % helper: contraction along 3rd index
    % -------------------------------------------------------------------------
        % function M = contr3(T,v)
        %     % T : p×q×r   , v : r×1   →  M : p×q   (matrix)
        %     M = reshape(permute(T,[1 2 3]),[],size(T,3)) * v;
        %     M = reshape(M,size(T,1),[]);
        % end
    % -------------------------------------------------------------------------
        % function y = contr3(T,v)
        % % T :  p×q×r ,   v : r×1   →   y : p×q
        % S  = permute(T,[1 2 3]);          % keep natural order
        % S  = reshape(S,[],size(S,3));     % (p*q) × r
        % y  = S * v;                       % (p*q) × 1
        % y  = reshape(y,size(T,1),[]);     % p × q
        % end
    
    N1_Gamma1 = zeros(Nm, Nm);
    N1_Gamma1_l = zeros(Nm, Nm);
    for j = 1:num_modes
        Gamma1_j = squeeze(G1(:, j, :));
        N1_Gamma1j = -q1(j) * Gamma1_j;
        for l = 1:num_modes
            Gamma1_l = G1(:, :, l);  % Slice of Γ₂ for mode l
            N1_Gamma1_tmp =  - q1(l) * (Gamma1_l );
            % N1_Gamma1_l = N1_Gamma1_l + N1_Gamma1_tmp;
            N1_Gamma1_l =  N1_Gamma1_tmp;
        end
        N1_Gamma1 = N1_Gamma1 +N1_Gamma1_l+N1_Gamma1j;
    end
   
    %===== Γ₂ Term: -q2_l * Γ₂^l * q2 (Equation 9a) =====
    N1_Gamma2 = zeros(num_modes, num_modes);
    N1_Gamma2_l= zeros(num_modes, num_modes);
    for j = 1:num_modes
        Gamma2_j = squeeze(G2(:, j, :));
        N1_Gamma2j = -q2(j) * Gamma2_j;

        for l = 1:num_modes
            Gamma2_l = G2(:, :, l);  % Slice of Γ₂ for mode l
            N1_Gamma2_tmp = - q2(l) * (Gamma2_l);
            % N1_Gamma2_l = N1_Gamma2_l + N1_Gamma2_tmp;
            N1_Gamma2_l =  N1_Gamma2_tmp;
        end
        N1_Gamma2 = N1_Gamma2 + N1_Gamma2_l + N1_Gamma2j ;
    end

    %===== Γ_g Term: +qξ_l * Γ_g^l * qξ (Equation 9a) =====
    N1_GammaG = zeros(num_modes, num_modes+1);
    N1_GammaG_l= zeros(num_modes, num_modes+1);
    for j = 1:num_modes +1
        GammaG_j = squeeze(Gg(:, j, :));
        N1_GammaGj = qxi(j) * GammaG_j;

        for l = 1:num_modes + 1  % Loop over num_modes + 1 for q_xi
            GammaG_l = Gg(:, :, l);

            N1_GammaG_tmp = qxi(l) * (GammaG_l );
            % N1_GammaG_l = N1_GammaG_l + N1_GammaG_tmp;
            N1_GammaG_l =  N1_GammaG_tmp;
        end
        N1_GammaG = N1_GammaG + N1_GammaG_l + N1_GammaGj ;
        
    end

    jacobian_nl_1 = [N1_Gamma1, N1_Gamma2, N1_GammaG, zeros(num_modes, num_aero_states2+3)];
    

    % 
    % Nxi = num_modes+1;
    % Nm = num_modes;

   
    %%====== N2 terms==========%%
    %===== Γ₂^⊤ Term: q2_l*(Γ₂^l)^⊤ * q1  =====
    N2_Gamma2T = zeros(num_modes, num_modes);
    for l = 1:num_modes
        Gamma2_l = G2(:, :, l);
        N2_Gamma2T = N2_Gamma2T + q2(l) * (Gamma2_l');
    end

    N2_q1_Gamma2T = zeros(num_modes, num_modes);
    for j = 1:num_modes
        Gamma2_jT = permute(G2,[3 2 1]); 
        Gamma2_jT_s = q1(j) *squeeze(Gamma2_jT(:, j, :));
        N2_q1_Gamma2T = N2_q1_Gamma2T+Gamma2_jT_s;
    end

    jacobian_nl_2 = [N2_Gamma2T, N2_q1_Gamma2T, zeros(num_modes, num_modes+1+num_aero_states2+3)];

    %%====== N_xi terms==========%%
    %===== Γ_xi Term: +qξ_l * Γ_xi^l * q1  =====
    N1_GammaXi = zeros(num_modes + 1, 1);  % Correct size for q_xi
    for l = 1:num_modes+ 1
        GammaXi_l = Gxi(:, :, l);  % Slice of Γ_xi for mode l
        N1_GammaXi = N1_GammaXi + qxi(l) * (GammaXi_l );
    end

    N1_q1_GammaXi = zeros(num_modes+1, num_modes+1);
    for j = 1:num_modes
        Gammaxi_j = squeeze(Gxi(:, j, :)); 
        Gammaxi_jT_s = q1(j) *Gammaxi_j;
        N1_q1_GammaXi = N1_q1_GammaXi + Gammaxi_jT_s;
    end

    jacobian_nl_xi = [N1_GammaXi, zeros(num_modes+1,num_modes), N1_q1_GammaXi, zeros(num_modes+1, num_aero_states2+3)];

    jacobian_nl = [jacobian_nl_1; jacobian_nl_2; jacobian_nl_xi; zeros(num_aero_states2+3,nx )];
    
    %% Need to also get u-input for dh/du in Nl but is linear
    % dqaendq0 = jacobian_nl * dqdq0;
    


   Nu = zeros(nx, nx+m);
    % Nu(:,1) = Bdelta;      % flap‑rate
    % Nu(:,2) = Bgust;       % gust
    Nu(idx.q1, nx+1) = par.gustSet.a*par.gustSet.t_inf * par.Dw;
    Nu(idx.qGam, nx+1) = par.Bw;
    
    Nu(idx.q1,nx+ 2:nx+1 + ctrlIdx1) = par.gustSet.a*par.Ddel;
    Nu(idx.q1, nx+2 + ctrlIdx1:nx+1 + ctrlIdx1 + ctrlIdx2) = par.scaleAero*par.Dddel;
    
    Nu(idx.qGam, nx+2:nx+1 + ctrlIdx1) = (1/par.t_inf)*par.Bdel;
    Nu(idx.qGam, nx+2 + ctrlIdx1:nx+1 + ctrlIdx1 + ctrlIdx2) =  par.Bddel;
    
end
% %BUILDNONLINEARJACOBIAN  returns a handle J = fJac(x,u,w)
% %
% % INPUT tensors (all column‑major, 1‑based):
% %   G1 ,G2   :  Γ₁, Γ₂  (Nm × Nm × Nm)
% %   Gg       :  Γ_g      (Nm × (Nm+1) × (Nm+1))
% %   Gxi      :  Γ_ξ      ((Nm+1) × Nm × (Nm+1))
% %   idx      :  struct with .q1 .q2 .qxi ranges inside full x vector
% %
% % OUTPUT
% %   fJac(x,u,w)  – returns sparse Jacobian ∂N/∂x   (size Nx×Nx)
% 
% Nm    = size(G1,1);   Nxi = Nm+1;
% 
% % pre‑compute tensor reshapes for pagemtimes -----------------------------
% G1r = reshape(G1 ,[Nm,Nm*Nm]);
% G2r = reshape(G2 ,[Nm,Nm*Nm]);
% GgL = reshape(permute(Gg,[1 3 2]),[Nm,Nxi*Nxi]);  % Γ_g * qξ , qξ*Γ_gᵀ
% GxiL= reshape(Gxi,[Nxi*Nxi,Nm]);                  % Γ_ξ * q1
% 
% fJac = @(x,u,w) localJac(x);
% 
%     function J = localJac(x)
%         % --- slice modal subsystems from state vector ---------------
%         q1  = x(idx.q1);
%         q2  = x(idx.q2);
%         qxi = x(idx.qxi);
% 
%         % --- build blocks via tensor contractions -------------------
%         N1  = -G1r * kron(q1,q1) ...
%               -G2r * kron(q2,q2) ...
%               +GgL * kron(qxi,qxi);
% 
%         N2  = -G2r.' * kron(q1,q2);          % Γ₂ᵀ block
%         Nxi = +GxiL.'* kron(q1,qxi);         % Γ_ξ block
% 
%         % assemble sparse (only top‑left 3Nm rows are non‑zero)
%         Nx  = numel(x);
%         J   = sparse(3*Nm+1,Nx);   % pre‑allocate
%         J(1:Nm,1:2*Nm)       = reshape(N1,[Nm,2*Nm]);
%         J(Nm+1:2*Nm,1:2*Nm)  = reshape(N2,[Nm,2*Nm]);
%         J(2*Nm+1:3*Nm+1,1:3*Nm+1) = ...
%                          [zeros(Nxi,2*Nm) reshape(Nxi,[Nxi,Nxi])];
%     end
% end
