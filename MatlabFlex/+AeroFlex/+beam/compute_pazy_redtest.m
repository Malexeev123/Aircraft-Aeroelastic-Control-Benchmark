
%% ------------------------------------------------------------------------

function [ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = compute_pazy_redtest( fem, aero,phi0, ...
    Mred, Kred, Nm, aoa, chain, w0, dubugPlot, rigid, NetworkPath)
% Notes:
% Will start with linear, will do cubic another day
% eigenvalues_rigid = load("eigenvalues_rigid.dat");
% eigenvectors_rigid = load("eigenvectors_rigid.dat");
%% Frame rotation Global to Local 
    % For full wing:
    idxMinL = min(fem.coordinates(:,2));
    idxMaxL = max(fem.coordinates(:,2));
    L = idxMaxL- idxMinL;
    discFact = 1; % for finer grid. Do 1 for same as SHARPy case file.
    nPts = discFact*fem.num_node;   % # of points for discretization
    sVals = linspace(idxMinL,idxMaxL,nPts).'; % 0 for half wing -L for full wing
    beamSegments(1).s0 = idxMinL;
    beamSegments(1).s1 = idxMaxL;
    keepDofs = fem.keepDofs;
    k1 = fem.structural_twist; 
    k2 = aero.twist;
    if k1 ~= k2
      error('somthing wrong double check')
    end 
    kappa0Vals = k2;
    num_elem = fem.num_elem;
    connects = fem.connectivities; 
    if size(kappa0Vals, 1) ~= size(connects, 1)
        kappa0Vals = kappa0Vals';
    end

    if isfield(aero,'sweep' )
        % disp('keke')
        twist0 = aero.twist;
        sweepVar = 1;
        if size(twist0, 1) ~= size(connects, 1)
            twist0 = twist0';
        end
    else
        sweepVar = 0;
    end

    twist_of_s = zeros(nPts,1);
    alpha_of_s = zeros(nPts,1);

    % Factor in Angle of Attack? Rotation about Y-axis (for all nodes)
    % disp(aoa)
    aoa_s = deg2rad(double(aoa))*ones(nPts, 1);

    for i = 1:num_elem
        if i>1 && connects(i,1) ~= connects(i-1, 2)
            % next wing connect
            nodeTip = connects(i-1,2);
            twistTip = kappa0Vals(i-1,2);
            twist_of_s(nodeTip) = twistTip;
            if sweepVar ==1
                alphaTip = aero.sweep(i-1,2);
                alpha_of_s(nodeTip)= alphaTip ;    
            end
        end
        conn_i = connects(i,:);
        node_Ai = conn_i(1);
        node_Bi = conn_i(3);

        twist_Ai = kappa0Vals(i,1);
        twist_Bi = kappa0Vals(i,3);
        twist_of_s(node_Ai) = twist_Ai;
        twist_of_s(node_Bi) = twist_Bi;

        if sweepVar ==1
            sweep_Ai = aero.sweep(i,1);
            sweep_Bi = aero.sweep(i,3);
            alpha_of_s(node_Ai) = sweep_Ai;
            alpha_of_s(node_Bi) = sweep_Bi;
        end
    end    

    xi_bar_array = zeros(nPts,4);
    % phi_omega = zeros(size(phi1,1)/2,Nm);
    kk =1;
    del_pts = 0; del_idx = [];
    for i=1:nPts
        sNow = sVals(i);
        alpha_s = alpha_of_s(i);  %  function for sweep vs. s
        theta_s = twist_of_s(i);  %  function for twist vs. s
        aoa_i = aoa_s(i);
        xi_bar_array(i,:) = combineSweepTwist(alpha_s, theta_s, aoa_i);

       dofRange= double((i-1)*6) + (4:6);
       dofOmega = double((kk-1)*3)+(1:3);
       idx_k = localDofIndex(keepDofs, dofRange);
       if isempty(idx_k)
           del_idx(end+1) = i;
           del_pts = del_pts+1; % to exclude clamped nodes in mode calc
           continue; 
       end
       kk=kk+1;
       % phi_omega(dofOmega,:) = phi1(idx_k,:);
    end

    xi_bar_array(del_idx,:) = []; % delete clamped nodes
    coord_del = fem.coordinates(del_idx, 2);
    idx_coord_del = find(sVals == coord_del);
    sVals(idx_coord_del) = [];
    xi_bar_v = xi_bar_array(:, 2:4);
    xi_bar_0 = xi_bar_array(:, 1);
    activePts = nPts-del_pts;
    Rlocal0 = zeros(activePts*3,3); % Euler
    Rglobal0 = Rlocal0;
    for i = 1:activePts  % Build Init rot field 
        xi_bar_v_i = xi_bar_v(i,:).';
        xi_bar_0_i = xi_bar_0(i);
        Rlmat = [-xi_bar_v_i, xi_bar_0_i*eye(3)+skew(xi_bar_v_i)];
        Rrmat = [-xi_bar_v_i, xi_bar_0_i*eye(3)-skew(xi_bar_v_i)].';
        Rs_idx0 = double((i-1)*3)+(1:3);
        Rlocal_tmp = Rlmat*Rrmat;
        Rlocal0(Rs_idx0,:) = Rlocal_tmp; % This is Local to Global
        RotPsiLoc(Rs_idx0, Rs_idx0)= Rlocal_tmp;
        Rglobal_tmp = inv(Rlocal_tmp);
        Rglobal0(Rs_idx0,:) = Rglobal_tmp; % This is Global to Local
        RotPsiGlob(Rs_idx0, Rs_idx0)= Rglobal_tmp;

    end

    % Test = Rlocal0(1:3, :)* [1;0;0]
    % Test2 = Rglobal0(1:3, :)* [1;0;0]
    phi1_discrete = zeros(size(phi0));
    Psi1_discrete = phi1_discrete;
    %% Build Discrete phi1 
    for j = 1:Nm
        for i = 1: activePts
            idx_R = double((i-1)*3)+(1:3);
            idx_0 = double((i-1)*6)+(1:6);
            % R_i3 = Rlocal0(idx_R,:);
            R_i3 = Rglobal0(idx_R,:);

            R_i6 = blkdiag(R_i3,R_i3);
            if Nm ==20

                % phi1_discrete(idx_0,j) = w0(j)*R_i6*phi0(idx_0,j);
                phi1_discrete(idx_0,j) = w0(j)*phi0(idx_0,j);
                % phi1_discrete(idx_0,j) = -w0(j)*phi0(idx_0,j);
                % phi1_discrete(idx_0,j) = phi0(idx_0,j);
                % phi1_discrete(idx_0,j) = R_i6*phi0(idx_0,j);
                % phi1_discrete(idx_0,j) = -phi0(idx_0,j);
                % phi1_discrete(idx_0,j) = w0(j)*phi0(idx_0,j);

            else
                if NetworkPath
                    % phi1_discrete(idx_0,j) = w0(j)*R_i6*phi0(idx_0,j);

                else
                    % phi1_discrete(idx_0,j) = w0(j)*phi0(idx_0,j);
                    % phi1_discrete(idx_0,j) = w0(j)*R_i6*phi0(idx_0,j);

                end
            end
            % if i >= (activePts/2+1)
            %     phi1_discrete(idx_0,j) = -phi1_discrete(idx_0,j);
            % end

            % phi1_rigid(idx_0,j) = eigenvalues_rigid(j)*R_i6*eigenvectors_rigid(idx_0,j);
        end

    end
    % phi1_discrete = phi0;

    % for k=1:Nm
    %     % simple forward difference along chain for slopes
    %     idx = chain(:,1:2);
    %     phi1(idx(:,1),k) = (phi0(idx(:,2),k)-phi0(idx(:,1),k))./ ...
    %                         vecnorm(fem.coordinates(idx(:,2),:)-fem.coordinates(idx(:,1),:),2,2);
    % end

    % phi1 = phi0;
    %% 1) figure out sub‐element lengths for chain Prog.
    chain_id = 1; % Start with 1 chain, use connects to determine val
    chain_end = [];
    chain_start = connects(1,1);
    for i = 1:num_elem
        if i>1 && connects(i,1) ~= connects(i-1, 2)
            chain_id = chain_id+1;
            chain_start(end+1) = connects(i,1);
            chain_end(end+1) = connects(i-1, 2);
        end
    end
    nModes = Nm;
    ndofs  = size(phi0,1);
    XYZ = fem.coordinates;
    Phi0Fun = cell(1,nModes);
    Phi1Fun = cell(1,nModes);
    seg_ids = 1; seg_ide = length(chain)/chain_id-1;
    % seg = struct([]);
    phi2_discrete = zeros(size(phi0));
    psi2_discrete = zeros(size(phi2_discrete));
    r_mid = zeros(activePts,3);

    %% Interpolation Phi0
    % Discrete Phi0,1 --> Cont. Phi0,1
    Nm_sub = zeros(Nm, Nm); 
    phi1_half = zeros(size(phi1_discrete));
    CorrTerms = struct('A1', {Nm_sub, Nm_sub, Nm_sub}, 'A2', {Nm_sub, Nm_sub, Nm_sub});
    phi1_loc_corrected = struct('Phi1_corr',{phi2_discrete, phi2_discrete, phi2_discrete} );
    phi2_loc_corrected= struct('Phi2_corr',{phi2_discrete, phi2_discrete, phi2_discrete} );
    for sIdx = 1:chain_id
        if sIdx ~= chain_id
            chain_s = chain(chain_start(sIdx):chain_end(sIdx));
        else % only for last chain
            chain_s = chain(chain_start(sIdx):end);

            
        end
        % disp('chain_s');disp(chain_s);
        nChain = length(chain_s);
        elemLen = zeros(1, nChain-1);
        if sIdx ~=chain_id
            for i=1:(nChain-1)
                nA = chain_s(i);
                nB = chain_s(i+1);
                elemLen(i)= norm( XYZ(nB,:) - XYZ(nA,:) );
                elemChainL(i, sIdx) = elemLen(i);
            end
        else
            for i=1:(nChain-1)
            nA = chain_s(i);
            nB = chain_s(i+1);
            elemLen(i)= norm( XYZ(nB,:) - XYZ(nA,:) );
            elemChainL(i, sIdx) = - elemLen(i);
            end
        end
        % disp(elemChainL)
        % E Matrix
        % Eall =zeros(fem.num_elem*6, 6);
        % % e1 = [1; 0; 0];
        % e1 = [0; 1; 0];
        % if isfield(fem,'structural_twist')
        %     kappa_all0 = fem.structural_twist; 
        % else
        %     kappa_all0 = zeros(3, fem.num_elem);
        % end
        % 
        % for i = 1:fem.num_elem
        %     elemL = elemLen(i);
        %     rowi = double((i-1)*6 )+ (1:6);
        %     kappa_0 = kappa_all0(:,i)/elemL;
        %     E_var = [e1;kappa_0];
        %     E6 = L1_operator(E_var);
        %     Eall(rowi, :) = E6;
        % end
        Eall = zeros(fem.num_elem*6,6);
        for i = 1:fem.num_elem
            % local frame for sub-element i  (node i -> i+1)
            nA = connects(i,1);  nB = connects(i,3);
            R3 = buildSegmentRotation2(fem, nA, nB, XYZ);   % local->global
            R6 = blkdiag(R3,R3);

            elemL = elemChainL(i,1);                         % length of sub
            kappa_loc = fem.structural_twist(:,i) / elemL;   % local curvature
            e1_loc = [1;0;0];
            % e1_loc = [0;1;0];

            E6_loc = L1_operator([e1_loc; kappa_loc]);       % operator in local coords
            % disp(E6_loc)
            E6_glb = R6 * E6_loc * R6.';                     % rotate into GLOBAL
            rowi = double((i-1)*6) + (1:6);
            % Eall(rowi,:) = E6_glb;
            Eall(rowi,:) = E6_loc;
        
            
        
        end
        % disp(Eall)
        %% Accumulate Phi2 and Psi2
        phi2_idxs = (chain_s(1)-1)*6+1-6*del_pts*(sIdx-1);
        phi2_idxe = chain_s(end-1)*6-6*del_pts*(sIdx-1);
        seg_int = double((sIdx-1)*8)+(1:8);

        doOmeg = 1; % to multiply Psi2 by 1/omega(j)
        % doOmeg = 0; % to multiply Psi2 by 1/omega(j)
        [phi2_discrete(phi2_idxs:phi2_idxe, :), r_mid(seg_int,:), psi2_discrete(phi2_idxs:phi2_idxe, :), phi1_half(phi2_idxs:phi2_idxe, :)] = accumulatePhi2_singleChain(fem, keepDofs, Kred, ...
            phi0, chain_s, phi1_discrete, w0, elemLen, Eall, doOmeg, Mred);

    end

    %%

%% Psi 1 Calc --> covert cont. (Can skip since dirac delta aka if s==si node then Psi(s) == Psi_si otherwise Psi == 0
for j = 1: Nm
      Psi1_discrete(:,j) = Mred*phi1_discrete(:, j); % This one for Comp
        % Psi1_discrete(:,j) = blkdiag(RotPsiGlob,RotPsiGlob)*Mred*phi1_discrete(:, j); % This one for Comp
        % Psi1_discrete(:,j) = blkdiag(RotPsiLoc,RotPsiLoc)*Mred*phi1_discrete(:, j); % This one for laptop
end

%% Test of M conversion
elemLen_sub1 = elemChainL(:,1);
elemLen_sub2 = elemChainL(:,2);
elemChainL = [elemLen_sub1;elemLen_sub2];


[alpha1_new, alpha2_new] = intrinsicAlphaNormalization( ...
             fem, keepDofs, phi1_discrete, phi2_discrete, Mred, Kred, nModes, ...
             r_mid, chain, elemChainL, Psi1_discrete, psi2_discrete, phi1_half);





% disp('alpha1_check'); disp(alpha1_new)

alpha1_new_symm = 0.5*(alpha1_new + alpha1_new.');
alpha2_new_symm = 0.5*(alpha2_new + alpha2_new.');

% Legacy Code:
[Valpha1, Dalpha1] = eig(alpha1_new_symm);
[Valpha2, Dalpha2] = eig(alpha2_new_symm);

% disp(rank(Valpha2))
alpha1_reconstructed = Valpha1*Dalpha1/Valpha1;
alpha2_reconstructed = (Valpha2*Dalpha2/Valpha2);

% alpha1_diag1 = diag(alpha1_new_symm);
alpha1_diag = diag(alpha1_reconstructed);
% alpha2_diag = diag(alpha2_reconstructed);
alpha2_diag = diag(alpha2_new_symm);

[alpha2_d, alpha2_idx]= find(alpha2_diag<1e-4);
% alpha2_diag(alpha2_d) =1;
alpha2_diag = abs(alpha2_diag);

phi1_corrected = zeros(size(phi1_discrete));
phi2_corrected = zeros(size(phi2_discrete));
psi1_corrected = zeros(size(Psi1_discrete));
psi2_corrected = zeros(size(psi2_discrete));
phi1_half_corrected =  zeros(size(phi1_discrete));
phi0_corrected = zeros(size(phi0));
% Version 1
% for j = 1:Nm
%     % phi1_corrected = phi1_discrete(:,j)/sqrt(alpha1_new(:,j));
%     phi1_corrected(:,j) = phi1_discrete(:,j)/sqrt(alpha1_diag(j));
%     % phi1_corrected2 = phi1_discrete/sqrt(alpha1_new);
%     phi1_half_corrected(:,j) = phi1_half(:,j)/sqrt(alpha1_diag(j));
%     % psi1_corrected = Psi1_discrete(:,j)/sqrt(alpha1_new(:,j));
%     psi1_corrected(:,j) = Psi1_discrete(:,j)/sqrt(alpha1_diag(j));
%     % psi1_corrected2 = Psi1_discrete/sqrt(alpha1_new);
% 
%     phi2_corrected(:,j) = phi2_discrete(:,j)/sqrt(alpha2_diag(j));
%     psi2_corrected(:,j) = psi2_discrete(:,j)/sqrt(alpha2_diag(j));
% 
% end

% Version 2
[alpha1_sc, alpha2_sc] = intrinsicAlphaNormalizationScalar( ...
             fem, keepDofs, phi1_discrete, phi2_discrete, Mred, Kred, nModes, ...
             r_mid, chain, elemChainL, Psi1_discrete, psi2_discrete, phi1_half);
% disp(size(alpha1_sc)); 
% disp(alpha1_sc); 
% disp(alpha2_sc); 

% error('pausing here')
% for i = 1:Nm
    for j = 1:Nm
        phi1_corrected(:,j) = phi1_discrete(:,j)/sqrt(alpha1_sc(j));
        phi1_half_corrected(:,j) = phi1_half(:,j)/sqrt(alpha1_sc(j));
        psi1_corrected(:,j) = Psi1_discrete(:,j)/sqrt(alpha1_sc(j));
    
        phi2_corrected(:,j) = phi2_discrete(:,j)/sqrt(alpha2_sc(j));
        psi2_corrected(:,j) = psi2_discrete(:,j)/sqrt(alpha2_sc(j));
        
        % Norm Phi 0 here
        phi0_corrected(:,j) = w0(j)*phi0(:,j)/sqrt(alpha1_sc(j));
    end
% end



phi1_half = phi1_half_corrected;
% phi2_discrete= (phi2_corrected);
% writematrix(phi2_corrected,'phi2_corrected.dat','Delimiter','tab');
% disp(conj(phi2_corrected))
% phi2_discrete= real(phi2_corrected);
phi2_discrete= (phi2_corrected);
% writematrix(phi1_corrected,'phi1_corrected.dat','Delimiter','tab');

phi1_discrete = phi1_corrected;
Psi1_discrete = psi1_corrected;
% psi2_discrete = real(psi2_corrected);
psi2_discrete = (psi2_corrected);
phi0 = phi0_corrected;

% 
% disp(phi0);
% disp(phi1_discrete);
% disp(phi2_discrete);
% % New Attempt:
% % --- recompute raw alphas once ---
% [alpha1_new, alpha2_new] = intrinsicAlphaNormalization( ...
%     fem, keepDofs, phi1_discrete, phi2_discrete, Mred, Kred, Nm, ...
%     r_mid, chain, elemChainL, Psi1_discrete, psi2_discrete, phi1_half);
% 
% 
% fprintf('cond(alpha1_new)=%.2e, cond(alpha2_new)=%.2e\n', cond(alpha1_new), cond(alpha2_new));
% fprintf('‖skew(alpha1_new)‖F=%.2e, ‖skew(alpha2_new)‖F=%.2e\n', ...
%         norm(alpha1_new - alpha1_new.','fro'), norm(alpha2_new - alpha2_new.','fro'));
% disp('diag(alpha1_new):'); disp(diag(alpha1_new).');
% disp('diagalpha2_newA2):'); disp(diag(alpha2_new).');
% 
% 
% % --- symmetrize ---
% A1 = 0.5*(alpha1_new + alpha1_new.');
% A2 = 0.5*(alpha2_new + alpha2_new.');
% 
% fprintf('cond(A1)=%.2e, cond(A2)=%.2e\n', cond(A1), cond(A2));
% fprintf('‖skew(A1)‖F=%.2e, ‖skew(A2)‖F=%.2e\n', ...
%         norm(A1 - A1.','fro'), norm(A2 - A2.','fro'));
% disp('diag(A1):'); disp(diag(A1).');
% disp('diag(A2):'); disp(diag(A2).');
% 
% 
% % Optional but recommended: flip column signs to make diag positive
% s1 = sign(diag(A1)); s1(s1==0)=1; S1=diag(s1); A1=S1*A1*S1; phi1_discrete=phi1_discrete*S1; Psi1_discrete=Psi1_discrete*S1;
% s2 = sign(diag(A2)); s2(s2==0)=1; S2=diag(s2); A2=S2*A2*S2; phi2_discrete=phi2_discrete*S2;  psi2_discrete=psi2_discrete*S2;
% 
% % Full whitening with ridge clamp (preserves dimension)
% eigcut = 1e-10;
% 
% [V1,D1] = eig(A1); lam1 = max(real(diag(D1)),0);
% lam1c = max(lam1, eigcut*max(lam1(lam1>0)));    % clamp tiny/neg
% T1_full = V1 * diag(1./sqrt(lam1c)) * V1.';     % 20x20
% 
% [V2,D2] = eig(A2); lam2 = max(real(diag(D2)),0);
% lam2c = max(lam2, eigcut*max(lam2(lam2>0)));
% T2_full = V2 * diag(1./sqrt(lam2c)) * V2.';     % 20x20
% 
% % Apply on columns (keep all arrays 96x20)
% phi1_discrete = phi1_discrete * T1_full;
% Psi1_discrete = Psi1_discrete * T1_full;
% phi1_half     = phi1_half     * T1_full;
% 
% phi2_discrete = phi2_discrete * T2_full;
% psi2_discrete = psi2_discrete * T2_full;

% % --- sign-fix so both are SPD-ish ---
% s1 = sign(diag(A1)); s1(s1==0)=1; S1 = diag(s1);  A1 = S1*A1*S1;  phi1_discrete = phi1_discrete*S1;  Psi1_discrete = Psi1_discrete*S1;
% s2 = sign(diag(A2)); s2(s2==0)=1; S2 = diag(s2);  A2 = S2*A2*S2;  phi2_discrete = phi2_discrete*S2;  psi2_discrete = psi2_discrete*S2;
% 
% % --- eigen-whitening with cutoff ---
% eigcut = 1e-10;  % relative cutoff
% % α1 -> I
% [V1,D1] = eig(A1); lam1 = max(diag(D1), 0);
% tol1 = eigcut*max(lam1);
% keep1 = lam1 > tol1;
% T1 = V1(:,keep1) * diag(1./sqrt(lam1(keep1)));
% 
% phi1_discrete  = phi1_discrete * T1;
% Psi1_discrete  = Psi1_discrete * T1;
% phi1_half      = phi1_half     * T1;   % keep φ1_half consistent if 
% 
% % α2 -> I  (because psi2 already has 1/ω)
% [V2,D2] = eig(A2); lam2 = max(diag(D2), 0);
% tol2 = eigcut*max(lam2);
% keep2 = lam2 > tol2;
% T2 = V2(:,keep2) * diag(1./sqrt(lam2(keep2)));
% 
% phi2_discrete  = phi2_discrete * T2;
% psi2_discrete  = psi2_discrete * T2;
% 
% 
% disp(size(phi1_discrete));
% disp(size(Psi1_discrete));
% 
% disp(size(phi2_discrete));
% disp(size(psi2_discrete));
% check alpha
[alpha1_check, alpha2_check] = intrinsicAlphaNormalization( ...
             fem, keepDofs, phi1_discrete, phi2_discrete, Mred, Kred, nModes, ...
             r_mid, chain, elemChainL, Psi1_discrete, psi2_discrete, phi1_half);

Inm = diag(eye(Nm));

% disp('alpha1_check'); disp(alpha1_check)
% 
% 
% disp('alpha2_check'); disp(alpha2_check)


% if diag(alpha1_check) ~= Inm
% disp(norm(alpha1_check))
% if norm(alpha1_check - eye(Nm),'fro') > 1e-8
if diag(alpha1_check) - diag(eye(Nm)) > 1e-8
    error('Alpha 1 diagonal is not Identity')
end


% if diag(alpha2_check) ~= Inm
% if norm(alpha2_check - eye(Nm),'fro') > 1e-8
if diag(alpha2_check) - diag(eye(Nm)) > 1e-8
% if norm(alpha2_check - diag(w0),'fro') > 1e-10
    error('Alpha 2 diagonal is not Identity')
end
%% Sigma Test
cfg.struct.damping_model = 'kv_sigma';
cfg.struct.zeta = 0.005;          % e.g. 0.5% critical across all 6 dofs
Sigma = rom_build_sigma_kv(phi1_discrete, chain, elemChainL, Eall, cfg);

% store for your ROM builder
ModeVars_discrete.Sigma = Sigma;
% disp(Sigma)
%%

% Node A phi modes
sL = fem.coordinates(2,2);
sR = fem.coordinates(end,2);
sA = fem.coordinates(1,2);
eLen = sR-sL;
phi1L = phi1_discrete(1:6, :);
phi1R = phi1_discrete(end-5:end, :);

% phi1L = phi1_discrete(43:48, :);
% phi1R = phi1_discrete(49:54, :);
% phi1_sA = (phi1L*(sR - sA)/elemLen(1)+ phi1R*(sA-sL)/elemLen(1));
phi1_sA = (phi1L*(sR - sA)/eLen+ phi1R*(sA-sL)/eLen);

phi2L = phi2_discrete(1:6, :);
% Rotation matrix to switch x and y axes

phi2R = phi2_discrete(end-5:end, :);
phi2_sA = (phi2L*(sR - sA)/eLen+ phi2R*(sA-sL)/eLen);
%% Eta_e Calculations
fext_nodes = fem.app_forces;
% fext_nodes(3,[9 10]) = -100;

% for i = 1:nNodes
%     thisF = fem.app_forces(1:3, i);
%     thisM = fem.app_forces(4:6, i);
%     dofRange = double((i-1)*6) + (1:6)
%     % app_forces(1:6, n) => a 6x1 giving [Fx, Fy, Fz, Mx, My, Mz]
%     fext_nodes(dofRange) = fem.app_forces(:, i);
%     fext_nodes(3,:) = -100;
% end
etaVec = projectExternalForces(phi1_discrete, fext_nodes, keepDofs);

% % 6) Build the “Gamma1, Gamma2” from an “intrinsic” perspective
disp('>>> Building gamma1,gamma2 from oriented cross-sections + compliance...');
%% Gamma Calculations
doCubic = 0; % Need to setup Cubic Interpol.
[Gamma1, Gamma2] = buildNonlinearGammas_intrinsic(elemChainL, fem, keepDofs, phi1_discrete, phi2_discrete, Mred, Kred, ...
    Psi1_discrete,psi2_discrete, chain, phi1_half, doCubic);



    for sIdx = 1:chain_id
        if sIdx ~= chain_id
            chain_s = chain(chain_start(sIdx):chain_end(sIdx));
        else % only for last chain
            chain_s = chain(chain_start(sIdx):end);
        end
        nChain = length(chain_s);
        elemLen = zeros(1, nChain-1);
        for i=1:(nChain-1)
            nA = chain_s(i);
            nB = chain_s(i+1);
            elemLen(i)= norm( XYZ(nB,:) - XYZ(nA,:) );
        end     
        % disp(elemLen)
        %% 2) Build Phi0 => piecewise polynomials
        for j=1:nModes
           nodeVals_j = zeros(6, nChain);
           for i=1:nChain
               nID= chain_s(i);
               dofRange= (nID-1)*6 + (1:6);
               idx_k = localDofIndex(keepDofs, dofRange);
               if isempty(idx_k), continue; end
               nodeVals_j(:,i)= phi0(idx_k,j);
           end
           seg_tmp = buildPiecewiseLin(nodeVals_j, elemLen, chain_s, XYZ);
           % seg(seg_ids:seg_ide) = seg_tmp;
           Phi0Fun{j}.seg(seg_ids:seg_ide) = seg_tmp;
        end
        % seg_ids = seg_ids+ seg_ide;
        % seg_ide = 2*seg_ide;
    % end  
        %% 3) Build Phi1 => piecewise polynomials
        for j=1:nModes
           nodeVals_j = zeros(6,nChain);
           for i=1:nChain
               nID= chain_s(i);
               dofRange= (nID-1)*6 + (1:6);
               idx_k = localDofIndex(keepDofs, dofRange);
               if isempty(idx_k), continue; end
               nodeVals_j(:,i)= phi1_discrete(idx_k,j);
           end
           seg_tmp2 = buildPiecewiseLin(nodeVals_j, elemLen, chain_s, XYZ);
           seg2(seg_ids:seg_ide) = seg_tmp2;
           Phi1Fun{j}.seg(seg_ids:seg_ide) = seg_tmp2;
        end
        seg_ids = seg_ids+ seg_ide;
        seg_ide = 2*seg_ide;

    end







%% Test Piecewise and Eval Funct
% Mode_plot_test = 2; % test mode to plot
for Mode_plot_test = 1:5
Nseg = 20; % number of intermediate pts. (min is 2 for sA and sB)
% Ntot = Nseg*(activePts)+activePts-2;
% sdasd = fem.coordinates;
% L_sub =  linspace(idxMinL, idxMaxL, activePts+1);
PhiVal = zeros(6*Nseg*activePts,1);
wing_s = []; 
Vx = []; Fx = [];
Vy = []; Fy = [];
Vz = []; Fz = [];
Wx = []; Mx = [];
Wy = []; My = [];
Wz= []; Mz = [];
    for iSeg = 1:activePts
        Phi1_mod = Phi1Fun{Mode_plot_test}.seg;
        srange_mod = Phi1_mod(1,iSeg);
        s_range = srange_mod.xrange;
        PhiVal_i = srange_mod.coefs;
        s0 = s_range(1);
        s1 = s_range(2);
        L_sub = linspace(s0,s1, Nseg);
        iSegsub = double((iSeg-1)*6*Nseg);
        idx_R = double((iSeg-1)*3)+(1:3);
        idx_p = double((iSeg-1))*6+(1:6);
        Rglobal0_s = Rglobal0(idx_R,:);
        for s = 1:length(L_sub)
            % disp('s'); disp(s)
            s_idx = double((s-1)*6)+iSegsub+(1:6);
            sLoc = L_sub(s);
            if sLoc == s0  
                PhiIter = PhiVal_i(:,1);
                PhiVal(s_idx,:) =PhiIter;
            elseif sLoc== s1
                PhiIter = PhiVal_i(:,2);
                PhiVal(s_idx,:) =PhiIter;
            else 
                PhiIter = Phi_eval(Rglobal0_s, s_range, PhiVal_i,sLoc);
                PhiVal(s_idx,:) =PhiIter;
            end
            r_mid_s = r_mid(iSeg);
            R6 = blkdiag(Rglobal0_s,Rglobal0_s);
            Phi2Iter = R6*phi2_discrete(idx_p,Mode_plot_test);
            Phi2Val(s_idx,:) = Phi2Iter;

            wing_s(end+1) = sLoc;
            Vx(end+1) = PhiIter(1); 
            Vy(end+1) = PhiIter(2);
            Vz(end+1) = PhiIter(3);
            Wx(end+1) = PhiIter(4);
            Wy(end+1) = PhiIter(5);
            Wz(end+1) = PhiIter(6);

            Fx(end+1) = Phi2Iter(1);
            Fy(end+1) = Phi2Iter(2);
            Fz(end+1)= Phi2Iter(3);
            Mx(end+1) =Phi2Iter(4);
            My(end+1) =Phi2Iter(5);
            Mz(end+1) = Phi2Iter(6);
        end
    end
    [wing_s, sortIdx] = sort(wing_s);
    for i = 1:length(wing_s)
        Vxp(i) = Vx(sortIdx(i));
        Vyp(i) = Vy(sortIdx(i));
        Vzp(i) = Vz(sortIdx(i));
        Wxp(i) = Wx(sortIdx(i));
        Wyp(i) = Wy(sortIdx(i));
        Wzp(i) = Wz(sortIdx(i));

        Fxp(i) = Fx(sortIdx(i));
        Fyp(i) =Fy(sortIdx(i));
        Fzp(i)= Fz(sortIdx(i));
        Mxp(i) = Mx(sortIdx(i));
        Myp(i) =My(sortIdx(i));
        Mzp(i) = Mz(sortIdx(i));
    end

    if dubugPlot ==1
    figure;
    subplot(2,1,1);
        plot(wing_s, Vxp,'r','LineWidth',1.2,'DisplayName','Vx');
        hold on;
        plot(wing_s, Vyp,'g','LineWidth',1.2,'DisplayName','Vy');
        plot(wing_s, Vzp,'b','LineWidth',1.2,'DisplayName','Vz');
        xlabel('sub-element index i');
        ylabel('Vel amplitude');
        title(sprintf('\\phi_1 Cont Full Winged Vel Components (Mode %d)', Mode_plot_test));
        legend('Location','best'); 
        grid on;


       subplot(2,1,2);
       plot(wing_s, Wxp,'r','LineWidth',1.2,'DisplayName','Wx');
       hold on;
       plot(wing_s, Wyp,'g','LineWidth',1.2,'DisplayName','Wy');
       plot(wing_s, Wzp,'b','LineWidth',1.2,'DisplayName','Wz');
       xlabel('sub-element index i');
       ylabel('Ang Vel amplitude');
       title(sprintf('\\phi_1 Ang Vel Components (Mode %d)', Mode_plot_test));
       legend('Location','best');
       grid on;

    figure;
    subplot(2,1,1);
       plot(wing_s, Fxp,'r','LineWidth',1.2,'DisplayName','Fx');
       hold on;
       plot(wing_s, Fyp,'g','LineWidth',1.2,'DisplayName','Fy');
       plot(wing_s, Fzp,'b','LineWidth',1.2,'DisplayName','Fz');
        xlabel('sub-element index i');
        ylabel('Vel amplitude');
        title(sprintf('\\phi_2 Full Winged Force Components (Mode %d)', Mode_plot_test));
        legend('Location','best'); 
        grid on;


       subplot(2,1,2);
       plot(wing_s, Mxp,'r','LineWidth',1.2,'DisplayName','Mx');
       hold on;
       plot(wing_s, Myp,'g','LineWidth',1.2,'DisplayName','My');
       plot(wing_s, Mzp,'b','LineWidth',1.2,'DisplayName','Mz');
       xlabel('sub-element index i');
       ylabel('Ang Vel amplitude');
       title(sprintf('\\phi_2 Moment Components (Mode %d)', Mode_plot_test));
       legend('Location','best');
       grid on;
    end
end
% Esteban:
    % Linear Int
    % % test
    % for l = 1:Nm
    %     jphi0_a = phi0(node_idx, l);
    %     jphi0_b = phi0(node_idx_next, l);
    % 
    %     for k = 1:activePts-1
    %         node_idx = double((k-1)*6)+(1:6);
    %         node_idx_next = double((k)*6)+(1:6);
    %         Rs_idx2 = double((k-1)*3)+(1:3);
    %         Rot_mat_global0_s = Rglobal0(Rs_idx2,:);
    %         sA = fem.coordinates(k+1, 2);
    %         sB = fem.coordinates(k+2, 2);
    %         sRange_i = struct('sA', sA, 'sB', sB);
    %         nSamplesPerSeg = 20;
    %         phi0_Fun = @(s) phi0toCont(Rot_mat_global0_s, sRange_i, jphi0_a, jphi0_b,s);
    %         for iS=0:nSamplesPerSeg
    %            xi = iS / nSamplesPerSeg;   % 0..1
    %            sLoc= ds* xi;              % local coordinate in [0..ds]
    %            % sLocD = sLoc - Acoord(2);
    %            phiVal = evalCubic6(myFun.seg(iSeg).coefs, sLoc);
    % 
    %            % The 2D or 3D location in space => param = Acoord + segVec*(xi).
    %            thisPoint = Acoord + segVec*xi;   % e.g. linear
    %            yVal= thisPoint(2);  % pick the Y coord
    % 
    %            Fy  = phiVal(2);     % for example, the second component => velocity in y
    %            allVal(end+1,:) = thisPoint;
    % 
    %            phiAll(SegIdx,iS+1) = phiVal;
    % 
    %            allY(end+1)= yVal;
    %            allFy(end+1)= Fy;
    %            allFx(end+1) =  phiVal(1);
    %            allFz(end+1) = phiVal(3);
    %            allMx(end+1) =  phiVal(4);
    %            allMy(end+1) = phiVal(5);
    %            allMz(end+1) =  phiVal(6);
    %         end
    %     end
    % end


% 
% 
% function phi0_Fun = phi0toCont(Rot_mat_global0_s, sRange_i, jphi0_a, jphi0_b, s)
%     sA = sRange_i.sA ;
%     sB= sRange_i.sB;
%     ds = sB-sA;
%     s_subleft = (sB-s)/ds;
%     s_subright = (s-sA)/ds;
%     phi0_Fun = Rot_mat_global0_s*(jphi0_a*s_subleft+jphi0_b*s_subright);
% end





%% Left off here



ModeVars_discrete.phi0_local = phi0;
ModeVars_discrete.phi1_local = phi1_discrete;
ModeVars_discrete.phi2_local = phi2_discrete;
ModeVars_discrete.etaVecExt = etaVec;
ModeVars_discrete.Gamma1 = Gamma1;
ModeVars_discrete.Gamma2 = Gamma2;
ModeVars_discrete.psi1_local = Psi1_discrete;
ModeVars_discrete.psi2_local = psi2_discrete;

ModeVars_continuous.phi0 = Phi0Fun;
ModeVars_continuous.phi1 = Phi1Fun;
% ModeVars_continuous.phi2 = [phi2FunTop; phi2FunBot];
% ModeVars_continuous.psi1 = [Psi1FunTop; Psi1FunBot];
% ModeVars_continuous.psi2 = [Psi2FunTop; Psi1FunBot];

Beam_Props.M = Mred;
Beam_Props.K = Kred;
Beam_Props.E = Eall;

% fem.keepDofs = keepDofs;
% [Xfull_time, orientation_time] = recoverFull3Dsolution( ...
%         fem, T_oa, phi0, q0Hist, dt, nodeCenters, warpingMatrices );
% debugPlotIntrinsicModes3D( Xfull_time, iTime, fem )


disp('All done.');


end % main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


function phi0_Fun = phi0toCont(Rot_mat_global0_s, sRange_i, jphi0_a, jphi0_b, s)
    sA = sRange_i.sA ;
    sB= sRange_i.sB;
    ds = sB-sA;
    s_subleft = (sB-s)/ds;
    s_subright = (s-sA)/ds;
    phi0_Fun = Rot_mat_global0_s*(jphi0_a*s_subleft+jphi0_b*s_subright);
end
function phi_val= Phi_eval(Rglobal0_s, s_range, PhiVal_i,s)
    sA = s_range(1) ;
    sB= s_range(2);
    ds = sB-sA;
    jphi0_a = PhiVal_i(:,1);
    jphi0_b = PhiVal_i(:,2);
    s_subleft = (sB-s)/ds;
    s_subright = (s-sA)/ds;
    R6 = blkdiag(Rglobal0_s,Rglobal0_s);
    phi_val = R6*(jphi0_a*s_subleft+jphi0_b*s_subright);
end
function phi_val= Phi2_eval(Rglobal0_s, s_range, PhiVal_i,s)
    sA = s_range(1) ;
    sB= s_range(2);
    ds = sB-sA;
    jphi0_a = PhiVal_i(:,1);
    jphi0_b = PhiVal_i(:,2);
    s_subleft = (sB-s)/ds;
    s_subright = (s-sA)/ds;
    R6 = blkdiag(Rglobal0_s,Rglobal0_s);
    phi_val = R6*(jphi0_a*s_subleft+jphi0_b*s_subright);
end

function R3 = buildElementRotation_sub(fem, eId, subId)
% BUILDELEMENTROTATION_SUB   build 3x3 local rotation for the subinterval
% If subId=1 => node1->node2, subId=2 => node2->node3
% We'll interpret fem.frame_of_reference_delta(:,:, eId) as the base orientation 
% and possibly add fem.structural_twist(eId, subId).

  % d = fem.frame_of_reference_delta(:,:, eId);  % e.g. 3x3 or 3-vector
  % % If it's 3x3, we assume it is a direct orientation. 
  % % Then do an extra structural twist for subId if you want.
  % 
  % R3 = d;  % if it is already 3x3
  % % Possibly do a small e.g. 
  R3 = eye(3);
  twistVal = fem.structural_twist(subId,eId); % or something
  if abs(twistVal)>1e-12
     axLocal= R3(:,1); % local x
     Rt = rotAxis(axLocal, twistVal);
     R3 = R3 * Rt; 
  end
end

function R = rotAxis(ax, theta)
%rotAxis builds 3×3 rotation about axis `ax` by angle `theta`
   ax = ax(:) / (norm(ax)+1e-16);
   ct = cos(theta);  st = sin(theta);
   omct = 1.0 - ct;
   x = ax(1);  y = ax(2);  z = ax(3);
   R = [ ct + x^2*omct,    x*y*omct - z*st,   x*z*omct + y*st
         y*x*omct + z*st,  ct + y^2*omct,     y*z*omct - x*st
         z*x*omct - y*st,  z*y*omct + x*st,   ct + z^2*omct ];
end







function [Gamma1, Gamma2, psi1_local, psi1_global] = buildNonlinearGammas_intrinsic(elemChainL, fem, keepDofs, phi1, phi2, Mred, Kred, psi1, ...
    psi2, chain, phi1_half, doCubic)
%  Improved approach:
%    Gamma1(j,k,l) = ∫ [ phi1_j^T * L1(phi1_k) * M phi1_l ] ds
%    Gamma2(j,k,l) = ∫ [ phi1_j^T * L2(phi2_k) *  C phi2_l ] ds
%  We'll do an element-by-element sum. For each 3-node element, we compute
%  or have the local cross-section Ssec => Csec=inv(Ssec). Then do
%  shape function with 2 sub-2node => trapezoid integration.
%
nm = size(phi1,2);
Gamma1 = zeros(nm,nm,nm);
Gamma2 = zeros(nm,nm,nm);

L1fun = @(x) L1_operator(x);
L2fun = @(x) L2_operator(x);

ne = fem.num_elem;


ndofs    = size(phi1,1);
% Cred = inv(Kred);

Nm = size(phi1, 2);             % number of modes
Nn = length(chain);            % number of nodes
Ne = Nn - 1;                    % number of elements

Cred = Kred\eye(size(Kred));

for i = 1:Nn
    chain_s = chain(i);
    idx_i     = double((chain_s-1)*6) + (1:6);
    % idx_ip1   = double(chain_s*6) + (1:6);   % node i+1
    % idx_mid   = idx_i;        % midpoint assumed between i and i+1

    idx_k = localDofIndex(keepDofs, idx_i);
    % idx_k2 = localDofIndex(keepDofs, idx_ip1);

    if isempty(idx_k)
            % phi1_i = zeros(6,1);
            % phi1_i_plus = phi1(idx_k2,jInd);
            % phi1_i_half(idx_i,jInd) = (phi1_i+phi1_i_plus)/2;
            % disp('skipped')
            continue; 

    % elseif isempty(idx_k2)
    %         % phi1_i = phi1(idx_k,jInd);
    %         % phi1_i_plus = zeros(6,1);
    %         % phi1_i_half(idx_i,jInd) = (phi1_i+phi1_i_plus)/2;
    %         continue; 
    % else
            % phi1_i = phi1(idx_k,jInd);
            % phi1_i_plus = phi1(idx_k2,jInd);
            % phi1_i_half(idx_i,jInd) = (phi1_i+phi1_i_plus)/2;

    end
    elem_s = elemChainL(i-1);

    for jInd=1:Nm
         for kInd=1:Nm
           for lInd=1:Nm
                phi1_j = phi1(idx_k,jInd);
                phi1_k = phi1(idx_k,kInd);
                psi1_l = psi1(idx_k, lInd);

                gamma1_local = phi1_j.' * L1fun(phi1_k) * (psi1_l);
                % gamma1_local = phi1_j.' * L1fun(phi1_k) * (psi1_l)*elem_s;
                Gamma1(jInd,kInd,lInd) = Gamma1(jInd,kInd,lInd) + gamma1_local;

                phi1_i_half =  phi1_half(idx_k, jInd);
                if doCubic == 1
                else
                    phi1_j_half = phi1_i_half;
                end
                phi2_k = phi2(idx_k,kInd);
                psi2_l = psi2(idx_k, lInd);

                gamma2_local= phi1_j_half.' * L2fun( phi2_k ) * psi2_l*elem_s;
                % gamma2_local= phi1_j_half.' * L2fun( phi2_k ) * psi2_l;
                Gamma2(jInd,kInd,lInd) = Gamma2(jInd,kInd,lInd) + gamma2_local;
           end
         end
    end
end



end


function [phi1_loc, phi2_loc] = transformModesToLocal( ...
           fem, keepDofs, phi1_global, phi2_global )
% transformModesToLocal
% Inputs:
%   fem            = your FE structure with .num_elem, .connectivities, etc.
%   keepDofs       = list of active dofs after condensation
%   phi1_global    = [ndofs x nModes], velocity-like modes in global coords
%   phi2_global    = [ndofs x nModes], force-like modes in global coords
%
% Outputs:
%   phi1_loc       = same dimension [ndofs x nModes], but each sub-block is
%                    expressed in local frames. 
%   phi2_loc       = likewise in local frames.


ndofs    = size(phi1_global,1);
nModes   = size(phi1_global,2);
phi1_loc = zeros(ndofs, nModes);
phi2_loc = zeros(ndofs, nModes);

ne = fem.num_elem;

for e = 1 : ne
    conn = fem.connectivities(e,:);
    dofA = (conn(1)-1)*6 + (1:6);
    dofB = (conn(3)-1)*6 + (1:6);
    dofC = (conn(2)-1)*6 + (1:6);

    subList = {[dofA,dofB], [dofB,dofC]};
    % subList = { [conn(1), conn(3)], [conn(3), conn(2)] };
    % We'll define sub-length ds=1 for each sub. Or we can use actual node distance:
    % ds ~ norm( xyzB-xyzA ), etc. We'll do actual:
    for subId=1:2
       dm = subList{subId};
       localInd = [];
       for dd=dm
           idxLocal = find(keepDofs==dd);
           if ~isempty(idxLocal)
               localInd=[localInd, idxLocal];
           end
       end
       % If fully present => 12 dofs. If clamp => partial
       if length(localInd)<12 
           continue; 
       end
       leftI  = localInd(1:6);
       rightI = localInd(7:12);

       R3_sub  = buildElementRotation_sub(fem, e, subId);
       R6_sub  = blkdiag(R3_sub, R3_sub);



       % transform each mode j:
       for j = 1:nModes
          % disp('mode num'); disp(j);
          x1_L = phi1_global(leftI, j);
          x1_R = phi1_global(rightI,j);
          x2_L = phi2_global(leftI, j);
          x2_R = phi2_global(rightI,j);

          x1_Lloc = R6_sub * x1_L;
          x1_Rloc = R6_sub * x1_R;
          x2_Lloc = R6_sub * x2_L;
          x2_Rloc = R6_sub * x2_R;

          % store back in phi1_loc, phi2_loc
          phi1_loc(leftI,j)  = x1_Lloc;
          phi1_loc(rightI,j) = x1_Rloc;
          phi2_loc(leftI,j)  = x2_Lloc;
          phi2_loc(rightI,j) = x2_Rloc;
          % disp(phi1_loc([leftI,rightI], j));
          sgrsgrs = norm(phi1_loc([leftI,rightI], j)- phi1_global([leftI,rightI], j));
          if sgrsgrs ~= 0
            % disp('not like us')
          end
       end
    end
end

end

function S = skew(u)
S= [ 0 -u(3) u(2);
     u(3) 0 -u(1);
    -u(2) u(1) 0 ];
end

function M6 = L1_operator(x6)
% x6 = [v(3); w(3)]
v = x6(1:3);
w = x6(4:6);
M6= zeros(6);
%   top-left = skew(-w), bottom-left=skew(-v), bottom-right= skew(-w)
% Why neg?
% M6(1:3,1:3)= skew(-w);
% M6(4:6,1:3)= skew(-v);
% M6(4:6,4:6)= skew(-w);
M6(1:3,1:3)= skew(w);
M6(4:6,1:3)= skew(v);
M6(4:6,4:6)= skew(w);
end

function M6 = L2_operator(x6)
% x6= [f(3); m(3)]
f = x6(1:3);
m = x6(4:6);
M6= zeros(6);
% top-right= skew(f), bottom-right= skew(m), symmetrical bottom-left= skew(f)
M6(1:3,4:6)= skew(f);
M6(4:6,4:6)= skew(m);
M6(4:6,1:3)= skew(f);
end

function etaVec = projectExternalForces(phi1, fext_nodes, keepDofs)
% phi1 = [nActiveDofs x nModes] of velocity modes
% fext_nodes = cell or array of node-based forces, e.g. fext_nodes{i} = [Fx,Fy,Fz,Mx,My,Mz]
% keepDofs = the list of condensed dofs
%
% Output:
%  etaVec = [nModes x 1], i.e. the total projection. 
%           If you want a time series, you can do it. 
%

nModes = size(phi1,2);
etaVec = zeros(nModes,1);

for iNode = 1:length(fext_nodes)
    % if node iNode has a 6-vector of external force+moment:
    F_i = fext_nodes(:,iNode);   % e.g. [Fx, Fy, Fz, Mx, My, Mz]
    if norm(F_i)==0
        continue;
    end

    dofRange = (iNode-1)*6 + (1:6);
    iLoc = findIndex(keepDofs, dofRange);
    if isempty(iLoc)
        continue; % means that node is omitted
    end
    for j = 1:nModes
       etaNew = phi1(iLoc,j)' * F_i(:);
       etaVec(j) = etaVec(j) + etaNew;
    end
end

end



function iLoc = findIndex(keepDofs, dofRange)
% FINDINDEX  Utility that finds, for each integer in dofRange,
%            where it appears in keepDofs. 
%
%  iLoc = findIndex(keepDofs, dofRange) 
%  returns a row vector of indices into keepDofs.

iLoc = zeros(size(dofRange));  % Pre-allocate
count = 0;
for ii = 1:length(dofRange)
    d = dofRange(ii);
    loc = find(keepDofs == d, 1, 'first');
    if ~isempty(loc)
        count = count + 1;
        iLoc(count) = loc;
    end
end
% Trim trailing zeros if some dofRange entry wasn’t found
iLoc = iLoc(1:count);
end



function [alpha1, alpha2] = intrinsicAlphaNormalization( ...
    fem, keepDofs, phi1_loc, phi2_loc, Mred, Kred, nModes, r_mid, chain, ...
    elemChainL, Psi1_discrete, psi2_discrete, phi1_half)

% We'll build alpha1, alpha2 => [nModes x nModes]

% elemLen_sub1 = elemLen(:,1);
% elemLen_sub2 = elemLen(:,2);
% elemLen = [elemLen_sub1;elemLen_sub2]


Nm = size(phi1_loc, 2);             % number of modes
Nn = length(chain);            % number of nodes
Ne = Nn - 1;                    % number of elements

Cred = Kred\eye(size(Kred));

alpha1 = zeros(Nm, Nm);
alpha2 = zeros(Nm, Nm);
phi1 = phi1_loc;
    for i = 1:Nn
        chain_s = chain(i);

        for j = 1:Nm
            for k = 1:Nm

                % Initialize accumulators
                A1_jk = 0;
                A2_jk = 0;
                % Indices for node i and midpoint i+1/2
                idx_i     = double((chain_s-1)*6) + (1:6);
                idx_ip1   = double(chain_s*6) + (1:6);   % node i+1
                % idx_mid   = idx_i;        % midpoint assumed between i and i+1


                idx_k = localDofIndex(keepDofs, idx_i);
                idx_k2 = localDofIndex(keepDofs, idx_ip1);

                if isempty(idx_k)

                    continue; 

                elseif isempty(idx_k2)

                    continue; 
                else
                    % elem_Ls = elemChainL(i-1);
                    elem_Ls = elemChainL(i);
                    % ii = max(1, min(i, numel(elemChainL)));
                    % elem_Ls = elemChainL(ii);


                    phi1_j_i = phi1_loc(idx_k, j); 
                    Psi1_k_i = Psi1_discrete(idx_k, k);
                    % First add local self-coupling
                    % A1_jk =  alpha1(j,k) + phi1_j_i' * Psi1_k_i;
                    % alpha1(j,k) = alpha1(j,k) + phi1_j_i' * Psi1_k_i;
                    alpha1(j,k) = alpha1(j,k) + phi1_j_i.' * Psi1_k_i;
                    % alpha1(j,k) = alpha1(j,k) + elem_Ls*phi1_j_i.' * Psi1_k_i;

                 

                    phi2_j_i = phi2_loc(idx_k, j);
                    Psi2_k_i = psi2_discrete(idx_k, k);
                    % testing 
                    phi1half_k_i = phi1_half(idx_k, k); % A calls for C*Phi1half 
                    phi2_k_i = phi2_loc(idx_k, k); % so CPhi2 half * Phi1half /Phi2 half
                    % test_var =  Psi2_k_i/phi2_k_i*phi1half_k_i;

                    alpha2(j,k) = alpha2(j,k) + elem_Ls*phi2_j_i.' * Psi2_k_i;
                    % alpha2(j,k) = alpha2(j,k) + phi2_j_i.' * Psi2_k_i;
                    % alpha2(j,k) = alpha2(j,k) + elem_Ls*phi2_j_i.' * test_var;

                end            


            end
        end
    end


end


function [alpha1, alpha2] = intrinsicAlphaNormalizationScalar( ...
    fem, keepDofs, phi1_loc, phi2_loc, Mred, Kred, nModes, r_mid, chain, ...
    elemChainL, Psi1_discrete, psi2_discrete, phi1_half)

% build alpha1, alpha2 => [nModes x nModes]

% elemLen_sub1 = elemLen(:,1);
% elemLen_sub2 = elemLen(:,2);
% elemLen = [elemLen_sub1;elemLen_sub2]


Nm = size(phi1_loc, 2);             % number of modes
Nn = length(chain);           % number of nodes
Ne = Nn - 1;                    % number of elements
% disp(chain)
Cred = Kred\eye(size(Kred));

alpha1 = zeros(Nm, 1);
alpha2 = zeros(Nm, 1);
phi1 = phi1_loc;
    for i = 1:Nn
        chain_s = chain(i);

        % for j = 1:Nm
            for k = 1:Nm

                % Initialize accumulators
                A1_jk = 0;
                A2_jk = 0;
                % Indices for node i and midpoint i+1/2
                idx_i     = double((chain_s-1)*6) + (1:6);
                idx_ip1   = double(chain_s*6) + (1:6);   % node i+1
                % idx_mid   = idx_i;        % midpoint assumed between i and i+1


                idx_k = localDofIndex(keepDofs, idx_i);
                idx_k2 = localDofIndex(keepDofs, idx_ip1);

                if isempty(idx_k)

                    continue; 

                elseif isempty(idx_k2)

                    continue; 
                else
                    % elem_Ls = elemChainL(i-1);
                    elem_Ls = elemChainL(i);
                    % ii = max(1, min(i, numel(elemChainL)));
                    % elem_Ls = elemChainL(ii);


                    phi1_j_i = phi1_loc(idx_k, k); 
                    Psi1_k_i = Psi1_discrete(idx_k, k);
                    % First add local self-coupling
                    % A1_jk =  alpha1(j,k) + phi1_j_i' * Psi1_k_i;
                    % alpha1(j,k) = alpha1(j,k) + phi1_j_i' * Psi1_k_i;
                    alpha1(k) = alpha1(k) + phi1_j_i.' * Psi1_k_i;
                    % alpha1(j,k) = alpha1(j,k) + elem_Ls*phi1_j_i.' * Psi1_k_i;

                 

                    phi2_j_i = phi2_loc(idx_k, k);
                    Psi2_k_i = psi2_discrete(idx_k, k);
                    
                    % testing 
                    phi1half_k_i = phi1_half(idx_k, k); % A calls for C*Phi1half 
                    phi2_k_i = phi2_loc(idx_k, k); % so CPhi2 half * Phi1half /Phi2 half
                    % test_var =  Psi2_k_i/phi2_k_i*phi1half_k_i;

                    alpha2(k) = alpha2(k) + elem_Ls*phi2_j_i.' * Psi2_k_i;
                    % alpha2(j,k) = alpha2(j,k) + phi2_j_i.' * Psi2_k_i;
                    % alpha2(j,k) = alpha2(j,k) + elem_Ls*phi2_j_i.' * test_var;

                end            


            end
        % end
    end


end


function [phi2, rMid, psi2, phi1_ihalf] = accumulatePhi2_singleChain( ...
    fem, keepDofs, Kred, phi0_global, chain, phi1_discrete, w0, elemLen, Emat, doOmeg, Mred)
%   Implementation of eq. (5.26) from Wang, for a single chain array
% 
%   We do piecewise-constant internal forces at each sub-mid i+1/2.
%   Summation of K*phi0 for all k>i (or k> i? depends on sign).
%
% if doOmeg ==1
%  disp('Doing Omeg')
% end
if chain(1) ~= 1
    chain = flip(chain)
end


numModes = size(phi0_global,2);
ndof     = size(phi0_global,1);
rMid = [];
psi2 = zeros((length(chain)-1)*6, numModes);
% Precompute fKAll
fKAll = zeros(ndof,numModes);
% for j=1:numModes
%     fKAll(:,j) = Kred*phi0_global(:,j);
% end
for j=1:numModes
    % fKAll(:,j) = -w0(j)^2*Mred*phi0_global(:,j);
    fKAll(:,j) = w0(j)^2*Mred*phi0_global(:,j);
end
% disp('size(fKAll)');
% disp(size(fKAll));

% I think you use phi0 not phi1 to calc the psi2, opposed to reference
phi1_discrete = phi0_global;

nChain = length(chain);
% Preallocate
phi2_loc = zeros(6*(nChain-1), numModes);
phi2   = zeros(6*(nChain-1), numModes);
phi1_ihalf = zeros(6*(nChain-1), numModes);
XYZ = fem.coordinates;  % [nNodes x 3], ensure row i => node i

% disp(chain)
sA = fem.coordinates(chain(1), 2);
sB = fem.coordinates(chain(end),2);
sign_ds = sign(sB - sA);          % +1 for increasing s, -1 for decreasing

for i=1:(nChain-1)
    nA = chain(i);
    nB = chain(i+1);
    % midpoint
    rMid_i = 0.5*( XYZ(nA,:) + XYZ(nB,:) );
    rMid(end+1,:) = rMid_i;
    % local rotation for sub-element i => from node nA -> nB
    R3_i = buildSegmentRotation2(fem, nA, nB, XYZ);  % see previous approach
    R6_i = blkdiag(R3_i, R3_i);

    rowIdx = (i-1)*6 + (1:6);
    elem_s = elemLen(i);
    % if i < Nn
    %     elem_s = elemLen(i);   % length of subelement [i -> i+1]
    % else
    %     continue;               % last node has no outgoing subelement
    % end
    for j=1:numModes
        forceSum_g = zeros(6,1);
        % Summation eq. (5.26) => k from i+1..end
        if chain(1) == 1 % need to change to auto but this is clamped node
            for k=(i):nChain
                nodeK = chain(k);
                dofRange = (nodeK-1)*6 + (1:6);
                idx_k = localDofIndex(keepDofs, dofRange);
                if isempty(idx_k), continue; end
                Fk_6 = fKAll(idx_k, j);

                forceSum_g(1:3) = forceSum_g(1:3) + Fk_6(1:3);
                rK = XYZ(nodeK,:);
                offset = (rK - rMid_i).';
                % offset = (rK - rMid_i)
                % ForceSec = Fk_6(1:3)
                crossTerm = cross(offset, Fk_6(1:3));
                % crossTerm = cross(offset, Fk_6(1:3).')
                forceSum_g(4:6) = forceSum_g(4:6) + Fk_6(4:6) + crossTerm;
            end 
                        % forceSum_g = -forceSum_g;

        else % Second wing start with clamped end
            for k=(i):nChain
            % for k=nChain:-1:i
                nodeK = chain(k);
                dofRange = (nodeK-1)*6 + (1:6);
                idx_k = localDofIndex(keepDofs, dofRange);
                if isempty(idx_k), continue; end
                Fk_6 = fKAll(idx_k, j);

                forceSum_g(1:3) = forceSum_g(1:3) + Fk_6(1:3);
                rK = XYZ(nodeK,:);
                offset = (rK - rMid_i).';
                crossTerm = cross(offset, Fk_6(1:3));
                forceSum_g(4:6) = forceSum_g(4:6) + Fk_6(4:6) + crossTerm;
                % forceSum_g(4:6) = forceSum_g(4:6) + Fk_6(4:6) - crossTerm;


            end
            forceSum_g = -forceSum_g;
        end

        % phi2_global_sub = -forceSum_g; 
        phi2_global_sub = forceSum_g; %this might need to be negative
        phi2(rowIdx,j)   = phi2_global_sub;
        % phi2_loc(rowIdx,j) = phi2_local_sub;

        dofRange2 = double((nA-1)*6) + (1:6);
        dofRange3 = double((nB-1)*6) + (1:6);
        idx_l = localDofIndex(keepDofs, dofRange2);
        idx_l2 = localDofIndex(keepDofs, dofRange3);
        kappa_sub = zeros(3,1);
        % if chain(1)~=1
        %     elem_s = -elem_s;
        % end
        if isfield(fem,'structural_twist')
            % per-sub estimate (twist about local x); adapt if you store per-node/element arrays:
            kappa_sub(1) = (fem.structural_twist(1,i)) / elem_s;
        end
        E6_i = L1_operator([1;0;0; kappa_sub]);   % [e1; kappa_sub]
        Emat(rowIdx,:) = E6_i;
        % Emat(rowIdx,:) = R6_i;
        if isempty(idx_l) 
            phi1_i = zeros(6,1);
            phi1_i_plus = phi1_discrete(idx_l2,j);
            % Psi2 calc
            % disp(Emat)

            psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ Emat(rowIdx,:)'*(phi1_i_plus+phi1_i)/2;
            % psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ (phi1_i_plus+phi1_i)/2;
        elseif isempty(idx_l2)
            phi1_i = phi1_discrete(idx_l,j);
            phi1_i_plus = zeros(6,1);
            % Psi2 calc
            psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ Emat(rowIdx,:)'*(phi1_i_plus+phi1_i)/2;
            % psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+(phi1_i_plus+phi1_i)/2;
            % continue; 
        else
            phi1_i = phi1_discrete(idx_l,j);
            phi1_i_plus = phi1_discrete(idx_l2,j);
            % Psi2 calc
            % psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ Emat(rowIdx,:)'*(phi1_i_plus+phi1_i)/2;
            psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ Emat(rowIdx,:)'*(phi1_i_plus+phi1_i)/2;
            % psi2(rowIdx,j) = -(phi1_i_plus-phi1_i)/elem_s+ (phi1_i_plus+phi1_i)/2;
        end
        % psi2 = -psi2;
           phi1_ihalf(rowIdx,j) =  (phi1_i_plus+phi1_i)/2;
        % if doOmeg
        %     psi2(rowIdx,j) = (1/w0(j))*psi2(rowIdx,j);
        % 
        % end
    end
end
end

function R3 = buildSegmentRotation2(fem, nodeA, nodeB, XYZ)
% BUILSEGMENTROTATION
%   Builds 3x3 local rotation from nodeA to nodeB.
%   If geometry is purely along global Y, we get a permutation matrix.

exdir = XYZ(nodeB,:) - XYZ(nodeA,:);
L = norm(exdir);
if L<1e-14
    R3=eye(3); 
    disp('skipping Segment Rot')
    return;
end
exdir= exdir(:)/L;
% pick trial ey
eyTrial = [0;1;0];
if abs(dot(exdir, eyTrial))>0.99
    eyTrial= [0;0;1];
end
ez= cross(exdir, eyTrial);
ez= ez/norm(ez);
ey= cross(ez, exdir);

R3_noTwist= [exdir, ey, ez];
% R3_noTwist = eye(3);
twVal= 0;
if isfield(fem,'structural_twist')
   % if there's an array, might do e.g. twVal= fem.structural_twist(nodeA) ...
end
if abs(twVal)>1e-12
   cT=cos(twVal); sT=sin(twVal);
   Rtw= [1  0   0;  0  cT -sT;  0  sT  cT];
   R3= R3_noTwist*Rtw;
else
   R3= R3_noTwist;
end

% check orthonormal:
if norm(R3'*R3 - eye(3),'fro')>1e-6
   warning('SegmentRotation not orthonormal?');
end
end


function R3 = buildSegmentRotation(fem, iElem, XYZ)
% buildSegmentRotation
%   Builds the 3x3 rotation from global -> local for the segment iElem->iElem+1
%   Usually the local x-axis is along the node i->i+1 direction, 
%   plus any structural twist from fem.structural_twist, etc.
% 
%   If you have 8 elements, iElem is 1..8, etc. But for a single chain, 
%   we can interpret iElem also as the "start node index".
%
%   For instance, define:
%       exdir = (XYZ(i+1,:) - XYZ(i,:)) / norm(...)
%       pick a trial ey, cross => ez, cross => ey, 
%       then apply twist from fem.structural_twist if you like.

exdir = ( XYZ(iElem+1,:) - XYZ(iElem,:) );
L = norm(exdir);
if L<1e-12
    R3 = eye(3); 
    return
end
exdir = exdir(:)/L;

% pick a trial "eyTrial" that isn't parallel to exdir
eyTrial = [0;1;0];
if abs(dot(exdir, eyTrial))>0.99
    eyTrial = [0;0;1];
end
ez = cross(exdir, eyTrial);
ez = ez / (norm(ez)+1e-16);
ey = cross(ez, exdir);

R3_noTwist = [exdir, ey, ez];  % 3x3

% If structural twist is stored in fem.structural_twist
% e.g. we define twist angle about x-axis:
twVal = 0;
if isfield(fem,'structural_twist')
   if length(fem.structural_twist)>= iElem
       twVal = fem.structural_twist(iElem);
   end
end
if abs(twVal) < 1e-12
    R3 = R3_noTwist;
else
    % rotate about local x by twVal
    cT = cos(twVal); sT= sin(twVal);
    Rtw = [1  0   0
           0  cT -sT
           0  sT  cT];
    R3 = R3_noTwist * Rtw;
end

% optionally check orthonormal
if norm(R3.'*R3 - eye(3),'fro')>1e-6
   warning('SegmentRotation not orthonormal?');
end

end



function c = cross(a,b)
% Overload if needed. Otherwise built-in cross is fine in MATLAB
c = [ a(2)*b(3) - a(3)*b(2);
      a(3)*b(1) - a(1)*b(3);
      a(1)*b(2) - a(2)*b(1) ];
end


function [psi1_local, psi1_global] = computeMomentumShapes(fem, keepDofs, Mred, phi1_global)
% COMPUTEMOMENTUMSHAPES   psi1 = M * phi1 in local coords.
% We'll do  node-by-node or sub-mid approach, similar to buildNonlinearGammas.

 nm = size(phi1_global,2);
 nd = size(phi1_global,1);
 psi1_local = zeros(nd,nm);

 ne = fem.num_elem;
 for e = 1:ne
    conn = fem.connectivities(e,:);
    subList= { [conn(1),conn(3)], [conn(3),conn(2)] };
    xyz  = fem.coordinates(conn,:);
    for subId=1:2
       nodes12 = subList{subId};
       dofA = (nodes12(1)-1)*6+(1:6);
       dofB = (nodes12(2)-1)*6+(1:6);
       subDofs = union(dofA,dofB);
       localInd= [];
       for dd=subDofs
         x = find(keepDofs==dd);
         localInd= [localInd, x];
       end
       if length(localInd)<12, continue; end
       leftI= localInd(1:6); rightI= localInd(7:12);

       R3= buildElementRotation_sub(fem,e,subId);
       R6= blkdiag(R3,R3);

       phi1L_glob = phi1_global(leftI,:);
       phi1R_glob = phi1_global(rightI,:);
       % local:
       phi1L_loc  = R6 * phi1L_glob;
       phi1R_loc  = R6 * phi1R_glob;

       % multiply by M sub-block => e.g. Mred(leftI,leftI)* phi1L
       M_L = Mred(leftI,  leftI);
       M_R = Mred(rightI, rightI);

       for j=1:nm
         p1L = M_L*phi1L_glob(:,j); 
         p1L_loc = R6*p1L;   % local
         % store
         psi1_local(leftI, j) = p1L_loc;
         % similarly for p1R
         p1R = M_R*phi1R_glob(:,j);
         p1R_loc= R6*p1R; 
         psi1_local(rightI,j)= p1R_loc;

         psi1_global(leftI, j) = p1L;
         psi1_global(rightI, j) = p1R;


       end
    end
 end
end



function plotNodes3D_V2(fem)
    figure; hold on; axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Node Layout');

    coords = fem.coordinates; % typically [num_node x 3] or [3 x num_node]
    if size(coords,2)==3 && size(coords,1)==fem.num_node
        % okay [node_i, (x,y,z)]
    elseif size(coords,1)==3 && size(coords,2)==fem.num_node
        coords = coords';
    end

    plot3(coords(:,1), coords(:,2), coords(:,3), 'ro','MarkerFaceColor','r');
    for i=1:fem.num_node
        text(coords(i,1), coords(i,2), coords(i,3), sprintf('Node %d', i));
    end
    plotElements3D(fem);
     hold off;
end

function plotElements3D(fem)
    % figure; hold on; axis equal;
    coords = fem.coordinates;
    if size(coords,2)==3, coords=coords'; end

    for e=1:fem.num_elem
        nodesE = fem.connectivities(e,:);  % e.g. [n1,n2,n3]
        xyz = coords(:,nodesE);
        % For 3-node beams, you might just plot lines from node1->node2->node3:
        plot3(xyz(1,[1,3,2]), xyz(2,[1,3,2]), xyz(3,[1,3,2]), '-b','LineWidth',1.5);
    end
end


function plotNodes3D(fem)
% plot the node coordinates with boundary colors
figure('Name','Pazy Nodes 3D','Color','white'); hold on; axis equal
XYZ = fem.coordinates;
bcs = fem.boundary_conditions;
for n=1:fem.num_node
   x= XYZ(n,1); y=XYZ(n,2); z=XYZ(n,3);
   switch bcs(n)
     case 1, clr='r';  % clamp
     case -1,clr='g';  % tip
     otherwise, clr='b';
   end
   plot3(x,y,z,'o','MarkerFaceColor',clr,'MarkerSize',5,'MarkerEdgeColor','k');
   text(x,y,z, sprintf('%d',n),'Color','k','FontSize',9);

end
% also draw lines for each element
ne = fem.num_elem;
for e=1:ne
   c= fem.connectivities(e,:);
   cXYZ= XYZ(c,:);
   dx = cXYZ(:,1);
   dy = cXYZ(:,2);
   dz = cXYZ(:,3);
   plot3(cXYZ([1,3, 2],1), cXYZ([1,3, 2],2), cXYZ([1,3, 2],3), '-k','LineWidth',1.2);
end
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Nodes, boundary conditions, and connectivity');
view(3);
grid on;
hold off;
end

function activeDofs = defineActiveDofs_pazyWing(fem, dofBC)
%DEFINEACTIVEDOFS_PAZYWING  Identify the master (active) DOFs for condensation
%
% Inputs:
%   numNodes    = total number of FE nodes in the Pazy wing mesh
%   dofsPerNode = number of DOFs per node (6 if 3D beam, 3 if 2D)
%   clampedNode = ID of the clamped (root) node (often node 0)
%   lumpNodes   = array of node IDs that have lumped masses or rods
%                 e.g. [5, 9, 12] if your model lumps are placed there
%
% Output:
%   activeDofs  = column vector listing the global DOFs that remain
%                 after condensation. (All others are 'omitted' or slaved.)
    % disp(fem)
    activeDofs = [];
    % disp(fem)
    dofPerNode = 6;

    for n = 1 : fem.num_node
        % If node is clamped, skip it entirely
        if ismember(n, dofBC)
            continue
        end

        % If node is among the lumps, keep all its DOFs
        if ismember(n, fem.lumped_mass_nodes)
            baseIndex = double((n-1)*dofPerNode); % (1-based "n", so subtract 1)
            theseDofs = baseIndex + (1 : dofPerNode);
            activeDofs = [activeDofs, theseDofs]; %#ok<AGROW>
        else
            % % Possibly keep only translations or skip entirely,
            % % % e.g. keep only translation DOFs for non-lump nodes:
            %   baseIndex = double((n-1)*dofPerNode);
            %   transDofs = baseIndex + [1 2 3]; % e.g. x,y,z only
            %   activeDofs = [activeDofs, transDofs];
            % 
            % Or skip them altogether:
            % do nothing
        end
    end

    % remove duplicates & sort:
    activeDofs = unique(activeDofs, 'sorted');

end


function [Xfull_time, orientation_time] = recoverFull3Dsolution( ...
        fem, T_oa, phi0, q0Hist, dt, nodeCenters, warpingMatrices )
%
%  This aims to do eq. (2.134) from Esteban or eq. (5.34) from Wang to 
%  reconstruct the omitted dofs or cross-section warp from the 
%  known skeleton solution q0(t).
%
%  Inputs:
%    fem          => your geometry (num_node, connectivities, etc.)
%    T_oa         => the [nOmitted x nActive] transformation or 
%                    some user approach for splitted DOFs
%    phi0         => the [nActive x nModes], the LNMs at active dofs
%    q0Hist       => [nModes x nTime], the time history of modal amplitude 
%                    for "displacement/rotation" or "q0(t)" in beam sense
%    dt           => time step
%    nodeCenters  => the p^i vectors or references for cross-sections
%    warpingMatrices => T_w^{ik}, etc. 
%
%  Output:
%    Xfull_time      => array storing the recovered 3D coordinates of 
%                       every node (both active and omitted) at each time step
%    orientation_time=> if you want local cross-section orientation, we store
%                       that, or you can do another approach
%

nTime = size(q0Hist,2);
nNodes = fem.num_node;
Xfull_time = cell(nTime,1);  % or use a 3D array [nNode x 3 x nTime]

% The main principle:
%   1) we get the “reference node position” r_a^k(t) from the beam solution 
%      i.e. from your NMROM integration, you got [q1(t), q2(t)] 
%      or something that also yields the master node positions & rotations 
%   2) For each omitted node i that belongs to cross-section k, do eq. (2.134)
%   3) Possibly do warping

for it = 1:nTime
    % get the current amplitude q0(t):
    q0_t = q0Hist(:,it);

    % get the master node positions from the “beam eqn.”:
    % you presumably store them in e.g. “ra(t), Rab(t)” 
    % for each master node. That depends on your structural solver.
    [rMaster, RMaster] = getSkeletonSolutionAtTime(it, fem, phi0, q0_t);
    % rMaster => [nMasterNodes x 3], RMaster => [3x3 x nMasterNodes]

    coords3d = zeros(nNodes, 3);
    for k = 1:nNodes
        % if k is a “master” node => coords3d(k,:) = rMaster(k,:)
        % if k is omitted => 
        if isMasterNode(k, fem) 
            coords3d(k,:) = rMaster(k,:); 
        else
            % find the “owner” or the cross-section center node that k is 
            % attached to => we call it node 'kC'. 
            kC = findCenterForOmittedNode(k, fem);
            % see eq. (2.134):
            p_i   = nodeCenters(k,:);  % e.g. the vector from node kC 
            Twik  = warpingMatrices{k}; 
            % T_w^{ik} is [3xN], so that T_w^{ik} * phi0^j => ...
            % sum over j => we do sum_j q0_t(j)* [T_w^{ik} * phi0(:,j)]
            warpVec = zeros(3,1);
            for j=1:length(q0_t)
                warpVec = warpVec + q0_t(j) * ( Twik * phi0(:,j) );
            end
            coords3d(k,:) = rMaster(kC,:) + ...
                ( RMaster(:,:,kC)* ( p_i(:) + warpVec ) ).';
        end
    end

    Xfull_time{it} = coords3d;
end
orientation_time = []; % or store the cross-section orientation if you like

end

function isM = isMasterNode(k, fem)
% e.g. if your keepDofs cover node k, it is a "master"? 
% or if boundary_conditions(k) is not omitted, etc. 
isM = true;  % or do a real logic
end

function debugPlotIntrinsicModes3D( Xfull_time, iTime, fem )
% Xfull_time => cell of nTime each with [nNode x 3], iTime => pick snapshot
coords = Xfull_time{iTime};
figure; hold on; axis equal
title(sprintf('Recovered 3D shape tstep=%d', iTime))
for e = 1:fem.num_elem
   c = fem.connectivities(e,:);
   plot3( coords([1,3, 2],1), coords([1,3, 2],2), coords([1,3, 2],3), '-o');
end
grid on;
end

function debugPlotMode3D(fem, phi0, modeNum, scale)
% DEBUGPLOTMODE3D  simple plot of the structure with mode-shape
%   phi0 => [ndofs x nmodes], modeNum => which mode to plot, scale => amplitude
figure('Name','ModePlot','Color','w');
hold on; axis equal; grid on;
coords = fem.coordinates;  % [nNodes x 3]
nNodes = fem.num_node;
dofPerNode=6;
phiAll = phi0;
thisMode= phiAll(:, modeNum);  % a 6*nNodes vector
% disp(fem.connectivities)

% Original
for e=1:fem.num_elem
    conn = fem.connectivities(e,:);
    xyz0 = coords(conn, :);  % e.g. Nx3
    plot3(xyz0([1,3, 2],1), xyz0([1,3, 2],2), xyz0([1,3, 2],3), '-k','LineWidth',1.2);
end

% Deformed
Deformed= zeros(nNodes,3);
for n=1:nNodes
    dofRange = double((n-1)*dofPerNode) + (1:3);
    disp3   = thisMode(dofRange)*scale;
    Deformed(n,:) = coords(n,:) + disp3.';
end
for e=1:fem.num_elem
    conn = fem.connectivities(e,:);
    xyzD = Deformed(conn,:);
    dx = xyzD(:,1);
    dy = xyzD(:,2);
    dz = xyzD(:,3);
    plot3(xyzD([1,3, 2],1), xyzD([1,3, 2],2), xyzD([1,3, 2],3), '-ro','LineWidth',1.5);
end
%  figure; hold on; axis equal;
%  coords = fem.coordinates;  % [n x 3]
%  nNodes = fem.num_node;
%  dofPerNode=6;
% 
%  thisMode= phi0(:, modeNum);
% 
%  for e=1:fem.num_elem
%     conn = fem.connectivities(e,:);  % e.g. [n1,n2,n3]
%     for c=1:(length(conn)-1)
%         nA = conn(c);
%         nB = conn(c+1);
%         baseA = coords(nA,:);
%         baseB = coords(nB,:);
%         dispA = phi0( (nA-1)*6+(1:3), modeNum )'*scale;
%         dispB = phi0( (nB-1)*6+(1:3), modeNum )'*scale;
%         plot3([baseA(1)+dispA(1), baseB(1)+dispB(1)], ...
%               [baseA(2)+dispA(2), baseB(2)+dispB(2)], ...
%               [baseA(3)+dispA(3), baseB(3)+dispB(3)], 'r-o','LineWidth',1.5);
%     end
% end
 view(3)

 title(sprintf('Mode #%d shape with scale=%.2f', modeNum,scale));
end

function debugPlotPhi2(phi2_loc, chain, modeNum, localOrGlobal, plotMoments)
% DEBUGPLOTPHI2  Debugging plot of the piecewise-constant force distribution
%  on a single chain. 
%
%  Inputs:
%    phi2_loc  => [6*nSub x nModes], piecewise constant force modes 
%                 in either local or global coords 
%    chain     => array of node IDs, e.g. [1,2,3,...] for plotting "span"
%    modeNum   => which mode index to plot
%    localOrGlobal => 'local' or 'global' (string) for plot labeling 
%    plotMoments => boolean or string indicating whether to plot M_x,y,z too
%
%  Plots F_x, F_y, F_z, (and optionally M_x, M_y, M_z) vs sub-element index.

if nargin<5
    plotMoments = true;
end
numModes  = size(phi2_loc,2);
if modeNum<1 || modeNum>numModes
    error('Invalid modeNum: must be between 1..%d', numModes);
end

nSub = length(chain)-1;   % number of sub-elements 
if 6*nSub ~= size(phi2_loc,1)
    error('Mismatch: chain of length %d => %d sub-e, but phi2_loc has %d rows.', ...
          length(chain), nSub, size(phi2_loc,1));
end

% if chain(1)>chain(end)
%     [B, Idx] = sort(chain, 'descend')
% end

% Extract the sub-block of phi2 for this mode
%   sub i => row block => (i-1)*6 + (1:6)
Fxyz = zeros(nSub,3);   % F_x, F_y, F_z for each sub
Mxyz = zeros(nSub,3);   % M_x, M_y, M_z for each sub
for i = 1:nSub
   rowIdx = (i-1)*6 + (1:6);
   phi2_i = phi2_loc(rowIdx, modeNum);  % [6x1]
   Fxyz(i,:) = phi2_i(1:3).';
   Mxyz(i,:) = phi2_i(4:6).';
end

% x-axis: sub-element index i or maybe midpoint's coordinate 
% For a simpler debug, just plot i=1..nSub
iVec = 1:nSub;   % or store the actual s_i midpoints if you want

figure('Name','DebugPlotPhi2','Color','w'); 
subplot(2,1,1);
plot(iVec, Fxyz(:,1),'r-o','LineWidth',1.2,'DisplayName','Fx');
hold on;
plot(iVec, Fxyz(:,2),'g-o','LineWidth',1.2,'DisplayName','Fy');
plot(iVec, Fxyz(:,3),'b-o','LineWidth',1.2,'DisplayName','Fz');
xlabel('sub-element index i');
ylabel('Force amplitude');
title(sprintf('\\phi_2 Force Components (Mode %d, %s coords)', modeNum, localOrGlobal));
legend('Location','best'); 
grid on;

if plotMoments
   subplot(2,1,2);
   plot(iVec, Mxyz(:,1),'r-s','LineWidth',1.2,'DisplayName','Mx');
   hold on;
   plot(iVec, Mxyz(:,2),'g-s','LineWidth',1.2,'DisplayName','My');
   plot(iVec, Mxyz(:,3),'b-s','LineWidth',1.2,'DisplayName','Mz');
   xlabel('sub-element index i');
   ylabel('Moment amplitude');
   title(sprintf('\\phi_2 Moment Components (Mode %d, %s coords)', modeNum, localOrGlobal));
   legend('Location','best');
   grid on;
end

end

function debugPlotPhi2Full(fem, phi2_loc, chain, modeNum, localOrGlobal, plotMoments, modeVal)
% DEBUGPLOTPHI2  Debugging plot of the piecewise-constant force distribution
%  on a single chain. 
%
%  Inputs:
%    phi2_loc  => [6*nSub x nModes], piecewise constant force modes 
%                 in either local or global coords 
%    chain     => array of node IDs, e.g. [1,2,3,...] for plotting "span"
%    modeNum   => which mode index to plot
%    localOrGlobal => 'local' or 'global' (string) for plot labeling 
%    plotMoments => boolean or string indicating whether to plot M_x,y,z too
%
%  Plots F_x, F_y, F_z, (and optionally M_x, M_y, M_z) vs sub-element index.

if nargin<5
    plotMoments = true;
end
numModes  = size(phi2_loc,2);
if modeNum<1 || modeNum>numModes
    error('Invalid modeNum: must be between 1..%d', numModes);
end
XYZ = fem.coordinates;
orig = XYZ(1,:);
del_idx = [];
for i = 2:fem.num_node

    if XYZ(i,2)<XYZ(i-1,2)
        del_idx(i) = i-1;
        continue;

    end 
    rmid(i-1,:) = .5*(XYZ(i,:)+XYZ(i-1,:));
    % also need last midpt
    if i == fem.num_node
        rmid(i,:) = .5*(XYZ(i,:)+orig);
    end
end
del = find(del_idx);
this =del_idx(del); 
rmid(this,:) = [];


[rMidSort, Idx] = sort(rmid(:,2));
[rMidSort2, Idx2] = sort(rMidSort(9:16), 'descend');

% if modeVal == 2
%     Idx = [Idx(1:8); Idx2];
% end
nSub = length(chain)-2;   % number of sub-elements 
if 6*nSub ~= size(phi2_loc,1)
    error('Mismatch: chain of length %d => %d sub-e, but phi2_loc has %d rows.', ...
          length(chain), nSub, size(phi2_loc,1));
end



% Extract the sub-block of phi2 for this mode
%   sub i => row block => (i-1)*6 + (1:6)
Fxyz = zeros(nSub,3);   % F_x, F_y, F_z for each sub
Mxyz = zeros(nSub,3);   % M_x, M_y, M_z for each sub
% for i = 1:nSub
%    rowIdx = (i-1)*6 + (1:6);
%    phi2_i = phi2_loc(rowIdx, modeNum);  % [6x1]
%    Fxyz(i,:) = phi2_i(1:3).';
%    Mxyz(i,:) = phi2_i(4:6).';
% end

for i = 1:nSub
   nIdx = Idx(i);
   rowIdx = (nIdx-1)*6 + (1:6);
   phi2_i = phi2_loc(rowIdx, modeNum);  % [6x1]
   Fxyz(i,:) = phi2_i(1:3).';
   Mxyz(i,:) = phi2_i(4:6).';
end

% x-axis: sub-element index i or maybe midpoint's coordinate 
% For a simpler debug, just plot i=1..nSub
iVec = 1:nSub;   % or store the actual s_i midpoints if you want
if modeVal == 2
    figure('Name','Phi2 Mode Shapes','Color','w'); 
    subplot(2,1,1);
    plot(rMidSort, Fxyz(:,1),'r-o','LineWidth',1.2,'DisplayName','Fx');
    hold on;
    plot(rMidSort, Fxyz(:,2),'g-o','LineWidth',1.2,'DisplayName','Fy');
    plot(rMidSort, Fxyz(:,3),'b-o','LineWidth',1.2,'DisplayName','Fz');
    xlabel('sub-element index i');
    ylabel('Force amplitude');
    title(sprintf('\\phi_2 Full Winged Force Components (Mode %d, %s coords)', modeNum, localOrGlobal));
    legend('Location','best'); 
    grid on;

    if plotMoments
       subplot(2,1,2);
       plot(rMidSort, Mxyz(:,1),'r-s','LineWidth',1.2,'DisplayName','Mx');
       hold on;
       plot(rMidSort, Mxyz(:,2),'g-s','LineWidth',1.2,'DisplayName','My');
       plot(rMidSort, Mxyz(:,3),'b-s','LineWidth',1.2,'DisplayName','Mz');
       xlabel('sub-element index i');
       ylabel('Moment amplitude');
       title(sprintf('\\phi_2 Moment Components (Mode %d, %s coords)', modeNum, localOrGlobal));
       legend('Location','best');
       grid on;
    end

elseif modeVal == 1
    figure('Name','Phi1 Mode Shapes','Color','w'); 
    subplot(2,1,1);
    plot(rMidSort, Fxyz(:,1),'r-o','LineWidth',1.2,'DisplayName','Fx');
    hold on;
    plot(rMidSort, Fxyz(:,2),'g-o','LineWidth',1.2,'DisplayName','Fy');
    plot(rMidSort, Fxyz(:,3),'b-o','LineWidth',1.2,'DisplayName','Fz');
    xlabel('sub-element index i');
    ylabel('Velocity amplitude');
    title(sprintf('\\phi_1 Full Winged Velocity Components (Mode %d, %s coords)', modeNum, localOrGlobal));
    legend('Location','best'); 
    grid on;

    if plotMoments
       subplot(2,1,2);
       plot(rMidSort, Mxyz(:,1),'r-s','LineWidth',1.2,'DisplayName','Mx');
       hold on;
       plot(rMidSort, Mxyz(:,2),'g-s','LineWidth',1.2,'DisplayName','My');
       plot(rMidSort, Mxyz(:,3),'b-s','LineWidth',1.2,'DisplayName','Mz');
       xlabel('sub-element index i');
       ylabel('Angular Vel amplitude');
       title(sprintf('\\phi_1 Ang Val Components (Mode %d, %s coords)', modeNum, localOrGlobal));
       legend('Location','best');
       grid on;
    end
else
    error('Incorrect ModeVal selected. Only for Phi1 and Phi2 validation')
end
end

function debugPlotPhi1Cubic(fem, chain, phi1Fun, modeNumber, modeVal, plotMoments)
% debugPlotPhi1Cubic
%  Evaluate the piecewise polynomial for the specified modeNumber
%  at many sample points, and then plot something, e.g. the "Fy" or
%  the "magnitude" vs. the y-coordinate in the real 3D geometry.

nSeg = length(chain)-1;
XYZ = fem.coordinates;
nSamplesPerSeg= 20;

% We'll store e.g. "arc array" or the actual Y positions
allY  = [];
allVal = [];
allFx= []; allFz = [];
allFy = [];
allMx =[]; allMy =[]; allMz=[];
arcSoFar= 0;

% Access the correct mode's piecewise struct
myFun = phi1Fun{modeNumber};
if modeVal == 1

    for iSeg = 1:nSeg
        s0 = myFun.seg(iSeg).xrange(1);
        s1 = myFun.seg(iSeg).xrange(2);
        ds = s1 - s0;

        nA = chain(iSeg);   nB = chain(iSeg+1);
        Acoord = XYZ(nA,:); Bcoord = XYZ(nB,:);
        segVec = (Bcoord - Acoord);   % direction in global coords
        SegIdx = double((iSeg-1)*6);
        SegIdx = SegIdx + (1:6);
        for iS=0:nSamplesPerSeg
           xi = iS / nSamplesPerSeg;   % 0..1
           sLoc= ds* xi;              % local coordinate in [0..ds]
           % sLocD = sLoc - Acoord(2);
           phiVal = evalCubic6(myFun.seg(iSeg).coefs, sLoc);

           % The 2D or 3D location in space => param = Acoord + segVec*(xi).
           thisPoint = Acoord + segVec*xi;   % e.g. linear
           yVal= thisPoint(2);  % pick the Y coord

           Fy  = phiVal(2);     % for example, the second component => velocity in y
           allVal(end+1,:) = thisPoint;

           phiAll(SegIdx,iS+1) = phiVal;

           allY(end+1)= yVal;
           allFy(end+1)= Fy;
           allFx(end+1) =  phiVal(1);
           allFz(end+1) = phiVal(3);
           allMx(end+1) =  phiVal(4);
           allMy(end+1) = phiVal(5);
           allMz(end+1) =  phiVal(6);
        end
    end



    % figure('Name','Phi1 Mode Shapes','Color','w'); 
    subplot(2,1,1);
    plot(allY, allFx,'r-o','LineWidth',1.2,'DisplayName','Fx');
    hold on;
    plot(allY, allFy,'g-o','LineWidth',1.2,'DisplayName','Fy');
    plot(allY, allFz,'b-o','LineWidth',1.2,'DisplayName','Fz');
    xlabel('sub-element index i');
    ylabel('Velocity amplitude');
    title(sprintf('\\phi_1 Full Winged Velocity Components (Mode %d)', modeNumber));
    legend('Location','best'); 
    grid on;

    if plotMoments
       subplot(2,1,2);
       plot(allY, allMx,'r-s','LineWidth',1.2,'DisplayName','Mx');
       hold on;
       plot(allY, allMy,'g-s','LineWidth',1.2,'DisplayName','My');
       plot(allY, allMz,'b-s','LineWidth',1.2,'DisplayName','Mz');
       xlabel('sub-element index i');
       ylabel('Angular Vel amplitude');
       title(sprintf('\\phi_1 Ang Val Components (Mode %d)', modeNumber));
       legend('Location','best');
       grid on;
    end

elseif modeVal == 2
    % for iSeg = 1:nSeg
    %     s0 = myFun.seg(iSeg).xrange(1);
    %     s1 = myFun.seg(iSeg).xrange(2);
    %     ds = s1 - s0;
    % 
    %     nA = chain(iSeg);   nB = chain(iSeg+1);
    %     Acoord = XYZ(nA,:); Bcoord = XYZ(nB,:);
    %     segVec = (Bcoord - Acoord);   % direction in global coords
    %     SegIdx = double((iSeg-1)*6);
    %     SegIdx = SegIdx + (1:6);
    %     for iS=0:nSamplesPerSeg
    %        xi = iS / nSamplesPerSeg;   % 0..1
    %        sLoc= ds* xi;              % local coordinate in [0..ds]
    %        phiVal = myFun.seg(iSeg).val;
    % 
    %        % The 2D or 3D location in space => param = Acoord + segVec*(xi).
    %        thisPoint = Acoord + segVec*xi;   % e.g. linear
    %        yVal= thisPoint(2);  % pick the Y coord
    % 
    %        Fy  = phiVal(2);     % for example, the second component => velocity in y
    %        allVal(end+1,:) = thisPoint;
    % 
    %        phiAll(SegIdx,iS+1) = phiVal;
    % 
    %        allY(end+1)= yVal;
    %        allFy(end+1)= Fy;
    %        allFx(end+1) =  phiVal(1);
    %        allFz(end+1) = phiVal(3);
    %        allMx(end+1) =  phiVal(4);
    %        allMy(end+1) = phiVal(5);
    %        allMz(end+1) =  phiVal(6);
    %     end
    % end
     for iSeg = 1:nSeg
        s0 = myFun.seg(iSeg).xrange(1);
        s1 = myFun.seg(iSeg).xrange(2);
        ds = s1 - s0;

        nA = chain(iSeg);   nB = chain(iSeg+1);
        Acoord = XYZ(nA,:); Bcoord = XYZ(nB,:);
        segVec = (Bcoord - Acoord);   % direction in global coords
        SegIdx = double((iSeg-1)*6);
        SegIdx = SegIdx + (1:6);
        for iS=0:nSamplesPerSeg
           xi = iS / nSamplesPerSeg;   % 0..1
           sLoc= ds* xi;              % local coordinate in [0..ds]
           phiVal = evalCubic6(myFun.seg(iSeg).coefs, sLoc);

           % The 2D or 3D location in space => param = Acoord + segVec*(xi).
           thisPoint = Acoord + segVec*xi;   % e.g. linear
           yVal= thisPoint(2);  % pick the Y coord

           Fy  = phiVal(2);     % for example, the second component => velocity in y
           allVal(end+1,:) = thisPoint;

           phiAll(SegIdx,iS+1) = phiVal;

           allY(end+1)= yVal;
           allFy(end+1)= Fy;
           allFx(end+1) =  phiVal(1);
           allFz(end+1) = phiVal(3);
           allMx(end+1) =  phiVal(4);
           allMy(end+1) = phiVal(5);
           allMz(end+1) =  phiVal(6);
        end
    end
    % figure('Name','Phi2 Mode Shapes','Color','w'); 
    subplot(2,1,1);
    plot(allY, allFx,'r-o','LineWidth',1.2,'DisplayName','Fx');
    hold on;
    plot(allY, allFy,'g-o','LineWidth',1.2,'DisplayName','Fy');
    plot(allY, allFz,'b-o','LineWidth',1.2,'DisplayName','Fz');
    xlabel('sub-element index i');
    ylabel('Force amplitude');
    title(sprintf('\\phi_2 Full Winged Force Components (Mode %d)', modeNumber));
    legend('Location','best'); 
    grid on;

    if plotMoments
       subplot(2,1,2);
       plot(allY, allMx,'r-s','LineWidth',1.2,'DisplayName','Mx');
       hold on;
       plot(allY, allMy,'g-s','LineWidth',1.2,'DisplayName','My');
       plot(allY, allMz,'b-s','LineWidth',1.2,'DisplayName','Mz');
       xlabel('sub-element index i');
       ylabel('Moment amplitude');
       title(sprintf('\\phi_2 Moment Components (Mode %d)', modeNumber));
       legend('Location','best');
       grid on;
    end

else
    error('Incorrect ModeVal selected. Only for Phi1 and Phi2 validation')
end

end

function debugPlotPhi0Cubic(fem, chain, phi1Fun, modeNumber)
% debugPlotPhi1Cubic
%  Evaluate the piecewise polynomial for the specified modeNumber
%  at many sample points, and then plot something, e.g. the "Fy" or
%  the "magnitude" vs. the y-coordinate in the real 3D geometry.

nSeg = length(chain)-1;
XYZ = fem.coordinates;
nSamplesPerSeg= 20;

% We'll store e.g. "arc array" or the actual Y positions
allY  = [];
allVal = [];
allFy = [];
arcSoFar= 0;

% Access the correct mode's piecewise struct
myFun = phi1Fun{modeNumber};

for iSeg = 1:nSeg
    s0 = myFun.seg(iSeg).xrange(1);
    s1 = myFun.seg(iSeg).xrange(2);
    ds = s1 - s0;

    nA = chain(iSeg);   nB = chain(iSeg+1);
    Acoord = XYZ(nA,:); Bcoord = XYZ(nB,:);
    segVec = (Bcoord - Acoord);   % direction in global coords
    SegIdx = double((iSeg-1)*6);
    SegIdx = SegIdx + (1:6);
    for iS=0:nSamplesPerSeg
       xi = iS / nSamplesPerSeg;   % 0..1
       sLoc= ds* xi;              % local coordinate in [0..ds]
       phiVal = evalCubic6(myFun.seg(iSeg).coefs, sLoc);

       % The 2D or 3D location in space => param = Acoord + segVec*(xi).
       thisPoint = Acoord + segVec*xi;   % e.g. linear
       yVal= thisPoint(2);  % pick the Y coord
       Fy  = phiVal(2);     % for example, the second component => velocity in y
       allVal(end+1,:) = thisPoint;

       phiAll(SegIdx,iS+1) = phiVal;

       allY(end+1)= yVal;
       allFy(end+1)= Fy;
    end
end

% figure('Name','Phi1Cubic debug'); 
% subplot(2,1,1);
% plot(allY, allFy,'-o','LineWidth',1.3);
% xlabel('Global Y coordinate (m)');
% ylabel('\phi_1^{(2)}(y)');  % the 2nd dof or your choice
% title(sprintf('Mode %d: dof #2 (Fy) vs. Y', modeNumber));
% grid on;

figure('Name','ModePlot','Color','w');
hold on; axis equal; grid on;
coords = fem.coordinates;  % [nNodes x 3]
nNodes = fem.num_node;
dofPerNode=6;
% phiAll = phi0;
% thisMode= phiAll(:, modeNum);  % a 6*nNodes vector
% disp(fem.connectivities)
scale =1;
% Original
for e=1:fem.num_elem
    conn = fem.connectivities(e,:);
    xyz0 = coords(conn, :);  % e.g. Nx3
    plot3(xyz0([1,3, 2],1), xyz0([1,3, 2],2), xyz0([1,3, 2],3), '-k','LineWidth',1.2);
end

% Deformed
Deformed= zeros(nNodes,3);
for n=1:nNodes
    dofRange = double((n-1)*dofPerNode) + (1:3);
    disp3   = phiAll(dofRange)*scale;
    Deformed(n,:) = coords(n,:) + disp3.';
end
for e=1:fem.num_elem
    conn = fem.connectivities(e,:);
    xyzD = Deformed(conn,:);
    dx = xyzD(:,1);
    dy = xyzD(:,2);
    dz = xyzD(:,3);
    plot3(xyzD([1,3, 2],1), xyzD([1,3, 2],2), xyzD([1,3, 2],3), '-ro','LineWidth',1.5);
end
% subplot(2,1,2); 
% figure('Name','Phi1Cubic debug'); 
% plot3(allVal, allFy,'-o','LineWidth',1.3);
% xlabel('Global Y coordinate (m)');
% ylabel('\phi_1^{(2)}(y)');  % the 2nd dof or your choice
% title(sprintf('Mode %d: dof #2 (Fy) vs. Y', modeNumber));

end

%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION #1: Build piecewise cubic polynomials with 4 BC for dofs 2->5, 3->6
function seg = buildPiecewiseLin(nodeVals_6_n, elemLen,chain, XYZ)
% nodeVals_6_n => [6 x nChain], each col => 6 dof values at node i
% We interpret the "bending" dof=2,3 with slope dof=5,6 => a 4-BC cubic:
%   y(0)= yA(2), y'(0)= yA(5), y(ds)= yB(2), y'(ds)= yB(5)  (like eq. 2.112).
% For dof=1,4 => linear.  For dof=5,6 => the derivative polynomial from dof=2,3.
%
% We produce a final .seg(i).coefs => [6 x 4]. The row d in [1..6] is the
% polynomial for that dof in the sub‐interval. 
%
nChain = size(nodeVals_6_n,2);
segArr= struct('xrange',[],'coefs',[]);
s0=XYZ(chain(1),2);
for i=1:(nChain-1)
   nID = chain(i);
   ds= elemLen(i);
   s1 =s0+ds;
   dofA = nodeVals_6_n(:, i);
   dofB = nodeVals_6_n(:, i+1);
   % R_id = double(i-1)*3+(1:3);
   % Rglobal0_s = Rglobal0(R_id,:);
   % First Lets do lin interp eq 2.112
   % phi0_Fun_s = @(s) Rglobal0_s*(dofA*((s1-s)/ds)+dofB*((s-s0)/ds));
   % Need an Eval funct --> storing only dofs
   phi0_Fun_coeff = [dofA  dofB];

   % coefs_6_4= zeros(6,4);  % each row => a0..a3
   % % DOF=1 => linear
   % [a0,a1] = linearBC(dofA(1), dofB(1), ds);
   % coefs_6_4(1,:) = [a0, a1, 0, 0];
   % % DOF=4 => linear
   % [b0,b1] = linearBC(dofA(4), dofB(4), ds);
   % coefs_6_4(4,:) = [b0, b1, 0, 0];
   % 
   % % DOF=2 => 4 BC with dof=5
   % %   y(0)= dofA(2),  y'(0)= dofA(5),
   % %   y(ds)= dofB(2), y'(ds)= dofB(5).
   % [a0,a1,a2,a3] = cubic4bc( dofA(2), dofA(5), dofB(2), dofB(5), ds );
   % coefs_6_4(2,:) = [a0,a1,a2,a3];
   % % slope dof=5 => derivative => d/ds of that polynomial => row(5)
   % % derivative => p'(s)= (a1) + 2(a2)*s + 3(a3)*s^2
   % % => but we want a polynomial of form c0 + c1*s + c2*s^2 + c3*s^3
   % % => c0= a1, c1= 2*a2, c2= 3*a3, c3= 0
   % coefs_6_4(5,1)= a1;     
   % coefs_6_4(5,2)= 2*a2;  
   % coefs_6_4(5,3)= 3*a3;  
   % coefs_6_4(5,4)= 0;
   % 
   % % DOF=3 => 4 BC with dof=6
   % [p0,p1,p2,p3] = cubic4bc( dofA(3), dofA(6), dofB(3), dofB(6), ds );
   % coefs_6_4(3,:) = [p0,p1,p2,p3];
   % coefs_6_4(6,1)= p1; 
   % coefs_6_4(6,2)= 2*p2; 
   % coefs_6_4(6,3)= 3*p3; 
   % coefs_6_4(6,4)= 0;

   segArr(i).xrange= [s0, s0+ds];
   segArr(i).coefs = phi0_Fun_coeff;
   s0= s0+ds;
end
seg= segArr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION #2: Solve for linear BC
function [a0,a1] = linearBC(valA, valB, ds)
a0= valA; 
a1= (valB - valA)/ds;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION #3: Solve the 4 boundary conditions for dof=2 => slope=5, etc.
%   y(0)= yA, y'(0)= yA_d, y(ds)= yB, y'(ds)= yB_d.
function [a0,a1,a2,a3] = cubic4bc(yA, dyA, yB, dyB, ds)
% We want polynomial p(s)= a0 + a1*s + a2*s^2 + a3*s^3
%   p(0)= yA => a0= yA
%   p'(0)= dyA => a1= dyA
%   p(ds)= yB => a0 + a1*ds + a2*ds^2 + a3*ds^3= yB
%   p'(ds)= dyB => a1 + 2*a2*ds + 3*a3*ds^2= dyB
% So we have:
%   a0= yA
%   a1= dyA
% Then we solve for a2,a3 from the last 2 equations:
%   yA + dyA*ds + a2*ds^2 + a3*ds^3= yB
%   dyA + 2*a2*ds + 3*a3*ds^2= dyB
% => rewrite:
%   a2*ds^2 + a3*ds^3= yB - (yA + dyA*ds)
%   2*a2*ds + 3*a3*ds^2= dyB - dyA
%
rhs1= yB - (yA + dyA*ds); 
rhs2= (dyB - dyA);
% Let us define:
% eq1 => a2*ds^2 + a3*ds^3= rhs1
% eq2 => 2*a2*ds + 3*a3*ds^2= rhs2
% We can solve quickly by standard 2x2 approach in (a2,a3).
% eq1 => a2= (rhs1 - a3*ds^3)/ ds^2
% then eq2 => 2*((rhs1 - a3*ds^3)/ds^2)*ds + 3*a3*ds^2= rhs2
% => 2*(rhs1 - a3*ds^3)/ds + 3*a3*ds^2= rhs2
% => 2*rhs1/ds - 2*a3*ds^2 + 3*a3*ds^2= rhs2
% => 2*rhs1/ds + a3*ds^2= rhs2
% => a3= (rhs2 - 2*rhs1/ds)/ ds^2
%
a3= (rhs2 - 2*rhs1/ds)/(ds^2);
a2= (rhs1 - a3*ds^3)/ (ds^2);
a0= yA;
a1= dyA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION #4: Evaluate polynomial coefs for [6 x 4] at local s
function val6 = evalCubicCoefs(coefs_6_4, sLoc)
val6= zeros(6,1);
for d=1:6
   a0= coefs_6_4(d,1);
   a1= coefs_6_4(d,2);
   a2= coefs_6_4(d,3);
   a3= coefs_6_4(d,4);
   val6(d)= a0 + a1*sLoc + a2*sLoc^2 + a3*sLoc^3;
end
end
function idxLocal = localDofIndex(keepDofs, dofSet)
% localDofIndex
%   given e.g. dofSet = [6*(k-1)+(1:6)], find them in keepDofs 
%   returns array of local indices 
idxLocal = [];
for dd = dofSet(:).'
   iPos = find( keepDofs==dd, 1, 'first' );
   if ~isempty(iPos)
       idxLocal(end+1) = iPos; %#ok<AGROW>
   end
end

end
function [Gamma1_val, Gamma2_val] = ...
      computeGammas_withPiecewise(Phi1Fun, Phi2Fun, M, K)
  nModes = length(Phi1Fun);  % same as length(Phi2Fun)
  Gamma1_val = zeros(nModes, nModes, nModes);
  Gamma2_val = zeros(nModes, nModes, nModes);

  L1fun = @(x) L1_operator(x);
  L2fun = @(x) L2_operator(x);
  C = inv(K);
  % For j=1..nModes, k=1..nModes, l=1..nModes
  for j=1:nModes
    nSeg_j = length(Phi1Fun{j}.seg);  % # of sub‐elements
    for k=1:nModes
      nSeg_k = length(Phi1Fun{k}.seg);
      for l=1:nModes
        nSeg_l = length(Phi1Fun{l}.seg);

        % We assume the sub‐element partitions match, or we unify them
        nSeg = nSeg_j;  % (but all your beams are the same chain)
        val1 = 0;  val2=0;
        for i=1:nSeg
          % get local polynomial coefs for j,k,l
          coefs_j = Phi1Fun{j}.seg(i).coefs; 
          coefs_k = Phi1Fun{k}.seg(i).coefs;
          coefs_l = Phi1Fun{l}.seg(i).coefs;
          % sub‐interval range:
          s0 = Phi1Fun{j}.seg(i).xrange(1);
          s1 = Phi1Fun{j}.seg(i).xrange(2);
          ds= (s1 - s0);

          % pick a small Gauss rule => say 2pt
          [gp, wgt] = gaussQuad_2pt();  % gp in [‐1,1], we map to [0,ds]
          for igp=1:length(gp)
             xi = gp(igp);   w= wgt(igp);
             sLoc= 0.5*(1+xi)*ds;   % map [‐1..1] => [0..ds]
             % Evaluate phi1_j,k,l at sLoc
             x1_j = evalCubicCoefs(coefs_j, sLoc);
             x1_k = evalCubicCoefs(coefs_k, sLoc);
             x1_l = evalCubicCoefs(coefs_l, sLoc);
             integrand1 = x1_j.' * L1fun(x1_k) * (M*x1_l);
             val1 = val1 + integrand1*w; 

             % Now for Gamma2 => we need phi1_j, phi2_k, phi2_l
             % i.e. x1_j, x2_k, x2_l.  So we also must evaluate phi2Fun
             % (which might be piecewise constant => easy).
             x2_k = Phi2Fun{k}.seg(i).val;  % piecewise constant
             x2_l = Phi2Fun{l}.seg(i).val;
             integrand2 = x1_j.' * L2fun(x2_k) * (C*x2_l);
             val2 = val2 + integrand2*w;
          end
          val1= val1*(ds*0.5);  % factor from [‐1..1] => length ds
          val2= val2*(ds*0.5);
        end
        Gamma1_val(j,k,l)= val1;
        Gamma2_val(j,k,l)= val2;
      end
    end
  end
end

function [xi,w] = gaussQuad_2pt()
  xi= [-1/sqrt(3); 1/sqrt(3)];
  w=  [1; 1];
end

function val6 = evalCubic6(coefs6x4, sLoc)
% EVALCUBIC6 evaluates a 6×1 cubic polynomial at coordinate sLoc
% 
% coefs6x4: a 6×4 array of the form:
%     coefs6x4(d,:) = [a0, a1, a2, a3]
% for d = 1..6 (one row per DOF).
% sLoc: the local coordinate along the sub-element 
%       (0 ≤ sLoc ≤ length_of_sub_element).
%
% val6: a 6×1 vector of evaluated values at sLoc.

    val6 = zeros(6,1);
    for d = 1:6
        % Extract the polynomial coefficients for DOF d
        a0 = coefs6x4(d,1);
        a1 = coefs6x4(d,2);
        a2 = coefs6x4(d,3);
        a3 = coefs6x4(d,4);
        % Evaluate a0 + a1*s + a2*s^2 + a3*s^3
        val6(d) = a0 + a1*sLoc + a2*sLoc^2 + a3*sLoc^3;
    end
end


function q_tot = combineSweepTwist(alpha, theta, aoa)
    % combineSweepTwist:
    %  Returns quaternion for "rotate by alpha about z, then by theta about x."
    %
    % q_sweep = [cos(alpha/2), 0, 0, sin(alpha/2)];
    % q_twist = [cos(theta/2), sin(theta/2), 0, 0];

    q_sweep = [cos(alpha/2),  0,            0,           sin(alpha/2)];
    q_twist = [cos(theta/2),  sin(theta/2), 0,           0           ];
    q_aoa = [cos(aoa/2), 0,  sin(aoa/2),           0           ];

    q_tot = quatMultiply(q_aoa.', q_sweep.', q_twist.');
end

function q_out = quatMultiply(qAOA, qA, qB)
    % Standard Hamilton product of quaternions qA, qB
    %   each in form [q0, q1, q2, q3].
    % sAOA = qAOA(1); vAOA = qAOA(2:4);
    sA = qA(1); vA = qA(2:4);
    sB = qB(1); vB = qB(2:4);

    sAB = sA*sB - dot(vA,vB);
    vAB = sA*vB + sB*vA + cross(vA,vB);
    qAB = [sAB; vAB];

    % final rotation about angle of attack (y-axis)
    % Step 2: Final product q_out = qAOA * qAB
    sAOA = qAOA(1); vAOA = qAOA(2:4);

    sOut = sAOA*sAB - dot(vAOA, vAB);
    vOut = sAOA*vAB + sAB*vAOA + cross(vAOA, vAB);
    q_out = [sOut; vOut];

    % Normalize for safety
    q_out = q_out / norm(q_out);
end

function Sigma = rom_build_sigma_kv(phi1_loc, chain, elemLen, Eall, cfg)
% phi1_loc : [ndof x Nm], in LOCAL coordinates (same frame as Eall)
% chain    : node ids of the active chain
% elemLen  : [nSub x 1] sub-element lengths (nSub = length(chain)-1)
% Eall     : stacked 6x6 E-matrices per sub-element: size (6*nSub x 6)
% cfg      : cfg.struct.damping_model, cfg.struct.zeta or cfg.struct.Cd_fun

Nm   = size(phi1_loc,2);
% nSub = length(chain)-1
nSub = length(chain)/2-1;
Sigma = zeros(Nm,Nm);

% Build a per-sub-element C_d (linearised KV) — three options:
% 1) user-supplied function handle Cd(sIdx) -> 6x6
% 2) uniform viscous ratio 'zeta' (scalar or 6x6) => Cd = 2*zeta*I6
% 3) fallback tiny damping
I6 = eye(6);
has_fun = isfield(cfg.struct,'Cd_fun') && ~isempty(cfg.struct.Cd_fun);

for e = 1:nSub
% for e = 1:nSub-1
    ds = elemLen(e);
    row = (e-1)*6 + (1:6);
    % disp(size(Eall));
    E6  = Eall(row,:);  % L1-operator matrix at this sub-element (6x6)

    if has_fun
        Cd = cfg.struct.Cd_fun(e);                   % 6x6
    elseif isfield(cfg.struct,'zeta') && ~isempty(cfg.struct.zeta)
        z  = cfg.struct.zeta;
        if isscalar(z), Cd = 2*z*I6; else, Cd = 2*z; end   % allow 6x6
    else
        Cd = 2*(5e-3)*I6;  % default 0.5% critical, uniformly
    end

    % Linearised KV D-matrices on this sub-element (primes≈0 inside e)
    D0 = -E6*Cd*E6.';            % 6x6
    D1 =  E6*Cd - Cd*E6.';       % 6x6 (skew)
    D2 =  Cd;                    % 6x6 (SPD)

    % Grab local dofs for the end nodes of this sub-element
    nA = chain(e);   nB = chain(e+1);
    idxA = (nA-1)*6 + (1:6);
    idxB = (nB-1)*6 + (1:6);

    % Map to kept dofs once outside if you prefer; here we assume phi1_loc already condensed.
    % mid-point value and slope for each mode
    PhiA = phi1_loc(idxA,:);     % 6 x Nm
    PhiB = phi1_loc(idxB,:);     % 6 x Nm
    PhiM = 0.5*(PhiA + PhiB);    % mid value
    G    = (PhiB - PhiA)/ds;     % constant slope over e

    % Accumulate Σ using the first-derivative form (integration by parts)
    % Σ += ds * [ PhiM' D0 PhiM + PhiM' D1 G - G' D2 G ]
    Sigma = Sigma ...
          + ds * ( PhiM.' * (D0*PhiM) ...
                 + PhiM.' * (D1*G)   ...
                 - G.'    * (D2*G) );
end

% enforce symmetry and small numerical clean-up
Sigma = 0.5*(Sigma + Sigma.');
end


