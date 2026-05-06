%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstration code that:
%   A. Twist is rotation about x-axis in rad 
%   B. Sweep is rotation about z-axis in rad
%   C. aoa is rotation about y-axis in rad
%  (1) Computes xi_bar(s) by solving eq. (2.36) or eq. (3.47) for a
%      baseline curvature kappa_0(s).
%  (2) Extracts xi_bar_v, the vector part.
%  (3) Builds phi_xi_j(s) from eq. (3.49)-(3.50) using 
%      U_t(xi_bar(s)) * phi_omega_j(s).
%  (4) Illustrates "iterative" usage: if you have a new kappa_0(s)
%      for a re-linearization step, you re-run step (1).
%
% The "xi_bar(s)" is NOT time-dependent in the usual sense: 
%    it's the baseline/undeformed (or equilibrium) rotation field. 
% But in a control or re-linearization loop, you can re-solve 
%    xi_bar(s) after each major iteration if the geometry changed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clamped node is assumed undeformed and identity quat [1 0 0 0]
function [phi_xi_modes, Gamma_xi_mat, Gamma_g_3D, phi_omega, xi_bar_array] =solve_baseline_quaternion_xi(fem, aero , ...
    phi1, Nm, M_global, aoa, phi0)

    
    % For full wing:
    idxMinL = min(fem.coordinates(:,2));
    idxMaxL = max(fem.coordinates(:,2));
    L = idxMaxL- idxMinL;
    discFact = 1; % for finer grid. Do 1 for same as SHARPy case file.
    nPts = discFact*fem.num_node;   % # of points for discretization
    sVals = linspace(idxMinL,idxMaxL,nPts).'; % 0 for half wing -L for full wing
    beamSegments(1).s0 = idxMinL;
    beamSegments(1).s1 = idxMaxL;
    % % For half wing
    % L = L/2;
    % sVals = linspace(0,L,nPts).'; % 0 for half wing -L for full wing
    gVec = [0; 0;-9.81];  % demonstration, in local coords if needed

    
    % shape nPts x 3, each row is [kappa_x, kappa_y, kappa_z].
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
        disp('keke')
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
    
    % Rotation about Y-axis (for all nodes)
    aoa_s = deg2rad(aoa)*ones(nPts, 1);
    


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
    phi_omega = zeros(size(phi1,1)/2,Nm);
    kk =1;
    del_pts = 0; del_idx = [];
    for i=1:nPts
        sNow = sVals(i);
        alpha_s = alpha_of_s(i);  % your function for sweep vs. s
        theta_s = twist_of_s(i);  % your function for twist vs. s
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
       phi_omega(dofOmega,:) = phi1(idx_k,:);
       % phi_omega(dofOmega,:) = phi0(idx_k,:);
    end
    xi_bar_array(del_idx,:) = []; % delete clamped nodes
    coord_del = fem.coordinates(del_idx, 2);
    idx_coord_del = find(sVals == coord_del);
    sVals(idx_coord_del) = [];

    xi_bar_v = xi_bar_array(:, 2:4);
    
    % This solves eq. (3.47) for an "undeformed" 
    
    xi_bar_prep = zeros((nPts-del_pts)*4,1);
    phi_xi_modes = zeros((nPts-del_pts)*4,Nm+1);
    for j=1:Nm
        for i=1:(nPts-del_pts)
            xiBar_pt = xi_bar_array(i,:).';   % 4x1
            dofCur = double((i-1)*3)+(1:3);
            w_j      = phi_omega(dofCur,j);  % 3x1 
            UT = halfTangentialOperator( xiBar_pt ); % eq. (3.49) or eq. (3.50)
            Idx_pxi = double((i-1)*4) + (1:4);
            xi_bar_prep(Idx_pxi) = xiBar_pt;
            phi_xi_modes(Idx_pxi,j+1) = UT * w_j;
        end
    end
    % writematrix(UT,'UT.dat','Delimiter','tab');
    % 
    % writematrix(phi_xi_modes,'phi_xi_mPre.dat','Delimiter','tab');
    % writematrix(xi_bar_prep,'xi_bar_prep.dat','Delimiter','tab');
% phi_xi_modes = 1/2*phi_xi_modes;
    phi_xi_modes(:,1) = xi_bar_prep; % Augmented approach
    % phi_xi_modes(:,1) = .675*xi_bar_prep; % Augmented approach
    % phi_xi_modes(:,1) = repmat([1;0;0;0], size(xi_bar_prep,1)/4,1); % Augmented approach
    % if trimAlg %fir trim we only need to update alpha
    %     Gamma_xi_mat =[]; Gamma_g_3D = []; 
    %     return;
    % end
    % phi_xi0(s)= xi_bar_array(s). phi_xi_j for j=1..Nm, and j=0 => phi_xi0= xi_bar.
    
% phi_xi_modes = 1/2*phi_xi_modes;

    % writematrix(phi_xi_modes,'phi_xi_modes.dat','Delimiter','tab');

    %% (D) Compute A_xi from eq. (3.53):
    %   [A_xi]_{ij} = sum_{m=1}^{n_b} int_{0}^{L_{b_m}} phi_{xi i}^T(s) * phi_{xi j}(s) ds

    A_xi = zeros(Nm+1, Nm+1);
    for i=1:(Nm+1)
        for j=1:(Nm+1)
            % integrate dot-product of phi_xi_i(s), phi_xi_j(s) over s
            A_xi(i,j) = integrate_phi_xi_dot( phi_xi_modes(:,i), phi_xi_modes(:,j), sVals, beamSegments);
        end
    end
    % writematrix(A_xi,'A_xi.dat');

    % disp('Matrix A_xi from eq. (3.53):');
    % disp(A_xi);

    Nm_plus_1 = Nm+1;
    
    % Build Gamma_xi:
    % size (Nm_plus_1, Nm_plus_1, Nm). 
    % For l=1..Nm, eq.:
    %   [Gamma_xi^l](i,j) = integral( phi_xi_i^T(s)* U( phi_omega_l(s) ) phi_xi_j(s) ds ).
    
% %%% OLD ATTEMPT FOR A 21x20x21
%    % Gamma_xi = zeros(Nm_plus_1, Nm_plus_1, Nm);
%     % Gamma_xi = zeros(Nm_plus_1, Nm, Nm);
%     Gamma_xi = zeros(Nm_plus_1, Nm, Nm+1);
%     Gamma_xi_mat = Gamma_xi;
%     % for l = 1:Nm+1
%     %     Gamma_xi(:,:,l) = compute_Gamma_xi_mode( l, ...
%     %         phi_xi_modes, phi_omega, sVals, beamSegments, Nm, Nm_plus_1 );
%     %     Gamma_xi_mat(:,:,l) = A_xi\Gamma_xi(:,:,l);
%     % end
    UOmeg_xi_mat = cell(Nm,1);
    Gamma_xi_cell = cell(Nm,1);
%     for l = 1:Nm+1
%     % for l = 1:Nm
%         Gamma_xi(:,:,l) = compute_Gamma_xi_mode( l, ...
%             phi_xi_modes, phi_omega, sVals, beamSegments, Nm, Nm_plus_1 );
% 
%         Gamma_xi_mat(:,:,l) = A_xi\Gamma_xi(:,:,l);
% 
%         % For debug:
%         UOmeg_xi_mat{l} = Gamma_xi(:,:,l);
%         Gamma_xi_cell{l} = Gamma_xi_mat(:,:,l);
% 
%         % disp(num2str(l))
%         % writematrix(Gamma_xi_mat(:,:,l),"Gamma_xi_mat"+ int2str(l)+".dat");
%         % writematrix(Gamma_xi(:,:,l),"UOmeg_xi_mat"+ int2str(l)+".dat");
% 
%     end

%%% NEW ATTEMPT FOR A 21x21x20
    % Gamma_xi = zeros(Nm_plus_1, Nm, Nm+1);
    Gamma_xi = zeros(Nm_plus_1, Nm_plus_1, Nm);
    Gamma_xi_mat = Gamma_xi;
    % for l = 1:Nm+1
    %     Gamma_xi(:,:,l) = compute_Gamma_xi_mode( l, ...
    %         phi_xi_modes, phi_omega, sVals, beamSegments, Nm, Nm_plus_1 );
    %     Gamma_xi_mat(:,:,l) = A_xi\Gamma_xi(:,:,l);
    % end
    UOmeg_xi_mat = cell(Nm_plus_1,1);
    Gamma_xi_cell = cell(Nm_plus_1,1);
    % for l = 1:Nm+1
    for l = 1:Nm
        Gamma_xi(:,:,l) = compute_Gamma_xi_mode( l, ...
            phi_xi_modes, phi_omega, sVals, beamSegments, Nm, Nm_plus_1 );
        
        Gamma_xi_mat(:,:,l) = A_xi\Gamma_xi(:,:,l);

        % For debug:
        UOmeg_xi_mat{l} = Gamma_xi(:,:,l);
        Gamma_xi_cell{l} = Gamma_xi_mat(:,:,l);

        % disp(num2str(l))
        % writematrix(Gamma_xi_mat(:,:,l),"Gamma_xi_mat"+ int2str(l)+".dat");
        % writematrix(Gamma_xi(:,:,l),"UOmeg_xi_mat"+ int2str(l)+".dat");
        
    end
    % writematrix(Gamma_xi_mat,'Gamma_xi_mat.dat');
    % disp(Gamma_xi_cell);  disp(size(Gamma_xi_cell)); 
    % writecell(Gamma_xi_cell,'Gamma_xi_cell.xls');
    % writecell(UOmeg_xi_mat,'UOmeg_xi_mat.xls');

    
    % ikoo = inv(A_xi);
    % dwd = ikoo*Gamma_xi2
    % Build Gamma_g:
    % eq. (3.42)-like expression. 
    
    % Gamma_g_2D = compute_Gamma_g( phi1, phi_xi_modes, sVals, beamSegments, Nm, gVec );
    % dimension => [Nm x Nm_plus_1]. Each row i => structural mode i, each column j => quaternion mode j.
    % Gamma_g_3D = zeros(Nm, Nm_plus_1, Nm);
    Gamma_g_3D = zeros(Nm, Nm_plus_1, Nm_plus_1);
    % M_global = Beam_Props.M;
    % Triple loop:
    % i in [1..Nm], j in [1..Nm+1], l in [1..Nm].
    for i = 1:Nm
        for j = 1:Nm_plus_1
            for l = 1:Nm+1
                val_ijl = integrate_Gamma_g_ijl( i, j, l, ...
                                 phi1, phi_xi_modes, sVals, beamSegments, gVec, M_global );
                Gamma_g_3D(i,j,l) = val_ijl;
            end
        end
    end
      


    % %%% (D) Illustrate an iterative approach for control:
 
    
    % new_kappa0Vals = kappa0Vals + 0.01;  % e.g. a small shift
    % xi_bar_array_new = integrate_xi_bar_ODE(sVals, new_kappa0Vals);
    % 
    % % and you'd build new phi_xi_j from that new baseline.
    % 
    % disp('Done building phi_xi_j(s) and Gamma');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function q_tot = combineSweepTwist(alpha, theta, aoa)
    % combineSweepTwist:
    %  Returns quaternion for "rotate by alpha about z, then by theta about x."
    %
  

    q_sweep = [cos(alpha/2),  0,            0,           sin(alpha/2)];
    q_twist = [cos(theta/2),  sin(theta/2), 0,           0           ];
    q_aoa = [cos(aoa/2), 0,  sin(aoa/2),           0           ];

    q_tot = quatMultiply(q_aoa, q_sweep, q_twist);
    % q_tot = quatMultiply(q_notessweep,q_aoa, q_twist);
end

function q_out = quatMultiply(qAOA, qA, qB)
    % Standard Hamilton product of quaternions qA, qB
    %   each in form [q0, q1, q2, q3].
    sAOA = qAOA(1); vAOA = qAOA(2:4);
    sA = qA(1); vA = qA(2:4);
    sB = qB(1); vB = qB(2:4);

    sAB = sA*sB - dot(vA,vB);
    vAB = sA*vB + sB*vA + cross(vA,vB);
    qAB = [sAB; vAB.'];

    % final rotation about angle of attack (y-axis)
    % Step 2: Final product q_out = qAOA * qAB
    sAOA = qAOA(1); vAOA = qAOA(2:4);
    
    sOut = sAOA*sAB - dot(vAOA, vAB);
    vOut = sAOA*vAB + sAB*vAOA + cross(vAOA, vAB);
    q_out = [sOut; vOut.'];

    % Normalize for safety
    q_out = q_out / norm(q_out);
end

function idxLocal = localDofIndex(keepDofs, dofSet)
idxLocal = [];
for dd = dofSet(:).'
   iPos = find( keepDofs==dd, 1, 'first' );
   if ~isempty(iPos)
       idxLocal(end+1) = iPos; 
   end
end
end

% function xi_bar_array = integrate_xi_bar_ODE(sVals, kappaVals)
% % Solve the ODE:
% %    d/ds xi_bar(s) = U( kappaVals(s) ) * xi_bar(s),
% % with boundary condition xi_bar(0)=[1,0,0,0].
% %
% 
%   nPts = length(sVals);
%   xi_bar_array = zeros(nPts,4);
%   % initial condition at s=0:
%   xi_bar_array(1,:) = [1,0,0,0];
%   for i=1:(nPts-1)
%     ds = sVals(i+1)-sVals(i);
%     kappa_mid = 0.5*(kappaVals(i,:) + kappaVals(i+1,:));  % average kappa
%     xi_now    = xi_bar_array(i,:).';
%     % U(kappa_mid) is 4x4 from eq. (2.36) style:
%     Umat      = bigU_of_kappa(kappa_mid);
%     % forward Euler:
%     xi_next   = xi_now + ds*( Umat*xi_now );
%     % normalize to avoid drift:
%     nq = norm(xi_next);
%     if nq<1e-12
%        xi_next = [1;0;0;0];
%     else
%        xi_next = xi_next./nq;
%     end
%     xi_bar_array(i+1,:) = xi_next.';
%   end
% end


% function U4 = bigU_of_kappa(kappa)
% % eq. (2.36) for the baseline:
% %   d/ds xi = 1/2 [ 0, -kappa^T; kappa, -~kappa ] xi
%   w=kappa(:);
%   U4= zeros(4,4);
%   U4(1,2:4)= -w.';
%   U4(2:4,1)=  w;
%   Wsk = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
%   U4(2:4,2:4) = -Wsk;
%   U4 = 0.5*U4;
% end


function UT = halfTangentialOperator( xibar )
% eq. (3.49) or eq. (3.50) from the reference:
%  U_t(xi_bar) = partial_x [ U(x)*xi_bar ] but the text states a closed form:
%  U_t( xi_bar ) = 1/2 [ -xi_bar_v^T;  xi_bar_0 I3 + tilde(xi_bar_v ) ]
%
% So we define:
  xi0 = xibar(1);
  xiv = xibar(2:4);
  UT = zeros(4,3);
  UT(1,1:3) = - (xiv.');
  skewv = [0 -xiv(3) xiv(2); xiv(3) 0 -xiv(1); -xiv(2) xiv(1) 0];
  UT(2:4,1:3) = xi0*eye(3) + skewv;
  UT = 0.5*UT;
end


function val = integrate_phi_xi_dot(phiA, phiB, sVals, beamSegments)
% integrate over the beam segments the dot-product (phiA^T phiB).
% phiA, phiB: [nPoints x 4],  sVals: [1 x nPoints].
% assume uniform spacing.
    val = 0;
    nPt = length(sVals);
    for i=1:length(beamSegments)
       s0 = beamSegments(i).s0;
       s1 = beamSegments(i).s1;
       % naive loop. 
       segmentVal=0;
       for k=1:(nPt-1)
         sMid = 0.5*(sVals(k)+sVals(k+1));
         if (sMid>=s0 && sMid<=s1)

          rowK    = 4*(k-1) + (1:4);
          rowKp1  = 4*(k)   + (1:4);
          % Extract the 4-dofs at node k
          a_k  = phiA(rowK);
          b_k  = phiB(rowK);
          dotVal = dot(a_k, b_k);
          % Extract the 4-dofs at node k+1
          a_kp1 = phiA(rowKp1);
          b_kp1 = phiB(rowKp1);
          dotVal2 = dot(a_kp1, b_kp1);

          ds = sVals(k+1) - sVals(k);
          segmentVal = segmentVal + 0.5*(dotVal + dotVal2)*ds;
           % dotVal = dot(phiA(k,:), phiB(k,:));
           % dotVal2= dot(phiA(k+1,:), phiB(k+1,:));
           % segmentVal = segmentVal + 0.5*(dotVal+dotVal2)*(sVals(k+1)-sVals(k));
         end
       end
       val = val + segmentVal;
    end
end


function Gxl = compute_Gamma_xi_mode(lMode, phi_xi_modes, phi_omega, sVals, beamSegments, Nm, Nm_plus_1)
% compute_Gamma_xi_mode:
%   Builds the matrix Gamma_xi^l from eq.:
%     [Gamma_xi^l]( i, j ) = < phi_xi_i, U( phi_omega_l ) phi_xi_j >_{L^2}.
%
%  assume:
%   - phi_xi_modes is [4*nNodes x (Nm_plus_1)] 
%   col i => phi_xi_i
%   - phi_omega is [3*nNodes x Nm], col lMode => phi_omega_l
%   - i, j in [1..(Nm+1)] => indexing the columns of phi_xi_modes
%   trapezoid or segment sum integral.
    %% Current Jan 2026
    % Gxl = zeros(Nm_plus_1, Nm);
    Gxl = zeros(Nm_plus_1, Nm_plus_1);

    % Extract the single rotational mode "phi_omega_l" from the big array:


    for kIdx = 1:Nm_plus_1
        for jIdx = 1:Nm_plus_1 % is the Nm or Nmplus 1
        % for jIdx = 1:Nm % is the Nm or Nmplus 1

        phi_omega_l_col = phi_omega(:, lMode);  % dimension [3*nNodes, 1]

            % integral over s: 
            val_ij = integrate_GammaXi_ij( ...
                phi_xi_modes(:, kIdx), ...  % entire column i
                phi_omega_l_col, ...
                phi_xi_modes(:, jIdx), ...
                sVals, beamSegments ...
                );
            Gxl(kIdx, jIdx) = val_ij;
            % Gxl(iIdx, jIdx-1) = val_ij;
        end
    end   

    %% Legacy 2025

    % Gxl = zeros(Nm_plus_1, Nm);
    % % Gxl = zeros(Nm_plus_1, Nm_plus_1);
    % 
    % % Extract the single rotational mode "phi_omega_l" from the big array:
    % 
    % 
    % for kIdx = 1:Nm_plus_1
    %     % for jIdx = 1:Nm_plus_1 % is the Nm or Nmplus 1
    %     for jIdx = 1:Nm % is the Nm or Nmplus 1
    %     % for jIdx = 1:Nm  % is the Nm or Nmplus 1
    %         phi_omega_j_col = phi_omega(:, jIdx);  % dimension [3*nNodes, 1]
    % 
    %         % integral over s: 
    %         val_ij = integrate_GammaXi_ij( ...
    %             phi_xi_modes(:, kIdx), ...  % entire column i
    %             phi_omega_j_col, ...
    %             phi_xi_modes(:, lMode), ...
    %             sVals, beamSegments ...
    %             );
    %         Gxl(kIdx, jIdx) = val_ij;
    %         % Gxl(iIdx, jIdx-1) = val_ij;
    %     end
    % end   
    %%
    % phi_omega_l_col = phi_omega(:, lMode);  % dimension [3*nNodes, 1]
    % 
    % for kIdx = 1:Nm_plus_1
    %     % for jIdx = 1:Nm_plus_1 % is the Nm or Nmplus 1
    %     % for jIdx = 1:Nm % is the Nm or Nmplus 1
    %     for jIdx = 1:Nm+1 % is the Nm or Nmplus 1
    %         % We'll define an integral over s: 
    %         val_ij = integrate_GammaXi_ij( ...
    %             phi_xi_modes(:, kIdx), ...  % entire column i
    %             phi_omega_l_col, ...
    %             phi_xi_modes(:, jIdx), ...
    %             sVals, beamSegments ...
    %             );
    %         Gxl(kIdx, jIdx) = val_ij;
    %         % Gxl(iIdx, jIdx-1) = val_ij;
    %     end
    % end
end


function val = integrate_GammaXi_ij(phi_xi_i, phi_omega_l, phi_xi_j, sVals, beamSegments)
% integrate for eq.:
%   \int phi_xi_i^T(s) * U( phi_omega_l(s) ) * phi_xi_j(s) ds

    val = 0;
    nNodes = length(sVals);
    for segID = 1:length(beamSegments)
        s0 = beamSegments(segID).s0;
        s1 = beamSegments(segID).s1;
        segVal = 0;
        for k=1:(nNodes-1)
            sMid = 0.5*(sVals(k) + sVals(k+1));
            if sMid>=s0 && sMid<=s1
               ds = sVals(k+1)-sVals(k);

               % Node k => extract the 4 dofs for phi_xi_i(k), phi_xi_j(k), 
               % and 3 dofs for phi_omega_l(k).
               rowXi_k    = 4*(k-1)+(1:4);
               rowOmega_k = 3*(k-1)+(1:3);

               a_i_k  = phi_xi_i(rowXi_k);
               a_j_k  = phi_xi_j(rowXi_k);
               w_k    = phi_omega_l(rowOmega_k);

               dotVal_k = a_i_k.' * bigU_of_omega(w_k) * a_j_k;

               % Node k+1 => similarly
               rowXi_kp1    = 4*(k)+(1:4);
               rowOmega_kp1 = 3*(k)+(1:3);

               a_i_kp1  = phi_xi_i(rowXi_kp1);
               a_j_kp1  = phi_xi_j(rowXi_kp1);
               w_kp1    = phi_omega_l(rowOmega_kp1);
               bigU = bigU_of_omega(w_kp1);
               dotVal_kp1 = a_i_kp1.' * bigU_of_omega(w_kp1) * a_j_kp1;

               % segVal = segVal + 0.5*(dotVal_k + dotVal_kp1)* ds;
               segVal = segVal + 0.5*(dotVal_k + dotVal_kp1)* ds;
            end
        end
        val = val + segVal;
    end
end

%%
% function Gamma_g_mat = compute_Gamma_g(phi_1, phi_xi_modes, sVals, beamSegments, Nm, gravVec)
% % eq. references:
% %   [Gamma_g]( i, j ) = < phi_{1 i}, mu*[I_3; ~r_cm]* T_phi(phi_xi_j, phi_xi_l)^T g >
% % but we keep it simpler and do a standard distribution approach. 
% % We'll produce a [Nm x (Nm+1)] matrix for demonstration, each entry integrated.
% %
% % In practice, you might have lumps or the local "r_cm" for each node. 
% % We'll just show the skeleton:
% 
%   nNodes = length(sVals);
%   Gamma_g_mat = zeros(Nm, Nm+1);
% 
%   % We'll define a (dummy) function that, for node k, returns Mv_k = mass at that node
%   % and r_cm_k = center of mass offset if you have that. For demonstration, we store it in 
%   % "massProps" or define in-line. We'll also define T_phi( a, b ) from eq. (32).
% 
%   for i = 1:Nm   % structural velocity mode i
%     for j = 1:(Nm+1)  % quaternion mode j
%       val_ij = integrate_Gamma_g_ij( i, j, phi_1, phi_xi_modes, sVals, beamSegments, gravVec );
%       Gamma_g_mat(i,j) = val_ij;
%     end
%   end
% end


% function val = integrate_Gamma_g_ij(iStruct, jQuat, phi_1, phi_xi_modes, sVals, beamSegments, gravVec)
% % integrand = phi_{1 iStruct}^T(s) * [mass*[I_3; ~r_cm]] * T_phi( phi_xi_j, phi_xi_l)^T * g
% % We'll define T_phi( a, b ) => eq. (32)
% % Here, we only pass one "jQuat" at a time; the reference also uses 'l' but in some eq. 
% % you might sum over l. We'll show the direct usage with a single j for simplicity.
% 
%   val = 0;
%   nNodes = length(sVals);
%   for segID=1:length(beamSegments)
%     s0 = beamSegments(segID).s0;
%     s1 = beamSegments(segID).s1;
%     segVal = 0;
%     for k=1:(nNodes-1)
%        sMid = 0.5*(sVals(k)+sVals(k+1));
%        if sMid>=s0 && sMid<=s1
%           ds = sVals(k+1)-sVals(k);
% 
%           % Extract dofs at node k for phi_1 iStruct => 6 dofs, 
%           %  but we only need the first 3 for translational or the full if you like
%           row_1_k = 6*(k-1)+(1:6);
%           a_1_k   = phi_1(row_1_k, iStruct);
% 
%           % The quaternion dofs for jQuat => 4 dofs => rowXi
%           rowXi_k = 4*(k-1)+(1:4);
%           q_xi_k  = phi_xi_modes(rowXi_k, jQuat);
% 
%           % Possibly we sum over another l, but let's keep it simple:
%           % e.g. T_phi( q_xi_k, q_xi_k ) or some eq. 
%           % The reference eq. might do T_phi( xi^{(a)}, xi^{(b)} ), or 
%           % might just do T_phi( phi_xi_j, phi_xi_l ) for a sum. 
%           % We'll do T_phi( q_xi_k, q_xi_k ) as a demonstration if the reference is eq. (32).
%           Tmat = T_phi_quat( q_xi_k, q_xi_k );
% 
%           % local mass or inertia => see if you define a function getMassMatrix(k)
%           Mv_local = eye(6); % for demonstration
%           % Typically something like:
%           %   Mv_local = [ mass_k*eye(3); mass_k*some ~r; ... ]
%           % We'll do a minimal approach: just use mass*g 
%           mass_k = 0.1; % example
% 
%           integrand_k = a_1_k(1:3).' * (mass_k * Tmat.' * gravVec);
% 
%           % node k+1 similarly if you want trapezoid
%           row_1_kp1= 6*(k)+(1:6);
%           a_1_kp1 = phi_1(row_1_kp1, iStruct);
%           rowXi_kp1= 4*(k)+(1:4);
%           q_xi_kp1= phi_xi_modes(rowXi_kp1, jQuat);
% 
%           Tmat_kp1   = T_phi_quat(q_xi_kp1, q_xi_kp1);
%           integrand_kp1 = a_1_kp1(1:3).' * (mass_k * Tmat_kp1.' * gravVec);
% 
%           segVal = segVal + 0.5*(integrand_k + integrand_kp1)* ds;
%        end
%     end
%     val = val + segVal;
%   end
% end
% 
% function T_AB = T_phi_quat(xiA, xiB)
% % eq. (32) style from the reference:
% %  T_phi(xi^{(a)}, xi^{(b)})= -2 (xi_v^(a))^T xi_v^(b)* I_3 + 2 xi_v^(a) (xi_v^(b))^T + ...
% %                             + 2 xi_0^(a) ~xi_v^(b).
% % Let xiA= [xi0_A; xiV_A], xiB= [xi0_B; xiV_B].
%   xi0A = xiA(1);
%   xvA  = xiA(2:4);
%   xi0B = xiB(1);
%   xvB  = xiB(2:4);
% 
%   dotVV = xvA.' * xvB;  % scalar
%   crossMat = skew3(xvB);
%   T_AB = (xi0A*xi0B - dotVV)*eye(3) + 2*(xvA*xvB.') + 2*xi0A*crossMat; 
% end
function T_AB = T_phi_quat(xiA, xiB)
  xi0A = xiA(1); xvA = xiA(2:4);
  xi0B = xiB(1); xvB = xiB(2:4);

  T_AB = (xi0A*xi0B - (xvA.'*xvB))*eye(3) ...
       + (xvA*xvB.' + xvB*xvA.') ...
       + xi0A*skew3(xvB) + xi0B*skew3(xvA);
end

function S = skew3(v)
  S=[0, -v(3), v(2);
     v(3), 0, -v(1);
    -v(2), v(1), 0];
end

function U4 = bigU_of_omega(omega3)
% eq. (2.36)-like operator for a 3-vector -> 4x4 block
  w = omega3(:);
  U4= zeros(4,4);
  U4(1,2:4)= -w.';
  U4(2:4,1)=  w;
  Wsk = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
  U4(2:4,2:4) = -Wsk;
  U4= 0.5*U4;
end

function val = integrate_Gamma_g_ijl(iStruct, jQuat, lQuat, ...
                             phi_1, phi_xi, sVals, beamSegments, gVec, M_global)
% integrate the expression from eq. (3.42)/(32)-type references. 
% Each node k => 
%   integrand ~ phi_{1 iStruct}(k)^T * M_local(...) * T_phi( phi_xi_j(k), phi_xi_l(k) )^T * gVec
% Summation or trapezoid across beamSegments and nodes.
%
% dimension logic:
%   iStruct => [1..Nm], picks out col iStruct from phi_1
%   jQuat   => [1..(Nm+1)], picks out col jQuat from phi_xi
%   lQuat   => [1..Nm], also from phi_xi (some references exclude the 0th mode, or not).
  val = 0;
  nNodes = length(sVals);

  for segID = 1:length(beamSegments)
    s0 = beamSegments(segID).s0;
    s1 = beamSegments(segID).s1;

    segVal = 0;
    for k = 1:(nNodes-1)
      sMid = 0.5*(sVals(k)+sVals(k+1));
      if (sMid >= s0 && sMid<=s1)
          ds = sVals(k+1)-sVals(k);
          
          % node k:
          row_1k  = 6*(k-1)+(1:6);   % 6 dofs for phi_1
          row_qjk = 4*(k-1)+(1:4);   % 4 dofs for phi_xi j
          row_qlk = 4*(k-1)+(1:4);   % 4 dofs for phi_xi l

          a_1_k  = phi_1(row_1k, iStruct);  % 6x1
          xi_j_k = phi_xi(row_qjk, jQuat);  % 4x1
          xi_l_k = phi_xi(row_qlk, lQuat);  % 4x1

          integrand_k = localGravityTerm(a_1_k, xi_j_k, xi_l_k, M_global, k, gVec);

          % node k+1:
          row_1kp1 = 6*k+(1:6);
          row_qjkp1= 4*k+(1:4);
          row_qlkp1= 4*k+(1:4);

          a_1_kp1  = phi_1(row_1kp1, iStruct);
          xi_j_kp1 = phi_xi(row_qjkp1, jQuat);
          xi_l_kp1 = phi_xi(row_qlkp1, lQuat);

          integrand_kp1 = localGravityTerm(a_1_kp1, xi_j_kp1, xi_l_kp1, M_global, k+1, gVec);

          segVal = segVal + 0.5*( integrand_k + integrand_kp1 )* ds;
      end
    end
    val = val + segVal;
  end
end


% function valk = localGravityTerm( phi1_i, xi_j, xi_l, mass_k, gVec )
% % eq. (32)-like:
% %   T_phi( xi_j, xi_l ) = -2 (xi_jv^T xi_lv) I_3 + 2 xi_jv xi_lv^T + 2 xi_j0 ~xi_lv
% % Then multiply by gVec, etc. 
% % We'll do a minimal approach:
%   xv_j = xi_j(2:4);
%   xv_l = xi_l(2:4);
%   xi0_j= xi_j(1);
%   dotVV = xv_j.'*xv_l;
% 
%   Tmat = -2*dotVV*eye(3) + 2*(xv_j*xv_l.') + 2*xi0_j*skew3(xv_l);
% 
%   % We assume phi1_i(1:3) => translational part
%   % so integrand ~ [phi1_i(1:3)]^T * (mass_k * Tmat^T ) * gVec
%   % or eq. (3.42) might have Tmat^T or Tmat. Check sign. We'll do Tmat^T for eq. 
%   Ttemp = Tmat.' * gVec; 
%   valk = (phi1_i(1:3).') * (mass_k * Ttemp);  
% end

function valk = localGravityTerm( phi1_i, xi_j, xi_l, M_global, nodeK, gVec )
% eq. (32)-like but now using the 6x3 block from M_global that maps 
% translational velocity => [linear+angular momentum].
%
% We'll do:
%   Mv_k = M_global( rowRange, colRange ), 
%   integrand = phi1_i^T * ( Mv_k * T_phi(...)^T*gVec )
% 
% T_phi( xi_j, xi_l ) = ...
% 
% Then we only take the "translational velocity" portion from phi1_i if needed?
% Actually, we can use all 6 dofs because Mv_k is 6x3. That yields a 6x1 multiplied by 6x3 => 3x1 => scalar.
%

  % 1) build T_phi matrix from eq. (32)
  Tmat = T_phi_quat( xi_j, xi_l );

  % 2) extract Mv_k from M_global for nodeK
  rowStart = 6*(nodeK-1)+1;  rowEnd = rowStart+5;
  colStart = rowStart;       colEnd = colStart+2; % first 3 columns
  Mv_k = M_global(rowStart:rowEnd, colStart:colEnd);  % 6x3

  % 3) "phi1_i" is 6x1 => the full velocity dofs [v(1..3), omega(1..3)].
  % multiply => phi1_i' * Mv_k => [1x3].
  % Then => Tmat^T*gVec => 3x1. 
  % final => scalar.

  temp31 = Mv_k * ( Tmat.' * gVec );    % [6x3]*[3x1] => [6x1]? Wait. 
  % Actually, Mv_k is [6x3], Tmat^T*gVec => [3x1], so product => [6x1].

  valk = (phi1_i.') * temp31;          % [1x6]*[6x1] => scalar
end

% 
