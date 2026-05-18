function [ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = compute_modal( ...
    fem, aero, phi0_keep, Mred, Kred, Nm, aoa, chain, w0, debugPlot, rigid, NetworkPath)
%COMPUTE_PAZY_REDTEST  Build intrinsic-style structural ROM quantities (clean, consistent, keepDofs-based)
%
% This version is a full rework focused on:
%   - strict frame consistency (global vs element-local)
%   - correct element-vs-node indexing (Phi2/CPhi2 are ELEMENT fields)
%   - no accidental overwrites (Phi1 is never replaced by Phi0)
%   - clean data flow and precomputed mappings for speed
%   - minimal assumptions about fem/keepDofs shapes
%
% Inputs (expected):
%   fem.keepDofs          : [nKeep x 1] global dof IDs kept after condensation
%   fem.coordinates       : [nNode x 3]
%   fem.num_node          : scalar
%   fem.num_elem          : scalar (optional; only used for kappa lookup)
%   fem.structural_twist  : optional, used to build kappa (curvature) per element
%   fem.app_forces        : optional, [6 x nNode] external nodal forces/moments
%
%   phi0_keep : [nKeep x Nm] displacement/rotation modal vectors (condensed DOFs)
%   w0        : [Nm x 1] natural frequencies (rad/s) consistent with phi0 ordering
%
% Outputs:
%   ModeVars_discrete: fields in keepDOF space + element fields stacked (6*nElemTotal x Nm)
%   ModeVars_continuous: piecewise-linear segment structs for phi0, phi1 (node-based 6-vectors)
%   Beam_Props: M,K,E (E stacked as 6*nElemTotal x 6 in LOCAL element frames)
%   phi1_sA, phi2_sA: convenience interpolated values at "center" node A (see end)

% ----------------------------
% 0) Input hygiene / defaults
% ----------------------------
if nargin < 12 || isempty(NetworkPath), NetworkPath = false; end %#ok<NASGU>
if nargin < 11 || isempty(rigid),       rigid = 0;           end %#ok<NASGU>
if nargin < 10 || isempty(debugPlot),   debugPlot = 0;       end
if nargin <  8 || isempty(chain),       error('chain must be provided'); end
if nargin <  9 || isempty(w0),          error('w0 must be provided');    end

keepDofs = fem.keepDofs(:);
nKeep    = numel(keepDofs);
Clamped = 1;
if size(phi0_keep,1) ~= nKeep
    error('phi0_keep must be size [numel(fem.keepDofs) x Nm]. Got %dx%d, expected %dxNm.', ...
        size(phi0_keep,1), size(phi0_keep,2), nKeep);
end

if size(phi0_keep,2) < Nm
    error('phi0_keep has only %d modes but Nm=%d was requested.', size(phi0_keep,2), Nm);
end

phi0_keep = phi0_keep(:,1:Nm);
w0        = w0(:);
if numel(w0) < Nm
    error('w0 must have at least Nm entries. Got %d, Nm=%d.', numel(w0), Nm);
end
w0 = w0(1:Nm);

XYZ   = fem.coordinates;
nNode = fem.num_node;

% Optional: apply a global AoA rotation (about global Y) if you truly want to
% reinterpret the structural vectors in an AoA-tilted global frame.
applyAoA = ~isempty(aoa) && ~isnan(aoa) && abs(double(aoa)) > 0;
if applyAoA
    R_aoa = rotY(deg2rad(double(aoa)));
else
    R_aoa = eye(3);
end
R6_aoa = blkdiag(R_aoa, R_aoa);

% -----------------------------------------
% 1) Build a fast globalDOF -> keepIndex map
% -----------------------------------------
keepMap = buildKeepMap(keepDofs);

% ----------------------------------------------------------
% 2) Split chain into monotone-span segments (robust default)
% ----------------------------------------------------------
chainSegments = splitChainBySpan(chain(:), XYZ(:,2));
nSeg = numel(chainSegments);

% ----------------------------------------------------------
% 3) Precompute node-wise 6-vectors from keep-space arrays
%    (Phi0, Phi1, Psi1, and nodal force vectors for Phi2)
% ----------------------------------------------------------

% Phi1 in keep-space: Phi1 = Omega * Phi0 (reference-consistent)
Phi1_keep = phi0_keep .* (ones(nKeep,1) * (w0(:).'));  % [nKeep x Nm]

% Psi1 in keep-space: Psi1 = M * Phi1
Psi1_keep = Mred * Phi1_keep;                          % [nKeep x Nm]

% Nodal generalized forces in keep-space: F_keep = -K * Phi0
F_keep    = -(Kred * phi0_keep);                       % [nKeep x Nm]
% F_keep    = (Kred * phi0_keep);                       % [nKeep x Nm]
% F_keep    = -w0^2.*(Mred * phi0_keep);                       % [nKeep x Nm]
% fKAll(:,j) = -w0(j)^2*Mred*phi0_global(:,j);

% Map all node-based quantities into explicit [6 x nNode x Nm] arrays.
Phi0_node = zeros(6, nNode, Nm);
Phi1_node = zeros(6, nNode, Nm);
Psi1_node = zeros(6, nNode, Nm);
F_node    = zeros(6, nNode, Nm);

for nid = 1:nNode
    dofs6 = double((nid-1)*6) + (1:6);
    idx6  = keepMapLookup(keepMap, dofs6); % [1x6], 0 if missing

    present = (idx6 > 0);
    if any(present)
        Phi0_node(present, nid, :) = reshape(phi0_keep(idx6(present),:).', 1, [], Nm);
        Phi1_node(present, nid, :) = reshape(Phi1_keep(idx6(present),:).', 1, [], Nm);
        Psi1_node(present, nid, :) = reshape(Psi1_keep(idx6(present),:).', 1, [], Nm);
        F_node(present,  nid, :)   = reshape(F_keep(idx6(present),:).',   1, [], Nm);

        % np = sum(present);
        % 
        % Phi0_node(present, nid, 1:Nm) = reshape(phi0_keep(idx6(present), 1:Nm), [np, 1, Nm]);
        % Phi1_node(present, nid, 1:Nm) = reshape(Phi1_keep(idx6(present), 1:Nm), [np, 1, Nm]);
        % Psi1_node(present, nid, 1:Nm) = reshape(Psi1_keep(idx6(present), 1:Nm), [np, 1, Nm]);
        % F_node(present,  nid, 1:Nm)   = reshape(F_keep(  idx6(present), 1:Nm), [np, 1, Nm]);

    end
    
    % Apply optional AoA rotation in 6D (both translational and rotational parts)
    if applyAoA
        for k = 1:Nm
            Phi0_node(:,nid,k) = R6_aoa * Phi0_node(:,nid,k);
            Phi1_node(:,nid,k) = R6_aoa * Phi1_node(:,nid,k);
            Psi1_node(:,nid,k) = R6_aoa * Psi1_node(:,nid,k);
            F_node(:,nid,k)    = R6_aoa * F_node(:,nid,k);
        end
    end
end
fprintf('Phi1_node size = [%d %d %d]\n', size(Phi1_node,1), size(Phi1_node,2), size(Phi1_node,3));
fprintf('Psi1_node size = [%d %d %d]\n', size(Psi1_node,1), size(Psi1_node,2), size(Psi1_node,3));

% ----------------------------------------------------------
% 4) Build element-wise fields (GEOMETRY ONCE, then fill per mode)
% ----------------------------------------------------------

% Count total elements across all segments (each segment: length-1)
nElemTotal = 0;
for s = 1:nSeg
    nElemTotal = nElemTotal + (numel(chainSegments{s}) - 1);
end

% Preallocate element-level containers
Phi2_g   = zeros(6, nElemTotal, Nm);
Phi2_l   = zeros(6, nElemTotal, Nm);
CPhi2_l  = zeros(6, nElemTotal, Nm);
Phi1m_l  = zeros(6, nElemTotal, Nm);
Psi1m_l  = zeros(6, nElemTotal, Nm);

elemLen  = zeros(nElemTotal, 1);
elemMid  = zeros(nElemTotal, 3);
elemNA   = zeros(nElemTotal, 1);
elemNB   = zeros(nElemTotal, 1);

Eall_loc = zeros(6*nElemTotal, 6);     % stacked [6x6] per element in LOCAL frame
R3_all   = zeros(3,3,nElemTotal);      % local->global rotation per element

% Also store segment->element-index mapping so we can fill per-mode later
segElemStart = zeros(nSeg,1);
segElemCount = zeros(nSeg,1);

% ---- PASS A: build geometry and per-element local operators ONCE ----
eGlobal = 0;
for s = 1:nSeg
    segNodes = chainSegments{s};
    nSegNodes = numel(segNodes);
    if nSegNodes < 2
        segElemStart(s) = eGlobal + 1;
        segElemCount(s) = 0;
        continue;
    end

    segElemStart(s) = eGlobal + 1;
    segElemCount(s) = nSegNodes - 1;

    for i = 1:(nSegNodes-1)
        nA = segNodes(i);
        nB = segNodes(i+1);

        eGlobal = eGlobal + 1;

        rA   = XYZ(nA,:);
        rB   = XYZ(nB,:);
        rMid = 0.5*(rA + rB);
        ds   = norm(rB - rA);
        if ds < 1e-14
            error('Zero/near-zero element length between nodes %d and %d.', nA, nB);
        end

        elemLen(eGlobal)   = ds;
        elemMid(eGlobal,:) = rMid;
        elemNA(eGlobal)    = nA;
        elemNB(eGlobal)    = nB;

        % local->global frame
        R3 = buildSegmentFrame(rA, rB);
        R3_all(:,:,eGlobal) = R3;

        % LOCAL E operator for this element
        kappa_loc = getKappaLocal(fem, eGlobal, ds);   % 3x1
        E6_loc    = L1_operator([1;0;0; kappa_loc]);   % 6x6 in LOCAL frame
        Eall_loc((eGlobal-1)*6 + (1:6), :) = E6_loc;
    end
end

if eGlobal ~= nElemTotal
    error('Internal error: counted nElemTotal=%d but assembled eGlobal=%d.', nElemTotal, eGlobal);
end

% ---- PASS B: fill mode-dependent element quantities WITHOUT changing element index ----
for s = 1:nSeg
    segNodes   = chainSegments{s};
    nSegNodes  = numel(segNodes);
    if nSegNodes < 2
        continue;
    end

    eStart = segElemStart(s);

    % segment node positions as 3xN
    segXYZ = XYZ(segNodes,:);
    rSeg   = segXYZ.';  % [3 x nSegNodes]

    for kMode = 1:Nm
        % nodal force/moment at segment nodes (GLOBAL)
        Fk = squeeze(F_node(1:3, segNodes, kMode));   % [3 x nSegNodes]
        Mk = squeeze(F_node(4:6, segNodes, kMode));   % [3 x nSegNodes]

        % M0k = M + r x F
        M0k = Mk + cross3xN(rSeg, Fk);

        % tail sums from node index p..end (inclusive)
        Ftail  = fliplr(cumsum(fliplr(Fk),  2));      % [3 x nSegNodes]
        M0tail = fliplr(cumsum(fliplr(M0k), 2));      % [3 x nSegNodes]

        % element loop: element i is between segNodes(i) and segNodes(i+1)
        for i = 1:(nSegNodes-1)
            e  = eStart + (i-1);
            ds = elemLen(e);

            nA = segNodes(i);
            nB = segNodes(i+1);

            rMid = elemMid(e,:).';
            R3   = R3_all(:,:,e);
            R6   = blkdiag(R3, R3);

            % --- Phi2 (global), use tail starting at node (i+1) ---
            % internal force cut at element midpoint: sum of nodal loads downstream
            if (i+1) <= nSegNodes
                Fint = Ftail(:, i+1);
                Mint = M0tail(:, i+1) - cross(rMid, Fint);
            else
                Fint = zeros(3,1);
                Mint = zeros(3,1);
            end
            Phi2_g(:, e, kMode) = [Fint; Mint];
            Phi2_l(:, e, kMode) = R6.' * Phi2_g(:, e, kMode);

            % --- CPhi2 (LOCAL) from Phi0 finite difference + E^T midpoint ---
            Phi0A_g = squeeze(Phi0_node(:, nA, kMode));
            Phi0B_g = squeeze(Phi0_node(:, nB, kMode));
            Phi0A_l = R6.' * Phi0A_g;
            Phi0B_l = R6.' * Phi0B_g;

            E6_loc = Eall_loc((e-1)*6 + (1:6), :);     % LOCAL 6x6
            CPhi2_l(:, e, kMode) = ...
                -(Phi0B_l - Phi0A_l)/ds + (E6_loc.') * (Phi0B_l + Phi0A_l)/2;

            % --- Midpoint Phi1, Psi1 (LOCAL) ---
            Phi1A_g = squeeze(Phi1_node(:, nA, kMode));
            Phi1B_g = squeeze(Phi1_node(:, nB, kMode));
            Psi1A_g = squeeze(Psi1_node(:, nA, kMode));
            Psi1B_g = squeeze(Psi1_node(:, nB, kMode));

            Phi1m_l(:, e, kMode) = 0.5*(R6.'*Phi1A_g + R6.'*Phi1B_g);
            Psi1m_l(:, e, kMode) = 0.5*(R6.'*Psi1A_g + R6.'*Psi1B_g);
        end
    end
end


%%
phiRed = reshape(Phi0_node, size(Phi0_node,1)* size(Phi0_node,2), []);
% disp((Phi0_node))
writematrix(phiRed,'phiRed.dat','Delimiter','tab')
% writematrix(phi0_keep,'phi0_keep.dat','Delimiter','tab')
% writematrix(Phi1_keep,'Phi1_keep.dat','Delimiter','tab')
% error('wait here')


% ----------------------------------------------------------
% 5) Modal alpha normalizations (reference-consistent scalars)
%    alpha1: keep-space kinetic pairing (Phi1^T * Psi1)
%    alpha2: element local pairing sum(ds * Phi2_l^T * CPhi2_l)
% ----------------------------------------------------------

MPhi1_node = Psi1_node;

[Alpha1, Alpha2] = intrinsic_integrals_alphas( ...
    chainSegments, segElemStart, elemLen, ...
    Phi1_node, MPhi1_node, ...
    Phi2_l, CPhi2_l, ...
    Nm, Clamped);
% alpha1 = (Alpha1);
% alpha1 = diag(abs((Alpha1)));
alpha1 = diag(abs((Alpha1)));
% alpha2 = Alpha2;
alpha2 = diag(Alpha2);


% alpha1 = zeros(Nm,1);
% alpha2 = zeros(Nm,1);
% [alpha1, alpha2] = intrinsicAlphaNormalizationScalar( ...
%     keepDofs, Phi1_keep, Phi2_l, chain, ...
%     elemLen, Psi1_keep, CPhi2_l, Phi1m_l);
%% ============================================================
%  ALPHAS (match FEMINAS python: integral_alphas)
%  alpha1 += Phi1(node j+1)^T * MPhi1(node j+1)
%  alpha2 += Phi2(elem j)^T * CPhi2(elem j) * ds
%  + (if ~Clamped) add extra root-node term once at (beam1, node0)
%% ============================================================

% alpha1 = zeros(Nm, Nm);
% alpha2 = zeros(Nm, Nm);
% 
% for iB = 1:NumBeams
%     nElemB = BeamSeg(iB).EnumNodes - 1;
% 
%     for e = 1:nElemB
%         n2 = e + 1;                          % python j+1 node
%         ds = BeamSeg(iB).NodeDL(e);          % python NodeDL[j]
% 
%         Phi1_n2  = sliceNode6xNm(Phi1{iB},  n2);
%         MPhi1_n2 = sliceNode6xNm(MPhi1{iB}, n2);
% 
%         Phi2_e   = sliceElem6xNm(Phi2{iB},  e);
%         CPhi2_e  = sliceElem6xNm(CPhi2{iB}, e);
% 
%         % alpha1: no ds factor in python
%         alpha1 = alpha1 + (Phi1_n2.' * MPhi1_n2);
% 
%         % alpha2: ds factor in python
%         alpha2 = alpha2 + (Phi2_e.' * CPhi2_e) * ds;
% 
%         % unclamped root correction: only once (i==0 && j==0 in python)
%         if ~Clamped && iB == 1 && e == 1
%             Phi1_root  = sliceNode6xNm(Phi1{1},  1);   % python node 0
%             MPhi1_root = sliceNode6xNm(MPhi1{1}, 1);
%             error('ye boi erroor')
%             alpha1 = alpha1 + (Phi1_root.' * MPhi1_root);
%         end
%     end
% end

% (optional) diagnostics you can keep:
% fprintf('||alpha1 - I||F = %.3e\n', norm(alpha1 - eye(Nm),'fro'));
% fprintf('sym(alpha1)=%.3e, sym(alpha2)=%.3e\n', ...
%     norm(alpha1-alpha1.','fro')/max(1,norm(alpha1,'fro')), ...
%     norm(alpha2-alpha2.','fro')/max(1,norm(alpha2,'fro')));
% [alpha1, alpha2] = assemble_alphas_ref( ...
%     chainSegments, segElemStart, elemLen, ...
%     Phi1_keep, Psi1_keep, Phi1m_l, Phi2_l, CPhi2_l, ...
%     Nm, Clamped);
% % alpha1 in keep-space (fast, robust)
% alpha1 = sum(Phi1_keep .* Psi1_keep, 1).';  % [Nm x 1]
% 
% % alpha2 in element-local form
% for kMode = 1:Nm
%     acc = 0;
%     for e = 1:nElemTotal
%         acc = acc + elemLen(e) * (Phi2_l(:,e,kMode).' * CPhi2_l(:,e,kMode));
%     end
%     alpha2(kMode) = acc;
% end
% disp('alpha1')
% disp(alpha1)
% 
% disp('alpha2')
% disp(alpha2)

% Safety: ensure positive scalars (common if sign conventions are correct)
% If you see negatives here, you have a sign/frame inconsistency upstream.
if any(alpha1 <= 0)
    warning('Some alpha1 <= 0 detected. Check w0, Phi0 scaling, and Mred definiteness.');
end
if any(alpha2 <= 0)
    warning('Some alpha2 <= 0 detected. Check Phi2/CPhi2 sign conventions and local frames.');
end

% ----------------------------------------------------------
% 6) Apply normalization scalings (match reference behavior)
%    - scale Phi0/Phi1/Psi1 and midpoint Phi1m/Psi1m by 1/sqrt(alpha1)
%    - scale Phi2 and CPhi2 by 1/sqrt(alpha2)
% ----------------------------------------------------------
s1 = 1 ./ sqrt(abs(alpha1));   % [Nm x 1]
s2 = 1 ./ sqrt(abs(alpha2));   % [Nm x 1]

% keep-space scalings (Phi0, Phi1, Psi1)
% phi0_keep_n = phi0_keep .* (ones(nKeep,1) * (s1(:).'));
% disp(size(phi0_keep)); disp(size(ones(nKeep,1) * (s1(:).')));
% disp(size(w0));
phi0_keep_n =phi0_keep .* ((ones(nKeep,1)) * (s1(:).'));
% phi0_keep_n =repmat(w0,1, nKeep)'.*phi0_keep .* ((ones(nKeep,1)) * (s1(:).'));
Phi1_keep_n = Phi1_keep .* (ones(nKeep,1) * (s1(:).'));
Psi1_keep_n = Psi1_keep .* (ones(nKeep,1) * (s1(:).'));

% node-based phi0/phi1/psi1 (for continuous output & debugging)
Phi0_node_n = Phi0_node;
Phi1_node_n = Phi1_node;
Psi1_node_n = Psi1_node;
for kMode = 1:Nm
    Phi0_node_n(:,:,kMode) = Phi0_node(:,:,kMode) * s1(kMode);
    % Phi0_node_n(:,:,kMode) = w0(kMode)*Phi0_node(:,:,kMode) * s1(kMode);
    Phi1_node_n(:,:,kMode) = Phi1_node(:,:,kMode) * s1(kMode);
    Psi1_node_n(:,:,kMode) = Psi1_node(:,:,kMode) * s1(kMode);
end

% element-based scalings (Phi2, CPhi2, plus midpoint fields)
Phi2_g_n  = Phi2_g;
Phi2_l_n  = Phi2_l;
CPhi2_l_n = CPhi2_l;
Phi1m_l_n = Phi1m_l;
Psi1m_l_n = Psi1m_l;

for kMode = 1:Nm
    Phi2_g_n(:, :, kMode)  = Phi2_g(:, :, kMode)  * s2(kMode);
    Phi2_l_n(:, :, kMode)  = Phi2_l(:, :, kMode)  * s2(kMode);
    CPhi2_l_n(:, :, kMode) = CPhi2_l(:, :, kMode) * s2(kMode);

    % Midpoint fields are “alpha1-normalized”, not alpha2
    Phi1m_l_n(:, :, kMode) = Phi1m_l(:, :, kMode) * s1(kMode);
    Psi1m_l_n(:, :, kMode) = Psi1m_l(:, :, kMode) * s1(kMode);
end

% ----------------------------------------------------------
% 7) Post-check alphas (optional but recommended)
% ----------------------------------------------------------
% [alpha1_check, alpha2_check] = assemble_alphas_ref( ...
%     chainSegments, segElemStart, elemLen, ...
%     Phi1_keep_n, Psi1_node_n, Phi1m_l_n, Phi2_l_n, CPhi2_l_n, ...
%     Nm, Clamped);
% [alpha1_check, alpha2_check] = intrinsicAlphaNormalizationScalar( ...
%     keepDofs, Phi1_keep_n, Phi2_l_n, chain, ...
%     elemLen, Psi1_node_n, CPhi2_l_n, Phi1m_l_n);

% [alpha1_check, alpha2_check] = intrinsic_integrals_alphas( ...
%     chainSegments, segElemStart, elemLen, ...
%     Phi1_node_n, Psi1_node_n, ...
%     Phi2_l_n, CPhi2_l_n, ...
%     Nm, Clamped);
alpha1_check = (Phi1_keep_n).' * (Mred * Phi1_keep_n); % [Nm x Nm]
alpha1_check = 0.5*(alpha1_check + alpha1_check.');

alpha2_check = zeros(Nm,Nm);
for j = 1:Nm
    for k = 1:Nm
        acc = 0;
        for e = 1:nElemTotal
            acc = acc + elemLen(e) * (Phi2_l_n(:,e,j).' * CPhi2_l_n(:,e,k));
        end
        alpha2_check(j,k) = acc;
    end
end
alpha2_check = 0.5*(alpha2_check + alpha2_check.');

if max(abs(diag(alpha1_check) - 1)) > 1e-6
    warning('alpha1_check diagonal not ~1. max|diag-1|=%.3e', max(abs(diag(alpha1_check)-1)));
end
if max(abs(diag(alpha2_check) - 1)) > 1e-6
    warning('alpha2_check diagonal not ~1. max|diag-1|=%.3e', max(abs(diag(alpha2_check)-1)));
end

% ----------------------------------------------------------
% 8) Build Gamma tensors efficiently (element-wise, local frames)
%   Gamma1(j,k,l) += ds * Phi1m_j^T * L1(Phi1m_k) * Psi1m_l
%   Gamma2(j,k,l) += ds * Phi1m_j^T * L2(Phi2_k)  * CPhi2_l
% ----------------------------------------------------------
Gamma1 = zeros(Nm,Nm,Nm);
Gamma2 = zeros(Nm,Nm,Nm);

for e = 1:nElemTotal
    ds = elemLen(e);

    % Assemble per-element matrices: [6 x Nm]
    Phi1m = squeeze(Phi1m_l_n(:,e,:));   % 6 x Nm
    Psi1m = squeeze(Psi1m_l_n(:,e,:));   % 6 x Nm
    Phi2  = squeeze(Phi2_l_n(:,e,:));    % 6 x Nm
    CPhi2 = squeeze(CPhi2_l_n(:,e,:));   % 6 x Nm

    % For each k, build L1/L2 and update entire (j,l) plane via matrix mult.
    for k = 1:Nm
        L1k = L1_operator(Phi1m(:,k));         % 6x6
        A1  = (Phi1m.' * L1k);                 % Nm x 6
        G1jl = A1 * Psi1m;                     % Nm x Nm  (rows j, cols l)
        Gamma1(:,k,:) = reshape(squeeze(Gamma1(:,k,:)) + ds*G1jl, [Nm,1,Nm]);

        L2k = L2_operator(Phi2(:,k));          % 6x6
        A2  = (Phi1m.' * L2k);                 % Nm x 6
        G2jl = A2 * CPhi2;                     % Nm x Nm
        Gamma2(:,k,:) = reshape(squeeze(Gamma2(:,k,:)) + ds*G2jl, [Nm,1,Nm]);
    end
end
% [Gamma1,Gamma2] = assemble_gammas_ref( ...
%     chainSegments, segElemStart, elemLen, ...
%      Phi1_keep_n, Psi1_node_n, Phi1m_l_n, Phi2_l_n, CPhi2_l_n, ...
%     Nm, Clamped);
%% ============================================================
%  GAMMAS (match FEMINAS python: integral_gammas, quadratic=0)
%  gamma1 += Phi1(node j+1)^T * L1(Phi1(node j+1)) * MPhi1(node j+1)
%  gamma2 += Phi1m(elem j)^T * L2(Phi2(elem j)) * CPhi2(elem j) * ds
%  + (if ~Clamped) add extra root-node gamma1 term once at (beam1,node0)
% %% ============================================================
% MPhi1_node_n = Psi1_node_n;
% 
% [Gamma1, Gamma2] = intrinsic_integrals_gammas( ...
%     chainSegments, segElemStart, elemLen, ...
%     Phi1_node_n, MPhi1_node_n, ...
%     Phi1m_l_n, Phi2_l_n, CPhi2_l_n, ...
%     Nm, Clamped);

% gamma1 = zeros(Nm, Nm, Nm);
% gamma2 = zeros(Nm, Nm, Nm);
% 
% for iB = 1:NumBeams
%     nElemB = BeamSeg(iB).EnumNodes - 1;
% 
%     for e = 1:nElemB
%         n2 = e + 1;
%         ds = BeamSeg(iB).NodeDL(e);
% 
%         % node-based slices (python uses j+1)
%         Phi1_n2  = sliceNode6xNm(Phi1{iB},  n2);     % 6 x Nm
%         MPhi1_n2 = sliceNode6xNm(MPhi1{iB}, n2);     % 6 x Nm
% 
%         % element-based slices (python uses j)
%         Phi1m_e  = sliceElem6xNm(Phi1m{iB}, e);      % 6 x Nm
%         Phi2_e   = sliceElem6xNm(Phi2{iB},  e);      % 6 x Nm
%         CPhi2_e  = sliceElem6xNm(CPhi2{iB}, e);      % 6 x Nm
% 
%         % -------- gamma1 accumulation (no ds factor) --------
%         % For each k2: T(k1,k3) = Phi1(:,k1)' * L1(Phi1(:,k2)) * MPhi1(:,k3)
%         for k2 = 1:Nm
%             L1 = L1fun(Phi1_n2(:,k2));               % 6 x 6
%             T  = Phi1_n2.' * L1 * MPhi1_n2;          % Nm x Nm  (k1 x k3)
%             gamma1(:,k2,:) = gamma1(:,k2,:) + reshape(T, [Nm 1 Nm]);
%         end
% 
%         % unclamped root correction: only once (i==0 && j==0 in python)
%         if ~Clamped && iB == 1 && e == 1
%             Phi1_root  = sliceNode6xNm(Phi1{1},  1);
%             MPhi1_root = sliceNode6xNm(MPhi1{1}, 1);
% 
%             for k2 = 1:Nm
%                 L1 = L1fun(Phi1_root(:,k2));
%                 T  = Phi1_root.' * L1 * MPhi1_root;
%                 gamma1(:,k2,:) = gamma1(:,k2,:) + reshape(T, [Nm 1 Nm]);
%             end
%         end
% 
%         % -------- gamma2 accumulation (WITH ds factor) --------
%         for k2 = 1:Nm
%             L2 = L2fun(Phi2_e(:,k2));                % 6 x 6
%             T  = Phi1m_e.' * L2 * CPhi2_e;           % Nm x Nm  (k1 x k3)
%             gamma2(:,k2,:) = gamma2(:,k2,:) + reshape(T, [Nm 1 Nm]) * ds;
%         end
%     end
% end

% ----------------------------------------------------------
% 9) External force projection eta (keepDofs-consistent)
% ----------------------------------------------------------
etaVec = zeros(Nm,1);
if isfield(fem,'app_forces') && ~isempty(fem.app_forces)
    fext_keep = zeros(nKeep,1);
    fext_nodes = fem.app_forces; % expected [6 x nNode]
    if size(fext_nodes,1) ~= 6 || size(fext_nodes,2) ~= nNode
        warning('fem.app_forces expected size [6 x fem.num_node]. Skipping etaVec.');
    else
        % Fill only kept DOFs
        for nid = 1:nNode
            dofs6 = double((nid-1)*6) + (1:6);
            idx6  = keepMapLookup(keepMap, dofs6);
            present = idx6>0;
            if any(present)
                f6 = fext_nodes(:,nid);
                if applyAoA, f6 = R6_aoa * f6; end
                fext_keep(idx6(present)) = fext_keep(idx6(present)) + f6(present);
            end
        end
        etaVec = Phi1_keep_n.' * fext_keep;
    end
end

% ----------------------------------------------------------
% 10) Kelvin-Voigt sigma (element-local consistent form)
% ----------------------------------------------------------
cfg.struct.damping_model = 'zeta';
cfg.struct.zeta = 0.005; % 0.5% critical, uniform (adjust as needed)
Sigma = rom_build_sigma_kv_elementlocal(Phi1_node_n, elemNA, elemNB, elemLen, R3_all, Eall_loc, cfg);
% Sigma = rom_build_sigma_kv(Phi1_node_n, chain, elemLen, Eall_loc, cfg);

% ----------------------------------------------------------
% 11) Build continuous piecewise-linear representations
% ----------------------------------------------------------
Phi0Fun = buildPiecewiseLinearFun(chainSegments, XYZ, Phi0_node_n);
Phi1Fun = buildPiecewiseLinearFun(chainSegments, XYZ, Phi1_node_n);

% ----------------------------------------------------------
% 12) Convenience: interpolate at “center node” sA (if root is clamped)
% ----------------------------------------------------------
phi1_sA = zeros(6,Nm);
phi2_sA = zeros(6,Nm);

if nNode >= 2
    sA = XYZ(1,2);
    sL = XYZ(2,2);
    sR = XYZ(end,2);
    eLen = (sR - sL);
    if abs(eLen) < 1e-14
        % fallback: just take node 2
        phi1_sA = squeeze(Phi1_node_n(:,2,:));
    else
        phi1L = squeeze(Phi1_node_n(:,2,:));      % 6 x Nm
        phi1R = squeeze(Phi1_node_n(:,end,:));    % 6 x Nm
        phi1_sA = phi1L * ((sR - sA)/eLen) + phi1R * ((sA - sL)/eLen);
    end

    % For Phi2, approximate using first element of first segment and last element
    if nElemTotal >= 2
        phi2L = squeeze(Phi2_g_n(:,1,:));              % 6 x Nm
        phi2R = squeeze(Phi2_g_n(:,nElemTotal,:));     % 6 x Nm
        if abs(eLen) < 1e-14
            phi2_sA = phi2L;
        else
            phi2_sA = phi2L * ((sR - sA)/eLen) + phi2R * ((sA - sL)/eLen);
        end
    end
end

% ----------------------------------------------------------
% 13) Pack outputs (clean, explicit)
% ----------------------------------------------------------
ModeVars_discrete = struct();

ModeVars_discrete.keepDofs     = keepDofs;
ModeVars_discrete.w0           = w0;

ModeVars_discrete.phi0_local    = phi0_keep_n;
% ModeVars_discrete.phi0_local    = phi0_keep;
ModeVars_discrete.phi1_local    = Phi1_keep_n;
ModeVars_discrete.phi2_local    = Psi1_keep_n;

% Element fields stacked as [6*nElemTotal x Nm] in LOCAL frames (recommended for ROM terms)
ModeVars_discrete.elemNA       = elemNA;
ModeVars_discrete.elemNB       = elemNB;
ModeVars_discrete.elemMid      = elemMid;
ModeVars_discrete.elemLen      = elemLen;
ModeVars_discrete.Eall_loc     = Eall_loc;   % stacked local E6 per element
ModeVars_discrete.R3_all       = R3_all;     % local->global per element

ModeVars_discrete.phi2_elem_g  = reshape(Phi2_g_n,  6*nElemTotal, Nm);
ModeVars_discrete.phi2_elem_l  = reshape(Phi2_l_n,  6*nElemTotal, Nm);
ModeVars_discrete.psi2_elem_l  = reshape(CPhi2_l_n, 6*nElemTotal, Nm);  % CPhi2 in local

ModeVars_discrete.alpha1       = alpha1;
ModeVars_discrete.alpha2       = alpha2;
ModeVars_discrete.alpha1_check = alpha1_check;
ModeVars_discrete.alpha2_check = alpha2_check;

% ModeVars_discrete.etaVecExt    = etaVec;
ModeVars_discrete.etaVecExt    = etaVec;
ModeVars_discrete.Sigma        = Sigma;

ModeVars_discrete.Gamma1       = Gamma1;
ModeVars_discrete.Gamma2       = Gamma2;

ModeVars_continuous = struct();
ModeVars_continuous.phi0 = Phi0Fun;
ModeVars_continuous.phi1 = Phi1Fun;

Beam_Props = struct();
Beam_Props.M = Mred;
Beam_Props.K = Kred;
Beam_Props.E = Eall_loc; % LOCAL element frames

% ----------------------------------------------------------
% 14) Optional debug plots (lightweight)
% ----------------------------------------------------------
if debugPlot
    try
        figure('Name','alpha checks','Color','w');
        subplot(2,1,1); plot(diag(alpha1_check),'-o'); grid on; ylabel('diag(\alpha_1)');
        subplot(2,1,2); plot(diag(alpha2_check),'-o'); grid on; ylabel('diag(\alpha_2)'); xlabel('mode');

        % quick Phi2 force components for mode 1 (in local frame)
        modePlot = 1;
        Fxyz = zeros(nElemTotal,3);
        for e = 1:nElemTotal
            Fxyz(e,:) = Phi2_l_n(1:3,e,modePlot).';
        end
        figure('Name','Phi2 local forces (mode 1)','Color','w');
        plot(Fxyz,'-o'); grid on; xlabel('element'); ylabel('force'); legend('Fx','Fy','Fz');
    catch ME
        warning('Debug plot failed: %s', ME.message);
    end
end

end % main function

% =========================================================================
%                               HELPERS
% =========================================================================

function keepMap = buildKeepMap(keepDofs)
% Map global dof id -> local keep index (0 if absent)
maxD = max(keepDofs);
keepMap = zeros(maxD,1,'int32');
keepMap(keepDofs) = int32(1:numel(keepDofs));
end

function idx6 = keepMapLookup(keepMap, dofs6)
idx6 = zeros(1,numel(dofs6),'int32');
for i = 1:numel(dofs6)
    d = dofs6(i);
    if d >= 1 && d <= numel(keepMap)
        idx6(i) = keepMap(d);
    else
        idx6(i) = 0;
    end
end
end

function segs = splitChainBySpan(chain, y)
% Split when span direction changes sign (robust for full-wing chains)
if isempty(chain)
    segs = {};
    return;
end
yc = y(chain);
d  = diff(yc);

% ignore zeros for sign detection
dz = d(abs(d) > 1e-12);
if isempty(dz)
    segs = {chain};
    return;
end
sgn0 = sign(dz(1));

breakIdx = [];
for i = 1:numel(d)
    if abs(d(i)) < 1e-12
        continue;
    end
    if sign(d(i)) ~= sgn0
        breakIdx = i;
        break;
    end
end

if isempty(breakIdx)
    segs = {chain};
else
    segs = {chain(1:breakIdx), chain(breakIdx+1:end)};
end
end

function R3 = buildSegmentFrame(rA, rB)
% local x-axis along tangent rB-rA; build right-handed frame
ex = (rB - rA).';
L = norm(ex);
ex = ex / (L + 1e-16);

eyTrial = [0;1;0];
if abs(dot(ex,eyTrial)) > 0.99
    eyTrial = [0;0;1];
end
ez = cross(ex, eyTrial);
ez = ez / (norm(ez) + 1e-16);
ey = cross(ez, ex);

R3 = [ex, ey, ez]; % columns are local axes in global coords
end

function kappa = getKappaLocal(fem, elemIdx, ds)
% Robust curvature extraction. Returns 3x1 in LOCAL element coordinates.
kappa = zeros(3,1);
if ~isfield(fem,'structural_twist') || isempty(fem.structural_twist)
    return;
end
st = fem.structural_twist;

% Accept several common shapes:
%   - [3 x num_elem] : intrinsic curvature integrated over element
%   - [num_elem x 3]
%   - [num_elem x 1] or [1 x num_elem] : scalar twist about x (convert)
if isnumeric(st)
    if ndims(st) == 2
        if size(st,1) == 3
            if size(st,2) >= elemIdx
                kappa = st(:,elemIdx) / ds;
            end
        elseif size(st,2) == 3
            if size(st,1) >= elemIdx
                kappa = st(elemIdx,:).' / ds;
            end
        elseif isvector(st)
            if numel(st) >= elemIdx
                kappa = [st(elemIdx)/ds; 0; 0];
            end
        end
    end
end
end

function M6 = L1_operator(x6)
% L1 operator for intrinsic beam (6x6)
v = x6(1:3);
w = x6(4:6);
M6 = zeros(6);
M6(1:3,1:3) = skew(w);
M6(4:6,1:3) = skew(v);
M6(4:6,4:6) = skew(w);
end

function M6 = L2_operator(x6)
% L2 operator for intrinsic beam (6x6)
f = x6(1:3);
m = x6(4:6);
M6 = zeros(6);
M6(1:3,4:6) = skew(f);
M6(4:6,1:3) = skew(f);
M6(4:6,4:6) = skew(m);
end

function S = skew(u)
S = [   0   -u(3)  u(2)
      u(3)    0   -u(1)
     -u(2)  u(1)    0 ];
end

function C = cross3xN(A, B)
% A,B are 3xN. Return C = A x B column-wise.
C = [ A(2,:).*B(3,:) - A(3,:).*B(2,:);
      A(3,:).*B(1,:) - A(1,:).*B(3,:);
      A(1,:).*B(2,:) - A(2,:).*B(1,:) ];
end

function R = rotY(th)
c = cos(th); s = sin(th);
R = [ c  0  s
      0  1  0
     -s  0  c ];
end

function Fun = buildPiecewiseLinearFun(chainSegments, XYZ, Phi_node)
% Build piecewise-linear segment structs for each mode:
% Fun{m}.seg(q) has fields:
%   .nodeA, .nodeB, .xrange [sA sB], .valA [6x1], .valB [6x1]
Nm = size(Phi_node,3);
Fun = cell(1,Nm);

% Pre-count segments
segList = [];
for s = 1:numel(chainSegments)
    segNodes = chainSegments{s};
    for i = 1:(numel(segNodes)-1)
        segList = [segList; segNodes(i), segNodes(i+1)]; %#ok<AGROW>
    end
end
nSegTotal = size(segList,1);

for m = 1:Nm
    segArr = repmat(struct('nodeA',[],'nodeB',[],'xrange',[],'valA',[],'valB',[]), nSegTotal, 1);
    for q = 1:nSegTotal
        nA = segList(q,1);
        nB = segList(q,2);
        segArr(q).nodeA  = nA;
        segArr(q).nodeB  = nB;
        segArr(q).xrange = [XYZ(nA,2), XYZ(nB,2)];
        segArr(q).valA   = Phi_node(:,nA,m);
        segArr(q).valB   = Phi_node(:,nB,m);
    end
    Fun{m}.seg = segArr;
end
end

function Sigma = rom_build_sigma_kv_elementlocal(Phi1_node, elemNA, elemNB, elemLen, R3_all, Eall_loc, cfg)
% Element-local Kelvin-Voigt sigma consistent with element local frames.
% Phi1_node: [6 x nNode x Nm] node velocities in GLOBAL frame already normalized.
% We rotate node values into each element-local frame using R6' (global->local).

Nm   = size(Phi1_node,3);
nEl  = numel(elemLen);
Sigma = zeros(Nm,Nm);

I6 = eye(6);
has_fun = isfield(cfg.struct,'Cd_fun') && ~isempty(cfg.struct.Cd_fun);

for e = 1:nEl
    ds = elemLen(e);
    if ds < 1e-14, continue; end

    R3 = R3_all(:,:,e);
    R6 = blkdiag(R3,R3);

    E6  = Eall_loc((e-1)*6+(1:6),:);

    if has_fun
        Cd = cfg.struct.Cd_fun(e);
    elseif isfield(cfg.struct,'zeta') && ~isempty(cfg.struct.zeta)
        z = cfg.struct.zeta;
        if isscalar(z), Cd = 2*z*I6; else, Cd = 2*z; end
    else
        Cd = 2*(5e-3)*I6;
    end

    D0 = -E6*Cd*E6.';        % 6x6
    D1 =  E6*Cd - Cd*E6.';   % 6x6
    D2 =  Cd;                % 6x6

    nA = elemNA(e);
    nB = elemNB(e);

    % Local node values for all modes: [6 x Nm]
    PhiA_l = zeros(6,Nm);
    PhiB_l = zeros(6,Nm);
    for k = 1:Nm
        PhiA_l(:,k) = R6.' * Phi1_node(:,nA,k);
        PhiB_l(:,k) = R6.' * Phi1_node(:,nB,k);
    end

    PhiM = 0.5*(PhiA_l + PhiB_l);
    G    = (PhiB_l - PhiA_l)/ds;

    Sigma = Sigma + ds * ( PhiM.'*(D0*PhiM) + PhiM.'*(D1*G) - G.'*(D2*G) );
end

Sigma = 0.5*(Sigma + Sigma.');
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

%% Gammas

function [Gamma1, Gamma2] = assemble_gammas_ref( ...
    chainSegments, segElemStart, elemLen, ...
    Phi1_node_l, MPhi1_node_l, Phi1m_l, Phi2_l, CPhi2_l, ...
    Nm, Clamped)


Gamma1 = zeros(Nm,Nm,Nm);   % Gamma1(k1,k2,k3)
Gamma2 = zeros(Nm,Nm,Nm);   % Gamma2(k1,k2,k3)

nSeg = numel(chainSegments);

% ---- root-node add-on for UNCLAMPED (matches: if not Clamped and i==0 and j==0) ----
if ~Clamped
    rootNode = chainSegments{1}(1);

    P  = reshape(Phi1_node_l(:, rootNode, 1:Nm),  6, Nm);  % 6xNm
    MP = reshape(MPhi1_node_l(:, rootNode, 1:Nm), 6, Nm);  % 6xNm


    for k2 = 1:Nm
        L1 = L1_operator(P(:,k2));          % 6x6  ( intrinsic.functions.L1fun)
        T  = P.' * (L1 * MP);               % Nm x Nm, rows=k1 cols=k3
        Gamma1(:,k2,:) = Gamma1(:,k2,:) + reshape(T, [Nm,1,Nm]);
    end
end

% ---- beam-wise loops (match Python beams and j+1 node usage) ----
for b = 1:nSeg
    segNodes = chainSegments{b};
    nN = numel(segNodes);
    if nN < 2, continue; end

    eStart = segElemStart(b);

    for j = 1:(nN-1)
        nB = segNodes(j+1);          % Python uses [j+1] node for alpha1/gamma1
        e  = eStart + (j-1);         % element index for alpha2/gamma2
        ds = elemLen(e);

        % --- alpha1/gamma1 node contribution (NO ds) ---
        P  = reshape(Phi1_node_l(:, nB, 1:Nm),  6, Nm);
        MP = reshape(MPhi1_node_l(:, nB, 1:Nm), 6, Nm);


        for k2 = 1:Nm
            L1 = L1_operator(P(:,k2));
            T  = P.' * (L1 * MP);            % (k1,k3)
            Gamma1(:,k2,:) = Gamma1(:,k2,:) + reshape(T, [Nm,1,Nm]);
        end

        % --- alpha2/gamma2 element contribution (WITH ds) ---
        Pm  = reshape(Phi1m_l(:, e, 1:Nm), 6, Nm);
        P2  = reshape(Phi2_l(:,  e, 1:Nm), 6, Nm);
        CP2 = reshape(CPhi2_l(:, e, 1:Nm), 6, Nm);


        for k2 = 1:Nm
            L2 = L2_operator(P2(:,k2));      %  intrinsic.functions.L2fun
            T  = Pm.' * (L2 * CP2);          % (k1,k3)
            Gamma2(:,k2,:) = Gamma2(:,k2,:) + reshape(T, [Nm,1,Nm]) * ds;
        end
    end
end

end
%% Alphas

function [Alpha1, Alpha2] = assemble_alphas_ref( ...
    chainSegments, segElemStart, elemLen, ...
    Phi1_node_l, MPhi1_node_l, Phi1m_l, Phi2_l, CPhi2_l, ...
    Nm, Clamped)

Alpha1 = zeros(Nm,Nm);
Alpha2 = zeros(Nm,Nm);


nSeg = numel(chainSegments);

% ---- root-node add-on for UNCLAMPED (matches: if not Clamped and i==0 and j==0) ----
if ~Clamped
    rootNode = chainSegments{1}(1);

    P  = reshape(Phi1_node_l(:, rootNode, 1:Nm),  6, Nm);  % 6xNm
    MP = reshape(MPhi1_node_l(:, rootNode, 1:Nm), 6, Nm);  % 6xNm

    Alpha1 = Alpha1 + (P.' * MP);

    
end

% ---- beam-wise loops (match Python beams and j+1 node usage) ----
for b = 1:nSeg
    segNodes = chainSegments{b};
    nN = numel(segNodes);
    if nN < 2, continue; end

    eStart = segElemStart(b);

    for j = 1:(nN-1)
        nB = segNodes(j+1);          % Python uses [j+1] node for alpha1/gamma1
        e  = eStart + (j-1);         % element index for alpha2/gamma2
        ds = elemLen(e);
        disp(size(Phi1_node_l))
        disp(size(nB))
        % --- alpha1/gamma1 node contribution (NO ds) ---
        P  = reshape(Phi1_node_l(:, nB, 1:Nm),  6, Nm);
        MP = reshape(MPhi1_node_l(:, nB, 1:Nm), 6, Nm);

        Alpha1 = Alpha1 + (P.' * MP);

        

        % --- alpha2/gamma2 element contribution (WITH ds) ---
        Pm  = reshape(Phi1m_l(:, e, 1:Nm), 6, Nm);
        P2  = reshape(Phi2_l(:,  e, 1:Nm), 6, Nm);
        CP2 = reshape(CPhi2_l(:, e, 1:Nm), 6, Nm);

        Alpha2 = Alpha2 + (P2.' * CP2) * ds;

    end
end

end

function [alpha1, alpha2] = intrinsicAlphaNormalizationScalar( ...
    keepDofs, phi1_loc, phi2_loc, chain, ...
    elemChainL, Psi1_discrete, psi2_discrete, phi1_half)

% We'll build alpha1, alpha2 => [nModes x nModes]

% elemLen_sub1 = elemLen(:,1);
% elemLen_sub2 = elemLen(:,2);
% elemLen = [elemLen_sub1;elemLen_sub2]


Nm = size(phi1_loc, 2);             % number of modes
Nn = length(chain);           % number of nodes
Ne = Nn - 1;                    % number of elements
disp(chain)

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
function [Alpha1, Alpha2] = intrinsic_integrals_alphas( ...
    chainSegments, segElemStart, elemLen, ...
    Phi1_node, MPhi1_node, ...
    Phi2_elem, CPhi2_elem, ...
    Nm, Clamped)
%INTRINSIC_INTEGRALS_ALPHAS  Exact MATLAB translation of FEMINAS integral_alphas
%
% Inputs:
%   chainSegments : {nSeg x 1} cell, each is vector of node IDs (global node indices)
%   segElemStart  : [nSeg x 1] starting global element index for each segment (1-based)
%   elemLen       : [nElemTotal x 1] element lengths (ds)
%   Phi1_node     : [6 x nNode x Nm]
%   MPhi1_node    : [6 x nNode x Nm]  (this corresponds to python MPhi1)
%   Phi2_elem     : [6 x nElemTotal x Nm]
%   CPhi2_elem    : [6 x nElemTotal x Nm]
%
% Outputs:
%   Alpha1, Alpha2 : [Nm x Nm]

nModesPhi = size(Phi1_node, 3);
if nModesPhi < Nm
    error(['intrinsic_integrals_alphas: Phi1_node has only %d mode slice(s) (size(Phi1_node,3)), ' ...
           'but Nm=%d was requested. Fix upstream node packing so Phi1_node is [6 x nNode x Nm].'], ...
           nModesPhi, Nm);
end

Alpha1 = zeros(Nm, Nm);
Alpha2 = zeros(Nm, Nm);

nSeg = numel(chainSegments);

% ---- Python: if not Clamped and i==0 and j==0 => add node 0 term once ----
if ~Clamped
    rootNode = chainSegments{1}(1);  % first node of first segment
    P  = reshape(Phi1_node(:, rootNode, 1:Nm),  6, Nm);
    MP = reshape(MPhi1_node(:, rootNode, 1:Nm), 6, Nm);
    Alpha1 = Alpha1 + (P.' * MP);
end

% ---- Beam/segment loop: j = 0..EnumNodes-2 -> MATLAB j = 1..(nN-1) ----
for b = 1:nSeg
    segNodes = chainSegments{b};
    nN = numel(segNodes);
    if nN < 2, continue; end

    eStart = segElemStart(b);

    for j = 1:(nN-1)
        nB = segNodes(j+1);      % python uses [j+1] node
        e  = eStart + (j-1);     % python uses element index [j]
        ds = elemLen(e);

        % alpha1 contribution (NO ds in python)
        P  = reshape(Phi1_node(:,  nB, 1:Nm),  6, Nm);
        MP = reshape(MPhi1_node(:, nB, 1:Nm),  6, Nm);
        Alpha1 = Alpha1 + (P.' * MP);

        % alpha2 contribution (WITH ds in python)
        P2  = reshape(Phi2_elem(:,  e, 1:Nm),  6, Nm);
        CP2 = reshape(CPhi2_elem(:, e, 1:Nm),  6, Nm);
        Alpha2 = Alpha2 + (P2.' * CP2) * ds;
    end
end

end
function [Gamma1, Gamma2] = intrinsic_integrals_gammas( ...
    chainSegments, segElemStart, elemLen, ...
    Phi1_node, MPhi1_node, ...
    Phi1m_elem, Phi2_elem, CPhi2_elem, ...
    Nm, Clamped)
%INTRINSIC_INTEGRALS_GAMMAS  Exact MATLAB translation of FEMINAS integral_gammas (quadratic=0)
%
% Uses L1_operator and L2_operator already defined in this file.
%
% Outputs:
%   Gamma1, Gamma2 : [Nm x Nm x Nm] with indexing Gamma(:,k2,:) matching python gamma(:,k2,:)

Gamma1 = zeros(Nm, Nm, Nm);
Gamma2 = zeros(Nm, Nm, Nm);

nSeg = numel(chainSegments);

% ---- Python: if not Clamped and i==0 and j==0 => add node 0 term once ----
if ~Clamped
    rootNode = chainSegments{1}(1);
    P  = reshape(Phi1_node(:,  rootNode, 1:Nm),  6, Nm);
    MP = reshape(MPhi1_node(:, rootNode, 1:Nm),  6, Nm);

    for k2 = 1:Nm
        L1 = L1_operator(P(:,k2));        % python intrinsic.functions.L1fun
        T  = P.' * (L1 * MP);             % Nm x Nm (rows k1, cols k3)
        Gamma1(:,k2,:) = Gamma1(:,k2,:) + reshape(T, [Nm, 1, Nm]);
    end
end

% ---- Beam/segment loop ----
for b = 1:nSeg
    segNodes = chainSegments{b};
    nN = numel(segNodes);
    if nN < 2, continue; end

    eStart = segElemStart(b);

    for j = 1:(nN-1)
        nB = segNodes(j+1);          % python Phi1[..., j+1, :]
        e  = eStart + (j-1);         % python Phi1m/Phi2/CPhi2[..., j, :]
        ds = elemLen(e);             % python NodeDL[j]

        % -------- gamma1 (node-based, NO ds) --------
        P  = reshape(Phi1_node(:,  nB, 1:Nm),  6, Nm);
        MP = reshape(MPhi1_node(:, nB, 1:Nm),  6, Nm);

        for k2 = 1:Nm
            L1 = L1_operator(P(:,k2));
            T  = P.' * (L1 * MP);
            Gamma1(:,k2,:) = Gamma1(:,k2,:) + reshape(T, [Nm, 1, Nm]);
        end

        % -------- gamma2 (element-based, WITH ds) --------
        Pm  = reshape(Phi1m_elem(:, e, 1:Nm),   6, Nm);  % Phi1m[i][k1][j,:]
        P2  = reshape(Phi2_elem(:,  e, 1:Nm),   6, Nm);  % Phi2[i][k2][j,:]
        CP2 = reshape(CPhi2_elem(:, e, 1:Nm),   6, Nm);  % CPhi2[i][k3][j,:]

        for k2 = 1:Nm
            L2 = L2_operator(P2(:,k2));         % python intrinsic.functions.L2fun
            T  = Pm.' * (L2 * CP2);             % Nm x Nm (rows k1, cols k3)
            Gamma2(:,k2,:) = Gamma2(:,k2,:) + reshape(T, [Nm, 1, Nm]) * ds;
        end
    end
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