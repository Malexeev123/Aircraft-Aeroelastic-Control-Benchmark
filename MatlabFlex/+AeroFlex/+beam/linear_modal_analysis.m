function [M,K,Vr,w] = linear_modal_analysis(fem, keepDofs)
% Assemble simple consistent M,K for uniform beam or reuse yours.
% Here we assume you already have Mfull,Kfull from your pipeline; stub:

[Mfull, Kfull] = your_assembly_function_for_uniform_beam(fem); %#ok<NASGU> 
% Replace ↑ with your current M,K assembly path

% reduce to active dofs
M = Mfull(keepDofs, keepDofs);
K = Kfull(keepDofs, keepDofs);

% eigen
[Vr, Lam] = eigs(K, M, 10, 'SM');
w = sqrt(diag(Lam));
end

function [Mfull, Kfull] = your_assembly_function_for_uniform_beam(fem)
% YOUR_ASSEMBLY_FUNCTION_FOR_UNIFORM_BEAM
% Assemble global mass (Mfull) and stiffness (Kfull) matrices for a
% uniform, straight Euler–Bernoulli 3D beam mesh with 6 DOF per node:
% [ux, uy, uz, rx, ry, rz].
%
% Inputs (fields expected in "fem"):
%   fem.nodes  : (nnode x 3) node coordinates [x y z]
%   fem.conn   : (nelem x 2)  1-based node indices for each 2-node element
%   fem.E      : Young's modulus           (scalar or nelem-vector)
%   fem.G      : Shear modulus             (scalar or nelem-vector)
%   fem.A      : Cross-section area        (scalar or nelem-vector)
%   fem.Iy     : Second moment about local y (scalar or nelem-vector)
%   fem.Iz     : Second moment about local z (scalar or nelem-vector)
%   fem.J      : Saint-Venant torsion const (scalar or nelem-vector)
%   fem.rho    : Density                   (scalar or nelem-vector)
%
% Output:
%   Mfull, Kfull : (6*nnode x 6*nnode) sparse, assembled global matrices
%
% Notes:
% - Euler–Bernoulli (no shear deformation); good for slender-beam tests.
% - Local DOF order per element (12 dof): 
%   [u1 v1 w1 rx1 ry1 rz1, u2 v2 w2 rx2 ry2 rz2] in the *local* frame
% - Each element is rotated into global with a robust triad built from its
%   end points (local x along element; y,z chosen orthonormal).
% - Apply BCs outside (e.g., via your keepDofs mask) just like you do now.

% nodes = fem.nodes;
nodes = fem.coordinates;
conn  = fem.connectivities;

ne = size(conn,1);
nn = size(nodes,1);
ndof = 6*nn;

% Expand per-element properties if scalars were given
E   = expand_prop(fem, 'E',   ne);
G   = expand_prop(fem, 'G',   ne);
A   = expand_prop(fem, 'A',   ne);
Iy  = expand_prop(fem, 'Iy',  ne);
Iz  = expand_prop(fem, 'Iz',  ne);
J   = expand_prop(fem, 'J',   ne);
rho = expand_prop(fem, 'rho', ne);

% Preallocate sparse triplets (rough upper bound: each elem contributes
% up to 12x12 nonzeros)
est = ne * 12 * 12;
I = zeros(est,1); Jc = zeros(est,1); Kv = zeros(est,1);
Im = zeros(est,1); Jm = zeros(est,1); Mv = zeros(est,1);
tK = 0; tM = 0;

for e = 1:ne
    n1 = conn(e,1);
    n2 = conn(e,2);
    x1 = nodes(n1,:).';
    x2 = nodes(n2,:).';
    Le = norm(x2 - x1);
    if Le <= eps
        error('Element %d has zero length.', e);
    end

    % Local material/section props for this element
    Ee   = E(e);  Ge = G(e);  Ae = A(e);
    Iye  = Iy(e); Ize = Iz(e); Je = J(e);
    rhoe = rho(e);

    % Local element matrices (12x12) in the element's own frame
    [Ke_loc, Me_loc] = elem_matrices_euler_bernoulli_3d(Ee, Ge, Ae, Iye, Ize, Je, rhoe, Le);

    % Build orthonormal local frame R: ex along element; ey, ez orthonormal
    R = local_frame(x1, x2);  % 3x3

    % Transform to global: T = blkdiag(R,R,R,R)
    T = blkdiag(R,R,R,R);
    Ke = T.' * Ke_loc * T;
    Me = T.' * Me_loc * T;

    % Assembly indices for the two nodes
    id1 = node_dofs(n1);
    id2 = node_dofs(n2);
    idx = [id1 id2];

    % Accumulate into sparse triplets
    [ii, jj, vv] = find(Ke);
    kN = numel(vv);
    I(tK+(1:kN))  = idx(ii);
    Jc(tK+(1:kN)) = idx(jj);
    Kv(tK+(1:kN)) = vv;
    tK = tK + kN;

    [ii, jj, vv] = find(Me);
    mN = numel(vv);
    Im(tM+(1:mN))  = idx(ii);
    Jm(tM+(1:mN))  = idx(jj);
    Mv(tM+(1:mN))  = vv;
    tM = tM + mN;
end

% Build sparse matrices
Kfull = sparse(I(1:tK),  Jc(1:tK), Kv(1:tK), ndof, ndof);
Mfull = sparse(Im(1:tM), Jm(1:tM), Mv(1:tM), ndof, ndof);

end

% ------------------------ helpers ------------------------

function arr = expand_prop(fem, name, ne)
    v = fem.(name);
    if isscalar(v), arr = repmat(v, ne, 1);
    elseif numel(v)==ne, arr = v(:);
    else, error('fem.%s must be scalar or length ne.', name);
    end
end

function dofs = node_dofs(n)
% 6 dofs per node
    base = 6*(n-1);
    dofs = base + (1:6);
end

function R = local_frame(x1, x2)
% Returns a 3x3 rotation matrix whose columns are the local basis
% [ex ey ez] expressed in global coordinates.
    ex = (x2 - x1) / norm(x2 - x1);
    % Choose a "not parallel" up vector to build a frame
    up = [0;0;1];
    if abs(dot(ex, up)) > 0.95
        up = [0;1;0];
    end
    ez = cross(ex, up); ez = ez / norm(ez);
    ey = cross(ez, ex); % already orthonormal
    R = [ex, ey, ez];
end

function [Ke, Me] = elem_matrices_euler_bernoulli_3d(E, G, A, Iy, Iz, J, rho, L)
% 12x12 local Euler–Bernoulli beam element (no shear deformation).
% DOF order (local): [u v w rx ry rz, u2 v2 w2 rx2 ry2 rz2]

Ke = zeros(12,12);
Me = zeros(12,12);

% --- Axial (ux) ---
Ka = E*A/L * [ 1 -1; -1  1 ];
Ma = rho*A*L/6 * [ 2  1;  1  2 ];
idx_u = [1, 7];
Ke(idx_u, idx_u) = Ke(idx_u, idx_u) + Ka;
Me(idx_u, idx_u) = Me(idx_u, idx_u) + Ma;

% --- Torsion (rx about local x) ---
Kt = G*J/L * [ 1 -1; -1  1 ];
% approximate consistent rotary inertia about x:
Jrho = rho*J; % (mass moment per unit length) * L later
Mt = Jrho*L/6 * [ 2  1;  1  2 ];
idx_rx = [4, 10];
Ke(idx_rx, idx_rx) = Ke(idx_rx, idx_rx) + Kt;
Me(idx_rx, idx_rx) = Me(idx_rx, idx_rx) + Mt;

% --- Bending in the x–y plane (deflection v, rotation rz), uses E*Iz ---
Kb = E*Iz/L^3 * [ 12,   6*L,  -12,   6*L;
                   6*L, 4*L^2, -6*L, 2*L^2;
                  -12,  -6*L,   12,  -6*L;
                   6*L, 2*L^2, -6*L, 4*L^2 ];
Mb = rho*A*L/420 * [ 156,   22*L,   54,  -13*L;
                     22*L,  4*L^2,  13*L, -3*L^2;
                      54,   13*L,  156,  -22*L;
                    -13*L, -3*L^2, -22*L,  4*L^2 ];
idx_v_rz = [2 6 8 12];
Ke(idx_v_rz, idx_v_rz) = Ke(idx_v_rz, idx_v_rz) + Kb;
Me(idx_v_rz, idx_v_rz) = Me(idx_v_rz, idx_v_rz) + Mb;

% --- Bending in the x–z plane (deflection w, rotation ry), uses E*Iy ---
Kb2 = E*Iy/L^3 * [ 12,   6*L,  -12,   6*L;
                    6*L, 4*L^2, -6*L, 2*L^2;
                   -12,  -6*L,   12,  -6*L;
                    6*L, 2*L^2, -6*L, 4*L^2 ];
Mb2 = rho*A*L/420 * [ 156,   22*L,   54,  -13*L;
                      22*L,  4*L^2,  13*L, -3*L^2;
                       54,   13*L,  156,  -22*L;
                     -13*L, -3*L^2, -22*L,  4*L^2 ];
idx_w_ry = [3 5 9 11];
Ke(idx_w_ry, idx_w_ry) = Ke(idx_w_ry, idx_w_ry) + Kb2;
Me(idx_w_ry, idx_w_ry) = Me(idx_w_ry, idx_w_ry) + Mb2;

% Symmetrize (guards against minor FP asymmetry)
Ke = (Ke + Ke.')/2;
Me = (Me + Me.')/2;

end
