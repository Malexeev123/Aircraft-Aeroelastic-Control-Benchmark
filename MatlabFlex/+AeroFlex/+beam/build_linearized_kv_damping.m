function Dlin = build_linearized_kv_damping(fem, keepDofs, chain, elemLen, Eall, opts)
% BUILD_LINEARIZED_KV_DAMPING  Artola §2.1.4 linearised Kelvin–Voigt damping
%
% Inputs
%   fem        :  FEM struct (coordinates, connectivities, num_node, num_elem, etc.)
%   keepDofs   : active DOFs (vector of indices into full 6*N)
%   chain      : node ordering along the wing/beam centerline
%   elemLen    : length per sub-element (size nSub = length(chain)-1)
%   Eall       : stacked 6x6 intrinsic E blocks per sub-element (rows follow your rowIdx layout)
%   opts       : struct with any of
%                .Cd_per_elem{e}  (6x6 per element)   OR
%                .eta_f, .eta_m   (scalars)           OR
%                .use_modal_zeta  (logical)           -- handled elsewhere
%
% Output
%   Dlin : struct containing
%          .D0_node   (ndof x ndof) node-space
%          .D1_sub    (6*nSub x ndof) multiplies S1*x1
%          .D2_sub    (6*nSub x ndof) multiplies S2*x1
%          .S1, .S2   (6*nSub x ndof) first/second derivative operators
%          .gather    function handle to gather sub->node (trap weights)
%          .nSub
%

ndof  = length(keepDofs);
nSub  = length(chain)-1;
D0_node = sparse(ndof, ndof);
D1_sub  = sparse(6*nSub, ndof);
D2_sub  = sparse(6*nSub, ndof);

% ---- Build spatial derivative operators S1 (sub->node) and S2 (sub->node)
[S1, S2] = AeroFlex.beam.build_spatial_derivative_operators(fem, keepDofs, chain, elemLen, Eall);

% ---- Build Cd per sub-element (local 6x6)
Cd_sub = cell(nSub,1);
for i = 1:nSub
    % Map sub-element to owning finite element 'e'
    e = AeroFlex.beam.sub_to_elem(fem, chain, i);
    if isfield(opts,'Cd_per_elem') && ~isempty(opts.Cd_per_elem)
        Cd = opts.Cd_per_elem{e};                 % full 6x6
    else
        % Material ETA fallback (uniform diagonal per block)
        eta_f = getfield(opts,'eta_f', 0); 
        eta_m = getfield(opts,'eta_m', 0); 
        Cd = blkdiag(eta_f*eye(3), eta_m*eye(3));
    end
    Cd_sub{i} = 0.5*(Cd+Cd.'); % enforce symmetry
end

% ---- Assemble D0_node and D1_sub/D2_sub (Artola 2.25)
rowBase = @(i) (i-1)*6 + (1:6);   % sub row block

for i = 1:nSub
    r = rowBase(i);
    E6 = Eall(r,:);        % 6x6 intrinsic block at sub i
    Cd = Cd_sub{i};        % 6x6

    % D2 is simply Cd, acts on x1'' -> we store as sub-operator
    D2_sub(r, :) = place_sub_to_node(D2_sub(r,:), i, chain, keepDofs, Cd);

    % Spatial derivatives of Cd and (Cd E^T) are neglected by default unless varying materially.
    % If you want, you can add centered differences across i-1,i,i+1 here.
    Cd_prime = zeros(6);       % default: piecewise-constant
    CdEt     = Cd * E6.';
    CdEt_prime = zeros(6);     % default

    % D1 = Cd' - Cd E^T + E Cd
    D1_local = Cd_prime - Cd*E6.' + E6*Cd;

    % D0 = -(Cd E^T)' - E Cd E^T  (neglecting (Cd E^T)' spatial derivative unless provided)
    D0_local = -(CdEt_prime) - E6*Cd*E6.';

    % Distribute D0_local to nodes with trapezoidal weights
    D0_node = D0_node + trapz_scatter_node(D0_local, i, chain, keepDofs);

    % Store D1_local acting on S1*x1 (so D1_sub * S1 * x1)
    D1_sub(r, :) = place_sub_to_node(D1_sub(r,:), i, chain, keepDofs, D1_local);
end

% Gatherer from sub -> node (trapezoidal)
gather = @(Zsub) gather_sub_to_node(Zsub, chain, keepDofs);

Dlin = struct('D0_node', D0_node, ...
              'D1_sub',  D1_sub, ...
              'D2_sub',  D2_sub, ...
              'S1',      S1, ...
              'S2',      S2, ...
              'gather',  gather, ...
              'nSub',    nSub);
end

% ---- helpers --------------------------------------------------------------

function Mplaced = place_sub_to_node(Mplaced, iSub, chain, keepDofs, B6)
% places a 6x6 block B6 into the sub-row r and node columns (i, i+1)
% so that Mplaced * x approximates B6 * (local quantity at sub i)
nA = chain(iSub);
nB = chain(iSub+1);
colsA = localIdx(keepDofs, (nA-1)*6 + (1:6));
colsB = localIdx(keepDofs, (nB-1)*6 + (1:6));
if ~isempty(colsA), Mplaced(:, colsA) = Mplaced(:, colsA) + 0.5*B6; end
if ~isempty(colsB), Mplaced(:, colsB) = Mplaced(:, colsB) + 0.5*B6; end
end

function D0_acc = trapz_scatter_node(D0_local, iSub, chain, keepDofs)
% scatter a 6x6 local node-multiplier onto global node x node with trap weights
ndof = length(keepDofs);
D0_acc = sparse(ndof, ndof);
nA = chain(iSub);
nB = chain(iSub+1);
ia = localIdx(keepDofs, (nA-1)*6 + (1:6));
ib = localIdx(keepDofs, (nB-1)*6 + (1:6));
wA = 0.5; wB = 0.5;
if ~isempty(ia), D0_acc(ia, ia) = D0_acc(ia, ia) + wA*D0_local; end
if ~isempty(ib), D0_acc(ib, ib) = D0_acc(ib, ib) + wB*D0_local; end
end

function I = localIdx(keepDofs, dofRange)
I = [];
for d = dofRange(:).'
    k = find(keepDofs==d, 1, 'first');
    if ~isempty(k), I(end+1) = k; end %#ok<AGROW>
end
end
function Znode = gather_sub_to_node(Zsub, chain, keepDofs)
% trapezoidal gather from sub-space (6*nSub) to node-space (ndof)
ndof = length(keepDofs);
nSub = size(Zsub,1)/6;
Znode = sparse(ndof, ndof);
rowBase = @(i) (i-1)*6 + (1:6);
for i=1:nSub
    r  = rowBase(i);
    nA = chain(i);
    nB = chain(i+1);
    ia = localIdx(keepDofs, (nA-1)*6 + (1:6));
    ib = localIdx(keepDofs, (nB-1)*6 + (1:6));
    if ~isempty(ia), Znode(ia, :) = Znode(ia, :) + 0.5*Zsub(r, :); end
    if ~isempty(ib), Znode(ib, :) = Znode(ib, :) + 0.5*Zsub(r, :); end
end
end

% function I = localIdx(keepDofs, dofRange)
% I = [];
% for d = dofRange(:).'
%     k = find(keepDofs==d, 1, 'first');
%     if ~isempty(k), I(end+1) = k; end
% end
% end