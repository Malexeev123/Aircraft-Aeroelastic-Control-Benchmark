function [fem, keepDofs, chain, elemLen, Eall] = make_uniform_beam(L,Nn,bc,E,A,Iy,Iz,GJ,rho)
% Straight beam along +y, 3-node Timoshenko-style elements (use your layout).
% bc: 'cantilever' | 'hinged'

% nodes
y = linspace(0,L,Nn).';
coords = [zeros(Nn,1), y, zeros(Nn,1)];

% simple chain
chain = (1:Nn)';

% elements: (n1, n3, n2) convention
Ne = Nn-1;
connect = zeros(Ne,3);
for e=1:Ne
    connect(e,:) = [e, e+1, e+1]; % stub: adjust to your 3-node convention if needed
end

% boundary conditions
bcs = zeros(Nn,1);
switch lower(bc)
    case 'cantilever'
        bcs(1) = 1;           % clamp node 1
    case 'hinged'
        bcs([1,end]) = -1;    % simply supported ends (you may need to adapt how you treat this)
    otherwise
        error('bc must be cantilever|hinged');
end

% pack fem
fem = struct();
fem.coordinates      = coords;
fem.connectivities   = connect;
fem.num_node         = Nn;
fem.num_elem         = Ne;
fem.boundary_conditions = bcs;
fem.structural_twist = zeros(Ne,3);  % no initial twist/curvature
fem.lumped_mass_nodes = (1:Nn);      % keep all nodes active

% active DOFs = all nodes except clamps for cantilever
if strcmpi(bc,'cantilever')
    root = 1;
    doNotKeep = ( (root-1)*6 + (1:6) );
    keep = setdiff(1:6*Nn, doNotKeep);
else
    keep = 1:6*Nn;
end
keepDofs = keep(:).';

% per-sub-element lengths
elemLen = diff(y);

% intrinsic E per sub (straight beam => E = [e1; 0]*L1 -> here e1=[1,0,0], no curvature)
Eall = zeros(6*Ne,6);
e1 = [1;0;0];
for i=1:Ne
    kappa0 = zeros(3,1);
    E6 = L1_operator([e1; kappa0]);
    r = (i-1)*6+(1:6);
    Eall(r,:) = E6;
end
end
function M6 = L1_operator(x6)
v = x6(1:3); w = x6(4:6);
M6 = zeros(6);
M6(1:3,1:3)= skew(w);
M6(4:6,1:3)= skew(v);
M6(4:6,4:6)= skew(w);
end

function S = skew(u)
S= [ 0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0 ];
end
