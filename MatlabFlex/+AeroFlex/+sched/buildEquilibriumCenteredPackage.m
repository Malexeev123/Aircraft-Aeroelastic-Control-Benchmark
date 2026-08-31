function sched = buildEquilibriumCenteredPackage(primitive,vertices,weights,xEq,uEq)
%BUILDEQUILIBRIUMCENTEREDPACKAGE Translate a static scheduled ROM about trim.
%   The returned absolute-state package is exactly equivalent to blending
%   the vertex vector fields at equal perturbations about their validated
%   equilibria. Vertex states use the convention x_i = xEq_i + delta_x.

arguments
    primitive (1,1) struct
    vertices (1,:) struct
    weights (:,1) double
    xEq (:,1) double
    uEq (:,1) double
end

weights = weights(:);
if numel(vertices) ~= numel(weights)
    error('AeroFlex:sched:CenteredVertexCount', ...
        'Vertex count (%d) differs from weight count (%d).', ...
        numel(vertices),numel(weights));
end
if abs(sum(weights)-1) > 1e-13 || any(weights < -1e-13)
    error('AeroFlex:sched:CenteredWeights', ...
        'Centered-package weights violate the interpolation policy.');
end

required = {'sched','xEq','uEq'};
for k = 1:numel(vertices)
    for j = 1:numel(required)
        if ~isfield(vertices(k),required{j})
            error('AeroFlex:sched:CenteredVertexField', ...
                'Vertex %d is missing %s.',k,required{j});
        end
    end
end

oneHot = find(abs(weights-1) <= 1e-14,1);
if ~isempty(oneHot) && nnz(abs(weights) > 1e-14) == 1
    sched = vertices(oneHot).sched;
    sched.x_eq = xEq;
    sched.u_eq = uEq;
    sched.equilibriumCentered = struct('enabled',false, ...
        'endpointExact',true,'vertexIds',vertexIds(vertices), ...
        'weights',weights.','rEq',evaluateVertex(vertices(oneHot)));
    return
end

sched = primitive;
nx = numel(xEq);
nu = numel(uEq);
rEq = zeros(nx,1);
AEq = zeros(nx);
BEq = zeros(nx,nu);
for k = 1:numel(vertices)
    vertex = vertices(k);
    pc = vertex.sched.parConst;
    pc.u_ctrl = vertex.uEq(:);
    [Nq,Nu] = AeroFlex.sim.nonlinearJacobian( ...
        vertex.xEq(:),vertex.sched.idx,pc);
    B = Nu(:,nx+2:nx+1+nu);
    rEq = rEq + weights(k)*evaluateVertex(vertex);
    AEq = AEq + weights(k)*(vertex.sched.L+Nq);
    BEq = BEq + weights(k)*B;
end

pcH = sched.parConst;
pcH.forces_0 = zeros(size(pcH.forces_0));
pcH.N_Thrust = zeros(size(pcH.N_Thrust));
pcH.u_ctrl = [];
JH = AeroFlex.sim.nonlinearJacobian(xEq,sched.idx,pcH);
Hxe = AeroFlex.sim.nonlinear_terms(xEq,pcH,sched.idx);
pcZero = sched.parConst;
pcZero.u_ctrl = [];
if isfield(pcZero,'affineOffset'), pcZero.affineOffset = []; end
nativeConstant = AeroFlex.sim.nonlinear_terms(zeros(nx,1),pcZero,sched.idx);

% Exact expansion of r + A*dx + H(dx,dx) + B*du in absolute coordinates.
sched.L = AEq-JH;
sched.parConst.affineOffset = ...
    rEq-AEq*xEq+Hxe-BEq*uEq-nativeConstant;
sched.x_eq = xEq;
sched.u_eq = uEq;
sched.parConst.u_ctrl = uEq;
sched.equilibriumCentered = struct('enabled',true, ...
    'endpointExact',false,'vertexIds',vertexIds(vertices), ...
    'weights',weights.','rEq',rEq,'AEq',AEq,'BEq',BEq, ...
    'form','absolute expansion of equilibrium-centered quadratic field', ...
    'recoveryVertices',buildRecoveryVertices(vertices,weights), ...
    'recoveryAnchorWrench',buildRecoveryAnchor(vertices,weights));
end

function out = buildRecoveryVertices(vertices,weights)
out = repmat(struct('weight',0,'xEq',[],'uEq',[],'Lq1',[], ...
    'parConst',struct(),'idx',struct(),'Pz',[],'Pr',[],'phi1_sA',[]), ...
    numel(vertices),1);
for k = 1:numel(vertices)
    v = vertices(k);
    out(k) = struct('weight',weights(k),'xEq',v.xEq(:), ...
        'uEq',v.uEq(:),'Lq1',v.sched.L(v.sched.idx.q1,:), ...
        'parConst',v.sched.parConst,'idx',v.sched.idx, ...
        'Pz',v.sched.beam.Pz,'Pr',v.sched.beam.Pr, ...
        'phi1_sA',v.sched.beam.red.phi1_sA);
end
end

function wrench = buildRecoveryAnchor(vertices,weights)
wrench = zeros(6,1);
for k = 1:numel(vertices)
    sched = vertices(k).sched;
    assert(isfield(sched,'p5') && isfield(sched.p5,'r2') && ...
        isfield(sched.p5.r2,'anchor') && ...
        isfield(sched.p5.r2.anchor,'wrench'), ...
        'AeroFlex:sched:CenteredRecoveryAnchor', ...
        'Vertex %d has no certified T2 recovery anchor.',k);
    sourceWrench = sched.p5.r2.anchor.wrench(:);
    assert(numel(sourceWrench)==6 && all(isfinite(sourceWrench)), ...
        'AeroFlex:sched:CenteredRecoveryAnchor', ...
        'Vertex %d has an invalid T2 recovery anchor.',k);
    wrench = wrench+weights(k)*sourceWrench;
end
end

function f = evaluateVertex(vertex)
pc = vertex.sched.parConst;
pc.u_ctrl = vertex.uEq(:);
f = vertex.sched.L*vertex.xEq(:) + ...
    AeroFlex.sim.nonlinear_terms(vertex.xEq(:),pc,vertex.sched.idx);
end

function ids = vertexIds(vertices)
ids = nan(1,numel(vertices));
for k = 1:numel(vertices)
    if isfield(vertices(k),'id'), ids(k) = vertices(k).id; end
end
end
