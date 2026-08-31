function wrench = recoverCenteredRootWrench(sched,x,u,gust)
%RECOVERCENTEREDROOTWRENCH Blend physical vertex reactions about equilibrium.

if nargin < 4 || isempty(gust)
    gust = [];
else
    gust = double(gust(:));
    assert(all(isfinite(gust)), ...
        'AeroFlex:sched:CenteredRecoveryGust', ...
        'The centered-recovery gust vector must be finite.');
end

if ~isfield(sched,'equilibriumCentered') || ...
        ~isfield(sched.equilibriumCentered,'enabled') || ...
        ~sched.equilibriumCentered.enabled
    error('AeroFlex:sched:NotCenteredPackage', ...
        'The scheduled package is not equilibrium centered.');
end
vertices = sched.equilibriumCentered.recoveryVertices;
anchored = isfield(sched.equilibriumCentered,'recoveryAnchorWrench');
if anchored
    anchor = sched.equilibriumCentered.recoveryAnchorWrench(:);
else
    anchor = zeros(6,1);
end
if numel(anchor)~=6 || any(~isfinite(anchor))
    error('AeroFlex:sched:CenteredRecoveryAnchorInvalid', ...
        'The absolute T2 reaction anchor must have six finite components.');
end
if isfield(sched,'physicalRecoveryReferenceState') && ...
        isfield(sched,'physicalRecoveryReferenceControl')
    referenceState = sched.physicalRecoveryReferenceState(:);
    referenceControl = sched.physicalRecoveryReferenceControl(:);
else
    referenceState = sched.x_eq(:);
    referenceControl = sched.u_eq(:);
end
assert(numel(referenceState)==numel(x) && ...
    numel(referenceControl)==numel(u), ...
    'AeroFlex:sched:CenteredRecoveryReferenceSize', ...
    'Physical recovery references must match the state and input sizes.');
dx = x(:)-referenceState;
du = u(:)-referenceControl;
wrench = anchor;
for k = 1:numel(vertices)
    v = vertices(k);
    xVertex = v.xEq+dx;
    uVertex = v.uEq+du;
    currentWrench = vertexReaction(v,xVertex,uVertex,gust);
    if anchored
        referenceWrench = vertexReaction(v,v.xEq,v.uEq,[]);
        increment = currentWrench-referenceWrench;
    else
        increment = currentWrench;
    end
    wrench = wrench+v.weight*increment;
end
end

function wrench = vertexReaction(vertex,state,control,gust)
pc = vertex.parConst;
pc.u_ctrl = control;
pc.gust = gust;
nonlinear = AeroFlex.sim.nonlinear_terms(state,pc,vertex.idx);
q1raw = vertex.Lq1*state+nonlinear(vertex.idx.q1);
wrench = AeroFlex.beam.recoverRootWrench( ...
    vertex.phi1_sA,vertex.Pr,vertex.Pr*q1raw);
end
