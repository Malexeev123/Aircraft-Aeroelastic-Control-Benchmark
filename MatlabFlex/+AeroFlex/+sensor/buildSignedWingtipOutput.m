function contract = buildSignedWingtipOutput(source)
%BUILDSIGNEDWINGTIPOUTPUT Build package-owned physical wingtip outputs.
%   The returned contract recovers the signed positive-span tip, signed
%   negative-span tip, and their symmetric-longitudinal mean in the A-frame.
%   Outputs are available as absolute positions and as deviations from the
%   supplied trim state. The recovery maps are owned by the active scheduled
%   package when available.

arguments
    source struct
end

required = ["beam","base","trim","p5"];
for field = required
    assert(isfield(source,field), ...
        "AeroFlex:SignedWingtipSource", ...
        "The signed-wingtip source is missing field %s.",field);
end
assert(isfield(source.trim,"states") && isfield(source.p5,"idx") && ...
    isfield(source.p5,"L"),"AeroFlex:SignedWingtipPackage", ...
    "The signed-wingtip source requires trim states, idx, and L.");

idx = source.p5.idx;
trimState = source.trim.states(:);
nx = numel(trimState);
assert(all(isfinite(trimState)) && max([idx.q2(:);idx.qxi(:)]) <= nx, ...
    "AeroFlex:SignedWingtipState", ...
    "The trim state or state ordering is invalid.");

[q2ToPhysical,qxiToPhysical,mapOwner] = localRecoveryMaps(source,idx);
[activeNodes,coordinates,positiveBlock,negativeBlock] = ...
    localPhysicalTips(source.beam,q2ToPhysical,qxiToPhysical);

measureAbsolute = @(state) localMeasure(state,idx,q2ToPhysical, ...
    qxiToPhysical,coordinates,positiveBlock,negativeBlock);
trimAbsolute = measureAbsolute(trimState);
measureTrimRelative = @(state) measureAbsolute(state)-trimAbsolute;

step = 1e-7*max(1,abs(trimState));
trimJacobian = localCentralJacobian(measureAbsolute,trimState,step);
halfStepJacobian = localCentralJacobian( ...
    measureAbsolute,trimState,0.5*step);
stepSensitivity = norm(trimJacobian-halfStepJacobian,"fro")/ ...
    max(norm(halfStepJacobian,"fro"),eps);

contract = struct();
contract.schemaVersion = "signed-wingtip-output-v1";
contract.ownerPolicy = "package_owned_signed_mirrored_tip";
contract.mapOwner = mapOwner;
contract.measureAbsolute = measureAbsolute;
contract.measureTrimRelative = measureTrimRelative;
contract.trimAbsoluteMeters = trimAbsolute;
contract.trimJacobian = trimJacobian;
contract.symmetricGradient = trimJacobian(3,:);
contract.gradientStep = step;
contract.gradientStepSensitivity = stepSensitivity;
contract.outputOrder = ["positive_tip_z","negative_tip_z", ...
    "symmetric_mean_z"];
contract.units = "m";
contract.positiveTip = struct("femNode",activeNodes(positiveBlock), ...
    "activeBlock",positiveBlock,"yMeters",coordinates(positiveBlock,2));
contract.negativeTip = struct("femNode",activeNodes(negativeBlock), ...
    "activeBlock",negativeBlock,"yMeters",coordinates(negativeBlock,2));
end

function [q2Map,qxiMap,owner] = localRecoveryMaps(source,idx)
chart = struct();
if isfield(source,"scheduleStateCoordinate") && ...
        isstruct(source.scheduleStateCoordinate)
    chart = source.scheduleStateCoordinate;
elseif isfield(source.p5,"scheduleStateCoordinate") && ...
        isstruct(source.p5.scheduleStateCoordinate)
    chart = source.p5.scheduleStateCoordinate;
end
if isfield(chart,"q2ToPhysical") && isfield(chart,"qxiToPhysical")
    q2Map = double(chart.q2ToPhysical);
    qxiMap = double(chart.qxiToPhysical);
    owner = "scheduleStateCoordinate";
else
    assert(localHasMember(source.p5,"beam") && ...
        localHasMember(source.p5,"base") && ...
        localHasMember(source.p5.beam,"red") && ...
        localHasMember(source.p5.beam.red,"ModeVars_discrete") && ...
        localHasMember(source.p5.beam.red.ModeVars_discrete,"phi0_local") && ...
        localHasMember(source.p5.base,"phi_xi_modes"), ...
        "AeroFlex:SignedWingtipRecovery", ...
        "The active package does not expose physical recovery maps.");
    omega = -double(source.p5.L(idx.q2,idx.q1));
    q2Map = -double( ...
        source.p5.beam.red.ModeVars_discrete.phi0_local)/omega;
    qxiMap = double(source.p5.base.phi_xi_modes);
    owner = "exact_package_recovery";
end
assert(size(q2Map,2) == numel(idx.q2) && ...
    size(qxiMap,2) == numel(idx.qxi) && ...
    all(isfinite(q2Map),"all") && all(isfinite(qxiMap),"all"), ...
    "AeroFlex:SignedWingtipRecovery", ...
    "The package-owned physical recovery maps are invalid.");
end

function [activeNodes,coordinates,positiveBlock,negativeBlock] = ...
        localPhysicalTips(beam,q2Map,qxiMap)
assert(localHasMember(beam,"fem") && ...
    localHasMember(beam.fem,"keepDofs") && ...
    localHasMember(beam.fem,"coordinates"), ...
    "AeroFlex:SignedWingtipFEM", ...
    "The retained FEM coordinates are required for physical tip ownership.");
activeNodes = unique(ceil(double(beam.fem.keepDofs(:))/6),"stable");
coordinates = double(beam.fem.coordinates(activeNodes,:));
assert(size(coordinates,2) >= 3 && all(isfinite(coordinates),"all") && ...
    size(q2Map,1) == 6*numel(activeNodes) && ...
    size(qxiMap,1) == 4*numel(activeNodes), ...
    "AeroFlex:SignedWingtipFEM", ...
    "The recovery rows do not match the retained physical FEM nodes.");
[~,positiveBlock] = max(coordinates(:,2));
[~,negativeBlock] = min(coordinates(:,2));
assert(positiveBlock ~= negativeBlock && ...
    coordinates(positiveBlock,2) > 0 && coordinates(negativeBlock,2) < 0, ...
    "AeroFlex:SignedWingtipFEM", ...
    "Distinct positive- and negative-span physical tips were not found.");
end

function present = localHasMember(value,name)
present = (isstruct(value) && isfield(value,name)) || ...
    (isobject(value) && isprop(value,name));
end

function output = localMeasure(state,idx,q2Map,qxiMap,coordinates, ...
        positiveBlock,negativeBlock)
state = state(:);
assert(numel(state) >= max([idx.q2(:);idx.qxi(:)]) && ...
    all(isfinite(state)),"AeroFlex:SignedWingtipState", ...
    "The physical-output state must be finite and correctly ordered.");
localConfiguration = q2Map*state(idx.q2);
quaternion = reshape(qxiMap*state(idx.qxi),4,[]).';
quaternionNorm = vecnorm(quaternion,2,2);
assert(all(isfinite(quaternionNorm)) && all(quaternionNorm > 0), ...
    "AeroFlex:SignedWingtipQuaternion", ...
    "The recovered physical quaternion contains a zero or nonfinite norm.");
quaternion = quaternion./quaternionNorm;
blocks = [positiveBlock,negativeBlock];
tipZ = zeros(2,1);
for index = 1:2
    block = blocks(index);
    translationRows = (block-1)*6+(1:3);
    translation = AeroFlex.core.T_phi_quat(quaternion(block,:).')* ...
        localConfiguration(translationRows);
    tipZ(index) = coordinates(block,3)+translation(3);
end
output = [tipZ;mean(tipZ)];
end

function jacobian = localCentralJacobian(functionHandle,state,step)
value = functionHandle(state);
jacobian = zeros(numel(value),numel(state));
for column = 1:numel(state)
    plus = state;
    minus = state;
    plus(column) = plus(column)+step(column);
    minus(column) = minus(column)-step(column);
    jacobian(:,column) = (functionHandle(plus)-functionHandle(minus))/ ...
        (2*step(column));
end
assert(all(isfinite(jacobian),"all"), ...
    "AeroFlex:SignedWingtipJacobian", ...
    "The signed-wingtip central-difference Jacobian is nonfinite.");
end
