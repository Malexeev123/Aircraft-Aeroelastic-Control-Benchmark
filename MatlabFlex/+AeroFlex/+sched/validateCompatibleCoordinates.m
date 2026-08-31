function val = validateCompatibleCoordinates(points,sourcePoints)
%VALIDATECOMPATIBLECOORDINATES Check transformed structural/recovery maps.

if nargin < 2 || isempty(sourcePoints)
    sourcePoints = points;
end
if numel(points) ~= numel(sourcePoints)
    error('validateCompatibleCoordinates:PointCount', ...
        'Transformed and source point counts must agree.');
end

val = struct('dPz',0,'dPr',0,'dPhi',0,'projectorError',0, ...
    'orthogonalityError',0,'recoveryEndpointError',0, ...
    'virtualPowerError',0,'compatible',true);

phiRef = points(1).beam.red.phi1_sA;
for i = 1:numel(points)
    target = points(i);
    source = sourcePoints(i);
    tr = target.compat.local_to_ref;
    T = tr.q1.T;
    S = tr.q1.Tinv;

    expectedPz = T*source.beam.Pz*S;
    expectedPr = T*source.beam.Pr*S;
    expectedPhi = source.beam.red.phi1_sA*S;
    val.dPz = max(val.dPz,norm(target.beam.Pz-expectedPz,'fro'));
    val.dPr = max(val.dPr,norm(target.beam.Pr-expectedPr,'fro'));
    val.dPhi = max(val.dPhi,norm(target.beam.red.phi1_sA-phiRef,'fro'));
    val.projectorError = max(val.projectorError, ...
        norm(target.beam.Pz^2-target.beam.Pz,'fro'));
    val.orthogonalityError = max(val.orthogonalityError, ...
        norm(T.'*T-eye(size(T)),'fro'));

    wrenchInput = [0.7;0;-1.1;0;0.3;0];
    reactionSource = source.beam.Pr*source.beam.red.phi1_sA.'*wrenchInput;
    reactionCommon = T*reactionSource;
    [wrenchSource,~] = AeroFlex.beam.recoverRootWrench( ...
        source.beam.red.phi1_sA,source.beam.Pr,reactionSource);
    [wrenchCommon,~] = AeroFlex.beam.recoverRootWrench( ...
        target.beam.red.phi1_sA,target.beam.Pr,reactionCommon);
    val.recoveryEndpointError = max(val.recoveryEndpointError, ...
        norm(wrenchCommon-wrenchSource)/max(1,norm(wrenchSource)));

    rateSource = sin((1:size(T,1)).'+0.1*i);
    rateCommon = T*rateSource;
    velocitySource = source.beam.red.phi1_sA*rateSource;
    velocityCommon = expectedPhi*rateCommon;
    powerScale = max(1,abs(wrenchSource.'*velocitySource));
    val.virtualPowerError = max(val.virtualPowerError, ...
        max(abs([wrenchSource.'*velocitySource- ...
        reactionSource.'*rateSource, ...
        wrenchCommon.'*velocityCommon-reactionCommon.'*rateCommon])) / ...
        powerScale);
end

tol = 1e-10;
val.compatible = max([val.dPz,val.dPr,val.projectorError, ...
    val.orthogonalityError,val.recoveryEndpointError, ...
    val.virtualPowerError]) <= tol;
end
