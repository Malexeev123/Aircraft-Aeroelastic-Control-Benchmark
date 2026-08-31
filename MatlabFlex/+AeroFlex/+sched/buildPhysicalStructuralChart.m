function [chart,cacheStats] = buildPhysicalStructuralChart( ...
        sources,weights,method,externalAerodynamicChart,options)
%BUILDPHYSICALSTRUCTURALCHART Build a source-only physical ROM chart.
%   The source packages are not modified. Reduced velocity coordinates use
%   q1_c=T1*q1_i, while their generalized-force dual uses
%   q2_c=inv(T1.')*q2_i. Quaternion equilibrium fields are interpolated in
%   the SO(3) tangent space and normalized after reconstruction.

arguments
    sources (1,:) struct
    weights (:,1) double
    method (1,1) string {mustBeMember(method, ...
        ["fixed_physical_svd","fixed_physical_biorthogonal", ...
        "cell_local_aligned","nearest_source"])}
    externalAerodynamicChart (1,1) logical = false
    options (1,1) struct = struct()
end

weights = weights(:);
nSource = numel(sources);
if numel(weights) ~= nSource || any(~isfinite(weights)) || ...
        abs(sum(weights)-1) > 1e-12
    error('AeroFlex:sched:PhysicalChartWeights', ...
        'Source count and finite partition-of-unity weights must agree.');
end

preparedRuntimeOwner = isfield(options,'preparedRuntimeOwner') && ...
    logical(options.preparedRuntimeOwner);
registry = struct();
if preparedRuntimeOwner
    requiredPrepared = {'preparedRuntimeProfile','preparedRuntimeChangeId', ...
        'preparedRuntimeRegistrySha256','registryPath'};
    missingPrepared = requiredPrepared(~isfield(options,requiredPrepared));
    assert(isempty(missingPrepared) && ...
        ismember(string(options.preparedRuntimeProfile), ...
            ["formal_case_b","formal_case_c"]) && ...
        string(options.preparedRuntimeChangeId)== ...
            "phase18c-v17a-casebc-prepared-scheduled-runtime-owner-v1", ...
        'AeroFlex:sched:PhysicalChartPreparedScope', ...
        ['Prepared physical-chart maps require the approved Case-B/C ', ...
         'profile, change identifier, and registry binding.']);
    registry = verifiedPreparedRegistry(string(options.registryPath), ...
        lower(string(options.preparedRuntimeRegistrySha256)));
end
cacheStats = struct('enabled',preparedRuntimeOwner,'hits',0, ...
    'misses',0,'h5Reads',0,'wallSeconds',0, ...
    'registrySha256',string.empty(1,0));
if preparedRuntimeOwner
    cacheStats.registrySha256 = ...
        lower(string(options.preparedRuntimeRegistrySha256));
end
maps = repmat(emptyMaps(),nSource,1);
for k = 1:nSource
    [maps(k),mapInfo] = extractMaps( ...
        sources(k),preparedRuntimeOwner,registry);
    cacheStats.hits = cacheStats.hits+double(mapInfo.hit);
    cacheStats.misses = cacheStats.misses+double(mapInfo.miss);
    cacheStats.h5Reads = cacheStats.h5Reads+mapInfo.h5Reads;
    cacheStats.wallSeconds = cacheStats.wallSeconds+mapInfo.wallSeconds;
end
oneHot = find(abs(weights-1) <= 1e-14,1);
isOneHot = ~isempty(oneHot) && nnz(abs(weights) > 1e-14) == 1;
if isOneHot
    referenceIndex = oneHot;
else
    [~,referenceIndex] = max(abs(weights));
end

H1 = selectBasis({maps.H1},weights,method,referenceIndex,isOneHot);
Hxi = selectBasis({maps.Hxi},weights,method,referenceIndex,isOneHot);
if externalAerodynamicChart
    % The atomic common-realization builder owns qGamma coordinates. Its
    % certified ambient E/R maps replace this legacy source-row chart.
    Hg = maps(referenceIndex).Hg;
else
    Hg = selectBasis({maps.Hg},weights,method,referenceIndex,isOneHot);
end

transforms = repmat(emptyTransform(),nSource,1);
H0 = zeros(size(maps(1).H0));
H2 = zeros(size(maps(1).H2));
maxCondition = 0;
maxBasisResidual = 0;
maxRoundTrip = 0;
maxPowerAbs = 0;
maxPowerRel = 0;
for k = 1:nSource
    if isOneHot
        T1 = eye(size(H1,2));
        Tx = eye(size(Hxi,2));
        Tg = eye(size(Hg,2));
        T2 = T1;
    else
        % The active intrinsic equations have an identity metric and reuse
        % Gamma2 in the q1 and q2 equations. An orthogonal physical-space
        % alignment is therefore required to retain the canonical ODE and
        % its nonlinear power cancellation without introducing a descriptor.
        if method == "fixed_physical_biorthogonal"
            T1 = H1\maps(k).H1;
            T2 = T1'\eye(size(T1));
        else
            T1 = AeroFlex.sched.procrustesLocalToRef( ...
                maps(k).H1,H1,'orthogonal');
            T2 = T1;
        end
        Tx = AeroFlex.sched.procrustesLocalToRef( ...
            maps(k).Hxi,Hxi,'orthogonal');
        if externalAerodynamicChart
            Tg = eye(size(maps(k).Hg,2));
        else
            Tg = Hg\maps(k).Hg;
        end
    end
    assertNonsingular(T1,'q1',k);
    assertNonsingular(Tx,'qxi',k);
    assertNonsingular(Tg,'qGamma',k);
    S1 = T1\eye(size(T1));
    S2 = S1;
    if method == "fixed_physical_biorthogonal" && ~isOneHot
        S2 = T1.';
    end
    Sx = Tx\eye(size(Tx));
    Sg = Tg\eye(size(Tg));
    T = blockTransform(sources(k).idx,T1,T2,Tx,Tg);
    Ti = blockTransform(sources(k).idx,S1,S2,Sx,Sg);
    transforms(k) = struct('q1',block(T1,S1), ...
        'q2',block(T2,S2),'qxi',block(Tx,Sx), ...
        'qGam',block(Tg,Sg),'T',T,'Tinv',Ti);
    H0 = H0+weights(k)*(maps(k).H0*S1);
    H2 = H2+weights(k)*(maps(k).H2*S2);
    maxCondition = max(maxCondition,max([cond(T1),cond(Tx),cond(Tg)]));
    basisResidual = norm(H1*T1-maps(k).H1,'fro')/ ...
        max(1,norm(maps(k).H1,'fro'));
    maxBasisResidual = max(maxBasisResidual,basisResidual);
    maxRoundTrip = max(maxRoundTrip,norm(Ti*T-eye(size(T)),'fro'));
    [powerAbs,powerRel] = virtualPowerErrors(T1,T2,k);
    maxPowerAbs = max(maxPowerAbs,powerAbs);
    maxPowerRel = max(maxPowerRel,powerRel);
end

physical = interpolateEquilibrium(maps,weights,H0,H1,H2,Hxi,referenceIndex);
chart = struct('schemaVersion',1,'method',char(method), ...
    'weights',weights.','referenceIndex',referenceIndex, ...
    'oneHot',isOneHot,'H0',H0,'H1',H1,'H2',H2, ...
    'Hxi',Hxi,'Hg',Hg,'transforms',transforms, ...
    'physicalEquilibrium',physical, ...
    'aerodynamicChartOwner',owner(externalAerodynamicChart), ...
    'metrics',struct('maxTransformCondition',maxCondition, ...
        'maxVelocityBasisResidual',maxBasisResidual, ...
        'maxStateRoundTripError',maxRoundTrip, ...
        'maxVirtualPowerAbsoluteError',maxPowerAbs, ...
        'maxVirtualPowerRelativeError',maxPowerRel));
end

function value = owner(externalAerodynamicChart)
if externalAerodynamicChart
    value = 'atomic_common_realization';
else
    value = 'legacy_physical_structural_chart';
end
end

function [maps,info] = extractMaps(source,usePreparedCache,registry)
startTime = tic;
info = struct('hit',false,'miss',false,'h5Reads',0,'wallSeconds',0);
required = {'L','idx','beam','base','x_eq'};
for k = 1:numel(required)
    if ~isfield(source,required{k})
        error('AeroFlex:sched:PhysicalChartSourceField', ...
            'Source package is missing %s.',required{k});
    end
end
disc = source.beam.red.ModeVars_discrete;
H0 = double(disc.phi0_local);
H1 = double(disc.phi1_local);
Hxi = double(source.base.phi_xi_modes);
xiBar = normalizeQuaternionRows(double(source.base.xi_bar));
omega = -source.L(source.idx.q2,source.idx.q1);
equilibriumState = source.x_eq(:);
if isfield(source,'physicalChartEquilibrium') && ...
        isfield(source.physicalChartEquilibrium,'state')
    equilibriumState = source.physicalChartEquilibrium.state(:);
end
q0 = -(omega\equilibriumState(source.idx.q2));
% phi2_local is a complex intrinsic eigenvector map and is not a certified
% real physical-force chart. At equilibrium, the accepted real physical
% configuration map from q2 is q2 -> q0=-Omega\q2 -> phi0_local*q0.
H2 = -H0/omega;
quaternion = normalizeQuaternionField(Hxi*equilibriumState(source.idx.qxi));
if usePreparedCache
    [maps,cacheInfo] = preparedSourceMaps(source,registry,H0,H1,H2,Hxi, ...
        xiBar,omega,equilibriumState,q0,quaternion);
    info.hit = cacheInfo.hit;
    info.miss = ~cacheInfo.hit;
    info.h5Reads = cacheInfo.h5Reads;
    info.wallSeconds = toc(startTime);
    return
end
Hg = readPythonMatrix(source.p5.artifactPath,'/krylov_projector');
info.miss = true;
info.h5Reads = 1;
maps = assembleMaps(H0,H1,H2,Hxi,Hg,xiBar,equilibriumState, ...
    source.idx,q0,quaternion);
info.wallSeconds = toc(startTime);
end

function maps = assembleMaps( ...
        H0,H1,H2,Hxi,Hg,xiBar,equilibriumState,idx,q0,quaternion)
maps = struct('H0',H0,'H1',H1,'H2',H2,'Hxi',Hxi,'Hg',Hg, ...
    'q0Eq',q0,'q1Eq',equilibriumState(idx.q1), ...
    'q2Eq',equilibriumState(idx.q2), ...
    'qxiEq',equilibriumState(idx.qxi), ...
    'quaternionEq',quaternion,'xiBar',xiBar, ...
    'yEq',H0*q0,'f2Eq',H2*equilibriumState(idx.q2));
end

function H = selectBasis(bases,weights,method,referenceIndex,isOneHot)
if isOneHot || method == "nearest_source"
    H = bases{referenceIndex};
    return
end
rankValue = size(bases{1},2);
if method == "fixed_physical_svd"
    aggregate = [];
    for k = 1:numel(bases)
        Q = orth(bases{k});
        aggregate = [aggregate,sqrt(abs(weights(k)))*Q]; %#ok<AGROW>
    end
    [U,~,~] = svd(aggregate,'econ');
    U = deterministicSigns(U(:,1:rankValue));
    % Retain the deterministic reference normalization. The SVD selects a
    % physical subspace; it does not redefine the accepted modal scale.
    H = U*(U.'*bases{referenceIndex});
    return
end
Href = bases{referenceIndex};
aggregate = zeros(size(Href));
for k = 1:numel(bases)
    T = AeroFlex.sched.procrustesLocalToRef( ...
        bases{k},Href,'orthogonal');
    aggregate = aggregate+weights(k)*(bases{k}*T.');
end
if rank(aggregate) < rankValue
    error('AeroFlex:sched:PhysicalChartRank', ...
        'The aligned cell-local physical basis lost rank.');
end
H = aggregate;
end

function H = deterministicSigns(H)
for k = 1:size(H,2)
    [~,row] = max(abs(H(:,k)));
    if H(row,k) < 0, H(:,k) = -H(:,k); end
end
end

function physical = interpolateEquilibrium(maps,w,H0,H1,H2,Hxi,referenceIndex)
yTarget = zeros(size(maps(1).yEq));
f2Target = zeros(size(maps(1).f2Eq));
q1Target = zeros(size(maps(1).H1,1),1);
for k = 1:numel(maps)
    yTarget = yTarget+w(k)*maps(k).yEq;
    f2Target = f2Target+w(k)*maps(k).f2Eq;
    q1Target = q1Target+w(k)*(maps(k).H1*maps(k).q1Eq);
end
qTarget = interpolateQuaternions([maps.quaternionEq],w,referenceIndex);
xiBarTarget = interpolateQuaternions([maps.xiBar],w,referenceIndex);
qVector = reshape(qTarget.',[],1);
oneHot = find(abs(w-1) <= 1e-14,1);
if ~isempty(oneHot) && nnz(abs(w) > 1e-14) == 1
    q0 = maps(oneHot).q0Eq;
    q1 = maps(oneHot).q1Eq;
    q2 = maps(oneHot).q2Eq;
    qxi = maps(oneHot).qxiEq;
else
    q0 = H0\yTarget;
    q1 = H1\q1Target;
    q2 = H2\f2Target;
    qxi = Hxi\qVector;
end
rawQuaternion = Hxi*qxi;
reconstructedQuaternion = normalizeQuaternionField(rawQuaternion);
physical = struct('translationRotationTarget',yTarget, ...
    'internalForceTarget',f2Target,'quaternionTarget',qTarget, ...
    'baselineQuaternion',xiBarTarget, ...
    'q0',q0,'q1',q1,'q2',q2,'qxi',qxi, ...
    'translationRotationResidual',norm(H0*q0-yTarget)/max(1,norm(yTarget)), ...
    'internalForceResidual',norm(H2*q2-f2Target)/max(1,norm(f2Target)), ...
    'orientationCoordinateResidual',norm(rawQuaternion-qVector)/max(1,norm(qVector)), ...
    'orientationResidual',norm(reconstructedQuaternion-qTarget,'fro')/ ...
        max(1,norm(qTarget,'fro')), ...
    'sourceQuaternionNormalizationCorrection', ...
        norm(rawQuaternion-reshape(reconstructedQuaternion.',[],1))/ ...
        max(1,norm(rawQuaternion)), ...
    'quaternionNormError',max(abs(vecnorm(qTarget,2,2)-1)));
end

function Q = interpolateQuaternions(allQ,w,referenceIndex)
nNode = size(allQ,1);
nSource = size(allQ,2)/4;
Q = zeros(nNode,4);
for node = 1:nNode
    qRef = allQ(node,4*(referenceIndex-1)+(1:4));
    tangent = zeros(3,1);
    for k = 1:nSource
        q = allQ(node,4*(k-1)+(1:4));
        if dot(q,qRef) < 0, q = -q; end
        tangent = tangent+w(k)*quaternionLog(quaternionMultiply( ...
            quaternionConjugate(qRef),q));
    end
    Q(node,:) = quaternionMultiply(qRef,quaternionExp(tangent));
end
Q = Q./vecnorm(Q,2,2);
end

function value = quaternionLog(q)
q = q/norm(q);
if q(1) < 0, q = -q; end
s = norm(q(2:4));
if s <= 10*eps
    value = 2*q(2:4).';
else
    value = (2*atan2(s,q(1))/s)*q(2:4).';
end
end

function q = quaternionExp(value)
angle = norm(value);
if angle <= 10*eps
    q = [1,0.5*value(:).'];
else
    q = [cos(angle/2),sin(angle/2)*value(:).'/angle];
end
q = q/norm(q);
end

function q = quaternionMultiply(a,b)
q = [a(1)*b(1)-dot(a(2:4),b(2:4)), ...
    a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
q = q/norm(q);
end

function q = quaternionConjugate(q)
q = [q(1),-q(2:4)];
end

function Q = normalizeQuaternionField(value)
Q = reshape(value,4,[]).';
norms = vecnorm(Q,2,2);
if any(norms < 1e-12)
    error('AeroFlex:sched:PhysicalChartQuaternion', ...
        'A source equilibrium contains a zero quaternion.');
end
Q = Q./norms;
end

function Q = normalizeQuaternionRows(Q)
norms = vecnorm(Q,2,2);
if any(norms < 1e-12)
    error('AeroFlex:sched:PhysicalChartQuaternion', ...
        'A source baseline contains a zero quaternion.');
end
Q = Q./norms;
end

function [absoluteError,relativeError] = virtualPowerErrors(T1,T2,seed)
absoluteError = 0;
relativeError = 0;
vectors = {sin((1:size(T1,2)).'+0.1*seed), ...
    cos((1:size(T1,2)).'-0.2*seed)};
rngState = rng;
cleanup = onCleanup(@() rng(rngState));
rng(1800+seed,'twister');
vectors(end+1:end+4) = {randn(size(T1,2),1),randn(size(T1,2),1), ...
    randn(size(T1,2),1),randn(size(T1,2),1)};
for k = 1:2:numel(vectors)
    dq = vectors{k};
    force = vectors{k+1};
    sourcePower = force.'*dq;
    commonPower = (T2*force).'*(T1*dq);
    errorValue = abs(sourcePower-commonPower);
    absoluteError = max(absoluteError,errorValue);
    relativeError = max(relativeError,errorValue/max(1,abs(sourcePower)));
end
clear cleanup
end

function T = blockTransform(idx,T1,T2,Tx,Tg)
nx = max([idx.q1(:);idx.q2(:);idx.qxi(:);idx.qGam(:);idx.chi(:)]);
T = eye(nx);
T(idx.q1,idx.q1) = T1;
T(idx.q2,idx.q2) = T2;
T(idx.qxi,idx.qxi) = Tx;
T(idx.qGam,idx.qGam) = Tg;
end

function value = block(T,Tinv)
value = struct('T',T,'Tinv',Tinv);
end

function assertNonsingular(T,name,index)
if any(~isfinite(T(:))) || rcond(T) < 1e-12
    error('AeroFlex:sched:PhysicalChartCondition', ...
        'Source %d %s transform is singular or ill-conditioned.',index,name);
end
end

function value = readPythonMatrix(path,dataset)
value = h5read(path,dataset);
if isstruct(value) && isfield(value,'r')
    value = double(value.r)+1i*double(value.i);
else
    value = double(value);
end
if ismatrix(value) && ~isvector(value), value = value.'; end
value = real(value);
end

function registry = verifiedPreparedRegistry(path,expectedHash)
persistent cachedPath cachedHash cachedStamp cachedRegistry
path = canonicalPath(resolveRepositoryPath(path));
stamp = fileStamp(path);
if ~isempty(cachedRegistry) && path==cachedPath && ...
        expectedHash==cachedHash && isequal(stamp,cachedStamp)
    registry = cachedRegistry;
    return
end
observedHash = fileHash(path);
assert(observedHash==expectedHash, ...
    'AeroFlex:sched:PhysicalChartPreparedRegistryHash', ...
    'The prepared physical-chart registry changed (actual %s; expected %s).', ...
    observedHash,expectedHash);
value = jsondecode(fileread(path));
required = {'registryVersion','sourceCount','passedCount', ...
    'productionPromotionApproved','sourceOnly','sources'};
missing = required(~isfield(value,required));
assert(isempty(missing) && ...
    string(value.registryVersion)== ...
        "phase18c-v17a-dimensional-steady-force-production-registry-v1" && ...
    logical(value.productionPromotionApproved) && logical(value.sourceOnly) && ...
    double(value.sourceCount)==29 && double(value.passedCount)==29 && ...
    numel(value.sources)==29, ...
    'AeroFlex:sched:PhysicalChartPreparedRegistryStatus', ...
    'The prepared physical-chart registry is not the promoted 29-source V17A registry.');
sourceIds = string({value.sources.sourceId});
assert(numel(unique(sourceIds))==29 && ...
    all(string({value.sources.status})=="PASS") && ...
    all([value.sources.sourceOnly]), ...
    'AeroFlex:sched:PhysicalChartPreparedRegistrySource', ...
    'The prepared physical-chart registry has invalid source records.');
registry = struct('path',path,'sha256',observedHash, ...
    'repositoryRoot',localRepositoryRoot(),'sources',value.sources);
cachedPath = path;
cachedHash = expectedHash;
cachedStamp = stamp;
cachedRegistry = registry;
end

function [maps,info] = preparedSourceMaps( ...
        source,registry,H0,H1,H2,Hxi,xiBar,omega,equilibriumState,q0,quaternion)
persistent keys values stamps
if isempty(keys)
    keys = strings(0,1);
    values = cell(0,1);
    stamps = cell(0,1);
end
assert(isfield(source,'p5') && isstruct(source.p5) && ...
    all(isfield(source.p5,{'artifactPath','packageHash'})), ...
    'AeroFlex:sched:PhysicalChartPreparedPackage', ...
    'Prepared physical-chart maps require an artifact path and package hash.');
sourceId = sourceIdentity(source);
rowIndex = find(string({registry.sources.sourceId})==sourceId);
assert(isscalar(rowIndex), ...
    'AeroFlex:sched:PhysicalChartPreparedSource', ...
    'Source %s is absent or duplicated in the V17A registry.',sourceId);
row = registry.sources(rowIndex);
assert(isfield(row,'selectedContract') && ...
    all(isfield(row.selectedContract,{'path','sha256'})), ...
    'AeroFlex:sched:PhysicalChartPreparedContract', ...
    'Source %s has no selected-contract hash binding.',sourceId);
artifactPath = canonicalPath(resolveRepositoryPath( ...
    string(source.p5.artifactPath)));
registeredPath = canonicalPath(resolveRepositoryPath( ...
    string(row.selectedContract.path)));
assert(strcmpi(artifactPath,registeredPath), ...
    'AeroFlex:sched:PhysicalChartPreparedContractPath', ...
    'Source %s does not use its registered selected contract.',sourceId);
expectedHash = lower(string(row.selectedContract.sha256));
stamp = fileStamp(artifactPath);
residentDigest = numericDigest({H0,H1,H2,Hxi,xiBar,omega, ...
    equilibriumState,q0,quaternion,double(source.idx.q1(:)), ...
    double(source.idx.q2(:)),double(source.idx.qxi(:))});
key = strjoin([sourceId,string(source.p5.packageHash),residentDigest, ...
    artifactPath,expectedHash,"/krylov_projector"],"|");
cacheIndex = find(keys==key,1);
if ~isempty(cacheIndex) && isequal(stamps{cacheIndex},stamp)
    maps = values{cacheIndex};
    info = struct('hit',true,'h5Reads',0);
    return
end
observedHash = fileHash(artifactPath);
assert(observedHash==expectedHash, ...
    'AeroFlex:sched:PhysicalChartPreparedContractHash', ...
    'Selected contract hash mismatch for source %s.',sourceId);
Hg = readPythonMatrix(artifactPath,'/krylov_projector');
maps = assembleMaps(H0,H1,H2,Hxi,Hg,xiBar,equilibriumState, ...
    source.idx,q0,quaternion);
if numel(keys)>=32
    keys(1) = [];
    values(1) = [];
    stamps(1) = [];
end
keys(end+1,1) = key;
values{end+1,1} = maps;
stamps{end+1,1} = stamp;
info = struct('hit',false,'h5Reads',1);
end

function value = sourceIdentity(source)
if isfield(source,'sourceContractId') && ...
        strlength(string(source.sourceContractId))>0
    value = string(source.sourceContractId);
elseif isfield(source,'name') && strlength(string(source.name))>0
    value = string(source.name);
else
    error('AeroFlex:sched:PhysicalChartPreparedSourceIdentity', ...
        'Prepared physical-chart maps require a source identity.');
end
assert(isscalar(value), ...
    'AeroFlex:sched:PhysicalChartPreparedSourceIdentity', ...
    'Prepared physical-chart source identity must be scalar.');
end

function value = numericDigest(parts)
engine = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
for index = 1:numel(parts)
    part = double(parts{index});
    assert(all(isfinite(part),'all') && isreal(part), ...
        'AeroFlex:sched:PhysicalChartPreparedMapFinite', ...
        'Prepared physical-chart resident maps must be finite and real.');
    dimensions = uint64(size(part));
    engine.update(typecast(dimensions(:),'int8'));
    engine.update(typecast(typecast(part(:),'uint8'),'int8'));
end
value = lower(string(reshape(dec2hex( ...
    typecast(engine.digest(),'uint8'),2).',1,[])));
end

function value = fileHash(path)
engine = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
fileId = fopen(path,'rb');
assert(fileId>=0,'AeroFlex:sched:PhysicalChartPreparedHashOpen', ...
    'Cannot read %s.',path);
cleanup = onCleanup(@() fclose(fileId));
while ~feof(fileId)
    bytes = fread(fileId,1024*1024,'*uint8');
    if isempty(bytes), break, end
    engine.update(typecast(bytes(:),'int8'));
end
value = lower(string(reshape(dec2hex( ...
    typecast(engine.digest(),'uint8'),2).',1,[])));
clear cleanup
end

function value = fileStamp(path)
listing = dir(path);
assert(isscalar(listing) && ~listing.isdir, ...
    'AeroFlex:sched:PhysicalChartPreparedFile', ...
    'Prepared physical-chart file is unavailable: %s',path);
value = struct('bytes',double(listing.bytes),'datenum',double(listing.datenum));
end

function path = resolveRepositoryPath(path)
path = string(path);
if ~isfile(path)
    path = fullfile(localRepositoryRoot(),path);
end
assert(isfile(path), ...
    'AeroFlex:sched:PhysicalChartPreparedFile', ...
    'Prepared physical-chart file is unavailable: %s',path);
end

function path = canonicalPath(path)
file = javaObject('java.io.File',char(path));
path = string(file.getCanonicalPath());
end

function root = localRepositoryRoot()
root = fileparts(fileparts(fileparts(fileparts( ...
    mfilename('fullpath')))));
end

function value = emptyMaps()
value = struct('H0',[],'H1',[],'H2',[],'Hxi',[],'Hg',[], ...
    'q0Eq',[],'q1Eq',[],'q2Eq',[],'qxiEq',[], ...
    'quaternionEq',[],'xiBar',[],'yEq',[],'f2Eq',[]);
end

function value = emptyTransform()
value = struct('q1',block([],[]),'q2',block([],[]), ...
    'qxi',block([],[]),'qGam',block([],[]),'T',[],'Tinv',[]);
end
