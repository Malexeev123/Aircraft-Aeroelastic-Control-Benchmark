function candidate = buildAtomicCommonRealization( ...
        sources,weights,chart,registryPath,ambientPath)
%BUILDATOMICCOMMONREALIZATION Build the locked source-only fixed10 tuple.
%   The full premodal source systems are embedded in the certified 1017-state
%   Option-A ambient space. Source sample times are first mapped through the
%   bilinear continuous-time representation and resampled at the query time
%   step. One blended biorthogonal V/W pair then produces A/B/C/D together.

arguments
    sources (1,:) struct
    weights (:,1) double
    chart (1,1) struct
    registryPath (1,1) string
    ambientPath (1,1) string
end

weights = weights(:);
if numel(sources) ~= numel(weights) || any(weights < -1e-13) || ...
        abs(sum(weights)-1) > 1e-12
    error('AeroFlex:sched:AtomicWeights', ...
        'Atomic realization requires nonnegative partition-of-unity weights.');
end
if ~isfile(registryPath) || ~isfile(ambientPath)
    error('AeroFlex:sched:AtomicContractMissing', ...
        'The source registry and ambient embedding must both exist.');
end

registry = jsondecode(fileread(registryPath));
ambient = jsondecode(fileread(ambientPath));
sourceIds = strings(numel(sources),1);
entries = cell(numel(sources),1);
for k = 1:numel(sources)
    sourceIds(k) = sourceIdAtCoordinates(registry,sources(k).mu);
    resolved = AeroFlex.sched.resolveSourceContractRegistry( ...
        registryPath,sourceIds(k));
    entries{k} = resolved.sources;
end

queryDt = min(arrayfun(@(p) double(p.parConst.dt),sources));
commonDimension = double(ambient.commonAmbientDimension);
reducedOrder = size(readPythonH5(entries{1}.sourceArtifactPath, ...
    '/krylov_projectors/V'),2);
sourceData = repmat(emptySource(),numel(sources),1);

for k = 1:numel(sources)
    artifact = string(entries{k}.sourceArtifactPath);
    [R,E,mapKind] = ambientMaps( ...
        sourceIds(k),artifact,ambient,ambientPath);
    V = readPythonH5(artifact,'/krylov_projectors/V');
    W = readPythonH5(artifact,'/krylov_projectors/W');
    if size(V,1) ~= size(E,2) || ~isequal(size(V),size(W)) || ...
            size(V,2) ~= reducedOrder
        error('AeroFlex:sched:AtomicProjectorDimensions', ...
            'Source %s has incompatible V/W dimensions.',sourceIds(k));
    end
    if norm(W'*V-eye(reducedOrder),'fro') > 1e-10
        error('AeroFlex:sched:AtomicBiorthogonality', ...
            'Source %s violates W''*V=I.',sourceIds(k));
    end
    A = readPythonH5(artifact,'/modal_pre_krylov/A');
    B = readPythonH5(artifact,'/modal_pre_krylov/B');
    C = readPythonH5(artifact,'/modal_pre_krylov/C');
    D = readPythonH5(artifact,'/modal_pre_krylov/D');
    sourceDt = double(h5readatt(artifact,'/modal_pre_krylov','dt'));
    [A,B,C,D] = resampleBilinear(A,B,C,D,sourceDt,queryDt);

    S1 = chart.transforms(k).q1.Tinv;
    T1 = chart.transforms(k).q1.T;
    inputTransform = blkdiag(S1,S1,1,eye(4));
    B = B*inputTransform;
    C = T1*C;
    D = T1*D*inputTransform;

    sourceData(k) = struct('id',sourceIds(k),'R',R,'E',E, ...
        'V',V,'W',W,'A',A,'B',B,'C',C,'D',D, ...
        'mapKind',mapKind,'sourceDt',sourceDt,'Q',eye(reducedOrder), ...
        'primal',zeros(commonDimension,reducedOrder), ...
        'dual',zeros(commonDimension,reducedOrder));
    sourceData(k).primal = E*V;
    sourceData(k).dual = R'*W;
end

% The largest barycentric source is the deterministic local reference.
[~,reference] = max(weights);
referencePrimal = sourceData(reference).primal;
for k = 1:numel(sources)
    if k == reference
        Q = eye(reducedOrder);
    else
        [left,~,right] = svd(sourceData(k).primal'*referencePrimal,'econ');
        Q = left*right';
    end
    sourceData(k).Q = Q;
    sourceData(k).primal = sourceData(k).primal*Q;
    sourceData(k).dual = sourceData(k).dual*Q;
end

Vq = zeros(commonDimension,reducedOrder);
Wraw = zeros(commonDimension,reducedOrder);
for k = 1:numel(sources)
    Vq = Vq+weights(k)*sourceData(k).primal;
    Wraw = Wraw+weights(k)*sourceData(k).dual;
end
metric = Wraw'*Vq;
if rcond(metric) < 1e-12
    error('AeroFlex:sched:AtomicBlendedMetric', ...
        'The blended V/W metric is singular or ill-conditioned.');
end
Wq = Wraw/metric';
biorthogonalityError = norm(Wq'*Vq-eye(reducedOrder),'fro');
if biorthogonalityError > 1e-10
    error('AeroFlex:sched:AtomicQueryBiorthogonality', ...
        'The normalized query pair violates Wq''*Vq=I.');
end

Ad = zeros(reducedOrder);
Bd = zeros(reducedOrder,size(sourceData(1).B,2));
Cd = zeros(size(sourceData(1).C,1),reducedOrder);
Dd = zeros(size(sourceData(1).D));
transforms = repmat(struct('T',[],'Tinv',[],'condition',NaN), ...
    numel(sources),1);
for k = 1:numel(sources)
    item = sourceData(k);
    restrictedQuery = item.R*Vq;
    dualAtSource = Wq'*item.E;
    Ad = Ad+weights(k)*(dualAtSource*item.A*restrictedQuery);
    Bd = Bd+weights(k)*(dualAtSource*item.B);
    Cd = Cd+weights(k)*(item.C*restrictedQuery);
    Dd = Dd+weights(k)*item.D;
    T = Wq'*(item.E*item.V);
    if rcond(T) < 1e-12
        error('AeroFlex:sched:AtomicSourceTransform', ...
            'Source %s cannot be mapped into the query realization.',item.id);
    end
    transforms(k) = struct('T',T,'Tinv',T\eye(reducedOrder), ...
        'condition',cond(T));
end
[A,B,C,D] = discreteToContinuous(Ad,Bd,Cd,Dd,queryDt);

candidate = struct();
candidate.enabled = true;
candidate.architecture = 'source_only_atomic_fixed10_option_a';
candidate.sourceIds = sourceIds;
candidate.weights = weights;
candidate.referenceSource = sourceIds(reference);
candidate.queryDt = queryDt;
candidate.A = A; candidate.B = B; candidate.C = C; candidate.D = D;
candidate.Ad = Ad; candidate.Bd = Bd; candidate.Cd = Cd; candidate.Dd = Dd;
candidate.V = Vq; candidate.W = Wq;
candidate.transforms = transforms;
candidate.metrics = struct( ...
    'biorthogonalityFrobenius',biorthogonalityError, ...
    'metricCondition',cond(metric), ...
    'minimumSourceTransformRcond',min(arrayfun(@(x)rcond(x.T),transforms)), ...
    'maximumSourceTransformCondition',max([transforms.condition]));
candidate.sourceMaps = arrayfun(@(x)struct( ...
    'sourceId',x.id,'kind',x.mapKind,'sourceDt',x.sourceDt),sourceData);
candidate.registryVersion = string(registry.registryVersion);
candidate.ambientEmbeddingSha256 = string(ambient.ambientEmbeddingSha256);
candidate.validationTruthUsedForConstruction = false;
candidate.noTuningAfterValidation = true;
end

function id = sourceIdAtCoordinates(registry,mu)
coordinates = reshape([registry.sources.coordinates],2,[]).';
match = find(max(abs(coordinates-double(mu(:).')),[],2)<=1e-12);
if numel(match) ~= 1
    error('AeroFlex:sched:AtomicSourceIdentity', ...
        'Coordinates [%g %g] resolve to %d registry sources.', ...
        mu(1),mu(2),numel(match));
end
id = string(registry.sources(match).sourceId);
end

function [R,E,kind] = ambientMaps(sourceId,artifact,ambient,ambientPath)
sourceMaps = ambient.sourceMaps;
match = find(string({sourceMaps.sourceId}) == sourceId,1);
if ~isempty(match)
    indices = double(sourceMaps(match).sourceToCommonIndices(:))+1;
    sourceDimension = double(sourceMaps(match).sourceAmbientDimension);
    commonDimension = double(sourceMaps(match).commonAmbientDimension);
    E = sparse(indices,1:sourceDimension,1,commonDimension,sourceDimension);
    R = E';
    kind = "certified_exact_semantic_injection";
    return
end

logPath = fullfile(fileparts(fileparts(artifact)),'pazy_krylov_ROM.log');
if ~isfile(logPath)
    error('AeroFlex:sched:AtomicDelayEvidence', ...
        'No accepted gust-station log exists for source %s.',sourceId);
end
text = fileread(logPath);
token = regexp(text,'Gust monitoring station domain:\s*\[([^\]]+)\]', ...
    'tokens','once');
if isempty(token)
    error('AeroFlex:sched:AtomicDelayEvidence', ...
        'The accepted log has no gust monitoring station domain for %s.',sourceId);
end
stations = sscanf(token{1},'%f').';
common = double(ambientCommonStations());
Rd = zeros(numel(stations),numel(common));
for row = 1:numel(stations)
    exact = find(common == stations(row),1);
    if ~isempty(exact)
        Rd(row,exact) = 1;
        continue
    end
    upper = find(common > stations(row),1);
    if isempty(upper) || upper == 1
        error('AeroFlex:sched:AtomicDelayExtrapolation', ...
            'Source %s requires prohibited delay-map extrapolation.',sourceId);
    end
    lower = upper-1;
    fraction = (stations(row)-common(lower))/(common(upper)-common(lower));
    Rd(row,lower) = 1-fraction;
    Rd(row,upper) = fraction;
end
Ed = Rd'/(Rd*Rd');
R = blkdiag(Rd,speye(832));
E = blkdiag(Ed,speye(832));
if norm(full(R*E-eye(size(R,1))),'fro') > 1e-11
    error('AeroFlex:sched:AtomicDelayRoundTrip', ...
        'Option-A delay map for %s is not a right inverse.',sourceId);
end
kind = "certified_option_a_finite_overlay";

    function values = ambientCommonStations()
        persistent cached
        if isempty(cached)
            delayPath = fullfile(fileparts(char(ambientPath)),'DELAY_STATE_MAP.json');
            delay = jsondecode(fileread(delayPath));
            cached = delay.commonStationCoordinatesM(:).';
        end
        values = cached;
    end
end

function [Aq,Bq,Cq,Dq] = resampleBilinear(A,B,C,D,sourceDt,queryDt)
if abs(sourceDt-queryDt) <= 16*eps(max(sourceDt,queryDt))
    Aq=A; Bq=B; Cq=C; Dq=D;
    return
end
[Ac,Bc,Cc,Dc] = discreteToContinuous(A,B,C,D,sourceDt);
omega = 2/queryDt;
M = omega*eye(size(Ac))-Ac;
Aq = (omega*eye(size(Ac))+Ac)/M;
Bq = sqrt(2*omega)*(M\Bc);
Cq = sqrt(2*omega)*(Cc/M);
Dq = Dc+Cc*(M\Bc);
end

function [Ac,Bc,Cc,Dc] = discreteToContinuous(A,B,C,D,dt)
omega = 2/dt;
M = A+eye(size(A));
if rcond(M) < 1e-14
    error('AeroFlex:sched:AtomicBilinearCondition', ...
        'A source has an ill-conditioned bilinear inverse.');
end
Ac = omega*(A-eye(size(A)))/M;
Bc = sqrt(2*omega)*(M\B);
Cc = sqrt(2*omega)*(C/M);
Dc = D-C*(M\B);
end

function value = readPythonH5(path,dataset)
value = h5read(path,dataset);
if isstruct(value) && isfield(value,'r')
    value = double(value.r)+1i*double(value.i);
else
    value = double(value);
end
if ismatrix(value) && ~isvector(value), value = value.'; end
end

function value = emptySource()
value = struct('id',"",'R',[],'E',[],'V',[],'W',[],'A',[],'B',[], ...
    'C',[],'D',[],'mapKind',"",'sourceDt',NaN,'Q',[], ...
    'primal',[],'dual',[]);
end
