function candidate = buildAlignedReducedRealization( ...
        sources,weights,chart,registryPath,ambientPath)
%BUILDALIGNEDREDUCEDREALIZATION Interpolate in the fixed node27 aero chart.
%   Every certified source tuple is converted to physical continuous time
%   and transported by x_ref = W_ref^H V_source x_source. The same map is
%   returned for equilibrium, recovery, and nonlinear-package transport.

arguments
    sources (1,:) struct
    weights (:,1) double
    chart (1,1) struct
    registryPath (1,1) string
    ambientPath (1,1) string
end
weights=weights(:);
if numel(sources)~=numel(weights) || any(weights < -1e-13) || ...
        abs(sum(weights)-1)>1e-12
    error('AeroFlex:sched:AlignedWeights', ...
        'Aligned realization requires nonnegative partition-of-unity weights.');
end
registry=jsondecode(fileread(registryPath));
ambient=jsondecode(fileread(ambientPath));
referenceId="node27";
reference=AeroFlex.sched.resolveSourceContractRegistry( ...
    registryPath,referenceId).sources;
[Rref,Eref]=commonMaps(referenceId,string(reference.sourceArtifactPath), ...
    ambient,ambientPath);
Vref=readPythonH5(reference.sourceArtifactPath,'/krylov_projectors/V');
Wref=readPythonH5(reference.sourceArtifactPath,'/krylov_projectors/W');
if norm(Wref'*Vref-eye(size(Vref,2)),'fro')>1e-10
    error('AeroFlex:sched:AlignedReferencePair', ...
        'The node27 reference violates W''*V=I.');
end
Vcommon=Eref*Vref;
Wcommon=Rref'*Wref;

sourceIds=strings(numel(sources),1);
data=repmat(emptySource(),numel(sources),1);
for k=1:numel(sources)
    sourceIds(k)=sourceIdAtCoordinates(registry,sources(k).mu);
    entry=AeroFlex.sched.resolveSourceContractRegistry( ...
        registryPath,sourceIds(k)).sources;
    artifact=string(entry.sourceArtifactPath);
    [~,E]=commonMaps(sourceIds(k),artifact,ambient,ambientPath);
    V=readPythonH5(artifact,'/krylov_projectors/V');
    W=readPythonH5(artifact,'/krylov_projectors/W');
    if norm(W'*V-eye(size(V,2)),'fro')>1e-10
        error('AeroFlex:sched:AlignedSourcePair', ...
            'Source %s violates W''*V=I.',sourceIds(k));
    end
    T=Wcommon'*(E*V);
    if rcond(T)<1e-12
        error('AeroFlex:sched:AlignedSourceTransform', ...
            'Source %s is singular in the node27 chart.',sourceIds(k));
    end
    Tinv=T\eye(size(T));
    A=readPythonH5(artifact,'/krylov_discrete/A');
    B=readPythonH5(artifact,'/krylov_discrete/B');
    C=readPythonH5(artifact,'/krylov_discrete/C');
    D=readPythonH5(artifact,'/krylov_discrete/D');
    dt=double(h5readatt(artifact,'/krylov_discrete','dt'));
    [A,B,C,D]=discreteToContinuous(A,B,C,D,dt);

    S1=chart.transforms(k).q1.Tinv;
    T1=chart.transforms(k).q1.T;
    inputTransform=blkdiag(S1,S1,1,eye(4));
    B=B*inputTransform;
    C=T1*C;
    D=T1*D*inputTransform;
    data(k)=struct('id',sourceIds(k),'A',T*A*Tinv, ...
        'B',T*B,'C',C*Tinv,'D',D,'dt',dt,'T',T,'Tinv',Tinv, ...
        'condition',cond(T));
end

A=weighted(data,'A',weights);
B=weighted(data,'B',weights);
C=weighted(data,'C',weights);
D=weighted(data,'D',weights);
transforms=arrayfun(@(x)struct('T',x.T,'Tinv',x.Tinv, ...
    'condition',x.condition),data);
candidate=struct('enabled',true, ...
    'architecture','source_only_aligned_reduced_node27_v1', ...
    'sourceIds',sourceIds,'weights',weights,'referenceSource',referenceId, ...
    'queryDt',min(arrayfun(@(p)double(p.parConst.dt),sources)), ...
    'A',A,'B',B,'C',C,'D',D,'Ad',[],'Bd',[],'Cd',[],'Dd',[], ...
    'V',Vcommon,'W',Wcommon,'transforms',transforms, ...
    'metrics',struct('biorthogonalityFrobenius', ...
        norm(Wcommon'*Vcommon-eye(size(Vcommon,2)),'fro'), ...
        'metricCondition',1, ...
        'minimumSourceTransformRcond',min(arrayfun(@(x)rcond(x.T),data)), ...
        'maximumSourceTransformCondition',max([data.condition])), ...
    'sourceMaps',arrayfun(@(x)struct('sourceId',x.id, ...
        'kind','fixed_node27_biorthogonal_overlap','sourceDt',x.dt),data), ...
    'registryVersion',string(registry.registryVersion), ...
    'ambientEmbeddingSha256',string(ambient.ambientEmbeddingSha256), ...
    'validationTruthUsedForConstruction',false, ...
    'noTuningAfterValidation',true);
end

function value=weighted(data,field,weights)
value=zeros(size(data(1).(field)));
for k=1:numel(data), value=value+weights(k)*data(k).(field); end
end

function id=sourceIdAtCoordinates(registry,mu)
coordinates=reshape([registry.sources.coordinates],2,[]).';
match=find(max(abs(coordinates-double(mu(:).')),[],2)<=1e-12);
if numel(match)~=1
    error('AeroFlex:sched:AlignedSourceIdentity', ...
        'Coordinates [%g %g] resolve to %d sources.',mu(1),mu(2),numel(match));
end
id=string(registry.sources(match).sourceId);
end

function [R,E]=commonMaps(sourceId,artifact,ambient,ambientPath)
maps=ambient.sourceMaps;
match=find(string({maps.sourceId})==sourceId,1);
if ~isempty(match)
    indices=double(maps(match).sourceToCommonIndices(:))+1;
    ns=double(maps(match).sourceAmbientDimension);
    nc=double(maps(match).commonAmbientDimension);
    E=sparse(indices,1:ns,1,nc,ns); R=E'; return
end
logPath=fullfile(fileparts(fileparts(artifact)),'pazy_krylov_ROM.log');
text=fileread(logPath);
token=regexp(text,'Gust monitoring station domain:\s*\[([^\]]+)\]', ...
    'tokens','once');
if isempty(token)
    error('AeroFlex:sched:AlignedDelayEvidence', ...
        'No accepted delay evidence exists for %s.',sourceId);
end
stations=sscanf(token{1},'%f').';
delay=jsondecode(fileread(fullfile(fileparts(char(ambientPath)), ...
    'DELAY_STATE_MAP.json')));
common=double(delay.commonStationCoordinatesM(:).');
Rd=zeros(numel(stations),numel(common));
for row=1:numel(stations)
    exact=find(common==stations(row),1);
    if ~isempty(exact), Rd(row,exact)=1; continue, end
    upper=find(common>stations(row),1);
    if isempty(upper) || upper==1
        error('AeroFlex:sched:AlignedDelayExtrapolation', ...
            'Source %s requires delay extrapolation.',sourceId);
    end
    lower=upper-1;
    fraction=(stations(row)-common(lower))/(common(upper)-common(lower));
    Rd(row,lower)=1-fraction; Rd(row,upper)=fraction;
end
Ed=Rd'/(Rd*Rd');
R=blkdiag(Rd,speye(832)); E=blkdiag(Ed,speye(832));
end

function [Ac,Bc,Cc,Dc]=discreteToContinuous(A,B,C,D,dt)
omega=2/dt; M=A+eye(size(A));
if rcond(M)<1e-14
    error('AeroFlex:sched:AlignedBilinearCondition', ...
        'A source has an ill-conditioned bilinear inverse.');
end
Ac=omega*(A-eye(size(A)))/M;
Bc=sqrt(2*omega)*(M\B);
Cc=sqrt(2*omega)*(C/M);
Dc=D-C*(M\B);
end

function value=readPythonH5(path,dataset)
value=h5read(path,dataset);
if isstruct(value) && isfield(value,'r')
    value=double(value.r)+1i*double(value.i);
else
    value=double(value);
end
if ismatrix(value) && ~isvector(value), value=value.'; end
end

function value=emptySource()
value=struct('id',"",'A',[],'B',[],'C',[],'D',[],'dt',NaN, ...
    'T',[],'Tinv',[],'condition',NaN);
end
