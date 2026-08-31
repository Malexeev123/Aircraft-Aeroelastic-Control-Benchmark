function candidate = buildFrequencyLimitedGlobalParallelQueryRealization( ...
        sources,weights,registryPath)
%BUILDFREQUENCYLIMITEDGLOBALPARALLELQUERYREALIZATION Build V17 diagnostic.
%   V17 forms the global parallel transfer function from the certified source
%   tuples, then uses a 35--46 Hz frequency-limited balanced realization.

arguments
    sources (1,:) struct
    weights (:,1) double
    registryPath (1,1) string
end

weights=weights(:);
if numel(sources)~=numel(weights) || any(weights < -1e-13) || ...
        abs(sum(weights)-1)>1e-12
    error('AeroFlex:sched:V17Weights', ...
        'V17 requires nonnegative partition-of-unity weights.');
end
if ~isfile(registryPath)
    error('AeroFlex:sched:V17Registry', ...
        'The V17 certified source registry is absent.');
end

registry=jsondecode(fileread(registryPath));
if string(registry.status)~="PASS" || ...
        string(registry.contractFamily)~="ACCEPTED_FIXED_NODE27_FIXED10"
    error('AeroFlex:sched:V17Registry', ...
        'V17 requires a passing fixed-node27/fixed-10 registry.');
end
allSources=loadRegistrySources(registry);
allIds=string({allSources.name});
activeIds=strings(numel(sources),1);
for k=1:numel(sources)
    activeIds(k)=sourceIdAtCoordinates(registry,sources(k).mu);
end
if numel(unique(activeIds))~=numel(activeIds) || ...
        ~all(ismember(activeIds,allIds))
    error('AeroFlex:sched:V17ActiveSources', ...
        'Active sources do not resolve uniquely in the V17 registry.');
end
globalWeights=zeros(numel(allSources),1);
for k=1:numel(activeIds)
    globalWeights(allIds==activeIds(k))=weights(k);
end

[Aglobal,Bglobal,Cglobal,Dglobal,fields,offsets]= ...
    assembleGlobalParallel(allSources,globalWeights);
[V,W,singularValues]=frequencyLimitedBasis(Aglobal,Bglobal,Cglobal);
overlap=W'*V;
if rcond(overlap)<1e-10
    error('AeroFlex:sched:V17BiorthogonalityCondition', ...
        'V17 frequency-limited biorthogonal overlap is ill-conditioned.');
end
W=W/overlap';
if norm(W'*V-eye(40),'fro')>1e-12
    error('AeroFlex:sched:V17Biorthogonality', ...
        'V17 frequency-limited projection is not biorthogonal.');
end

candidateA=real(W'*Aglobal*V);
candidateB=real(W'*Bglobal);
candidateC=real(Cglobal*V);
candidateD=real(Dglobal);
candidateFields=projectFields(fields,V,W);
oneHotTransportCondition=NaN;
oneHotInverseResidual=NaN;

if nnz(globalWeights>1e-14)==1
    sourceIndex=find(globalWeights>1e-14,1);
    rows=offsets(sourceIndex)+(1:40);
    transport=V(rows,:);
    if rcond(transport)<1e-10
        error('AeroFlex:sched:V17IdentityTransport', ...
            'The V17 one-hot source transport is singular.');
    end
    oneHotTransportCondition=cond(transport);
    oneHotInverseResidual=norm(transport*(W(rows,:).')-eye(40),'fro');
    candidateA=real(transport*candidateA/transport);
    candidateB=real(transport*candidateB);
    candidateC=real(candidateC/transport);
    candidateFields=transportOneHotFields(candidateFields,transport);
end

candidate=struct('enabled',true, ...
    'architecture','frequency_limited_global_parallel_transfer_function_diagnostic_v17', ...
    'sourceIds',activeIds,'weights',weights, ...
    'referenceSource',"GLOBAL_DIRECT_SUM_FIXED_COORDINATE", ...
    'A',candidateA,'B',candidateB, ...
    'C',candidateC,'D',candidateD, ...
    'qGammaEquilibrium',candidateFields.qGammaEquilibrium, ...
    'qGammaAffine',candidateFields.qGammaAffine, ...
    'Bchi',candidateFields.Bchi,'rootC',candidateFields.rootC, ...
    'rootD',candidateFields.rootD,'V',V,'W',W, ...
    'metrics',struct('frequencyBandHz',[35 46], ...
        'quadratureNodeCount',16,'biorthogonalityFrobenius', ...
        norm(W'*V-eye(40),'fro'),'hsv40',singularValues(40), ...
        'hsv41',singularValues(min(41,numel(singularValues))), ...
        'oneHotTransportCondition',oneHotTransportCondition, ...
        'oneHotInverseResidual',oneHotInverseResidual), ...
    'registrySha256',char(fileHash(registryPath)), ...
    'fieldAssemblyOrder','global-parallel-frequency-limited-project', ...
    'validationTruthUsedForConstruction',false, ...
    'noReducedInverseTransport',true);
end

function sources=loadRegistrySources(registry)
repo=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
sources=repmat(struct('name',"",'package',struct()),1,numel(registry.sources));
for k=1:numel(registry.sources)
    entry=registry.sources(k);
    path=fullfile(repo,char(entry.acceptedMatlabPackage.path));
    if ~isfile(path) || fileHash(path)~=string(entry.acceptedMatlabPackage.sha256)
        error('AeroFlex:sched:V17SourceHash', ...
            'Certified package for %s is absent or stale.',entry.sourceId);
    end
    data=load(path);
    if isfield(data,'p5'), package=data.p5;
    elseif isfield(data,'package'), package=data.package;
    else, error('AeroFlex:sched:V17PackageSchema', ...
            'Certified package for %s has no p5/package field.',entry.sourceId);
    end
    required={'trim','cfg','beam','aero','base'};
    if ~all(isfield(data,required))
        error('AeroFlex:sched:V17PackageSchema', ...
            'Certified package for %s lacks complete plant data.',entry.sourceId);
    end
    package=AeroFlex.sched.installP5R2T2Anchor( ...
        package,data.trim,"SOURCE_TRAINING_T2_ANCHOR");
    sources(k)=struct('name',string(entry.sourceId),'package',package);
end
end

function [Aglobal,Bglobal,Cglobal,Dglobal,fields,offsets]= ...
        assembleGlobalParallel(sources,weights)
count=numel(sources); tuples=cell(1,count); offsets=40*(0:count-1);
for k=1:count, tuples{k}=tupleFromPackage(sources(k).package); end
aBlocks=cellfun(@(item)item.A,tuples,'UniformOutput',false);
bBlocks=cellfun(@(item)item.B,tuples,'UniformOutput',false);
Aglobal=blkdiag(aBlocks{:});
Bglobal=vertcat(bBlocks{:});
Cglobal=[]; Dglobal=zeros(size(tuples{1}.D));
fields=struct('qGammaEquilibrium',[],'qGammaAffine',[], ...
    'Bchi',[],'rootC',[],'rootD',zeros(size(sources(1).package.p5.r2.rootD)));
for k=1:count
    weight=weights(k); tuple=tuples{k}; package=sources(k).package;
    Cglobal=[Cglobal weight*tuple.C]; %#ok<AGROW>
    Dglobal=Dglobal+weight*tuple.D;
    idx=package.idx;
    fields.qGammaEquilibrium=[fields.qGammaEquilibrium;package.x_eq(idx.qGam)];
    fields.qGammaAffine=[fields.qGammaAffine; ...
        package.parConst.affineOffset(idx.qGam)];
    fields.Bchi=[fields.Bchi;package.base.FM.Bchi];
    fields.rootC=[fields.rootC weight*package.p5.r2.rootC];
    fields.rootD=fields.rootD+weight*package.p5.r2.rootD;
end
end

function tuple=tupleFromPackage(package)
fm=package.aero.forceMap; modeCount=numel(package.idx.q1);
tuple=struct('A',fm.A_Gamma, ...
    'B',[fm.B0 fm.B1 fm.Bw fm.B_delta fm.B_ddelta], ...
    'C',fm.C_Gamma, ...
    'D',[fm.D0 fm.D1 fm.Dw fm.D_delta fm.D_ddelta]);
if ~isequal(size(tuple.A),[40 40]) || size(tuple.B,1)~=40 || ...
        size(tuple.C,2)~=40 || size(tuple.B,2)~=2*modeCount+5 || ...
        ~isequal(size(tuple.D),[size(tuple.C,1) size(tuple.B,2)])
    error('AeroFlex:sched:V17TupleDimensions', ...
        'Certified source tuple violates the fixed-10 port contract.');
end
end

function [V,W,singularValues]=frequencyLimitedBasis(A,B,C)
[nodes,quadratureWeights]=gaussLegendre(16);
frequencies=40.5+5.5*nodes; dOmega=2*pi*5.5*quadratureWeights;
stateCount=size(A,1); X=[]; Y=[];
for frequencyIndex=1:numel(frequencies)
    omega=2*pi*frequencies(frequencyIndex);
    scale=sqrt(dOmega(frequencyIndex)/(2*pi));
    resolvent=1i*omega*eye(stateCount)-A;
    inputs=scale*(resolvent\B);
    outputs=scale*(C/resolvent);
    X=[X sqrt(2)*real(inputs) sqrt(2)*imag(inputs)]; %#ok<AGROW>
    Y=[Y;sqrt(2)*real(outputs);sqrt(2)*imag(outputs)]; %#ok<AGROW>
end
cross=Y*X;
imaginaryRatio=norm(imag(cross),'fro')/max(1,norm(real(cross),'fro'));
if imaginaryRatio>1e-10
    error('AeroFlex:sched:V17ComplexPairing', ...
        'Conjugate quadrature pairing leaves a non-negligible imaginary cross Gramian.');
end
[U,S,right]=svd(real(cross),'econ'); singularValues=diag(S);
if numel(singularValues)<40 || singularValues(40)<=eps(singularValues(1))
    error('AeroFlex:sched:V17Rank', ...
        'The frequency-limited global parallel realization has rank below 40.');
end
sqrtInverse=diag(1./sqrt(singularValues(1:40)));
V=real(X*right(:,1:40)*sqrtInverse);
W=real(Y'*U(:,1:40)*sqrtInverse);
end

function fields=projectFields(fields,V,W)
fields.qGammaEquilibrium=real(W'*fields.qGammaEquilibrium);
fields.qGammaAffine=real(W'*fields.qGammaAffine);
fields.Bchi=real(W'*fields.Bchi);
fields.rootC=real(fields.rootC*V);
end

function fields=transportOneHotFields(fields,transport)
fields.qGammaEquilibrium=real(transport*fields.qGammaEquilibrium);
fields.qGammaAffine=real(transport*fields.qGammaAffine);
fields.Bchi=real(transport*fields.Bchi);
fields.rootC=real(fields.rootC/transport);
end

function [nodes,weights]=gaussLegendre(count)
index=(1:count-1).'; beta=index./sqrt(4*index.^2-1);
[vectors,values]=eig(diag(beta,1)+diag(beta,-1));
[nodes,order]=sort(diag(values)); weights=2*vectors(1,order).^2;
end

function digest=fileHash(path)
engine=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
fid=fopen(path,'rb');
assert(fid>=0,'AeroFlex:sched:V17HashOpen','Cannot read %s.',path);
cleanup=onCleanup(@() fclose(fid));
while ~feof(fid)
    buffer=fread(fid,1024*1024,'*uint8');
    if isempty(buffer), break; end
    engine.update(typecast(buffer(:),'int8'));
end
bytes=typecast(engine.digest(),'uint8');
digest=lower(string(reshape(dec2hex(bytes,2).',1,[])));
clear cleanup
end

function id=sourceIdAtCoordinates(registry,mu)
coordinates=reshape([registry.sources.coordinates],2,[]).';
match=find(max(abs(coordinates-double(mu(:).')),[],2)<=1e-12);
if numel(match)~=1
    error('AeroFlex:sched:V17SourceIdentity', ...
        'Coordinates [%g %g] resolve to %d sources.',mu(1),mu(2),numel(match));
end
id=string(registry.sources(match).sourceId);
end
