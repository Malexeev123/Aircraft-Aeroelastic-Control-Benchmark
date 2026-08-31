function candidate = buildFixedBasisQueryRealization( ...
        sources,weights,chart,registryPath)
%BUILDFIXEDBASISQUERYREALIZATION Reduce an interpolated fixed-node27 operator.
%   Source full operators are first converted to physical continuous time,
%   interpolated in the certified 1017-state fixed-node27 realization, and
%   discretized at a source-only query sample time. The accepted ordered
%   observability-side MIMO rational-Arnoldi reduction is then rerun.

arguments
    sources (1,:) struct
    weights (:,1) double
    chart (1,1) struct
    registryPath (1,1) string
end

weights=weights(:);
if numel(sources)~=numel(weights) || any(weights < -1e-13) || ...
        abs(sum(weights)-1)>1e-12
    error('AeroFlex:sched:FixedBasisWeights', ...
        'Fixed-basis query reduction requires convex weights.');
end
if ~isfile(registryPath)
    error('AeroFlex:sched:FixedBasisRegistryMissing', ...
        'The fixed-node27 registry is unavailable: %s.',registryPath);
end
registry=jsondecode(fileread(registryPath));
if string(registry.status)~="PASS" || ...
        string(registry.contractFamily)~="ACCEPTED_FIXED_NODE27_FIXED10"
    error('AeroFlex:sched:FixedBasisRegistryFamily', ...
        'The registry is not the certified fixed-node27 family.');
end

sourceIds=strings(numel(sources),1);
entries=cell(numel(sources),1);
sourceDt=zeros(numel(sources),1);
for k=1:numel(sources)
    sourceIds(k)=sourceIdAtCoordinates(registry,sources(k).mu);
    entry=registry.sources(string({registry.sources.sourceId})==sourceIds(k));
    if numel(entry)~=1 || string(entry.status)~="PASS"
        error('AeroFlex:sched:FixedBasisSourceResolution', ...
            'Source %s does not resolve to one certified contract.',sourceIds(k));
    end
    artifact=string(entry.selectedContract.path);
    if ~isfile(artifact) || fileHash(artifact)~=string(entry.selectedContract.sha256)
        error('AeroFlex:sched:FixedBasisSourceHash', ...
            'Source %s fixed-basis artifact is absent or stale.',sourceIds(k));
    end
    entries{k}=entry;
    sourceDt(k)=double(h5readatt(artifact, ...
        '/fixed_node27_premodal_common_discrete','dt'));
end

% The weighted source sample time is locked before development execution.
% It preserves exact-node identity and uses no held-out target quantity.
queryDt=weights.'*sourceDt;
commonDimension=double(registry.sources(1).dimensions.commonPremodal);
inputCount=double(registry.sources(1).dimensions.inputs);
outputCount=double(registry.sources(1).dimensions.outputs);
Ac=zeros(commonDimension); Bc=zeros(commonDimension,inputCount);
Cc=zeros(outputCount,commonDimension); Dc=zeros(outputCount,inputCount);
sourceV=cell(numel(sources),1);
for k=1:numel(sources)
    artifact=string(entries{k}.selectedContract.path);
    Ad=readH5(artifact,'/fixed_node27_premodal_common_discrete/A');
    Bd=readH5(artifact,'/fixed_node27_premodal_common_discrete/B');
    Cd=readH5(artifact,'/fixed_node27_premodal_common_discrete/C');
    Dd=readH5(artifact,'/fixed_node27_premodal_common_discrete/D');
    [Aci,Bci,Cci,Dci]=discreteToContinuous(Ad,Bd,Cd,Dd,sourceDt(k));

    % Structural trial and force-test maps use the already certified
    % physical chart. Gust, chi, and the two control/rate pairs retain their
    % accepted port ordering and units.
    S1=chart.transforms(k).q1.Tinv;
    T1=chart.transforms(k).q1.T;
    inputTransform=blkdiag(S1,S1,1,eye(4));
    Bci=Bci*inputTransform;
    Cci=T1*Cci;
    Dci=T1*Dci*inputTransform;
    Ac=Ac+weights(k)*Aci;
    Bc=Bc+weights(k)*Bci;
    Cc=Cc+weights(k)*Cci;
    Dc=Dc+weights(k)*Dci;
    sourceV{k}=readH5(artifact,'/projectors/V');
end
[Ad,Bd,Cd,Dd]=continuousToDiscrete(Ac,Bc,Cc,Dc,queryDt);
[Vq,reduction]=acceptedObservabilityReduction(Ad,Bd,Cd,Dd,4);
Wq=Vq;
[A,B,C,D]=discreteToContinuous( ...
    reduction.A,reduction.B,reduction.C,reduction.D,queryDt);

reducedOrder=size(Vq,2);
transforms=repmat(struct('T',[],'Tinv',[],'condition',NaN),numel(sources),1);
for k=1:numel(sources)
    T=Vq'*sourceV{k};
    if rcond(T)<1e-12
        error('AeroFlex:sched:FixedBasisStateTransport', ...
            'Source %s has singular transport into the query gauge.',sourceIds(k));
    end
    transforms(k)=struct('T',T,'Tinv',T\eye(reducedOrder), ...
        'condition',cond(T));
end

candidate=struct('enabled',true, ...
    'architecture','fixed_basis_full_operator_query_reduction_v1', ...
    'sourceIds',sourceIds,'weights',weights,'queryDt',queryDt, ...
    'sourceDt',sourceDt,'A',A,'B',B,'C',C,'D',D, ...
    'Ad',reduction.A,'Bd',reduction.B,'Cd',reduction.C,'Dd',reduction.D, ...
    'V',Vq,'W',Wq,'transforms',transforms, ...
    'referenceSource',sourceIds(find(weights==max(weights),1)), ...
    'metrics',struct( ...
        'biorthogonalityFrobenius',norm(Wq'*Vq-eye(reducedOrder),'fro'), ...
        'queryReductionRank',rank(Vq), ...
        'queryReductionCondition',cond(Vq), ...
        'minimumSourceTransformRcond',min(arrayfun(@(x)rcond(x.T),transforms)), ...
        'maximumSourceTransformCondition',max([transforms.condition]), ...
        'fullOperatorNorm',norm(Ac,'fro')), ...
    'registrySha256',char(fileHash(registryPath)), ...
    'physicalTimeConvention','continuous physical seconds', ...
    'queryDtRule','convex weighted source dt', ...
    'validationTruthUsedForConstruction',false, ...
    'noTuningAfterValidation',true);
end

function id=sourceIdAtCoordinates(registry,mu)
coordinates=reshape([registry.sources.coordinates],2,[]).';
match=find(max(abs(coordinates-double(mu(:).')),[],2)<=1e-12);
if numel(match)~=1
    error('AeroFlex:sched:FixedBasisSourceIdentity', ...
        'Coordinates [%g %g] resolve to %d fixed-basis sources.', ...
        mu(1),mu(2),numel(match));
end
id=string(registry.sources(match).sourceId);
end

function value=readH5(path,dataset)
value=h5read(path,dataset);
if isstruct(value) && isfield(value,'r')
    value=double(value.r)+1i*double(value.i);
else
    value=double(value);
end
% h5py persists matrices in C order; MATLAB exposes their dimensions in
% reverse order through h5read.
if ismatrix(value) && ~isvector(value), value=value.'; end
if ~isnumeric(value) || any(~isfinite(value),'all')
    error('AeroFlex:sched:FixedBasisDataset', ...
        'Dataset %s in %s is nonnumeric or nonfinite.',dataset,path);
end
end

function [V,reduced]=acceptedObservabilityReduction(A,B,C,D,r)
% This is the accepted SHARPy 2.4 ordered MGS algorithm at z=1. The
% observability seeds are the ten ordered columns of C'.
stateCount=size(A,1); seedCount=size(C,1);
targetCount=r*seedCount;
M=eye(stateCount)-A;
if rcond(M)<1e-14
    error('AeroFlex:sched:FixedBasisKrylovShift', ...
        'The query z=1 Krylov shift is ill-conditioned.');
end
V=zeros(stateCount,targetCount);
current=M'\C';
column=0;
for moment=1:r
    for seed=1:seedCount
        vector=current(:,seed);
        for previous=1:column
            vector=vector-V(:,previous)*(V(:,previous)'*vector);
        end
        magnitude=norm(vector);
        if ~isfinite(magnitude) || magnitude<=1e-12
            error('AeroFlex:sched:FixedBasisKrylovDeflation', ...
                'Query Krylov column %d deflated (norm %.3e).',column+1,magnitude);
        end
        column=column+1;
        V(:,column)=vector/magnitude;
    end
    if moment<r
        current=M'\V(:,column-seedCount+(1:seedCount));
    end
end
if column~=targetCount || norm(V'*V-eye(targetCount),'fro')>1e-10
    error('AeroFlex:sched:FixedBasisKrylovOrthogonality', ...
        'The query reduction did not construct the accepted 40-state basis.');
end
reduced=struct('A',V'*A*V,'B',V'*B,'C',C*V,'D',D);
end

function [Ac,Bc,Cc,Dc]=discreteToContinuous(A,B,C,D,dt)
omega=2/dt;
M=A+eye(size(A));
if rcond(M)<1e-14
    error('AeroFlex:sched:FixedBasisBilinearInverse', ...
        'A source bilinear inverse is ill-conditioned.');
end
Ac=omega*(A-eye(size(A)))/M;
Bc=sqrt(2*omega)*(M\B);
Cc=sqrt(2*omega)*(C/M);
Dc=D-C*(M\B);
end

function [A,B,C,D]=continuousToDiscrete(Ac,Bc,Cc,Dc,dt)
omega=2/dt;
M=omega*eye(size(Ac))-Ac;
if rcond(M)<1e-14
    error('AeroFlex:sched:FixedBasisBilinearForward', ...
        'The query bilinear map is ill-conditioned.');
end
A=(omega*eye(size(Ac))+Ac)/M;
B=sqrt(2*omega)*(M\Bc);
C=sqrt(2*omega)*(Cc/M);
D=Dc+Cc*(M\Bc);
end

function digest=fileHash(path)
fileId=fopen(path,'rb');
if fileId<0
    error('AeroFlex:sched:FixedBasisHashRead', ...
        'Unable to open fixed-basis artifact: %s.',path);
end
cleanup=onCleanup(@()fclose(fileId));
data=fread(fileId,Inf,'*uint8');
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(data,'int8'));
bytes=typecast(engine.digest(),'uint8');
digest=lower(string(reshape(dec2hex(bytes,2).',1,[])));
end
