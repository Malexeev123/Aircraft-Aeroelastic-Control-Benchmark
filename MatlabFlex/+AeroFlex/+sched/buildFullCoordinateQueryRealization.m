function candidate = buildFullCoordinateQueryRealization( ...
        sources,weights,chart,registryPath,fieldRoot,options)
%BUILDFULLCOORDINATEQUERYREALIZATION Lift, interpolate, and project all fields.
%   Every qGamma-dependent field is assembled in the shared 1017-state
%   fixed-node27 realization. A single source-only biorthogonal query pair
%   then projects dynamics, equilibrium, affine, chi, and recovery fields.

arguments
    sources (1,:) struct
    weights (:,1) double
    chart (1,1) struct
    registryPath (1,1) string
    fieldRoot (1,1) string
    options (1,1) struct = struct()
end

weights=weights(:);
if numel(sources)~=numel(weights) || any(weights < -1e-13) || ...
        abs(sum(weights)-1)>1e-12
    error('AeroFlex:sched:FullCoordinateWeights', ...
        'Full-coordinate construction requires convex weights.');
end
if ~isfile(registryPath) || ~isfolder(fieldRoot)
    error('AeroFlex:sched:FullCoordinateInputs', ...
        'The fixed-node27 registry and full-field sidecars are required.');
end
registry=jsondecode(fileread(registryPath));
if string(registry.status)~="PASS" || ...
        string(registry.contractFamily)~="ACCEPTED_FIXED_NODE27_FIXED10"
    error('AeroFlex:sched:FullCoordinateFamily', ...
        'The registry is not the certified fixed-node27 family.');
end
cacheImmutableMetadata=isfield(options,'cacheImmutableMetadata') && ...
    logical(options.cacheImmutableMetadata);
cacheImmutableSourceData=isfield(options,'cacheImmutableSourceData') && ...
    logical(options.cacheImmutableSourceData);
preparedRuntimeOwner=isfield(options,'preparedRuntimeOwner') && ...
    logical(options.preparedRuntimeOwner);
projectedTustinCondensation=isfield(options, ...
    'projectedTustinCondensation') && ...
    logical(options.projectedTustinCondensation);
if preparedRuntimeOwner
    requiredPrepared={'preparedRuntimeProfile','preparedRuntimeChangeId'};
    missingPrepared=requiredPrepared(~isfield(options,requiredPrepared));
    if ~isempty(missingPrepared)
        error('AeroFlex:sched:PreparedRuntimeScope', ...
            'Prepared runtime options are missing: %s.', ...
            strjoin(missingPrepared,', '));
    end
    preparedProfile=string(options.preparedRuntimeProfile);
    preparedChangeId=string(options.preparedRuntimeChangeId);
    if ~isscalar(preparedProfile) || ...
            ~ismember(preparedProfile,["formal_case_b","formal_case_c"]) || ...
            ~isscalar(preparedChangeId) || preparedChangeId~= ...
            "phase18c-v17a-casebc-prepared-scheduled-runtime-owner-v1" || ...
            ~cacheImmutableSourceData
        error('AeroFlex:sched:PreparedRuntimeScope', ...
            ['Prepared source operators require the approved Case-B/C ', ...
             'profile, change identifier, and immutable source cache.']);
    end
else
    preparedProfile="disabled";
    preparedChangeId="";
end
if projectedTustinCondensation && ~preparedRuntimeOwner
    error('AeroFlex:sched:ProjectedTustinScope', ...
        ['Projected Tustin condensation requires the approved prepared ', ...
         'Case-B/C runtime owner.']);
end

count=numel(sources); ids=strings(count,1); entries=cell(count,1);
sourceDt=zeros(count,1); sourceData=repmat(emptySource(),count,1);
for k=1:count
    ids(k)=sourceIdAtCoordinates(registry,sources(k).mu);
    entry=registry.sources(string({registry.sources.sourceId})==ids(k));
    if numel(entry)~=1 || string(entry.status)~="PASS"
        error('AeroFlex:sched:FullCoordinateResolution', ...
            'Source %s does not resolve to one certified contract.',ids(k));
    end
    artifact=resolveRepositoryPath(string(entry.selectedContract.path));
    sidecar=fullfile(fieldRoot,ids(k)+".h5");
    metadata=verifiedImmutableMetadata(artifact,sidecar, ...
        string(entry.selectedContract.sha256),ids(k), ...
        cacheImmutableMetadata || cacheImmutableSourceData);
    entries{k}=entry;
    sourceDt(k)=metadata.dt;
end

physicsNormalized=isfield(options,'physicsNormalized') && ...
    logical(options.physicsNormalized);
orderedSchur=isfield(options,'orderedSchur') && logical(options.orderedSchur);
invariantCluster=isfield(options,'invariantCluster') && ...
    logical(options.invariantCluster);
jointCoupling=isfield(options,'jointCoupling') && logical(options.jointCoupling);
jointDiscrete=isfield(options,'jointDiscrete') && logical(options.jointDiscrete);
jointDiscreteCouplingOnly=isfield(options,'jointDiscreteCouplingOnly') && ...
    logical(options.jointDiscreteCouplingOnly);
fullCoupledCommonCoordinate=isfield(options,'fullCoupledCommonCoordinate') && ...
    logical(options.fullCoupledCommonCoordinate);
branchPreservingStructural=isfield(options,'branchPreservingStructural') && ...
    logical(options.branchPreservingStructural);
fullCoupledPartitionPreserving=isfield(options,'fullCoupledPartitionPreserving') && ...
    logical(options.fullCoupledPartitionPreserving);
unpartitionedCoupledInterface=isfield(options,'unpartitionedCoupledInterface') && ...
    logical(options.unpartitionedCoupledInterface);
branchAwareCrossCoupling=isfield(options,'branchAwareCrossCoupling') && ...
    logical(options.branchAwareCrossCoupling);
fullCoupledOrderedSchur=isfield(options,'fullCoupledOrderedSchur') && ...
    logical(options.fullCoupledOrderedSchur);
selectiveHighBandCoupled=isfield(options,'selectiveHighBandCoupled') && ...
    logical(options.selectiveHighBandCoupled);
multiBandCoupled=isfield(options,'multiBandCoupled') && logical(options.multiBandCoupled);
if orderedSchur && ~physicsNormalized
    error('AeroFlex:sched:OrderedSchurNormalization', ...
        'Ordered-Schur construction requires physics-normalized operators.');
end
if (jointCoupling || jointDiscrete || jointDiscreteCouplingOnly || ...
        fullCoupledCommonCoordinate || branchPreservingStructural || ...
        fullCoupledPartitionPreserving || unpartitionedCoupledInterface || ...
        branchAwareCrossCoupling || fullCoupledOrderedSchur || selectiveHighBandCoupled || multiBandCoupled) && ...
        (~physicsNormalized || orderedSchur || invariantCluster)
    error('AeroFlex:sched:JointCouplingOptions', ...
        ['Joint structural-aerodynamic assembly requires the direct ', ...
        'physics-normalized path without Schur or invariant overrides.']);
end
if nnz([jointCoupling,jointDiscrete,jointDiscreteCouplingOnly, ...
        fullCoupledCommonCoordinate,branchPreservingStructural, ...
        fullCoupledPartitionPreserving,unpartitionedCoupledInterface, ...
        branchAwareCrossCoupling,fullCoupledOrderedSchur,selectiveHighBandCoupled,multiBandCoupled])>1
    error('AeroFlex:sched:JointCouplingOptions', ...
        'Only one joint-coupling construction may be active at a time.');
end
if physicsNormalized
    if ~isfield(options,'query') || ~isequal(size(options.query),[1 2]) || ...
            options.query(1)<=0
        error('AeroFlex:sched:PhysicsNormalizedQuery', ...
            'Physics-normalized interpolation requires [U alpha] with U>0.');
    end
    sourceSpeed=arrayfun(@(item)double(item.mu(1)),sources).';
    sourceChord=sourceDt.*sourceSpeed;
    referenceChord=0.025;
    if max(abs(sourceChord-referenceChord))>1e-12
        error('AeroFlex:sched:PhysicsNormalizedTimeIdentity', ...
            'Source dt*U does not reproduce the accepted 0.025 m chord.');
    end
    queryDt=referenceChord/double(options.query(1));
else
    queryDt=weights.'*sourceDt;
    referenceChord=NaN;
end
commonDimension=double(registry.sources(1).dimensions.commonPremodal);
inputCount=double(registry.sources(1).dimensions.inputs);
outputCount=double(registry.sources(1).dimensions.outputs);
sourceCacheStats=struct('enabled',cacheImmutableSourceData, ...
    'preparedRuntimeOwner',preparedRuntimeOwner, ...
    'preparedRuntimeProfile',preparedProfile,'hits',0,'misses',0, ...
    'preparedHits',0,'preparedMisses',0,'rawH5Reads',0, ...
    'preparationSeconds',0,'wallSeconds',0);

for k=1:count
    artifact=resolveRepositoryPath(string(entries{k}.selectedContract.path));
    sidecar=fullfile(fieldRoot,ids(k)+".h5");
    [raw,rawInfo]=verifiedImmutableSourceData(artifact,sidecar, ...
        string(entries{k}.selectedContract.sha256),ids(k), ...
        cacheImmutableSourceData,preparedRuntimeOwner);
    sourceCacheStats.hits=sourceCacheStats.hits+double(rawInfo.hit);
    sourceCacheStats.misses=sourceCacheStats.misses+double(~rawInfo.hit);
    sourceCacheStats.rawH5Reads=sourceCacheStats.rawH5Reads+rawInfo.rawH5Reads;
    sourceCacheStats.preparedHits=sourceCacheStats.preparedHits+ ...
        double(rawInfo.preparedHit);
    sourceCacheStats.preparedMisses=sourceCacheStats.preparedMisses+ ...
        double(rawInfo.preparedMiss);
    sourceCacheStats.preparationSeconds=sourceCacheStats.preparationSeconds+ ...
        rawInfo.preparationSeconds;
    sourceCacheStats.wallSeconds=sourceCacheStats.wallSeconds+rawInfo.wallSeconds;
    Ad=raw.Ad; Bd=raw.Bd; Cd=raw.Cd; Dd=raw.Dd;
    if preparedRuntimeOwner
        Ac=raw.prepared.Ac; Bc=raw.prepared.Bc;
        Cc=raw.prepared.Cc; Dc=raw.prepared.Dc;
    else
        [Ac,Bc,Cc,Dc]=discreteToContinuous( ...
            Ad,Bd,Cd,Dd,sourceDt(k));
    end
    S1=chart.transforms(k).q1.Tinv; T1=chart.transforms(k).q1.T;
    T2=chart.transforms(k).q2.T;
    inputTransform=blkdiag(S1,S1,1,eye(4));
    Bc=Bc*inputTransform; Cc=T1*Cc; Dc=T1*Dc*inputTransform;
    joint=struct();
    if jointCoupling
        Nm=numel(sources(k).idx.q1);
        omegaSource=-T2*sources(k).L(sources(k).idx.q2,sources(k).idx.q1)*S1;
        if rcond(omegaSource)<1e-12
            error('AeroFlex:sched:JointCouplingOmega', ...
                'Source %s has an ill-conditioned physical Omega block.',ids(k));
        end
        joint.q1q1=sources(k).parConst.scaleA*sources(k).parConst.t_inf* ...
            Dc(:,Nm+(1:Nm));
        joint.q1q2=omegaSource-sources(k).parConst.scaleA* ...
            (Dc(:,1:Nm)/omegaSource);
        joint.q1Gamma=sources(k).parConst.scaleA*sources(k).parConst.t_inf*Cc;
        joint.gammaQ1=Bc(:,Nm+(1:Nm));
        joint.gammaQ2=-(Bc(:,1:Nm)/omegaSource)/sources(k).parConst.t_inf;
        joint.gammaGamma=Ac/sources(k).parConst.t_inf;
    end
    if jointDiscrete || jointDiscreteCouplingOnly || fullCoupledCommonCoordinate || ...
            branchPreservingStructural || fullCoupledPartitionPreserving || ...
            unpartitionedCoupledInterface || branchAwareCrossCoupling || fullCoupledOrderedSchur || ...
            selectiveHighBandCoupled || multiBandCoupled
        Nm=numel(sources(k).idx.q1);
        S2=chart.transforms(k).q2.Tinv;
        sourceL=sources(k).L;
        sourceIdx=sources(k).idx;
        subsystem=zeros(2*Nm+commonDimension);
        subsystem(1:Nm,1:Nm)=T1*sourceL(sourceIdx.q1,sourceIdx.q1)*S1;
        subsystem(1:Nm,Nm+(1:Nm))=T1*sourceL(sourceIdx.q1,sourceIdx.q2)*S2;
        subsystem(Nm+(1:Nm),1:Nm)= ...
            T2*sourceL(sourceIdx.q2,sourceIdx.q1)*S1;
        subsystem(Nm+(1:Nm),Nm+(1:Nm))= ...
            T2*sourceL(sourceIdx.q2,sourceIdx.q2)*S2;
        gamma=2*Nm+(1:commonDimension);
        subsystem(1:Nm,gamma)=sources(k).parConst.scaleA*Cc;
        subsystem(gamma,1:Nm)=Bc(:,Nm+(1:Nm));
        omegaSource=-T2*sourceL(sourceIdx.q2,sourceIdx.q1)*S1;
        if rcond(omegaSource)<1e-12
            error('AeroFlex:sched:JointDiscreteOmega', ...
                'Source %s has an ill-conditioned physical Omega block.',ids(k));
        end
        subsystem(gamma,Nm+(1:Nm))= ...
            -(Bc(:,1:Nm)/omegaSource)/sources(k).parConst.t_inf;
        subsystem(gamma,gamma)=Ac/sources(k).parConst.t_inf;
        joint=struct('normalizedSubsystem',sourceDt(k)*subsystem);
    end
    reducedA=raw.reducedA;
    if orderedSchur
        [schurGauge,schurBlocks]=orderedRealSchurGauge(reducedA);
    else
        schurGauge=eye(size(reducedA)); schurBlocks=[];
    end
    if physicsNormalized
        % Remove only the exact inverse-Tustin rate factors. Dependence on
        % alpha and the source equilibrium remains in the residual fields.
        Ac=sourceDt(k)*Ac;
        Bc=sqrt(sourceDt(k))*Bc;
        Cc=sqrt(sourceDt(k))*Cc;
    end
    recoveries=struct();
    for kind=["nodal","root"]
        if preparedRuntimeOwner
            recoveryC=raw.prepared.(kind+"C");
            recoveryD=raw.prepared.(kind+"D");
        else
            recoveryCd=raw.(kind+"Cd");
            recoveryDd=raw.(kind+"Dd");
            [recoveryC,recoveryD]=outputDiscreteToContinuous( ...
                Ad,Bd,recoveryCd,recoveryDd,sourceDt(k));
        end
        recoveries.(kind)=struct('C',recoveryC, ...
            'D',recoveryD*inputTransform);
    end
    sourceData(k)=struct('id',ids(k),'V',raw.V, ...
        'W',raw.W,'A',Ac,'B',Bc,'C',Cc,'D',Dc, ...
        'reducedA',reducedA, ...
        'equilibrium',zeros(commonDimension,1), ...
        'affine',zeros(commonDimension,1), ...
        'Bchi',raw.Bchi, ...
        'nodalC',recoveries.nodal.C,'nodalD',recoveries.nodal.D, ...
        'rootC',recoveries.root.C,'rootD',recoveries.root.D, ...
        'Q',eye(40),'schurGauge',schurGauge,'schurBlocks',schurBlocks, ...
        'primal',zeros(commonDimension,40),'dual',zeros(commonDimension,40), ...
        'joint',joint);
    sourceData(k).primal=sourceData(k).V*schurGauge;
    sourceData(k).dual=sourceData(k).W*schurGauge;
    % The accepted package owns the centered qGamma state and affine term.
    % Lift them with the source primal basis so one-hot projection recovers
    % the protected source values exactly; sidecar equilibrium fields are
    % retained only as diagnostic provenance because they are pre-centering
    % quantities for some source contracts.
    sourceData(k).equilibrium=sourceData(k).V* ...
        sources(k).x_eq(sources(k).idx.qGam);
    sourceData(k).affine=sourceData(k).V* ...
        sources(k).parConst.affineOffset(sources(k).idx.qGam);
end

% Align trial and test pairs by the same source-only unitary gauge, then
% normalize the blended test basis without inferring it from the trial basis.
[~,reference]=max(weights); referencePrimal=sourceData(reference).primal;
for k=1:count
    if orderedSchur
        if ~isequal([sourceData(k).schurBlocks.dimension], ...
                [sourceData(reference).schurBlocks.dimension])
            error('AeroFlex:sched:OrderedSchurBlockMismatch', ...
                'Ordered real-Schur block dimensions differ across active sources.');
        end
        Q=eye(size(referencePrimal,2));
        for blockIndex=1:numel(sourceData(k).schurBlocks)
            columns=sourceData(k).schurBlocks(blockIndex).columns;
            [left,~,right]=svd( ...
                sourceData(k).primal(:,columns)'*referencePrimal(:,columns), ...
                'econ');
            Q(columns,columns)=left*right';
        end
    elseif k==reference
        Q=eye(size(referencePrimal,2));
    else
        [left,~,right]=svd(sourceData(k).primal'*referencePrimal,'econ');
        Q=left*right';
    end
    sourceData(k).Q=Q;
    sourceData(k).primal=sourceData(k).primal*Q;
    sourceData(k).dual=sourceData(k).dual*Q;
end
Vq=zeros(commonDimension,40); Wraw=zeros(commonDimension,40);
for k=1:count
    Vq=Vq+weights(k)*sourceData(k).primal;
    Wraw=Wraw+weights(k)*sourceData(k).dual;
end
if orderedSchur
    % Return the spaces to the reference source's public fixed10 gauge.
    % This preserves one-hot V/W identity without changing either subspace.
    referenceGauge=sourceData(reference).schurGauge;
    Vq=Vq*referenceGauge'; Wraw=Wraw*referenceGauge';
end
metric=Wraw'*Vq;
if rcond(metric)<1e-12
    error('AeroFlex:sched:FullCoordinateMetric', ...
        'The blended query trial/test metric is singular.');
end
Wq=Wraw/metric';
if norm(Wq'*Vq-eye(40),'fro')>1e-10
    error('AeroFlex:sched:FullCoordinateBiorthogonality', ...
        'The query pair violates Wq''*Vq=I.');
end

Ac=zeros(commonDimension); Bc=zeros(commonDimension,inputCount);
Cc=zeros(outputCount,commonDimension); Dc=zeros(outputCount,inputCount);
xeq=zeros(commonDimension,1); affine=xeq; Bchi=zeros(commonDimension,3);
nodalC=zeros(size(sourceData(1).nodalC)); nodalD=zeros(size(sourceData(1).nodalD));
rootC=zeros(size(sourceData(1).rootC)); rootD=zeros(size(sourceData(1).rootD));
structuralCount=size(chart.transforms(1).q1.T,1);
jointFull=emptyJointCoupling();
if jointCoupling
    Nm=size(sourceData(1).joint.q1q1,1);
    jointFull=struct('q1q1',zeros(Nm),'q1q2',zeros(Nm), ...
        'q1Gamma',zeros(Nm,commonDimension), ...
        'gammaQ1',zeros(commonDimension,Nm), ...
        'gammaQ2',zeros(commonDimension,Nm), ...
        'gammaGamma',zeros(commonDimension));
end
jointDiscreteFull=[];
if jointDiscrete || jointDiscreteCouplingOnly || fullCoupledCommonCoordinate || ...
        unpartitionedCoupledInterface || branchAwareCrossCoupling
    jointDiscreteFull=zeros(2*structuralCount+commonDimension);
end
for k=1:count
    item=sourceData(k); weight=weights(k);
    Ac=Ac+weight*item.A; Bc=Bc+weight*item.B;
    Cc=Cc+weight*item.C; Dc=Dc+weight*item.D;
    xeq=xeq+weight*item.equilibrium;
    affine=affine+weight*item.affine;
    Bchi=Bchi+weight*item.Bchi;
    nodalC=nodalC+weight*item.nodalC; nodalD=nodalD+weight*item.nodalD;
    rootC=rootC+weight*item.rootC; rootD=rootD+weight*item.rootD;
    if jointCoupling
        jointFull.q1q1=jointFull.q1q1+weight*sourceDt(k)*item.joint.q1q1;
        jointFull.q1q2=jointFull.q1q2+weight*sourceDt(k)*item.joint.q1q2;
        jointFull.q1Gamma=jointFull.q1Gamma+weight*sourceDt(k)*item.joint.q1Gamma;
        jointFull.gammaQ1=jointFull.gammaQ1+weight*sourceDt(k)*item.joint.gammaQ1;
        jointFull.gammaQ2=jointFull.gammaQ2+weight*sourceDt(k)*item.joint.gammaQ2;
        jointFull.gammaGamma=jointFull.gammaGamma+weight*sourceDt(k)*item.joint.gammaGamma;
    end
    if jointDiscrete || jointDiscreteCouplingOnly || fullCoupledCommonCoordinate
        jointDiscreteFull=jointDiscreteFull+weight*item.joint.normalizedSubsystem;
    end
end
if physicsNormalized
    Ac=Ac/queryDt;
    Bc=Bc/sqrt(queryDt);
    Cc=Cc/sqrt(queryDt);
end
if jointCoupling
    jointFull.q1q1=jointFull.q1q1/queryDt;
    jointFull.q1q2=jointFull.q1q2/queryDt;
    jointFull.q1Gamma=jointFull.q1Gamma/queryDt;
    jointFull.gammaQ1=jointFull.gammaQ1/queryDt;
    jointFull.gammaQ2=jointFull.gammaQ2/queryDt;
    jointFull.gammaGamma=jointFull.gammaGamma/queryDt;
end
if projectedTustinCondensation
    [Ad,Bd,Cd,Dd,projectedTustin]=projectedContinuousToDiscrete( ...
        Ac,Bc,Cc,Dc,Vq,Wq,queryDt);
    AdFull=[]; BdFull=[]; CdFull=[]; DdFull=[];
else
    [AdFull,BdFull,CdFull,DdFull]=continuousToDiscrete( ...
        Ac,Bc,Cc,Dc,queryDt);
    Ad=Wq'*AdFull*Vq; Bd=Wq'*BdFull; Cd=CdFull*Vq; Dd=DdFull;
    projectedTustin=struct();
end
[A,B,C,D]=discreteToContinuous(Ad,Bd,Cd,Dd,queryDt);
if invariantCluster
    A=zeros(40);
    for k=1:count
        alignedA=sourceData(k).Q'*sourceData(k).reducedA*sourceData(k).Q;
        A=A+weights(k)*sourceDt(k)*alignedA;
    end
    A=A/queryDt;
    [Ad,Bd,Cd,Dd]=continuousToDiscrete(A,B,C,D,queryDt);
end
if multiBandCoupled
    jointSubsystem=multiBandCoupledSubsystem(sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart);
elseif selectiveHighBandCoupled
    jointSubsystem=selectiveHighBandCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart);
elseif fullCoupledOrderedSchur
    [jointSubsystem,coordinateTransport]=orderedSchurCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart);
elseif unpartitionedCoupledInterface || branchAwareCrossCoupling
    [jointSubsystem,coordinateTransport]=unpartitionedCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart);
elseif branchPreservingStructural || fullCoupledPartitionPreserving
    oneHot=find(abs(weights-1)<=1e-14,1);
    if ~isempty(oneHot) && nnz(weights>1e-14)==1
        sourceIdx=sources(oneHot).idx;
        coupled=[sourceIdx.q1 sourceIdx.q2 sourceIdx.qGam];
        jointSubsystem=sources(oneHot).L(coupled,coupled);
    else
        jointSubsystem=branchPreservingSubsystem(sourceData,weights,Vq,Wq, ...
            sourceDt,queryDt,fullCoupledPartitionPreserving);
    end
elseif jointDiscrete || fullCoupledCommonCoordinate
    oneHot=find(abs(weights-1)<=1e-14,1);
    if ~isempty(oneHot) && nnz(weights>1e-14)==1
        sourceIdx=sources(oneHot).idx;
        coupled=[sourceIdx.q1 sourceIdx.q2 sourceIdx.qGam];
        jointSubsystem=sources(oneHot).L(coupled,coupled);
    else
        Nm=(size(jointDiscreteFull,1)-commonDimension)/2;
        if Nm~=floor(Nm)
            error('AeroFlex:sched:JointDiscreteDimensions', ...
                'The joint structural-aerodynamic subsystem has invalid dimensions.');
        end
        fullDiscrete=continuousToDiscreteA(jointDiscreteFull);
        trial=blkdiag(eye(Nm),eye(Nm),Vq);
        test=blkdiag(eye(Nm),eye(Nm),Wq);
        projectedDiscrete=test'*fullDiscrete*trial;
        jointSubsystem=discreteToContinuousA(projectedDiscrete)/queryDt;
    end
end
if jointDiscreteCouplingOnly
    oneHot=find(abs(weights-1)<=1e-14,1);
    if ~isempty(oneHot) && nnz(weights>1e-14)==1
        sourceIdx=sources(oneHot).idx;
        coupled=[sourceIdx.q1 sourceIdx.q2 sourceIdx.qGam];
        jointSubsystem=sources(oneHot).L(coupled,coupled);
    else
        Nm=(size(jointDiscreteFull,1)-commonDimension)/2;
        if Nm~=floor(Nm)
            error('AeroFlex:sched:JointDiscreteDimensions', ...
                'The joint structural-aerodynamic subsystem has invalid dimensions.');
        end
        fullDiscrete=continuousToDiscreteA(jointDiscreteFull);
        trial=blkdiag(eye(Nm),eye(Nm),Vq);
        test=blkdiag(eye(Nm),eye(Nm),Wq);
        projectedDiscrete=test'*fullDiscrete*trial;
        jointSubsystem=discreteToContinuousA(projectedDiscrete)/queryDt;
    end
end
if fullCoupledCommonCoordinate || selectiveHighBandCoupled || multiBandCoupled
    core = forceMapCoreFromCoupledSubsystem(jointSubsystem,sources,weights);
    A = core.A;
    B(:,1:(2*core.modeCount)) = [core.B0 core.B1];
    C = core.C;
    D(:,1:(2*core.modeCount)) = [core.D0 core.D1];
end
if fullCoupledOrderedSchur
    physicalSubsystem=coordinateTransport.coupledInternalToPhysical* ...
        jointSubsystem*coordinateTransport.coupledPhysicalToInternal;
    core=forceMapCoreFromCoupledSubsystem(physicalSubsystem,sources,weights);
    A=core.A;
    B(:,1:(2*core.modeCount))=[core.B0 core.B1];
    C=core.C;
    D(:,1:(2*core.modeCount))=[core.D0 core.D1];
end
if projectedTustinCondensation
    [nodalCdV,nodalDd]=projectedContinuousOutputToDiscrete( ...
        nodalC,nodalD,projectedTustin,queryDt);
    [rootCdV,rootDd]=projectedContinuousOutputToDiscrete( ...
        rootC,rootD,projectedTustin,queryDt);
    nodalCd=[]; rootCd=[];
else
    [nodalCd,nodalDd]=outputContinuousToDiscrete( ...
        Ac,Bc,nodalC,nodalD,queryDt);
    [rootCd,rootDd]=outputContinuousToDiscrete( ...
        Ac,Bc,rootC,rootD,queryDt);
    nodalCdV=nodalCd*Vq;
    rootCdV=rootCd*Vq;
end
[nodalCq,nodalDq]=outputDiscreteToContinuous( ...
    Ad,Bd,nodalCdV,nodalDd,queryDt);
[rootCq,rootDq]=outputDiscreteToContinuous( ...
    Ad,Bd,rootCdV,rootDd,queryDt);

transforms=repmat(struct('T',[],'Tinv',[],'condition',NaN),count,1);
for k=1:count
    T=Wq'*sourceData(k).V;
    transforms(k)=struct('T',T,'Tinv',pinv(T), ...
        'condition',cond(T));
end
architecture='full_coordinate_atomic_lift_interpolate_project_v1';
if physicsNormalized
    architecture='physics_normalized_blockwise_full_operator_v1';
end
if orderedSchur
    architecture='ordered_real_schur_biorthogonal_v1';
end
if invariantCluster
    architecture='invariant_cluster_biorthogonal_v1';
end
if jointCoupling
    architecture='joint_structural_aerodynamic_lift_interpolate_project_v1';
end
if jointDiscrete
    architecture='joint_structural_aerodynamic_lift_interpolate_project_v2';
end
if jointDiscreteCouplingOnly
    architecture='mixed_diagonal_offdiagonal_joint_coupling_lift_interpolate_project_v3';
end
if fullCoupledCommonCoordinate
    architecture='full_coupled_common_coordinate_realization_diagnostic_v4';
end
if branchPreservingStructural
    architecture='branch_preserving_structural_common_realization_diagnostic_v6';
end
if fullCoupledPartitionPreserving
    architecture='full_coupled_partition_preserving_branch_realization_diagnostic_v7';
end
if unpartitionedCoupledInterface
    architecture='unpartitioned_coupled_interface_realization_diagnostic_v8';
end
if branchAwareCrossCoupling
    architecture='branch_aware_cross_coupling_realization_diagnostic_v9';
end
if fullCoupledOrderedSchur
    architecture='full_coupled_ordered_schur_invariant_subspace_diagnostic_v13';
end
if selectiveHighBandCoupled
    architecture='selective_high_band_coupled_invariant_subspace_diagnostic_v14';
end
if multiBandCoupled
    architecture='multiband_coupled_invariant_subspace_diagnostic_v15';
end
candidate=struct('enabled',true, ...
    'architecture',architecture, ...
    'sourceIds',ids,'weights',weights,'referenceSource',ids(reference), ...
    'queryDt',queryDt,'sourceDt',sourceDt,'A',A,'B',B,'C',C,'D',D, ...
    'Ad',Ad,'Bd',Bd,'Cd',Cd,'Dd',Dd,'V',Vq,'W',Wq, ...
    'qGammaEquilibrium',Wq'*xeq,'qGammaAffine',Wq'*affine, ...
    'Bchi',Wq'*Bchi,'nodalC',nodalCq,'nodalD',nodalDq, ...
    'rootC',rootCq,'rootD',rootDq,'transforms',transforms, ...
    'metrics',struct('biorthogonalityFrobenius', ...
    norm(Wq'*Vq-eye(40),'fro'),'metricCondition',cond(metric), ...
    'queryReductionRank',rank(Vq),'queryReductionCondition',cond(Vq)), ...
    'registrySha256',char(fileHash(registryPath)), ...
    'queryDtRule',ternary(physicsNormalized, ...
        'accepted reference chord divided by query airspeed', ...
        'convex weighted source dt'), ...
    'physicsNormalized',physicsNormalized, ...
    'orderedSchur',orderedSchur, ...
    'invariantCluster',invariantCluster, ...
    'jointCoupling',jointCoupling, ...
    'jointDiscrete',jointDiscrete, ...
    'jointDiscreteCouplingOnly',jointDiscreteCouplingOnly, ...
    'fullCoupledCommonCoordinate',fullCoupledCommonCoordinate, ...
    'branchPreservingStructural',branchPreservingStructural, ...
    'fullCoupledPartitionPreserving',fullCoupledPartitionPreserving, ...
    'unpartitionedCoupledInterface',unpartitionedCoupledInterface, ...
    'branchAwareCrossCoupling',branchAwareCrossCoupling, ...
    'fullCoupledOrderedSchur',fullCoupledOrderedSchur, ...
    'selectiveHighBandCoupled',selectiveHighBandCoupled, ...
    'multiBandCoupled',multiBandCoupled, ...
    'referenceChordM',referenceChord, ...
    'scalingRule',ternary(physicsNormalized, ...
        'Abar=dt*A; Bbar=sqrt(dt)*B; Cbar=sqrt(dt)*C; Dbar=D', ...
        'direct continuous tuple interpolation'), ...
    'fieldAssemblyOrder','lift-interpolate-project', ...
    'immutableSourceCache',sourceCacheStats, ...
    'preparedRuntimeOwner',struct('enabled',preparedRuntimeOwner, ...
        'profile',preparedProfile,'changeId',preparedChangeId, ...
        'projectedTustinCondensation',projectedTustinCondensation, ...
        'fullDiscreteDiagnosticRetained',~projectedTustinCondensation, ...
        'sourceInvariantFields', ...
        {{'Ac','Bc','Cc','Dc','nodalC','nodalD','rootC','rootD'}}, ...
        'queryDependentFieldsRebuilt',true), ...
    'validationTruthUsedForConstruction',false,'noReducedInverseTransport',true);
candidate.full=struct('A',Ac,'B',Bc,'C',Cc,'D',Dc, ...
    'Ad',AdFull,'Bd',BdFull,'Cd',CdFull,'Dd',DdFull, ...
    'qGammaEquilibrium',xeq,'qGammaAffine',affine,'Bchi',Bchi, ...
    'nodalC',nodalC,'nodalD',nodalD,'nodalCd',nodalCd,'nodalDd',nodalDd, ...
    'rootC',rootC,'rootD',rootD,'rootCd',rootCd,'rootDd',rootDd);
if jointCoupling
    candidate.jointCoupling=struct('q1q1',jointFull.q1q1, ...
        'q1q2',jointFull.q1q2,'q1Gamma',jointFull.q1Gamma*Vq, ...
        'gammaQ1',Wq'*jointFull.gammaQ1, ...
        'gammaQ2',Wq'*jointFull.gammaQ2, ...
        'gammaGamma',Wq'*jointFull.gammaGamma*Vq, ...
        'assemblyRule', ...
        'direct source q1/q2/qGamma blocks; lift-interpolate-project');
end
if jointDiscrete
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'assemblyRule', ...
        'joint physics-normalized discrete lift-interpolate-project');
end
if jointDiscreteCouplingOnly
    candidate.jointCoupling=struct('couplingOnlySubsystem',jointSubsystem, ...
        'assemblyRule',['physics-normalized discrete structural-aerodynamic ', ...
        'off-diagonal lift-interpolate-project with retained diagonal blocks']);
end
if fullCoupledCommonCoordinate
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'forceMapCore',core, ...
        'assemblyRule',[ ...
        'single physics-normalized discrete q1/q2/qGamma tuple; ', ...
        'derive installed A/B0/B1/C/D0/D1 from that tuple before package assembly']);
end
if fullCoupledOrderedSchur
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'coordinateTransport',coordinateTransport,'forceMapCore',core, ...
        'architecture',architecture, ...
        'assemblyRule',['source-only ordered real-Schur q1/q2/qGamma ', ...
        'invariant realization; fixed-10 projection with physical recovery']);
end
if selectiveHighBandCoupled
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'forceMapCore',core,'architecture',architecture, ...
        'assemblyRule',['source-only high-band coupled invariant-subspace ', ...
        'transport with locked physical complementary block retention']);
end
if multiBandCoupled
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'forceMapCore',core,'architecture',architecture, ...
        'assemblyRule','source-only joint 6/14/43 Hz invariant-plane transport with closest complement map');
end
if branchPreservingStructural
    core=forceMapCoreFromCoupledSubsystem(jointSubsystem,sources,weights);
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'forceMapCore',core,'assemblyRule',[ ...
        'source-only 43 Hz real q1/q2 branch-plane alignment; ', ...
        'consistent qGamma coupling transport; normalized blend and fixed-10 projection']);
end
if fullCoupledPartitionPreserving
    core=forceMapCoreFromCoupledSubsystem(jointSubsystem,sources,weights);
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'forceMapCore',core,'assemblyRule',[ ...
        'source-only q1/q2/qGamma real branch-plane alignment; ', ...
        'partition-preserving normalized blend and fixed-10 projection']);
end
if unpartitionedCoupledInterface
    candidate.jointCoupling=struct('subsystem',jointSubsystem, ...
        'coordinateTransport',coordinateTransport,'assemblyRule',[ ...
        'source-only unpartitioned q1/q2/qGamma branch coordinates; ', ...
        'one common internal realization with physical public recovery']);
end
if branchAwareCrossCoupling
    physical=coordinateTransport.coupledInternalToPhysical*jointSubsystem* ...
        coordinateTransport.coupledPhysicalToInternal;
    candidate.jointCoupling=struct('couplingOnlySubsystem',physical, ...
        'assemblyRule',['source-projected high-branch common realization; ', ...
        'retain physical q1/q2/qGamma diagonal blocks']);
end
end

function subsystem=multiBandCoupledSubsystem(sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart)
count=numel(sourceData); dimension=2*numel(sources(1).idx.q1)+size(Vq,2);
oneHot=find(abs(weights-1)<=1e-14,1);
if ~isempty(oneHot) && nnz(weights>1e-14)==1
    idx=sources(oneHot).idx; subsystem=sources(oneHot).L([idx.q1 idx.q2 idx.qGam],[idx.q1 idx.q2 idx.qGam]); return
end
[~,reference]=max(weights); referenceBasis=multiBandCoupledBasis(sources(reference),sourceData(reference),chart.transforms(reference),Wq);
normalized=zeros(dimension);
for k=1:count
    modeCount=numel(sources(k).idx.q1); trial=blkdiag(eye(modeCount),eye(modeCount),Vq); test=blkdiag(eye(modeCount),eye(modeCount),Wq);
    physical=test'*sourceData(k).joint.normalizedSubsystem*trial/sourceDt(k);
    sourceBasis=multiBandCoupledBasis(sources(k),sourceData(k),chart.transforms(k),Wq);
    map=minimalSubspaceMap(sourceBasis,referenceBasis);
    normalized=normalized+weights(k)*sourceDt(k)*(map*physical*map.');
end
subsystem=discreteToContinuousA(continuousToDiscreteA(normalized))/queryDt;
end

function basis=multiBandCoupledBasis(source,sourceData,transform,Wq)
limits=[5.5 8.5;12 16;40 46]; idx=source.idx; coupled=[idx.q1 idx.q2 idx.qGam];
pc=source.parConst; pc.u_ctrl=source.u_eq(:); if isfield(pc,'gust'),pc.gust=zeros(size(pc.gust));end
pc.N_Thrust=zeros(numel(idx.q1),1); [Nq,~]=AeroFlex.sim.nonlinearJacobian(source.x_eq(:),idx,pc);
J=source.Ldyn+Nq; J(idx.q1,:)=source.Ldyn(idx.q1,:)+source.beam.Pz*Nq(idx.q1,:);
coordinate=blkdiag(transform.q1.T,transform.q2.T,Wq'*sourceData.V); [vectors,values,left]=eig(coordinate*J(coupled,coupled)/coordinate,'vector'); frequency=imag(values)/(2*pi);
basis=[]; structural=1:(2*numel(idx.q1));
for b=1:size(limits,1)
    candidates=find(frequency>0 & frequency>=limits(b,1) & frequency<=limits(b,2)); assert(~isempty(candidates),'AeroFlex:sched:MultiBandMode'); score=zeros(size(candidates));
    for k=1:numel(candidates), c=candidates(k); overlap=left(:,c)'*vectors(:,c); dual=left(:,c)/conj(overlap); e=abs(conj(dual).*vectors(:,c)); score(k)=sum(e(structural))/max(sum(e),eps); end
    [~,k]=max(score); mode=vectors(:,candidates(k)); basis=[basis real(mode) imag(mode)]; %#ok<AGROW>
end
[basis,~]=qr(basis,0); assert(size(basis,2)==6,'AeroFlex:sched:MultiBandBasis');
end

function map=minimalSubspaceMap(sourceBasis,referenceBasis)
[left,~,right]=svd(referenceBasis.'*sourceBasis,'econ'); sourceBasis=sourceBasis*right*left.';
sourceComplement=null(sourceBasis.'); referenceComplement=null(referenceBasis.');
[left,~,right]=svd(referenceComplement.'*sourceComplement,'econ'); map=referenceBasis*sourceBasis.'+referenceComplement*(left*right.')*sourceComplement.';
assert(norm(map.'*map-eye(size(map)),'fro')<=1e-10,'AeroFlex:sched:MultiBandMap');
end

function subsystem=selectiveHighBandCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart)
count=numel(sourceData); gammaCount=size(Vq,2);
fullDimension=size(sourceData(1).joint.normalizedSubsystem,1);
structuralCount=(fullDimension-size(Vq,1))/2;
if structuralCount~=floor(structuralCount) || structuralCount<=0
    error('AeroFlex:sched:SelectiveHighBandDimensions', ...
        'The coupled source subsystem has invalid dimensions.');
end
dimension=2*structuralCount+gammaCount;
oneHot=find(abs(weights-1)<=1e-14,1);
if ~isempty(oneHot) && nnz(weights>1e-14)==1
    idx=sources(oneHot).idx; coupled=[idx.q1 idx.q2 idx.qGam];
    subsystem=sources(oneHot).L(coupled,coupled);
    return
end
[~,reference]=max(weights);
referencePlane=selectiveHighBandPlane(sources(reference),sourceData(reference), ...
    chart.transforms(reference),Wq);
projector=referencePlane*referencePlane.'; complement=eye(dimension)-projector;
physicalBlend=zeros(dimension); transported=zeros(dimension);
for k=1:count
    trial=blkdiag(eye(structuralCount),eye(structuralCount),Vq);
    test=blkdiag(eye(structuralCount),eye(structuralCount),Wq);
    physical=test'*sourceData(k).joint.normalizedSubsystem*trial/sourceDt(k);
    sourcePlane=selectiveHighBandPlane(sources(k),sourceData(k), ...
        chart.transforms(k),Wq);
    map=minimalPlaneMap(sourcePlane,referencePlane);
    aligned=map*physical*map.';
    physicalBlend=physicalBlend+weights(k)*sourceDt(k)*physical;
    transported=transported+weights(k)*sourceDt(k)*( ...
        projector*aligned*projector+projector*aligned*complement+ ...
        complement*aligned*projector);
end
normalized=complement*physicalBlend*complement+transported;
discrete=continuousToDiscreteA(normalized);
subsystem=discreteToContinuousA(discrete)/queryDt;
end

function map=minimalPlaneMap(sourcePlane,referencePlane)
[left,~,right]=svd(referencePlane.'*sourcePlane,'econ');
sourcePlane=sourcePlane*right*left.';
sourceComplement=null(sourcePlane.');
referenceComplement=null(referencePlane.');
[left,~,right]=svd(referenceComplement.'*sourceComplement,'econ');
complementMap=left*right.';
map=referencePlane*sourcePlane.'+ ...
    referenceComplement*complementMap*sourceComplement.';
if norm(map.'*map-eye(size(map)),'fro')>1e-10 || ...
        norm(map*sourcePlane-referencePlane,'fro')>1e-10
    error('AeroFlex:sched:SelectiveHighBandMap', ...
        'The high-band map is not a valid closest orthogonal transport.');
end
end

function plane=selectiveHighBandPlane(source,sourceData,transform,Wq)
idx=source.idx; coupled=[idx.q1 idx.q2 idx.qGam];
pc=source.parConst; pc.u_ctrl=source.u_eq(:);
if isfield(pc,'gust'), pc.gust=zeros(size(pc.gust)); end
pc.N_Thrust=zeros(numel(idx.q1),1);
[Nq,~]=AeroFlex.sim.nonlinearJacobian(source.x_eq(:),idx,pc);
J=source.Ldyn+Nq;
J(idx.q1,:)=source.Ldyn(idx.q1,:)+source.beam.Pz*Nq(idx.q1,:);
gammaMap=Wq'*sourceData.V;
coordinate=blkdiag(transform.q1.T,transform.q2.T,gammaMap);
[vectors,values,left]=eig(coordinate*J(coupled,coupled)/coordinate,'vector');
frequency=imag(values)/(2*pi); participation=zeros(numel(values),1);
structural=[1:numel(idx.q1),numel(idx.q1)+(1:numel(idx.q2))];
for k=1:numel(values)
    overlap=left(:,k)'*vectors(:,k);
    if abs(overlap)<=eps, continue; end
    dual=left(:,k)/conj(overlap); energy=abs(conj(dual).*vectors(:,k));
    participation(k)=sum(energy(structural))/max(sum(energy),eps);
end
candidates=find(frequency>0 & frequency>=40 & frequency<=46 & participation>=0.2);
if isempty(candidates)
    error('AeroFlex:sched:SelectiveHighBandMode', ...
        'Source %s has no eligible 40-46 Hz coupled branch.',source.name);
end
[~,position]=max(participation(candidates)); mode=vectors(:,candidates(position));
plane=orth([real(mode),imag(mode)]);
if size(plane,2)~=2
    error('AeroFlex:sched:SelectiveHighBandPlane', ...
        'The selected high-band branch is not a real two-plane.');
end
end

function [subsystem,transport]=orderedSchurCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart)
count=numel(sourceData); gammaCount=size(Vq,2);
fullDimension=size(sourceData(1).joint.normalizedSubsystem,1);
structuralCount=(fullDimension-size(Vq,1))/2;
if structuralCount~=floor(structuralCount) || structuralCount<=0
    error('AeroFlex:sched:OrderedCoupledDimensions', ...
        'The coupled source subsystem has invalid dimensions.');
end
dimension=2*structuralCount+gammaCount;
oneHot=find(abs(weights-1)<=1e-14,1);
if ~isempty(oneHot) && nnz(weights>1e-14)==1
    sourceIdx=sources(oneHot).idx; coupled=[sourceIdx.q1 sourceIdx.q2 sourceIdx.qGam];
    physical=sources(oneHot).L(coupled,coupled);
    coordinate=orderedProjectedCoupledGauge(sources(oneHot),sourceData(oneHot), ...
        chart.transforms(oneHot),Wq,coupled);
    subsystem=coordinate*physical*coordinate.';
    transport=struct('coupledPhysicalToInternal',coordinate, ...
        'coupledInternalToPhysical',coordinate.', ...
        'orthogonalityError',norm(coordinate*coordinate.'-eye(dimension),'fro'), ...
        'dimension',dimension,'method','ordered_real_schur_full_coupled');
    return
end
internal=zeros(dimension); transforms=zeros(dimension,dimension,count);
for k=1:count
    trial=blkdiag(eye(structuralCount),eye(structuralCount),Vq);
    test=blkdiag(eye(structuralCount),eye(structuralCount),Wq);
    physical=test'*sourceData(k).joint.normalizedSubsystem*trial/sourceDt(k);
    coordinate=orderedProjectedCoupledGauge(sources(k),sourceData(k), ...
        chart.transforms(k),Wq,[]);
    transforms(:,:,k)=coordinate;
    internal=internal+weights(k)*sourceDt(k)*(coordinate*physical*coordinate.');
end
[left,~,right]=svd(sum(transforms.*reshape(weights,1,1,[]),3),'econ');
physicalToInternal=left*right';
discrete=continuousToDiscreteA(internal);
subsystem=discreteToContinuousA(discrete)/queryDt;
transport=struct('coupledPhysicalToInternal',physicalToInternal, ...
    'coupledInternalToPhysical',physicalToInternal.', ...
    'orthogonalityError',norm(physicalToInternal*physicalToInternal.'- ...
    eye(dimension),'fro'),'dimension',dimension,'method','ordered_real_schur_full_coupled');
end

function coordinate=orderedProjectedCoupledGauge(source,sourceData,transform,Wq,coupled)
idx=source.idx;
if isempty(coupled), coupled=[idx.q1 idx.q2 idx.qGam]; end
pc=source.parConst; pc.u_ctrl=source.u_eq(:);
if isfield(pc,'gust'), pc.gust=zeros(size(pc.gust)); end
pc.N_Thrust=zeros(numel(idx.q1),1);
[Nq,~]=AeroFlex.sim.nonlinearJacobian(source.x_eq(:),idx,pc);
J=source.Ldyn+Nq;
J(idx.q1,:)=source.Ldyn(idx.q1,:)+source.beam.Pz*Nq(idx.q1,:);
gammaMap=Wq'*sourceData.V;
S=blkdiag(transform.q1.T,transform.q2.T,gammaMap);
Jcommon=S*J(coupled,coupled)/S;
[gauge,~]=orderedRealSchurGauge(Jcommon);
coordinate=gauge.';
if norm(coordinate*coordinate.'-eye(size(coordinate)),'fro')>1e-10
    error('AeroFlex:sched:OrderedCoupledOrthogonality', ...
        'The ordered full-coupled gauge is not orthogonal.');
end
end

function [subsystem,transport]=unpartitionedCoupledSubsystem( ...
        sources,sourceData,weights,Vq,Wq,sourceDt,queryDt,chart)
count=numel(sourceData); gammaCount=size(Vq,2);
fullDimension=size(sourceData(1).joint.normalizedSubsystem,1);
structuralCount=(fullDimension-size(Vq,1))/2;
if structuralCount~=floor(structuralCount) || structuralCount<=0
    error('AeroFlex:sched:UnpartitionedDimensions', ...
        'The coupled source subsystem has invalid dimensions.');
end
dimension=2*structuralCount+gammaCount;
oneHot=find(abs(weights-1)<=1e-14,1);
if ~isempty(oneHot) && nnz(weights>1e-14)==1
    sourceIdx=sources(oneHot).idx;
    coupled=[sourceIdx.q1 sourceIdx.q2 sourceIdx.qGam];
    physical=sources(oneHot).L(coupled,coupled);
    physicalToInternal=unpartitionedProjectedBranchCoordinates( ...
        sources(oneHot),coupled);
    subsystem=physicalToInternal*physical*physicalToInternal.';
    transport=struct('coupledPhysicalToInternal',physicalToInternal, ...
        'coupledInternalToPhysical',physicalToInternal.', ...
        'orthogonalityError',norm(physicalToInternal*physicalToInternal.'- ...
        eye(dimension),'fro'),'dimension',dimension);
    return
end
internal=zeros(dimension); transforms=zeros(dimension,dimension,count);
for k=1:count
    trial=blkdiag(eye(structuralCount),eye(structuralCount),Vq);
    test=blkdiag(eye(structuralCount),eye(structuralCount),Wq);
    physical=test'*sourceData(k).joint.normalizedSubsystem*trial/sourceDt(k);
    mode=unpartitionedProjectedBranchMode(sources(k));
    transform=chart.transforms(k);
    mode=[transform.q1.T*mode(sources(k).idx.q1); ...
        transform.q2.T*mode(sources(k).idx.q2); ...
        Wq'*sourceData(k).V*mode(sources(k).idx.qGam)];
    coordinate=unpartitionedPlaneCoordinates(mode);
    transforms(:,:,k)=coordinate;
    internal=internal+weights(k)*sourceDt(k)*(coordinate*physical*coordinate.');
end
[left,~,right]=svd(sum(transforms.*reshape(weights,1,1,[]),3),'econ');
physicalToInternal=left*right';
discrete=continuousToDiscreteA(internal);
subsystem=discreteToContinuousA(discrete)/queryDt;
transport=struct('coupledPhysicalToInternal',physicalToInternal, ...
    'coupledInternalToPhysical',physicalToInternal.', ...
    'orthogonalityError',norm(physicalToInternal*physicalToInternal.'- ...
    eye(dimension),'fro'),'dimension',dimension);
end

function transform=unpartitionedBranchCoordinates(subsystem)
[vectors,values]=eig(subsystem,'vector');
frequency=imag(values)/(2*pi);
candidates=find(frequency>0 & frequency>=40 & frequency<=46);
if isempty(candidates)
    error('AeroFlex:sched:UnpartitionedBranch', ...
        'The source coupled subsystem has no positive 40-46 Hz branch.');
end
[~,position]=min(abs(frequency(candidates)-43.25));
mode=vectors(:,candidates(position));
transform=unpartitionedPlaneCoordinates(mode);
end

function transform=unpartitionedProjectedBranchCoordinates(source,coupled)
mode=unpartitionedProjectedBranchMode(source);
transform=unpartitionedPlaneCoordinates(mode(coupled));
end

function mode=unpartitionedProjectedBranchMode(source)
idx=source.idx;
pc=source.parConst;
pc.u_ctrl=source.u_eq(:);
if isfield(pc,'gust'), pc.gust=zeros(size(pc.gust)); end
pc.N_Thrust=zeros(numel(idx.q1),1);
[Nq,~]=AeroFlex.sim.nonlinearJacobian(source.x_eq(:),idx,pc);
J=source.Ldyn+Nq;
J(idx.q1,:)=source.Ldyn(idx.q1,:)+source.beam.Pz*Nq(idx.q1,:);
[vectors,values,left]=eig(J,'vector');
frequency=imag(values)/(2*pi); participation=zeros(numel(values),1);
for modeIndex=1:numel(values)
    overlap=left(:,modeIndex)'*vectors(:,modeIndex);
    if abs(overlap)<=eps, continue; end
    weight=left(:,modeIndex)/conj(overlap);
    element=abs(conj(weight).*vectors(:,modeIndex));
    participation(modeIndex)=sum(element([idx.q1 idx.q2]))/max(sum(element),eps);
end
candidates=find(frequency>0 & frequency>=30 & frequency<=55 & participation>=0.2);
if isempty(candidates)
    error('AeroFlex:sched:UnpartitionedProjectedBranch', ...
        'The projected source Jacobian has no eligible structural branch.');
end
[~,position]=max(participation(candidates));
mode=vectors(:,candidates(position));
end

function transform=unpartitionedPlaneCoordinates(mode)
anchor=find(abs(mode)==max(abs(mode)),1,'first');
phase=exp(-1i*angle(mode(anchor))); mode=mode*phase;
if real(mode(anchor))<0, mode=-mode; end
plane=orth([real(mode),imag(mode)]);
if size(plane,2)~=2
    error('AeroFlex:sched:UnpartitionedPlane', ...
        'The selected branch does not provide a real two-plane.');
end
dimension=numel(mode); basis=plane;
for index=1:dimension
    direction=zeros(dimension,1); direction(index)=1;
    direction=direction-basis*(basis.'*direction); magnitude=norm(direction);
    if magnitude>1e-10, basis(:,end+1)=direction/magnitude; end %#ok<AGROW>
    if size(basis,2)==dimension, break; end
end
if size(basis,2)~=dimension
    error('AeroFlex:sched:UnpartitionedCompletion', ...
        'Cannot complete the V8 internal coupled coordinate basis.');
end
transform=basis.';
end

function subsystem=branchPreservingSubsystem(sourceData,weights,Vq,Wq,sourceDt,queryDt,alignGamma)
count=numel(sourceData); gammaCount=size(Vq,2);
fullDimension=size(sourceData(1).joint.normalizedSubsystem,1);
structuralCount=(fullDimension-size(Vq,1))/2;
if structuralCount~=floor(structuralCount) || structuralCount<=0
    error('AeroFlex:sched:BranchPreservingDimensions', ...
        'The coupled source subsystem has invalid structural dimensions.');
end
dimension=2*structuralCount+gammaCount;
reduced=zeros(dimension,dimension,count);
for k=1:count
    trial=blkdiag(eye(structuralCount),eye(structuralCount),Vq);
    test=blkdiag(eye(structuralCount),eye(structuralCount),Wq);
    reduced(:,:,k)=test'*sourceData(k).joint.normalizedSubsystem*trial/sourceDt(k);
end
[~,reference]=max(weights);
[referencePlane,referenceGammaPlane]=branchPlanes(reduced(:,:,reference));
normalized=zeros(dimension);
for k=1:count
    [plane,gammaPlane]=branchPlanes(reduced(:,:,k));
    R=planeMap(plane,referencePlane);
    if alignGamma
        gammaMap=planeMap(gammaPlane,referenceGammaPlane);
    else
        gammaMap=eye(gammaCount);
    end
    transform=blkdiag(R,R,gammaMap);
    normalized=normalized+weights(k)*sourceDt(k)* ...
        (transform*reduced(:,:,k)*transform.');
end
discrete=continuousToDiscreteA(normalized);
subsystem=discreteToContinuousA(discrete)/queryDt;
end

function [plane,gammaPlane]=branchPlanes(subsystem)
[vectors,values]=eig(subsystem,'vector');
frequency=imag(values)/(2*pi);
candidates=find(frequency>0 & frequency>=40 & frequency<=46);
if isempty(candidates)
    error('AeroFlex:sched:BranchPreservingMode', ...
        'The source coupled subsystem has no positive 40-46 Hz branch.');
end
[~,position]=min(abs(frequency(candidates)-43.25));
structuralCount=(size(subsystem,1)-40)/2;
q1=vectors(1:structuralCount,candidates(position));
plane=orth([real(q1),imag(q1)]);
if size(plane,2)~=2
    error('AeroFlex:sched:BranchPreservingPlane', ...
        'The selected structural branch does not provide a real two-plane.');
end
gamma=vectors(2*structuralCount+(1:40),candidates(position));
gammaPlane=orth([real(gamma),imag(gamma)]);
if size(gammaPlane,2)~=2
    error('AeroFlex:sched:BranchPreservingGammaPlane', ...
        'The selected branch does not provide a real qGamma two-plane.');
end
end

function map=planeMap(sourcePlane,referencePlane)
sourceComplement=null(sourcePlane.');
referenceComplement=null(referencePlane.');
map=[referencePlane,referenceComplement]*[sourcePlane,sourceComplement].';
if norm(map.'*map-eye(size(map)),'fro')>1e-10
    error('AeroFlex:sched:BranchPreservingMap', ...
        'The structural branch alignment is not orthogonal.');
end
end

function value=ternary(condition,left,right)
if condition,value=left;else,value=right;end
end

function value=emptySource()
value=struct('id',"",'V',[],'W',[],'A',[],'B',[],'C',[],'D',[], ...
    'reducedA',[], ...
    'equilibrium',[],'affine',[],'Bchi',[],'nodalC',[],'nodalD',[], ...
    'rootC',[],'rootD',[],'Q',[],'schurGauge',[],'schurBlocks',[], ...
    'primal',[],'dual',[],'joint',struct());
end
function value=emptyJointCoupling()
value=struct('q1q1',[],'q1q2',[],'q1Gamma',[], ...
    'gammaQ1',[],'gammaQ2',[],'gammaGamma',[]);
end
function core=forceMapCoreFromCoupledSubsystem(subsystem,sources,weights)
modeCount=numel(sources(1).idx.q1);
gammaCount=size(subsystem,1)-2*modeCount;
if gammaCount<=0 || ~isequal(size(subsystem),[2*modeCount+gammaCount, ...
        2*modeCount+gammaCount]) || any(~isfinite(subsystem),'all')
    error('AeroFlex:sched:FullCoupledCoreDimensions', ...
        'The full coupled subsystem is invalid.');
end
scaleA=0; tInf=0;
for k=1:numel(sources)
    scaleA=scaleA+weights(k)*double(sources(k).parConst.scaleA);
    tInf=tInf+weights(k)*double(sources(k).parConst.t_inf);
end
if ~isfinite(scaleA) || ~isfinite(tInf) || scaleA<=0 || tInf<=0
    error('AeroFlex:sched:FullCoupledCoreScale', ...
        'The query aerodynamic scale or convective time is invalid.');
end
q1=1:modeCount; q2=modeCount+(1:modeCount);
qGamma=2*modeCount+(1:gammaCount);
omega=-subsystem(q2,q1);
if rcond(omega)<1e-12
    error('AeroFlex:sched:FullCoupledCoreOmega', ...
        'The transported structural Omega block is ill-conditioned.');
end
core=struct('A',subsystem(qGamma,qGamma)*tInf, ...
    'B0',-subsystem(qGamma,q2)*omega*tInf, ...
    'B1',subsystem(qGamma,q1), ...
    'C',subsystem(q1,qGamma)/scaleA, ...
    'D0',(omega-subsystem(q1,q2))*omega/scaleA, ...
    'D1',subsystem(q1,q1)/(scaleA*tInf), ...
    'Omega',omega,'scaleA',scaleA,'tInf',tInf, ...
    'modeCount',modeCount,'gammaCount',gammaCount);
reconstructed=subsystem;
reconstructed(q1,q1)=scaleA*tInf*core.D1;
reconstructed(q1,q2)=omega-scaleA*core.D0/omega;
reconstructed(q1,qGamma)=scaleA*core.C;
reconstructed(qGamma,q1)=core.B1;
reconstructed(qGamma,q2)=-core.B0/omega/tInf;
reconstructed(qGamma,qGamma)=core.A/tInf;
core.coupledOperatorRelativeResidual=norm(reconstructed-subsystem,'fro')/ ...
    max(1,norm(subsystem,'fro'));
if core.coupledOperatorRelativeResidual>1e-10
    error('AeroFlex:sched:FullCoupledCoreClosure', ...
        'The force-map core does not reconstruct the transported coupled operator.');
end
end
function [gauge,blocks]=orderedRealSchurGauge(A)
relativeImaginary=norm(imag(A),'fro')/max(norm(A,'fro'),eps);
if relativeImaginary>1e-12
    error('AeroFlex:sched:OrderedSchurComplexSource', ...
        'Source realization has relative imaginary part %.3e.',relativeImaginary);
end
[rawGauge,T]=schur(real(A),'real');
n=size(T,1); start=1; raw=struct('columns',{},'dimension',{}, ...
    'frequencyHz',{},'growth',{},'originalIndex',{});
blockIndex=0; threshold=100*eps(max(1,norm(T,'fro')));
while start<=n
    blockIndex=blockIndex+1; dimension=1;
    if start<n && abs(T(start+1,start))>threshold, dimension=2; end
    columns=start:(start+dimension-1); values=eig(T(columns,columns));
    raw(blockIndex)=struct('columns',columns,'dimension',dimension, ...
        'frequencyHz',max(abs(imag(values)))/(2*pi), ...
        'growth',max(real(values)),'originalIndex',blockIndex);
    start=start+dimension;
end
keys=[[raw.frequencyHz].',[raw.growth].',[raw.originalIndex].'];
[~,order]=sortrows(keys,[1 2 3]); permutation=[]; blocks=raw(order);
for k=1:numel(order)
    permutation=[permutation,raw(order(k)).columns]; %#ok<AGROW>
    first=numel(permutation)-raw(order(k)).dimension+1;
    blocks(k).columns=first:numel(permutation);
end
gauge=rawGauge(:,permutation);
if norm(gauge'*gauge-eye(n),'fro')>1e-10
    error('AeroFlex:sched:OrderedSchurOrthogonality', ...
        'Ordered real-Schur gauge is not orthogonal.');
end
end
function id=sourceIdAtCoordinates(registry,mu)
coordinates=reshape([registry.sources.coordinates],2,[]).';
match=find(max(abs(coordinates-double(mu(:).')),[],2)<=1e-12);
if numel(match)~=1
    error('AeroFlex:sched:FullCoordinateSourceIdentity', ...
        'Coordinates [%g %g] resolve to %d sources.',mu(1),mu(2),numel(match));
end
id=string(registry.sources(match).sourceId);
end
function metadata=verifiedImmutableMetadata(artifact,sidecar,expectedHash,id,useCache)
% Verify immutable source metadata once per unchanged file pair in a session.
persistent cache
if isempty(cache), cache=struct('key',{},'artifactStamp',{}, ...
        'sidecarStamp',{},'dt',{}); end
artifact=string(artifact); sidecar=string(sidecar); expectedHash=string(expectedHash);
if ~isfile(artifact) || ~isfile(sidecar)
    error('AeroFlex:sched:FullCoordinateSourceHash', ...
        'Source %s fixed-node27 contract or sidecar is absent.',id);
end
artifactStamp=fileStamp(artifact); sidecarStamp=fileStamp(sidecar);
key=canonicalPath(artifact)+"|"+expectedHash+"|"+canonicalPath(sidecar)+"|"+id;
if useCache
    hit=find(string({cache.key})==key,1);
    if ~isempty(hit) && isequal(cache(hit).artifactStamp,artifactStamp) && ...
            isequal(cache(hit).sidecarStamp,sidecarStamp)
        metadata=struct('dt',cache(hit).dt, ...
            'artifactStamp',artifactStamp,'sidecarStamp',sidecarStamp); return
    end
end
if fileHash(artifact)~=expectedHash
    error('AeroFlex:sched:FullCoordinateSourceHash', ...
        'Source %s fixed-node27 contract is stale.',id);
end
if string(h5readatt(sidecar,'/','source_id'))~=id || ...
        string(h5readatt(sidecar,'/','contract_family'))~= ...
        "FULL_COORDINATE_FIELD_SIDECAR"
    error('AeroFlex:sched:FullCoordinateSidecar', ...
        'Source %s full-coordinate sidecar is invalid.',id);
end
metadata=struct('dt',double(h5readatt(artifact, ...
    '/fixed_node27_premodal_common_discrete','dt')), ...
    'artifactStamp',artifactStamp,'sidecarStamp',sidecarStamp);
if useCache
    cache(end+1)=struct('key',key,'artifactStamp',artifactStamp, ...
        'sidecarStamp',sidecarStamp,'dt',metadata.dt);
end
end
function stamp=fileStamp(path)
info=dir(path); stamp=struct('bytes',double(info.bytes), ...
    'datenum',double(info.datenum));
end
function path=canonicalPath(path)
path=string(char(java.io.File(char(path)).getCanonicalPath()));
end
function [raw,info]=verifiedImmutableSourceData( ...
        artifact,sidecar,expectedHash,id,useCache,prepareRuntime)
% Read immutable source arrays once per unchanged source pair in a session.
persistent cache
startTime=tic;
info=struct('hit',false,'preparedHit',false,'preparedMiss',false, ...
    'rawH5Reads',0,'preparationSeconds',0,'wallSeconds',0);
if isempty(cache), cache=struct('key',{},'artifactStamp',{}, ...
        'sidecarStamp',{},'raw',{}); end
metadata=verifiedImmutableMetadata(artifact,sidecar,expectedHash,id,useCache);
key=canonicalPath(artifact)+"|"+string(expectedHash)+"|"+ ...
    canonicalPath(sidecar)+"|"+string(id);
if useCache
    hit=find(string({cache.key})==key,1);
    if ~isempty(hit) && isequal(cache(hit).artifactStamp,metadata.artifactStamp) && ...
            isequal(cache(hit).sidecarStamp,metadata.sidecarStamp)
        raw=cache(hit).raw; info.hit=true;
        if prepareRuntime
            if isfield(raw,'prepared')
                info.preparedHit=true;
            else
                preparationTimer=tic;
                raw.prepared=prepareImmutableSourceData(raw,metadata.dt);
                info.preparationSeconds=toc(preparationTimer);
                info.preparedMiss=true;
                cache(hit).raw=raw;
            end
        end
        info.wallSeconds=toc(startTime); return
    end
end
raw=struct( ...
    'Ad',readH5(artifact,'/fixed_node27_premodal_common_discrete/A'), ...
    'Bd',readH5(artifact,'/fixed_node27_premodal_common_discrete/B'), ...
    'Cd',readH5(artifact,'/fixed_node27_premodal_common_discrete/C'), ...
    'Dd',readH5(artifact,'/fixed_node27_premodal_common_discrete/D'), ...
    'reducedA',readH5(artifact,'/continuous/A'), ...
    'V',readH5(artifact,'/projectors/V'), ...
    'W',readH5(artifact,'/projectors/W'), ...
    'Bchi',readH5(sidecar,'/ports/Bchi'), ...
    'nodalCd',readH5(sidecar,'/recovery/nodal_C_discrete'), ...
    'nodalDd',readH5(sidecar,'/recovery/nodal_D_discrete'), ...
    'rootCd',readH5(sidecar,'/recovery/root_C_discrete'), ...
    'rootDd',readH5(sidecar,'/recovery/root_D_discrete'));
info.rawH5Reads=12;
if prepareRuntime
    preparationTimer=tic;
    raw.prepared=prepareImmutableSourceData(raw,metadata.dt);
    info.preparationSeconds=toc(preparationTimer);
    info.preparedMiss=true;
end
if useCache
    maximumEntries=14;
    if numel(cache)>=maximumEntries, cache(1)=[]; end
    cache(end+1)=struct('key',key,'artifactStamp',metadata.artifactStamp, ...
        'sidecarStamp',metadata.sidecarStamp,'raw',raw);
end
info.wallSeconds=toc(startTime);
end
function prepared=prepareImmutableSourceData(raw,dt)
% Convert only immutable source arrays; query charts remain online-owned.
[Ac,Bc,Cc,Dc]=discreteToContinuous( ...
    raw.Ad,raw.Bd,raw.Cd,raw.Dd,dt);
[nodalC,nodalD]=outputDiscreteToContinuous( ...
    raw.Ad,raw.Bd,raw.nodalCd,raw.nodalDd,dt);
[rootC,rootD]=outputDiscreteToContinuous( ...
    raw.Ad,raw.Bd,raw.rootCd,raw.rootDd,dt);
prepared=struct('Ac',Ac,'Bc',Bc,'Cc',Cc,'Dc',Dc, ...
    'nodalC',nodalC,'nodalD',nodalD, ...
    'rootC',rootC,'rootD',rootD);
end
function value=readH5(path,dataset)
value=h5read(path,dataset);
if isstruct(value) && isfield(value,'r')
    value=double(value.r)+1i*double(value.i);
else
    value=double(value);
end
if ismatrix(value) && ~isvector(value),value=value.';end
if any(~isfinite(value),'all')
    error('AeroFlex:sched:FullCoordinateFinite', ...
        'Dataset %s in %s is nonfinite.',dataset,path);
end
end
function [Ac,Bc,Cc,Dc]=discreteToContinuous(A,B,C,D,dt)
omega=2/dt; M=A+eye(size(A));
if rcond(M)<1e-14,error('AeroFlex:sched:FullCoordinateTustin','Singular inverse Tustin map.');end
Ac=omega*(A-eye(size(A)))/M; Bc=sqrt(2*omega)*(M\B);
Cc=sqrt(2*omega)*(C/M); Dc=D-C*(M\B);
end
function [A,B,C,D]=continuousToDiscrete(Ac,Bc,Cc,Dc,dt)
omega=2/dt; M=omega*eye(size(Ac))-Ac;
if rcond(M)<1e-14,error('AeroFlex:sched:FullCoordinateTustin','Singular forward Tustin map.');end
A=(omega*eye(size(Ac))+Ac)/M; B=sqrt(2*omega)*(M\Bc);
C=sqrt(2*omega)*(Cc/M); D=Dc+Cc*(M\Bc);
end
function [A,B,C,D,prepared]=projectedContinuousToDiscrete( ...
        Ac,Bc,Cc,Dc,V,W,dt)
omega=2/dt; identity=eye(size(Ac)); M=omega*identity-Ac;
if rcond(M)<1e-14
    error('AeroFlex:sched:FullCoordinateTustin', ...
        'Singular projected forward Tustin map.');
end
factor=decomposition(M,'lu');
XV=factor\V; XB=factor\Bc; scale=sqrt(2*omega);
A=W'*((omega*identity+Ac)*XV);
B=scale*(W'*XB);
C=scale*(Cc*XV);
D=Dc+Cc*XB;
prepared=struct('XV',XV,'XB',XB);
end
function [C,D]=projectedContinuousOutputToDiscrete(Cc,Dc,prepared,dt)
scale=sqrt(4/dt);
C=scale*(Cc*prepared.XV);
D=Dc+Cc*prepared.XB;
end
function A=continuousToDiscreteA(Ac)
M=eye(size(Ac))-0.5*Ac;
if rcond(M)<1e-14
    error('AeroFlex:sched:JointDiscreteTustin', ...
        'The normalized joint forward Tustin map is singular.');
end
A=(eye(size(Ac))+0.5*Ac)/M;
end
function Ac=discreteToContinuousA(A)
M=A+eye(size(A));
if rcond(M)<1e-14
    error('AeroFlex:sched:JointDiscreteTustin', ...
        'The normalized joint inverse Tustin map is singular.');
end
Ac=2*(A-eye(size(A)))/M;
end
function [Cc,Dc]=outputDiscreteToContinuous(A,B,C,D,dt)
omega=2/dt; M=A+eye(size(A));
Cc=sqrt(2*omega)*(C/M); Dc=D-C*(M\B);
end
function [C,D]=outputContinuousToDiscrete(Ac,Bc,Cc,Dc,dt)
omega=2/dt; M=omega*eye(size(Ac))-Ac;
C=sqrt(2*omega)*(Cc/M); D=Dc+Cc*(M\Bc);
end
function path=resolveRepositoryPath(path)
if isfile(path),return,end
repositoryRoot=fileparts(fileparts(fileparts(fileparts( ...
    mfilename('fullpath')))));
candidate=fullfile(repositoryRoot,path);
if isfile(candidate),path=string(candidate);end
end
function digest=fileHash(path)
f=fopen(path,'rb');c=onCleanup(@()fclose(f));data=fread(f,Inf,'*uint8');
engine=java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(data,'int8'));bytes=typecast(engine.digest(),'uint8');
digest=lower(string(reshape(dec2hex(bytes,2).',1,[])));clear c
end
