function sched = evalLibrary(ROMlibIn, mu, cfgLibrary)
%EVALLIBRARY Evaluate/interpolate the scheduled ROM at mu = [U, alpha_deg].
%
% Returned fields are ready to apply to ROMIntegrator/SimRunner/PlantRunTime:
%   sched.L, sched.parConst, sched.idx, sched.beam, sched.base, sched.aero
%
% Stage 1 performs direct componentwise interpolation of all fields consumed
% by nonlinear_terms and the linear IMEX solve.  If a compatible-coordinate
% library has been precomputed, the transformed matrices/fields are used.

if nargin < 3, cfgLibrary = struct(); end
ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
requireCompatible = true;
if isfield(cfgLibrary,'requireCompatible')
    requireCompatible = logical(cfgLibrary.requireCompatible);
end

if requireCompatible
    if ~isfield(ROMlib,'compatibleCoordinates') || ~ROMlib.compatibleCoordinates
        error('evalLibrary:NotCompatible', ...
            'This library has not been finalized into common reduced coordinates.');
    end
end

[w, ids, info] = AeroFlex.sched.interpWeights(ROMlib, mu, cfgLibrary);
P = ROMlib.points(ids);
registryResolution = struct();
if isfield(cfgLibrary,'sourceContractRegistryPath') && ...
        strlength(string(cfgLibrary.sourceContractRegistryPath))>0
    sourceIds = strings(1,numel(P));
    for k = 1:numel(P)
        if isfield(P(k),'sourceContractId') && ...
                strlength(string(P(k).sourceContractId))>0
            sourceIds(k) = string(P(k).sourceContractId);
        else
            sourceIds(k) = string(P(k).name);
        end
    end
    registryResolution = AeroFlex.sched.resolveSourceContractRegistry( ...
        cfgLibrary.sourceContractRegistryPath,sourceIds);
end
if numel(P)>1 && isfield(cfgLibrary, ...
        'activeCellRecoveryCertificateRegistryPath') && ...
        strlength(string( ...
        cfgLibrary.activeCellRecoveryCertificateRegistryPath))>0
    [ambientContract,activeCertificate] = ...
        AeroFlex.sched.resolveActiveCellPhysicalRecoveryCertificate( ...
        cfgLibrary.activeCellRecoveryCertificateRegistryPath, ...
        cfgLibrary.caseView,string({P.name}),w,registryResolution);
    ROMlib.ambientContract = ambientContract;
    ROMlib.physicalRecoveryCertificates = activeCertificate;
end

% Default-disabled Case-B path: an explicitly approved, hash-bound stencil
% may install the existing full-coordinate query package.  This branch is
% intentionally before the componentwise active-cell certificate check:
% the complete dynamic/recovery tuple replaces that incompatible path rather
% than bypassing its rejection for an arbitrary cell.
[fullCoordinateEnabled,fullCoordinateOptions] = ...
    resolveFullCoordinateRuntimeCandidate(cfgLibrary,P,mu);
if fullCoordinateEnabled
    [sched,~] = AeroFlex.sched.buildPhysicalChartPackage(P,w,mu, ...
        "fixed_physical_svd",struct(),fullCoordinateOptions);
    sched.mu = double(mu(:).');
    sched.weights = w;
    sched.pointIds = ids;
    sched.pointMu = reshape([P.mu],numel(sched.mu),[]).';
    sched.info = info;
    sched.method = ROMlib.method;
    sched.muNames = ROMlib.muNames;
    sched.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    sched.fullCoordinateRuntimeCandidate = struct( ...
        'enabled',true,'sourceIds',string({P.name}), ...
        'query',sched.mu,'registryPath',string(fullCoordinateOptions.registryPath), ...
        'fieldRoot',string(fullCoordinateOptions.fieldRoot));
    sourceDt = arrayfun(@(p) double(p.parConst.dt),P(:));
    sched.time = struct('convention','physical_seconds', ...
        'policy','full_coordinate_query_package', ...
        'sourceDt',sourceDt,'runtimeDt',sched.parConst.dt);
    sched.eigL = eig(full(sched.L));
    AeroFlex.sched.assertPackageConsistency(sched);
    return
end

% A legacy compatibleCoordinates or library-wide Boolean proves only
% reduced-state alignment. Multi-source evaluation requires a hash-bound
% certificate for the exact active cell.
[recoveryCompatible,recoveryCertificate] = ...
    verifyPhysicalRecoveryCertificate(ROMlib,P,ids,info,registryResolution);
if numel(P) > 1 && ~recoveryCompatible
    sourceText = strjoin(string({P.name}),', ');
    error('evalLibrary:IncompatibleBeamMaps', ...
        ['Active ROM cell rejected: no verified common physical recovery.\n', ...
         'query=[%.12g %.12g], cell=%s, sources={%s}.\n', ...
         'certificate=%s, first incompatible field=%s.\n', ...
         'Exact-node evaluation remains available. Audit artifact: %s'], ...
        mu(1),mu(2),recoveryCertificate.cellId,sourceText, ...
        recoveryCertificate.state,recoveryCertificate.firstField, ...
        recoveryCertificate.auditPath);
end
% Store the physical coordinates of the active interpolation vertices.
% This is useful for checking that the scheduler is using the expected
% U-alpha simplex, not merely the expected point indices.
schedPointMu = NaN(numel(P), numel(mu(:).'));
for k = 1:numel(P)
    if isfield(P(k), 'mu') && ~isempty(P(k).mu)
        schedPointMu(k,1:numel(P(k).mu(:).')) = P(k).mu(:).';
    end
end


sched = struct();
sched.mu = double(mu(:).');
sched.weights = w;
sched.pointIds = ids;
sched.pointMu = schedPointMu;
sched.info = info;
if ~isempty(fieldnames(registryResolution))
    sched.sourceContractRegistry = registryResolution;
end
sched.method = ROMlib.method;
sched.muNames = ROMlib.muNames;
sched.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));

% Linear operator.
vals = cell(numel(P),1);
for k = 1:numel(P), vals{k} = P(k).L; end
sched.L = AeroFlex.sched.lincombNumeric(vals,w);

% Nonlinear bundle: Gamma tensors, force maps, steady loads, scalings.
sched.parConst = AeroFlex.sched.interpParConst(P,w);

% Source sample times are provenance from the discrete-to-continuous
% conversion. They are never interpolated. Use an explicit configured
% physical step when supplied; otherwise use the smallest active-vertex step.
sourceDt = arrayfun(@(p) double(p.parConst.dt), P(:));
if any(~isfinite(sourceDt) | sourceDt <= 0)
    error('evalLibrary:InvalidSourceDt', ...
        'Every active source point must provide a positive finite dt.');
end
if isfield(cfgLibrary,'runtimeDt') && ~isempty(cfgLibrary.runtimeDt)
    runtimeDt = double(cfgLibrary.runtimeDt);
    timePolicy = 'configured_physical_dt';
else
    runtimeDt = min(sourceDt);
    timePolicy = 'minimum_active_source_dt';
end
if ~isscalar(runtimeDt) || ~isfinite(runtimeDt) || runtimeDt <= 0
    error('evalLibrary:InvalidRuntimeDt', ...
        'The scheduled runtime dt must be a positive finite scalar.');
end
sched.parConst.dt = runtimeDt;
sched.time = struct('convention','physical_seconds', ...
    'policy',timePolicy,'sourceDt',sourceDt,'runtimeDt',runtimeDt);

% The constraint projectors share the common q1 coordinates. Physical root
% observation maps are parameter-dependent outputs and are interpolated in
% those coordinates; equality with the reference map is not required.
sched.idx = P(1).idx;
sched.beam = P(1).beam;
sched.aero = P(1).aero;
sched.base = P(1).base;

% Beam maps are copied, not interpolated. Check that the active vertices agree.
sched.validation.beamMap = struct('dPz',0,'dPr',0,'dPhi',0);


for k = 2:numel(P)
    if isfield(P(k).beam,'Pz') && isfield(P(1).beam,'Pz')
        assertSameMapSize(P(k).beam.Pz,P(1).beam.Pz,'beam.Pz',P(k),P(1), ...
            sched.mu,sched.info);
        sched.validation.beamMap.dPz = max(sched.validation.beamMap.dPz, ...
            norm(P(k).beam.Pz - P(1).beam.Pz,'fro'));
    end

    if isfield(P(k).beam,'Pr') && isfield(P(1).beam,'Pr')
        assertSameMapSize(P(k).beam.Pr,P(1).beam.Pr,'beam.Pr',P(k),P(1), ...
            sched.mu,sched.info);
        sched.validation.beamMap.dPr = max(sched.validation.beamMap.dPr, ...
            norm(P(k).beam.Pr - P(1).beam.Pr,'fro'));
    end

    if isfield(P(k).beam,'red') && isfield(P(k).beam.red,'phi1_sA') && ...
            isfield(P(1).beam,'red') && isfield(P(1).beam.red,'phi1_sA')
        sched.validation.beamMap.dPhi = max(sched.validation.beamMap.dPhi, ...
            norm(P(k).beam.red.phi1_sA - P(1).beam.red.phi1_sA,'fro'));
    end
end

if isfield(cfgLibrary,'debug') && cfgLibrary.debug
    fprintf('[evalLibrary] beam map diffs over active simplex: dPz=%.3e dPr=%.3e dPhi=%.3e\n', ...
        sched.validation.beamMap.dPz, ...
        sched.validation.beamMap.dPr, ...
        sched.validation.beamMap.dPhi);

    fprintf('[evalLibrary] active vertices:\n');
    for k = 1:numel(ids)
        fprintf('  id=%4d  w=%+.6f  mu=[%+.6f %+.6f]\n', ...
            ids(k), w(k), sched.pointMu(k,1), sched.pointMu(k,2));
    end
end

tolBeamMap = 1e-8;
if isfield(cfgLibrary,'beamMapTol') && ~isempty(cfgLibrary.beamMapTol)
    tolBeamMap = cfgLibrary.beamMapTol;
end

allowBadBeamMaps = false;
if isfield(cfgLibrary,'allowIncompatibleBeamMaps')
    allowBadBeamMaps = logical(cfgLibrary.allowIncompatibleBeamMaps);
end

badBeamMap = sched.validation.beamMap.dPz  > tolBeamMap || ...
             sched.validation.beamMap.dPr  > tolBeamMap;

if badBeamMap && ~allowBadBeamMaps
    error('evalLibrary:IncompatibleBeamMaps', ...
          ['Active ROM simplex is not in a common beam basis.\n', ...
           'dPz=%.3e, dPr=%.3e, dPhi=%.3e. ', ...
           'Direct interpolation is not valid.'], ...
           sched.validation.beamMap.dPz, ...
           sched.validation.beamMap.dPr, ...
           sched.validation.beamMap.dPhi);
end

if isfield(sched.beam,'red')
    recoveryFields = {'phi1_sA','phi2_sA','phi_sA'};
    for j = 1:numel(recoveryFields)
        field = recoveryFields{j};
        if isfield(P(1).beam.red,field)
            vals = arrayfun(@(point){point.beam.red.(field)},P);
            sched.beam.red.(field) = ...
                AeroFlex.sched.lincombNumeric(vals,w);
        end
    end
end

% Baseline fields that vary with alpha.
if isfield(P(1).base,'Gamma_xi')
    vals = arrayfun(@(p){p.base.Gamma_xi},P);
    sched.base.Gamma_xi = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).base,'Gamma_g')
    vals = arrayfun(@(p){p.base.Gamma_g},P);
    sched.base.Gamma_g = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).base,'xi_bar')
    vals = arrayfun(@(p){p.base.xi_bar},P);
    sched.base.xi_bar = AeroFlex.sched.lincombNumeric(vals,w);
end

% Trim/equilibrium offsets.
vals = cell(numel(P),1);
for k = 1:numel(P)
    vals{k} = P(k).x_eq;
end
sched.x_eq = AeroFlex.sched.lincombNumeric(vals,w);

vals = cell(numel(P),1);
for k = 1:numel(P)
    vals{k} = P(k).u_eq;
end
sched.u_eq = AeroFlex.sched.lincombNumeric(vals,w);

if isscalar(P)
    sched = installExactNodeRecovery(sched,P);
elseif recoveryCompatible
    sched = installActiveCellRecovery(sched,P,w);
end

sched.trim = P(1).trim;
if isfield(P(1).trim,'alphaDeg')
    vals = arrayfun(@(p){p.trim.alphaDeg},P);
    sched.trim.alphaDeg = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'deltaDeg')
    vals = arrayfun(@(p){p.trim.deltaDeg},P);
    sched.trim.deltaDeg = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'deltaElev')
    vals = arrayfun(@(p){p.trim.deltaElev},P);
    sched.trim.deltaElev = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'thrust')
    vals = arrayfun(@(p){p.trim.thrust},P);
    sched.trim.thrust = AeroFlex.sched.lincombNumeric(vals,w);
end
sched.trim.states = sched.x_eq;
sched.compatibleCoordinates = isfield(ROMlib,'compatibleCoordinates') && ...
                              logical(ROMlib.compatibleCoordinates);

% Keep scheduled forceMap fields consistent with the interpolated parConst.
% This is mainly for runtime diagnostics and code paths that inspect aero.forceMap.
if ~isfield(sched,'aero') || ~isstruct(sched.aero)
    sched.aero = struct();
end
if ~isfield(sched.aero,'forceMap') || ~isstruct(sched.aero.forceMap)
    sched.aero.forceMap = struct();
end

sched.aero.forceMap.Bw       = sched.parConst.Bw;
sched.aero.forceMap.Dw       = sched.parConst.Dw;
sched.aero.forceMap.B_delta  = sched.parConst.Bdel;
sched.aero.forceMap.D_delta  = sched.parConst.Ddel;
sched.aero.forceMap.B_ddelta = sched.parConst.Bddel;
sched.aero.forceMap.D_ddelta = sched.parConst.Dddel;

% Diagnostics.
sched.eigL = eig(full(sched.L));
AeroFlex.sched.assertPackageConsistency(sched);
end

function [enabled,atomicOptions] = resolveFullCoordinateRuntimeCandidate( ...
        cfgLibrary,points,query)
% Resolve only an explicit, source-set-bound full-coordinate runtime stencil.
% Absence of this audit/promoted selector leaves the standard fail-closed
% multi-source path unchanged.
enabled = false;
atomicOptions = struct();
if ~isfield(cfgLibrary,'fullCoordinateRuntimeCandidate') || ...
        ~isstruct(cfgLibrary.fullCoordinateRuntimeCandidate)
    return
end
candidate = cfgLibrary.fullCoordinateRuntimeCandidate;
if ~isfield(candidate,'enabled') || ~isscalar(candidate.enabled) || ...
        ~logical(candidate.enabled)
    return
end
required = {'registryPath','fieldRoot'};
missing = required(~isfield(candidate,required));
if ~isempty(missing)
    error('evalLibrary:FullCoordinateCandidateConfig', ...
        'Full-coordinate runtime candidate is missing: %s.', ...
        strjoin(missing,', '));
end
tolerance = 1e-12;
if isfield(candidate,'queryTolerance') && ~isempty(candidate.queryTolerance)
    tolerance = double(candidate.queryTolerance);
end
if ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance < 0
    error('evalLibrary:FullCoordinateCandidateConfig', ...
        'Candidate queryTolerance must be a finite nonnegative scalar.');
end
sourceIds = sort(string({points.name}).');
matchedStencil = false;
if isfield(candidate,'stencils') && ~isempty(candidate.stencils)
    stencils = candidate.stencils;
    if ~isstruct(stencils) || ~all(isfield(stencils,{'query','sourceIds'}))
        error('evalLibrary:FullCoordinateCandidateConfig', ...
            'Candidate stencils must provide query and sourceIds fields.');
    end
    match = false(1,numel(stencils));
    for k = 1:numel(stencils)
        stencilQuery = double(stencils(k).query(:).');
        stencilSources = sort(string(stencils(k).sourceIds(:)));
        if isequal(size(stencilQuery),[1,2]) && ...
                all(isfinite(stencilQuery)) && ...
                norm(query(:).'-stencilQuery,inf) <= tolerance && ...
                isequal(stencilSources,sourceIds)
            match(k) = true;
        end
    end
    if nnz(match) > 1
        error('evalLibrary:FullCoordinateCandidateConfig', ...
            'Candidate exact stencils are not unique at query=[%.12g %.12g].', ...
            query(1),query(2));
    end
    matchedStencil = nnz(match) == 1;
end

matchedCenterline = false;
if ~matchedStencil && isfield(candidate,'centerline') && ...
        ~isempty(candidate.centerline)
    centerline = candidate.centerline;
    centerlineRequired = {'speedRangeMps','alphaOffsetDeg','sourceSetKeys'};
    if ~isstruct(centerline) || ~isscalar(centerline)
        error('evalLibrary:FullCoordinateCandidateConfig', ...
            'Candidate centerline must be a scalar structure.');
    end
    centerlineMissing = centerlineRequired(~isfield(centerline,centerlineRequired));
    if ~isempty(centerlineMissing)
        error('evalLibrary:FullCoordinateCandidateConfig', ...
            'Candidate centerline must provide speedRangeMps, alphaOffsetDeg, and sourceSetKeys.');
    end
    speedRange = double(centerline.speedRangeMps(:).');
    alphaOffset = double(centerline.alphaOffsetDeg);
    centerlineTolerance = tolerance;
    if isfield(centerline,'tolerance') && ~isempty(centerline.tolerance)
        centerlineTolerance = double(centerline.tolerance);
    end
    sourceSetKeys = string(centerline.sourceSetKeys(:));
    if ~isequal(size(speedRange),[1,2]) || any(~isfinite(speedRange)) || ...
            speedRange(1) > speedRange(2) || ~isscalar(alphaOffset) || ...
            ~isfinite(alphaOffset) || ~isscalar(centerlineTolerance) || ...
            ~isfinite(centerlineTolerance) || centerlineTolerance < 0 || ...
            isempty(sourceSetKeys) || any(strlength(sourceSetKeys)==0) || ...
            numel(unique(sourceSetKeys)) ~= numel(sourceSetKeys)
        error('evalLibrary:FullCoordinateCandidateConfig', ...
            'Candidate centerline has invalid bounds, tolerance, or source-set keys.');
    end
    query = double(query(:).');
    activeKey = strjoin(sourceIds,'|');
    matchedCenterline = query(1) >= speedRange(1)-centerlineTolerance && ...
        query(1) <= speedRange(2)+centerlineTolerance && ...
        abs(query(2)-(alphaOffset-query(1))) <= centerlineTolerance && ...
        nnz(sourceSetKeys == activeKey) == 1;
end
if ~(matchedStencil || matchedCenterline)
    error('evalLibrary:FullCoordinateCandidateStencil', ...
        ['Full-coordinate runtime candidate is enabled, but query=[%.12g %.12g] ', ...
         'and its active source set are outside its approved exact/centerline selector.'], ...
        query(1),query(2));
end
if ~isfile(candidate.registryPath) || ~isfolder(candidate.fieldRoot)
    error('evalLibrary:FullCoordinateCandidatePath', ...
        'Candidate registryPath and fieldRoot must exist.');
end
atomicOptions = struct('enabled',true, ...
    'architecture',"full_coordinate_atomic_lift_interpolate_project", ...
    'registryPath',string(candidate.registryPath), ...
    'fieldRoot',string(candidate.fieldRoot));
if isfield(candidate,'cacheImmutableMetadata')
    atomicOptions.cacheImmutableMetadata = logical(candidate.cacheImmutableMetadata);
end
if isfield(candidate,'cacheImmutableSourceData')
    atomicOptions.cacheImmutableSourceData = logical(candidate.cacheImmutableSourceData);
end
preparedFields = {'preparedRuntimeOwner','preparedRuntimeProfile', ...
    'preparedRuntimeChangeId','preparedRuntimeRegistrySha256', ...
    'projectedTustinCondensation'};
for fieldIndex = 1:numel(preparedFields)
    fieldName = preparedFields{fieldIndex};
    if isfield(candidate,fieldName)
        atomicOptions.(fieldName) = candidate.(fieldName);
    end
end
enabled = true;
end

function sched = installActiveCellRecovery(sched,points,weights)
% Recovery vertices retain their own raw-L/Pr/phi ownership. Query weights
% multiply the source-local vertex weights, while the physical recentering
% origin is constructed from the same complete source fields.
vertexCount = 0;
for k = 1:numel(points)
    if isfield(points(k),'equilibriumCentered') && ...
            isfield(points(k).equilibriumCentered,'recoveryVertices')
        vertexCount = vertexCount+ ...
            numel(points(k).equilibriumCentered.recoveryVertices);
    end
end

vertices = repmat(struct('weight',0,'xEq',[],'uEq',[],'Lq1',[], ...
    'parConst',struct(),'idx',struct(),'Pz',[],'Pr',[], ...
    'phi1_sA',[]),vertexCount,1);
vertexIndex = 0;
referenceState = zeros(size(points(1).x_eq(:)));
referenceControl = zeros(size(points(1).u_eq(:)));
recoveryAnchorWrench = zeros(6,1);
for k = 1:numel(points)
    point = points(k);
    if ~isfield(point,'equilibriumCentered') || ...
            ~logical(point.equilibriumCentered.enabled) || ...
            ~isfield(point.equilibriumCentered,'recoveryVertices') || ...
            isempty(point.equilibriumCentered.recoveryVertices)
        error('evalLibrary:IncompatibleBeamMaps', ...
            'Source %s has no certified equilibrium-centered recovery.', ...
            string(point.name));
    end
    sourceVertices = point.equilibriumCentered.recoveryVertices;
    sourceWeight = [sourceVertices.weight].';
    if abs(sum(sourceWeight)-1)>1e-12
        error('evalLibrary:IncompatibleBeamMaps', ...
            'Source %s recovery-vertex weights do not sum to one.', ...
            string(point.name));
    end
    localState = zeros(size(referenceState));
    localControl = zeros(size(referenceControl));
    for j = 1:numel(sourceVertices)
        vertex = sourceVertices(j);
        vertex.weight = weights(k)*vertex.weight;
        vertexIndex = vertexIndex+1;
        vertices(vertexIndex) = vertex;
        localState = localState+sourceWeight(j)*vertex.xEq(:);
        localControl = localControl+sourceWeight(j)*vertex.uEq(:);
    end
    if isfield(point,'physicalRecoveryReferenceState') && ...
            isfield(point,'physicalRecoveryReferenceControl')
        localState = point.physicalRecoveryReferenceState(:);
        localControl = point.physicalRecoveryReferenceControl(:);
    end
    referenceState = referenceState+weights(k)*localState;
    referenceControl = referenceControl+weights(k)*localControl;
    if ~isfield(point,'p5') || ~isfield(point.p5,'r2') || ...
            ~isfield(point.p5.r2,'anchor') || ...
            ~isfield(point.p5.r2.anchor,'wrench')
        error('evalLibrary:IncompatibleBeamMaps', ...
            'Source %s has no certified T2 recovery anchor.', ...
            string(point.name));
    end
    sourceAnchor = point.p5.r2.anchor.wrench(:);
    if numel(sourceAnchor)~=6 || any(~isfinite(sourceAnchor))
        error('evalLibrary:IncompatibleBeamMaps', ...
            'Source %s has an invalid T2 recovery anchor.', ...
            string(point.name));
    end
    recoveryAnchorWrench = recoveryAnchorWrench+weights(k)*sourceAnchor;
end

if abs(sum([vertices.weight])-1)>1e-12
    error('evalLibrary:IncompatibleBeamMaps', ...
        'Active-cell recovery weights do not sum to one.');
end
sched.equilibriumCentered = struct( ...
    'enabled',true,'endpointExact',false, ...
    'vertexIds',{cellstr(string({points.name}))}, ...
    'weights',weights(:), ...
    'rEq',zeros(size(sched.L,1),1), ...
    'form','ACTIVE_CELL_COMMON_COORDINATE_RECOVERY_V1', ...
    'recoveryVertices',vertices, ...
    'recoveryAnchorWrench',recoveryAnchorWrench);
sched.physicalRecoveryReferenceState = referenceState;
sched.physicalRecoveryReferenceControl = referenceControl;
sched.recoveryConstruction = struct( ...
    'path','A_COMMON_COORDINATE_FIELD_CONSTRUCTION', ...
    'sourceIds',{cellstr(string({points.name}))}, ...
    'weights',weights(:),'atomicWithQueryPackage',true);
end

function sched = installExactNodeRecovery(sched,points)
% Exact-node evaluation must retain the source-owned physical-recovery
% contract.  The query package above already records query provenance; do
% not reconstruct or translate a one-source recovery reference here.
point = points(1);
assert(isfield(point,'equilibriumCentered'), ...
    'evalLibrary:ExactNodeRecoveryContract', ...
    'Exact source %s omits its equilibrium-centered recovery contract.', ...
    string(point.name));
assert(isstruct(point.equilibriumCentered) && ...
    isfield(point.equilibriumCentered,'enabled') && ...
    logical(point.equilibriumCentered.enabled), ...
    'evalLibrary:ExactNodeRecoveryContract', ...
    'Exact source %s has no enabled equilibrium-centered recovery.', ...
    string(point.name));
sched.equilibriumCentered = point.equilibriumCentered;
if ~isfield(sched.equilibriumCentered,'recoveryAnchorWrench')
    assert(isfield(point,'p5') && isfield(point.p5,'r2') && ...
        isfield(point.p5.r2,'anchor') && ...
        isfield(point.p5.r2.anchor,'wrench'), ...
        'evalLibrary:ExactNodeRecoveryAnchor', ...
        'Exact source %s has no certified T2 recovery anchor.', ...
        string(point.name));
    anchor = point.p5.r2.anchor.wrench(:);
    assert(numel(anchor) == 6 && all(isfinite(anchor)), ...
        'evalLibrary:ExactNodeRecoveryAnchor', ...
        'Exact source %s has an invalid certified T2 recovery anchor.', ...
        string(point.name));
    sched.equilibriumCentered.recoveryAnchorWrench = anchor;
end
if isfield(point,'physicalRecoveryReferenceState')
    sched.physicalRecoveryReferenceState = ...
        point.physicalRecoveryReferenceState;
end
if isfield(point,'physicalRecoveryReferenceControl')
    sched.physicalRecoveryReferenceControl = ...
        point.physicalRecoveryReferenceControl;
end
if isfield(point,'recoveryConstruction')
    sched.recoveryConstruction = point.recoveryConstruction;
end
end

function [verified,details] = verifyPhysicalRecoveryCertificate(ROMlib,P,ids,info,registryResolution)
details = struct('state','MISSING','cellId',cellIdentifier(info), ...
    'firstField','physicalRecoveryCertificate', ...
    'auditPath',defaultAuditPath());
verified = isscalar(P);
if verified
    details.state = 'EXACT_NODE_NOT_REQUIRED';
    return;
end
if ~isfield(ROMlib,'physicalRecoveryCertificates') || ...
        isempty(ROMlib.physicalRecoveryCertificates)
    if isfield(ROMlib,'physicalRecoveryCompatible') && ...
            isscalar(ROMlib.physicalRecoveryCompatible) && ...
            logical(ROMlib.physicalRecoveryCompatible)
        details.state = 'TRUE_BUT_STALE_LIBRARY_WIDE_FLAG';
    end
    return;
end

certificates = ROMlib.physicalRecoveryCertificates;
sourceIds = string({P.name});
match = false(1,numel(certificates));
for k = 1:numel(certificates)
    certificate = certificates(k);
    if ~isfield(certificate,'sourceIds'), continue; end
    match(k) = sameStringSet(string(certificate.sourceIds),sourceIds) && ...
        certificateCellMatches(certificate,info);
end
index = find(match,1);
if isempty(index)
    details.state = 'MISSING_FOR_ACTIVE_CELL';
    return;
end
certificate = certificates(index);
details.state = 'FALSE';
if isfield(certificate,'auditPath') && strlength(string(certificate.auditPath))>0
    details.auditPath = char(string(certificate.auditPath));
end
if isfield(certificate,'firstIncompatibleField') && ...
        strlength(string(certificate.firstIncompatibleField))>0
    details.firstField = char(string(certificate.firstIncompatibleField));
end
requiredLogical = {'passed','semanticAmbientPassed', ...
    'physicalRecoveryPassed','virtualPowerPassed','petrovGalerkinPassed'};
for k = 1:numel(requiredLogical)
    field = requiredLogical{k};
    if ~isfield(certificate,field) || ~isscalar(certificate.(field)) || ...
            ~logical(certificate.(field))
        return;
    end
end
requiredText = {'auditHash','chartHash','semanticMapHash','ambientMapHash', ...
    'delayArchitectureHash','delayOverlayHash','evalLibraryCodeHash', ...
    'interpWeightsCodeHash', ...
    'physicalRecoveryAuditHash','virtualPowerAuditHash', ...
    'petrovGalerkinAuditHash','certificateTimestamp','certificateSchemaVersion'};
for k = 1:numel(requiredText)
    field = requiredText{k};
    if ~isfield(certificate,field) || strlength(string(certificate.(field)))==0
        details.state = 'TRUE_BUT_STALE_MISSING_HASH';
        return;
    end
end
if ~isfield(ROMlib,'ambientContract') || ...
        ~isstruct(ROMlib.ambientContract) || ...
        ~isscalar(ROMlib.ambientContract)
    details.state = 'TRUE_BUT_STALE_LIBRARY_AMBIENT_CONTRACT_MISSING';
    details.firstField = 'ambientContract';
    return;
end
contractFields = {'auditHash','chartHash','semanticMapHash','ambientMapHash', ...
    'delayArchitectureHash','delayOverlayHash','evalLibraryCodeHash', ...
    'interpWeightsCodeHash', ...
    'physicalRecoveryAuditHash','virtualPowerAuditHash', ...
    'petrovGalerkinAuditHash','certificateSchemaVersion'};
for k = 1:numel(contractFields)
    field = contractFields{k};
    if ~isfield(ROMlib.ambientContract,field) || ...
            strlength(string(ROMlib.ambientContract.(field)))==0 || ...
            string(certificate.(field))~= ...
            string(ROMlib.ambientContract.(field))
        details.state = 'TRUE_BUT_STALE_AMBIENT_HASH_MISMATCH';
        details.firstField = field;
        return;
    end
end
currentEvalHash = fileSha256(which('AeroFlex.sched.evalLibrary'));
currentWeightsHash = fileSha256(which('AeroFlex.sched.interpWeights'));
if string(certificate.evalLibraryCodeHash) ~= currentEvalHash
    details.state = 'TRUE_BUT_STALE_CODE_HASH_MISMATCH';
    details.firstField = 'evalLibraryCodeHash';
    return;
end
if string(certificate.interpWeightsCodeHash) ~= currentWeightsHash
    details.state = 'TRUE_BUT_STALE_CODE_HASH_MISMATCH';
    details.firstField = 'interpWeightsCodeHash';
    return;
end
if ~isfield(ROMlib.ambientContract,'delayOverlayPath') || ...
        strlength(string(ROMlib.ambientContract.delayOverlayPath))==0
    details.state = 'TRUE_BUT_STALE_MISSING_HASH';
    details.firstField = 'delayOverlayPath';
    return;
end
overlayPath = char(string(ROMlib.ambientContract.delayOverlayPath));
if ~isfile(overlayPath) || fileSha256(overlayPath) ~= ...
        string(certificate.delayOverlayHash)
    details.state = 'TRUE_BUT_STALE_DELAY_OVERLAY_HASH_MISMATCH';
    details.firstField = 'delayOverlayHash';
    return;
end
if ~isfield(certificate,'sourcePackageHashes') || ...
        numel(certificate.sourcePackageHashes)~=numel(sourceIds)
    details.state = 'TRUE_BUT_STALE_MISSING_SOURCE_HASHES';
    return;
end
certificateIds = string(certificate.sourceIds);
certificateHashes = string(certificate.sourcePackageHashes);
for k = 1:numel(ids)
    sourceIndex = find(certificateIds==sourceIds(k),1);
    if ~isempty(fieldnames(registryResolution))
        registryIndex = find(registryResolution.sourceIds==sourceIds(k),1);
        currentHash = registryResolution.sourceContractHashes(registryIndex);
    else
        currentHash = pointPackageHash(ROMlib.points(ids(k)));
    end
    if isempty(sourceIndex) || strlength(currentHash)==0 || ...
            certificateHashes(sourceIndex)~=currentHash
        details.state = 'TRUE_BUT_STALE_SOURCE_HASH_MISMATCH';
        details.firstField = 'sourcePackageHash';
        return;
    end
end
verified = true;
details.state = 'TRUE_AND_VERIFIED';
details.firstField = 'none';
end

function value = fileSha256(path)
if isempty(path) || ~isfile(path)
    value = "";
    return;
end
fileId = fopen(path,'rb');
if fileId < 0
    value = "";
    return;
end
cleanup = onCleanup(@() fclose(fileId));
data = fread(fileId,Inf,'*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(data,'int8'));
bytes = typecast(engine.digest(),'uint8');
value = lower(string(reshape(dec2hex(bytes,2).',1,[])));
end

function pass = certificateCellMatches(certificate,info)
pass = true;
if isfield(certificate,'cellId')
    pass = string(certificate.cellId)==string(cellIdentifier(info));
end
end

function pass = sameStringSet(left,right)
left = sort(left(:)); right = sort(right(:));
pass = numel(left)==numel(right) && all(left==right);
end

function value = pointPackageHash(point)
value = "";
fields = {'physicalRecoverySourceHash','packageHash','package_hash'};
for k = 1:numel(fields)
    if isfield(point,fields{k}) && strlength(string(point.(fields{k})))>0
        value = string(point.(fields{k}));
        return;
    end
end
if isfield(point,'provenance') && isstruct(point.provenance) && ...
        isfield(point.provenance,'packageHash')
    value = string(point.provenance.packageHash);
end
end

function value = cellIdentifier(info)
if isfield(info,'simplex') && isscalar(info.simplex) && isfinite(info.simplex)
    value = char(string(info.simplex));
elseif isfield(info,'mode')
    value = char(string(info.mode));
else
    value = 'UNKNOWN';
end
end

function value = defaultAuditPath()
value = ['context/audits/phase18b-interpolation-attribution/', ...
    'EVAL_LIBRARY_DIAGNOSTIC.json'];
end

function assertSameMapSize(actual,reference,fieldName,actualPoint, ...
        referencePoint,mu,info)
if isequal(size(actual),size(reference)), return; end
error('evalLibrary:IncompatibleBeamMaps', ...
    ['Physical recovery dimension mismatch at query=[%.12g %.12g], ', ...
     'cell=%s: %s is %dx%d for %s but %dx%d for %s. ', ...
     'Exact-node evaluation remains available. Audit artifact: %s'], ...
    mu(1),mu(2),cellIdentifier(info),fieldName,size(actual,1),size(actual,2), ...
    char(string(actualPoint.name)),size(reference,1),size(reference,2), ...
    char(string(referencePoint.name)),defaultAuditPath());
end
