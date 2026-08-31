function assertPackageConsistency(sched)
%ASSERTPACKAGECONSISTENCY Validate one atomic scheduled ROM package.

required = {'mu','weights','pointIds','pointMu','L','idx','parConst', ...
    'beam','base','aero','x_eq','u_eq','time','compatibleCoordinates'};
for k = 1:numel(required)
    if ~isfield(sched,required{k})
        error('assertPackageConsistency:MissingField', ...
            'Scheduled package is missing field "%s".',required{k});
    end
end

if ~sched.compatibleCoordinates
    error('assertPackageConsistency:Coordinates', ...
        'Scheduled package is not in declared compatible coordinates.');
end
if abs(sum(sched.weights)-1) > 1e-13 || any(sched.weights < -1e-13)
    error('assertPackageConsistency:Weights', ...
        'Interpolation weights violate the Phase-18 acceptance budget.');
end
if numel(sched.weights) ~= numel(sched.pointIds) || ...
        size(sched.pointMu,1) ~= numel(sched.pointIds)
    error('assertPackageConsistency:Provenance', ...
        'Source-node IDs, coordinates, and weights are inconsistent.');
end

nx = size(sched.L,1);
if size(sched.L,2) ~= nx || numel(sched.x_eq) ~= nx
    error('assertPackageConsistency:StateSize', ...
        'L and x_eq do not share one state dimension.');
end
blocks = {'q1','q2','qxi','qGam','chi'};
allIdx = [];
for k = 1:numel(blocks)
    if ~isfield(sched.idx,blocks{k})
        error('assertPackageConsistency:StateOrdering', ...
            'State index block "%s" is missing.',blocks{k});
    end
    allIdx = [allIdx sched.idx.(blocks{k})]; %#ok<AGROW>
end
if ~isequal(allIdx,1:nx)
    error('assertPackageConsistency:StateOrdering', ...
        'State index blocks are not a complete ordered partition.');
end

dt = sched.parConst.dt;
if ~isscalar(dt) || ~isfinite(dt) || dt <= 0 || dt ~= sched.time.runtimeDt
    error('assertPackageConsistency:Timestep', ...
        'Package runtime dt is invalid or inconsistent with time provenance.');
end
if ~isfield(sched.beam,'Pz') || ~isfield(sched.beam,'Pr') || ...
        ~isfield(sched.beam,'red')
    error('assertPackageConsistency:StructuralMaps', ...
        'Fixed structural projection/recovery maps are incomplete.');
end
if isfield(sched,'equilibriumCentered') && ...
        isfield(sched.equilibriumCentered,'enabled') && ...
        logical(sched.equilibriumCentered.enabled)
    centered = sched.equilibriumCentered;
    if ~isfield(centered,'recoveryVertices') || ...
            isempty(centered.recoveryVertices) || ...
            ~isfield(centered,'recoveryAnchorWrench') || ...
            numel(centered.recoveryAnchorWrench)~=6 || ...
            any(~isfinite(centered.recoveryAnchorWrench))
        error('assertPackageConsistency:PhysicalRecoveryAnchor', ...
            ['Enabled centered recovery requires vertices and one finite ', ...
            'six-component absolute T2 anchor.']);
    end
end

mapPairs = {'Bw','Bw';'Dw','Dw';'Bdel','B_delta';'Ddel','D_delta'; ...
    'Bddel','B_ddelta';'Dddel','D_ddelta'};
for k = 1:size(mapPairs,1)
    parName = mapPairs{k,1};
    aeroName = mapPairs{k,2};
    if ~isequaln(sched.parConst.(parName),sched.aero.forceMap.(aeroName))
        error('assertPackageConsistency:ForceMaps', ...
            'parConst.%s and aero.forceMap.%s differ.',parName,aeroName);
    end
end

numericFields = {'L','x_eq','u_eq'};
for k = 1:numel(numericFields)
    if any(~isfinite(sched.(numericFields{k})(:)))
        error('assertPackageConsistency:NonFinite', ...
            'Scheduled field "%s" contains NaN/Inf.',numericFields{k});
    end
end
end
