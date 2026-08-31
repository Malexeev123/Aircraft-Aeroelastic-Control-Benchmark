function [ambient,certificate] = ...
        resolveActiveCellPhysicalRecoveryCertificate( ...
        registryPath,caseView,sourceIds,weights,sourceRegistryResolution)
%RESOLVEACTIVECELLPHYSICALRECOVERYCERTIFICATE Resolve one governed cell.

arguments
    registryPath (1,1) string
    caseView (1,1) struct
    sourceIds (1,:) string
    weights (:,1) double
    sourceRegistryResolution (1,1) struct
end
if ~isfile(registryPath)
    error('resolveActiveCellPhysicalRecoveryCertificate:MissingRegistry', ...
        'Active-cell recovery registry is missing: %s',registryPath);
end
registry=jsondecode(fileread(registryPath));
required={'status','certificateSchemaVersion','versionedRegistrySha256', ...
    'caseViewRegistrySha256','ambientContract','certificates'};
requireFields(registry,required,'active-cell registry');
if string(registry.status)~="LOCKED_FOR_EXECUTION"
    error('resolveActiveCellPhysicalRecoveryCertificate:NotLocked', ...
        'Active-cell recovery registry status is %s.',registry.status);
end
if ~isfield(sourceRegistryResolution,'registryHash') || ...
        string(sourceRegistryResolution.registryHash)~= ...
        string(registry.versionedRegistrySha256)
    error('resolveActiveCellPhysicalRecoveryCertificate:ParentHash', ...
        'Source and active-cell recovery registries have different hashes.');
end
requireFields(caseView,{'caseManifestSha256'},'case view');
caseHash=string(caseView.caseManifestSha256);
rows=registry.certificates;
match=false(1,numel(rows));
for k=1:numel(rows)
    allowed=string(rows(k).caseViewHashes);
    match(k)=any(allowed==caseHash) && ...
        sameStringSet(string(rows(k).sourceIds),sourceIds);
end
index=find(match);
if numel(index)~=1
    error('resolveActiveCellPhysicalRecoveryCertificate:MissingCell', ...
        'Case view %s resolves to %d active-cell certificates.',caseHash,numel(index));
end
row=rows(index);
caseIndex=find(string(row.caseViewHashes)==caseHash,1);
expected=caseWeightAt(row.caseWeights,caseIndex,numel(row.caseViewHashes));
if numel(expected)~=numel(weights) || ...
        norm(expected-weights(:),inf)>1e-12
    error('resolveActiveCellPhysicalRecoveryCertificate:WeightMismatch', ...
        'Runtime weights differ from the locked case-view contract.');
end
if ~logical(row.passed) || string(row.status)~="CERTIFIED"
    error('resolveActiveCellPhysicalRecoveryCertificate:NotCertified', ...
        'Active stencil %s is not certified.',row.stencilId);
end
ambient=registry.ambientContract;
certificate=row.certificate;
certificate.cellId=string(caseView.expectedInterpolationMode);
certificate.sourceIds=sourceIds;
certificate.sourcePackageHashes= ...
    sourceRegistryResolution.sourceContractHashes;
end

function value=caseWeightAt(caseWeights,caseIndex,caseCount)
if iscell(caseWeights)
    value=double(caseWeights{caseIndex}(:));
elseif caseCount==1
    value=double(caseWeights(:));
elseif size(caseWeights,1)==caseCount
    value=double(caseWeights(caseIndex,:).');
elseif size(caseWeights,2)==caseCount
    value=double(caseWeights(:,caseIndex));
else
    error('resolveActiveCellPhysicalRecoveryCertificate:WeightSchema', ...
        'The locked case-weight array cannot be indexed by case view.');
end
end

function pass=sameStringSet(left,right)
left=sort(left(:)); right=sort(right(:));
pass=numel(left)==numel(right) && all(left==right);
end

function requireFields(value,fields,label)
for k=1:numel(fields)
    if ~isfield(value,fields{k})
        error('resolveActiveCellPhysicalRecoveryCertificate:MissingField', ...
            '%s is missing %s.',label,fields{k});
    end
end
end
