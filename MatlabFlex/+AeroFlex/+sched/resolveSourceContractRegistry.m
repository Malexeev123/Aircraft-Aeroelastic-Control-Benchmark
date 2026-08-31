function resolved = resolveSourceContractRegistry(registryInput,sourceIds)
%RESOLVESOURCECONTRACTREGISTRY Resolve fixed10 source contracts fail-closed.

arguments
    registryInput
    sourceIds (1,:) string
end

if isstruct(registryInput)
    registry = registryInput;
    registryPath = "IN_MEMORY";
    registryHash = "IN_MEMORY";
else
    registryPath = string(registryInput);
    if ~isscalar(registryPath) || strlength(registryPath)==0 || ...
            ~isfile(registryPath)
        error('resolveSourceContractRegistry:MissingRegistry', ...
            'The source-contract registry does not exist: %s',registryPath);
    end
    registry = jsondecode(fileread(registryPath));
    registryHash = fileHash(registryPath);
end
requiredTop = {'status','registryVersion','requiredSourceCount','sources', ...
    'aliases','dependencyVersion'};
requireFields(registry,requiredTop,'registry');
if string(registry.status)~="LOCKED"
    error('resolveSourceContractRegistry:NotLocked', ...
        'The source-contract registry status is %s, not LOCKED.', ...
        string(registry.status));
end
supportedVersions = ["phase18b-fixed10-source-contract-registry-v1", ...
    "phase18b-fixed-node27-source-contract-registry-v1"];
if ~isscalar(string(registry.registryVersion)) || ...
        ~any(string(registry.registryVersion)==supportedVersions)
    error('resolveSourceContractRegistry:UnsupportedRegistryVersion', ...
        'Unsupported source-contract registry version: %s.', ...
        string(registry.registryVersion));
end
sources = registry.sources;
if numel(sources)~=double(registry.requiredSourceCount)
    error('resolveSourceContractRegistry:CountMismatch', ...
        'Registry source count %d does not match declared count %d.', ...
        numel(sources),double(registry.requiredSourceCount));
end
allIds = string({sources.sourceId});
if numel(unique(allIds))~=numel(allIds)
    error('resolveSourceContractRegistry:DuplicateIdentity', ...
        'The registry contains duplicate source identities.');
end
if ~isempty(fieldnames(registry.aliases))
    error('resolveSourceContractRegistry:AmbiguousAlias', ...
        'Source aliases are prohibited in the locked Phase-18B registry.');
end

root = repositoryRoot();
rows = sources(ones(1,numel(sourceIds)));
hashes = strings(1,numel(sourceIds));
for k = 1:numel(sourceIds)
    index = find(allIds==sourceIds(k));
    if isempty(index)
        error('resolveSourceContractRegistry:MissingSource', ...
            'No registry entry exists for source %s.',sourceIds(k));
    end
    if numel(index)~=1
        error('resolveSourceContractRegistry:DuplicateIdentity', ...
            'Source %s resolves ambiguously.',sourceIds(k));
    end
    row = sources(index);
    requireFields(row,{'fixedOrder','sourceArtifactPath', ...
        'sourceArtifactSha256','fullTupleHashes','reducedTupleHashes', ...
        'VHash','WHash','semanticHash','ambientMapHash','delayOverlayHash', ...
        'equilibriumHash','affineHash','recoveryHash','rootReactionHash', ...
        'tensorHash','PzHash','PrHash','rawLHash','LdynHash','T2Hash'}, ...
        sprintf('source %s',sourceIds(k)));
    if double(row.fixedOrder)~=10
        error('resolveSourceContractRegistry:WrongFixedOrder', ...
            'Source %s has fixed order %g; fixed10 is required.', ...
            sourceIds(k),double(row.fixedOrder));
    end
    requiredHashes = {'VHash','WHash','semanticHash','ambientMapHash', ...
        'delayOverlayHash','equilibriumHash','affineHash','recoveryHash', ...
        'rootReactionHash','tensorHash','PzHash','PrHash','rawLHash', ...
        'LdynHash','T2Hash'};
    for j = 1:numel(requiredHashes)
        field = requiredHashes{j};
        if strlength(string(row.(field)))==0
            error('resolveSourceContractRegistry:MissingHash', ...
                'Source %s is missing %s.',sourceIds(k),field);
        end
    end
    requireTupleHashes(row.fullTupleHashes,sourceIds(k),'full');
    requireTupleHashes(row.reducedTupleHashes,sourceIds(k),'reduced');
    artifactPath = string(row.sourceArtifactPath);
    if ~isfile(artifactPath)
        artifactPath = fullfile(root,artifactPath);
    end
    if ~isfile(artifactPath)
        error('resolveSourceContractRegistry:MissingArtifact', ...
            'Source artifact is missing for %s: %s',sourceIds(k),artifactPath);
    end
    observed = fileHash(artifactPath);
    if observed~=string(row.sourceArtifactSha256)
        error('resolveSourceContractRegistry:StaleArtifactHash', ...
            'Source artifact hash mismatch for %s.',sourceIds(k));
    end
    rows(k) = row;
    hashes(k) = string(row.sourceArtifactSha256);
end
resolved = struct('registryPath',registryPath,'registryHash',registryHash, ...
    'registryVersion',string(registry.registryVersion), ...
    'sourceIds',sourceIds,'sourceContractHashes',hashes,'sources',rows, ...
    'dependencyVersion',registry.dependencyVersion);
end

function requireTupleHashes(value,sourceId,label)
requireFields(value,{'A','B','C','D'},sprintf('%s tuple for %s',label,sourceId));
for field = {'A','B','C','D'}
    if strlength(string(value.(field{1})))==0
        error('resolveSourceContractRegistry:MissingHash', ...
            'Source %s is missing the %s %s hash.',sourceId,label,field{1});
    end
end
end

function requireFields(value,fields,label)
for k = 1:numel(fields)
    if ~isfield(value,fields{k})
        error('resolveSourceContractRegistry:MissingField', ...
            '%s is missing required field %s.',label,fields{k});
    end
end
end

function root = repositoryRoot()
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end

function value = fileHash(path)
fileId = fopen(path,'rb');
if fileId<0
    error('resolveSourceContractRegistry:HashRead', ...
        'Unable to open source-contract artifact: %s',path);
end
cleanup = onCleanup(@() fclose(fileId));
data = fread(fileId,Inf,'*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(typecast(data,'int8'));
bytes = typecast(engine.digest(),'uint8');
value = lower(string(reshape(dec2hex(bytes,2).',1,[])));
end
