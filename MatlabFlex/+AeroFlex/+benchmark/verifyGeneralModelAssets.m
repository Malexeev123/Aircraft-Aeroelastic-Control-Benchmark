function status = verifyGeneralModelAssets(repositoryRoot,options)
%VERIFYGENERALMODELASSETS Verify the shared wing/coupled model-data payload.

arguments
    repositoryRoot (1,1) string
    options.PrintSummary (1,1) logical = false
end

manifestPath = fullfile(repositoryRoot,"MatlabFlex","configs", ...
    "benchmark","pazy_general_model_assets_v1.json");
assert(isfile(manifestPath),"AeroFlex:GeneralModelAssetManifest", ...
    "General model-asset manifest is missing: %s",manifestPath);
manifest = jsondecode(fileread(manifestPath));
assert(string(manifest.schemaVersion)=="pazy-general-model-assets-v1" && ...
    isfield(manifest,"assets") && ~isempty(manifest.assets), ...
    "AeroFlex:GeneralModelAssetManifestSchema", ...
    "General model-asset manifest has an unsupported schema.");

assets = manifest.assets(:);
records = repmat(struct("path","","expectedSha256","", ...
    "actualSha256","","exists",false,"passed",false),numel(assets),1);
for index = 1:numel(assets)
    relativePath = string(assets(index).path);
    expectedHash = lower(string(assets(index).sha256));
    absolutePath = fullfile(repositoryRoot,localNativePath(relativePath));
    exists = isfile(absolutePath);
    actualHash = "";
    if exists
        actualHash = localFileHash(absolutePath);
    end
    records(index) = struct("path",relativePath, ...
        "expectedSha256",expectedHash,"actualSha256",actualHash, ...
        "exists",exists,"passed",exists && actualHash==expectedHash);
end

status = struct( ...
    "schemaVersion","pazy-general-model-asset-status-v1", ...
    "manifestPath",string(manifestPath), ...
    "caseName",string(manifest.caseName), ...
    "assetCount",numel(records), ...
    "passed",all([records.passed]), ...
    "records",records);
if options.PrintSummary
    fprintf("  Shared model data    : %s (%d files)\n", ...
        localPassFail(status.passed),status.assetCount);
end
end

function value = localFileHash(path)
file = fopen(path,"rb");
assert(file>=0,"AeroFlex:GeneralModelAssetRead", ...
    "Cannot read model asset for hashing: %s",path);
cleanup = onCleanup(@()fclose(file));
engine = javaMethod("getInstance", ...
    "java.security.MessageDigest","SHA-256");
while ~feof(file)
    bytes = fread(file,1024*1024,"*uint8");
    if isempty(bytes), break, end
    engine.update(typecast(bytes(:),"int8"));
end
value = lower(string(reshape(dec2hex(typecast( ...
    engine.digest(),"uint8"),2).',1,[])));
clear cleanup
end

function value = localNativePath(value)
value = replace(string(value),["/","\"],filesep);
end

function value = localPassFail(flag)
if flag, value = "PASS"; else, value = "NOT READY"; end
end
