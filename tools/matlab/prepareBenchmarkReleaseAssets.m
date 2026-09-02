function status = prepareBenchmarkReleaseAssets(options)
%PREPAREBENCHMARKRELEASEASSETS Verify or stage locked benchmark data.
%   STATUS = PREPAREBENCHMARKRELEASEASSETS() checks every selected runtime
%   asset against its accepted SHA-256 without copying files.
%
%   STATUS = PREPAREBENCHMARKRELEASEASSETS(Action="stage", ...
%   DestinationRoot=FOLDER) copies the verified bytes into FOLDER while
%   preserving repository-relative paths. Existing matching files are reused;
%   an existing mismatch fails closed and is never overwritten.

arguments
    options.Action (1,1) string {mustBeMember(options.Action, ...
        ["check","stage"])} = "check"
    options.DestinationRoot (1,1) string = ""
    options.ProjectInfo (1,1) struct = struct()
    options.PrintSummary (1,1) logical = true
end

if isempty(fieldnames(options.ProjectInfo))
    project = setupProject(ValidateEntryPoints=true);
else
    project = options.ProjectInfo;
end
root = string(project.repositoryRoot);
manifestPath = fullfile(root,"MatlabFlex","configs","benchmark", ...
    "phase18c_v17a_release_asset_manifest_v1.json");
manifest = jsondecode(fileread(manifestPath));
assert(string(manifest.status)=="LOCKED_SOURCE_EQUIVALENT_SELECTION", ...
    "AeroFlex:ReleaseAssetsManifest", ...
    "The release-asset manifest is not locked for staging.");

records = localDirectRecords(manifest.directAssets);
records = [records;localV17aSourceRecords(root,manifest.sourceGroups.v17aSources)];
records = [records;localRecordGroup(root, ...
    manifest.sourceGroups.fullCoordinateFields,"path","sha256")];
records = [records;localRecordGroup(root, ...
    manifest.sourceGroups.scheduledReciprocalBank, ...
    "runtimePath","runtimeSha256")];
records = localUniqueRecords(records);

for index = 1:numel(records)
    source = fullfile(root,localNativePath(records(index).path));
    records(index).exists = isfile(source);
    records(index).observedSha256 = "";
    records(index).passed = false;
    if records(index).exists
        records(index).observedSha256 = localFileHash(source);
        records(index).passed = ...
            records(index).observedSha256==records(index).expectedSha256;
    end
end

if ~all([records.passed])
    failed = records(~[records.passed]);
    details = strings(numel(failed),1);
    for index = 1:numel(failed)
        if failed(index).exists
            details(index) = failed(index).path + " (hash mismatch)";
        else
            details(index) = failed(index).path + " (missing)";
        end
    end
    error("AeroFlex:ReleaseAssetsVerification", ...
        "Selected release assets did not verify:\n%s", ...
        strjoin(details,newline));
end

stagedCount = 0;
reusedCount = 0;
if options.Action=="stage"
    destination = localValidateDestination(root,options.DestinationRoot);
    for index = 1:numel(records)
        source = fullfile(root,localNativePath(records(index).path));
        target = fullfile(destination,localNativePath(records(index).path));
        folder = string(fileparts(target));
        if ~isfolder(folder), mkdir(folder); end
        if isfile(target)
            assert(localFileHash(target)==records(index).expectedSha256, ...
                "AeroFlex:ReleaseAssetsExistingMismatch", ...
                "Refusing to overwrite a changed staged asset: %s",target);
            reusedCount = reusedCount+1;
        else
            [copied,message] = copyfile(source,target);
            assert(copied,"AeroFlex:ReleaseAssetsCopy", ...
                "Unable to stage %s: %s",source,message);
            assert(localFileHash(target)==records(index).expectedSha256, ...
                "AeroFlex:ReleaseAssetsPostCopyHash", ...
                "The staged asset hash changed during copying: %s",target);
            stagedCount = stagedCount+1;
        end
    end
    localWriteInventory(destination,manifestPath,records);
else
    destination = "";
end

status = struct( ...
    "schemaVersion","pazy-benchmark-release-assets-status-v1", ...
    "action",options.Action,"passed",all([records.passed]), ...
    "assetCount",numel(records), ...
    "totalBytes",localTotalBytes(root,records), ...
    "stagedCount",stagedCount,"reusedCount",reusedCount, ...
    "destinationRoot",destination,"records",records);

if options.PrintSummary
    fprintf("\nPazy benchmark release assets\n");
    fprintf("  Action          : %s\n",options.Action);
    fprintf("  Verified assets : %d\n",status.assetCount);
    fprintf("  Total size      : %.1f MiB\n",status.totalBytes/2^20);
    if options.Action=="stage"
        fprintf("  Newly staged    : %d\n",stagedCount);
        fprintf("  Reused          : %d\n",reusedCount);
        fprintf("  Destination     : %s\n",destination);
    end
    fprintf("  Status          : PASS\n\n");
end
end

function records = localDirectRecords(values)
records = repmat(localRecord(),numel(values),1);
for index = 1:numel(values)
    records(index) = localRecord(string(values(index).path), ...
        string(values(index).sha256),string(values(index).role));
end
end

function records = localV17aSourceRecords(root,group)
manifestPath = fullfile(root,localNativePath(string(group.manifestPath)));
assert(localFileHash(manifestPath)==string(group.manifestSha256), ...
    "AeroFlex:ReleaseAssetsSourceManifest", ...
    "The V17A source manifest changed.");
manifest = jsondecode(fileread(manifestPath));
assert(numel(manifest.sources)==double(group.expectedSourceCount), ...
    "AeroFlex:ReleaseAssetsSourceCount", ...
    "The V17A source count differs from the release contract.");
pathFields = string(group.pathFields(:));
assert(~isempty(pathFields),"AeroFlex:ReleaseAssetsSourceFields", ...
    "The V17A source group declares no dynamic asset fields.");
records = repmat(localRecord(), ...
    numel(pathFields)*numel(manifest.sources)+1,1);
records(1) = localRecord(string(group.manifestPath), ...
    string(group.manifestSha256),"v17a_source_registry");
cursor = 1;
for index = 1:numel(manifest.sources)
    source = manifest.sources(index);
    for fieldIndex = 1:numel(pathFields)
        field = pathFields(fieldIndex);
        assert(isfield(source,field) && ...
            all(isfield(source.(field),["path","sha256"])), ...
            "AeroFlex:ReleaseAssetsSourceField", ...
            "Source %s has no complete %s asset record.", ...
            string(source.sourceId),field);
        cursor = cursor+1;
        records(cursor) = localRecord( ...
            string(source.(field).path), ...
            string(source.(field).sha256), ...
            "v17a_source_"+field+":"+string(source.sourceId));
    end
end
end

function records = localRecordGroup(root,group,pathField,hashField)
manifestPath = fullfile(root,localNativePath(string(group.manifestPath)));
assert(localFileHash(manifestPath)==string(group.manifestSha256), ...
    "AeroFlex:ReleaseAssetsGroupManifest", ...
    "A release-asset group manifest changed: %s",manifestPath);
manifest = jsondecode(fileread(manifestPath));
values = manifest.(char(group.recordField));
assert(numel(values)==double(group.expectedSourceCount), ...
    "AeroFlex:ReleaseAssetsGroupCount", ...
    "A release-asset group count differs from its contract.");
records = repmat(localRecord(),numel(values)+1,1);
records(1) = localRecord(string(group.manifestPath), ...
    string(group.manifestSha256),"transitive_manifest");
for index = 1:numel(values)
    records(index+1) = localRecord(string(values(index).(pathField)), ...
        string(values(index).(hashField)), ...
        "transitive_asset:"+string(values(index).sourceId));
end
end

function record = localRecord(path,hash,role)
if nargin==0
    path = ""; hash = ""; role = "";
end
record = struct("path",replace(string(path),"\","/"), ...
    "expectedSha256",lower(string(hash)),"role",string(role), ...
    "exists",false,"observedSha256","","passed",false);
end

function records = localUniqueRecords(records)
paths = string({records.path}).';
[uniquePaths,first] = unique(paths,"stable");
for path = uniquePaths.'
    indices = paths==path;
    hashes = unique(string({records(indices).expectedSha256}));
    assert(isscalar(hashes),"AeroFlex:ReleaseAssetsDuplicateHash", ...
        "A selected asset has conflicting expected hashes: %s",path);
end
records = records(first);
end

function destination = localValidateDestination(root,value)
assert(strlength(value)>0,"AeroFlex:ReleaseAssetsDestination", ...
    "DestinationRoot is required when Action is stage.");
destination = string(value);
if ~isfolder(destination), mkdir(destination); end
sourceCanonical = string(java.io.File(char(root)).getCanonicalPath());
destinationCanonical = string(java.io.File(char(destination)).getCanonicalPath());
assert(destinationCanonical~=sourceCanonical && ...
    ~startsWith(sourceCanonical,destinationCanonical+filesep), ...
    "AeroFlex:ReleaseAssetsDestinationScope", ...
    "The staging destination must not contain or equal the source repository.");
destination = destinationCanonical;
end

function localWriteInventory(destination,manifestPath,records)
inventory = struct("schemaVersion","pazy-benchmark-asset-inventory-v1", ...
    "generatedUtc",string(datetime("now","TimeZone","UTC", ...
        "Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
    "selectionManifestSha256",localFileHash(manifestPath), ...
    "assetCount",numel(records),"records",records);
path = fullfile(destination,"BENCHMARK_ASSET_INVENTORY.json");
file = fopen(path,"w");
assert(file>=0,"AeroFlex:ReleaseAssetsInventoryWrite", ...
    "Cannot write staged asset inventory: %s",path);
cleanup = onCleanup(@()fclose(file));
fprintf(file,"%s\n",jsonencode(inventory,"PrettyPrint",true));
clear cleanup
end

function bytes = localTotalBytes(root,records)
bytes = 0;
for index = 1:numel(records)
    info = dir(fullfile(root,localNativePath(records(index).path)));
    bytes = bytes+info.bytes;
end
end

function value = localNativePath(value)
value = replace(string(value),["/","\"],filesep);
end

function value = localFileHash(path)
file = fopen(path,"rb");
assert(file>=0,"AeroFlex:ReleaseAssetsHashRead", ...
    "Cannot read selected asset: %s",path);
cleanup = onCleanup(@()fclose(file));
engine = javaMethod("getInstance","java.security.MessageDigest","SHA-256");
while ~feof(file)
    bytes = fread(file,1024*1024,"*uint8");
    if isempty(bytes), break, end
    engine.update(typecast(bytes(:),"int8"));
end
value = lower(string(reshape(dec2hex( ...
    typecast(engine.digest(),"uint8"),2).',1,[])));
clear cleanup
end
