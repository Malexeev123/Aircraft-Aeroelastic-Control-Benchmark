function status = prepareBenchmarkReleasePackage(options)
%PREPAREBENCHMARKRELEASEPACKAGE Classify and stage a public source package.
%   STATUS = PREPAREBENCHMARKRELEASEPACKAGE() computes the executable MATLAB
%   dependency closure, verifies locked runtime assets, and reports the public
%   package boundary without copying files.
%
%   STATUS = PREPAREBENCHMARKRELEASEPACKAGE(Action="stage", ...
%   DestinationRoot=FOLDER) copies only the classified source, documentation,
%   examples, validation entries, and verified runtime assets. Existing files
%   are reused only when their bytes match; changed files are never overwritten.

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

entryPoints = localEntryPoints(root);
[dependencyFiles,products] = matlab.codetools.requiredFilesAndProducts( ...
    cellstr(entryPoints));
dependencyFiles = localRepositoryFiles(root,string(dependencyFiles(:)));
publicFiles = localUniqueFiles([dependencyFiles;localPublicFiles(root)]);
records = localRecords(root,publicFiles);

assets = prepareBenchmarkReleaseAssets(Action="check", ...
    ProjectInfo=project,PrintSummary=false);
licenseFiles = localLicenseFiles(root);
licensePresent = ~isempty(licenseFiles);

destination = "";
stagedCount = 0;
reusedCount = 0;
if options.Action=="stage"
    assert(licensePresent,"AeroFlex:ReleasePackageLicense", ...
        ["A public package requires a reviewed repository license. ", ...
         "Add the selected LICENSE file before staging."]);
    destination = localValidateDestination(root,options.DestinationRoot);
    [stagedCount,reusedCount] = localStageRecords(root,destination,records);
    stagedAssets = prepareBenchmarkReleaseAssets(Action="stage", ...
        DestinationRoot=destination,ProjectInfo=project,PrintSummary=false);
    localWriteInventory(destination,records,stagedAssets,products);
end

status = struct( ...
    "schemaVersion","pazy-public-release-package-status-v1", ...
    "action",options.Action, ...
    "passed",all([records.passed]) && assets.passed, ...
    "publicationReady",all([records.passed]) && assets.passed && ...
        licensePresent, ...
    "sourceFileCount",numel(records), ...
    "runtimeAssetCount",assets.assetCount, ...
    "dependencyProductNames",string({products.Name}).', ...
    "licensePresent",licensePresent, ...
    "licenseFiles",licenseFiles, ...
    "stagedCount",stagedCount, ...
    "reusedCount",reusedCount, ...
    "destinationRoot",destination, ...
    "records",records, ...
    "excludedCategories",localExcludedCategories());

if options.PrintSummary
    fprintf("\nPazy public release package\n");
    fprintf("  Required source files : %d\n",status.sourceFileCount);
    fprintf("  Verified runtime data : %d\n",status.runtimeAssetCount);
    fprintf("  Repository license    : %s\n",localYesNo(licensePresent));
    fprintf("  Source/data integrity : %s\n",localPassFail(status.passed));
    fprintf("  Publication ready     : %s\n", ...
        localPassFail(status.publicationReady));
    if options.Action=="stage"
        fprintf("  Newly staged          : %d\n",stagedCount);
        fprintf("  Reused                : %d\n",reusedCount);
        fprintf("  Destination           : %s\n",destination);
    end
    fprintf("\n");
end
end

function files = localEntryPoints(root)
files = [ ...
    fullfile(root,"Run_Pazy_Benchmark.m")
    fullfile(root,"setupProject.m")
    fullfile(root,"MatlabFlex","sim_init.m")
    fullfile(root,"MatlabFlex","sim_run.m")
    fullfile(root,"MatlabFlex","runBenchmarkCase.m")
    fullfile(root,"MatlabFlex","runPazyModelWorkflow.m")
    fullfile(root,"MatlabFlex","verifyBenchmarkInstallation.m")
    fullfile(root,"MatlabFlex","+AeroFlex","+benchmark","+runtime", ...
        "run_phase18c_v17a_formal_casea_freeflight_attitude_hold_v1.m")
    fullfile(root,"MatlabFlex","+AeroFlex","+benchmark","+runtime", ...
        "run_phase18c_v17a_caseb_integrated_profile_v1.m")
    fullfile(root,"tools","matlab","buildBenchmarkTools.m")
    fullfile(root,"tools","matlab", ...
        "buildPhase18CompiledControlKernels.m")
    fullfile(root,"tools","matlab", ...
        "buildPhase18ScheduledCompiledControlKernel.m")
    fullfile(root,"tools","matlab", ...
        "buildPhase18ScheduledCompiledHorizonKernels.m")
    fullfile(root,"tools","matlab", ...
        "buildPhase18ScheduledCompiledValueHorizonKernels.m")
    fullfile(root,"tools","matlab", ...
        "buildPhase18ScheduledCompiledCausalRolloutKernels.m")
    fullfile(root,"tools","matlab","prepareBenchmarkReleaseAssets.m")
    fullfile(root,"tools","matlab","prepareBenchmarkReleasePackage.m")];
kernelNames = [ ...
    "fixedReciprocalIntervalKernelAudit.m"
    "buildFixedReciprocalIntervalPacketAudit.m"
    "scheduledReciprocalIntervalKernelAudit.m"
    "buildScheduledReciprocalIntervalPacketAudit.m"
    "buildScheduledReciprocalHorizonPacketAudit.m"
    "scheduledReciprocalHorizonKernelCoreAudit.m"
    "scheduledReciprocalEstimatorHorizonKernelAudit.m"
    "scheduledReciprocalControllerHorizonKernelAudit.m"
    "scheduledReciprocalReducedTangentHorizonCoreAudit.m"
    "scheduledReciprocalEstimatorReducedTangentHorizonAudit.m"
    "scheduledReciprocalControllerReducedTangentHorizonAudit.m"
    "scheduledReciprocalValueHorizonCoreAudit.m"
    "scheduledReciprocalEstimatorValueHorizonAudit.m"
    "scheduledReciprocalControllerValueHorizonAudit.m"
    "scheduledReciprocalCausalRolloutCoreAudit.m"
    "scheduledReciprocalEstimatorCausalRolloutAudit.m"
    "scheduledReciprocalControllerCausalRolloutAudit.m"];
files = [files;fullfile(root,"MatlabFlex","+AeroFlex","+ctrl", ...
    kernelNames)];
assert(all(isfile(files)),"AeroFlex:ReleasePackageEntryPoint", ...
    "One or more declared public entry points are missing.");
end

function files = localPublicFiles(root)
files = strings(0,1);
exact = [ ...
    "README.md"
    "AUTHORS.md"
    "CITATION.cff"
    "THIRD_PARTY_NOTICES.md"
    "tools/release/public.gitattributes"
    "tools/release/public.gitignore"
    "tests/README.md"
    "tests/Run_Linear_Validation.m"
    "tests/Run_Extended_Validation.m"
    "tests/test_native_tools.m"
    "MatlabFlex/configs/benchmark/phase18c_v17a_release_asset_manifest_v1.json"
    "MatlabFlex/configs/benchmark/native-build-fixtures/README.md"];
for path = exact.'
    candidate = fullfile(root,localNativePath(path));
    if isfile(candidate), files(end+1,1) = candidate; end %#ok<AGROW>
end

patterns = ["docs/*.md","references/*.bib","references/*.json", ...
    "references/*.md","MatlabFlex/configs/*.m", ...
    "MatlabFlex/configs/benchmark/*.m", ...
    "MatlabFlex/configs/benchmark/*.json"];
for pattern = patterns
    listing = dir(fullfile(root,localNativePath(pattern)));
    listing = listing(~[listing.isdir]);
    files = [files;string(fullfile({listing.folder},{listing.name})).']; ...
        %#ok<AGROW>
end

licenseListing = [dir(fullfile(root,"LICENSE*")); ...
    dir(fullfile(root,"COPYING*"))];
licenseListing = licenseListing(~[licenseListing.isdir]);
files = [files;string(fullfile( ...
    {licenseListing.folder},{licenseListing.name})).'];
end

function files = localLicenseFiles(root)
listing = [dir(fullfile(root,"LICENSE*"));dir(fullfile(root,"COPYING*"))];
listing = listing(~[listing.isdir]);
files = string({listing.name}).';
end

function files = localRepositoryFiles(root,files)
prefix = lower(localNormalize(root)+"/");
normalized = lower(localNormalize(files));
files = files(startsWith(normalized,prefix));
end

function files = localUniqueFiles(files)
files = files(strlength(files)>0 & isfile(files));
[~,first] = unique(lower(localNormalize(files)),"stable");
files = files(first);
end

function records = localRecords(root,files)
records = repmat(struct("path","","sha256","", ...
    "bytes",0,"passed",false),numel(files),1);
for index = 1:numel(files)
    relative = extractAfter(localNormalize(files(index)), ...
        strlength(localNormalize(root))+1);
    if relative=="tools/release/public.gitignore"
        relative = ".gitignore";
    elseif relative=="tools/release/public.gitattributes"
        relative = ".gitattributes";
    end
    info = dir(files(index));
    records(index) = struct("path",relative, ...
        "sha256",localFileHash(files(index)),"bytes",info.bytes, ...
        "passed",true);
end
[~,order] = sort(lower(string({records.path})));
records = records(order);
end

function [staged,reused] = localStageRecords(root,destination,records)
staged = 0;
reused = 0;
for index = 1:numel(records)
    source = localRecordSource(root,records(index).path);
    target = fullfile(destination,localNativePath(records(index).path));
    folder = string(fileparts(target));
    if ~isfolder(folder), mkdir(folder); end
    if isfile(target)
        assert(localFileHash(target)==records(index).sha256, ...
            "AeroFlex:ReleasePackageExistingMismatch", ...
            "Refusing to overwrite a changed staged file: %s",target);
        reused = reused+1;
        continue
    end
    [copied,message] = copyfile(source,target);
    assert(copied,"AeroFlex:ReleasePackageCopy", ...
        "Unable to stage %s: %s",source,message);
    assert(localFileHash(target)==records(index).sha256, ...
        "AeroFlex:ReleasePackagePostCopyHash", ...
        "The staged file hash changed during copying: %s",target);
    staged = staged+1;
end
end

function source = localRecordSource(root,publicPath)
if publicPath==".gitignore"
    source = fullfile(root,"tools","release","public.gitignore");
elseif publicPath==".gitattributes"
    source = fullfile(root,"tools","release","public.gitattributes");
else
    source = fullfile(root,localNativePath(publicPath));
end
end

function destination = localValidateDestination(root,value)
assert(strlength(value)>0,"AeroFlex:ReleasePackageDestination", ...
    "DestinationRoot is required when Action is stage.");
destination = string(value);
if ~isfolder(destination), mkdir(destination); end
sourceCanonical = string(java.io.File(char(root)).getCanonicalPath());
destinationCanonical = string(java.io.File(char(destination)).getCanonicalPath());
assert(destinationCanonical~=sourceCanonical && ...
    ~startsWith(sourceCanonical,destinationCanonical+filesep), ...
    "AeroFlex:ReleasePackageDestinationScope", ...
    "The staging destination must not contain or equal the source repository.");
destination = destinationCanonical;
end

function localWriteInventory(destination,records,assets,products)
inventory = struct( ...
    "schemaVersion","pazy-public-release-inventory-v1", ...
    "generatedUtc",string(datetime("now","TimeZone","UTC", ...
        "Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'")), ...
    "sourceFiles",records, ...
    "runtimeAssets",assets.records, ...
    "matlabProducts",string({products.Name}).', ...
    "excludedCategories",localExcludedCategories());
path = fullfile(destination,"PUBLIC_RELEASE_INVENTORY.json");
file = fopen(path,"w");
assert(file>=0,"AeroFlex:ReleasePackageInventoryWrite", ...
    "Cannot write public release inventory: %s",path);
cleanup = onCleanup(@()fclose(file));
fprintf(file,"%s\n",jsonencode(inventory,"PrettyPrint",true));
clear cleanup
end

function values = localExcludedCategories()
values = [ ...
    "private development context and conversations"
    "campaign audit scripts, checkpoints, logs, and rejected candidates"
    "ordinary simulation histories and local result directories"
    "MATLAB MEX binaries, codegen products, and first-call caches"
    "Python bytecode, virtual environments, and test caches"
    "editor, operating-system, credential, and machine-local state"
    "publisher PDFs without verified redistribution permission"
    "obsolete standalone launchers outside the documented interface"];
end

function value = localFileHash(path)
file = fopen(path,"rb");
assert(file>=0,"AeroFlex:ReleasePackageHashRead", ...
    "Cannot read selected release file: %s",path);
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

function value = localNormalize(value)
value = replace(string(value),"\","/");
end

function value = localNativePath(value)
value = replace(string(value),["/","\"],filesep);
end

function value = localPassFail(flag)
if flag, value = "PASS"; else, value = "NOT READY"; end
end

function value = localYesNo(flag)
if flag, value = "present"; else, value = "not selected"; end
end
