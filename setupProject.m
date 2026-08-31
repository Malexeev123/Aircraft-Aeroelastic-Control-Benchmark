function projectInfo = setupProject(options)
%SETUPPROJECT Configure the AeroFlex repository without a MATLAB Project.
%
%   projectInfo = setupProject()
%   projectInfo = setupProject(ValidateEntryPoints=true)
%   projectInfo = setupProject(ChangeCurrentFolder=true)
%
% This function configures the MATLAB environment required to inspect,
% analyze, and execute the AeroFlex source tree. It deliberately does not:
%
%   - require a MATLAB Project to be open;
%   - open a .prj file;
%   - modify MATLAB Project membership;
%   - recursively add package folders;
%   - run trim or simulation;
%   - run dependency analysis;
%   - clear the workspace;
%   - close figures.
%
% Only the parent folder of MATLAB namespaces is added to the path:
%
%   <repository>/MatlabFlex
%
% This correctly resolves namespaces such as:
%
%   AeroFlex.sim.ROMIntegrator
%   RigidBody.methods.paramsRigid_PazyUAV
%
% without directly adding +AeroFlex, +RigidBody, or their child namespaces.
    arguments
        options.ValidateEntryPoints (1,1) logical = true
        options.ChangeCurrentFolder (1,1) logical = false
    end

    %% Determine repository root from this file
    setupFile = string(mfilename("fullpath"));

    assert(strlength(setupFile) > 0, ...
        "RepositorySetup:UnknownSetupLocation", ...
        "MATLAB could not determine the location of setupProject.m.");

    setupFolder = string(fileparts(setupFile));
    repositoryRoot = localFindRepositoryRoot(setupFolder);

    matlabRoot = fullfile(repositoryRoot, "MatlabFlex");
    matlabConfigsRoot = fullfile(matlabRoot, "configs");
    matlabToolsRoot = fullfile(repositoryRoot, "tools", "matlab");

    assert(isfolder(matlabRoot), ...
        "RepositorySetup:MissingMatlabFlex", ...
        "The expected MatlabFlex folder was not found: %s", ...
        matlabRoot);

    %% Add only ordinary path roots
    %
    % Do not use genpath(MatlabFlex). It would independently add namespace
    % folders such as +AeroFlex and +RigidBody, which is incorrect.

    localRemoveRepositoryPaths(repositoryRoot);

    localAddPathIfNeeded(repositoryRoot, "-end");
    localAddPathIfNeeded(matlabConfigsRoot, "-end");

    if isfolder(matlabToolsRoot)
        localAddPathIfNeeded(matlabToolsRoot, "-end");
    end

    localAddPathIfNeeded(matlabRoot, "-begin");

    compiledControl = localConfigureCompiledControlCache( ...
        matlabRoot,string(version("-release")),string(computer("arch")));
    scheduledCompiledControl = localConfigureScheduledCompiledControlCache( ...
        matlabRoot,string(version("-release")),string(computer("arch")));

    rehash path;

    %% Optionally make repository root the current folder

    if options.ChangeCurrentFolder
        cd(repositoryRoot);
    end

    %% Inspect important entry points without executing them

    resolvableEntryPoints = [
        "nominalConfig"
        "sim_init"
        "sim_run"
        "RigidBody.methods.paramsRigid_PazyUAV"
    ];

    resolvedEntryPoints = strings(numel(resolvableEntryPoints), 1);
    missingEntryPoints = strings(0, 1);

    for entryIndex = 1:numel(resolvableEntryPoints)
        entryPoint = resolvableEntryPoints(entryIndex);
        resolvedPath = string(which(entryPoint));

        resolvedEntryPoints(entryIndex) = resolvedPath;

        if strlength(resolvedPath) == 0
            missingEntryPoints(end + 1, 1) = entryPoint; %#ok<AGROW>
        end
    end

    expectedEntryPoints = [
        fullfile(matlabConfigsRoot, "nominalConfig.m")
        fullfile(matlabRoot, "sim_init.m")
        fullfile(matlabRoot, "sim_run.m")
        fullfile(matlabRoot, "+RigidBody", "+methods", "paramsRigid_PazyUAV.m")
    ];
    pathAuthorityMessages = strings(0, 1);

    for entryIndex = 1:numel(resolvableEntryPoints)
        allMatches = string(which(resolvableEntryPoints(entryIndex), "-all"));
        allMatches(allMatches == "") = [];

        if numel(allMatches) ~= 1 || ...
                localNormalizePath(allMatches(1)) ~= ...
                localNormalizePath(expectedEntryPoints(entryIndex))
            pathAuthorityMessages(end + 1, 1) = sprintf( ... %#ok<AGROW>
                '%s must resolve uniquely to %s; found: %s', ...
                resolvableEntryPoints(entryIndex), ...
                expectedEntryPoints(entryIndex), ...
                strjoin(allMatches, ', '));
        end
    end

    %% Locate TrimRBwFlex independently of namespace qualification
    %
    % TrimRBwFlex may be:
    %   - an ordinary unqualified function;
    %   - inside +CoupledTrim;
    %   - inside a deeper RigidBody namespace.
    %
    % Context generation only needs to know where the source file exists.
    % It does not require that "TrimRBwFlex" resolve as an unqualified name.

    trimListing = dir(fullfile( ...
        matlabRoot, ...
        "**", ...
        "TrimRBwFlex.m"));

    trimFiles = strings(numel(trimListing), 1);

    for trimIndex = 1:numel(trimListing)
        trimFiles(trimIndex) = string(fullfile( ...
            trimListing(trimIndex).folder, ...
            trimListing(trimIndex).name));
    end

    %% Optional strict validation

    if options.ValidateEntryPoints
        validationMessages = strings(0, 1);

        if ~isempty(missingEntryPoints)
            validationMessages(end + 1, 1) = ...
                "Unresolvable MATLAB entry points:"; %#ok<AGROW>

            for entryIndex = 1:numel(missingEntryPoints)
                validationMessages(end + 1, 1) = ...
                    "  " + missingEntryPoints(entryIndex); %#ok<AGROW>
            end
        end

        if isempty(trimFiles)
            validationMessages(end + 1, 1) = ...
                "TrimRBwFlex.m was not found beneath MatlabFlex."; %#ok<AGROW>
        elseif numel(trimFiles) > 1
            validationMessages(end + 1, 1) = ...
                "Multiple TrimRBwFlex.m files were found:"; %#ok<AGROW>

            for trimIndex = 1:numel(trimFiles)
                validationMessages(end + 1, 1) = ...
                    "  " + trimFiles(trimIndex); %#ok<AGROW>
            end
        end

        validationMessages = [validationMessages; pathAuthorityMessages];

        if ~isempty(validationMessages)
            error("RepositorySetup:ValidationFailed", ...
                "%s", ...
                strjoin(validationMessages, newline));
        end
    end

    %% Read repository revision

    gitRevision = localGitRevision(repositoryRoot);

    %% Build returned information

    projectInfo = struct();

    projectInfo.repositoryRoot = repositoryRoot;
    projectInfo.matlabRoot = matlabRoot;
    projectInfo.matlabConfigsRoot = matlabConfigsRoot;
    projectInfo.matlabToolsRoot = matlabToolsRoot;

    projectInfo.matlabRelease = string(version("-release"));
    projectInfo.matlabVersion = string(version);
    projectInfo.computer = string(computer);
    projectInfo.timestamp = datetime("now");

    projectInfo.gitRevision = gitRevision;

    projectInfo.entryPointNames = resolvableEntryPoints;
    projectInfo.resolvedEntryPoints = resolvedEntryPoints;
    projectInfo.missingEntryPoints = missingEntryPoints;
    projectInfo.trimFiles = trimFiles;

    projectInfo.matlabProjectRequired = false;
    projectInfo.matlabProjectOpen = ~isempty(matlab.project.rootProject);
    projectInfo.compiledControl = compiledControl;
    projectInfo.scheduledCompiledControl = scheduledCompiledControl;

    %% Console summary

    fprintf("\n");
    fprintf("AeroFlex repository environment configured.\n");
    fprintf("  Repository root : %s\n", repositoryRoot);
    fprintf("  MATLAB root     : %s\n", matlabRoot);
    fprintf("  MATLAB release  : %s\n", projectInfo.matlabRelease);
    fprintf("  Git revision    : %s\n", gitRevision);
    fprintf("  Project required: no\n");
    fprintf("  Compiled control: %s\n",compiledControl.message);
    fprintf("  Scheduled compiled control: %s\n", ...
        scheduledCompiledControl.message);

    fprintf("\nEntry-point resolution:\n");

    for entryIndex = 1:numel(resolvableEntryPoints)
        if strlength(resolvedEntryPoints(entryIndex)) == 0
            fprintf("  %-45s : NOT FOUND\n", ...
                resolvableEntryPoints(entryIndex));
        else
            fprintf("  %-45s : %s\n", ...
                resolvableEntryPoints(entryIndex), ...
                resolvedEntryPoints(entryIndex));
        end
    end

    fprintf("\nTrimRBwFlex source files found: %d\n", numel(trimFiles));

    for trimIndex = 1:numel(trimFiles)
        fprintf("  %s\n", trimFiles(trimIndex));
    end

    if ~isempty(missingEntryPoints)
        fprintf(2, ...
            "\nWarning: %d optional entry points are currently unresolved.\n", ...
            numel(missingEntryPoints));

        fprintf(2, ...
            "Context generation can still proceed.\n");
    end
end

function status = localConfigureCompiledControlCache( ...
        matlabRoot,releaseName,architecture)
%LOCALCONFIGURECOMPILEDCONTROLCACHE Select only a compatible verified cache.
    status = struct('active',false,'path',"", ...
        'message',"BUILD REQUIRED (run buildPhase18CompiledControlKernels)");
    sourcePaths = [ ...
        string(fullfile(matlabRoot,'+AeroFlex','+ctrl', ...
            'fixedReciprocalIntervalKernelAudit.m')); ...
        string(fullfile(matlabRoot,'+AeroFlex','+ctrl', ...
            'buildFixedReciprocalIntervalPacketAudit.m')); ...
        string(fullfile(matlabRoot,'configs','benchmark', ...
            'native-build-fixtures', ...
            'phase18c_v17a_fixed_interval_build_fixture_v1.mat'))];
    if ~all(isfile(sourcePaths))
        status.message = ...
            "FALLBACK ACTIVE; compiled-kernel source files are unavailable";
        localPrintCompiledControlStatus(status.message);
        return
    end
    sourceHashes = strings(size(sourcePaths));
    for sourceIndex = 1:numel(sourcePaths)
        sourceHashes(sourceIndex) = localSetupFileHash(sourcePaths(sourceIndex));
    end
    currentSourceSignature = localSetupTextHash(strjoin(sourceHashes,"|"));
    cacheRoot = fullfile(matlabRoot,'cache','compiled_control', ...
        releaseName,architecture);
    if ~isfolder(cacheRoot)
        localPrintCompiledControlStatus(status.message);
        return
    end
    manifests = dir(fullfile(cacheRoot,'**', ...
        'PHASE18_COMPILED_CONTROL_MANIFEST.json'));
    for index = numel(manifests):-1:1
        manifestPath = fullfile(manifests(index).folder,manifests(index).name);
        try
            manifest = jsondecode(fileread(manifestPath));
            binaryPath = fullfile(fileparts(manifestPath), ...
                ['AeroFlex_ctrl_fixedReciprocalIntervalKernelProduction_mex.',mexext]);
            compatible = logical(manifest.passed) && ...
                string(manifest.matlabRelease)==releaseName && ...
                string(manifest.architecture)==architecture && ...
                string(manifest.mexExtension)==string(mexext) && ...
                isfield(manifest,'sourceSignature') && ...
                string(manifest.sourceSignature)==currentSourceSignature && ...
                isfield(manifest,'functionName') && ...
                string(manifest.functionName)== ...
                    "AeroFlex_ctrl_fixedReciprocalIntervalKernelProduction_mex" && ...
                isfile(binaryPath) && string(manifest.binarySha256)== ...
                    localSetupFileHash(binaryPath);
            if compatible
                localAddPathIfNeeded(fileparts(manifestPath),'-begin');
                status.active = true;
                status.path = string(fileparts(manifestPath));
                status.message = "ACTIVE (verified local MEX cache)";
                localPrintCompiledControlStatus(status.message);
                return
            end
        catch
            % Ignore incomplete or stale cache entries and keep searching.
        end
    end
    status.message = ...
        "FALLBACK ACTIVE; cached binary is stale or incompatible (rebuild recommended)";
    localPrintCompiledControlStatus(status.message);
end

function localPrintCompiledControlStatus(message)
    persistent previousMessage
    if isempty(previousMessage) || string(previousMessage)~=string(message)
        fprintf('Compiled control acceleration: %s\n',message);
        previousMessage = string(message);
    end
end

function status = localConfigureScheduledCompiledControlCache( ...
        matlabRoot,releaseName,architecture)
%LOCALCONFIGURESCHEDULEDCOMPILEDCONTROLCACHE Select verified scheduled MEX.
    status = struct('active',false,'path',"", ...
        'message', ...
        "BUILD REQUIRED (run buildPhase18ScheduledCompiledControlKernel)");
    sourcePaths = [ ...
        string(fullfile(matlabRoot,'+AeroFlex','+ctrl', ...
            'scheduledReciprocalIntervalKernelAudit.m')); ...
        string(fullfile(matlabRoot,'+AeroFlex','+ctrl', ...
            'buildScheduledReciprocalIntervalPacketAudit.m')); ...
        string(fullfile(matlabRoot,'configs','benchmark', ...
            'native-build-fixtures', ...
            'phase18c_v17a_scheduled_build_fixture_v1.mat'))];
    if ~all(isfile(sourcePaths))
        status.message = ...
            "FALLBACK ACTIVE; scheduled compiled sources are unavailable";
        return
    end
    sourceHashes = strings(size(sourcePaths));
    for sourceIndex = 1:numel(sourcePaths)
        sourceHashes(sourceIndex) = localSetupFileHash(sourcePaths(sourceIndex));
    end
    currentSourceSignature = localSetupTextHash(strjoin(sourceHashes,"|"));
    cacheRoot = fullfile(matlabRoot,'cache','compiled_control_scheduled', ...
        releaseName,architecture);
    if ~isfolder(cacheRoot), return, end
    manifests = dir(fullfile(cacheRoot,'**', ...
        'PHASE18_SCHEDULED_COMPILED_CONTROL_MANIFEST.json'));
    functionName = ...
        "AeroFlex_ctrl_scheduledReciprocalIntervalKernelAudit_mex";
    for index = numel(manifests):-1:1
        manifestPath = fullfile(manifests(index).folder,manifests(index).name);
        try
            manifest = jsondecode(fileread(manifestPath));
            binaryPath = fullfile(fileparts(manifestPath), ...
                char(functionName+"."+mexext));
            compatible = logical(manifest.passed) && ...
                string(manifest.matlabRelease)==releaseName && ...
                string(manifest.architecture)==architecture && ...
                string(manifest.mexExtension)==string(mexext) && ...
                string(manifest.sourceSignature)==currentSourceSignature && ...
                string(manifest.functionName)==functionName && ...
                isfile(binaryPath) && string(manifest.binarySha256)== ...
                    localSetupFileHash(binaryPath);
            if compatible
                localAddPathIfNeeded(fileparts(manifestPath),'-begin');
                status.active = true;
                status.path = string(fileparts(manifestPath));
                status.message = "ACTIVE (verified local MEX cache)";
                return
            end
        catch
            % Ignore incomplete or stale cache entries and keep searching.
        end
    end
    status.message = ...
        "FALLBACK ACTIVE; scheduled cached binary is stale or incompatible";
end

function value = localSetupFileHash(pathIn)
    file = fopen(pathIn,'rb');
    assert(file>=0,'RepositorySetup:CompiledControlHashRead', ...
        'Cannot read %s for hashing.',pathIn);
    cleanup = onCleanup(@() fclose(file));
    bytes = fread(file,inf,'*uint8');
    engine = java.security.MessageDigest.getInstance('SHA-256');
    engine.update(bytes);
    clear cleanup
    digest = typecast(engine.digest(),'uint8');
    value = lower(string(sprintf('%02x',digest)));
end

function value = localSetupTextHash(text)
    engine = java.security.MessageDigest.getInstance('SHA-256');
    engine.update(uint8(char(text)));
    digest = typecast(engine.digest(),'uint8');
    value = lower(string(sprintf('%02x',digest)));
end

function repositoryRoot = localFindRepositoryRoot(startFolder)
%LOCALFINDREPOSITORYROOT Ascend until the Git repository root is found.

    candidate = string(startFolder);

    while true
        gitDirectory = fullfile(candidate, ".git");
        publicLauncher = fullfile(candidate, "Run_Pazy_Benchmark.m");
        matlabFolder = fullfile(candidate, "MatlabFlex");

        isRepository = ...
            isfolder(gitDirectory) || ...
            isfile(gitDirectory) || ...
            (isfile(publicLauncher) && isfolder(matlabFolder));

        if isRepository
            repositoryRoot = candidate;
            return;
        end

        parent = string(fileparts(candidate));

        if parent == candidate || strlength(parent) == 0
            break;
        end

        candidate = parent;
    end

    error("RepositorySetup:RepositoryRootNotFound", ...
        ["Could not locate the repository root while ascending from:" ...
         newline + startFolder]);
end

function localAddPathIfNeeded(folder, position)
%LOCALADDPATHIFNEEDED Add a folder only when it is not already on path.

    folder = string(folder);

    currentPathEntries = string(split(path, pathsep));
    currentPathEntries(currentPathEntries == "") = [];

    normalizedFolder = localNormalizePath(folder);
    normalizedEntries = arrayfun( ...
        @localNormalizePath, ...
        currentPathEntries);

    if ~any(normalizedEntries == normalizedFolder)
        addpath(folder, position);
    end
end

function localRemoveRepositoryPaths(repositoryRoot)
%LOCALREMOVEREPOSITORYPATHS Remove stale paths beneath this checkout.

    entries = string(split(path, pathsep));
    entries(entries == "") = [];
    repositoryRootNormalized = localNormalizePath(repositoryRoot);
    repositoryPrefix = repositoryRootNormalized + "/";

    for entryIndex = 1:numel(entries)
        normalizedEntry = localNormalizePath(entries(entryIndex));
        if normalizedEntry == repositoryRootNormalized || ...
                startsWith(normalizedEntry, repositoryPrefix)
            rmpath(entries(entryIndex));
        end
    end
end

function normalizedPath = localNormalizePath(inputPath)
%LOCALNORMALIZEPATH Normalize separators and case for comparison.

    normalizedPath = lower(replace(string(inputPath), "\", "/"));

    while endsWith(normalizedPath, "/")
        normalizedPath = extractBefore( ...
            normalizedPath, ...
            strlength(normalizedPath));
    end
end

function revision = localGitRevision(repositoryRoot)
%LOCALGITREVISION Read Git HEAD without invoking a shell from a UNC folder.

    revision = "unavailable";
    gitFolder = fullfile(repositoryRoot, ".git");
    headFile = fullfile(gitFolder, "HEAD");
    if ~isfile(headFile)
        return;
    end

    headValue = strtrim(string(fileread(headFile)));
    if startsWith(headValue, "ref: ")
        refName = extractAfter(headValue, "ref: ");
        refFile = fullfile(gitFolder, replace(refName, "/", filesep));
        if isfile(refFile)
            headValue = strtrim(string(fileread(refFile)));
        else
            packedRefsFile = fullfile(gitFolder, "packed-refs");
            if ~isfile(packedRefsFile)
                return;
            end
            packedLines = splitlines(string(fileread(packedRefsFile)));
            match = packedLines(endsWith(packedLines, " " + refName));
            if isempty(match)
                return;
            end
            headValue = extractBefore(match(1), " ");
        end
    end

    if ~isempty(regexp(headValue, "^[0-9a-fA-F]{40}$", "once"))
        revision = lower(headValue);
    end
end


% function projectInfo = setupProject()
% %SETUPPROJECT Configure and validate the Pazy aeroelastic-control project.
% %
% % projectInfo = setupProject()
% %
% % The setup is intentionally conservative. It does not clear the workspace,
% % close figures, reset persistent state, modify preferences, or run a
% % simulation.
% thisFile = mfilename("fullpath");
% repositoryRoot = fileparts(thisFile);
% if strlength(repositoryRoot) == 0 || ~isfolder(repositoryRoot)
%     error("ProjectSetup:InvalidRoot", ...
%         "Unable to determine a valid repository root.");
% end
% candidatePaths = [
%     string(repositoryRoot)
%     fullfile(repositoryRoot, "tests")
%     fullfile(repositoryRoot, "tests", "unit")
%     fullfile(repositoryRoot, "tests", "integration")
%     fullfile(repositoryRoot, "tests", "regression")
%     fullfile(repositoryRoot, "tools", "matlab")
%     ];
% addedPaths = strings(0, 1);
% for pathIndex = 1:numel(candidatePaths)
%     candidatePath = candidatePaths(pathIndex);
%     if isfolder(candidatePath)
%         addpath(candidatePath);
%         addedPaths(end + 1, 1) = candidatePath; %#ok<AGROW>
%     end
% end
% requiredEntryPoints = [
%     "sim_init"
%     "sim_run"
%     "RigidBody.methods.CoupledTrim.TrimRBwFlex"
%     "nominalConfig"
%     ];
% missingEntryPoints = strings(0, 1);
% for entryIndex = 1:numel(requiredEntryPoints)
%     entryName = requiredEntryPoints(entryIndex);
%     if isempty(which(entryName))
%         missingEntryPoints(end + 1, 1) = entryName; %#ok<AGROW>
%     end
% end
% if ~isempty(missingEntryPoints)
%     error("ProjectSetup:MissingEntryPoint", ...
%         "Required entry points are missing from the MATLAB path:%s%s", ...
%         newline, strjoin(missingEntryPoints, newline));
% end
% gitRevision = localGitRevision(repositoryRoot);
% projectInfo = struct( ...
%     "repositoryRoot", string(repositoryRoot), ...
%     "matlabRelease", string(version("-release")), ...
%     "matlabVersion", string(version), ...
%     "computer", string(computer), ...
%     "workingDirectory", string(pwd), ...
%     "addedPaths", addedPaths, ...
%     "gitRevision", gitRevision, ...
%     "timestamp", datetime("now"));
% fprintf("Pazy aeroelastic-control project configured.\n");
% fprintf(" Repository root : %s\n", projectInfo.repositoryRoot);
% fprintf(" MATLAB release : %s\n", projectInfo.matlabRelease);
% fprintf(" Git revision : %s\n", projectInfo.gitRevision);
% fprintf(" Added paths : %d\n", numel(projectInfo.addedPaths));
% end
% function revision = localGitRevision(repositoryRoot)
% %LOCALGITREVISION Return the current Git revision when Windows Git can
% %resolve the WSL project. Git metadata generation does not fail project setup.
% command = sprintf('git -C "%s" rev-parse HEAD', repositoryRoot);
% [status, output] = system(command);
% if status == 0
%     revision = strtrim(string(output));
% else
%     revision = "unavailable-use-node-context-generator";
% end
% end
