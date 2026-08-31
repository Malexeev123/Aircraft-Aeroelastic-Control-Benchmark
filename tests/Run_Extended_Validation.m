%% Pazy benchmark: extended source validation
% Validate the accepted V17 source packages directly. The suite checks source
% hashes, trim ownership, nonlinear residuals, and flexible/aerodynamic poles,
% then saves a checkpoint after every source.

if ~exist("validationSettings","var") || ~isstruct(validationSettings)
    validationSettings = struct();
end

if ~isfield(validationSettings,"mode")
    validationSettings.mode = "preflight";  % Set to "full" for execution.
end
if ~isfield(validationSettings,"sourceIndices")
    validationSettings.sourceIndices = 1:29;
end

repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
setupProject(ValidateEntryPoints=true,ChangeCurrentFolder=false);

registryPath = fullfile(repositoryRoot,"context","audits", ...
    "phase18b-kkt-path-closure","static-feedback-source-placement-v1", ...
    "V17_PRODUCTION_SOURCE_REGISTRY_V1.json");
expectedRegistryHash = ...
    "2702cd0dcdae15c13f6bcbc2c59d8b4f899054c63a5df607dcacd0d282ac473f";

assert(isfile(registryPath),"PazyValidation:MissingRegistry", ...
    "The locked V17 source registry is unavailable.");
assert(fileSha256(registryPath) == expectedRegistryHash, ...
    "PazyValidation:RegistryHashMismatch", ...
    "The V17 registry hash does not match the locked production value.");

registry = jsondecode(fileread(registryPath));
assert(registry.sourceCount == 29 && numel(registry.sources) == 29, ...
    "PazyValidation:RegistryCountMismatch", ...
    "The V17 registry must contain exactly 29 accepted sources.");

sourceIndices = double(validationSettings.sourceIndices(:).');
assert(all(isfinite(sourceIndices)) && ...
       all(sourceIndices == fix(sourceIndices)) && ...
       all(sourceIndices >= 1 & sourceIndices <= registry.sourceCount), ...
    "PazyValidation:InvalidSourceIndex", ...
    "Source indices must be integers in the range 1:29.");

outputRoot = fullfile(repositoryRoot,"results","validation", ...
    "extended-source-linearization");
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

checkpointPath = fullfile(outputRoot,"source_validation_checkpoint.mat");
jsonPath = fullfile(outputRoot,"source_validation_summary.json");
plotPath = fullfile(outputRoot,"source_validation_trends.png");

fprintf("\nPazy Aeroelastic Control Benchmark\n");
fprintf("Extended V17 source trim and linearization validation\n");
fprintf("  Mode       : %s\n",validationSettings.mode);
fprintf("  Sources    : %s\n",mat2str(sourceIndices));
fprintf("  Runtime    : typically under two minutes after setup\n");
fprintf("  Checkpoint : %s\n",checkpointPath);

switch lower(string(validationSettings.mode))
    case "preflight"
        extendedValidationReport = struct( ...
            "schemaVersion",1, ...
            "status","PREFLIGHT_PASS", ...
            "mode","preflight", ...
            "sourceIndices",sourceIndices, ...
            "registryPath",string(registryPath), ...
            "registrySha256",expectedRegistryHash, ...
            "outputRoot",string(outputRoot));
        fprintf("  Status     : PREFLIGHT PASS (numerical stages not executed)\n\n");

    case "full"
        extendedValidationReport = initializeReport( ...
            registry,sourceIndices,registryPath,expectedRegistryHash);
        suiteTimer = tic;

        for localIndex = 1:numel(sourceIndices)
            sourceIndex = sourceIndices(localIndex);
            sourceTimer = tic;
            sourceRecord = registry.sources(sourceIndex);

            try
                row = analyzeSource(repositoryRoot,sourceRecord);
                row.elapsedSeconds = toc(sourceTimer);
                fprintf("  %2d/%2d %-12s PASS  residual=%9.3e maxRe=%+9.3e\n", ...
                    localIndex,numel(sourceIndices),row.sourceId, ...
                    row.propagatedResidualNorm,row.maximumPoleRealPerSecond);
            catch exception
                row = emptyRow();
                row.sourceId = string(sourceRecord.sourceId);
                row.coordinates = double(sourceRecord.coordinates(:).');
                row.status = "FAIL";
                row.message = string(exception.message);
                row.elapsedSeconds = toc(sourceTimer);
                fprintf("  %2d/%2d %-12s FAIL  %s\n", ...
                    localIndex,numel(sourceIndices),row.sourceId,row.message);
            end

            extendedValidationReport.rows(localIndex) = row;
            extendedValidationReport.completed(localIndex) = true;
            extendedValidationReport.updated = string(datetime("now"));
            save(checkpointPath,"extendedValidationReport","-v7");
        end

        extendedValidationReport.totalElapsedSeconds = toc(suiteTimer);
        extendedValidationReport.allPassed = ...
            all(extendedValidationReport.completed) && ...
            all([extendedValidationReport.rows.status] == "PASS");
        extendedValidationReport.status = passFail( ...
            extendedValidationReport.allPassed);

        save(checkpointPath,"extendedValidationReport","-v7");
        writeSummary(jsonPath,extendedValidationReport);
        figureHandle = plotSummary(extendedValidationReport.rows);
        exportgraphics(figureHandle,plotPath,"Resolution",180);
        close(figureHandle);

        extendedValidationReport.checkpointPath = string(checkpointPath);
        extendedValidationReport.jsonPath = string(jsonPath);
        extendedValidationReport.plotPath = string(plotPath);

        fprintf("  Elapsed    : %.3f s\n", ...
            extendedValidationReport.totalElapsedSeconds);
        fprintf("  Status     : %s\n",extendedValidationReport.status);
        fprintf("  Plot       : %s\n\n",plotPath);

    otherwise
        error("PazyValidation:UnknownExtendedMode", ...
            "validationSettings.mode must be 'preflight' or 'full'.");
end


function report = initializeReport( ...
    registry,sourceIndices,registryPath,registryHash)
%INITIALIZEREPORT Create a deterministic checkpoint-compatible result.

report = struct( ...
    "schemaVersion",1, ...
    "status","RUNNING", ...
    "mode","full", ...
    "sourceIndices",sourceIndices, ...
    "registryPath",string(registryPath), ...
    "registrySha256",registryHash, ...
    "registryVersion",string(registry.registryVersion), ...
    "completed",false(1,numel(sourceIndices)), ...
    "rows",repmat(emptyRow(),1,numel(sourceIndices)), ...
    "allPassed",false, ...
    "totalElapsedSeconds",NaN, ...
    "updated",string(datetime("now")));
end


function row = analyzeSource(repositoryRoot,sourceRecord)
%ANALYZESOURCE Evaluate one immutable accepted package at its own node.

relativePath = string(sourceRecord.acceptedMatlabPackage.path);
packagePath = fullfile(repositoryRoot,relativePath);
assert(isfile(packagePath),"PazyValidation:MissingSourcePackage", ...
    "Accepted source package not found: %s",packagePath);

expectedHash = lower(string(sourceRecord.acceptedMatlabPackage.sha256));
actualHash = fileSha256(packagePath);
assert(actualHash == expectedHash,"PazyValidation:SourceHashMismatch", ...
    "Accepted source package hash mismatch: %s",relativePath);

data = load(packagePath,"p5","trim");
assert(isfield(data,"p5") && isfield(data,"trim"), ...
    "PazyValidation:IncompleteSourcePackage", ...
    "The accepted source MAT file must contain p5 and trim.");

package = data.p5;
trim = data.trim;
requiredFields = ["L","Ldyn","x_eq","u_eq","idx","beam","parConst"];
assert(all(isfield(package,requiredFields)), ...
    "PazyValidation:IncompleteAcceptedPackage", ...
    "The accepted exact-source package is missing a required runtime field.");
assert(size(package.L,1) == numel(package.x_eq) && ...
       isequal(size(package.L),size(package.Ldyn)), ...
    "PazyValidation:InvalidAcceptedDimensions", ...
    "The accepted exact-source operators and equilibrium state disagree.");

assert(isfield(trim,"thrust") && isfinite(trim.thrust) && trim.thrust >= 0, ...
    "PazyValidation:NegativeThrust", ...
    "Every source validation requires finite nonnegative thrust.");

state = package.x_eq(:);
pc = package.parConst;
pc.u_ctrl = package.u_eq(:);

if isfield(pc,"gust")
    pc.gust = zeros(size(pc.gust));
end
if isfield(pc,"N_Thrust")
    pc.N_Thrust = zeros(size(pc.N_Thrust));
end

nonlinear = AeroFlex.sim.nonlinear_terms(state,pc,package.idx);
rawResidual = package.L*state+nonlinear;
propagatedResidual = rawResidual;
propagatedResidual(package.idx.q1) = ...
    package.beam.Pz*rawResidual(package.idx.q1);

[stateJacobian,~] = AeroFlex.sim.nonlinearJacobian( ...
    state,package.idx,pc);
linearOperator = package.Ldyn+stateJacobian;
linearOperator(package.idx.q1,:) = package.Ldyn(package.idx.q1,:)+ ...
    package.beam.Pz*stateJacobian(package.idx.q1,:);

poles = eig(full(linearOperator));
finite = all(isfinite([state;rawResidual;propagatedResidual;poles]));
trimConverged = isfield(trim,"converged") && logical(trim.converged);

row = emptyRow();
row.sourceId = string(sourceRecord.sourceId);
row.coordinates = double(sourceRecord.coordinates(:).');
row.packagePath = relativePath;
row.packageSha256 = actualHash;
row.trimConverged = trimConverged;
row.thrustNewtons = double(trim.thrust);
row.rawResidualNorm = norm(rawResidual);
row.propagatedResidualNorm = norm(propagatedResidual);
row.maximumPoleRealPerSecond = max(real(poles));
row.maximumPoleFrequencyHz = max(abs(imag(poles)))/(2*pi);
row.unstablePoleCount = nnz(real(poles) > 0);
row.poles = poles;
row.status = passFail(trimConverged && finite);

if row.status ~= "PASS"
    row.message = "Source trim or linearization did not satisfy finite/converged checks.";
end
end


function row = emptyRow()
%EMPTYROW Keep the checkpoint schema stable across source failures.

row = struct( ...
    "sourceId","", ...
    "coordinates",[NaN NaN], ...
    "packagePath","", ...
    "packageSha256","", ...
    "trimConverged",false, ...
    "thrustNewtons",NaN, ...
    "rawResidualNorm",NaN, ...
    "propagatedResidualNorm",NaN, ...
    "maximumPoleRealPerSecond",NaN, ...
    "maximumPoleFrequencyHz",NaN, ...
    "unstablePoleCount",NaN, ...
    "poles",complex(zeros(0,1)), ...
    "elapsedSeconds",NaN, ...
    "status","NOT_RUN", ...
    "message","");
end


function figureHandle = plotSummary(rows)
%PLOTSUMMARY Plot physical-node trends without opening figures in the loop.

coordinates = vertcat(rows.coordinates);
residual = [rows.propagatedResidualNorm];
growth = [rows.maximumPoleRealPerSecond];
frequency = [rows.maximumPoleFrequencyHz];
elapsed = [rows.elapsedSeconds];

figureHandle = figure("Visible","off","Color","w", ...
    "Name","Pazy V17 source validation");
layout = tiledlayout(figureHandle,2,2,"TileSpacing","compact", ...
    "Padding","compact");
title(layout,"V17 exact-source validation trends");

ax = nexttile(layout,1);
scatter(ax,coordinates(:,1),coordinates(:,2),55,log10(max(residual,eps)), ...
    "filled");
grid(ax,"on"); colorbar(ax);
xlabel(ax,"Airspeed (m/s)"); ylabel(ax,"Angle of attack (deg)");
title(ax,"log_{10} propagated trim residual");

ax = nexttile(layout,2);
scatter(ax,coordinates(:,1),growth,45,coordinates(:,2),"filled");
grid(ax,"on"); colorbar(ax);
xlabel(ax,"Airspeed (m/s)"); ylabel(ax,"Maximum Re(\lambda) (s^{-1})");
title(ax,"Flexible/aerodynamic spectral abscissa");

ax = nexttile(layout,3);
scatter(ax,coordinates(:,1),frequency,45,coordinates(:,2),"filled");
grid(ax,"on"); colorbar(ax);
xlabel(ax,"Airspeed (m/s)"); ylabel(ax,"Maximum modal frequency (Hz)");
title(ax,"Resolved frequency extent");

ax = nexttile(layout,4);
bar(ax,elapsed);
grid(ax,"on"); xlabel(ax,"Source index"); ylabel(ax,"Elapsed time (s)");
title(ax,"Per-source validation time");
end


function writeSummary(path,report)
%WRITESUMMARY Store scalar rows in JSON; complete poles remain in MAT.

summary = report;
summary.rows = rmfield(summary.rows,"poles");

fileId = fopen(path,"w");
assert(fileId >= 0,"PazyValidation:JsonOpenFailed", ...
    "Unable to open JSON output: %s",path);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,"%s\n",jsonencode(summary,PrettyPrint=true));
clear cleanup
end


function value = fileSha256(path)
%FILESHA256 Compute a lowercase SHA-256 digest from binary file bytes.

fileId = fopen(path,"r");
assert(fileId >= 0,"PazyValidation:HashOpenFailed", ...
    "Unable to open file for hashing: %s",path);
cleanup = onCleanup(@() fclose(fileId));
bytes = fread(fileId,Inf,"*uint8");
clear cleanup

digest = java.security.MessageDigest.getInstance("SHA-256");
digest.update(bytes);
value = lower(join(string(dec2hex(typecast(digest.digest(),"uint8"),2)),""));
end


function status = passFail(passed)
%PASSFAIL Return stable saved/console status text.

status = "FAIL";
if passed
    status = "PASS";
end
end
