function [summary,plan] = runBenchmarkCase(caseId,options)
%RUNBENCHMARKCASE Plan or execute a Pazy benchmark scenario.
%   [SUMMARY,PLAN] = RUNBENCHMARKCASE(CASEID) executes a qualified case.
%   Use Execute=false to inspect the complete resolved plan without changing
%   files, creating figures, or starting MATLAB simulation work.
%
%   Formal Case A members are qualified. Case B uses the retained production
%   production scheduler and accelerated runtime. Its qualification state is
%   reported independently from execution availability.

arguments
    caseId (1,1) string
    options.Execute (1,1) logical = true
    options.DurationSeconds (1,1) double = nan
    options.ControllerMode (1,1) string = "full"
    options.DiagnosticsLevel (1,1) double = 1
    options.FiguresVisible (1,1) logical = false
    options.SavePlots (1,1) logical = true
    options.PublicationMode (1,1) logical = false
    options.OutputRoot (1,1) string = ""
    options.RunId (1,1) string = "auto"
    options.NativeKernelPolicy (1,1) string = "auto"
    options.AllowUnqualified (1,1) logical = false
    options.InitialFlightCondition (1,1) struct = struct()
    options.TerminalFlightCondition (1,1) struct = struct()
end

% Resolve the immutable case intent before any setup or file creation.
repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
definition = AeroFlex.benchmark.caseDefinition(caseId);
definition.initialFlightCondition = localFormalCondition( ...
    options.InitialFlightCondition,definition.initialFlightCondition, ...
    "initial");
definition.terminalFlightCondition = localFormalCondition( ...
    options.TerminalFlightCondition,definition.terminalFlightCondition, ...
    "terminal");
settings = AeroFlex.benchmark.simulationOptions( ...
    Execute=options.Execute, ...
    DurationSeconds=options.DurationSeconds, ...
    ControllerMode=options.ControllerMode, ...
    DiagnosticsLevel=options.DiagnosticsLevel, ...
    FiguresVisible=options.FiguresVisible, ...
    SavePlots=options.SavePlots, ...
    PublicationMode=options.PublicationMode, ...
    OutputRoot=options.OutputRoot, ...
    RunId=options.RunId, ...
    NativeKernelPolicy=options.NativeKernelPolicy, ...
    AllowUnqualified=options.AllowUnqualified);
plan = AeroFlex.benchmark.resolvePlan(definition,settings,repositoryRoot);
summary = struct("status","PLAN_ONLY","caseId",plan.caseId, ...
    "qualificationStatus",plan.qualificationStatus, ...
    "executed",false,"plan",plan);
if ~settings.execute
    return
end

% Execution gates remain separate from plan inspection. Case C and unsupported
% custom configurations continue to fail closed.
assert(plan.executionAllowed,"AeroFlex:BenchmarkCaseUnavailable", ...
    "%s",plan.unavailableReason);
assert(plan.controllerMode=="full" || ...
    ~ismember(plan.caseId,["B1","B2"]), ...
    "AeroFlex:BenchmarkControllerMode", ...
    "The retained Case-B profile requires ControllerMode='full'.");

entryTimer = tic;
projectInfo = setupProject(ValidateEntryPoints=true);
localCheckNativeKernelPolicy(plan,projectInfo);
layout = AeroFlex.benchmark.outputLayout(plan,Create=true);
plan.outputRoot = layout.runRoot;

oldVisibility = get(groot,"DefaultFigureVisible");
initialFigures = findall(groot,"Type","figure");
if plan.figuresVisible
    set(groot,"DefaultFigureVisible","on");
else
    set(groot,"DefaultFigureVisible","off");
end
figureCleanup = onCleanup(@()localRestoreFigures( ...
    oldVisibility,initialFigures));

% Only the numerical owner is charged to simulation wall time. Environment
% setup, hashing, plotting, and serialization are reported separately.
AeroFlex.benchmark.printBanner(plan,"start");
preparationSeconds = toc(entryTimer);
wallTimer = tic;
executionError = [];
try
    summary = AeroFlex.benchmark.internal.executePlan(plan,layout.runRoot);
catch exception
    executionError = exception;
    artifact = AeroFlex.benchmark.internal.recoverExecutionArtifact( ...
        layout.runRoot);
    if localRecoverableValidationFailure(plan,exception,artifact)
        summary = artifact.summary;
        summary.status = "VALIDATION_PENDING_"+string(summary.status);
        summary.executionError = struct("identifier", ...
            string(exception.identifier),"message",string(exception.message));
    else
        summary = struct("status","ERROR","caseId",plan.caseId, ...
            "executionError",struct("identifier", ...
            string(exception.identifier),"message",string(exception.message)));
    end
end
wallSeconds = toc(wallTimer);

% Preserve useful validation evidence even when the runner deliberately
% fails its final qualification assertion. All other errors are rethrown
% after the failure manifest and any recoverable products are written.
plan.executed = true;
plan.preparationSeconds = preparationSeconds;
summary.entryPlan = plan;
summary.wallSeconds = wallSeconds;
summary.preparationSeconds = preparationSeconds;
summary.outputLayout = layout;
products = AeroFlex.benchmark.finalizeResults( ...
    summary,plan,layout,wallSeconds);
summary.products = products;
summary.onlineWallSeconds = products.metrics.onlineWallSeconds;
summary.onlineRealTimeFactor = products.metrics.wallRealTimeFactor;
summary.totalEntryWallSeconds = toc(entryTimer);
manifest = AeroFlex.benchmark.writeManifest( ...
    plan,layout,summary,projectInfo);
summary.manifest = manifest;
AeroFlex.benchmark.printBanner(plan,"complete",struct( ...
    "status",localSummaryStatus(summary), ...
    "onlineWallSeconds",summary.onlineWallSeconds, ...
    "runnerExecutionWallSeconds",wallSeconds, ...
    "preparationSeconds",preparationSeconds, ...
    "totalEntryWallSeconds",summary.totalEntryWallSeconds, ...
    "runRoot",layout.runRoot));
clear figureCleanup
if ~isempty(executionError)
    artifact = AeroFlex.benchmark.internal.recoverExecutionArtifact( ...
        layout.runRoot);
    if ~localRecoverableValidationFailure(plan,executionError,artifact)
        rethrow(executionError)
    end
end
end

function condition = localFormalCondition(requested,frozen,label)
% Named cases accept explicit conditions only when they restate the manifest.
if isempty(fieldnames(requested))
    condition = frozen;
    return
end
required = ["airspeedMps","angleOfAttackDeg"];
assert(all(isfield(requested,required)) && ...
    all(isfield(frozen,required)), ...
    "AeroFlex:BenchmarkFlightCondition", ...
    "The %s flight condition requires airspeedMps and angleOfAttackDeg.", ...
    label);
condition = struct("airspeedMps",double(requested.airspeedMps), ...
    "angleOfAttackDeg",double(requested.angleOfAttackDeg));
assert(isscalar(condition.airspeedMps) && ...
    isfinite(condition.airspeedMps) && condition.airspeedMps>0 && ...
    isscalar(condition.angleOfAttackDeg) && ...
    isfinite(condition.angleOfAttackDeg), ...
    "AeroFlex:BenchmarkFlightCondition", ...
    "The %s flight condition must contain finite scalar values.",label);
tolerance = 1e-12;
assert(abs(condition.airspeedMps-frozen.airspeedMps)<=tolerance && ...
    abs(condition.angleOfAttackDeg- ...
    frozen.angleOfAttackDeg)<=tolerance, ...
    "AeroFlex:BenchmarkFrozenFlightCondition", ...
    "The named case has a frozen %s condition of U_inf=%.12g m/s and " + ...
    "alpha=%.12g deg. Use a custom case for a different condition.", ...
    label,frozen.airspeedMps,frozen.angleOfAttackDeg);
end

function localCheckNativeKernelPolicy(plan,projectInfo)
if plan.nativeKernelPolicy~="required"
    return
end


if ismember(plan.caseId,["B1","B2"])
    available = projectInfo.scheduledCompiledControl.active;
    owner = "scheduled Case-B";
else
    available = projectInfo.compiledControl.active;
    owner = "fixed-source Case-A";
end
assert(available,"AeroFlex:BenchmarkNativeKernelRequired", ...
    "NativeKernelPolicy='required', but the %s kernels are unavailable.", ...
    owner);
end

function status = localSummaryStatus(summary)
if isfield(summary,"status")
    status = string(summary.status);
elseif isfield(summary,"passed") && summary.passed
    status = "PASS";
else
    status = "COMPLETED";
end
end

function recoverable = localRecoverableValidationFailure( ...
        plan,exception,artifact)
recoverable = ismember(plan.caseId,["B1","B2"]) && artifact.available && ...
    string(exception.identifier)=="Phase18C:CaseBProfile";
end

function localRestoreFigures(oldVisibility,initialFigures)
currentFigures = findall(groot,"Type","figure");
newFigures = setdiff(currentFigures,initialFigures);
if ~isempty(newFigures)
    delete(newFigures);
end
set(groot,"DefaultFigureVisible",oldVisibility);
end
