function [result,plan] = runPazyModelWorkflow(options)
%RUNPAZYMODELWORKFLOW Prepare or execute the general Pazy model workflow.
%   This entry point retains the established sim_init/sim_run path for
%   wing-only and coupled-full development cases. It is separate from the
%   formally qualified A/B benchmark facade.
%
%   Use Execute=false to inspect the setup request without reading SHARPy
%   products, constructing the model, or creating output files.

arguments
    options.BodyCase (1,1) string ...
        {mustBeMember(options.BodyCase,["wingOnly","coupledFull"])} = ...
        "wingOnly"
    options.SimulationMode (1,1) string ...
        {mustBeMember(options.SimulationMode,["openloop","nmhe_nmpc"])} = ...
        "openloop"
    options.GustEnabled (1,1) logical = true
    options.Execute (1,1) logical = false
    options.PrepareSetup (1,1) logical = true
    options.SetupDirectory (1,1) string = ""
    options.SharpyRoot (1,1) string = ""
    options.CaseName (1,1) string = "pazy_krylov_ROM"
    options.RunId (1,1) string = "auto"
    options.Debug (1,1) logical = false
    options.FiguresVisible (1,1) logical = false
end

repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
sharpyRoot = options.SharpyRoot;
if strlength(sharpyRoot)==0
    sharpyRoot = fullfile(repositoryRoot,"TestBenchPazy");
end
runId = options.RunId;
if runId=="auto"
    runId = lower(options.BodyCase)+"_"+options.SimulationMode+"_" + ...
        string(datetime("now","TimeZone","UTC", ...
        "Format","yyyyMMdd'T'HHmmss'Z'"));
end

plan = struct( ...
    "schemaVersion","pazy-general-model-workflow-v1", ...
    "bodyCase",options.BodyCase, ...
    "simulationMode",options.SimulationMode, ...
    "gustEnabled",options.GustEnabled, ...
    "caseName",options.CaseName, ...
    "sharpyRoot",sharpyRoot, ...
    "prepareSetup",options.PrepareSetup, ...
    "setupDirectory",options.SetupDirectory, ...
    "runId",runId, ...
    "rateProjectionEnabled",options.BodyCase=="coupledFull", ...
    "figuresVisible",options.FiguresVisible, ...
    "formalBenchmarkCase",false, ...
    "executed",false);
result = struct("status","PLAN_ONLY","success",false, ...
    "setup",struct(),"history",struct(),"plan",plan);
if ~options.Execute
    return
end

% Setup generation and simulation are deliberately separate. A prepared
% setup can be reused without charging model construction to online timing.
workflowTimer = tic;
setupProject(ValidateEntryPoints=true);
setupDirectory = options.SetupDirectory;
setup = struct();
setupTimer = tic;
if options.PrepareSetup
    [setup,setupPassed] = sim_init(sharpyRoot, ...
        'case_name',options.CaseName, ...
        'body_case',options.BodyCase, ...
        'sim_case',options.SimulationMode, ...
        'runner','PlantROM', ...
        'gustOn',options.GustEnabled, ...
        'debug',options.Debug, ...
        'overwrite',false, ...
        'date_only_runs',false, ...
        'run_id',runId);
    assert(setupPassed,"AeroFlex:ModelWorkflowSetup", ...
        "The requested model setup did not complete successfully.");
    setupDirectory = string(setup.paths.run_dir);
else
    assert(strlength(setupDirectory)>0 && isfolder(setupDirectory), ...
        "AeroFlex:ModelWorkflowSetupDirectory", ...
        "Provide an existing SetupDirectory when PrepareSetup=false.");
end
preparationSeconds = toc(setupTimer);

oldVisibility = get(groot,"DefaultFigureVisible");
initialFigures = findall(groot,"Type","figure");
set(groot,"DefaultFigureVisible",localVisibility(options.FiguresVisible));
figureCleanup = onCleanup(@()localRestoreFigures( ...
    oldVisibility,initialFigures));
simulationTimer = tic;
[simulationPassed,history] = sim_run(options.CaseName,options.BodyCase, ...
    'setup_dir',setupDirectory, ...
    'date_only_runs',false, ...
    'overwrite',false);
simulationAndPostSeconds = toc(simulationTimer);
assert(simulationPassed,"AeroFlex:ModelWorkflowSimulation", ...
    "The requested model workflow did not complete successfully.");

plan.setupDirectory = setupDirectory;
plan.executed = true;
outputs = localWorkflowOutputs(history);
timing = localWorkflowTiming(history,preparationSeconds, ...
    simulationAndPostSeconds,toc(workflowTimer));
metrics = localWorkflowMetrics(outputs);
result = struct("status","PASS_COMPLETED_MODEL_WORKFLOW", ...
    "success",true,"setup",setup,"history",history,"plan",plan, ...
    "timing",timing,"metrics",metrics,"outputs",outputs);

% Keep the compact run summary beside the established sim_run products.
% The full state history is already owned by the serialized run bundle.
workflowSummary = rmfield(result,["setup","history"]);
summaryPath = fullfile(outputs.runRoot,"model_workflow_summary.mat");
save(summaryPath,"workflowSummary","-v7");
result.outputs.summaryMat = string(summaryPath);
clear figureCleanup
end

function value = localVisibility(isVisible)
if isVisible
    value = "on";
else
    value = "off";
end
end

function localRestoreFigures(oldVisibility,initialFigures)
currentFigures = findall(groot,"Type","figure");
newFigures = setdiff(currentFigures,initialFigures);
if ~isempty(newFigures)
    delete(newFigures);
end
set(groot,"DefaultFigureVisible",oldVisibility);
end

function outputs = localWorkflowOutputs(history)
assert(isfield(history,"cfg") && isfield(history.cfg,"paths") && ...
    isfield(history.cfg.paths,"run_dir"), ...
    "AeroFlex:ModelWorkflowOutputRoot", ...
    "The completed run did not report its output directory.");
paths = history.cfg.paths;
runRoot = string(paths.run_dir);
plotsRoot = localPathField(paths,"plots",fullfile(runRoot,"plots"));
postMat = localPathField(paths,"for_matlab", ...
    fullfile(runRoot,"for_matlab"));
postMat = fullfile(postMat,"post_out.mat");
plotFiles = strings(0,1);
if isfolder(plotsRoot)
    listing = [dir(fullfile(plotsRoot,"*.png")); ...
        dir(fullfile(plotsRoot,"*.pdf"))];
    plotFiles = string(fullfile({listing.folder},{listing.name})).';
end
outputs = struct("runRoot",runRoot,"plotsRoot",string(plotsRoot), ...
    "postProcessMat",string(postMat),"plotFiles",plotFiles, ...
    "summaryMat","");
end

function timing = localWorkflowTiming(history,preparationSeconds, ...
        simulationAndPostSeconds,totalEntrySeconds)
simulatedSeconds = nan;
if isfield(history,"t") && ~isempty(history.t)
    simulatedSeconds = max(history.t)-min(history.t);
end
timing = struct( ...
    "preparationSeconds",preparationSeconds, ...
    "simulationAndPostProcessingSeconds",simulationAndPostSeconds, ...
    "totalEntrySeconds",totalEntrySeconds, ...
    "simulatedSeconds",simulatedSeconds, ...
    "conservativeRealTimeFactor", ...
        simulatedSeconds/max(simulationAndPostSeconds,eps), ...
    "onlineComponentTiming",localComponentTiming(history));
end

function distributions = localComponentTiming(history)
distributions = struct();
if ~isfield(history,"log") || ~isfield(history.log,"diag")
    return
end
diagnostics = history.log.diag;
fields = ["scheduleTime","senseTime","mheTime","mpcTime", ...
    "allocationTime","fusionTime","actuatorTime","plantTime"];
labels = ["scheduling","sensing","nMhe","nMpc", ...
    "allocation","fusion","actuator","plant"];
for index = 1:numel(fields)
    if ~isfield(diagnostics,fields(index))
        continue
    end
    values = diagnostics.(fields(index));
    values = values(isfinite(values) & values>=0);
    if isempty(values)
        continue
    end
    distributions.(labels(index)) = struct( ...
        "meanSeconds",mean(values), ...
        "p95Seconds",prctile(values,95), ...
        "maximumSeconds",max(values), ...
        "sampleCount",numel(values));
end
end

function metrics = localWorkflowMetrics(outputs)
metrics = struct();
if ~isfile(outputs.postProcessMat)
    return
end
saved = load(outputs.postProcessMat,"out");
if isfield(saved,"out") && isfield(saved.out,"metrics")
    metrics = saved.out.metrics;
end
end

function value = localPathField(paths,name,fallback)
if isfield(paths,name) && ~isempty(paths.(name))
    value = string(paths.(name));
else
    value = string(fallback);
end
end
