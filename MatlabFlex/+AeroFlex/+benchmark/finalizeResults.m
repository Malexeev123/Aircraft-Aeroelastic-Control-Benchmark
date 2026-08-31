function products = finalizeResults(summary,plan,layout,wallSeconds)
%FINALIZERESULTS Write standardized metrics and benchmark trend plots.
%   This function consumes saved runner output after integration has ended.
%   It does not alter the plant, estimator, controller, or accepted commands.

arguments
    summary (1,1) struct
    plan (1,1) struct
    layout (1,1) struct
    wallSeconds (1,1) double {mustBeNonnegative,mustBeFinite}
end

% A runner writes the authoritative history before its final qualification
% assertion. Recovering that file lets failed experimental runs remain useful
% without changing their status or rerunning the plant.
artifact = AeroFlex.benchmark.internal.recoverExecutionArtifact(layout.runRoot);
if artifact.available
    runnerSummary = artifact.summary;
    history = artifact.history;
else
    runnerSummary = summary;
    history = struct();
end

% Post-processing is intentionally outside the online simulation timer.
metricTimer = tic;
metrics = localMetrics(runnerSummary,history,plan,wallSeconds);
metrics.offlineTiming.metricAssemblySeconds = toc(metricTimer);
jsonPath = fullfile(layout.metrics,"metrics.json");
csvPath = fullfile(layout.metrics,"metrics.csv");
matPath = fullfile(layout.data,"standardized_results.mat");
plotPaths = strings(0,1);
plotTimer = tic;
if plan.savePlots && artifact.available && isfield(history,"log")
    plotPaths = localPlotHistory(history,plan,layout);
end
metrics.offlineTiming.plottingSeconds = toc(plotTimer);

% The first-pass duration is recorded before the self-describing metrics file
% is rewritten with that duration. MAT output remains v7 on the WSL path.
serializationTimer = tic;
[names,values,units] = localMetricRows(metrics);
writetable(table(names,values,units),csvPath);
save(matPath,"metrics","runnerSummary","-v7");
metrics.offlineTiming.serializationFirstPassSeconds = ...
    toc(serializationTimer);
localWriteJson(jsonPath,metrics);
save(matPath,"metrics","runnerSummary","-v7");
products = struct("metrics",metrics,"metricsJson",jsonPath, ...
    "metricsCsv",csvPath,"dataMat",matPath, ...
    "runnerArtifact",artifact.path,"plots",plotPaths);
end

function metrics = localMetrics(summary,history,plan,wallSeconds)
metrics = struct( ...
    "schemaVersion","pazy-benchmark-metrics-v1", ...
    "caseId",plan.caseId, ...
    "qualificationStatus",plan.qualificationStatus, ...
    "simulatedSeconds",plan.durationSeconds, ...
    "wallSeconds",wallSeconds, ...
    "onlineWallSeconds",wallSeconds, ...
    "runnerExecutionWallSeconds",wallSeconds, ...
    "preparationSeconds",localField(plan,"preparationSeconds"), ...
    "wallRealTimeFactor",plan.durationSeconds/max(wallSeconds,eps), ...
    "onlineTiming",struct(),"offlineTiming",struct(),"solver",struct(), ...
    "tracking",struct(),"estimation",struct(), ...
    "actuator",struct(),"loads",struct(),"constraints",struct());
metrics.solver.nMheSuccessRate = localNested(summary,["mhe","successRate"]);
metrics.solver.nMpcSuccessRate = localNested(summary,["mpc","successRate"]);
metrics.solver.nMpcFallbackCount = localField(summary,"mpcFallbackCount");
metrics.constraints.minimumThrustNewton = ...
    localField(summary,"minimumThrustNewton");
metrics.constraints.maximumSourceDomainRatio = ...
    localNested(summary,["realizedSourceDomain","maximumRatio"]);
metrics.constraints.sourceDomainRejectedIntervals = ...
    localNested(summary,["realizedSourceDomain","rejectedIntervalCount"]);

if ~isfield(history,"log") || ~isstruct(history.log)
    return
end
log = history.log;
if isfield(log,"controlDiagnostics")
    diagnostics = log.controlDiagnostics;
    count = localValidCount(diagnostics);
    timingFields = ["scheduleTime","senseTime","mheTime","mpcTime", ...
        "allocationTime","fusionTime","actuatorTime","plantTime"];
    labels = ["scheduling","sensing","nMhe","nMpc", ...
        "allocation","fusion","actuator","plant"];
    for index = 1:numel(timingFields)
        if isfield(diagnostics,timingFields(index))
            values = diagnostics.(timingFields(index));
            values = values(1:min(count,numel(values)));
            metrics.onlineTiming.(labels(index)) = localDistribution(values);
        end
    end
    onlineWallSeconds = localOnlineWallSeconds(diagnostics,count);
    if isfinite(onlineWallSeconds)
        metrics.onlineWallSeconds = onlineWallSeconds;
        metrics.wallSeconds = onlineWallSeconds;
        metrics.wallRealTimeFactor = ...
            plan.durationSeconds/max(onlineWallSeconds,eps);
    elseif isfield(summary,"wallRealTimeFactor")
        metrics.wallRealTimeFactor = summary.wallRealTimeFactor;
    end
    if all(isfield(diagnostics,["mheAccepted","mheAttempted"]))
        attempted = logical(diagnostics.mheAttempted(1:count));
        accepted = logical(diagnostics.mheAccepted(1:count));
        if any(attempted)
            metrics.solver.nMheSuccessRate = mean(accepted(attempted));
        else
            metrics.solver.nMheSuccessRate = NaN;
        end
    end
    if all(isfield(diagnostics,["mpcAccepted","mpcAttempted"]))
        attempted = logical(diagnostics.mpcAttempted(1:count));
        accepted = logical(diagnostics.mpcAccepted(1:count));
        if any(attempted)
            metrics.solver.nMpcSuccessRate = mean(accepted(attempted));
            metrics.solver.nMpcFallbackCount = ...
                nnz(~accepted(attempted));
        else
            metrics.solver.nMpcSuccessRate = NaN;
            metrics.solver.nMpcFallbackCount = 0;
        end
    end
end

if all(isfield(log,["uOuterRigidState","uOuterRigidReference"]))
    count = min(size(log.uOuterRigidState,2), ...
        size(log.uOuterRigidReference,2));
    state = log.uOuterRigidState(:,1:count);
    reference = log.uOuterRigidReference(:,1:count);
    speedError = vecnorm(state(1:3,:),2,1)- ...
        vecnorm(reference(1:3,:),2,1);
    metrics.tracking.speedErrorMetersPerSecond = ...
        localErrorMetrics(speedError);
    metrics.tracking.pitchErrorDegrees = localErrorMetrics( ...
        rad2deg(state(5,:)-reference(5,:)));
end
if isfield(log,"uOuterRigidPosition")
    altitude = -log.uOuterRigidPosition(3,:);
    metrics.tracking.altitudeDeviationMeters = localErrorMetrics(altitude);
end
if all(isfield(log,["estimatorTruthGust","wHorizon"]))
    count = min(numel(log.estimatorTruthGust),size(log.wHorizon,2));
    gustError = log.wHorizon(end,1:count)- ...
        log.estimatorTruthGust(1:count);
    metrics.estimation.gustErrorMetersPerSecond = ...
        localErrorMetrics(gustError);
end
if isfield(log,"uActuatorEndpoint")
    endpoint = log.uActuatorEndpoint;
    half = floor(size(endpoint,1)/2);
    metrics.actuator.positionPeakRadians = ...
        max(abs(endpoint(1:half,:)),[],"all","omitnan");
    metrics.actuator.ratePeakRadiansPerSecond = ...
        max(abs(endpoint(half+1:2*half,:)),[],"all","omitnan");
end
if isfield(log,"loads") && isstruct(log.loads)
    if isfield(log.loads,"Clamp6")
        clamp = log.loads.Clamp6;
        metrics.loads.rootFzNewton = localSignalMetrics(clamp(3,:));
        metrics.loads.rootMyNewtonMeter = localSignalMetrics(clamp(5,:));
    end
    if isfield(log.loads,"Fthrust_B")
        metrics.constraints.minimumThrustNewton = ...
            min(log.loads.Fthrust_B(1,:),[],"all","omitnan");
    end
end
if isfield(log,"physicalWingtip") && ...
        isfield(log.physicalWingtip,"symmetricMeanAbsoluteMeters")
    tip = log.physicalWingtip.symmetricMeanAbsoluteMeters;
    metrics.loads.wingtipMeters = localSignalMetrics(tip-tip(1));
end
end

function paths = localPlotHistory(history,plan,layout)
log = history.log;
assert(isfield(log,"controlDiagnostics"), ...
    "AeroFlex:BenchmarkPlotDiagnostics", ...
    "The runner artifact has no control diagnostics to plot.");
diagnostics = log.controlDiagnostics;
count = localValidCount(diagnostics);
time = diagnostics.tCtrl(1:count);
handle = figure("Visible",localVisibility(plan.figuresVisible), ...
    "Color","w","Position",[20 20 1600 1050], ...
    "Name","Pazy benchmark trends");
cleanup = onCleanup(@()localCloseFigure(handle));
tiles = tiledlayout(5,3,"TileSpacing","compact","Padding","compact");

localPlotRigid(time,log,count);
localPlotEstimation(time,history,count);
localPlotCommands(time,log,count);
localPlotLoads(time,log,count);
localPlotSolvers(time,diagnostics,count);
localPlotTiming(time,diagnostics,count);
localFormatSingletonHistory(handle,time);
xlabel(tiles,"Time (s)");
if ~plan.publicationMode
    title(tiles,sprintf("Pazy %s benchmark trends (%s)", ...
        plan.caseId,plan.qualificationStatus),"Interpreter","none");
end
path = fullfile(layout.publicationPlots, ...
    lower(plan.caseId)+"_benchmark_trends.png");
exportgraphics(handle,path,"Resolution",220);
paths = [string(path);localExportIndividualPlots( ...
    handle,plan,layout)];
clear cleanup
end

function paths = localExportIndividualPlots(handle,plan,layout)
names = ["airspeed","pitch","altitude","gust_estimate", ...
    "state_estimation_error","wingtip","actuator_positions", ...
    "actuator_rates","thrust","root_fz","root_my", ...
    "source_domain_gate","solver_status","primary_timing", ...
    "handoff_timing"];
axesHandles = findall(handle,"Type","axes");
tileNumbers = nan(size(axesHandles));
for index = 1:numel(axesHandles)
    tileNumbers(index) = double(axesHandles(index).Layout.Tile);
end
[tileNumbers,order] = sort(tileNumbers);
axesHandles = axesHandles(order);
valid = isfinite(tileNumbers) & tileNumbers>=1 & ...
    tileNumbers<=numel(names);
axesHandles = axesHandles(valid);
tileNumbers = tileNumbers(valid);

paths = strings(0,1);
for index = 1:numel(axesHandles)
    sourceAxis = axesHandles(index);
    figureHandle = figure("Visible",localVisibility(plan.figuresVisible), ...
        "Color","w","Position",[40 40 760 520], ...
        "Name","Pazy benchmark result");
    cleanup = onCleanup(@()localCloseFigure(figureHandle));
    axisHandle = copyobj(sourceAxis,figureHandle);
    axisHandle.Units = "normalized";
    axisHandle.Position = [0.14 0.15 0.82 0.80];
    xlabel(axisHandle,"Time (s)");
    localCopyLegend(handle,sourceAxis,axisHandle);
    if plan.publicationMode
        title(axisHandle,"");
        axisHandle.FontName = "Arial";
        axisHandle.FontSize = 9;
    end
    stem = lower(plan.caseId)+"_"+names(tileNumbers(index));
    pngPath = fullfile(layout.publicationPlots,stem+".png");
    exportgraphics(figureHandle,pngPath,"Resolution",300);
    paths(end+1,1) = string(pngPath); %#ok<AGROW>
    if plan.publicationMode
        pdfPath = fullfile(layout.publicationPlots,stem+".pdf");
        exportgraphics(figureHandle,pdfPath,"ContentType","vector");
        paths(end+1,1) = string(pdfPath); %#ok<AGROW>
    end
    clear cleanup
end
end

function localCopyLegend(sourceFigure,sourceAxis,targetAxis)
legends = findall(sourceFigure,"Type","legend");
for index = 1:numel(legends)
    candidate = legends(index);
    if isprop(candidate,"Axes") && candidate.Axes==sourceAxis
        legend(targetAxis,string(candidate.String), ...
            "Location",candidate.Location);
        return
    end
end
end

function localPlotRigid(time,log,count)
nexttile
if all(isfield(log,["uOuterRigidState","uOuterRigidReference"]))
    state = log.uOuterRigidState(:,1:count);
    reference = log.uOuterRigidReference(:,1:count);
    plot(time,vecnorm(state(1:3,:),2,1),"k","LineWidth",1.2); hold on
    plot(time,vecnorm(reference(1:3,:),2,1),"--","LineWidth",1.2);
    legend(["Achieved","Reference"],"Location","best");
end
ylabel("Speed (m/s)"); title("Airspeed"); grid on
nexttile
if all(isfield(log,["uOuterRigidState","uOuterRigidReference"]))
    plot(time,rad2deg(log.uOuterRigidState(5,1:count)),"LineWidth",1.2); hold on
    plot(time,rad2deg(log.uOuterRigidReference(5,1:count)),"--", ...
        "LineWidth",1.2);
    legend(["Achieved","Reference"],"Location","best");
end
ylabel("Pitch (deg)"); title("Attitude"); grid on
nexttile
if isfield(log,"uOuterRigidPosition")
    plot(time,-log.uOuterRigidPosition(3,1:count),"LineWidth",1.2);
end
ylabel("Altitude deviation (m)"); title("Altitude"); grid on
end

function localPlotEstimation(time,history,count)
log = history.log;
nexttile
if all(isfield(log,["estimatorTruthGust","wHorizon"]))
    stairs(time,log.estimatorTruthGust(1:count),"k--","LineWidth",1.1); hold on
    stairs(time,log.wHorizon(end,1:count),"LineWidth",1.1);
    legend(["Truth","Estimate"],"Location","best");
end
ylabel("Gust (m/s)"); title("Gust estimation"); grid on
nexttile
if isfield(log,"xhat") && isfield(history,"x")
    stateCount = min([count,size(log.xhat,2),size(history.x,2)]);
    stateRows = min(size(log.xhat,1),size(history.x,1));
    errorInf = max(abs(log.xhat(1:stateRows,1:stateCount)- ...
        history.x(1:stateRows,1:stateCount)),[],1);
    semilogy(time(1:stateCount),max(errorInf,eps),"LineWidth",1.1);
end
ylabel("Error infinity norm"); title("State estimation"); grid on
nexttile
if isfield(log,"physicalWingtip") && ...
        isfield(log.physicalWingtip,"symmetricMeanAbsoluteMeters")
    tip = log.physicalWingtip.symmetricMeanAbsoluteMeters;
    tTip = linspace(time(1),time(end),numel(tip));
    plot(tTip,tip-tip(1),"LineWidth",1.1);
end
ylabel("Wingtip excursion (m)"); title("Flexible response"); grid on
end

function localPlotCommands(time,log,count)
nexttile
if isfield(log,"uActuatorEndpoint")
    endpoint = log.uActuatorEndpoint(:,1:count);
    half = floor(size(endpoint,1)/2);
    localStairsSeries(time,endpoint(1:half,:));
end
ylabel("Position (rad)"); title("Actuator positions"); grid on
nexttile
if isfield(log,"uActuatorEndpoint")
    endpoint = log.uActuatorEndpoint(:,1:count);
    half = floor(size(endpoint,1)/2);
    localStairsSeries(time,endpoint(half+1:2*half,:));
end
ylabel("Rate (rad/s)"); title("Actuator rates"); grid on
nexttile
if isfield(log,"loads") && isfield(log.loads,"Fthrust_B")
    stairs(time,log.loads.Fthrust_B(1,1:count),"LineWidth",1.1); hold on
    yline(0,"k--","T = 0");
end
ylabel("Thrust (N)"); title("Thrust"); grid on
end

function localPlotLoads(time,log,count)
nexttile
if isfield(log,"loads") && isfield(log.loads,"Clamp6")
    plot(time,log.loads.Clamp6(3,1:count),"LineWidth",1.1);
end
ylabel("F_z (N)"); title("Wing-root force"); grid on
nexttile
if isfield(log,"loads") && isfield(log.loads,"Clamp6")
    plot(time,log.loads.Clamp6(5,1:count),"LineWidth",1.1);
end
ylabel("M_y (N m)"); title("Wing-root bending"); grid on
nexttile
if isfield(log,"reciprocalSourceDomainRatio")
    ratio = log.reciprocalSourceDomainRatio(:,1:count);
    semilogy(time,max(ratio,[],1,"omitnan"),"LineWidth",1.1); hold on
    yline(1,"k--","Gate");
end
ylabel("Ratio (-)"); title("Source-domain gate"); grid on
end

function localPlotSolvers(time,diagnostics,count)
nexttile
if all(isfield(diagnostics,["mheAccepted","mpcAccepted"]))
    stairs(time,double(diagnostics.mheAccepted(1:count)),"LineWidth",1); hold on
    stairs(time,double(diagnostics.mpcAccepted(1:count)),"--","LineWidth",1);
    ylim([-0.05,1.05]); legend(["nMHE","nMPC"],"Location","best");
end
ylabel("Accepted (-)"); title("Solver status"); grid on
end

function localPlotTiming(time,diagnostics,count)
nexttile
fields = ["plantTime","mheTime","mpcTime"];
labels = ["Plant","nMHE","nMPC"];
plotted = strings(0,1);
for index = 1:numel(fields)
    if isfield(diagnostics,fields(index))
        values = diagnostics.(fields(index));
        plotCount = min([count,numel(time),numel(values)]);
        stairs(time(1:plotCount),values(1:plotCount),"LineWidth",1); hold on
        plotted(end+1,1) = labels(index); %#ok<AGROW>
    end
end
ylabel("Online time (s)"); title("Primary online costs");
if ~isempty(plotted), legend(plotted,"Location","best"); end
grid on
nexttile
fields = ["scheduleTime","fusionTime","actuatorTime"];
labels = ["Schedule","Fusion","Actuator"];
plotted = strings(0,1);
for index = 1:numel(fields)
    if isfield(diagnostics,fields(index))
        values = diagnostics.(fields(index));
        plotCount = min([count,numel(time),numel(values)]);
        stairs(time(1:plotCount),values(1:plotCount),"LineWidth",1); hold on
        plotted(end+1,1) = labels(index); %#ok<AGROW>
    end
end
ylabel("Online time (s)"); title("Handoff costs");
if ~isempty(plotted), legend(plotted,"Location","best"); end
grid on
end

function count = localValidCount(diagnostics)
if isfield(diagnostics,"iCtrl")
    count = double(diagnostics.iCtrl);
else
    count = numel(diagnostics.tCtrl);
end
count = min(count,numel(diagnostics.tCtrl));
end

function result = localDistribution(values)
values = values(isfinite(values));
if isempty(values)
    result = struct("meanSeconds",nan,"p95Seconds",nan,"maxSeconds",nan);
else
    result = struct("meanSeconds",mean(values), ...
        "p95Seconds",prctile(values,95),"maxSeconds",max(values));
end
end

function seconds = localOnlineWallSeconds(diagnostics,count)
seconds = nan;
if isfield(diagnostics,"totalCtrlTime")
    control = diagnostics.totalCtrlTime(1:min(count, ...
        numel(diagnostics.totalCtrlTime)));
    if isfield(diagnostics,"plantTime")
        plant = diagnostics.plantTime(1:min(count, ...
            numel(diagnostics.plantTime)));
        sampleCount = min(numel(control),numel(plant));
        values = control(1:sampleCount)+plant(1:sampleCount);
    else
        values = control;
    end
else
    fields = ["scheduleTime","senseTime","mheTime","mpcTime", ...
        "allocationTime","fusionTime","actuatorTime","plantTime"];
    values = zeros(1,count);
    available = false;
    for field = fields
        if isfield(diagnostics,field)
            component = diagnostics.(field);
            component = component(1:min(count,numel(component)));
            values(1:numel(component)) = ...
                values(1:numel(component))+component;
            available = true;
        end
    end
    if ~available
        return
    end
end
values = values(isfinite(values));
if ~isempty(values)
    seconds = sum(values);
end
end

function result = localErrorMetrics(values)
values = values(isfinite(values));
if isempty(values)
    result = struct("rms",nan,"peakAbsolute",nan,"terminal",nan);
else
    result = struct("rms",rms(values), ...
        "peakAbsolute",max(abs(values)),"terminal",values(end));
end
end

function result = localSignalMetrics(values)
values = values(isfinite(values));
if isempty(values)
    result = struct("rms",nan,"peakAbsolute",nan);
else
    result = struct("rms",rms(values),"peakAbsolute",max(abs(values)));
end
end

function value = localNested(source,path)
value = nan;
for name = path
    if ~isstruct(source) || ~isfield(source,name), return, end
    source = source.(name);
end
if isnumeric(source) || islogical(source), value = double(source); end
end

function value = localField(source,name)
value = nan;
if isfield(source,name) && (isnumeric(source.(name)) || islogical(source.(name)))
    value = double(source.(name));
end
end

function [names,values,units] = localMetricRows(metrics)
rows = {
    "online_wall_seconds",metrics.onlineWallSeconds,"s";
    "runner_execution_wall_seconds",metrics.runnerExecutionWallSeconds,"s";
    "wall_real_time_factor",metrics.wallRealTimeFactor,"-";
    "nmhe_success_rate",metrics.solver.nMheSuccessRate,"fraction";
    "nmpc_success_rate",metrics.solver.nMpcSuccessRate,"fraction";
    "nmpc_fallback_count",metrics.solver.nMpcFallbackCount,"count";
    "minimum_thrust",metrics.constraints.minimumThrustNewton,"N";
    "maximum_source_domain_ratio", ...
        metrics.constraints.maximumSourceDomainRatio,"-"};
names = string(rows(:,1));
values = cell2mat(rows(:,2));
units = string(rows(:,3));
end

function localWriteJson(path,value)
file = fopen(path,"w");
assert(file>=0,"AeroFlex:BenchmarkMetricsWrite","Cannot write %s.",path);
cleanup = onCleanup(@()fclose(file));
fprintf(file,"%s\n",jsonencode(value,"PrettyPrint",true));
clear cleanup
end

function value = localVisibility(visible)
if visible, value = "on"; else, value = "off"; end
end

function localCloseFigure(handle)
if isgraphics(handle), close(handle); end
end

function localFormatSingletonHistory(handle,time)
if ~isscalar(time)
    return
end

% A production-entry smoke contains only one held control interval. Mark its
% samples explicitly; otherwise line objects are valid but visually empty.
graphics = findall(handle);
for index = 1:numel(graphics)
    item = graphics(index);
    if isprop(item,"XData") && isscalar(item.XData) && ...
            isprop(item,"Marker")
        item.Marker = "o";
        item.MarkerSize = 4;
    end
end
span = 0.01;
axesHandles = findall(handle,"Type","axes");
for index = 1:numel(axesHandles)
    xlim(axesHandles(index),[max(0,time(1)-span),time(1)+span]);
end
end

function localStairsSeries(time,signals)
% Plot each held signal explicitly so a one-sample, multi-channel history is
% not interpreted as a single vector with incompatible X/Y lengths.
for row = 1:size(signals,1)
    stairs(time,signals(row,:),"LineWidth",1); hold on
end
end
