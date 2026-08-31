%% Pazy benchmark: SHARPy wingtip comparison
% Reproduce the accepted trim-relative wingtip comparison at 40 m/s and
% 1 degree angle of attack. The supplied data preserve the unchanged SHARPy
% reference trace and the matched MATLAB wing-only response.

if ~exist("comparisonSettings","var") || ~isstruct(comparisonSettings)
    comparisonSettings = struct();
end

comparisonSettings = applyDefaults(comparisonSettings);
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
setupProject(ValidateEntryPoints=true,ChangeCurrentFolder=false);

dataPath = fullfile(repositoryRoot,"results","validation", ...
    "sharpy-wingtip-comparison","data","wingtip_comparison.csv");
summaryPath = fullfile(repositoryRoot,"results","validation", ...
    "sharpy-wingtip-comparison","wingtip_comparison_summary.json");

assert(isfile(dataPath),"PazyValidation:MissingWingtipComparisonData", ...
    "The supplied wingtip comparison data are missing: %s",dataPath);
assert(isfile(summaryPath),"PazyValidation:MissingWingtipComparisonSummary", ...
    "The supplied wingtip comparison summary is missing: %s",summaryPath);

data = readtable(dataPath,"VariableNamingRule","preserve");
accepted = jsondecode(fileread(summaryPath));
metrics = calculateMetrics(data);
checks = validateComparison(data,metrics,accepted);

runStamp = string(datetime("now","Format","yyyyMMdd'T'HHmmss"));
outputRoot = fullfile(repositoryRoot,comparisonSettings.outputRoot,runStamp);
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

comparisonReport = struct();
comparisonReport.schemaVersion = 1;
comparisonReport.status = passFail(all(struct2array(checks)));
comparisonReport.flightCondition = accepted.flightCondition;
comparisonReport.gust = accepted.gust;
comparisonReport.comparisonFrame = accepted.comparisonFrame;
comparisonReport.projectionPolicy = accepted.projectionPolicy;
comparisonReport.metrics = metrics;
comparisonReport.checks = checks;
comparisonReport.sourceDataSha256 = accepted.sourceDataSha256;

dataOutputPath = fullfile(outputRoot,"wingtip_comparison_data.mat");
jsonOutputPath = fullfile(outputRoot,"wingtip_comparison_summary.json");
plotOutputPath = fullfile(outputRoot,"wingtip_comparison.png");

save(dataOutputPath,"data","comparisonReport","-v7");
writeJson(jsonOutputPath,comparisonReport);

figureHandle = plotComparison(data,comparisonSettings);
exportgraphics(figureHandle,plotOutputPath,"Resolution", ...
    comparisonSettings.resolutionDpi);
if comparisonSettings.closeFigureAfterSave
    close(figureHandle);
end

comparisonReport.outputRoot = string(outputRoot);
comparisonReport.dataPath = string(dataOutputPath);
comparisonReport.jsonPath = string(jsonOutputPath);
comparisonReport.plotPath = string(plotOutputPath);

fprintf("\nPazy Aeroelastic Control Benchmark\n");
fprintf("SHARPy/MATLAB trim-relative wingtip comparison\n");
fprintf("  Peak ratio       : %.9f\n",metrics.peakRatioMatlabToSharpy);
fprintf("  RMS error         : %.9f m\n",metrics.rmsErrorMeters);
fprintf("  Peak error        : %.9f m\n",metrics.peakAbsoluteErrorMeters);
fprintf("  Status            : %s\n",comparisonReport.status);
fprintf("  Plot              : %s\n\n",plotOutputPath);


function settings = applyDefaults(settings)
%APPLYDEFAULTS Fill settings not supplied by the caller.

defaults = struct( ...
    "outputRoot",fullfile("results","validation", ...
        "sharpy-wingtip-comparison","runs"), ...
    "publicationMode",false, ...
    "visible",false, ...
    "closeFigureAfterSave",true, ...
    "resolutionDpi",220);

names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(settings,name) || isempty(settings.(name))
        settings.(name) = defaults.(name);
    end
end

settings.outputRoot = string(settings.outputRoot);
settings.publicationMode = logical(settings.publicationMode);
settings.visible = logical(settings.visible);
settings.closeFigureAfterSave = logical(settings.closeFigureAfterSave);
settings.resolutionDpi = double(settings.resolutionDpi);
end


function metrics = calculateMetrics(data)
%CALCULATEMETRICS Evaluate the symmetric trim-relative comparison.

sharpy = double(data.sharpy_mean_tip_m(:));
matlab = double(data.matlab_symmetric_tip_m(:));
error = matlab-sharpy;

metrics = struct();
metrics.sharpyPeakMeters = max(sharpy);
metrics.matlabPeakMeters = max(matlab);
metrics.peakRatioMatlabToSharpy = ...
    metrics.matlabPeakMeters/metrics.sharpyPeakMeters;
metrics.rmsErrorMeters = sqrt(mean(error.^2));
metrics.peakAbsoluteErrorMeters = max(abs(error));
end


function checks = validateComparison(data,metrics,accepted)
%VALIDATECOMPARISON Verify data integrity and accepted numerical metrics.

expected = accepted.metrics;
tolerance = accepted.integrityTolerance;
time = double(data.time_s(:));
numericData = data{:,vartype("numeric")};

checks = struct( ...
    "expectedSampleCount",height(data)==accepted.sampleCount, ...
    "strictlyIncreasingTime",all(diff(time)>0), ...
    "expectedTimeRange", ...
        abs(time(1)-accepted.timeRangeSeconds(1))<=tolerance.timeSeconds && ...
        abs(time(end)-accepted.timeRangeSeconds(2))<=tolerance.timeSeconds, ...
    "finiteData",all(isfinite(numericData),"all"), ...
    "sharpyPeakMatch",abs(metrics.sharpyPeakMeters- ...
        expected.sharpyPeakMeters)<=tolerance.metricAbsolute, ...
    "matlabPeakMatch",abs(metrics.matlabPeakMeters- ...
        expected.matlabPeakMeters)<=tolerance.metricAbsolute, ...
    "peakRatioMatch",abs(metrics.peakRatioMatlabToSharpy- ...
        expected.peakRatioMatlabToSharpy)<=tolerance.ratioAbsolute, ...
    "rmsErrorMatch",abs(metrics.rmsErrorMeters- ...
        expected.rmsErrorMeters)<=tolerance.metricAbsolute, ...
    "peakErrorMatch",abs(metrics.peakAbsoluteErrorMeters- ...
        expected.peakAbsoluteErrorMeters)<=tolerance.metricAbsolute);
end


function figureHandle = plotComparison(data,settings)
%PLOTCOMPARISON Plot matched symmetric wingtip displacement histories.

visibility = "off";
if settings.visible
    visibility = "on";
end

figureHandle = figure("Visible",visibility,"Color","w", ...
    "Name","SHARPy wingtip comparison");
axesHandle = axes(figureHandle);

plot(axesHandle,data.time_s,data.sharpy_mean_tip_m, ...
    "k-","LineWidth",1.45,"DisplayName","SHARPy");
hold(axesHandle,"on");
plot(axesHandle,data.time_s,data.matlab_symmetric_tip_m, ...
    "Color",[0 0.4470 0.7410],"LineWidth",1.35, ...
    "DisplayName","MATLAB reduced-order model");
xline(axesHandle,0.25,"--","Gust end", ...
    "Color",[0.45 0.45 0.45],"HandleVisibility","off");

grid(axesHandle,"on");
box(axesHandle,"on");
xlabel(axesHandle,"Time (s)");
ylabel(axesHandle,"Trim-relative symmetric wingtip displacement (m)");
legend(axesHandle,"Location","best");
if ~settings.publicationMode
    title(axesHandle,"40 m/s, 1 deg wing-only gust comparison");
end
end


function status = passFail(passed)
%PASSFAIL Return an unambiguous validation status.

status = "FAIL";
if passed
    status = "PASS";
end
end


function writeJson(path,value)
%WRITEJSON Write a human-readable UTF-8 summary.

fileId = fopen(path,"w");
assert(fileId>=0,"PazyValidation:JsonOpenFailed", ...
    "Unable to open JSON output: %s",path);
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s\n",jsonencode(value,PrettyPrint=true));
clear cleanup
end
