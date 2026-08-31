%% Pazy benchmark: standalone linear validation
% Inspect the accepted compact longitudinal model without running a nonlinear
% benchmark case. The response plots are diagnostics; formal case acceptance
% remains owned by the frozen A/B case manifests.

if ~exist("validationSettings","var") || ~isstruct(validationSettings)
    validationSettings = struct();
end

validationSettings = applyDefaults(validationSettings);
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));

setupProject(ValidateEntryPoints=true,ChangeCurrentFolder=false);

if strlength(validationSettings.artifactPath) == 0
    validationSettings.artifactPath = fullfile(repositoryRoot, ...
        "context","audits","phase18c-control-validation","lqr", ...
        "v17a-casea-corrected-model-rigid-residual-allocation-v1", ...
        "PHASE18C_V17A_CASEA_CORRECTED_NONWING_LQR_BUNDLE_V1.mat");
end

assert(isfile(validationSettings.artifactPath), ...
    "PazyValidation:MissingLinearArtifact", ...
    "Linear validation artifact not found: %s", ...
    validationSettings.artifactPath);

fprintf("\nPazy Aeroelastic Control Benchmark\n");
fprintf("Standalone longitudinal linear validation\n");
fprintf("  Artifact : %s\n",validationSettings.artifactPath);

loaded = load(validationSettings.artifactPath,"bundle");
assert(isfield(loaded,"bundle") && ...
       isfield(loaded.bundle,"directLqrComparator"), ...
    "PazyValidation:InvalidLinearArtifact", ...
    "The selected artifact does not contain the accepted LQR comparator.");

model = loaded.bundle.directLqrComparator;
sampleTime = loaded.bundle.sampleTime.outerSeconds;

A = double(model.A);
B = double(model.B);
K = double(model.gain);
Acl = A-B*K;

stateNames = string(model.stateOrder(:));
inputNames = string(model.inputOrder(:));

assert(isequal(size(A),[6 6]) && isequal(size(B),[6 3]), ...
    "PazyValidation:UnexpectedLinearDimensions", ...
    "Expected a 6-state, 3-input compact longitudinal model.");
assert(isequal(size(K),[3 6]) && all(isfinite([A(:);B(:);K(:)])), ...
    "PazyValidation:InvalidLinearModel", ...
    "The compact model or feedback gain is nonfinite or dimensionally invalid.");

%% Pole and damping summary

openPoles = eig(A);
closedPoles = eig(Acl);
recordedRadius = double(model.maximumPoleMagnitude);
closedRadius = max(abs(closedPoles));
radiusError = abs(closedRadius-recordedRadius);

equivalentPoles = log(closedPoles)/sampleTime;
naturalFrequencyHz = abs(equivalentPoles)/(2*pi);
dampingRatio = -real(equivalentPoles)./max(abs(equivalentPoles),eps);

%% Frequency response

frequencyHz = logspace(log10(0.05),log10(10),500).';
outputNames = ["airspeed perturbation","pitch attitude","pitch rate", ...
    "symmetric actuator"];
outputRows = [1 3 4 6];
C = eye(size(A,1));
C = C(outputRows,:);
D = zeros(numel(outputRows),size(B,2));

openResponse = discreteFrequencyResponse(A,B,C,D,sampleTime,frequencyHz);
closedResponse = discreteFrequencyResponse( ...
    Acl,B,C,D,sampleTime,frequencyHz);

%% Small physical input step

stepTime = (0:sampleTime:10).';
stepInput = zeros(size(B,2),1);
stepInput(1) = deg2rad(0.05);

openTrace = propagateStep(A,B,stepInput,numel(stepTime));
closedTrace = propagateStep(Acl,B,stepInput,numel(stepTime));

%% Acceptance and stored products

checks = struct( ...
    "finitePoles",all(isfinite([openPoles;closedPoles])), ...
    "stableClosedComparator",closedRadius < 1, ...
    "recordedPoleRadiusMatch",radiusError <= 1e-10, ...
    "finiteFrequencyResponse", ...
        all(isfinite([openResponse(:);closedResponse(:)])), ...
    "finiteStepResponse",all(isfinite([openTrace(:);closedTrace(:)])));

validationReport = struct();
validationReport.schemaVersion = 1;
validationReport.status = passFail(all(struct2array(checks)));
validationReport.artifactPath = string(validationSettings.artifactPath);
validationReport.sampleTimeSeconds = sampleTime;
validationReport.stateOrder = stateNames;
validationReport.inputOrder = inputNames;
validationReport.openLoopPoles = openPoles;
validationReport.closedLoopPoles = closedPoles;
validationReport.closedLoopRadius = closedRadius;
validationReport.recordedClosedLoopRadius = recordedRadius;
validationReport.closedLoopRadiusError = radiusError;
validationReport.closedNaturalFrequencyHz = naturalFrequencyHz;
validationReport.closedDampingRatio = dampingRatio;
validationReport.checks = checks;

runStamp = string(datetime("now","Format","yyyyMMdd'T'HHmmss"));
outputRoot = fullfile(repositoryRoot,validationSettings.outputRoot,runStamp);
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

dataPath = fullfile(outputRoot,"linear_validation_data.mat");
jsonPath = fullfile(outputRoot,"linear_validation_summary.json");
plotPath = fullfile(outputRoot,"linear_validation_trends.png");

save(dataPath,"validationReport","frequencyHz","openResponse", ...
    "closedResponse","stepTime","openTrace","closedTrace", ...
    "stepInput","outputNames","-v7");

jsonSummary = validationReport;
jsonSummary = rmfield(jsonSummary,["openLoopPoles","closedLoopPoles"]);
jsonSummary.openLoopPoleReal = real(openPoles);
jsonSummary.openLoopPoleImaginary = imag(openPoles);
jsonSummary.closedLoopPoleReal = real(closedPoles);
jsonSummary.closedLoopPoleImaginary = imag(closedPoles);
writeJson(jsonPath,jsonSummary);

figureHandle = plotLinearValidation(openPoles,closedPoles,frequencyHz, ...
    openResponse,closedResponse,stepTime,openTrace,closedTrace, ...
    outputRows,validationSettings.visible);
exportgraphics(figureHandle,plotPath,"Resolution",180);

if validationSettings.closeFigureAfterSave
    close(figureHandle);
end

validationReport.outputRoot = string(outputRoot);
validationReport.dataPath = string(dataPath);
validationReport.jsonPath = string(jsonPath);
validationReport.plotPath = string(plotPath);

fprintf("  Closed pole radius : %.9f\n",closedRadius);
fprintf("  Status             : %s\n",validationReport.status);
fprintf("  Plot               : %s\n\n",plotPath);


function settings = applyDefaults(settings)
%APPLYDEFAULTS Fill only values not supplied in the user-settings structure.

defaults = struct( ...
    "artifactPath","", ...
    "outputRoot",fullfile("results","validation","linear-response"), ...
    "visible",false, ...
    "closeFigureAfterSave",true);

names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(settings,name) || isempty(settings.(name))
        settings.(name) = defaults.(name);
    end
end

settings.artifactPath = string(settings.artifactPath);
settings.outputRoot = string(settings.outputRoot);
settings.visible = logical(settings.visible);
settings.closeFigureAfterSave = logical(settings.closeFigureAfterSave);
end


function status = passFail(passed)
%PASSFAIL Keep console and saved status labels unambiguous.

status = "FAIL";
if passed
    status = "PASS";
end
end


function response = discreteFrequencyResponse(A,B,C,D,dt,frequencyHz)
%DISCRETEFREQUENCYRESPONSE Evaluate C(zI-A)^(-1)B+D on the unit circle.

response = complex(zeros(size(C,1),size(B,2),numel(frequencyHz)));
identity = eye(size(A));

for index = 1:numel(frequencyHz)
    z = exp(1i*2*pi*frequencyHz(index)*dt);
    response(:,:,index) = C*((z*identity-A)\B)+D;
end
end


function trace = propagateStep(A,B,input,count)
%PROPAGATESTEP Propagate a held small-signal input from zero initial state.

trace = zeros(size(A,1),count);
for index = 1:count-1
    trace(:,index+1) = A*trace(:,index)+B*input;
end
end


function figureHandle = plotLinearValidation( ...
    openPoles,closedPoles,frequencyHz,openResponse,closedResponse, ...
    stepTime,openTrace,closedTrace,outputRows,visible)
%PLOTLINEARVALIDATION Create the deterministic publication-style summary.

visibility = "off";
if visible
    visibility = "on";
end

figureHandle = figure("Visible",visibility,"Color","w", ...
    "Name","Pazy linear validation");
layout = tiledlayout(figureHandle,2,2,"TileSpacing","compact", ...
    "Padding","compact");
title(layout,"Pazy compact longitudinal validation");

ax = nexttile(layout,1);
plot(ax,real(openPoles),imag(openPoles),"o","DisplayName","Open loop");
hold(ax,"on");
plot(ax,real(closedPoles),imag(closedPoles),"x", ...
    "LineWidth",1.3,"DisplayName","Closed loop");
angleGrid = linspace(0,2*pi,361);
plot(ax,cos(angleGrid),sin(angleGrid),"k--","DisplayName","Unit circle");
axis(ax,"equal"); grid(ax,"on");
xlabel(ax,"Real(z)"); ylabel(ax,"Imaginary(z)");
title(ax,"Discrete poles"); legend(ax,"Location","best");

ax = nexttile(layout,2);
pitchIndex = 2;
wingInput = 1;
semilogx(ax,frequencyHz,20*log10(max( ...
    squeeze(abs(openResponse(pitchIndex,wingInput,:))),eps)), ...
    "DisplayName","Open loop");
hold(ax,"on");
semilogx(ax,frequencyHz,20*log10(max( ...
    squeeze(abs(closedResponse(pitchIndex,wingInput,:))),eps)), ...
    "LineWidth",1.2,"DisplayName","Closed loop");
grid(ax,"on"); xlabel(ax,"Frequency (Hz)");
ylabel(ax,"Magnitude (dB re rad/rad)");
title(ax,"Wing-to-pitch frequency response");
legend(ax,"Location","best");

ax = nexttile(layout,3);
plot(ax,stepTime,rad2deg(closedTrace(outputRows(2),:)), ...
    "LineWidth",1.2,"DisplayName","Closed-loop pitch");
grid(ax,"on"); xlabel(ax,"Time (s)"); ylabel(ax,"Pitch (deg)");
title(ax,"0.05 deg symmetric-wing step");
legend(ax,"Location","best");

ax = nexttile(layout,4);
yyaxis(ax,"left");
plot(ax,stepTime,rad2deg(closedTrace(outputRows(3),:)), ...
    "DisplayName","Pitch rate");
ylabel(ax,"Pitch rate (deg/s)");
yyaxis(ax,"right");
plot(ax,stepTime,rad2deg(closedTrace(outputRows(4),:)), ...
    "DisplayName","Symmetric actuator");
grid(ax,"on"); xlabel(ax,"Time (s)");
ylabel(ax,"Actuator deflection (deg)");
title(ax,"Closed-loop rate and actuator response");
legend(ax,"Location","best");
end


function writeJson(path,value)
%WRITEJSON Write a human-readable UTF-8 summary without changing numerics.

fileId = fopen(path,"w");
assert(fileId >= 0,"PazyValidation:JsonOpenFailed", ...
    "Unable to open JSON output: %s",path);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,"%s\n",jsonencode(value,PrettyPrint=true));
clear cleanup
end
