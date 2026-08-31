function setup = prepareProductionModelSetup(options)
%PREPAREPRODUCTIONMODELSETUP Serialize a qualified source-trim setup.
%   This path avoids recomputing a coupled nonlinear trim during ordinary
%   benchmark use. Full trim reconstruction remains part of the explicit
%   library-regeneration and validation workflow.

arguments
    options.RepositoryRoot (1,1) string
    options.SharpyRoot (1,1) string
    options.CaseName (1,1) string = "pazy_krylov_ROM"
    options.BodyCase (1,1) string {mustBeMember(options.BodyCase, ...
        ["wingOnly","coupledFull"])} = "coupledFull"
    options.SimulationMode (1,1) string {mustBeMember( ...
        options.SimulationMode,["openloop","nmhe_nmpc"])} = "openloop"
    options.GustEnabled (1,1) logical = false
    options.DurationSeconds (1,1) double = nan
    options.OperatingPoint (1,2) double {mustBeFinite} = [15,10]
    options.RunId (1,1) string
end

assert(options.OperatingPoint(1)>0, ...
    "AeroFlex:ProductionSetupSpeed","Airspeed must be positive.");
assert(strlength(options.RunId)>0 && isempty(regexp(options.RunId, ...
    "[^A-Za-z0-9_-]","once")),"AeroFlex:ProductionSetupRunId", ...
    "RunId may contain only letters, numbers, underscores, and hyphens.");

preparationTimer = tic;
[library,source] = AeroFlex.benchmark.loadProductionScheduledLibrary( ...
    options.RepositoryRoot,options.OperatingPoint);
cfg = source.cfg;
beam = source.beam;
aero = source.aero;
base = source.base;
trim = source.trim;
idx = source.p5.idx;
x0 = trim.states(:);

cfg.case_name = char(options.CaseName);
cfg.case = char(options.CaseName);
cfg.flight.U_inf = options.OperatingPoint(1);
cfg.flight.aoa_deg = options.OperatingPoint(2);
assert(isfield(source.p5,'parConst') && ...
    isfield(source.p5.parConst,'dt') && ...
    isscalar(source.p5.parConst.dt) && ...
    isfinite(source.p5.parConst.dt) && source.p5.parConst.dt>0, ...
    "AeroFlex:ProductionSetupPackageTimeStep", ...
    "The locked runtime package must own a positive finite timestep.");
cfg.sim.dt = double(source.p5.parConst.dt);
durationSeconds = options.DurationSeconds;
assert(isnan(durationSeconds) || ...
    (isfinite(durationSeconds) && durationSeconds>0), ...
    "AeroFlex:ProductionSetupDuration", ...
    "DurationSeconds must be positive and finite when specified.");
if isnan(durationSeconds)
    durationSeconds = cfg.sim.t_end;
end
cfg.sim.t_end = durationSeconds;
cfg.sim.bodyCase = lower(options.BodyCase);
cfg.sim.body_case = lower(options.BodyCase);
cfg.gust.on = options.GustEnabled;
cfg.trim.useRateProjection = options.BodyCase=="coupledFull";
% Both production body cases install their locked source package. Wing-only
% keeps the package frozen and RateProject off; only coupled-full may use the
% qualified moving schedule.
cfg.library.enable = true;
cfg.library.method = "locked-coordinate-aligned-source-library";
cfg.library.requireCompatible = true;
cfg.library.requireExactNode = true;
cfg.library.noExtrapolate = true;
cfg.library.interpTol = 1e-10;
cfg.library.schedulerQueryPoint = options.OperatingPoint;
cfg.library.updateMode = "frozenTrim";
cfg.library.scheduleAlphaMode = "trimHold";
if source.dynamicSchedulingSupported
    cfg.library.fullCoordinateRuntimeCandidate = ...
        source.fullCoordinateRuntimeCandidate;
else
    % The full-coordinate realization is qualified only on the Case-B
    % centerline. Exact wing-only/Case-A sources retain their raw locked
    % coordinates and must not enter that selector.
    cfg.library.fullCoordinateRuntimeCandidate = struct("enabled",false);
end
cfg.library.frozenPackage = source.p5;
cfg.sharedModelWorkflow = struct( ...
    "enabled",true, ...
    "entryPoint","runPazyModelWorkflow", ...
    "exactPackagePrediction",true, ...
    "actuatorCommandFrame","trim_relative_increment", ...
    "sourceId",string(source.trim.sourceId), ...
    "packageSha256",string(source.packageSha256));
cfg.operatingPoint = struct( ...
    "requestedTarget",struct("U_inf",options.OperatingPoint(1), ...
        "alpha_deg",options.OperatingPoint(2),"rho",cfg.flight.rho), ...
    "initialRuntime",struct("U_inf",options.OperatingPoint(1), ...
        "alpha_deg",options.OperatingPoint(2),"rho",cfg.flight.rho), ...
    "schedulerQuery",options.OperatingPoint, ...
    "selection",struct("mode","qualifiedSourceTrim", ...
        "exact",true,"sourceNodeName",char(trim.sourceId), ...
        "weights",1,"runtimeUpdateFrozen",true));

bodyFolder = localBodyFolder(options.BodyCase);
runRoot = fullfile(options.SharpyRoot,"sim_setup",options.CaseName, ...
    bodyFolder,options.RunId);
assert(~isfolder(runRoot),"AeroFlex:ProductionSetupExists", ...
    "The setup directory already exists: %s.",runRoot);
paths = localCreateSetupDirectories(runRoot);
cfg.paths.run_dir = runRoot;
cfg.paths.sharpy.root = options.SharpyRoot;

time = localTimeGrid(durationSeconds,cfg.sim.dt);
gust = localGust(aero,time,options.GustEnabled);
surfaceCount = localSurfaceCount(cfg,source.p5.u_eq);
deflection = zeros(numel(time),surfaceCount);
deflection_rate = zeros(numel(time),surfaceCount);
inputCount = size(aero.ROM_dsc.B,2);
gustColumn = inputCount;
input = zeros(numel(time),inputCount);
input(:,gustColumn) = gust;

sim_config = struct( ...
    "cfg",cfg,"beam",beam,"aero",aero,"base",base, ...
    "trim",trim,"idx",idx,"x0",x0,"sysd",aero.ROM_dsc, ...
    "ROMlib",library,"paths",paths,"case_name",options.CaseName, ...
    "body_case",options.BodyCase,"initialization",struct( ...
        "success",true,"trimConverged",logical(trim.converged), ...
        "policy","qualified_source_trim", ...
        "sourceId",trim.sourceId, ...
        "sourcePackage",source.packagePath, ...
        "sourcePackageSha256",source.packageSha256, ...
        "operatingPoint",cfg.operatingPoint));
state_space_system = aero.ROM_dsc;
input_settings = struct("Ts",aero.ROM_dsc.Ts,"Nc",surfaceCount, ...
    "m",inputCount,"gust_col",gustColumn);
t = time;
U = input;
% The closed-loop runner uses the measurement-array row count to recover
% the requested endpoint grid. Measurements are generated from the plant at
% runtime, so retain zero columns while preserving one row per time sample.
y = zeros(numel(time),0);
tip_deflection = [];

save(fullfile(paths.for_matlab,"sim_bundle.mat"), ...
    "sim_config","state_space_system","input_settings","t","U","y", ...
    "tip_deflection","deflection","deflection_rate","-v7");
save(fullfile(paths.for_matlab,"aero_bundle.mat"),"aero","-v7");
save(fullfile(paths.for_matlab,"beam_bundle.mat"),"beam","-v7");
save(fullfile(paths.for_matlab,"base_bundle.mat"),"base","-v7");
L = source.p5.L;
blk = source.p5.blk;
save(fullfile(paths.for_matlab,"L_bundle.mat"),"L","blk","-v7");

run_settings = struct( ...
    "sim_case",char(options.SimulationMode),"runner","PlantROM", ...
    "modelFcn","AeroFlex.sim.ROMIntegrator", ...
    "sensorFcn","AeroFlex.sensor.WingVelSensor", ...
    "estimatorFcn","AeroFlex.ctrl.nMHEv2", ...
    "controllerFcn","AeroFlex.ctrl.nMPC", ...
    "ctrl_Ts",cfg.ctrl.Ts,"Nc",cfg.ctrl.Nc,"Ne",cfg.ctrl.Ne, ...
    "x0",x0);
save(fullfile(paths.for_matlab,"run_settings.mat"), ...
    "run_settings","-v7");
roots = struct("sharpy_root",options.SharpyRoot, ...
    "sim_setup_root",fullfile(options.SharpyRoot,"sim_setup"), ...
    "sim_run_root",fullfile(options.SharpyRoot,"sim_run"), ...
    "matlab_root",fullfile(options.RepositoryRoot,"MatlabFlex"));
save(fullfile(paths.for_matlab,"roots.mat"),"roots","-v7");
save(fullfile(paths.run_dir,"roots.mat"),"roots","-v7");
localWriteExchange(paths.for_sharpy,time,gust,deflection, ...
    deflection_rate);

setup = struct("paths",paths,"cfg",cfg,"trim",trim, ...
    "initialization",sim_config.initialization, ...
    "preparationSeconds",toc(preparationTimer), ...
    "libraryPointCount",numel(library.points), ...
    "dynamicSchedulingSupported",source.dynamicSchedulingSupported);
end

function paths = localCreateSetupDirectories(runRoot)
paths = struct("run_dir",runRoot, ...
    "from_sharpy",fullfile(runRoot,"from_sharpy"), ...
    "for_matlab",fullfile(runRoot,"for_matlab"), ...
    "for_sharpy",fullfile(runRoot,"for_sharpy"), ...
    "plots",fullfile(runRoot,"plots"), ...
    "logs",fullfile(runRoot,"logs"));
for path = string(struct2cell(paths)).'
    if path~=runRoot
        mkdir(path);
    else
        mkdir(runRoot);
    end
end
end

function folder = localBodyFolder(bodyCase)
if bodyCase=="wingOnly"
    folder = "wing_only";
else
    folder = "coupled_full";
end
end

function time = localTimeGrid(duration,dt)
stepCount = round(duration/dt);
assert(abs(stepCount*dt-duration)<=100*eps(max(1,duration)), ...
    "AeroFlex:ProductionSetupGrid", ...
    "DurationSeconds must align with the source plant step %.17g s.",dt);
time = (0:stepCount).'*dt;
end

function gust = localGust(aero,time,enabled)
gust = zeros(size(time));
if ~enabled
    return
end
source = aero.gust_input(:);
if isempty(source)
    return
end
count = min(numel(source),numel(gust));
gust(1:count) = source(1:count);
if count<numel(gust)
    gust(count+1:end) = source(count);
end
end

function count = localSurfaceCount(cfg,equilibriumInput)
if isfield(cfg,"ctrl") && isfield(cfg.ctrl,"n_surf")
    count = double(cfg.ctrl.n_surf);
else
    count = numel(equilibriumInput)/2;
end
assert(isscalar(count) && isfinite(count) && count>=1 && ...
    count==round(count),"AeroFlex:ProductionSetupSurfaces", ...
    "Unable to resolve the physical surface count.");
end

function localWriteExchange(root,time,gust,deflection,deflectionRate)
path = fullfile(root,"sim_for_sharpy.h5");
h5create(path,"/time",size(time),"Datatype","double");
h5create(path,"/gust",size(gust),"Datatype","double");
h5create(path,"/deflection",size(deflection),"Datatype","double");
h5create(path,"/deflection_rate",size(deflectionRate), ...
    "Datatype","double");
h5write(path,"/time",time);
h5write(path,"/gust",gust);
h5write(path,"/deflection",deflection);
h5write(path,"/deflection_rate",deflectionRate);
end
