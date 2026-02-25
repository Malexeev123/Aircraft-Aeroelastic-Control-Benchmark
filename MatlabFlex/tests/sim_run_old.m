function success = sim_run(case_name, body_case, varargin)
%SIM_RUN  Execute a simulation using the previously prepared sim_setup bundle.
%
% Usage (auto-pick latest setup):
%   sharpy_root = '/home/maxal/Research/TestBench';      % SHARPy root
%   matlab_root = '/home/maxal/Research/MatlabFlex';     % MATLAB project root
%   success = sim_run(sharpy_root, matlab_root, 'pazy_ROM', 'wingOnly');
%
% Usage (explicit setup dir):
%   success = sim_run(sharpy_root, matlab_root, 'pazy_ROM','wingOnly', ...
%                     'setup_dir', '/.../sim_setup/pazy_ROM/wing_only/20251110');
%
% Name-Value:
%   'setup_dir'       : explicit sim_setup run directory to read from
%   'date_only_runs'  : true  -> sim_run/<date>, false -> <date_time> (default: true)
%   'overwrite'       : whether to reuse an existing run dir (default: true)
%
% Notes
%   • Reads sim_bundle.mat, run_settings.mat, base/beam/aero bundles, sim_for_sharpy.h5
%   • Supports runner families: 'SimRunner' and 'PlantROM'
%   • Supports sim_case: 'openloop' and 'nmhe_nmpc'
%   • Minimal outputs now; post-processing will follow later

% ------------------- parse args -------------------
p = inputParser;
p.addParameter('setup_dir','',@(s)ischar(s)||isstring(s));
p.addParameter('date_only_runs', true, @(b)islogical(b)||isnumeric(b));
p.addParameter('overwrite',      true, @(b)islogical(b)||isnumeric(b));
p.parse(varargin{:});
setup_dir      = char(p.Results.setup_dir);
date_only_runs = logical(p.Results.date_only_runs);
overwrite      = logical(p.Results.overwrite);

% ---- locate setup_dir and roots without explicit args ----
% Option A: user provided an explicit setup_dir via name-value (honored as-is)
if isempty(setup_dir)
    % Option B: find latest sim_setup by using a saved roots.mat (preferred)
    % Try the most recent sim_setup/<case>/<body>/... by scanning roots.sim_setup_root
    roots = [];
    % 1) If a previous run left a roots.mat in the working dir, load it
    if exist('roots.mat','file')==2
        Sroots = load('roots.mat'); roots = Sroots.roots;
    end
    % 2) If not, try to load roots from the most recent sim_setup for this case/body
    if isempty(roots)
        % We can guess roots by looking up any sim_setup that contains our case/body:
        % If there’s a recent setup, it has for_matlab/roots.mat — search upward from current folder
        % (fast path) if you usually call sim_run from MatlabFlex/configs:
        here = fileparts(mfilename('fullpath'));   % .../MatlabFlex/configs
        research_root = fileparts(here);           % .../MatlabFlex
        % Try sibling 'TestBench/sim_setup' if that’s your layout; otherwise we’ll fall back below.
        candidate_sharpy = fullfile(fileparts(research_root),'TestBench');
        candidate_sim_setup = fullfile(candidate_sharpy,'sim_setup');

        if exist(candidate_sim_setup,'dir')==7
            setup_dir = find_latest_setup(candidate_sim_setup, case_name, sanitize_bodycase(body_case));
            if ~isempty(setup_dir)
                fm = fullfile(setup_dir,'for_matlab');
                if exist(fullfile(fm,'roots.mat'),'file')==2
                    Sroots = load(fullfile(fm,'roots.mat')); roots = Sroots.roots;
                end
            end
        end
    end
    % 3) If still empty, last fallback: we will derive roots from sim_bundle once setup_dir is known
    if isempty(setup_dir)
        % If we already found roots, we can now find the latest setup cleanly
        if ~isempty(roots) && isfield(roots,'sim_setup_root')
            setup_dir = find_latest_setup(roots.sim_setup_root, case_name, sanitize_bodycase(body_case));
        end
    end
end

% If setup_dir still empty, error out with a helpful message
assert(~isempty(setup_dir), 'Could not locate a sim_setup directory. Consider passing ''setup_dir'', or ensure roots.mat exists from sim_init.');

% Load sim_bundle to derive roots if we didn’t have them yet
fm_setup = fullfile(setup_dir,'for_matlab');
assert(exist(fullfile(fm_setup,'sim_bundle.mat'),'file')==2, 'sim_bundle.mat not found in setup for_matlab.');

Ssim_setup = load(fullfile(fm_setup,'sim_bundle.mat'));    % has sim_config.paths.run_dir
if isempty(roots)
    % Reconstruct sharpy_root from run_dir = <sharpy_root>/sim_setup/<case>/<body>/<date>
    run_dir = Ssim_setup.sim_config.paths.run_dir;
    p1 = fileparts(run_dir);    % .../sim_setup/<case>/<body>
    p2 = fileparts(p1);         % .../sim_setup/<case>
    p3 = fileparts(p2);         % .../sim_setup
    sharpy_root = fileparts(p3);
    roots.sharpy_root    = sharpy_root;
    roots.sim_setup_root = p3;
    roots.sim_run_root   = fullfile(sharpy_root,'sim_run');

    % matlab_root from this file’s location
    roots.matlab_root    = fileparts(fileparts(mfilename('fullpath')));  % .../MatlabFlex
else
    sharpy_root = roots.sharpy_root;
end

% With roots determined, create sim_run/<case>/<body>/<date>
sim_run_root = fullfile(roots.sharpy_root,'sim_run');
[run_dir, paths] = create_run_dirs(sim_run_root, case_name, sanitize_bodycase(body_case), overwrite, date_only_runs);

% Put project on path
addpath(genpath(roots.matlab_root));

% 
% 
% % ------------------- paths & MATLAB path -------------------
% addpath(genpath(fullfile(matlab_root)));  % project code
% sim_setup_root = fullfile(sharpy_root, 'sim_setup');
% sim_run_root   = fullfile(sharpy_root, 'sim_run');
% 
% body_case_dir  = sanitize_bodycase(body_case);
% 
% if isempty(setup_dir)
%     setup_dir = find_latest_setup(sim_setup_root, case_name, body_case_dir);
%     assert(~isempty(setup_dir), 'No sim_setup folder found for %s/%s.', case_name, body_case_dir);
% end
% 
% % Create run folder under sim_run/<case>/<body>/<date>
% [run_dir, paths] = create_run_dirs(sim_run_root, case_name, body_case_dir, overwrite, date_only_runs);

% ------------------- global console log -------------------
if ~exist(paths.logs,'dir'), mkdir(paths.logs); end
log_name = sprintf('run_%s_%s_%s.log', case_name, body_case_dir, datestr(now,'yyyymmdd'));
console_log = fullfile(paths.logs, log_name);
diary off; diary(console_log); diary on;
cleanupLog = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('=== sim_run started %s ===\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf('setup_dir: %s\n', setup_dir);
fprintf('run_dir  : %s\n', run_dir);

% ------------------- load setup artifacts -------------------
fm = fullfile(setup_dir,'for_matlab');
fs = fullfile(setup_dir,'for_sharpy');
assert(exist(fm,'dir')==7 && exist(fs,'dir')==7, 'Setup dir missing for_matlab/for_sharpy.');

% Core bundles
Ssim = load(fullfile(fm,'sim_bundle.mat'));         % sim_config, state_space_system, input_settings, ...
Saer = load(fullfile(fm,'aero_bundle.mat'));        % aero
Sbem = load(fullfile(fm,'beam_bundle.mat'));        % beam
Sbas = load(fullfile(fm,'base_bundle.mat'));        % base
Srun = load(fullfile(fm,'run_settings.mat'));       % run_settings

% Optional L
if exist(fullfile(fm,'L_bundle.mat'),'file')
    SL = load(fullfile(fm,'L_bundle.mat')); L = SL.L; blk = SL.blk; %#ok<NASGU>
else
    L = []; blk = [];
end

% Gust/time H5
h5 = fullfile(fs,'sim_for_sharpy.h5');
assert(exist(h5,'file')==2, 'Missing sim_for_sharpy.h5 in setup.');
t_vec       = h5read(h5,'/time');                   % nsamp×1
gust_series = h5read(h5,'/gust');                   % nsamp×1
deflection  = h5read(h5,'/deflection');             % nsamp×Nc
defl_rate   = h5read(h5,'/deflection_rate');        % nsamp×Nc

% Rehydrate objects/config
cfg  = Ssim.sim_config.cfg;          %#ok<NASGU>
aero = Saer.aero;
beam = Sbem.beam;
base = Sbas.base;
trim = Ssim.sim_config.trim;
idx  = Ssim.sim_config.idx;
x0   = Ssim.sim_config.x0;

% Make sure gust pipe is consistent for SimRunner path
if isprop(aero, 'gust_input')
    aero.gust_input = gust_series(:);   % guarantee column
end

% Runner preferences
run_settings = Srun.run_settings;
sim_case     = lower(char(run_settings.sim_case));
runner       = lower(char(run_settings.runner));

% ------------------- prep outputs & metadata -------------------
copyfile(setup_dir, fullfile(paths.run_dir,'_setup_snapshot'));  % helpful provenance

% Save a small meta for the run
meta_fn = fullfile(paths.run_dir,'meta_run.yaml');
fid = fopen(meta_fn,'w');
if fid~=-1
    fprintf(fid,'case_name: %s\n', Ssim.sim_config.case_name);
    fprintf(fid,'body_case: %s\n', body_case_dir);
    fprintf(fid,'setup_dir: %s\n', setup_dir);
    fprintf(fid,'created:   %s\n', datestr(now,'yyyy-mm-dd'));
    fprintf(fid,'sim_case:  %s\n', sim_case);
    fprintf(fid,'runner:    %s\n', runner);
    fclose(fid);
end

% ------------------- RUN DISPATCH ----------------------------------------
results = struct(); hadError = false; MEout = [];

try
    switch sim_case

        % ================= OPENLOOP =================
        case 'openloop'
            switch runner
                case 'simrunner'
                    % Fresh SimRunner (don’t reuse any stale handles from MAT)
                    simObj = AeroFlex.sim.SimRunner(cfg, beam, aero, base, trim);
                    cfg.sim.storeSens = true;               % mirror setup
                    u_ctrl = zeros(cfg.ctrl.n_surf*cfg.ctrl.var_per,1); % open-loop
                    [t, x, log, sensEq] = simObj.run(x0, Ssim.sim_config.cfg.sim.t_end, u_ctrl, cfg.sim.storeSens);

                case 'plantrom'
                    % Plant integrator path (ROMIntegrator + PlantRunTime)
                    modelFcn      = str2func(run_settings.modelFcn);
                    sensorFcn     = str2func(run_settings.sensorFcn);
                    estimatorFcn  = str2func(run_settings.estimatorFcn);
                    controllerFcn = @(cfg_in) []; %#ok<NASGU>   % not used in open-loop

                    % Build model; many projects pass cfg, beam, aero, base
                    rom = modelFcn(cfg, beam, aero, base);
                    plant = AeroFlex.sim.PlantRunTime(cfg, beam, aero, base, x0);
                    % Open-loop propagate with gust only (no control)
                    Ts = run_settings.ctrl_Ts;
                    Nt = floor(Ssim.sim_config.cfg.sim.t_end / Ts);
                    t = (0:Nt)*Ts; x = zeros(length(x0), Nt+1); x(:,1) = x0;
                    log = struct(); sensEq = [];
                    for k = 1:Nt
                        gk = gust_series( min(k+1, numel(gust_series)) );
                        uk = zeros(cfg.ctrl.n_surf*cfg.ctrl.var_per, 1);
                        x(:,k+1) = plant.step(uk, gk, "ROM");
                    end

                otherwise
                    error('Unknown runner="%s".', runner);
            end

        % ================= NMHE/NMPC =================
        case 'nmhe_nmpc'
            switch runner
                case 'simrunner'
                    % Let SimRunner own the inner loop (if your SimRunner supports control hooks)
                    % Here we just call it; user’s SimRunner implementation applies control inside.
                    simObj = AeroFlex.sim.SimRunner(cfg, beam, aero, base, trim);
                    cfg.sim.storeSens = false;
                    [t, x, log, sensEq] = simObj.run(x0, Ssim.sim_config.cfg.sim.t_end, cfg.sim.storeSens);

                case 'plantrom'
                    % Full loop: Plant + Estimator + Controller (like your example)
                    modelFcn      = str2func(run_settings.modelFcn);
                    sensorFcn     = str2func(run_settings.sensorFcn);
                    estimatorFcn  = str2func(run_settings.estimatorFcn);
                    controllerFcn = str2func(run_settings.controllerFcn);

                    % Instantiate
                    rom    = modelFcn(cfg, beam, aero, base);                     %#ok<NASU> 
                    sensor = sensorFcn(beam, cfg);                                %#ok<NASU>
                    plant  = AeroFlex.sim.PlantRunTime(cfg, beam, aero, base, x0);
                    est    = estimatorFcn(cfg);
                    ctrl   = controllerFcn(cfg);

                    % Time grid
                    Ts = run_settings.ctrl_Ts;
                    Nt = floor(Ssim.sim_config.cfg.sim.t_end / Ts);
                    t  = (0:Nt-1)*Ts;

                    % Logs
                    xhat_hist = zeros(length(x0), Nt);
                    Uhist     = zeros(cfg.ctrl.n_surf*cfg.ctrl.var_per, Nt);
                    wHorizon  = zeros(run_settings.Ne, Nt);
                    trueGust  = zeros(1, Nt);

                    xk = x0; u_prev = zeros(cfg.ctrl.n_surf*cfg.ctrl.var_per,1);
                    for k = 1:Nt
                        % Measurements (your PlantRunTime API returns z_k, t_k)
                        [z_k, t_k] = plant.readSensors(cfg, []);   %#ok<ASGLU>

                        % True gust from prepared profile
                        gk = gust_series( min(round(t_k/cfg.sim.dt)+1, numel(gust_series)) );
                        trueGust(k) = gk;

                        % Estimate
                        [xhat, what, ~] = est.estimate(z_k, u_prev, t_k);

                        % Control
                        [uk, Uinfo] = ctrl.computeControl(xhat, gk, u_prev, t_k); %#ok<ASGLU>
                        u_prev = uk;

                        % Propagate plant
                        xk = plant.step(uk, gk, "ROM");

                        % Log
                        xhat_hist(:,k) = xhat;
                        Uhist(:,k)     = uk;
                        if isnumeric(what)
                            wHorizon(1:min(end,numel(what)),k) = what(1:min(end,numel(what)));
                        end
                    end

                    % Pack outputs to match open-loop fields
                    x = []; sensEq = [];
                    log = struct('xhat',xhat_hist,'U',Uhist,'wHorizon',wHorizon,'wTrue',trueGust);

                otherwise
                    error('Unknown runner="%s".', runner);
            end

        otherwise
            error('Unknown sim_case="%s".', sim_case);
    end

    % Store results
    results.t = exist('t','var') * t;             %#ok<NASGU,NBRAK>
    results.x = exist('x','var') * x;             %#ok<NASGU,NBRAK>
    results.log = exist('log','var') * log;       %#ok<NASGU,NBRAK>
    results.sensEq = exist('sensEq','var') * sensEq; %#ok<NASGU,NBRAK>

catch ME
    hadError = true; MEout = ME;
    fprintf(2,'[sim_run] ERROR: %s\n', ME.message);
    fprintf(2,'%s\n', getReport(ME,'extended','hyperlinks','off'));
end

% ------------------- Save run bundle -------------------
try
    save(fullfile(paths.for_matlab,'run_bundle.mat'), ...
         'results','cfg','aero','beam','base','trim','idx','run_settings','-v7.3');
catch ME
    warning('[sim_run] Could not save run_bundle: %s', ME.message);
end

% Also drop a minimal H5 with time & (if present) a primary state trace
try
    h5out = fullfile(paths.for_sharpy,'sim_run_out.h5');
    if exist(h5out,'file'), delete(h5out); end
    if exist('t','var') && ~isempty(t)
        h5create(h5out,'/time',[numel(t) 1]); h5write(h5out,'/time',t(:));
    end
    if exist('x','var') && ~isempty(x)
        h5create(h5out,'/x_shape',[2 1]); h5write(h5out,'/x_shape',[size(x,1); size(x,2)]);
        % write first modal velocity row if it exists
        row = min(size(x,1), 1);
        h5create(h5out,'/x_row1',[1 size(x,2)]); h5write(h5out,'/x_row1',x(row,:));
    end
catch ME
    warning('[sim_run] H5 export skipped: %s', ME.message);
end

% ------------------- finalize log -------------------
fprintf('=== sim_run finished %s (status: %s) ===\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), ternary(~hadError,'OK','ERROR'));
try, copyfile(console_log, fullfile(paths.run_dir,'console.log')); catch, end

success = ~hadError;
if hadError
    % rethrow(MEout);  % uncomment if you want hard failure
end
end
% --------------- helpers ----------------
function [run_dir, paths] = create_run_dirs(sim_root, case_name, body_case_dir, overwrite, date_only_runs)
    if ~exist(sim_root,'dir'), mkdir(sim_root); end
    case_dir = fullfile(sim_root, case_name); if ~exist(case_dir,'dir'), mkdir(case_dir); end
    body_dir = fullfile(case_dir, body_case_dir); if ~exist(body_dir,'dir'), mkdir(body_dir); end
    if date_only_runs
        stamp = datestr(now,'yyyymmdd');
    else
        stamp = datestr(now,'yyyymmdd_HHMMSS');
    end
    run_dir = fullfile(body_dir, stamp);
    if exist(run_dir,'dir')
        if ~overwrite, error('Run dir exists: %s', run_dir); end
    else
        mkdir(run_dir);
    end
    paths.run_dir     = run_dir;
    paths.from_sharpy = must_mkdir(fullfile(run_dir,'from_sharpy'));
    paths.for_matlab  = must_mkdir(fullfile(run_dir,'for_matlab'));
    paths.for_sharpy  = must_mkdir(fullfile(run_dir,'for_sharpy'));
    paths.plots       = must_mkdir(fullfile(run_dir,'plots'));
    paths.logs        = must_mkdir(fullfile(run_dir,'logs'));
end

function d = must_mkdir(d), if ~exist(d,'dir'), mkdir(d); end, end

function s = sanitize_bodycase(s)
    s = lower(char(s)); s = strrep(s,' ','_');
    if strcmp(s,'wingonly'),    s = 'wing_only';    end
    if strcmp(s,'coupledfull'), s = 'coupled_full'; end
end

function setup_dir = find_latest_setup(sim_setup_root, case_name, body_case_dir)
    base = fullfile(sim_setup_root, case_name, body_case_dir);
    if exist(base,'dir')~=7, setup_dir = ''; return; end
    dd = dir(base); dd = dd([dd.isdir]);
    % keep only folders named yyyymmdd or yyyymmdd_HHMMSS
    names = setdiff({dd.name},{'.','..'});
    if isempty(names), setup_dir=''; return; end
    % sort by datenum
    [~,ix] = sort([dd.datenum],'descend');
    setup_dir = fullfile(base, dd(ix(1)).name);
end

function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end
