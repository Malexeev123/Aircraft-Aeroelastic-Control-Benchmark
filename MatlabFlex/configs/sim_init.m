function [sim_config, success] = sim_init(sharpy_root_in, varargin)
%SIM_INIT  Prepare a timestamped sim_setup under the SHARPy tree, load ROM,
%          merge config with defaults, run open-loop test, and export
%          artifacts for MATLAB and SHARPy HIL.
%
% FILE LOCATION (assumed):
%   Research/MatlabFlex/configs/sim_init.m    <-- this file
% DIRECTORY LAYOUT:
%   Research/
%   ├─ MatlabFlex/    (MATLAB project/sources; this file lives here)
%   └─ TestBench/     (SHARPy cases/output/<case>/...  ==> sharpy_root)
%
% REQUIRED POSITIONAL ARG
%   sharpy_root_in : path to SHARPy "TestBench" folder (root of cases/output)
%
% NAME-VALUE OPTIONS
%   'case_name'     : default 'pazy_ROM'
%   'sim_root'      : default '<sharpy_root>/sim_setup'    % << you asked this
%   'debug'         : logical, default true
%   'overwrite'     : logical, default false
%   'research_root' : optional, defaults to parent of matlab_root
%   'matlab_root'   : optional, defaults to folder auto-detected from this file
%
% OUTPUTS
%   sim_config : struct with cfg/beam/aero/base/trim/idx/x0 and metadata
%   success    : true on success
%
% FOLDER LAYOUT CREATED
%   sim_root/<case>/<YYYYMMDD_HHMMSS>/
%     ├─ from_sharpy/    (copy of linss.h5 & sim params if available)
%     ├─ for_matlab/     (MAT bundle for your MATLAB controllers/sims)
%     ├─ for_sharpy/     (CSV/HDF5 time-series to drive SHARPy plant)
%     ├─ plots/          (open-loop figures)
%     └─ logs/           (diags / text logs)

%% ------------------------------------------------------------
%  0) Parse inputs and set flags
% -------------------------------------------------------------
% auto-detect matlab_root from this file's location (.../MatlabFlex/configs)
this_file = mfilename('fullpath');
this_dir  = fileparts(this_file);        % .../MatlabFlex/configs
matlab_root_autodetect = fileparts(this_dir);  % .../MatlabFlex

p = inputParser;
p.addRequired('sharpy_root_in', @(s)ischar(s) || isstring(s));
p.addParameter('case_name', 'pazy_krylov_ROM', @(s)ischar(s) || isstring(s));
p.addParameter('sim_root', '', @(s)ischar(s) || isstring(s));
p.addParameter('debug', true, @(b)islogical(b) || isnumeric(b));
p.addParameter('overwrite', false, @(b)islogical(b) || isnumeric(b));
p.addParameter('research_root','', @(s)ischar(s) || isstring(s));
p.addParameter('matlab_root','',   @(s)ischar(s) || isstring(s));
p.addParameter('body_case', 'wingOnly', @(s)ischar(s) || isstring(s));
p.addParameter('sim_case', 'openloop', @(s)ischar(s) || isstring(s));
% p.addParameter('runner', 'SimRunner', @(s)ischar(s) || isstring(s));
p.addParameter('runner', 'PlantROM', @(s)ischar(s) || isstring(s));
p.parse(sharpy_root_in, varargin{:});
opt = p.Results;

% Resolve roots (SHARPy is required positional)
sharpy_root = strip_trailing_filesep(char(opt.sharpy_root_in));
if isempty(opt.matlab_root)
    matlab_root = matlab_root_autodetect;
else
    matlab_root = strip_trailing_filesep(char(opt.matlab_root));
end
if isempty(opt.research_root)
    research_root = fileparts(matlab_root);  % parent of MatlabFlex
else
    research_root = strip_trailing_filesep(char(opt.research_root));
end
% Default sim_root: ***UNDER SHARPY ROOT*** (your request)
if isempty(opt.sim_root)
    sim_root = fullfile(sharpy_root, 'sim_setup');
else
    sim_root = strip_trailing_filesep(char(opt.sim_root));
end

case_name = char(opt.case_name);
debug     = logical(opt.debug);
overwrite = logical(opt.overwrite);

% Put MatlabFlex project on path
if exist(matlab_root,'dir')
    addpath(genpath(matlab_root));
    if debug, fprintf('[sim_init] Added MatlabFlex to path: %s\n', matlab_root); end
else
    warning('[sim_init] matlab_root not found: %s', matlab_root);
end

% Try to add SHARPy MATLAB helpers (if present)
try
    addpath(genpath(fullfile(sharpy_root, 'example_notebooks', 'UDP_control', 'matlab_functions')));
    if debug, fprintf('[sim_init] Added SHARPy MATLAB helpers.\n'); end
catch
    % ok if absent — fallback to raw HDF5 reader
end

fprintf('[sim_init] Roots:\n  research_root = %s\n  matlab_root   = %s\n  sharpy_root   = %s\n', ...
        research_root, matlab_root, sharpy_root);
fprintf('[sim_init] sim_root (write target) = %s\n', sim_root);
fprintf('[sim_init] CWD: %s\n', pwd);

roots = struct();
roots.matlab_root   = matlab_root;
roots.sharpy_root   = sharpy_root;                     % you already have this in sim_init
roots.sim_setup_root= fullfile(sharpy_root,'sim_setup');
roots.sim_run_root  = fullfile(sharpy_root,'sim_run');



success = false;

cfg.case_name = case_name;  % make it explicit
cfg.NetworkPath = getOr(cfg,'NetworkPath',0); % preserve if already set

cfg.paths.sharpy.root         = sharpy_root;
cfg.paths.sharpy.case_root    = fullfile(sharpy_root, 'output', case_name);
cfg.paths.sharpy.savedata     = fullfile(sharpy_root, 'output', case_name, 'savedata');
cfg.paths.sharpy.linear       = fullfile(sharpy_root, 'output', case_name, 'linear_results');
body_case_str = char(opt.body_case);


% (Optional) sanity check, fail early if savedata missing
if ~exist(cfg.paths.sharpy.savedata,'dir')
    error('SHARPy savedata folder not found: %s', cfg.paths.sharpy.savedata);
end



%% ------------------------------------------------------------
%  1) Create timestamped run folder ***UNDER SHARPY ROOT***
% -------------------------------------------------------------
[run_dir, paths] = create_run_dirs(sim_root, case_name, body_case_str, overwrite);
fprintf('[sim_init] Run folder: %s\n', run_dir);

% Also point at canonical SHARPy output/linear_results (for mirror writes)
sharpy_output_linear = fullfile(sharpy_root, 'output', case_name, 'linear_results');
if ~exist(sharpy_output_linear,'dir')
    % mkdir(sharpy_output_linear);  % ensure exists for txt mirrors
    error('Cannot find SHARPy linear_results file for case: %s', case_name)
end

% Save a copy next to the other MATLAB bundles (and at run root for convenience)
% save(fullfile(paths.for_matlab,'roots.mat'),'roots','-v7.3');
save(fullfile(paths.for_matlab,'roots.mat'),'roots');
% save(fullfile(paths.run_dir,'roots.mat'),'roots','-v7.3');
save(fullfile(paths.run_dir,'roots.mat'),'roots');

%% ===== Console log capture (start) =======================================
% One log per run; captures ALL Command Window output (fprintf, warnings, errors).
if ~exist(paths.logs,'dir'), mkdir(paths.logs); end
log_name    = sprintf('sim_%s_%s.log', case_name, sanitize_bodycase(body_case_str),datestr(now,'yyyymmdd'));
console_log = fullfile(paths.logs, log_name);

% If any prior diary is active, turn it off to avoid switching files mid-run.
diary off;
diary(console_log);
diary on;

% Ensure we always close the diary even if an error occurs.
diary_cleanup__do_not_delete = onCleanup(@() diary('off'));

% Nice header
fprintf('=== sim_init started at %s ===\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf('Run folder : %s\n', paths.run_dir);
fprintf('Case       : %s\n', case_name);
% fprintf('Body-case  : %s\n', exist('opt','var') && isfield(opt,'body_case') ? opt.body_case : 'wingOnly');
% Decide body-case string (robust even if opt doesn't exist/is empty)
if exist('opt','var') && isstruct(opt) && isfield(opt,'body_case') && ~isempty(opt.body_case)
    body_case_str = char(opt.body_case);   % ensure char for fprintf
else
    body_case_str = 'wingOnly';
end
fprintf('Body-case  : %s\n', body_case_str);

fprintf('MATLAB     : %s (PID %d)\n', version, feature('getpid'));
warning('on','backtrace');  % include backtraces in the log


hadError = false; MEout = [];
try

%% ------------------------------------------------------------
%  2) Load configs (SHARPy simulation_parameters.mat if available),
%     then merge with nominal defaults (defaults fill missing only)
% -------------------------------------------------------------
cfg_default = nominalConfig();   % your existing default config

sharpy_param_mat = fullfile(sharpy_root, 'output', case_name, 'linear_results', 'simulation_parameters.mat');
cfg_loaded = struct();
if exist(sharpy_param_mat, 'file')
    S = load(sharpy_param_mat);
    % disp(fieldnames(S))

    cfg_loaded = map_sharpy_sim_params_to_cfg(S);
    cfg_loaded.flight.aoa_deg = double(S.aoa_deg);
    cfg_loaded.flight.U_inf = S.u_inf;
    % disp(cfg_loaded.flight.U_inf);
    % disp(fieldnames(cfg_loaded))
    fprintf('[sim_init] Loaded SHARPy simulation_parameters.mat\n');
    copyfile(sharpy_param_mat, fullfile(paths.from_sharpy, 'simulation_parameters.mat'));
else
    fprintf('[sim_init] No SHARPy simulation_parameters.mat found. Proceeding with defaults+local.\n');
end

% cfg.flight.aoa_deg = cfg.aoa_deg;
cfg = merge_cfg(cfg, cfg_default);

cfg = merge_cfg(cfg, cfg_loaded);
% q_inf = 0.5 * cfg.flight.rho * cfg.flight.U_inf^2; 
%     % disp(cfg.flight.U_inf)
%     % disp(q_inf);
% cfg.flight.aeroscale = q_inf*cfg.flight.b_ref^2; 
% cfg.flight.t_scale = cfg.flight.b_ref/cfg.flight.U_inf; 
% cfg.flight.Fscale = sqrt(cfg.ctrl.Ts/(cfg.flight.aeroscale*cfg.flight.t_scale)); 
% disp(cfg.flight.Fscale)
cfg = ensure_nominal_extras(cfg);
validate_cfg(cfg);
if debug, show_cfg_summary(cfg); end


    cfg.flight.b_ref = 0.05;  % [m] reference semi‑chord (0.1 / 2)
    % cfg.flight.b_ref = 0.1;  % [m] reference chord 
    % cfg.flight.b_ref = .55/2;  % [m] reference semi‑chord (0.1 / 4)
    % cfg.flight.b_ref = .55;  % [m] reference semi‑chord (0.1 / 4)
    % cfg.flight.aoa_deg = 1; % in degrees
    % cfg.flight.S_ref = cfg.struct.L *cfg.flight.b_ref; % [m^2] Since it is rectangular b*c

q_inf = 0.5*cfg.flight.rho*cfg.flight.U_inf^2;
% cfg.flight.a = q_inf*cfg.flight.b_ref^2*cfg.sim.dt;
cfg.flight.a = q_inf*cfg.flight.b_ref^2;
% cfg.flight.a = q_inf*(cfg.flight.b_ref/2)^2;
cfg.flight.t_inf = cfg.flight.b_ref/cfg.flight.U_inf;
% cfg.flight.t_inf = cfg.flight.b_ref/cfg.flight.U_inf/cfg.sim.dt;
% cfg.flight.t_inf = 2*cfg.flight.b_ref/cfg.flight.U_inf;
% % 
cfg.force_map.scale =  cfg.flight.a;
% disp(cfg.flight.a);
cfg.t_map =  cfg.flight.t_inf;

% cfg.force_map.scale =  1;
% cfg.t_map =  1;

% cfg.flight.a = 1/cfg.flight.a;
% cfg.flight.t_inf =1/cfg.flight.t_inf;
% cfg.flight.t_inf =1/cfg.flight.t_inf*cfg.sim.dt;

cfg.flight.a = 1;
cfg.flight.t_inf = 1;


% disp(cfg.flight.aoa_deg)

%% A. Pazy Wing Init

%% ------------------------------------------------------------
%  3) Aerodynamic ROM
% -------------------------------------------------------------

% Build the aerodynamic ROM from SHARPy files:
%   - cases/<case>/<case>.aero.h5
%   - output/<case>/savedata/<case>.linss.h5
%   - output/<case>/savedata/<case>.data.h5
%   - output/<case>/savedata/*.rom.h5        (Krylov ROM, discrete)
%   - output/<case>/savedata/*aerorob*.h5    (projection V, optional)
fprintf('[sim_init] Building AeroROM from SHARPy savedata...\n');
aero = AeroFlex.aero.AeroROM(cfg);

cfg.ctrl.Ts = aero.ROM_dsc.Ts;
cfg.sim.dt = aero.ROM_dsc.Ts;
q_inf = 0.5 * cfg.flight.rho * cfg.flight.U_inf^2; 
    % disp(cfg.flight.U_inf)
    % disp(q_inf);
cfg.flight.aeroscale = q_inf*cfg.flight.b_ref^2; 
cfg.flight.t_scale = cfg.flight.b_ref/cfg.flight.U_inf; 
% cfg.flight.t_scale = 1;
cfg.flight.Fscale = 1; 
% cfg.flight.Fscale = sqrt(cfg.ctrl.Ts/(cfg.flight.aeroscale*cfg.flight.t_scale)); 
% cfg.flight.Fscale = 1/(cfg.flight.aeroscale); 
cfg.flight.Fscale = (cfg.flight.aeroscale); 
% disp(cfg.flight.Fscale)

% Quick diagnostics
n = size(aero.ROM_dsc.A,1); m = size(aero.ROM_dsc.B,2); p = size(aero.ROM_dsc.C,1);
fprintf('[sim_init] AeroROM OK: n=%d, m=%d, p=%d, Ts=%.6g\n', n, m, p, aero.ROM_dsc.Ts);
% disp(fieldnames(aero))
% Persist a small bundle for MATLAB-side inspection/reuse
save(fullfile(paths.for_matlab,'aero_bundle.mat'), 'aero');
cfg.Ts = aero.ROM_dsc.Ts;
% Plot a minimal aero-only impulse check to verify stability
if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'debugPlots') && cfg.plotOpts.debugPlots
    try
        sysd_a = aero.ROM_dsc;
        t_imp = (0:400)'*sysd_a.Ts;              % ~ few hundred steps
        u_imp = zeros(numel(t_imp), size(sysd_a.B,2)); 
        u_imp(1,end) = 1;                        % impulse on last input (often gust)
        y_imp = lsim(sysd_a, u_imp, t_imp);
        fA = figure('Color','w','Position',[100 100 900 360]);
        plot(t_imp, y_imp(:,1), 'LineWidth', 1.25); grid on
        xlabel('t [s]'); ylabel('y_1'); title('AeroROM discrete impulse (y_1)')
        saveas(fA, fullfile(paths.plots,'aero_rom_impulse_y1.png'));
        close(fA);
    catch ME
        warning('[sim_init] AeroROM impulse plot skipped: %s', ME.message);
    end
end
% S = load(fullfile(paths.for_matlab,'aero_bundle.mat'),'aero'); A=S.aero; fprintf('AeroROM OK: n=%d, m=%d, p=%d, Ts=%.6g\n',size(A.ROM_dsc.A,1),size(A.ROM_dsc.B,2),size(A.ROM_dsc.C,1),A.ROM_dsc.Ts);
% disp(fieldnames(S.aero))


%% ===== STRUCTURAL BEAM (SHARPy-aware) ====================================
% Required cfg fields for BeamModel (set earlier in sim_init):
%   cfg.case_name
%   cfg.paths.sharpy.root
%   cfg.paths.sharpy.cases     = <sharpy_root>/cases/<case>
%   cfg.paths.sharpy.savedata  = <sharpy_root>/output/<case>/savedata
%   cfg.paths.sharpy.beam_modal= <sharpy_root>/output/<case>/beam_modal_analysis  

fprintf('[sim_init] Building BeamModel from SHARPy cases/savedata...\n');
beam = AeroFlex.beam.BeamModel(cfg);
% Vkeep = load("eigenvectors.dat");
% beam.phi0 = Vkeep;
% Quick diagnostics
fprintf('[sim_init] BeamModel OK: nDofs=%d, nFlex=%d, Nm=%d, first ω=%.3f rad/s\n', ...
        beam.nDofs, beam.nFlex, beam.Nm, beam.Omega(1,1));

% Test
% beam.Gamma1 = zeros(beam.Nm, beam.Nm,beam.Nm);
% beam.Gamma2 = zeros(beam.Nm, beam.Nm,beam.Nm);

% Persist for MATLAB-side use
save(fullfile(paths.for_matlab,'beam_bundle.mat'), 'beam');
% S = load(fullfile(paths.for_matlab,'beam_bundle.mat'),'beam'); B=S.beam;
% fprintf('Beam φ0 size = %dx%d, Ω(1)=%.3f rad/s\n', size(B.phi0,1), size(B.phi0,2), B.Omega(1,1));
% disp(fieldnames(S.beam))

%% Test
% beam.Gamma1 = zeros(beam.Nm, beam.Nm,beam.Nm);
% beam.Gamma2 = zeros(beam.Nm, beam.Nm,beam.Nm);
% cfg.tests.struct.enable = true;
% cfg.tests.struct.kind   = 'tip-impulse'; % 'free' | 'tip-step'
% % cfg.tests.struct.kind   = 'free'; % 'free' | 'tip-step'
% cfg.tests.struct.dir    = 'z';
% cfg.tests.struct.Nm     = 12;
% cfg.tests.struct.zeta   = 0.01;          % 1% critical per mode (or vector)
% cfg.tests.struct.F0     = 2.0;           % N
% cfg.tests.struct.Timp   = 0.02;          % s
% cfg.tests.struct.tFinal = 1.0;           % s
% cfg.tests.struct.dt     = 1e-3;          % s
% cfg.tests.struct.tip_node = [];          % auto-pick farthest x
% cfg.tests.struct.ic.mode = 1;          

% disp(beam.fem.coordinates)

% ---- structural-only test harness (defaults) ----
if ~isfield(cfg,'tests'); cfg.tests = struct(); end
if ~isfield(cfg.tests,'struct'); cfg.tests.struct = struct(); end
if ~isfield(cfg.tests.struct,'ic'); cfg.tests.struct.ic = struct(); end
if ~isfield(cfg.tests.struct,'force'); cfg.tests.struct.force = struct(); end
if ~isfield(cfg.tests.struct,'trim'); cfg.tests.struct.trim = struct(); end
S = cfg.tests.struct;
S.enable      = getfieldwithdefault(S,'enable',      true);   % set true to run
S.T           = getfieldwithdefault(S,'T',           5.0);     % [s] total time
S.dt          = getfieldwithdefault(S,'dt',          1e-3);    % output step
S.ic.q0       = getfieldwithdefault(S.ic,'q0',       zeros(beam.Nm,1));
S.ic.q1       = getfieldwithdefault(S.ic,'q1',       zeros(beam.Nm,1));
S.ic.q2       = getfieldwithdefault(S.ic,'q2',       zeros(beam.Nm,1));
S.ic.kick_idx = getfieldwithdefault(S.ic,'kick_idx', 1);       % which mode to kick
S.ic.kick_amp = getfieldwithdefault(S.ic,'kick_amp', 0);    % small initial q1 1e-3
S.use_nonlinear= getfieldwithdefault(S,'use_nonlinear', true); % run NL branch too
S.force.type  = getfieldwithdefault(S.force,'type',  'impulse');  % 'none'|'impulse'|'sine'
S.force.A     = getfieldwithdefault(S.force,'A',     1);       % modal force amplitude
S.force.w     = getfieldwithdefault(S.force,'w',     0);       % [rad/s] for 'sine'
S.force.doDamp = cfg.struct.damping;
% S.force.doDamp = false; % testing no damped response

cfg.tests.struct = S;
clear S



% ---- run structural-only modal test  ----
if cfg.tests.struct.enable
    run_struct_only_test(beam, cfg.tests.struct);
end

% 
% if isfield(cfg, 'tests') && isfield(cfg.tests, 'struct') && cfg.tests.struct.enable
%     T = cfg.tests.struct;   % short alias
%     defaults = struct('kind','free', 'dir','z', 'Nm',min(20,size(beam.phi0,2)), ...
%                       'zeta',0.0, 'F0',5.0, 'Timp',0.01, ...
%                       'tFinal',1.0, 'dt',1e-3, 'tip_node',[]);
%     T = merge_cfg(T, defaults);   % simple helper you may already have
%     % disp(size(beam.Sigma));
%     res_struct = run_struct_only_test(beam, beam.fem, beam.keepDofs, struct('test',T));
%     % Optionally stash:
%     sim.results.struct_only = res_struct;
% end

%% Test
% force_lin_cur = aero.DataMatrix.forces_aero_beam_dof;
% force_lin_t =  aero.DataMatrix.Kforces*aero.DataMatrix.forces_aero;
% % disp(size(force_lin_cur));
% % disp(size(force_lin_t));
% % disp(force_lin_t(end-7:end));
% 
% 
% force_lin_tst1 = beam.phi0.'*beam.Mglobal*force_lin_t(1:96);
% force_lin_tst = beam.phi0.'*force_lin_t(1:96);
% force_lin_tst2 = beam.phi1.'*force_lin_t(1:96);
% force_lin_tst3 = beam.phi1.'*beam.Mglobal*force_lin_t(1:96);
% % disp(force_lin_cur)
% % disp(force_lin_tst)
% % disp(force_lin_tst1)
% % disp(force_lin_tst2)
% % disp(force_lin_tst3)
% % err = norm(force_lin_tst-force_lin_cur)
% % err = norm(force_lin_tst-force_lin_cur)
% % err = norm(force_lin_tst-force_lin_cur)
% % err = norm(force_lin_tst-force_lin_cur)
% aero.DataMatrix.forces_aero_beam_dof= force_lin_tst;
% aero.DataMatrix.forces_aero_beam_dof= aero.DataMatrix.forces_aero_beam_dof / cfg.flight.aeroscale;

%% ===== BASE / ξ–COUPLING ======================================
fprintf('[sim_init] Building BASE (ξ-coupling)...\n');

base = AeroFlex.core.build_xi_coupling(beam, aero, cfg);

% Save alongside the others (for MATLAB-side inspection/reuse)
save(fullfile(paths.for_matlab,'base_bundle.mat'), 'base');

% Quick diagnostics
fprintf('[sim_init] BASE OK: φξ %dx%d, Dχ %dx%d, Bχ %dx%d\n', ...
    size(base.phi_xi_modes,1), size(base.phi_xi_modes,2), ...
    size(base.FM.Dchi,1), size(base.FM.Dchi,2), ...
    size(base.FM.Bchi,1), size(base.FM.Bchi,2));
% 
% S = load(fullfile(paths.for_matlab,'base_bundle.mat'),'base'); B0=S.base;
% fprintf('φξ: %dx%d | Dχ: %dx%d | Bχ: %dx%d\n', size(B0.phi_xi_modes,1),size(B0.phi_xi_modes,2), ...
%         size(B0.FM.Dchi,1),size(B0.FM.Dchi,2), size(B0.FM.Bchi,1),size(B0.FM.Bchi,2));
% disp(fieldnames(S.base))







% % Controllability (numeric rank)
% Co = ctrb(A,B);
% rankCo = rank(Co);
% if debug
%     fprintf('[sim_init] ctrb rank = %d (n=%d)\n', rankCo, size(A,1));
% end


%% ------------------------------------------------------------
%  6) Trim (steady states) + indices + runner
% ------------------------------------------------------------

if ~exist('opt','var') || ~isfield(opt,'body_case') || isempty(opt.body_case)
    opt.body_case = 'wingOnly';  % default
end
fprintf('Initial AOA = %g degrees\n',cfg.flight.aoa_deg);
fprintf('Beginning steady-state (trim) for case: %s\n', opt.body_case);
% disp(cfg.flight.aoa_deg) % in degrees

trim = struct(); x0 = []; idx = struct(); simRunner = [];
cfg.sim.storeSens = true;
% trim = struct(); x0 = []; idx = struct(); 
% % base = []; aero = []; beam = []; 
% simRunner = [];
cfg.SwTest =0;
% cfg.sim.storeSens = true;
% optional logfile
% try
%     if exist(paths.logs,'dir')
%         trim_log = fullfile(paths.logs, sprintf('trim_%s.txt', lower(opt.body_case)));
%         diary(trim_log); diary on;
%     end
% catch, end

try
    switch lower(opt.body_case)
        case 'wingonly'
            trim = AeroFlex.sim.trimAircraft(cfg, beam, aero, base);
            x0   = trim.states;
            if ~isfield(trim,'aoa'), trim.aoa = NaN; end
            if ~isfield(trim,'u_ctrl'), trim.u_ctrl = zeros(getOr(cfg.ctrl,'Nc',6),1); end
            fprintf('[trim] wingOnly: |x0|=%d, aoa=%.4g deg\n', numel(x0), trim.aoa);

        case 'coupledfull'
            % Rigid-body + flexible coupled trim
            trim = RigidBody.methods.CoupledTrim.TrimRBwFlex(cfg, beam, aero, base);
            x0   = trim.states;
            if ~isfield(trim,'aoa'), trim.aoa = NaN; end
            if ~isfield(trim,'u_ctrl'), trim.u_ctrl = zeros(getOr(cfg.ctrl,'Nc',6),1); end
            fprintf('[trim] coupledFull: |x0|=%d, aoa=%.4g deg\n', numel(x0), trim.aoa);

        otherwise
            error('Invalid opt.body_case="%s". Use "wingOnly" or "coupledFull".', opt.body_case);
    end
catch ME
    warning('[trim] Trim failed: %s', ME.message);
    trim = struct('states', [], 'aoa', NaN, 'u_ctrl', []);
    x0   = [];  % leave empty; you can still run ROM-only tests below
end

cfg.sim.storeSens = false;

% Build index struct (for downstream controllers/sim)
try
    idx = AeroFlex.core.buildIndexStruct(beam.Nm, aero.Na);
catch ME
    warning('[idx] Failed to build index struct: %s', ME.message);
    idx = struct();
end

% Optional SimRunner
try
    simRunner = AeroFlex.sim.SimRunner(cfg, beam, aero, base, trim);
catch ME
    warning('[SimRunner] Skipped (%s). Continue with ROM export.', ME.message);
end

% % stop logging for trim
% try, diary off; catch, end

%% ------------------------------------------------------------
%  7) Linear ROM handles
% ------------------------------------------------------------
sysd = aero.ROM_dsc;                  % discrete ROM (from SHARPy *.rom.h5)
A = sysd.A; B = sysd.B; C = sysd.C; D = sysd.D; Ts = sysd.Ts;

%% ------------------------------------------------------------
%  8) Timebase + inputs (fix off-by-one)
%     Make the *length of time, gust, deflection* all identical.
% ------------------------------------------------------------
% Number of samples including t=0
nsamp = floor(cfg.sim.t_end / Ts) + 1;      % e.g., t_end=1.0, Ts=0.01 -> 101 samples
t     = (0:nsamp-1)' * Ts;                  % nsamp×1

% Gust column = last input by convention
m = size(B,2);
gust_col = m;

% Fit gust signal to nsamp
ug = aero.gust_input(:);
if isempty(ug), ug = zeros(nsamp,1); end
if numel(ug) < nsamp
    ug = [ug; ug(end) * ones(nsamp - numel(ug), 1)];
elseif numel(ug) > nsamp
    ug = ug(1:nsamp);
end

% Full input matrix U (control + gust)
U = zeros(nsamp, m);
U(:,gust_col) = ug;

% Control-surface command placeholders (zeros for open-loop)
Nc = getOr(cfg.ctrl,'n_surf',2);
deflection       = zeros(nsamp, Nc);
deflection_rate  = zeros(nsamp, Nc);

%% ------------------------------------------------------------
%  9) Optional: quick open-loop linear response (ROM only)
% ------------------------------------------------------------
y = []; tip_deflection = [];
try
    sysd_ol = ss(A,B,C,D,Ts);
    y = lsim(sysd_ol, U, t);
    if ~isempty(y)
        tip_deflection = y(:,1);          % pick first output as a proxy
    end
catch ME
    warning('[openloop] Linear ROM lsim skipped: %s', ME.message);
end

%%

cfg.tests.struct.ic.qxi = zeros(beam.Nm+1,1);
cfg.tests.struct.trim = x0;
if cfg.tests.struct.enable
    run_struct_and_base_test(beam, base, cfg.tests.struct);
end

%% ------------------------------------------------------------
% 10) Export “for_sharpy” artifacts (HIL plant driver)
%     sim_for_sharpy.h5: /time, /gust, /deflection, /deflection_rate
%     meta.yaml: where this run lives, what case, created date
% ------------------------------------------------------------
h5fn = fullfile(paths.for_sharpy,'sim_for_sharpy.h5');
if exist(h5fn,'file'), delete(h5fn); end

% Always create datasets with EXACT sizes you will write
h5create(h5fn, '/time',            [nsamp 1], 'Datatype','double');
h5create(h5fn, '/gust',            [nsamp 1], 'Datatype','double');
h5create(h5fn, '/deflection',      [nsamp Nc], 'Datatype','double');
h5create(h5fn, '/deflection_rate', [nsamp Nc], 'Datatype','double');

h5write(h5fn, '/time',            t);
h5write(h5fn, '/gust',            U(:,gust_col));
h5write(h5fn, '/deflection',      deflection);
h5write(h5fn, '/deflection_rate', deflection_rate);

% % Minimal meta.yaml (no dependency on external YAML writer)
% meta_case     = case_name;
% meta_linss    = fullfile(paths.from_sharpy, [case_name,'.linss.h5']);
% meta_created  = datestr(now,'yyyy-mm-dd');
% meta_bodycase = lower(opt.body_case);
% yaml_text = sprintf([ ...
%     "case_name: %s\n", ...
%     "linss_path: %s\n", ...
%     "created: %s\n", ...
%     "body_case: %s\n", ...
%     "note: Open-loop export. Fill deflection later for closed-loop HIL.\n" ...
%     ], meta_case, meta_linss, meta_created, meta_bodycase);
% fid = fopen(fullfile(paths.run_dir,'meta.yaml'),'w');
% if fid~=-1, fwrite(fid, yaml_text); fclose(fid); end

% yaml_fn = fullfile(paths.run_dir,'meta.yaml');
% fid = fopen(yaml_fn,'w');  assert(fid ~= -1);
% cleanupObj = onCleanup(@() fclose(fid)); 
% fprintf(fid, 'case_name: %s\n', case_name);
% fprintf(fid, 'linss_path: %s\n', fullfile(paths.from_sharpy,[case_name '.linss.h5']));
% fprintf(fid, 'created: %s\n', datestr(now,'yyyy-mm-dd'));
% fprintf(fid, 'body_case: %s\n', lower(opt.body_case));
% fprintf(fid, 'note: Open-loop export. Fill deflection later for closed-loop HIL.\n');
% --- compute body_case string safely + sanitize for consistency with folder name
if exist('opt','var') && isstruct(opt) && isfield(opt,'body_case') && ~isempty(opt.body_case)
    body_case_str = char(opt.body_case);
else
    body_case_str = 'wingOnly';
end
meta_bodycase = sanitize_bodycase(body_case_str);  % 'wing_only' or 'coupled_full'

% --- Write meta.yaml
yaml_fn = fullfile(paths.run_dir,'meta.yaml');
fid = fopen(yaml_fn,'w');  assert(fid ~= -1);
cleanupObj = onCleanup(@() fclose(fid));   % ensures file closes on error

fprintf(fid, 'case_name: %s\n', case_name);
fprintf(fid, 'linss_path: %s\n', fullfile(paths.from_sharpy,[case_name '.linss.h5']));
fprintf(fid, 'created: %s\n', datestr(now,'yyyy-mm-dd'));
fprintf(fid, 'body_case: %s\n', meta_bodycase);
fprintf(fid, 'note: Open-loop export. Fill deflection later for closed-loop HIL.\n');
catch ME
    hadError = true; MEout = ME;
    fprintf(2, '\n[sim_init] ERROR: %s\n', ME.message);
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
end
%% ------------------------------------------------------------
% 11) Save consolidated MATLAB bundle (for your controllers/sims)
% ------------------------------------------------------------
state_space_system = sysd;    % keep same naming as SHARPy MATLAB examples
input_settings     = struct('Ts',Ts,'Nc',Nc,'m',m,'gust_col',gust_col);

sim_config.cfg    = cfg;
sim_config.beam   = beam;
sim_config.aero   = aero;
sim_config.base   = base;
sim_config.trim   = trim;
sim_config.idx    = idx;
sim_config.x0     = x0;
sim_config.sysd   = sysd;
sim_config.paths  = paths;
sim_config.case_name = case_name;
sim_config.body_case = opt.body_case;

save(fullfile(paths.for_matlab,'sim_bundle.mat'), ...
    'sim_config','state_space_system','input_settings','t','U','y','tip_deflection','deflection','deflection_rate');

% S = load(fullfile(paths.for_matlab,'sim_bundle.mat'),'sim_config'); 
% disp(fieldnames(S.sim_config))

%% ===== Runner settings for sim_run =======================================

fprintf('Generating run files for %s runner and assembling L matrix\n', opt.runner);
% Choose defaults once here; you can edit these before calling sim_run.m
run_settings = struct();

% Which overall experiment to run later: 'openloop' or 'nmhe_nmpc'
run_settings.sim_case   = opt.sim_case;   % openloop or 'nmhe_nmpc'

% Which plant integrator family to use later: 'SimRunner' or 'PlantROM'
run_settings.runner     = opt.runner;  % 'SimRunner' or 'PlantROM'

% Controller/estimator model choices (only used by nmhe_nmpc + PlantROM)
run_settings.modelFcn       = 'AeroFlex.sim.ROMIntegrator';
run_settings.sensorFcn      = 'AeroFlex.sensor.WingVelSensor';
run_settings.estimatorFcn   = 'AeroFlex.ctrl.nMHEv2';
run_settings.controllerFcn  = 'AeroFlex.ctrl.nMPC';

% Control sampling; if you already set in cfg, mirror here for robustness
if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'Ts') && ~isempty(cfg.ctrl.Ts)
    run_settings.ctrl_Ts = cfg.ctrl.Ts;
else
    run_settings.ctrl_Ts = cfg.sim.dt;  % fallback
end
GLAConfig;
% Number of control/estimation horizon steps if applicable
run_settings.Nc = getOr(cfg.ctrl,'Nc', 6);
run_settings.Ne = getOr(cfg.ctrl,'Ne', 6);

% Persist trim initial state for convenience
run_settings.x0 = sim_config.x0;

% Save L-matrix if available (guarded)
try
    [L, blk] = AeroFlex.core.assemble_L_matrix(cfg, beam, aero, base);
    save(fullfile(paths.for_matlab,'L_bundle.mat'),'L','blk');
catch ME
    warning('[sim_init] assemble_L_matrix skipped: %s', ME.message);
end

% Save runner settings alongside other for_matlab artifacts
save(fullfile(paths.for_matlab,'run_settings.mat'),'run_settings');


%% ------------------------------------------------------------
% 12) Final OK
% ------------------------------------------------------------
%% ===== Console log capture (finalize) ====================================
status_str = ternary(~hadError, 'OK', 'ERROR');
fprintf('=== sim_init finished at %s (status: %s) ===\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), status_str);

% Duplicate the log in run root for convenience:
try
    copyfile(console_log, fullfile(paths.run_dir, 'console.log'));
catch MEc
    fprintf(2, '[log] Could not copy console log to run root: %s\n', MEc.message);
end

% If you want a copy next to the MATLAB bundles too:
try
    copyfile(console_log, fullfile(paths.for_matlab, 'console.log'));
catch, end

% (diary will be closed by onCleanup, but we close explicitly as well)
diary off;

% Set success and rethrow if you want failures to stop CI:
if hadError
    success = false;
    % rethrow(MEout); % uncomment if you want the error to bubble up
else
    success = true;
end


fprintf('[sim_init] Completed. Run folder: %s\n', paths.run_dir);

% input(here)

end
% ------------- small helpers -------------
function val = getOr(S, f, dflt)
    if nargin<3, dflt = []; end
    if isstruct(S)
        if isfield(S,f) && ~isempty(S.(f)), val = S.(f); else, val = dflt; end
    else
        val = dflt;
    end
end

%% ================================
%  Local helpers (single-file)
%  ===============================

function s = strip_trailing_filesep(s)
    s = char(s);
    if ~isempty(s) && s(end) == filesep, s = s(1:end-1); end
end

function [run_dir, paths] = create_run_dirs(sim_root, case_name, body_case, overwrite, date_only_runs)
    if ~exist(sim_root,'dir'), mkdir(sim_root); end
    case_dir = fullfile(sim_root, case_name);
    if ~exist(case_dir,'dir'), mkdir(case_dir); end

    % Normalize folder name for body_case
    body_case = sanitize_bodycase(body_case);
    body_dir  = fullfile(case_dir, body_case);
    if ~exist(body_dir,'dir'), mkdir(body_dir); end

    % --- choose stamp style ---
    % if date_only_runs
        stamp  = datestr(now,'yyyymmdd');           % e.g. 20251110
        run_dir = fullfile(body_dir, stamp);

        if exist(run_dir,'dir')
            if overwrite
                warning('[create_run_dirs] Reusing existing %s (overwrite=true). Files may be overwritten.', run_dir);
            else
                fprintf('[create_run_dirs] Using existing daily folder: %s\n', run_dir);
            end
        else
            mkdir(run_dir);
        end
    % else
    %     stamp  = datestr(now,'yyyymmdd_HHMMSS');    % e.g. 20251110_161245
    %     run_dir = fullfile(body_dir, stamp);
    %     if exist(run_dir,'dir')
    %         if ~overwrite
    %             error('Run folder already exists: %s (pass ''overwrite'',true or enable date_only_runs)', run_dir);
    %         end
    %     else
    %         mkdir(run_dir);
    %     end
    % end

    % Subfolders (under the new run_dir)
    paths.run_dir     = run_dir;
    paths.from_sharpy = must_mkdir(fullfile(run_dir,'from_sharpy'));
    paths.for_matlab  = must_mkdir(fullfile(run_dir,'for_matlab'));
    paths.for_sharpy  = must_mkdir(fullfile(run_dir,'for_sharpy'));
    paths.plots       = must_mkdir(fullfile(run_dir,'plots'));
    paths.logs        = must_mkdir(fullfile(run_dir,'logs'));
end

function d = must_mkdir(d)
    if ~exist(d,'dir'), mkdir(d); end
end

function s = sanitize_bodycase(s)
    s = lower(char(s));
    s = strrep(s,' ','_');
    % Map your canonical names to nicer folder names
    if strcmp(s,'wingonly'),    s = 'wing_only';    end
    if strcmp(s,'coupledfull'), s = 'coupled_full'; end
end

% function [run_dir, paths] = create_run_dirs(sim_root, case_name, overwrite)
%     if ~exist(sim_root,'dir'), mkdir(sim_root); end
%     case_dir = fullfile(sim_root, case_name);
%     if ~exist(case_dir,'dir'), mkdir(case_dir); end
%     % stamp = datestr(now,'yyyymmdd_HHMMSS');
%     stamp = datestr(now,'yyyymmdd');
% 
%     run_dir = fullfile(case_dir, stamp);
%     if exist(run_dir,'dir')
%         if ~overwrite
%             error('Run folder already exists: %s (use ''overwrite'',true)', run_dir);
%         end
%     else
%         mkdir(run_dir);
%     end
%     paths.run_dir     = run_dir;
%     paths.from_sharpy = must_mkdir(fullfile(run_dir,'from_sharpy'));
%     paths.for_matlab  = must_mkdir(fullfile(run_dir,'for_matlab'));
%     paths.for_sharpy  = must_mkdir(fullfile(run_dir,'for_sharpy'));
%     paths.plots       = must_mkdir(fullfile(run_dir,'plots'));
%     paths.logs        = must_mkdir(fullfile(run_dir,'logs'));
% end
% 
% function d = must_mkdir(d)
%     if ~exist(d,'dir'), mkdir(d); end
% end


function cfg_out = merge_cfg(defaults, loaded)
    cfg_out = defaults;
    f = fieldnames(loaded);
    for k = 1:numel(f)
        key = f{k};
        if isstruct(getfield(loaded,key)) && isfield(cfg_out,key) && isstruct(cfg_out.(key))
            cfg_out.(key) = merge_cfg(cfg_out.(key), loaded.(key));
        else
            cfg_out.(key) = loaded.(key);
        end
    end
end

function validate_cfg(cfg)
    % Define required high-level fields and nested ones
    req_top = {'Nm','Na','sim','flight','ctrl','gust','trim','struct'};
    require_fields(cfg, req_top, 'cfg');

    require_fields(cfg.sim,    {'dt','t_end'},           'cfg.sim');
    require_fields(cfg.flight, {'U_inf','rho','b_ref'},  'cfg.flight');
    require_fields(cfg.ctrl,   {'n_surf'},               'cfg.ctrl');

    % Derived sanity checks
    assert(cfg.Nm>0 && cfg.Na>=0, 'cfg.Nm and cfg.Na must be positive (Na>=0)');
    assert(cfg.sim.dt>0 && cfg.sim.t_end>0, 'cfg.sim.{dt,t_end} must be > 0');
    assert(cfg.flight.U_inf>0, 'cfg.flight.U_inf must be > 0');
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function require_fields(S, names, rootname)
    for i=1:numel(names)
        n = names{i};
        if ~isfield(S, n)
            error('Missing required field "%s.%s"', rootname, n);
        end
    end
end

% function v = getOr(S, name, default)
%     if isfield(S, name), v = S.(name); else, v = default; end
% end


function cfg = ensure_nominal_extras(cfg)
%ENSURE_NOMINAL_EXTRAS
% Fill *missing* fields with sensible nominal defaults and fix dimensions.
% Never overwrites fields that already exist in cfg.

    % ---------- STRUCT ----------
    if ~isfield(cfg,'struct') || ~isstruct(cfg.struct)
        cfg.struct = struct();
    end
    if ~isfield(cfg.struct,'L'), cfg.struct.L = 1.1; end  % [m] span for S_ref

    % ---------- FLIGHT ----------
    if ~isfield(cfg,'flight') || ~isstruct(cfg.flight)
        cfg.flight = struct();
    end
    must_have(cfg.flight,'U_inf','cfg.flight.U_inf (free-stream speed) missing');
    must_have(cfg.flight,'rho',  'cfg.flight.rho (air density) missing');
    must_have(cfg.flight,'b_ref','cfg.flight.b_ref (semi-chord) missing');

    if ~isfield(cfg.flight,'aoa_deg'), cfg.flight.aoa_deg = 1; end
    if ~isfield(cfg.flight,'S_ref')
        cfg.flight.S_ref = cfg.struct.L * cfg.flight.b_ref; % rectangular planform
    end

    % Derived non-dimensionalization knobs (keep nominal values unless user/SHARPy gave them)
    if ~isfield(cfg.flight,'a'),      cfg.flight.a = 1; end
    if ~isfield(cfg.flight,'t_inf'),  cfg.flight.t_inf = 1; end

    % Also keep a handy dynamic pressure for other derived fields (not stored if you don't want)
    q_inf = 0.5 * cfg.flight.rho * cfg.flight.U_inf^2; 
    % disp(cfg.flight.U_inf)
    % disp(q_inf);
    % if ~isfield(cfg.flight.aeroscale), cfg.flight.aeroscale = q_inf*cfg.flight.b_ref^2; end
    % if ~isfield(cfg.flight.t_scale), cfg.flight.t_scale = cfg.flight.b_ref/cfg.flight.U_inf; end
    % if ~isfield(cfg.flight.Fscale), cfg.flight.Fscale = sqrt(cfg.ctrl.Ts/(cfg.flight.aeroscale*cfg.flight.t_scale)); end
    
    
    % disp(cfg.flight.Fscale)
    % ---------- TRIM ----------
    if ~isfield(cfg,'trim') || ~isstruct(cfg.trim)
        cfg.trim = struct();
    end
    if ~isfield(cfg.trim,'do'),         cfg.trim.do        = true;     end
    if ~isfield(cfg.trim,'tolSteady'),  cfg.trim.tolSteady = 1e-6;     end
    if ~isfield(cfg.trim,'tolNewton'),  cfg.trim.tolNewton = 1e-4;     end
    if ~isfield(cfg.trim,'alphaDeg'),   cfg.trim.alphaDeg  = 0;        end
    if ~isfield(cfg.trim,'deltaDeg'),   cfg.trim.deltaDeg  = 0;        end
    if ~isfield(cfg.trim,'itersInt'),   cfg.trim.itersInt  = 10000;    end
    if ~isfield(cfg.trim,'dtQS'),       cfg.trim.dtQS      = 1e-4;     end

    % Ensure thrust length matches Nm
    if ~isfield(cfg,'Nm'), error('cfg.Nm missing before ensuring trim.thrust'); end
    if ~isfield(cfg.trim,'thrust') || ~isvector(cfg.trim.thrust) || numel(cfg.trim.thrust)~=cfg.Nm
        cfg.trim.thrust = zeros(cfg.Nm,1);
    else
        % force column vector shape
        cfg.trim.thrust = cfg.trim.thrust(:);
    end

    % ---------- AERO ----------
    if ~isfield(cfg,'aero') || ~isstruct(cfg.aero)
        cfg.aero = struct();
    end
    if ~isfield(cfg.aero,'discrete'), cfg.aero.discrete = false; end

    % ---------- GUST ----------
    if ~isfield(cfg,'gust') || ~isstruct(cfg.gust)
        cfg.gust = struct();
    end
    if ~isfield(cfg.gust,'on'),            cfg.gust.on            = true;   end
    if ~isfield(cfg.gust,'gust_length'),   cfg.gust.gust_length   = 10.0;   end
    if ~isfield(cfg.gust,'gust_intensity'),cfg.gust.gust_intensity= 0.02;   end
    if ~isfield(cfg.gust,'gust_offset'),   cfg.gust.gust_offset   = 0;      end

    % Only set these if absent, using dimensional forms from your nominal config
    if ~isfield(cfg.gust,'a')
        cfg.gust.a = (0.5*cfg.flight.rho*cfg.flight.U_inf^2) * cfg.flight.b_ref^2;
    end
    if ~isfield(cfg.gust,'t_inf')
        cfg.gust.t_inf = cfg.flight.b_ref / cfg.flight.U_inf;
    end

    % ---------- RIGID (mirror gust scales by default) ----------
    if ~isfield(cfg,'rigid') || ~isstruct(cfg.rigid)
        cfg.rigid = struct();
    end
    if ~isfield(cfg.rigid,'a'),     cfg.rigid.a     = cfg.gust.a;     end
    if ~isfield(cfg.rigid,'t_inf'), cfg.rigid.t_inf = cfg.gust.t_inf; end

    % ---------- CONTROL-SURFACE META ----------
    if ~isfield(cfg,'ctrl') || ~isstruct(cfg.ctrl)
        cfg.ctrl = struct();
    end
    if ~isfield(cfg.ctrl,'n_surf'),   cfg.ctrl.n_surf = 2; end % validate_cfg already checks n_surf presence
    if ~isfield(cfg.ctrl,'var_per'),  cfg.ctrl.var_per = 2; end
    if ~isfield(cfg.ctrl,'Nc'),       cfg.ctrl.Nc      = 6; end
    if ~isfield(cfg.ctrl,'Ne'),       cfg.ctrl.Ne      = 6; end
end

function must_have(S,field,msg)
    if ~isfield(S,field)
        error('%s', msg);
    end
end

function show_cfg_summary(cfg)
    fprintf('---- cfg summary ----\n');
    fprintf('Nm=%d, Na=%d\n', cfg.Nm, cfg.Na);
    fprintf('dt=%.6g s, t_end=%.3f s\n', cfg.sim.dt, cfg.sim.t_end);
    fprintf('U_inf=%.3f m/s, rho=%.3f kg/m^3, b_ref=%.4f m\n', ...
        cfg.flight.U_inf, cfg.flight.rho, cfg.flight.b_ref);
    fprintf('gust: L=%.3f sec, I=%.4f\n', getOr(cfg.gust,'gust_length',NaN), getOr(cfg.gust,'gust_intensity',NaN));
    fprintf('ctrl.n_surf=%d\n', cfg.ctrl.n_surf);
    fprintf('---------------------\n');
end

function [A,B,C,D,Ts,eta_ref] = read_linss_h5_fallback(h5file)
    % Minimal HDF5 fallback reader that matches SHARPy linss export
    A = h5read(h5file,'/A'); B = h5read(h5file,'/B');
    C = h5read(h5file,'/C'); D = h5read(h5file,'/D');
    Ts = h5read(h5file,'/dt'); 
    if numel(Ts)==1, Ts = double(Ts); else, Ts = double(Ts(1)); end
    eta_ref = [];
    try
        eta_ref = h5read(h5file,'/eta_ref');
    catch
        % ok if absent
    end
end

function [ok,msg] = sanity_check_linss(A,B,C,D,Ts)
    ok = true; msg = '';
    n = size(A,1);
    if ~ismatrix(A) || size(A,2)~=n, ok=false; msg='A must be square'; return; end
    if size(B,1)~=n || size(C,2)~=n, ok=false; msg='B/C dimension mismatch'; return; end
    if size(D,1)~=size(C,1) || size(D,2)~=size(B,2), ok=false; msg='D size mismatch'; return; end
    if ~isscalar(Ts) || ~(Ts>0), ok=false; msg='Ts must be positive scalar'; return; end
end

function input_settings = set_input_parameters(u_inf, Ts, num_modes, num_aero_states, ...
    rigid_body_motions, simulation_time, num_control_surfaces, n_nodes, control_input_start, gust_input_start)
    % Keep the same field names SHARPy helpers expect downstream
    input_settings = struct();
    input_settings.u_inf = u_inf;
    input_settings.Ts = Ts;
    input_settings.num_modes = num_modes;
    input_settings.num_aero_states = num_aero_states;
    input_settings.rigid_body_motions = logical(rigid_body_motions);
    input_settings.simulation_time = simulation_time;
    input_settings.num_control_surfaces = num_control_surfaces;
    input_settings.n_nodes = n_nodes;
    input_settings.control_input_start = control_input_start;
    input_settings.gust_input_start = gust_input_start;
end

function gust = get_1minuscosine_gust_input(gust_length, gust_intensity, Ts, u_inf, sim_time)
    % Time-domain 1-cos discrete gust
    t = (0:round(sim_time/Ts))' * Ts;
    s = t * u_inf; % distance along flight path
    gust = struct();
    gust.time = t;
    gust.wgust = zeros(size(t));
    L = gust_length;
    Uds = gust_intensity * u_inf; % scale by U_inf
    for k=1:numel(t)
        if s(k)<=L
            gust.wgust(k) = Uds/2*(1 - cos(pi*s(k)/L));
        else
            gust.wgust(k) = 0;
        end
    end
end

function ug = local_one_minus_cos(t, L, Ts, u_inf, sim_time, intensity)
    s = t * u_inf;
    Uds = intensity * u_inf;
    ug = zeros(size(t));
    in = s<=L;
    ug(in) = Uds/2*(1 - cos(pi*s(in)/L));
end

function cfg2 = map_sharpy_sim_params_to_cfg(S)
    % Map SHARPy's simulation_parameters.mat to our cfg fields, when present
    cfg2 = struct();
    safe_copy = @(dst,src,name) setfield(dst, name, src.(name));
    try
        if isfield(S,'u_inf'),   cfg2.flight.U_inf = S.u_inf; end
        if isfield(S,'rho'),     cfg2.flight.rho   = S.rho;   end
        if isfield(S,'b_ref'),   cfg2.flight.b_ref = S.b_ref; end
        if isfield(S,'simulation_time'), cfg2.sim.t_end = S.simulation_time; end
        if isfield(S,'num_modes'), cfg2.Nm = double(S.num_modes); end
        if isfield(S,'num_aero_states'), cfg2.Na = double(S.num_aero_states); end
        if isfield(S,'gust_length'),    cfg2.gust.gust_length = S.gust_length; end
        if isfield(S,'gust_intensity'), cfg2.gust.gust_intensity = S.gust_intensity; end
        if isfield(S,'num_control_surfaces'), cfg2.ctrl.n_surf = double(S.num_control_surfaces); end
        if isfield(S,'dt'), cfg2.sim.dt = S.dt; end
        if isfield(S,'n_nodes'), cfg2.n_nodes = uint8(S.n_nodes); end
    catch
        % Best-effort mapping
    end
end

function writelines_struct_yaml(meta, fn)
    % Tiny YAML-like dumper (safe for small meta)
    f = fieldnames(meta);
    lines = strings(numel(f),1);
    for i=1:numel(f)
        key = f{i}; val = meta.(key);
        if ischar(val) || isstring(val)
            v = string(val);
        elseif isnumeric(val) && isscalar(val)
            v = num2str(val);
        else
            v = jsonencode(val);
        end
        lines(i) = sprintf('%s: %s', key, v);
    end
    fid = fopen(fn,'w'); 
    for i=1:numel(lines), fprintf(fid,'%s\n', lines(i)); end
    fclose(fid);
end



function run_struct_only_test(beam, S)
%STRUCTURAL-ONLY MODAL TEST (Artola-style)
% Runs linear and nonlinear ROMs, reconstructs nodal outputs, and plots.

Nm   = beam.Nm;
Om   = beam.Omega;                 % diag(omega_i)
% Sig  = beam.Sigma(:);              % Nm x 1 damping (Kelvin-Voigt-projected) – already built in beam
Sig  = beam.Sigma;              % Nm x Nm damping (Kelvin-Voigt-projected) – already built in beam

if ~S.force.doDamp
    Sig = zeros(Nm,Nm);
end
% disp(S.use_nonlinear);
% disp(isfield(beam,'Gamma1'));
% disp(isfield(beam,'Gamma2'));
% disp(~isempty(beam.Gamma1));
% disp(~isempty(beam.Gamma2));
% haveNL = S.use_nonlinear && isfield(beam,'Gamma1') && isfield(beam,'Gamma2') ...
%          && ~isempty(beam.Gamma1) && ~isempty(beam.Gamma2)
haveNL = S.use_nonlinear && ~isempty(beam.Gamma1) && ~isempty(beam.Gamma2);
% -------------------------------------------------------------------------
% Initial conditions (small modal kick in q1 by default)
q0_0 = S.ic.q0(:); 
q1_0 = S.ic.q1(:); 
q2_0 = S.ic.q2(:);
if ~any(q1_0) && S.ic.kick_amp~=0
    k = max(1, min(Nm, S.ic.kick_idx));
    q1_0(k) = S.ic.kick_amp;
end
%% trim val test
% q1_0 = S.trim(1:Nm);
% q2_0 = S.trim(Nm+1:2*Nm);
%%
x0_lin = [q1_0; q2_0];     % state = [q1;q2]
x0_nl  = x0_lin;

disp('starting Structure Test')
% -------------------------------------------------------------------------
% External modal forcing η_e(t) (acts in the q1-equation)
eta_fun = make_eta_fun(S.force, Nm);

% -------------------------------------------------------------------------
% Time span
T  = S.T;
tq = 0:S.dt:T;
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

% ===== Linear system (Artola 3.33, without Γ terms) ======================
f_lin = @(t,x) rhs_linear(t,x,Om,Sig,eta_fun);
[tl, xl] = ode45(f_lin, [tq(1) tq(end)], x0_lin, opts);
xl = interp1(tl, xl, tq, 'linear');  % shape Nt x (2Nm)
disp('Linear Structure Test Done')

% ===== Nonlinear system ===================================================
if haveNL
    f_nl = @(t,x) rhs_full(t,x,Om,Sig,beam.Gamma1,beam.Gamma2,eta_fun);
    [tn, xn] = ode45(f_nl, [tq(1) tq(end)], x0_nl, opts);
    xn = interp1(tn, xn, tq, 'linear');
else
    xn = [];
end

% -------------------------------------------------------------------------
% Recover q0 by integrating  q̇0 = q1  (component-wise)
q1_lin = xl(:, 1:Nm).';                      % Nm x Nt
q2_lin = xl(:, Nm+1:end).';
% q0_lin = cumtrapz(tq, q1_lin, 2);            % Nm x Nt
q0_lin = -Om\q2_lin;            % Nm x Nt Using Artola Def

if ~isempty(xn)
    q1_nl = xn(:, 1:Nm).';
    q2_nl = xn(:, Nm+1:end).';
    % q0_nl = cumtrapz(tq, q1_nl, 2);
    q0_nl = -Om\q2_nl;            % Nm x Nt Using Artola Def
end

% -------------------------------------------------------------------------
% Nodal post-processing (tip selection + z deflection using phi0)
tip = pick_tip_node(beam.fem);               % robust: max(Y)
row_phi0_tipz = dof_row_in_phi0(beam.keepDofs, tip, 3); % (u,v,w,rx,ry,rz) => w=3
if row_phi0_tipz==0
    warning('Tip z row was eliminated; picking nearest available flexible z dof.');
    row_phi0_tipz = nearest_flexible_z(beam.keepDofs); % fallback
end
%% Test eigvect
% disp(row_phi0_tipz);
% phi0_tst = load("eigenvectors.dat");
% disp(size(phi0_tst));
% phi0z = phi0_tst(row_phi0_tipz, :);         % 1 x Nm

% dlmwrite('phi0Mat', beam.phi0,'');
% writematrix(beam.phi0,'phi0Mat.dat','Delimiter','tab');

phi0z = beam.phi0(row_phi0_tipz, :);         % 1 x Nm
% disp(size(phi0z));
z_tip_lin = (phi0z * q0_lin).';              % Nt x 1
if ~isempty(xn)
    z_tip_nl  = (phi0z * q0_nl).';
end

% -------------------------------------------------------------------------
% Plots
figure('Name','Struct-only: tip deflection (linear vs nonlinear)','Color','w');
plot(tq, z_tip_lin*100, 'LineWidth', 1.8); hold on;
if ~isempty(xn), plot(tq, z_tip_nl*100, 'LineWidth', 1.8); end
grid on; xlabel('t [s]'); ylabel('z_{tip} [% of m]');   % (raw meters *100 = “percent of m”)
ttl = 'Structural-only: tip z(t) from \phi_0 (integrated q_0)';
if ~isempty(xn), legend('Linear','Full (nonlinear)'); else, legend('Linear'); end
title(ttl);

% q0, q1, q2 panels (show first up-to-4 modes)
idx_show = 1:min(4,Nm);
figure('Name','Modal coordinates (q0,q1,q2)','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on;  title('q_0(t)'); grid on; xlabel('t [s]');
for j = idx_show, plot(tq, q0_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q0_nl(j,:), '--', 'LineWidth', 1.0); end
end
% legend(make_leg(idx_show,~isempty(xn)));
legend(make_leg(idx_show,isempty(xn)));

nexttile; hold on; title('q_1(t)'); grid on; xlabel('t [s]');
for j = idx_show, plot(tq, q1_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q1_nl(j,:), '--', 'LineWidth', 1.0); end
end

nexttile; hold on; title('q_2(t)'); grid on; xlabel('t [s]');
for j = idx_show, plot(tq, q2_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q2_nl(j,:), '--', 'LineWidth', 1.0); end
end

% (Optional) also expose strain/velocity fields via phi1/phi2 if you want:
% x1_lin(s,t) = beam.phi1 * q1_lin(:,t);  ( 6*Ne x 1 per t ) – similar for x2.

end

% ================== helpers ==================

function eta_fun = make_eta_fun(force, Nm)
% modal forcing η_e(t) in the q1-equation
switch lower(force.type)
    case 'none'
        eta_fun = @(t) zeros(Nm,1);
    case 'impulse'
        % narrow pulse at t≈0 (numerically integrable)
        A = force.A;
        eta_fun = @(t) (t<1e-3) * (A/1e-3) * ones(Nm,1);
    case 'sine'
        A = force.A; w = force.w;
        eta_fun = @(t) A*sin(w*t) * ones(Nm,1);
    otherwise
        eta_fun = @(t) zeros(Nm,1);
end
end

function dx = rhs_linear(t, x, Om, Sig, eta_fun)
% function dx = rhs_linear(~, x, Om, Sig, eta_fun)
Nm = size(Om,1);
q1 = x(1:Nm);
q2 = x(Nm+1:end);

% dx1 = Om*q2 + Sig*q1 + eta_fun(0);     % Σ is used as diagonal (modal) damping
dx1 = Om*q2 + Sig*q1 + eta_fun(t);     % Σ is used as diagonal (modal) damping
dx2 = -Om*q1;
dx  = [dx1; dx2];
end

function dx = rhs_full(t, x, Om, Sig, G1, G2, eta_fun)
% Full structural ROM: linear part + quadratic contractions
Nm = size(Om,1);
q1 = x(1:Nm);
q2 = x(Nm+1:end);

% Quadratic/bilinear terms:
% g11 = contractGamma(G1, q1, q1);  % i: sum_{jk} G1(i,j,k) q1_j q1_k
% g22 = contractGamma(G2, q2, q2);  % i: sum_{jk} G2(i,j,k) q2_j q2_k
% g21 = contractGamma(G2, q2, q1);  % i: sum_{jk} G2(i,j,k) q2_j q1_k  (bilinear)

%%
% ---------- Γ₁ & Γ₂  contractions  (vectorised) ---------------------------
%    −Σ_l  q1_l Γ₁(:,:,l) q1    = - Σ_l  diag(q1) *
%                                  vec(Gamma1)^T * kron(q1,I)
% Use pagemtimes (R2023b) :
G1q1 = squeeze(pagemtimes(G1, 'none', q1, 'none'));  % Nm×Nm×1  <- Σ over 3rd dim
N1_G1 = -(G1q1(:,:,1) * q1);
% N(idx.q1) = N(idx.q1) + N1_G1;

G2q2 = squeeze(pagemtimes(G2, 'none', q2, 'none'));
N1_G2 = -(G2q2(:,:,1) * q2);
% N(idx.q2) = N(idx.q2) + N1_G2;

% % ---------- Γ_g  term -----------------------------------------------------
% Ggqx = squeeze(pagemtimes(par.Gamma_g, 'none', q_xi, 'none'));
% N1_Gg =  (Ggqx(:,:,1) * q_xi);
% % N(idx.q2) = N(idx.q2) + N1_Gg;
% 
% % ---------- Γ_xi  term  (Nm+1×1) -----------------------------------------
% Gxiq1 = squeeze(pagemtimes(par.Gamma_xi, 'none', q1, 'none'));
% N_qxi =  (Gxiq1(:,:,1) * q_xi);

% ---------- Γ2^T q1  (Eq 9b) ----------------------------------------------
GT = squeeze(pagemtimes(G2, 'transpose', q1, 'none'));
N_GT =  GT(:,:,1)*q2;     % note: transpose on first arg
%%

% dx1 = Om*q2 + Sig*q1 - g11 - g22 + eta_fun(t);
dx1 = Om*q2 + Sig*q1 + N1_G1 + N1_G2 + eta_fun(t);
% dx2 = -Om*q1 + g21;
dx2 = -Om*q1 + N_GT;
dx  = [dx1; dx2];
end

function run_struct_and_base_test(beam, base, S)
%STRUCTURAL-ONLY MODAL TEST (Artola-style)
% Runs linear and nonlinear ROMs, reconstructs nodal outputs, and plots.

Nm   = beam.Nm;
Om   = beam.Omega;                 % diag(omega_i)
% Sig  = beam.Sigma(:);              % Nm x 1 damping (Kelvin-Voigt-projected) – already built in beam
Sig  = beam.Sigma;              % Nm x Nm damping (Kelvin-Voigt-projected) – already built in beam

if ~S.force.doDamp
    Sig = zeros(Nm,Nm);
end
% disp(S.use_nonlinear);
% disp(isfield(beam,'Gamma1'));
% disp(isfield(beam,'Gamma2'));
% disp(~isempty(beam.Gamma1));
% disp(~isempty(beam.Gamma2));
% haveNL = S.use_nonlinear && isfield(beam,'Gamma1') && isfield(beam,'Gamma2') ...
%          && ~isempty(beam.Gamma1) && ~isempty(beam.Gamma2)
haveNL = S.use_nonlinear && ~isempty(beam.Gamma1) && ~isempty(beam.Gamma2);
% -------------------------------------------------------------------------
% Initial conditions (small modal kick in q1 by default)
q0_0 = S.ic.q0(:); 
q1_0 = S.ic.q1(:); 
q2_0 = S.ic.q2(:);
qxi_0 = S.ic.qxi(:);
if ~any(q1_0) && S.ic.kick_amp~=0
    k = max(1, min(Nm, S.ic.kick_idx));
    q1_0(k) = S.ic.kick_amp;
end
%% trim val test
% q1_0 = S.trim(1:Nm);
% q2_0 = S.trim(Nm+1:2*Nm);
qxi_0 = S.trim(2*Nm+1:3*Nm+1);

quat_ini = base.phi_xi_modes*qxi_0;
qxi_oned = zeros(size(qxi_0)); qxi_oned(1) = 1;
quat_ini_oned = base.phi_xi_modes*qxi_oned;
qxi_0  = qxi_oned;
% qxi_0 =qxi_oned;
% disp(norm(qxi_0))
%%
x0_lin = [q1_0; q2_0; qxi_0];     % state = [q1;q2]
x0_nl  = x0_lin;

disp('starting Structure-Base Test')
% -------------------------------------------------------------------------
% External modal forcing η_e(t) (acts in the q1-equation)
eta_fun = make_eta_fun(S.force, Nm);

% -------------------------------------------------------------------------
% Time span
T  = S.T;
tq = 0:S.dt:T;
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

% % ===== Linear system (Artola 3.33, without Γ terms) ======================
% f_lin = @(t,x) rhs_linear(t,x,Om,Sig,eta_fun);
% [tl, xl] = ode45(f_lin, [tq(1) tq(end)], x0_lin, opts);
% xl = interp1(tl, xl, tq, 'linear');  % shape Nt x (2Nm)
% disp('Linear Structure Test Done')

% ===== Nonlinear system ===================================================
if haveNL
    f_nl = @(t,x) rhs_full_base(t,x,Om,Sig,beam.Gamma1,beam.Gamma2, base.Gamma_g, base.Gamma_xi, eta_fun);
    [tn, xn] = ode45(f_nl, [tq(1) tq(end)], x0_nl, opts);
    xn = interp1(tn, xn, tq, 'linear');
else
    xn = [];
end

% -------------------------------------------------------------------------
% % Recover q0 by integrating  q̇0 = q1  (component-wise)
% q1_lin = xl(:, 1:Nm).';                      % Nm x Nt
% q2_lin = xl(:, Nm+1:end).';
% % q0_lin = cumtrapz(tq, q1_lin, 2);            % Nm x Nt
% q0_lin = -Om\q2_lin;            % Nm x Nt Using Artola Def

% if ~isempty(xn)
    q1_nl = xn(:, 1:Nm).';
    q2_nl = xn(:, Nm+1:2*Nm).';
    % q0_nl = cumtrapz(tq, q1_nl, 2);
    q0_nl = -Om\q2_nl;            % Nm x Nt Using Artola Def
    qxi_nl = xn(:, 2*Nm+1:3*Nm+1).';
% end

% -------------------------------------------------------------------------
% Nodal post-processing (tip selection + z deflection using phi0)
tip = pick_tip_node(beam.fem);               % robust: max(Y)
row_phi0_tipz = dof_row_in_phi0(beam.keepDofs, tip, 3); % (u,v,w,rx,ry,rz) => w=3 for tip
if row_phi0_tipz==0
    warning('Tip z row was eliminated; picking nearest available flexible z dof.');
    row_phi0_tipz = nearest_flexible_z(beam.keepDofs); % fallback
end
%% Test eigvect
% disp(row_phi0_tipz);
% phi0_tst = load("eigenvectors.dat");
% disp(size(phi0_tst));
% phi0z = phi0_tst(row_phi0_tipz, :);         % 1 x Nm

% dlmwrite('phi0Mat', beam.phi0,'');
% writematrix(beam.phi0,'phi0Mat.dat','Delimiter','tab');

phi0z = beam.phi0(row_phi0_tipz, :);         % 1 x Nm
% disp(size(phi0z));
% z_tip_lin = (phi0z * q0_lin).';              % Nt x 1
if ~isempty(xn)
    z_tip_nl  = (phi0z * q0_nl).';
end
% disp(size(base.phi_xi_modes));
% disp(size(qxi_nl));
node_quat = tip;
% node_quat = tip-1;
quat_dof = double((node_quat - 1)*4) + (1:4);
quat_modes = base.phi_xi_modes(quat_dof,:);
% quat_full = (base.phi_xi_modes * qxi_nl);
quat_full = (quat_modes * qxi_nl);
eul_full = zeros(3,size(quat_full,2));
for i = 1:size(quat_full,2)
    eul_full(:,i) = quat2eul(quat_full(:,i).').';
end
eul_full = eul_full.';
% disp(size(quat_full));
% disp(quat_full);

% -------------------------------------------------------------------------
% Plots
figure('Name','Struct-only: tip deflection (linear vs nonlinear)','Color','w');
% plot(tq, z_tip_lin*100, 'LineWidth', 1.8); hold on;
if ~isempty(xn), plot(tq, z_tip_nl*100, 'LineWidth', 1.8); end
grid on; xlabel('t [s]'); ylabel('z_{tip} [% of m]');   % (raw meters *100 = “percent of m”)
ttl = 'Structural-only: tip z(t) from \phi_0 (integrated q_0)';
if ~isempty(xn), legend('Linear','Full (nonlinear)'); else, legend('Linear'); end
title(ttl);

% q0, q1, q2 panels (show first up-to-4 modes)


% disp(size(q0_nl));
% disp(q0_nl);

idx_show = 1:min(4,Nm);
figure('Name','Modal coordinates (q0,q1,q2)','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
nexttile; hold on;  title('q_0(t)'); grid on; xlabel('t [s]');
% for j = idx_show, plot(tq, q0_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q0_nl(j,:), 'LineWidth', 1.0); end
end
legend(make_leg(idx_show,~isempty(xn)));

nexttile; hold on; title('q_1(t)'); grid on; xlabel('t [s]');
% for j = idx_show, plot(tq, q1_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q1_nl(j,:), '--', 'LineWidth', 1.0); end
end

nexttile; hold on; title('q_2(t)'); grid on; xlabel('t [s]');
% for j = idx_show, plot(tq, q2_lin(j,:), 'LineWidth', 1.2); end
if ~isempty(xn)
    for j = idx_show, plot(tq, q2_nl(j,:), '--', 'LineWidth', 1.0); end
end

figure('Name','Base-only: Quat','Color','w');
% plot(tq, z_tip_lin*100, 'LineWidth', 1.8); hold on;
if ~isempty(xn), plot(tq, eul_full, 'LineWidth', 1.8); end
grid on; xlabel('t [s]'); ylabel('eul_{tip} [% of m]');   % (raw meters *100 = “percent of m”)
legend('phi', 'theta', 'psi');
ttl = 'Structural-only: Quat Eul at tip';
title(ttl);

% (Optional) also expose strain/velocity fields via phi1/phi2 if you want:
% x1_lin(s,t) = beam.phi1 * q1_lin(:,t);  ( 6*Ne x 1 per t ) – similar for x2.

end

function dx = rhs_full_base(t, x, Om, Sig, G1, G2, Gg, Gxi, eta_fun)
% Full structural ROM: linear part + quadratic contractions
Nm = size(Om,1);
q1 = x(1:Nm);
q2 = x(Nm+1:2*Nm);
q_xi = x(2*Nm+1: 3*Nm+1);
% Quadratic/bilinear terms:
% g11 = contractGamma(G1, q1, q1);  % i: sum_{jk} G1(i,j,k) q1_j q1_k
% g22 = contractGamma(G2, q2, q2);  % i: sum_{jk} G2(i,j,k) q2_j q2_k
% g21 = contractGamma(G2, q2, q1);  % i: sum_{jk} G2(i,j,k) q2_j q1_k  (bilinear)

%%
% ---------- Γ₁ & Γ₂  contractions  (vectorised) ---------------------------
%    −Σ_l  q1_l Γ₁(:,:,l) q1    = - Σ_l  diag(q1) *
%                                  vec(Gamma1)^T * kron(q1,I)
% Use pagemtimes (R2023b) :
G1q1 = squeeze(pagemtimes(G1, 'none', q1, 'none'));  % Nm×Nm×1  <- Σ over 3rd dim
N1_G1 = -(G1q1(:,:,1) * q1);
% N(idx.q1) = N(idx.q1) + N1_G1;
G2q2 = squeeze(pagemtimes(G2, 'none', q2, 'none'));
N1_G2 = -(G2q2(:,:,1) * q2);
% N(idx.q2) = N(idx.q2) + N1_G2;

% % ---------- Γ_g  term -----------------------------------------------------
Ggqx = squeeze(pagemtimes(Gg, 'none', q_xi, 'none'));
% Ggqx = squeeze(pagemtimes(permute(Gg, [1 3 2]), 'none', q_xi, 'none'));
N1_Gg =  (Ggqx(:,:,1) * q_xi);
% N(idx.q2) = N(idx.q2) + N1_Gg;
% 
% % ---------- Γ_xi  term  (Nm+1×1) -----------------------------------------
% disp(size(Gxi))
% Gxiq1 = squeeze(pagemtimes(Gxi, 'none', q1, 'none'));
Gxiq1 = squeeze(pagemtimes(Gxi, 'none', q_xi, 'none'));
% N_qxi =  (Gxiq1(:,:,1) * q_xi);
N_qxi =  (Gxiq1(:,:,1) * q1);
% disp(size(Gg))
% disp(size(q_xi))

% Gxiq1 = squeeze(pagemtimes(permute(Gxi, [1 3 2]), 'none', q_xi, 'none'));
% N_qxi =  (Gxiq1(:,:,1) * q1);

% ---------- Γ2^T q1  (Eq 9b) ----------------------------------------------
GT = squeeze(pagemtimes(G2, 'transpose', q1, 'none'));
N_GT =  GT(:,:,1)*q2;     % note: transpose on first arg
%%

% dx1 = Om*q2 + Sig*q1 - g11 - g22 + eta_fun(t);
dx1 = Om*q2 + Sig*q1 + N1_G1 + N1_G2 + N1_Gg + eta_fun(t);
% dx2 = -Om*q1 + g21;
dx2 = -Om*q1 + N_GT;
dx3 = N_qxi;


dx  = [dx1; dx2; dx3];
end

function y = contractGamma(G, a, b)
% y_i = Σ_{j,k} G(i,j,k) a_j b_k
% Accepts G as Nm×Nm×Nm or tries to adapt common permutations.
sz = size(G);
Nm = numel(a);
y  = zeros(Nm,1);
if isequal(sz, [Nm Nm Nm])
    for i = 1:Nm
        Gi = squeeze(G(i,:,:));
        y(i) = a.' * (Gi * b);
    end
else
    % try to reshape/permute if user stores as Nm×Nm×Nm in another order
    if numel(G)==Nm^3
        G = reshape(G, [Nm Nm Nm]);
        for i = 1:Nm
            Gi = squeeze(G(i,:,:));
            y(i) = a.' * (Gi * b);
        end
    else
        error('Gamma tensor has unexpected size.');
    end
end

% ---------- extract modal sub‑vectors (no find/ismember in loop) ---------
q1      = S(idx.q1 );            % Nm
q2      = S(idx.q2 );            % Nm
q_xi    = S(idx.qxi);            % Nm+1
% q_Gamma = S(idx.qGam);           %#ok<NASGU>  % Na2 (unused in open loop)
% chi_hat = S(idx.chi );           %#ok<NASGU>  % 3   (unused in open loop)

Nm = numel(idx.q1);
Nx = numel(S);
N_terms = zeros(Nx,1);


end

function tip = pick_tip_node(fem)
% pick node with maximum +Y (consistent with a half-span wing along +Y)
[~, tip] = max(fem.coordinates(:,2));
end

function row = dof_row_in_phi0(keepDofs, node, dof)
% map node,dof -> row index in phi0 (flexible-only slice)
fullRow = (node-1)*6 + dof;     % (u,v,w,rx,ry,rz) = (1..6)
[tf, row] = ismember(fullRow, keepDofs);
if ~tf, row = 0; end
end

function row = nearest_flexible_z(keepDofs)
% pick first available w-dof in keepDofs as fallback
zRows = keepDofs(mod(keepDofs-1,6)==2);   % DOF 3 -> (index-1) mod 6 == 2
if isempty(zRows)
    row = 1;
else
    firstFull = zRows(1);
    [~, row] = ismember(firstFull, keepDofs);
end
end

function L = make_leg(idx_show, hasNL)
tags = compose('mode %d', idx_show);
if hasNL, L = [tags; strcat(tags,' (NL)')]; L = L(:).'; else, L = tags; end
end


function v = getfieldwithdefault(s, f, default)
    if isfield(s,f); v = s.(f); else; v = default; end
end