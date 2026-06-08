function ROMlib = buildLibraryFromSetups(setupDirs, varargin)
%BUILDLIBRARYFROMSETUPS Build a ROM library from completed sim_init folders.
%
% setupDirs must point to folders containing for_matlab/sim_bundle.mat,
% aero_bundle.mat, beam_bundle.mat and base_bundle.mat.  This function does
% not run SHARPy; it packages existing MATLAB assemblies into the scheduler
% database.  Use it after SHARPy and sim_init have been run at all grid
% points.

p = inputParser;
p.addParameter('library_name','U_alpha_library',@(s)ischar(s)||isstring(s));
p.addParameter('save_path','',@(s)ischar(s)||isstring(s));
p.addParameter('method','stage1_direct',@(s)ischar(s)||isstring(s));
p.addParameter('body_case','',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

if ischar(setupDirs) || isstring(setupDirs)
    setupDirs = cellstr(setupDirs);
end

points = repmat(struct(),0,1);
for i = 1:numel(setupDirs)
    sd = char(setupDirs{i});
    fm = fullfile(sd,'for_matlab');
    assert(exist(fullfile(fm,'sim_bundle.mat'),'file')==2, 'Missing sim_bundle.mat in %s.', fm);
    assert(exist(fullfile(fm,'aero_bundle.mat'),'file')==2, 'Missing aero_bundle.mat in %s.', fm);
    assert(exist(fullfile(fm,'beam_bundle.mat'),'file')==2, 'Missing beam_bundle.mat in %s.', fm);
    assert(exist(fullfile(fm,'base_bundle.mat'),'file')==2, 'Missing base_bundle.mat in %s.', fm);

    Ssim = load(fullfile(fm,'sim_bundle.mat'));
    Saer = load(fullfile(fm,'aero_bundle.mat'));
    Sbem = load(fullfile(fm,'beam_bundle.mat'));
    Sbas = load(fullfile(fm,'base_bundle.mat'));

    cfg = Ssim.sim_config.cfg;
    trim = Ssim.sim_config.trim;
    if isempty(opt.body_case) && isfield(Ssim.sim_config,'body_case')
        body_case = Ssim.sim_config.body_case;
    else
        body_case = opt.body_case;
    end

    points(i,1) = AeroFlex.sched.extractPoint(cfg, Sbem.beam, Saer.aero, Sbas.base, trim, ...
        'source_dir', sd, 'body_case', body_case); %#ok<AGROW>
end

ROMlib = struct();
ROMlib.name = char(opt.library_name);
ROMlib.method = char(opt.method);
ROMlib.muNames = {'U_inf','alpha_deg'};
ROMlib.points = points;
ROMlib.mu = reshape([points.mu],2,[]).';
ROMlib.created = datestr(now,'yyyy-mm-dd HH:MM:SS');
ROMlib.version = 1;
ROMlib.notes = ['Stage-1 direct interpolation of L and nonlinear parConst fields. ', ...
                'All Gamma tensors and force maps are interpolated; no coordinate ', ...
                'transform is applied unless finalizeCompatibleLibrary is called.'];

AeroFlex.sched.validateLibrary(ROMlib);

if ~isempty(opt.save_path)
    save_path = char(opt.save_path);
    outDir = fileparts(save_path);
    if ~isempty(outDir) && exist(outDir,'dir')~=7, mkdir(outDir); end
    save(save_path,'ROMlib','-v7.3');
end
end
