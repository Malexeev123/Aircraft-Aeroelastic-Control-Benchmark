function cfg = defaultLibraryConfig(cfg)
%DEFAULTLIBRARYCONFIG Fill missing ROM-library scheduling options.
%
% The scheduler uses a local library of MATLAB-assembled ROM points.  Each
% point is generated at a SHARPy/MATLAB flight condition and stores the
% quantities consumed by ROMIntegrator, SimRunner, PlantRunTime and trim.
%
% Required library coordinates for the first implementation:
%   mu = [U_inf, alpha_deg]
%
% The library is intentionally kept independent of nMHE/nMPC.  The control
% horizon uses the frozen scheduled ROM supplied at the current major step.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

if ~isfield(cfg,'library') || ~isstruct(cfg.library)
    cfg.library = struct();
end

L = cfg.library;

if ~isfield(L,'enable'),          L.enable = false; end
if ~isfield(L,'method'),          L.method = 'stage1_direct'; end
if ~isfield(L,'muNames'),         L.muNames = {'U_inf','alpha_deg'}; end
if ~isfield(L,'noExtrapolate'),   L.noExtrapolate = true; end
if ~isfield(L,'updateMode'),      L.updateMode = 'perPlantStep'; end
if ~isfield(L,'freezeInOptimizer'), L.freezeInOptimizer = true; end
if ~isfield(L,'interpTol'),       L.interpTol = 1e-10; end
if ~isfield(L,'debug'),           L.debug = false; end

if ~isfield(L,'root') || isempty(L.root)
    if isfield(cfg,'paths') && isfield(cfg.paths,'sharpy') && isfield(cfg.paths.sharpy,'root')
        L.root = fullfile(cfg.paths.sharpy.root,'rom_library');
    else
        L.root = fullfile(pwd,'rom_library');
    end
end

if ~isfield(L,'name') || isempty(L.name)
    if isfield(cfg,'case_name')
        L.name = sprintf('%s_U_alpha_library', cfg.case_name);
    else
        L.name = 'U_alpha_library';
    end
end

if ~isfield(L,'path') || isempty(L.path)
    L.path = fullfile(L.root,L.name,'library_stage1_direct.mat');
end

if ~isfield(L,'U_grid') || isempty(L.U_grid)
    if isfield(cfg,'flight') && isfield(cfg.flight,'U_inf')
        L.U_grid = cfg.flight.U_inf * [0.90 1.00 1.10];
    else
        L.U_grid = [36 40 44];
    end
end

if ~isfield(L,'alpha_grid_deg') || isempty(L.alpha_grid_deg)
    if isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
        L.alpha_grid_deg = cfg.flight.aoa_deg + [-4 -2 0 2 4];
    else
        L.alpha_grid_deg = [-4 -2 0 2 4];
    end
end

if ~isfield(L,'stateTransform')
    L.stateTransform = struct();
end
if ~isfield(L.stateTransform,'enable'), L.stateTransform.enable = false; end
if ~isfield(L.stateTransform,'blocks'), L.stateTransform.blocks = {'qGam'}; end

% Clamp-load reference correction. Set r_A_CG_B from CG to wing clamp point A.
if ~isfield(L,'loads')
    L.loads = struct();
end
if ~isfield(L.loads,'r_A_CG_B'), L.loads.r_A_CG_B = [0;0;0]; end
if ~isfield(L.loads,'shiftWingMomentToCG'), L.loads.shiftWingMomentToCG = true; end

cfg.library = L;
end
