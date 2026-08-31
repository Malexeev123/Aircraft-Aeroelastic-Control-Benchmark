function [setupDir, forMatlabDir] = resolveSetupDir(setupDirIn)
%RESOLVESETUPDIR Locate the setup folder that contains for_matlab/sim_bundle.mat.
%
% Accepted inputs:
%   1) Exact setup folder:
%        .../wing_only/pt_U040_alpha_m01
%
%   2) Exact MATLAB artifact folder:
%        .../wing_only/pt_U040_alpha_m01/for_matlab
%
%   3) Parent folder containing setup folders:
%        .../wing_only
%
% Point-folder naming convention:
%        pt_U040_alpha_m01
%        pt_U040_alpha_p00
%        pt_U040_alpha_p02
%
% Legacy date-folder convention remains supported:
%        20260610
%        20260610_153012
%
% Search priority for parent-folder input:
%   1. deterministic point folders: pt_U*_alpha_*
%   2. legacy date folders:        yyyymmdd or yyyymmdd_HHMMSS
%   3. any other folder containing for_matlab/sim_bundle.mat

    setupDirIn = char(setupDirIn);

    if exist(setupDirIn,'dir') ~= 7
        error('resolveSetupDir:MissingDir', ...
            'Setup directory does not exist:\n  %s', setupDirIn);
    end

    % ---------------------------------------------------------------------
    % Case 1: input already points to the for_matlab artifact folder.
    % ---------------------------------------------------------------------
    [parentDir, leaf] = fileparts(setupDirIn);

    if strcmpi(leaf, 'for_matlab')
        setupDir     = parentDir;
        forMatlabDir = setupDirIn;

        requireSimBundle(forMatlabDir);
        return
    end

    % ---------------------------------------------------------------------
    % Case 2: input is an exact setup folder.
    % ---------------------------------------------------------------------
    fm = fullfile(setupDirIn, 'for_matlab');

    if exist(fullfile(fm,'sim_bundle.mat'),'file') == 2
        setupDir     = setupDirIn;
        forMatlabDir = fm;
        return
    end

    % ---------------------------------------------------------------------
    % Case 3: input is a parent folder containing point/date subfolders.
    % ---------------------------------------------------------------------
    hits = dir(fullfile(setupDirIn, '*', 'for_matlab', 'sim_bundle.mat'));

    if isempty(hits)
        error('resolveSetupDir:MissingBundle', ...
            ['Could not find sim_bundle.mat under:\n  %s\n\n', ...
             'Expected either:\n  %s\n\n', ...
             'or one level below it, for example:\n  %s'], ...
             setupDirIn, ...
             fullfile(setupDirIn,'for_matlab','sim_bundle.mat'), ...
             fullfile(setupDirIn,'pt_U040_alpha_m01','for_matlab','sim_bundle.mat'));
    end

    % Rank candidates so point-named library folders win over old date folders.
    nHits = numel(hits);
    rank  = zeros(1,nHits);
    dnum  = zeros(1,nHits);

    for i = 1:nHits
        fm_i    = hits(i).folder;
        setup_i = fileparts(fm_i);
        [~, runLeaf] = fileparts(setup_i);

        if isPointRunFolder(runLeaf)
            rank(i) = 3;
        elseif isDateRunFolder(runLeaf)
            rank(i) = 2;
        else
            rank(i) = 1;
        end

        dnum(i) = hits(i).datenum;
    end

    % Highest rank first, newest within the same class.
    score = rank(:)*1e10 + dnum(:);
    [~, ix] = max(score);

    forMatlabDir = hits(ix).folder;
    setupDir     = fileparts(forMatlabDir);

    if nHits > 1
        warning('[resolveSetupDir] Parent input matched %d setup folders. Selected:\n  %s', ...
            nHits, setupDir);
    else
        warning('[resolveSetupDir] Resolved nested setup folder:\n  input : %s\n  using : %s', ...
            setupDirIn, setupDir);
    end
end

function requireSimBundle(forMatlabDir)
%REQUIRESIMBUNDLE Guard exact for_matlab inputs.

    fn = fullfile(forMatlabDir, 'sim_bundle.mat');

    if exist(fn,'file') ~= 2
        error('resolveSetupDir:MissingBundle', ...
            'Missing sim_bundle.mat in:\n  %s', forMatlabDir);
    end
end

function tf = isPointRunFolder(name)
%ISPOINTRUNFOLDER Match deterministic ROM-library point folders.
%
% Examples:
%   pt_U040_alpha_m01
%   pt_U040_alpha_p00
%   pt_U040p5_alpha_m01p5

    name = char(name);

    expr = '^pt_U[0-9]+(p[0-9]+)?_alpha_[pm][0-9]+(p[0-9]+)?$';
    tf = ~isempty(regexp(name, expr, 'once'));
end

function tf = isDateRunFolder(name)
%ISDATERUNFOLDER Match legacy date-based setup folders.
%
% Examples:
%   20260610
%   20260610_153012

    name = char(name);

    exprDateOnly = '^[0-9]{8}$';
    exprDateTime = '^[0-9]{8}_[0-9]{6}$';

    tf = ~isempty(regexp(name, exprDateOnly, 'once')) || ...
         ~isempty(regexp(name, exprDateTime, 'once'));
end