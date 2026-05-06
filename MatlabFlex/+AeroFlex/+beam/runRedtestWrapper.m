function red = runRedtestWrapper(fem, aero, phi0, Nm, cfg, obj)
%RUNREDTESTWRAPPER  Call legacy compute_pazy_redtest with all inputs set.
%   • Prefers obj.Kglobal/Mglobal already prepared by BeamModel
%   • Otherwise loads from SHARPy (beam_modal_analysis -> data.h5 fallback)

% ----- inputs & defaults --------------------------------------------------
aoa_deg   = cfg.flight.aoa_deg;
% disp(aoa_deg)
debugPlot = isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level >= 2;
aoa       = aoa_deg; % legacy wants degrees

% ----- pick Mred/Kred -----------------------------------------------------
Kred = []; Mred = [];
if isfield(obj,'Kglobal') && ~isempty(obj.Kglobal) && isfield(obj,'Mglobal') && ~isempty(obj.Mglobal)
    Kred = obj.Kglobal;  Mred = obj.Mglobal;
else
    % try beam_modal_analysis first
    [~, ~, ~, beam_modal_dir] = resolve_sharpy_beam_roots(cfg);
    [Kred, Mred] = try_load_KM_from_beam_modal(beam_modal_dir, cfg.struct.rigid_body_modes);

    % fallback: savedata data.h5
    if isempty(Kred) || isempty(Mred)
        data_h5 = fullfile(cfg.paths.sharpy.savedata, [cfg.case_name, '.data.h5']);
        [Kred, Mred] = try_load_KM_from_data_h5(data_h5, cfg.struct.rigid_body_modes);
    end
end

% sanity
assert(~isempty(Kred) && ~isempty(Mred), 'runRedtestWrapper: cannot obtain K/M');

% ----- eigen-freq vector (needed for derivative shapes) -------------------
Omega = obj.Omega; w0 = diag(Omega);    % keep your legacy signature (compute_pazy_redtest expects it)

% ----- chain for phi2 accumulation ---------------------------------------
chain = AeroFlex.core.buildChainFrom3NodeElements(fem);
% disp('fem.coordinates');
% disp(fem.coordinates);
% %%
% % ----- New Test ------------------------------------------
% [ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = AeroFlex.beam.compute_instrinsic_modes( ...
%       fem, aero, phi0, Mred, Kred, Nm, aoa, chain, w0, debugPlot, cfg.struct.rigid_body_modes, cfg.NetworkPath);
% 
% error('Proceed with legacy')

% ----- call the legacy function ------------------------------------------
[ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = AeroFlex.beam.compute_pazy_redtest( ...
      fem, aero, phi0, Mred, Kred, Nm, aoa, chain, w0, debugPlot, cfg.struct.rigid_body_modes, cfg.NetworkPath);

% [ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = AeroFlex.beam.compute_modal( ...
%       fem, aero, phi0, Mred, Kred, Nm, aoa, chain, w0, debugPlot, cfg.struct.rigid_body_modes, cfg.NetworkPath);
%% Test
 % phi1_sA_load = load("eigenvectors_rigid.dat");   
 % phi1_sA = phi1_sA_load(97:102,:);
            %%

red = struct('ModeVars_continuous', ModeVars_continuous, ...
             'ModeVars_discrete',  ModeVars_discrete, ...
             'Beam_Props',         Beam_Props, ...
             'phi1_sA',            phi1_sA, ...
             'phi2_sA',            phi2_sA);
end

%% ------------- helpers duplicated from BeamModel locals -----------------

function [root, cases_dir, savedata_dir, beam_modal_dir] = resolve_sharpy_beam_roots(cfg)
    assert(isfield(cfg,'paths') && isfield(cfg.paths,'sharpy') && isfield(cfg,'case_name'), ...
        'cfg.paths.sharpy and cfg.case_name are required.');
    root = cfg.paths.sharpy.root;
    case_name = cfg.case_name;
    if isfield(cfg.paths.sharpy,'cases') && ~isempty(cfg.paths.sharpy.cases)
        cases_dir = cfg.paths.sharpy.cases;
    else
        cases_dir = fullfile(root,'cases',case_name);
    end
    if isfield(cfg.paths.sharpy,'savedata') && ~isempty(cfg.paths.sharpy.savedata)
        savedata_dir = cfg.paths.sharpy.savedata;
    else
        savedata_dir = fullfile(root,'output',case_name,'savedata');
    end
    if isfield(cfg.paths.sharpy,'beam_modal') && ~isempty(cfg.paths.sharpy.beam_modal)
        beam_modal_dir = cfg.paths.sharpy.beam_modal;
    else
        beam_modal_dir = fullfile(root,'output',case_name,'beam_modal_analysis');
    end
end

function [K,M] = try_load_KM_from_beam_modal(beam_modal_dir, use_rigid)
    K=[]; M=[];
    if ~isfolder(beam_modal_dir), return; end
    if use_rigid
        Kfile = fullfile(beam_modal_dir,'Kglobal_rigid.dat');
        Mfile = fullfile(beam_modal_dir,'Mglobal_rigid.dat');
        if ~isfile(Kfile), Kfile = fullfile(beam_modal_dir,'Kglobal.dat'); end
        if ~isfile(Mfile), Mfile = fullfile(beam_modal_dir,'Mglobal.dat'); end
    else
        Kfile = fullfile(beam_modal_dir,'Kglobal.dat');
        Mfile = fullfile(beam_modal_dir,'Mglobal.dat');
    end
    if isfile(Kfile) && isfile(Mfile)
        K = load(Kfile); M = load(Mfile);
    end
end

function [K,M] = try_load_KM_from_data_h5(data_h5, use_rigid)
    K=[]; M=[];
    info = h5info(data_h5);
    [K,~] = try_read_first(data_h5, {'/data/structural/Kglobal','/data/structural/K_global','/Kglobal','/K_global'});
    [M,~] = try_read_first(data_h5, {'/data/structural/Mglobal','/data/structural/M_global','/Mglobal','/M_global'});
    if isempty(K) || isempty(M)
        if isempty(K), [K,~] = discover_square(data_h5, info, 'K'); end
        if isempty(M), [M,~] = discover_square(data_h5, info, 'M'); end
    end
    if use_rigid
        [Kr,~] = try_read_first(data_h5, {'/data/structural/Kglobal_rigid','/data/structural/K_rigid'});
        [Mr,~] = try_read_first(data_h5, {'/data/structural/Mglobal_rigid','/data/structural/M_rigid'});
        if ~isempty(Kr) && ~isempty(Mr), K=Kr; M=Mr; end
    end
end

function [X, path] = try_read_first(h5fn, paths)
    X = []; path = '';
    for i=1:numel(paths)
        p = paths{i};
        if has_dataset(h5fn, p)
            X = h5read(h5fn, p); path = p; return;
        end
    end
end

function tf = has_dataset(h5fn, path)
    tf=false; try, h5info(h5fn, path); tf=true; catch, tf=false; end
end

function [X, path] = discover_square(h5fn, info, keyChar)
    X=[]; path='';
    for d=1:numel(info.Datasets)
        ds = info.Datasets(d);
        if contains(lower(ds.Name), lower(keyChar))
            sz = fliplr(ds.Dataspace.Size);
            if numel(sz)==2 && sz(1)==sz(2) && sz(1)>=6
                trial = h5read(h5fn, ['/' ds.Name]);
                if isfloat(trial), X=trial; path=['/' ds.Name]; return; end
            end
        end
    end
    groups = info.Groups;
    while ~isempty(groups)
        g = groups(1); groups(1) = [];
        for d=1:numel(g.Datasets)
            ds = g.Datasets(d);
            if contains(lower(ds.Name), lower(keyChar))
                fullp = [g.Name '/' ds.Name];
                sz = fliplr(ds.Dataspace.Size);
                if numel(sz)==2 && sz(1)==sz(2) && sz(1)>=6
                    trial = h5read(h5fn, fullp);
                    if isfloat(trial), X=trial; path=fullp; return; end
                end
            end
        end
        groups = [groups; g.Groups]; 
    end
end

% function red = runRedtestWrapper(fem, aero, phi0, Nm, cfg, obj)
% %RUNREDTESTWRAPPER  Call legacy compute_pazy_redtest with all inputs set.
% %   • Builds Mred & Kred from section props if they are not on disk
% %   • Generates chain array from FEM connectivity
% %   • Uses cfg.trim.aoa_deg (default 1 deg)
% %   • Disables debug plots unless cfg.debug.level ≥ 2
% 
% % ----- defaults ----------------------------------------------------------
% aoa_deg   = cfg.flight.aoa_deg;                  % 1 deg by default
% debugPlot = cfg.debug.level >= 2;
% aoa       = aoa_deg;                           % legacy wants degrees
% 
% Omega = obj.Omega;
% w0 = diag(Omega);
% Mred = obj.Mglobal; 
% Kred = obj.Kglobal; 
% % Mred = obj.Mhat; 
% % Kred = obj.Khat; 
% 
% % ----- reduced matrices (build once, save under /cache) ------------------
% if Nm == 20
% redCache = fullfile(cfg.paths.cache,'MKred.mat');
% else
%     redCache = fullfile(cfg.paths.cache,cfg.case_name,'MKred.mat');
% end
% if isfile(redCache)
%     S = load(redCache,'Mred','Kred');
%     Mred = S.Mred; Kred = S.Kred;
% else
%     % [Kred,Mred] = buildReducedKM(fem); % user helper
%     save(redCache,'Kred','Mred','-v7.3');
% end
% % Mred = obj.Mhat; 
% % Kred = obj.Khat;
% % Mred = eye(size(Mred)); 
% % Kred = obj.Khat;
% % ----- eigen-freq vector (needed for derivative shapes) ------------------
% % w0 = sqrt(diag(phi0.'*Kred*phi0));   % ω = √(φᵀKφ) for mass-norm φ0
% 
% % ----- chain for phi2 accumulation --------------------------------------
% chain = AeroFlex.core.buildChainFrom3NodeElements(fem);
% 
% % if cfg.struct.rigid_body_modes 
% % 
% % end
% 
% % ----- call the legacy function -----------------------------------------
% [ModeVars_continuous, ModeVars_discrete, Beam_Props, phi1_sA, phi2_sA] = AeroFlex.beam.compute_pazy_redtest( ...
%           fem, aero, phi0, Mred, Kred, Nm, aoa, chain, w0, debugPlot, cfg.struct.rigid_body_modes, cfg.NetworkPath);
% 
% %% Temp 
% % rigid_phi = load('eigenvectors_rigid.dat');
% % % rigid_phi = load('data/eigenvectors_rigid');
% % phi1_sA = rigid_phi(97:102,:);
% %%
% red = struct('ModeVars_continuous', ModeVars_continuous, 'ModeVars_discrete',ModeVars_discrete, 'Beam_Props', Beam_Props, 'phi1_sA', phi1_sA, ...
%     'phi2_sA', phi2_sA);
% end
