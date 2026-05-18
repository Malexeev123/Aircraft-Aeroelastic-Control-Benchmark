
classdef AeroROM
    % AEROROM  Build aerodynamic ROM using SHARPy case + output trees.
    %
    % Reads:
    %   - cases/<case>/<case>.aero.h5              (aero mesh/settings)
    %   - output/<case>/savedata/<case>.linss.h5   (linearisation vectors)
    %   - output/<case>/savedata/<case>.data.h5    (coupling matrices)
    %   - output/<case>/savedata/*.rom.h5          (Krylov ROM, DISCRETE)
    %   - output/<case>/savedata/*aerorob*.h5      (projection V/gain; optional)
    %
    % Produces:
    %   ROM_cts, ROM_dsc, forceMap, DataMatrix (with KrylovSharpyContSS), gust_input,
    %   aeroMesh and bookkeeping paths.
    %
    % Notes:
    %   * No cache: deterministic re-build from SHARPy files.
    %   * cfg.aero.discrete == true -> puts discrete model in ROM_cts (legacy naming).

    % properties (SetAccess = private)
    properties 
        Na
        Nm
        case_name
        aeroMesh
        ROM_cts
        ROM_dsc
        forceMap
        DataMatrix
        gust_input
        NetworkPath

        % Path bookkeeping (resolved on construction)
        sharpy_root
        cases_dir      % .../cases/<case>
        savedata_dir   % .../output/<case>/savedata
        save_pmor_dir
        aero_h5        % cases/<case>/<case>.aero.h5
        linss_h5       % savedata/<case>.linss.h5
        data_h5        % savedata/<case>.data.h5
        rom_h5         % savedata/<something>.rom.h5 (pattern match)
        vproj_h5       % savedata/*aerorob*.h5 (optional)
    end

    methods
        function obj = AeroROM(cfg)
            % Sizes and IDs
            obj.Na          = cfg.Na;
            obj.Nm          = cfg.Nm;
            obj.case_name   = cfg.case_name;
            obj.NetworkPath = getfield(cfg, 'NetworkPath'); %#ok<GFLD>
            % disp(obj.case_name)
            % Resolve SHARPy paths
            [obj.sharpy_root, obj.cases_dir, obj.savedata_dir, obj.save_pmor_dir] = ...
                AeroFlex.aero.AeroROM.resolve_sharpy_roots(cfg);

            % Compose expected file paths (with robust fallbacks)
            obj.aero_h5  = AeroFlex.aero.AeroROM.find_file_with_fallback( ...
                fullfile(obj.cases_dir, [obj.case_name,'.aero.h5']), ...
                {obj.cases_dir}, '*.aero.h5', true, ...
                'Aero mesh (*.aero.h5) not found under cases');

            obj.linss_h5 = AeroFlex.aero.AeroROM.find_file_with_fallback( ...
                fullfile(obj.savedata_dir, [obj.case_name,'.linss.h5']), ...
                {obj.savedata_dir}, '*.linss.h5', true, ...
                '.linss.h5 not found under savedata');

            obj.data_h5  = AeroFlex.aero.AeroROM.find_file_with_fallback( ...
                fullfile(obj.savedata_dir, [obj.case_name,'.data.h5']), ...
                {obj.savedata_dir}, '*.data.h5', true, ...
                '.data.h5 not found under savedata');

            % ROM now lives in SAVEDATA. Support a few common names.
            obj.rom_h5   = AeroFlex.aero.AeroROM.find_file_with_fallback( ...
                fullfile(obj.savedata_dir, [obj.case_name,'.rom.h5']), ...
                {obj.savedata_dir}, {'*.rom.h5','*krylov*.rom.h5'}, true, ...
                'Krylov ROM (*.rom.h5) not found under savedata');

            % Projection V often written next to ROM in savedata (optional).
            obj.vproj_h5 = AeroFlex.aero.AeroROM.find_file_with_fallback( ...
                fullfile(obj.save_pmor_dir, [obj.case_name,'_krylov_aerorob.h5']), ...
                {obj.save_pmor_dir}, {'*krylov*_aerorob*.h5','*aerorob*.h5'}, false, ...
                'Projection V file not found (will continue without it)');

            % 1) Read aero mesh (cases/<case>/<case>.aero.h5)
            obj.aeroMesh = AeroFlex.aero.readAEROFromH5(obj.aero_h5);

            % 2) Build ROM & collect matrices/vectors (ROM from SAVEDATA)
            [ROM_krylov_SS, forces_aero_beam_dof, mode_shape, ...
             ROM_krylov_Disc, forces_aero, forces_structure, DataMatrix, etas] = ...
                AeroFlex.aero.AeroDiscToCont(obj, cfg); % <- no rom_h5 arg; uses obj.rom_h5
            % 3) Choose representation
            if isfield(cfg,'aero') && isfield(cfg.aero,'discrete') && cfg.aero.discrete
                obj.ROM_cts = ROM_krylov_Disc;
                obj.ROM_dsc = ROM_krylov_Disc;
            else
                if isfield(DataMatrix,'KrylovSharpyContSS') && ~isempty(DataMatrix.KrylovSharpyContSS)
                    obj.ROM_cts = DataMatrix.KrylovSharpyContSS;
                else
                    obj.ROM_cts = ROM_krylov_SS;
                end
                obj.ROM_dsc = ROM_krylov_Disc;
            end

            % 4) Force map
            obj.forceMap = AeroFlex.aero.buildForceMap( ...
                obj.ROM_cts, obj.Na, obj.Nm, ...
                getfield(cfg.aero,'discrete'), ...
                getfield(cfg.ctrl,'n_surf'), ...
                getfield(cfg.ctrl,'var_per'), ...
                getfield(cfg.gust,'on'), ...
                1, ... %cfg.force_map.scale
                obj.savedata_dir);
            cfg.sim.dt = obj.ROM_dsc.Ts;
            % 5) Gust input 
            if cfg.gust.on
                gust = AeroFlex.aero.get_1minuscosine_gust_input( ...
                    cfg.gust.gust_length, ...
                    cfg.gust.gust_offset, ...
                    cfg.gust.gust_intensity, ...
                    cfg.sim.dt, ...
                    cfg.flight.U_inf, ...
                    cfg.sim.t_end);
                % obj.gust_input = gust(:,2)/cfg.force_map.scale;
                % obj.gust_input = 0.525 * gust(:,2);
                % obj.gust_input = 0.515 * gust(:,2);
                % obj.gust_input = 0.5 * gust(:,2);
                % obj.gust_input =  10*gust(:,2);
                % obj.gust_input =  5*gust(:,2);
                obj.gust_input =  gust(:,2);
            else
                obj.gust_input = zeros(round(cfg.sim.t_end/cfg.sim.dt),1);
            end

            % 6) Stash DataMatrix (+ extras)
            % DataMatrix.forces_aero_beam_dof = forces_aero_beam_dof/cfg.force_map.scale;
            DataMatrix.forces_aero_beam_dof = forces_aero_beam_dof/cfg.force_map.scale;
            DataMatrix.forces_aero          = forces_aero;
            DataMatrix.forces_structure     = forces_structure;
            DataMatrix.etas_sharpy          = etas;
            obj.DataMatrix = DataMatrix;
            % disp(fieldnames(DataMatrix));

            
            
            if isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level>=1
                fprintf('[AeroROM] Built: n=%d, m=%d, p=%d, Ts=%.6g\n', ...
                    size(obj.ROM_dsc.A,1), size(obj.ROM_dsc.B,2), size(obj.ROM_dsc.C,1), obj.ROM_dsc.Ts);
            end
        end
    end

    methods(Static, Access=private)
        function [root, cases_dir, savedata_dir, save_pmor_dir] = resolve_sharpy_roots(cfg)
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

            save_pmor_dir = fullfile(root,'output',case_name,'save_pmor_data');
        end

        function fn = find_file_with_fallback(primary_guess, search_roots, patterns, hard, errmsg)
            if nargin<5, errmsg = 'file not found'; end
            if nargin<4, hard = true; end
            if ischar(patterns) || isstring(patterns), patterns = {char(patterns)}; end

            if ~isempty(primary_guess) && isfile(primary_guess)
                fn = primary_guess; return;
            end

            for r = 1:numel(search_roots)
                root = search_roots{r};
                if ~isfolder(root), continue; end
                for p = 1:numel(patterns)
                    D = dir(fullfile(root, patterns{p}));
                    if ~isempty(D)
                        fn = fullfile(D(1).folder, D(1).name);
                        return;
                    end
                end
            end

            if hard
                error('%s (patterns %s not found in: %s)', errmsg, strjoin(patterns,', '), strjoin(search_roots,', '));
            else
                fn = '';
            end
        end
    end
end



% classdef AeroROM
%     properties (SetAccess = private)
%         Na
%         Nm, case_name
%         aeroMesh           % struct from readAEROFromH5
%         ROM_cts            % continuous‑time SS
%         ROM_dsc            % discrete‑time SS
%         forceMap           % struct with D0, D1, B0, B1, CGamma, ...
%         cacheFile
%         dataPath
%         aeroFile
%         DataMatrix
%         gust_input
%         NetworkPath
%     end
%     methods
%         function obj = AeroROM(cfg)
%             obj.Na       = cfg.Na;
%             obj.Nm       = cfg.Nm;
%             obj.NetworkPath = cfg.NetworkPath;
%             if obj.Nm ==20
%                 obj.cacheFile = fullfile(cfg.paths.cache,'rom.mat');
%                 obj.dataPath  = cfg.paths.data;
%                 obj.aeroFile = "pazy_ROM.aero.h5";
%                 obj.case_name = 'pazy_ROM';
%             else
%                 cfg.case_name = append('pazy_ROM_',string(cfg.Nm),'M');
%                 obj.case_name = cfg.case_name;
%                 obj.aeroFile =append(obj.case_name,'.aero.h5');
% 
%                 if not(isfolder(fullfile(cfg.paths.cache,obj.case_name)))
%                     mkdir(fullfile(cfg.paths.cache,obj.case_name))
%                 end            
%                 if not(isfolder(fullfile(cfg.paths.data,obj.case_name)))
%                     mkdir(fullfile(cfg.paths.data,obj.case_name))
%                 end 
%                 obj.cacheFile= fullfile(cfg.paths.cache,obj.case_name,'rom.mat');
%                 obj.dataPath  = fullfile(cfg.paths.data,obj.case_name);
% 
%             end
% 
% 
% 
% 
%             obj = obj.loadOrCompute(cfg);
% 
%             % Testing the Projection for Trim of aero states
%             aeroZ = null(obj.DataMatrix.aero_gamma_state.project_V.r);
%             aeroPz = aeroZ/(aeroZ.'*aeroZ)*aeroZ.';
% 
%             aeroZr = orth(obj.DataMatrix.aero_gamma_state.project_V.r.');
%             aeroPr = aeroZr/(aeroZr.'*aeroZr)*aeroZr.';
%             obj.DataMatrix.aeroPzT = obj.DataMatrix.aero_gamma_state.project_V.r*aeroPz*obj.DataMatrix.aero_gamma_state.project_V.r.';
%             obj.DataMatrix.aeroPrT = obj.DataMatrix.aero_gamma_state.project_V.r*aeroPr*obj.DataMatrix.aero_gamma_state.project_V.r.';
% 
% 
%         end
%         function obj = loadOrCompute(obj,cfg)
%             if isfile(obj.cacheFile)
%                 S = load(obj.cacheFile);
%                 % obj = S.obj;
%                 obj = obj.assign(S.obj); 
%                 return
%             end
%             % -- read aero mesh ------------------------------------------------
%             h5_file = fullfile(obj.dataPath, "pazy_h5");
%             h5_lin = fullfile(obj.dataPath, "pazy_lin");
%             aeroLoc = fullfile(h5_file, obj.aeroFile);
%             obj.aeroMesh = AeroFlex.aero.readAEROFromH5(aeroLoc);
%             % -- build ROM -----------------------------------------------------
%             [ROM_krylov_SS, forces_aero_beam_dof, ~, ROM_krylov_Disc, forces_aero, forces_structure, obj.DataMatrix, etas_sharpy] = AeroFlex.aero.AeroDiscToCont(h5_lin, obj.case_name, obj.NetworkPath);
% 
% 
%             if cfg.aero.discrete
%                 obj.ROM_cts = ROM_krylov_Disc;
%                 % obj.ROM_cts  = obj.DataMatrix.KrylovSharpyContSS;
% 
%                 % obj.ROM_dsc  = c2d(obj.ROM_cts, obj.ROM_cts.Ts,'tustin')
%             else
%                  % obj.ROM_dsc = ROM_krylov_Disc;
%                  % obj.ROM_cts  = d2c(obj.ROM_dsc,'tustin');
%                  obj.ROM_cts  = obj.DataMatrix.KrylovSharpyContSS;
%                 % obj.ROM_cts = ROM_krylov_SS;
% 
%             end
% 
%             obj.forceMap = AeroFlex.aero.buildForceMap(obj.ROM_cts,obj.Na, obj.Nm, cfg.aero.discrete, cfg.ctrl.n_surf, cfg.ctrl.var_per, cfg.gust.on, obj.dataPath);
% 
%             % -- Get gust input ------------------------------------------------
%             if cfg.gust.on
%                 disp('Gust Enabled')
%                 gust = AeroFlex.aero.get_1minuscosine_gust_input(...
%                     cfg.gust.gust_length, ...
%                     cfg.gust.gust_offset, ...
%                     cfg.gust.gust_intensity, ...
%                     cfg.sim.dt, ...
%                     cfg.flight.U_inf, ...
%                     cfg.sim.t_end);
%                 % % obj.gust_input = 10*gust(:,2);
%                 % obj.gust_input = gust(:,2);
%                 % obj.gust_input = cfg.flight.a*gust(:,2);
%                 % obj.gust_input = .5*gust(:,2);
%                 obj.gust_input = .525*gust(:,2);
%                 % obj.gust_input = .4082*gust(:,2);
%                 % gust_input = zeros(simulation_time/dt,1);
%             else
%                 disp('Gust Disabled')
%                 obj.gust_input = zeros(cfg.sim.t_end/cfg.sim.dt,1);
% 
%             end
% 
%             obj.DataMatrix.forces_aero_beam_dof = forces_aero_beam_dof;
%             obj.DataMatrix.forces_aero = forces_aero;
%             obj.DataMatrix.forces_structure = forces_structure;
%             obj.DataMatrix.etas_sharpy = etas_sharpy;
% 
% 
% 
% 
%             % obj.DataMatrix = DataMatrix;
%             % switch nargout('AeroDiscToCont')
%             %     case 1
%             %         obj.ROM_cts = AeroDiscToCont(h5_lin);
%             %         obj.forceMap = buildForceMap(obj.ROM_cts,obj.Na);
%             %         obj.ROM_dsc  = c2d(obj.ROM_cts, obj.ROM_cts.Ts,'tustin');
%             %     otherwise
%             %         [obj.ROM_cts,out2] = AeroDiscToCont(h5_lin);
%             %         if isstruct(out2) && isfield(out2,'D1')
%             %             obj.forceMap = out2;
%             %             obj.ROM_dsc  = c2d(obj.ROM_cts, obj.ROM_cts.Ts,'tustin');
%             %         else
%             %             obj.ROM_dsc = out2;
%             %             obj.forceMap = buildForceMap(obj.ROM_cts,obj.Na);
%             %         end
%             % end
%             save(obj.cacheFile,'obj','-v7.3');
%         end
%         function obj = assign(obj,S)
%             f = fieldnames(S);
%             for k=1:numel(f), if isprop(obj,f{k}), obj.(f{k}) = S.(f{k}); end, end
%         end
%     end
% end
