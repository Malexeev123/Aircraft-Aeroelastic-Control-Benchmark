function [ROM_krylov_SS, forces_aero_beam_dof, mode_shape, ROM_krylov_Disc, forces_aero, forces_structure, DataMatrix, etas] = AeroDiscToCont(aeroObj, cfg)
% AERODISCTOCONT  Build continuous/discrete ROM and collect SHARPy matrices.
%
% INPUTS
%   aeroObj  : instance of AeroFlex.aero.AeroROM (has resolved paths & filenames)
%   cfg      : configuration struct (for options)
%
% OUTPUTS
%   ROM_krylov_SS     : continuous-time ROM (via Tustin d2c)
%   ROM_krylov_Disc   : discrete-time ROM (as exported by SHARPy)
%   forces_*          : linearisation vectors from <case>.linss.h5
%   DataMatrix        : coupling matrices from <case>.data.h5
%                       + projection V (if found)
%                       + KrylovSharpyContSS (continuous via analytic Tustin)
%   etas              : struct with eta, eta_dot from linss
%
% FILES READ (ALL under savedata except aero mesh):
%   - aeroObj.rom_h5        (…/savedata/*.rom.h5)
%   - aeroObj.linss_h5      (…/savedata/<case>.linss.h5)
%   - aeroObj.data_h5       (…/savedata/<case>.data.h5)
%   - aeroObj.vproj_h5 (opt) (…/savedata/*aerorob*.h5)

    debug = (isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level>=1);

    %% 1) Read the discrete ROM from *.rom.h5 (SAVEDATA)
    [A_d, B_d, C_d, D_d, Ts] = h5_read_ssrom(aeroObj.rom_h5);
    ROM_krylov_Disc = ss(A_d, B_d, C_d, D_d, Ts);
    if debug
        fprintf('[AeroDiscToCont] ROM_dsc: n=%d, m=%d, p=%d, Ts=%.6g\n', size(A_d,1), size(B_d,2), size(C_d,1), Ts);
    end

    % dt_nd = ROM_krylov_Disc.Ts;
    % t_ref = cfg.t_map;
    % 
    % dt_phys = dt_nd * t_ref;    
    % % dt_phys = dt_nd /t_ref;    
    % 
    % ROM_krylov_Disc.Ts = dt_phys;

    %% 2) Continuous versions (two ways)
    try
        ROM_krylov_SS = d2c(ROM_krylov_Disc, d2cOptions('Method','tustin'));
    catch
        ROM_krylov_SS = tustin_ct(ROM_krylov_Disc);
    end
    KrylovSharpyContSS = tustin_ct(ROM_krylov_Disc); % keep your validator

    %% 3) Linearisation vectors from <case>.linss.h5
    [forces_aero_beam_dof, mode_shape, forces_aero, forces_structure, etas] = ...
        read_linearisation_vectors(aeroObj.linss_h5);

    %% 4) Coupling matrices from <case>.data.h5
    DataMatrix = read_coupling_matrices(aeroObj.data_h5);

    %% 5) Projection V (optional; savedata)
    if ~isempty(aeroObj.vproj_h5) && isfile(aeroObj.vproj_h5)
        try
            Vgain = h5read(aeroObj.vproj_h5,'/V/gain');
            DataMatrix.aero_gamma_state.project_V = struct('r', Vgain);
        catch ME
            warning('[AeroDiscToCont] Could not read V/gain from %s: %s', aeroObj.vproj_h5, ME.message);
            DataMatrix.aero_gamma_state.project_V = struct('r', []);
        end
    else
        DataMatrix.aero_gamma_state.project_V = struct('r', []);
    end

    %% 6) Stash CT validator into DataMatrix
    DataMatrix.KrylovSharpyContSS = KrylovSharpyContSS;

    %% 7) Optional check
    if debug
        del = ROM_krylov_SS.A - KrylovSharpyContSS.A;
        del(abs(del) < 1e-10) = 0;
        fprintf('[AeroDiscToCont] CT: %d A-entries differ beyond 1e-10.\n', nnz(del));
    end
end


%% ============================ Local helpers ==============================

function [A,B,C,D,Ts] = h5_read_ssrom(fn)
    assert(isfile(fn), 'ROM file not found: %s', fn);
    A = h5read_complex(fn,'/ssrom/A');
    B = h5read_complex(fn,'/ssrom/B');
    C = h5read_complex(fn,'/ssrom/C');
    D = h5read_complex(fn,'/ssrom/D');
    Ts = h5read(fn,'/ssrom/dt'); Ts = double(Ts(1));
    if size(A,1) ~= size(B,1), A = A.'; B = B.'; end
    if size(C,1) ~= size(D,1), C = C.'; D = D.'; end
end

function X = h5read_complex(fn, path)
    X = h5read(fn, path);
    if isstruct(X) && isfield(X,'r')  % SHARPy complex compound
        X = X.r;
    end
end

function sysc = tustin_ct(sysd)
    n = size(sysd.A,1);
    Ts = sysd.Ts;
    omega0 = 2/Ts;
    I = eye(n);
    A = sysd.A; B = sysd.B; C = sysd.C; D = sysd.D;
    M = (I + A);
    syscA = omega0 * (A - I) / M;
    syscB = sqrt(2*omega0) * (M \ B);
    syscC = sqrt(2*omega0) * (C / M);
    syscD = D - C / M * B;
    sysc = ss(syscA, syscB, syscC, syscD);
end

function [f_abd, modes, f_aero, f_struct, etas] = read_linearisation_vectors(linss_h5)
    assert(isfile(linss_h5), 'linss file not found: %s', linss_h5);
    f_abd   = h5read(linss_h5, '/linearisation_vectors/forces_aero_beam_dof');
    modes   = h5read(linss_h5, '/linearisation_vectors/mode_shapes');
    f_aero  = h5read(linss_h5, '/linearisation_vectors/forces_aero');
    f_struct= h5read(linss_h5, '/linearisation_vectors/forces_struct');
    eta     = h5read(linss_h5, '/linearisation_vectors/eta');
    eta_dot = h5read(linss_h5, '/linearisation_vectors/eta_dot');
    etas = struct('eta',eta,'eta_dot',eta_dot);
end

function DM = read_coupling_matrices(data_h5)
    assert(isfile(data_h5), 'data file not found: %s', data_h5);
    DM.Crr       = h5read(data_h5, '/data/linear/linear_system/Crr').';
    DM.Crs       = h5read(data_h5, '/data/linear/linear_system/Crs');
    DM.Csr       = h5read(data_h5, '/data/linear/linear_system/Csr');
    DM.Kdisp     = h5read(data_h5, '/data/linear/linear_system/Kdisp');
    DM.Kdisp_vel = h5read(data_h5, '/data/linear/linear_system/Kdisp_vel');
    DM.Kforces   = h5read(data_h5, '/data/linear/linear_system/Kforces').';
    DM.Krs       = h5read(data_h5, '/data/linear/linear_system/Krs');
    DM.Kss       = h5read(data_h5, '/data/linear/linear_system/Kss');
    DM.Kvel_disp = h5read(data_h5, '/data/linear/linear_system/Kvel_disp');
    DM.Kvel_vel  = h5read(data_h5, '/data/linear/linear_system/Kvel_vel');

    DM.B_to_vertex = h5read(data_h5, '/data/linear/linear_system/uvlm/B_to_vertex_forces');
    DM.C_to_vertex = h5read(data_h5, '/data/linear/linear_system/uvlm/C_to_vertex_forces');
    DM.D_to_vertex = h5read(data_h5, '/data/linear/linear_system/uvlm/D_to_vertex_forces');

    DM.q0start   = h5read(data_h5, '/data/linear/tsstruct0/q');
    DM.q1start   = h5read(data_h5, '/data/linear/tsstruct0/dqdt');

    DM.aero_gamma_state.gamma0     = h5read(data_h5, '/data/linear/tsaero0/gamma/_as_array');
    DM.aero_gamma_state.gamma      = h5read(data_h5, '/data/linear/tsaero0/gamma/_as_array');
    DM.aero_gamma_state.gamma_dot  = h5read(data_h5, '/data/linear/tsaero0/gamma_dot/_as_array');
    DM.aero_gamma_state.gamma_star = h5read(data_h5, '/data/linear/tsaero0/gamma_star/_as_array');
    DM.aero_gamma_state.gamma_star0= h5read(data_h5, '/data/linear/tsaero0/gamma_star/_as_array');
end


% function [ROM_krylov_SS, forces_aero_beam_dof, mode_shape, ROM_krylov_Disc, forces_aero, forces_structure, DataMatrix, etas] = AeroDiscToCont(h5_lin, case_name, NetworkPath)
%     %, KrylovSharpyContSS
%     % Need to change to actual file path, works for now;
%     warning("off");
%     if NetworkPath 
%         route_directory = "\\malexeev\users\maxal\Documents\JupyterNotebooks\ControlsProj\";
%     else
%         route_directory = "/Users/maxal/Documents/JupyterNotebooks/ControlsProj";
%     end
%     [ROM_krylov_SS, ROM_krylov_Disc, KrylovSharpyContSS, forces_aero_beam_dof, mode_shape,DataMatrix, etas] = AeroDiscToContF(route_directory, case_name);
%     delta = ROM_krylov_SS.A-ROM_krylov_Disc.A;
%     contValidation = ROM_krylov_SS.A - KrylovSharpyContSS.A;
%     contValidation = contValidation.*(contValidation>1e-10 | contValidation<-1e-10);
%     if ~NetworkPath
%         save('AeroSS', "ROM_krylov_SS")
%         save('KrylovSharpyContSS', "KrylovSharpyContSS")
%     end
%     %%
%     function [ROM_krylov_SS, state_space_system, KrylovSharpyContSS, forces_aero_beam_dof, mode_shape,DataMatrix, etas] = AeroDiscToContF(route_directory, case_name)
%         addpath(strcat(route_directory,'example_notebooks/UDP_control/matlab_functions/'));
%         success = false;
%         %% Define parameters;
%         % case_name = 'pazy_ROM';
% 
%         % fold_name = 'pazy_ROM_krylov';        
%         fold_name = append(case_name,'_krylov');
%         output_folder = strcat(route_directory,'/output/',case_name);
%         output_folder_linear = strcat(output_folder, '/linear_results/');
%         complex_matrices = true;
% 
%         %% Reade from SHARPy generated state space model
%         krylovFolder = strcat(fold_name, '.rom.h5');
%         [state_space_system, forces_aero_beam_dof, mode_shape, forces_aero, forces_structure, etas] = read_SHARPy_state_space_system2(...
%                                                 krylovFolder, ...
%                                                 strcat(output_folder, '/savedata/'), ...
%                                                 complex_matrices, ...
%                                                 case_name);
%         opts = d2cOptions('Method','tustin');%,'PrewarpFrequency',20
%         ROM_krylov_SS = d2c(state_space_system, opts);
%         KrylovSharpyContSS = tustin(state_space_system);
% 
%         Asd = state_space_system.A;
%         Asd_size= size(Asd);
%         Bsd = state_space_system.B;
%         Csd = size(state_space_system.C);
%         Dsd = size(state_space_system.D);
% 
% 
%         Co = ctrb(Asd,Bsd);
%         % rank(Co)
% 
%         DataMatrix = readCoupleMats(strcat(output_folder, '/savedata/'), case_name);
%         project_V = readProjectionKrylov(strcat(output_folder, '/save_pmor_data/'), case_name);
%         DataMatrix.aero_gamma_state.project_V = project_V;
%         DataMatrix.KrylovSharpyContSS = KrylovSharpyContSS;
%         % %% Load simulation settings and store in struct
%         % load(strcat(output_folder_linear, 'simulation_parameters.mat'));
%         % n_nodes = uint8(n_nodes);
%         % rigid_body_motions = false;
% 
%         % input_settings =  set_input_parameters(u_inf, ...
%         %                                       state_space_system.Ts, ...
%         %                                       num_modes, ...
%         %                                       num_aero_states, ...
%         %                                       rigid_body_motions, ...
%         %                                       simulation_time, ...
%         %                                       num_control_surfaces, ...
%         %                                       n_nodes, ...
%         %                                       control_input_start, ...
%         %                                       gust_input_start);
% 
%         % %% Remove unused input and outputs from the generated ROM in SHARPy and link deflection and its rate    
%         % state_space_system = adjust_state_space_system(state_space_system, input_settings);
%         % sadasd = size(state_space_system.A)
%         % Asd = state_space_system.A;
%         % Bsd = state_space_system.B;
%         % Csd = size(state_space_system.C);
%         % Dsd = size(state_space_system.D);
%         % 
% 
% 
% 
%         success = true;
%         disp(success)
% 
%     end
% 
%     function project_V = readProjectionKrylov(folder, case_name)
%         % fold_nameROM = 'pazy_ROM_krylov_aerorob';
%         fold_nameROM = append(case_name,'_krylov_aerorob');
%         folderROM = strcat(fold_nameROM, '.h5');
%         forces_file_path = strcat(folder, folderROM);
%         project_V = h5read(forces_file_path, '/V/gain');
% 
%     end
%     function DataMatrix = readCoupleMats(folder, case_name)
%         % fold_nameROM = 'pazy_ROM';
%         fold_nameROM = case_name;
%         folderROM = strcat(fold_nameROM, '.data.h5');
%         forces_file_path = strcat(folder, folderROM);
%         Crr = h5read(forces_file_path, '/data/linear/linear_system/Crr');
%         Crr = Crr.';
%         Crs = h5read(forces_file_path, '/data/linear/linear_system/Crs');
%         Csr = h5read(forces_file_path, '/data/linear/linear_system/Csr');
%         Kdisp = h5read(forces_file_path, '/data/linear/linear_system/Kdisp');
%         Kdisp_vel = h5read(forces_file_path, '/data/linear/linear_system/Kdisp_vel');
%         Kforces = h5read(forces_file_path, '/data/linear/linear_system/Kforces');
%         Kforces = Kforces.';
%         Krs = h5read(forces_file_path, '/data/linear/linear_system/Krs');
%         Kss = h5read(forces_file_path, '/data/linear/linear_system/Kss');
%         Kvel_disp = h5read(forces_file_path, '/data/linear/linear_system/Kvel_disp');
%         Kvel_vel = h5read(forces_file_path, '/data/linear/linear_system/Kvel_vel');
% 
%         B_to_vertex = h5read(forces_file_path, '/data/linear/linear_system/uvlm/B_to_vertex_forces');
%         C_to_vertex = h5read(forces_file_path, '/data/linear/linear_system/uvlm/C_to_vertex_forces');
%         D_to_vertex = h5read(forces_file_path, '/data/linear/linear_system/uvlm/D_to_vertex_forces');
% 
%         q0start = h5read(forces_file_path, '/data/linear/tsstruct0/q');
%         q1start = h5read(forces_file_path, '/data/linear/tsstruct0/dqdt');
% 
% 
%         % gamma0 = h5read(forces_file_path, '/data/aero/timestep_info/00000/gamma/_as_array');
%         % gamma = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma/_as_array');
%         % gamma_dot = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma_dot/_as_array');
%         % gamma_star = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma_star/_as_array');
%         % gamma_star0 = h5read(forces_file_path, '/data/aero/timestep_info/00000/gamma_star/_as_array');
% 
%         gamma0 = h5read(forces_file_path, '/data/linear/tsaero0/gamma/_as_array');
%         gamma = h5read(forces_file_path, '/data/linear/tsaero0/gamma/_as_array');
%         gamma_dot = h5read(forces_file_path, '/data/linear/tsaero0/gamma_dot/_as_array');
%         gamma_star = h5read(forces_file_path, '/data/linear/tsaero0/gamma_star/_as_array');
%         gamma_star0 = h5read(forces_file_path, '/data/linear/tsaero0/gamma_star/_as_array');
% 
%         aero_gamma_state = struct('gamma0',gamma0, 'gamma', gamma, 'gamma_dot',gamma_dot, 'gamma_star', gamma_star, 'gamma_star0',gamma_star0);
% 
%         DataMatrix = struct('Crr',Crr,'Crs', Crs,'Csr', Csr,'Kdisp', Kdisp,'Kdisp_vel',Kdisp_vel, ...
%             'Kforces',Kforces,'Krs',Krs,'Kss',Kss,'Kvel_disp',Kvel_disp, ...
%             'Kvel_vel',Kvel_vel, 'q0start', q0start, 'q1start', q1start, 'aero_gamma_state', aero_gamma_state, ...
%             'B_to_vertex', B_to_vertex, 'C_to_vertex', C_to_vertex, 'D_to_vertex', D_to_vertex);
% 
%     end
% 
%     function [state_space_system, forces_aero_beam_dof, mode_shape, forces_aero, forces_structure, etas] = read_SHARPy_state_space_system2(file_name_sharpy, folder, complex, case_name)
%     %read_SHARPy_state_space_system Gets sate space model from SHARPy output
%     absolute_file_path = strcat(folder, file_name_sharpy);
% 
%     fold_nameROM = case_name;
%     folderROM = strcat(fold_nameROM, '.linss.h5');
%     forces_file_path = strcat(folder, folderROM);
%     forces_aero_beam_dof = h5read(forces_file_path, '/linearisation_vectors/forces_aero_beam_dof');
%     mode_shape = h5read(forces_file_path, '/linearisation_vectors/mode_shapes');
%     forces_aero = h5read(forces_file_path, '/linearisation_vectors/forces_aero');
%     forces_structure = h5read(forces_file_path, '/linearisation_vectors/forces_struct');
%     eta = h5read(forces_file_path, '/linearisation_vectors/eta');
%     eta_dot = h5read(forces_file_path, '/linearisation_vectors/eta_dot');
%     etas.eta = eta; etas.eta_dot = eta_dot;
%     % Reference Point and Timestep
%     %eta_ref = h5read(absolute_file_path, '/linearisation_vectors/eta');
%     dt = h5read(absolute_file_path, '/ssrom/dt');
%     % Read state space matrices 
%     if complex
%         A = h5read(absolute_file_path, '/ssrom/A').r; 
%         B = h5read(absolute_file_path, '/ssrom/B').r; 
%         C = h5read(absolute_file_path, '/ssrom/C').r;
%         D = h5read(absolute_file_path, '/ssrom/D');
%     else
%         A = h5read(absolute_file_path, '/ssrom/A'); 
%         B = h5read(absolute_file_path, '/ssrom/B'); 
%         C = h5read(absolute_file_path, '/ssrom/C');
%         D = h5read(absolute_file_path, '/ssrom/D');
%     end
% 
%     % Sometimes matrices are exported transposed
%     if size(A, 1) ~= size(B, 1)
%        A = transpose(A); 
%        B = transpose(B); 
%     end
%     if size(C, 1) ~= size(D, 1)
%        C = transpose(C); 
%        D = transpose(D); 
%     end
% 
%     % Assemble final discrete state space system
%     state_space_system = ss(A, B, C, D, dt);
% 
%     end 
% 
%     function [KrylovSharpyContSS] = tustin(ssDisc)
%         omeg0 = 2/ssDisc.Ts;
%         num_states = length(ssDisc.A);
%         Abar = omeg0*(ssDisc.A-eye(num_states))/(eye(num_states)+ssDisc.A);
% 
%         Bbar = (sqrt(2*omeg0)*inv(eye(num_states)+ssDisc.A))*ssDisc.B;
%         Cbar = sqrt(2*omeg0)*ssDisc.C/(eye(num_states)+ssDisc.A);
%         Dbar = ssDisc.D-ssDisc.C/(eye(num_states)+ssDisc.A)*ssDisc.B;
%         KrylovSharpyContSS = ss(Abar, Bbar, Cbar, Dbar);
% 
%     end
% end
% 
