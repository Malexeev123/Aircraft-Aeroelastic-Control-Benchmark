function [ROM_krylov_SS, forces_aero_beam_dof, mode_shape, ROM_krylov_Disc, forces_aero, forces_structure, DataMatrix, etas] = AeroDiscToCont(h5_lin)
    %, KrylovSharpyContSS
    
    route_directory = "/Users/maxal/Documents/JupyterNotebooks/ControlsProj";
    [ROM_krylov_SS, ROM_krylov_Disc, KrylovSharpyContSS, forces_aero_beam_dof, mode_shape,DataMatrix, etas] = AeroDiscToContF(route_directory);
    delta = ROM_krylov_SS.A-ROM_krylov_Disc.A;
    contValidation = ROM_krylov_SS.A - KrylovSharpyContSS.A;
    contValidation = contValidation.*(contValidation>1e-10 | contValidation<-1e-10);
    save('AeroSS', "ROM_krylov_SS")
    save('KrylovSharpyContSS', "KrylovSharpyContSS")

    %%
    function [ROM_krylov_SS, state_space_system, KrylovSharpyContSS, forces_aero_beam_dof, mode_shape,DataMatrix, etas] = AeroDiscToContF(route_directory)
        addpath(strcat(route_directory,'example_notebooks/UDP_control/matlab_functions/'));
        success = false;
        %% Define parameters;
        case_name = 'pazy_ROM';
        fold_name = 'pazy_ROM_krylov';
        output_folder = strcat(route_directory,'/output/',case_name);
        output_folder_linear = strcat(output_folder, '/linear_results/');
        complex_matrices = true;
        
        %% Reade from SHARPy generated state space model
        krylovFolder = strcat(fold_name, '.rom.h5');
        [state_space_system, forces_aero_beam_dof, mode_shape, forces_aero, forces_structure, etas] = read_SHARPy_state_space_system2(...
                                                krylovFolder, ...
                                                strcat(output_folder, '/savedata/'), ...
                                                complex_matrices...
                                                );
        opts = d2cOptions('Method','tustin');%,'PrewarpFrequency',20
        ROM_krylov_SS = d2c(state_space_system, opts);
        KrylovSharpyContSS = tustin(state_space_system);
        
        Asd = state_space_system.A;
        Asd_size= size(Asd);
        Bsd = state_space_system.B;
        Csd = size(state_space_system.C);
        Dsd = size(state_space_system.D);
    
    
        Co = ctrb(Asd,Bsd);
        % rank(Co)

        DataMatrix = readCoupleMats(strcat(output_folder, '/savedata/'));
        project_V = readProjectionKrylov(strcat(output_folder, '/save_pmor_data/'));
        DataMatrix.aero_gamma_state.project_V = project_V;
        % %% Load simulation settings and store in struct
        % load(strcat(output_folder_linear, 'simulation_parameters.mat'));
        % n_nodes = uint8(n_nodes);
        % rigid_body_motions = false;
        
        % input_settings =  set_input_parameters(u_inf, ...
        %                                       state_space_system.Ts, ...
        %                                       num_modes, ...
        %                                       num_aero_states, ...
        %                                       rigid_body_motions, ...
        %                                       simulation_time, ...
        %                                       num_control_surfaces, ...
        %                                       n_nodes, ...
        %                                       control_input_start, ...
        %                                       gust_input_start);
        
        % %% Remove unused input and outputs from the generated ROM in SHARPy and link deflection and its rate    
        % state_space_system = adjust_state_space_system(state_space_system, input_settings);
        % sadasd = size(state_space_system.A)
        % Asd = state_space_system.A;
        % Bsd = state_space_system.B;
        % Csd = size(state_space_system.C);
        % Dsd = size(state_space_system.D);
        % 
        
       
    
        success = true;
        disp(success)
    
    end

    function project_V = readProjectionKrylov(folder)
        fold_nameROM = 'pazy_ROM_krylov_aerorob';
        folderROM = strcat(fold_nameROM, '.h5');
        forces_file_path = strcat(folder, folderROM);
        project_V = h5read(forces_file_path, '/V/gain');
        
    end
    function DataMatrix = readCoupleMats(folder)
        fold_nameROM = 'pazy_ROM';
        folderROM = strcat(fold_nameROM, '.data.h5');
        forces_file_path = strcat(folder, folderROM);
        Crr = h5read(forces_file_path, '/data/linear/linear_system/Crr');
        Crr = Crr.';
        Crs = h5read(forces_file_path, '/data/linear/linear_system/Crs');
        Csr = h5read(forces_file_path, '/data/linear/linear_system/Csr');
        Kdisp = h5read(forces_file_path, '/data/linear/linear_system/Kdisp');
        Kdisp_vel = h5read(forces_file_path, '/data/linear/linear_system/Kdisp_vel');
        Kforces = h5read(forces_file_path, '/data/linear/linear_system/Kforces');
        Kforces = Kforces.';
        Krs = h5read(forces_file_path, '/data/linear/linear_system/Krs');
        Kss = h5read(forces_file_path, '/data/linear/linear_system/Kss');
        Kvel_disp = h5read(forces_file_path, '/data/linear/linear_system/Kvel_disp');
        Kvel_vel = h5read(forces_file_path, '/data/linear/linear_system/Kvel_vel');
        
        q0start = h5read(forces_file_path, '/data/linear/tsstruct0/q');
        q1start = h5read(forces_file_path, '/data/linear/tsstruct0/dqdt');


        % gamma0 = h5read(forces_file_path, '/data/aero/timestep_info/00000/gamma/_as_array');
        % gamma = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma/_as_array');
        % gamma_dot = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma_dot/_as_array');
        % gamma_star = h5read(forces_file_path, '/data/aero/timestep_info/00001/gamma_star/_as_array');
        % gamma_star0 = h5read(forces_file_path, '/data/aero/timestep_info/00000/gamma_star/_as_array');

        gamma0 = h5read(forces_file_path, '/data/linear/tsaero0/gamma/_as_array');
        gamma = h5read(forces_file_path, '/data/linear/tsaero0/gamma/_as_array');
        gamma_dot = h5read(forces_file_path, '/data/linear/tsaero0/gamma_dot/_as_array');
        gamma_star = h5read(forces_file_path, '/data/linear/tsaero0/gamma_star/_as_array');
        gamma_star0 = h5read(forces_file_path, '/data/linear/tsaero0/gamma_star/_as_array');
        
        aero_gamma_state = struct('gamma0',gamma0, 'gamma', gamma, 'gamma_dot',gamma_dot, 'gamma_star', gamma_star, 'gamma_star0',gamma_star0);

        DataMatrix = struct('Crr',Crr,'Crs', Crs,'Csr', Csr,'Kdisp', Kdisp,'Kdisp_vel',Kdisp_vel, ...
            'Kforces',Kforces,'Krs',Krs,'Kss',Kss,'Kvel_disp',Kvel_disp, ...
            'Kvel_vel',Kvel_vel, 'q0start', q0start, 'q1start', q1start, 'aero_gamma_state', aero_gamma_state);
        
    end

    function [state_space_system, forces_aero_beam_dof, mode_shape, forces_aero, forces_structure, etas] = read_SHARPy_state_space_system2(file_name_sharpy, folder, complex)
    %read_SHARPy_state_space_system Gets sate space model from SHARPy output
    absolute_file_path = strcat(folder, file_name_sharpy);

    fold_nameROM = 'pazy_ROM';
    folderROM = strcat(fold_nameROM, '.linss.h5');
    forces_file_path = strcat(folder, folderROM);
    forces_aero_beam_dof = h5read(forces_file_path, '/linearisation_vectors/forces_aero_beam_dof');
    mode_shape = h5read(forces_file_path, '/linearisation_vectors/mode_shapes');
    forces_aero = h5read(forces_file_path, '/linearisation_vectors/forces_aero');
    forces_structure = h5read(forces_file_path, '/linearisation_vectors/forces_struct');
    eta = h5read(forces_file_path, '/linearisation_vectors/eta');
    eta_dot = h5read(forces_file_path, '/linearisation_vectors/eta_dot');
    etas.eta = eta; etas.eta_dot = eta_dot;
    % Reference Point and Timestep
    %eta_ref = h5read(absolute_file_path, '/linearisation_vectors/eta');
    dt = h5read(absolute_file_path, '/ssrom/dt');
    % Read state space matrices 
    if complex
        A = h5read(absolute_file_path, '/ssrom/A').r; 
        B = h5read(absolute_file_path, '/ssrom/B').r; 
        C = h5read(absolute_file_path, '/ssrom/C').r;
        D = h5read(absolute_file_path, '/ssrom/D');
    else
        A = h5read(absolute_file_path, '/ssrom/A'); 
        B = h5read(absolute_file_path, '/ssrom/B'); 
        C = h5read(absolute_file_path, '/ssrom/C');
        D = h5read(absolute_file_path, '/ssrom/D');
    end
    
    % Sometimes matrices are exported transposed
    if size(A, 1) ~= size(B, 1)
       A = transpose(A); 
       B = transpose(B); 
    end
    if size(C, 1) ~= size(D, 1)
       C = transpose(C); 
       D = transpose(D); 
    end
    
    % Assemble final discrete state space system
    state_space_system = ss(A, B, C, D, dt);
    
    end 

    function [KrylovSharpyContSS] = tustin(ssDisc)
        omeg0 = 2/ssDisc.Ts;
        num_states = length(ssDisc.A);
        Abar = omeg0*(ssDisc.A-eye(num_states))/(eye(num_states)+ssDisc.A);

        Bbar = (sqrt(2*omeg0)*inv(eye(num_states)+ssDisc.A))*ssDisc.B;
        Cbar = sqrt(2*omeg0)*ssDisc.C/(eye(num_states)+ssDisc.A);
        Dbar = ssDisc.D-ssDisc.C/(eye(num_states)+ssDisc.A)*ssDisc.B;
        KrylovSharpyContSS = ss(Abar, Bbar, Cbar, Dbar);
    
    end
end

