function base = build_xi_coupling(beam, aero, cfg)
%BUILD_XI_COUPLING  Construct quaternion strain basis and ROM couplings (no cache).
%
% Inputs
%   beam : BeamModel   (needs: fem, aero, phi0, phi1, Mglobal, red.phi1_sA)
%   aero : AeroROM     (needs: DataMatrix with Csr, K*, B_to_vertex, aero_gamma_state.project_V.r)
%   cfg  : struct      (needs: Nm, ctrl.n_surf, ctrl.var_per, flight.aoa_deg; optional: struct.n_rigid_nodes)
%
% Outputs (struct)
%   base.phi_xi_modes : (Nq×Nm)  quaternion/strain modal matrix
%   base.Gamma_xi     : cell{Nm}  quadratic tensors from beam Γ
%   base.Gamma_g      : (Nm×1)    gravity projection
%   base.xi_bar       : (Nq×1)    static quaternion strain offset
%   base.FM           : struct with χ-coupling maps
%       • Dchi : Nm×3      = φ0ᵀ * Csrᵀ(:,end-2:end)
%       • Bchi : Na×3      = V * (B_to_vertex_trim * Kas(:,end-2:end))
%       • X    : 3×3       = orientation at clamped node (from φ1 slice & ξ̄)
%   base.phiXi_sA     : 4×Nm      quaternion slice at the clamped node
%
% Notes
%   • No SHARPy file I/O here — everything comes from beam/aero objects.
%   • No caching — caller (sim_init) is responsible for saving the bundle.


    % -------- guard: required fields --------
    requireField(beam,'phi1','beam.phi1 missing (run redtest in BeamModel first).');
    requireField(beam,'phi0','beam.phi0 missing.');
    requireField(beam,'Mglobal','beam.Mglobal missing.');
    requireField(beam,'red','beam.red missing.');
    requireField(beam,'red.phi1_sA','beam.red.phi1_sA missing (from redtest).');

    requireField(aero,'DataMatrix','aero.DataMatrix missing.');
    requireField(aero.DataMatrix,'Csr','aero.DataMatrix.Csr missing.');
    requireField(aero.DataMatrix,'Kdisp','aero.DataMatrix.Kdisp missing.');
    requireField(aero.DataMatrix,'Kdisp_vel','aero.DataMatrix.Kdisp_vel missing.');
    requireField(aero.DataMatrix,'Kvel_disp','aero.DataMatrix.Kvel_disp missing.');
    requireField(aero.DataMatrix,'Kvel_vel','aero.DataMatrix.Kvel_vel missing.');
    requireField(aero.DataMatrix,'B_to_vertex','aero.DataMatrix.B_to_vertex missing.');
    requireField(aero.DataMatrix,'aero_gamma_state','aero.DataMatrix.aero_gamma_state missing.');
    requireField(aero.DataMatrix.aero_gamma_state,'project_V','aero.DataMatrix.aero_gamma_state.project_V missing.');
    requireField(aero.DataMatrix.aero_gamma_state.project_V,'r','aero.DataMatrix.aero_gamma_state.project_V.r missing.');

    requireField(cfg,'Nm','cfg.Nm missing.');
    requireField(cfg,'ctrl','cfg.ctrl missing.');
    requireField(cfg.ctrl,'n_surf','cfg.ctrl.n_surf missing.');
    requireField(cfg.ctrl,'var_per','cfg.ctrl.var_per missing.');
    requireField(cfg,'flight','cfg.flight missing.');
    requireField(cfg.flight,'aoa_deg','cfg.flight.aoa_deg missing.');

    % -------- 1) Baseline quaternion ξ pack (uses your existing solver) --------
    [phi_xi_modes, Gamma_xi_mat, Gamma_g_3D, ~, xi_bar_array] = ...
        AeroFlex.core.solve_baseline_quaternion_xi( ...
            beam.fem, beam.aero, beam.phi1, beam.Nm, beam.Mglobal, cfg.flight.aoa_deg, beam.phi0);

    base.phi_xi_modes = phi_xi_modes;
    base.Gamma_xi     = Gamma_xi_mat;
    base.Gamma_g      = Gamma_g_3D;
    base.xi_bar       = xi_bar_array;

    % base.Gamma_xi = zeros(size(base.Gamma_xi));
    % base.Gamma_g = zeros(size(base.Gamma_g));
    % -------- 2) χ-force mappings (FM) ----------------------------------------
    FM = struct();

    % (a) Dchi = φ0ᵀ * Csrᵀ(:, end-2:end)  — 
    CsrT    = aero.DataMatrix.Csr.';                     
    FM.Dchi = beam.phi0.' * CsrT(:, end-2:end);          % Nm×3
    % FM.Dchi = cfg.force_map.scale*beam.phi0.' * CsrT(:, end-2:end);          % Nm×3
    
    % (b) Assemble Kas from SHARPy’s K blocks
    if isfield(cfg,'struct') && isfield(cfg.struct,'n_rigid_nodes') && ~isempty(cfg.struct.n_rigid_nodes)
        num_rigid_node = cfg.struct.n_rigid_nodes;
    else
        num_rigid_node = 9;  % legacy default in your script
    end

    Kdisp     = aero.DataMatrix.Kdisp(1:end-num_rigid_node,:).';
    Kdisp_vel = aero.DataMatrix.Kdisp_vel.';     % full
    Kvel_disp = aero.DataMatrix.Kvel_disp(1:end-num_rigid_node,:).';
    Kvel_vel  = aero.DataMatrix.Kvel_vel.';      % full

    Kas = [Kdisp, Kdisp_vel; Kvel_disp, Kvel_vel];      % (Γ+χ)×(Γ+χ)

    % (c) Select χ columns ( last 3 columns)
    Kas_xi = Kas(:, end-2:end);                          % (Γ+χ)×3

    % (d) Project B_to_vertex to remove CS+gust channels, then to Krylov Γ subspace
    u_ext_idx   = cfg.ctrl.n_surf * cfg.ctrl.var_per + 1;   % +1 gust
    B_project   = aero.DataMatrix.B_to_vertex(1:end-u_ext_idx,:).';  % Γ×(#panels)
    Bchi_unproj = B_project * Kas_xi;                        % Γ×3

    V_project = aero.DataMatrix.aero_gamma_state.project_V.r; % Na×Γ
    FM.Bchi   = V_project.r * Bchi_unproj;                     % Na×3

    % (e) Orientation at clamped node
    [FM.X, FM.chi_bar ]= AeroFlex.core.compute_orientation_matrix(beam.red.phi1_sA, base.xi_bar);
    
    base.FM = FM;

    % -------- 3) Quaternion slice at the clamped node -------------------------
    coords  = beam.fem.coordinates;  % assume y-axis is index 2
    sL = coords(2,2); sR = coords(end,2); sA = coords(1,2);
    elemLen  = (sR - sL);

    phixiL = base.phi_xi_modes(1:4, :);               % left quat rows
    phixiR = base.phi_xi_modes(end-3:end, :);         % right quat rows
    base.phiXi_sA = phixiL * (sR - sA)/elemLen + phixiR * (sA - sL)/elemLen;  % 4×Nm

    % -------- 4) Optional debug prints ---------------------------------------
    if dbg(cfg)
        fprintf('[build_xi_coupling] Nm=%d | size(φξ)=%dx%d | Dχ=%dx%d | Bχ=%dx%d\n', ...
            cfg.Nm, size(base.phi_xi_modes,1), size(base.phi_xi_modes,2), ...
            size(base.FM.Dchi,1), size(base.FM.Dchi,2), ...
            size(base.FM.Bchi,1), size(base.FM.Bchi,2));
    end
end

% ------------------------ helpers ------------------------
% ------------------------ helpers ------------------------
function requireField(obj, path, msg)
    if ~hasPath(obj, path)
        error('%s', msg);
    end
    % Optional: also enforce not-empty (comment out if you only care about existence)
    val = getPath(obj, path);
    if isempty(val)
        error('%s is empty.', path);
    end
end

function tf = hasPath(obj, path)
    toks = strsplit(path, '.');
    cur = obj;
    for k = 1:numel(toks)
        t = toks{k};
        if isstruct(cur)
            if ~isfield(cur, t), tf = false; return; end
        else
            % object (handle or value)
            if ~(isprop(cur, t) || isfield(struct(cur), t))
                tf = false; return;
            end
        end
        cur = cur.(t);
    end
    tf = true;
end

function val = getPath(obj, path)
    toks = strsplit(path, '.');
    cur = obj;
    for k = 1:numel(toks)
        cur = cur.(toks{k});
    end
    val = cur;
end

function tf = dbg(cfg)
    tf = isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level>=1;
end

% function base = build_xi_coupling(beam, aero, cfg)
% %BASELINEMODES  Build modal strain basis ξ and nonlinear tensors.
% %
% % Inputs
% %   beam : BeamModel   (has phi0, Gamma1, Gamma2, eta_e)
% %   aero : AeroROM     (only used for panel z, y in Γg ⟵ gravity term)
% %   cfg  : config      (aoa, rho, U_inf in case needed later)
% %
% % Outputs (packed in struct)
% %   phi_xi   : (Nm+1) × Nm    augmented displacement/strain modes
% %   Gamma_xi : cell{Nm}       quadratic tensors ψᵀΓψ from beam.Gamma*
% %   Gamma_g  : Nm × 1         gravity projection
% %   xi_bar   : (Nm+1) × 1     static strain offset
% if cfg.Nm ==20
%     cacheFile = fullfile(cfg.paths.cache,'base.mat');
% else
%     cacheFile = fullfile(cfg.paths.cache, cfg.case_name,'base.mat');
% 
% end
% % cacheFileAero = fullfile(cfg.paths.cache,'rom.mat');
% if isfile(cacheFile)
%                 S = load(cacheFile);
%                 % obj = S.obj;
%                 base = S.base; 
% 
% else
% 
%     [phi_xi_modes, Gamma_xi_mat, Gamma_g_3D, phi_omega, xi_bar_array] =AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, beam.aero , ...
%         beam.phi1, beam.Nm, beam.Mglobal, cfg.flight.aoa_deg);
% 
%     base.phi_xi_modes = phi_xi_modes;
%     base.Gamma_xi = Gamma_xi_mat;
%     base.Gamma_g = Gamma_g_3D;
%     base.xi_bar = xi_bar_array;
%     save(cacheFile,'base','-v7.3');
% end
% Csr = aero.DataMatrix.Csr.';
% FM.Dchi = beam.phi0.'*Csr(:,end-2:end);
% % Dchi = FM.Dchi;
% % Build Kas 
% num_rigid_node = 9;
% 
% Kdisp = aero.DataMatrix.Kdisp(1:end-num_rigid_node,:).';
% Kdisp_vel = aero.DataMatrix.Kdisp_vel.';
% Kvel_disp = aero.DataMatrix.Kvel_disp(1:end-num_rigid_node,:).';
% Kvel_vel = aero.DataMatrix.Kvel_vel.';
% Kas = [Kdisp, Kdisp_vel;
%     Kvel_disp, Kvel_vel];
% Kas_xi = Kas(:, end-2:end);
% % if cfg.gust.on
%     u_ext_idx = cfg.ctrl.n_surf*cfg.ctrl.var_per+1;
% % else
%     % u_ext_idx = cfg.ctrl.n_surf*cfg.ctrl.var_per+1;
% 
% % end
% B_project =  aero.DataMatrix.B_to_vertex(1:end-u_ext_idx,:).';
% Bchi_unproj = B_project*Kas_xi;
% V_project = aero.DataMatrix.aero_gamma_state.project_V.r;
% FM.Bchi = V_project*Bchi_unproj;
% % 
% % FM.Bchi = zeros(size(FM.Bchi));
% % FM.Dchi = zeros(size(FM.Dchi));
% 
% 
% % Bchi = FM.Bchi;
% base.FM = FM;
% 
% 
% X_mat = AeroFlex.core.compute_orientation_matrix(beam.red.phi1_sA, base.xi_bar);
% base.FM.X = X_mat;
% 
% 
% 
% sL = beam.fem.coordinates(2,2);
% sR = beam.fem.coordinates(end,2);
% sA = beam.fem.coordinates(1,2);
% elemLen = (sR-sL);
% phixiL = base.phi_xi_modes(1:4, :); % first 4 quat
% phixiR = base.phi_xi_modes(end-3:end, :); % last 4 quat
% phiXi_sA = (phixiL*(sR - sA)/elemLen+ phixiR*(sA-sL)/elemLen);
% base.phiXi_sA = phiXi_sA;
% save(cacheFile,'-append','FM');
% % Kas_dot =  [Kvel_disp, Kvel_vel];
% 
% 
% % Nm = beam.Nm;
% % 
% % % --- φξ  (augment flexible modes with a “1” row for chord‑wise strain) ---
% % phi_xi         = [beam.phi0 ; zeros(1,Nm)];
% % phi_xi(end, :) = 1;           % constant axial component (Artola Eq. 34)
% % 
% % % --- Γξ  (stack Γ1 & Γ2 tensors along 3rd dim for fast access) ----------
% % Gamma_xi = cell(1,Nm);
% % for k = 1:Nm
% %     Gamma_xi{k} = beam.Gamma1{k} + beam.Gamma2{k};  % sample combination
% % end
% % 
% % % --- Γg  (gravity) ------------------------------------------------------
% % rho = cfg.flight.rho;  g = 9.80665;
% % dz  = aero.aeroMesh.z;   % vector Na×1 of panel z‑coords in local frame
% % % Simple strip formula: ∫ φ0 ρg z ds
% % Gamma_g = rho*g * (beam.phi0.' * dz);
% % 
% % % --- ξ̄  (static offset) -------------------------------------------------
% % xi_bar = zeros(Nm+1,1);     % zero because trimmed reference is undeformed
% % 
% % base = struct('phi_xi',phi_xi,...
% %               'Gamma_xi',{Gamma_xi},...
% %               'Gamma_g',Gamma_g,...
% %               'xi_bar',xi_bar);
% %   Uses spanwise strip theory: Δα = (z p − y r) / U∞
% %   Inputs
% %       beam : BeamModel   (phi0 only)
% %       aero : AeroROM     (aeroMesh with bound panel coords)
% %       cfg  : config      (U_inf, b_ref)
% %
% %   Outputs
% %       Dchi : Nm × 3
% %       Bchi : Na × 3
% 
% % U   = cfg.flight.U_inf;
% % bref = cfg.flight.b_ref;
% % rho  = cfg.flight.rho;          % for magnitude scaling
% % % CLα  = 2*pi;                    % thin‑airfoil slope
% % 
% % % --- geometry vectors ---------------------------------------------------
% % Z = aero.aeroMesh.z(:);         % Na × 1 panel z
% % Y = aero.aeroMesh.y(:);         % Na × 1 panel y
% % X = aero.aeroMesh.x(:);         % Na × 1 panel x  (for q down‑wash)
% % 
% % Na = length(Z);  Nm = beam.Nm;
% % Bchi = zeros(Na,3);
% % Bchi(:,1) =  Z / bref;          % ∂Γ̇/∂p
% % Bchi(:,3) = -Y / bref;          % ∂Γ̇/∂r
% % % q‑column is zero for lift‑line; add X‑term if pitch‑rate wash wanted
% % 
% % % --- Dchi via modal projection of strip lift ----------------------------
% % phi0 = beam.phi0_flex;          % 96 × Nm (translation DOFs only)
% % dy   = diff(Y); dy(end+1)=dy(end); dy = abs(dy);   % span strip width
% % c    = aero.aeroMesh.c(:);      % local chord
% % % LiftPerStrip = rho*U^2 .* c .* CLα;   % scalar×vector Na×1
% % 
% % Dchi = zeros(Nm,3);
% % for i = 1:Nm
% %     Dchi(i,1) =  sum( LiftPerStrip .* Z .* phi0(:,i) .* dy ) / bref;
% %     Dchi(i,3) = -sum( LiftPerStrip .* Y .* phi0(:,i) .* dy ) / bref;
% %     % Dchi(:,2) pitch‑rate term negligible for straight wing
% % end
% end
