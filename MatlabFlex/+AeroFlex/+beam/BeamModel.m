classdef BeamModel
    % BEAMMODEL  SHARPy-aware structural model and modal extractor.
    %
    % Reads:
    %   - cases/<case>/<case>.fem.h5                (FEM mesh & BCs)
    %   - cases/<case>/<case>.aero.h5               (for consistency with legacy redtest)
    %   - output/<case>/beam_modal_analysis/        (Kglobal.dat, Mglobal.dat, *_rigid.dat)  [preferred]
    %   - output/<case>/savedata/<case>.data.h5     (fallback: auto-discover K/M datasets)
    %
    % Produces:
    %   fem, aero, Kglobal, Mglobal, modal set (phi0, Omega, etc.),
    %   red tensors via runRedtestWrapper (phi1, phi2, Gamma1, Gamma2, eta_e),
    %   projectors Pz/Pr for phi1 and phi2, Sigma damping vector (zeros unless cfg.struct.damping=true).
    %
    % Required cfg fields (set in sim_init):
    %   cfg.Nm
    %   cfg.case_name
    %   cfg.paths.sharpy.root
    %   cfg.paths.sharpy.cases     = <root>/cases/<case>
    %   cfg.paths.sharpy.savedata  = <root>/output/<case>/savedata
    %   cfg.paths.sharpy.beam_modal= <root>/output/<case>/beam_modal_analysis   (optional; will be inferred if absent)
    %
    % Notes:
    %   * No cache is required. If you want to keep a cache, you can re-enable obj.cacheFile usage.
    %   * Robust HDF5 traversal is used to locate K/M in <case>.data.h5 if the DATs are missing.

    properties
        % sizes & ids
        Nm
        case_name

        % Mesh & aero support
        femFile
        aeroFile
        fem
        aero

        % DOF bookkeeping
        nDofs
        nRig
        nFlex
        keepDofs

        % Global structural matrices
        Kglobal
        Mglobal

        % Modal data
        phi0        % Nm mass-normalised displacement modes (flex rows if retain_rigid=false)
        phi1
        phi2
        Omega       % diag(omega_i)
        Khat
        Mhat
        eigVals
        red 

        % “rigid” six-DOF slice used in projectors (as produced by redtest)
        phi_sA
        Pr, Pz
        Pr_phi2, Pz_phi2
        Zva, Zvr

        % tensors / offsets
        Gamma1
        Gamma2
        eta_e

        % options
        Sigma       % damping vector
        cfg

        % SHARPy paths (resolved on construction)
        sharpy_root
        cases_dir
        savedata_dir
        beam_modal_dir

        % fallback data file
        data_h5
    end

    methods
        function obj = BeamModel(cfg)
            % save("cfg", cfg, '-mat')
            % ---- sizes/ids ----
            obj.Nm        = cfg.Nm;
            obj.case_name = cfg.case_name;
            obj.cfg       = cfg;

            % ---- resolve SHARPy folders ----
            [obj.sharpy_root, obj.cases_dir, obj.savedata_dir, obj.beam_modal_dir] = resolve_sharpy_beam_roots(cfg);

            % expected filenames in cases/
            obj.femFile  = fullfile(obj.cases_dir,  [obj.case_name, '.fem.h5']);
            obj.aeroFile = fullfile(obj.cases_dir,  [obj.case_name, '.aero.h5']);

            % fallback data.h5 for K/M search
            obj.data_h5  = fullfile(obj.savedata_dir, [obj.case_name, '.data.h5']);

            % ---- 1) load mesh (and aero for redtest consistency) ----
            obj = obj.loadMesh();

            % ---- 2) load/import modal set (K/M -> eig -> mass-normalise -> red tensors) ----
            obj = obj.buildModalSet(cfg);
        end

        %% ----------------------- mesh loading -------------------------
        function obj = loadMesh(obj)
            assert(isfile(obj.femFile),  'FEM file not found: %s', obj.femFile);
            assert(isfile(obj.aeroFile), 'AERO file not found: %s', obj.aeroFile);

            obj.fem  = AeroFlex.beam.readFEMFromH5(obj.femFile);
            obj.aero = AeroFlex.aero.readAEROFromH5(obj.aeroFile);

            % DOF bookkeeping from boundary conditions
            dofPerNode = 6;
            allDofs    = 1:(obj.fem.num_node * dofPerNode);

            bcs = obj.fem.boundary_conditions;  % 1=clamped, etc.
            clampedNodes = find(bcs == 1);
            removeDofs = [];
            for cN = clampedNodes(:).'
                these = (cN-1)*6 + (1:6);
                removeDofs = [removeDofs, these]; %#ok<AGROW>
            end
            obj.keepDofs = setdiff(allDofs, removeDofs);
            obj.fem.keepDofs = obj.keepDofs;

            obj.nDofs = numel(allDofs);
            obj.nFlex = numel(obj.keepDofs);
            obj.nRig  = obj.nDofs - obj.nFlex;
        end

        %% ------------------- modal import / build ----------------------
        function obj = buildModalSet(obj, cfg)
            % 2.A) try beam_modal_analysis DATs first
            [K, M, sourceTag] = try_load_KM_from_beam_modal(obj.beam_modal_dir, cfg.struct.rigid_body_modes);
   
            % 2.B) fallback: auto-discover K/M datasets in <case>.data.h5
            if isempty(K) || isempty(M)
                assert(isfile(obj.data_h5), 'Fallback data file not found: %s', obj.data_h5);
                [K, M, sourceTag] = try_load_KM_from_data_h5(obj.data_h5, cfg.struct.rigid_body_modes);
            end

            obj.Kglobal = K;
            obj.Mglobal = M;

            if isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level>=1
                fprintf('[BeamModel] Using K/M from: %s\n', sourceTag);
            end

            % 2.C) solve generalised eigenproblem K v = lambda M v
            % [V, Lam] = eig((K+K')/2, (M+M')/2);  % enforce symmetry numerically
            % writematrix(V,'VFull.dat','Delimiter','tab');

            [V, Lam] = eig(K,M);  
            % writematrix(Vt,'VFullUnsymm.dat','Delimiter','tab');

            % Ks = 0.5*(K+K'); Ms = 0.5*(M+M');
            % symK = norm(K-K.','fro')/max(1,norm(K,'fro'));
            % symM = norm(M-M.','fro')/max(1,norm(M,'fro'));
            % if symK < 1e-10 && symM < 1e-10
            %     [V,Lam] = eig(K,M);
            % else
            %     warning('[BeamModel] Using symmetrized (K,M) due to symmetry defects.');
            %     [V,Lam] = eig(Ks,Ms);
            % end


            [w_all, ix] = sort(sqrt(max(eps, real(diag(Lam)))));
            V = V(:, ix);

            % pick Nm
            w  = w_all(1:obj.Nm);
            Vr = V(:,   1:obj.Nm);

            % 2.D) mass-normalise φ0 (column-wise)
            Mred = Vr.' * M * Vr;
            Mred = 0.5*(Mred + Mred.');  % numerical symmetry
            
            
            % obj.Kglobal = K;
            % disp('size(obj.Mglobal)')           
            % disp(size(obj.Mglobal))
            % 
            % obj.Mglobal = Mred;
            % disp('size(obj.Mglobal)')
            % disp(size(obj.Mglobal))

            % diagnostics
            symK = norm(K-K.','fro')/max(1,norm(K,'fro'));
            symM = norm(M-M.','fro')/max(1,norm(M,'fro'));
            minEigM = min(eig(0.5*(M+M.')));
            if isfield(cfg,'debug') && isfield(cfg.debug,'level') && cfg.debug.level>=1
                fprintf('[BeamModel] sym(K)=%.3e, sym(M)=%.3e, minEig(Msym)=%.3e\n', symK, symM, minEigM);
            end
            
            % robust M-orthonormalization
            [Rchol,p] = chol(Mred,'upper');
            if p ~= 0
                warning('[BeamModel] Mred not SPD in chol; falling back to diagonal scaling.');
                scale = 1 ./ sqrt(max(eps, diag(Mred)));
                Vmass = Vr .* scale.';
            else
                Vmass = Vr / Rchol;  % Vr * inv(Rchol) => Vmass' M Vmass = I
            end
            % Vmass = Vr;
            % writematrix(Vmass,'V_Curr.dat','Delimiter','tab');
            % writematrix(-Vmass,'V_CurrNeg.dat','Delimiter','tab');
            % deterministic sign convention to avoid "negative tip" confusion:
            % make the largest-magnitude component in each mode positive.
            % negCount = 0;
            % for k = 1:size(Vmass,2)
            %     [~,imax] = max(abs(Vmass(:,k)));
            %     if Vmass(imax,k) < 0
            %         % disp('to neg')
            %          negCount = negCount+1;
            %          Vmass(:,k) = -Vmass(:,k);
            %         % Vmass(:,k) = -Vmass(:,k);
            %     end
            % end
            % if negCount > size(Vmass,2)/2
            %     Vmass= -Vmass;
            % end


            

            % [~,idMax] = max(abs(Vmass))
            % if Vmass(idMax) < 0
            %     disp('Neg')
            %     Vmass = -Vmass;
            % end


            % scale = 1 ./ sqrt(diag(Mred));
            % Vmass = Vr .* scale.';      % mass-normalised
            % % Vmass = Vr;
            % writematrix(Vmass,'V_Curr.dat','Delimiter','tab');
            % writematrix(-Vmass,'V_CurrNeg.dat','Delimiter','tab');
            % Vmass = -Vmass;
            % disp('Vmass norm'); disp(norm(Vmass));


            obj.Omega = diag(w);
            % obj.Mhat  = eye(obj.Nm);
            % obj.Khat  = diag(w.^2);
            obj.eigVals = w.^2;

            % 2.E) drop rigid rows unless asked to retain
            if ~cfg.struct.retain_rigid
                % To do, re-implement Guyan reduction, just use SHARPy
                % input for now
                % Vkeep = Vmass(obj.keepDofs, :);
                Vkeep = Vmass;

            else
                Vkeep = Vmass;
            end
            % Temp
            % numN = size(Vkeep,1);
            % halfN = numN/2+1;
            % Vkeep(halfN,:) = -Vkeep(halfN,:);

            obj.phi0 = Vkeep;
       
            % 2.F) Reduced-order geometric tensors & higher derivatives
            if isempty(obj.Gamma1) || isempty(obj.Gamma2) || isempty(obj.phi1) || isempty(obj.phi2)
                red = AeroFlex.beam.runRedtestWrapper(obj.fem, obj.aero, obj.phi0, obj.Nm, cfg, obj);
                obj.phi0   = red.ModeVars_discrete.phi0_local;
                obj.phi1   = red.ModeVars_discrete.phi1_local;
                obj.phi2   = red.ModeVars_discrete.phi2_local;
                obj.Gamma1 = red.ModeVars_discrete.Gamma1;
                obj.Gamma2 = red.ModeVars_discrete.Gamma2;
                obj.eta_e  = red.ModeVars_discrete.etaVecExt;
                obj.phi_sA = red.phi1_sA;   % clamped-node rigid slice that redtest exposes
                obj.red = red;
            end
            
            % 2.G) Damping (placeholder)
            if isfield(cfg.struct,'damping') && cfg.struct.damping
                obj.Sigma = red.ModeVars_discrete.Sigma;
            else
                obj.Sigma = zeros(obj.Nm,1);
            end

            % 2.H) Projectors (phi1-based)
            obj.Zva = null(obj.phi_sA);
            obj.Pz  = obj.Zva/(obj.Zva.'*obj.Zva) * obj.Zva.';
            obj.Zvr = orth(obj.phi_sA.');
            obj.Pr  = obj.Zvr/(obj.Zvr.'*obj.Zvr) * obj.Zvr.';

            % 2.I) Projectors (phi2-based)
            Zva2    = null(red.phi2_sA);
            obj.Pz_phi2 = Zva2/(Zva2.'*Zva2) * Zva2.';
            Zvr2    = orth(red.phi2_sA.');
            obj.Pr_phi2 = Zvr2/(Zvr2.'*Zvr2) * Zvr2.';
        end

        %% ----------------------- utility ------------------------------
        function obj = assign(obj,S)
            f = fieldnames(S);
            for k = 1:numel(f)
                if isprop(obj,f{k}) && ~isempty(S.(f{k}))
                    obj.(f{k}) = S.(f{k});
                end
            end
        end
    end
end

%% ====================== Local helpers (file I/O) =========================

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

function [K,M,tag] = try_load_KM_from_beam_modal(beam_modal_dir, use_rigid)
    K = []; M = []; tag = '';
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
        K = load(Kfile);
        M = load(Mfile);
        tag = sprintf('beam_modal_analysis (%s, %s)', get_last(Kfile), get_last(Mfile));
    end
end

function [K,M,tag] = try_load_KM_from_data_h5(data_h5, use_rigid)
    info = h5info(data_h5);
    % search likely names first
    candK = {'/data/structural/Kglobal','/data/structural/K_global','/Kglobal','/K_global'};
    candM = {'/data/structural/Mglobal','/data/structural/M_global','/Mglobal','/M_global'};

    [K, kpath] = try_read_first(data_h5, candK);
    [M, mpath] = try_read_first(data_h5, candM);

    if isempty(K) || isempty(M)
        % brute-force discovery: first square float dataset whose name contains 'K' / 'M'
        if isempty(K), [K,kpath] = discover_square(data_h5, info, 'K'); end
        if isempty(M), [M,mpath] = discover_square(data_h5, info, 'M'); end
    end

    if isempty(K) || isempty(M)
        error('Could not locate K/M in %s. Please provide DATs under beam_modal_analysis or ensure datasets exist.', data_h5);
    end

    tag = sprintf('data.h5 (%s , %s)', kpath, mpath);

    % if "rigid" requested and there exist *_rigid datasets, prefer them
    if use_rigid
        [Kr, krp] = try_read_first(data_h5, {'/data/structural/Kglobal_rigid','/data/structural/K_rigid'});
        [Mr, mrp] = try_read_first(data_h5, {'/data/structural/Mglobal_rigid','/data/structural/M_rigid'});
        if ~isempty(Kr) && ~isempty(Mr)
            K = Kr; M = Mr; tag = sprintf('%s [rigid override: %s , %s]', tag, krp, mrp);
        end
    end
end

function [X, path] = try_read_first(h5fn, paths)
    X = []; path = '';
    for i=1:numel(paths)
        p = paths{i};
        if has_dataset(h5fn, p)
            X = h5read(h5fn, p);
            path = p;
            return;
        end
    end
end

function tf = has_dataset(h5fn, path)
    tf = false;
    try
        info = h5info(h5fn, path); %#ok<NASGU>
        tf = true;
    catch
        tf = false;
    end
end

function [X, path] = discover_square(h5fn, info, keyChar)
    % Depth-first search for a square, >= 6x6, floating dataset whose name contains keyChar (e.g., 'K' or 'M')
    X = []; path = '';
    stack = info.Groups;
    % also inspect datasets at root
    for d = 1:numel(info.Datasets)
        ds = info.Datasets(d);
        if contains(lower(ds.Name), lower(keyChar))
            sz = fliplr(ds.Dataspace.Size);
            if numel(sz)==2 && sz(1)==sz(2) && sz(1)>=6
                trial = h5read(h5fn, ['/' ds.Name]);
                if isfloat(trial)
                    X = trial; path = ['/' ds.Name]; return;
                end
            end
        end
    end
    % recurse groups
    groups = info.Groups;
    while ~isempty(groups)
        g = groups(1); groups(1) = [];
        % datasets in group
        for d=1:numel(g.Datasets)
            ds = g.Datasets(d);
            if contains(lower(ds.Name), lower(keyChar))
                fullp = [g.Name '/' ds.Name];
                sz = fliplr(ds.Dataspace.Size);
                if numel(sz)==2 && sz(1)==sz(2) && sz(1)>=6
                    trial = h5read(h5fn, fullp);
                    if isfloat(trial)
                        X = trial; path = fullp; return;
                    end
                end
            end
        end
        % enqueue subgroups
        groups = [groups; g.Groups]; %#ok<AGROW>
    end
end

function s = get_last(p)
    [~,s,ext] = fileparts(p); s = [s ext];
end
%%
function [Vcanon, meta] = canonicalise_internal_modes(obj, Vc, M, opts)
% CANONICALISE_INTERNAL_MODES
% Make internal (flexible) modes consistent:
%   • expand weighting Wx/Wy/Wz to ndofs×ndofs
%   • compute per-mode dominant global axis (x/y/z) via Gram energies
%   • flip each mode so tip translation along that dominant axis is positive
%   • (re)mass-normalise columns defensively
%
% Inputs
%   obj  : BeamModel instance (must have obj.fem.coordinates, obj.keepDofs)
%   Vc   : ndofs×Nm modal matrix (displacement-type, i.e. φ0)
%   M    : ndofs×ndofs reduced mass matrix (same space as Vc)
%   opts : (optional) struct with fields:
%          .tip_side   ∈ {'+','-','auto'}  default 'auto'
%          .weights    struct with fields Wx,Wy,Wz (6-vec | 6×6 | nd×1 | nd×nd)
%          .reorth     logical, default true (re-mass-normalise)
%
% Outputs
%   Vcanon : canonicalised modes (same size as Vc)
%   meta   : struct with fields:
%            .span_axis, .axis_order, .tip_node, .energies (Nm×3),
%            .flips (1×Nm), .diag_after (Nm×1)

    arguments
        obj
        Vc double
        M  double
        opts.tip_side char {mustBeMember(opts.tip_side,{'auto','+','-'})} = 'auto'
        opts.weights struct = struct()
        opts.reorth logical = true
    end

    % --- sanity ---
    nd = size(M,1);
    assert(size(Vc,1)==nd, 'Vc and M dimension mismatch.');
    Nm = size(Vc,2);
    dpn = 6;  % DOFs per node
    assert(mod(nd,dpn)==0, 'ndofs must be multiple of 6.');
    nActN = nd/dpn;

    % --- geometry & span axis (no "range" dependency) ---
    fem = obj.fem;
    XYZ = fem.coordinates;
    if size(XYZ,2) ~= 3, XYZ = XYZ.'; end
    assert(size(XYZ,2)==3, 'fem.coordinates must be [num_node×3] or [3×num_node].');
    if size(XYZ,1) ~= (numel(obj.keepDofs)/dpn)
        % condensed DOFs: coordinates may include omitted nodes; that’s fine
        % we only need span detection & tip selection
    end

    mins   = min(XYZ,[],1);
    maxs   = max(XYZ,[],1);
    ranges = maxs - mins;     % 1×3
    [~, span_axis] = max(ranges);  % 1=X,2=Y,3=Z
    % order the axes by decreasing range (span, chord, thickness)
    [~, order] = sort(ranges,'descend');  % e.g., [2 1 3] for Pazy (Y biggest)
    meta.span_axis = span_axis;
    meta.axis_order = order;

    % --- default per-node weights for x/y/z translations only (ux,uy,uz,rx,ry,rz) ---
    % Allow overrides via opts.weights (6-vector, 6×6, nd×1 or nd×nd).
    defaultWx = [1 0 0 0 0 0].';
    defaultWy = [0 1 0 0 0 0].';
    defaultWz = [0 0 1 0 0 0].';
    Wx = pick_weight(opts, 'Wx', defaultWx);
    Wy = pick_weight(opts, 'Wy', defaultWy);
    Wz = pick_weight(opts, 'Wz', defaultWz);

    % --- expand Wx/Wy/Wz to nd×nd diagonal (or pass-through if already nd×nd) ---
    Wx = make_W_nd(Wx, nd, nActN, dpn);
    Wy = make_W_nd(Wy, nd, nActN, dpn);
    Wz = make_W_nd(Wz, nd, nActN, dpn);

    % --- symmetricise M (defensive) ---
    M  = 0.5*(M+M.');

    % --- mass re-normalise columns (defensive; preserves span) ---
    H  = Vc.'*(M*Vc);
    dH = real(diag(H));
    dH(dH<=eps) = 1;              % avoid divide by 0
    Vc = Vc ./ sqrt(dH.');        % column-wise scaling

    % --- compute per-axis Gram energies for each mode (diagonals only) ---
    Gx = Vc.' * (M * (Wx * Vc));  Gx = 0.5*(Gx+Gx.');
    Gy = Vc.' * (M * (Wy * Vc));  Gy = 0.5*(Gy+Gy.');
    Gz = Vc.' * (M * (Wz * Vc));  Gz = 0.5*(Gz+Gz.');
    ex = real(diag(Gx));  ey = real(diag(Gy));  ez = real(diag(Gz));
    E  = [ex ey ez];              % Nm×3
    [~, domAxis] = max(E,[],2);   % 1=x, 2=y, 3=z per mode
    meta.energies = E;

    % --- choose a tip node robustly (honours keepDofs) ---
    [tipNode, tipSign] = pick_tip_node(fem, span_axis, opts.tip_side);
    meta.tip_node = tipNode;
    meta.tip_side = tipSign;

    % --- ensure we can access the tip’s translational DOFs in the condensed set ---
    keep = obj.keepDofs(:).';
    tip_dofs_full = (tipNode-1)*dpn + (1:3); % ux,uy,uz at tip
    tip_loc = localDofIndex(keep, tip_dofs_full);
    if numel(tip_loc) < 1
        % walk inward until we find a node that survived condensation
        tipNodeBak = tipNode;
        step = sign(tipSign); if step==0, step = 1; end
        nodeList = nearest_nodes_along_span(XYZ, tipNode, span_axis, step);
        found = false;
        for k = 1:numel(nodeList)
            cand = nodeList(k);
            cand_dofs = (cand-1)*dpn + (1:3);
            cand_loc  = localDofIndex(keep, cand_dofs);
            if ~isempty(cand_loc)
                tipNode = cand; tip_loc = cand_loc; found = true; break;
            end
        end
        if ~found
            error('canonicalise_internal_modes:tipDOFsNotKept', ...
                'Could not find a near-tip node present in keepDofs.');
        end
        meta.tip_node_fallback = tipNode;
    end

    % --- flip modes so tip translation along dominant axis is positive ---
    flips = ones(1,Nm);
    for j = 1:Nm
        axis_j = domAxis(j);  % 1=x,2=y,3=z
        % tip_loc are local indices for ux,uy,uz (subset of [1..nd])
        if numel(tip_loc) < axis_j
            % if rotation-only kept at this node (unlikely for Pazy), skip
            continue;
        end
        tip_tr = Vc(tip_loc(axis_j), j);  % scalar
        if tip_tr < 0
            Vc(:,j) = -Vc(:,j);
            flips(j) = -1;
        end
    end

    % --- optional re-orth / re-normalise after flips ---
    if opts.reorth
        H  = Vc.'*(M*Vc);
        % enforce near-orthonormal by symmetric whitening (cheap, stable)
        [Vh, Dh] = eig(0.5*(H+H.'));
        lam = real(diag(Dh));
        lam(lam<=eps) = 1;
        T = Vh * diag(1./sqrt(lam)) * Vh.';   % Nm×Nm
        Vc = Vc * T;
        H2 = Vc.'*(M*Vc);
        meta.diag_after = real(diag(H2));
    else
        meta.diag_after = real(diag(Vc.'*(M*Vc)));
    end

    Vcanon = Vc;
    meta.flips = flips;
end

% ============================== helpers ===============================

function W = pick_weight(opts, field, default6)
    if isfield(opts,'weights') && isfield(opts.weights, field)
        W = opts.weights.(field);
    else
        W = default6;
    end
end

function Wnd = make_W_nd(Waxis, nd, nActN, dpn)
% Expand a per-node weight (6-vec or 6×6) into a nd×nd diagonal matrix.
% If already nd×1 or nd×nd, pass through appropriately.
    if isvector(Waxis) && numel(Waxis)==dpn
        Wdof = repmat(Waxis(:), nActN, 1);
        Wnd  = spdiags(Wdof, 0, nd, nd);
    elseif isequal(size(Waxis), [dpn, dpn])
        Wdof = repmat(diag(Waxis), nActN, 1);
        Wnd  = spdiags(Wdof, 0, nd, nd);
    elseif isvector(Waxis) && numel(Waxis)==nd
        Wnd  = spdiags(Waxis(:), 0, nd, nd);
    elseif isequal(size(Waxis), [nd, nd])
        Wnd  = sparse(Waxis);
    else
        error('make_W_nd:badSize', ...
             'Unsupported weight size %s for nd=%d.', mat2str(size(Waxis)), nd);
    end
end

function [tipNode, tipSign] = pick_tip_node(fem, span_axis, side)
% Choose a tip node on +span or -span side (or auto = farther one).
    XYZ = fem.coordinates; if size(XYZ,2)~=3, XYZ = XYZ.'; end
    [~, iPlus]  = max(XYZ(:,span_axis));  % farthest +span
    [~, iMinus] = min(XYZ(:,span_axis));  % farthest -span
    switch side
        case '+', tipNode = iPlus;  tipSign = +1;
        case '-', tipNode = iMinus; tipSign = -1;
        otherwise
            % pick farther one from the root along span
            if abs(XYZ(iPlus,span_axis)) >= abs(XYZ(iMinus,span_axis))
                tipNode = iPlus;  tipSign = +1;
            else
                tipNode = iMinus; tipSign = -1;
            end
    end
    % If boundary_conditions available, avoid clamped node
    if isfield(fem,'boundary_conditions')
        if fem.boundary_conditions(tipNode)==1
            % if we accidentally chose a clamp, switch side
            if tipSign>0, tipNode = iMinus; tipSign = -1;
            else,         tipNode = iPlus;  tipSign = +1;
            end
        end
    end
end

function idxLocal = localDofIndex(keepDofs, dofSet)
% Map a set of full-space DOF ids to local condensed indices.
    idxLocal = [];
    for dd = dofSet(:).'
        k = find(keepDofs==dd, 1, 'first');
        if ~isempty(k), idxLocal(end+1) = k; end %#ok<AGROW>
    end
end

function nodeList = nearest_nodes_along_span(XYZ, tipNode, span_axis, step)
% Return nodes ordered moving from tipNode inward along span direction.
    y = XYZ(:,span_axis);
    [~, orderIdx] = sort(y,'ascend');
    if step>0
        % from most negative to most positive; start near tipNode and go backwards
        [~, pos] = ismember(tipNode, orderIdx);
        nodeList = fliplr(orderIdx(1:pos));   % toward root
    else
        [~, pos] = ismember(tipNode, orderIdx);
        nodeList = orderIdx(pos:end);         % toward root in the other direction
    end
end


function clusters = cluster_by_gap(w, relTol)
% Group indices of w where consecutive gaps are <= relTol
    if nargin<2, relTol = 0.02; end
    w = w(:);
    Nm = numel(w);
    clusters = {};
    i = 1;
    while i <= Nm
        j = i;
        while j < Nm
            gap = abs(w(j+1)-w(j))/max(w(j),1e-12);
            if gap <= relTol, j = j+1; else, break; end
        end
        clusters{end+1} = i:j; %#ok<AGROW>
        i = j+1;
    end
end

function W = make_W_axis_weighted(fem, axisSel, spanAxis, pow)
% Build a sparse diagonal selector on TRANSLATIONAL DOFs of `axisSel`
% and weight by (|spanCoord|/max)^pow to emphasise tips.
% axisSel, spanAxis ∈ {1(x),2(y),3(z)}
    ndof = fem.num_node*6;
    W = spalloc(ndof, ndof, fem.num_node); % only translational diag entries
    XYZ = fem.coordinates; if size(XYZ,2)~=3, XYZ=XYZ.'; end
    span = XYZ(:, spanAxis);
    smax = max(abs(span)); if smax<1e-12, smax=1; end
    weights = (abs(span)/smax).^pow;   % [nNode x 1]

    for n = 1:fem.num_node
        base = (n-1)*6;
        d = base + axisSel;  % 1..3 translation along axisSel
        W(d,d) = weights(n);
    end
end

function V = enforce_mode_signs(V, M, fem, tipNode, outAxis, Wz, Wy, spanIdx)
% For each mode, choose a meaningful DOF at the right tip and make it positive.
% If a mode has larger out-of-plane energy -> use outAxis at tip,
% else use in-plane axis automatically.

    % decide in-plane axis from span & out-of-plane
    allAxes = [1 2 3];
    inAxis = setdiff(allAxes, [spanIdx, outAxis]);
    if isempty(inAxis), inAxis = 1; end
    inAxis = inAxis(1);

    ndof = size(V,1);
    dofPerNode = 6;
    tip_out = (tipNode-1)*dofPerNode + outAxis;
    tip_in  = (tipNode-1)*dofPerNode + inAxis;

    Nm = size(V,2);
    for j = 1:Nm
        vj = V(:,j);
        Ez = vj.' * (M * (Wz * vj));
        Ey = vj.' * (M * (Wy * vj));

        if Ez >= Ey
            s = sign(V(tip_out, j));
            if s==0, s = sign(sum(abs(V(:,j)).*(diag(Wz)~=0))); end
        else
            s = sign(V(tip_in, j));
            if s==0, s = sign(sum(abs(V(:,j)).*(diag(Wy)~=0))); end
        end
        if s==0, s = 1; end
        V(:,j) = V(:,j) * s;
    end
end


% classdef BeamModel
%     %BEAMMODEL  Structural mesh reader and (later) modal builder.
%     %
%     % Public usage pattern:
%     %   cfg  = nominalConfig;
%     %   beam = AeroFlex.beam.BeamModel(cfg);
%     %   nodes = beam.fem.coordinates;   % access mesh
%     %
%     % Expected fields in cfg:
%     %   cfg.Nm           – number of modes requested
%     %   cfg.paths.cache  – folder for heavy cache *.mat files
%     % ------------------------------------------------------------------
% 
%     properties %(SetAccess = private)
%         Nm, case_name                              % # modes requested
%         % ── mesh
%         femFile                         % .fem.h5 filename
%         fem                             % struct with mesh data
%         aeroFile
%         aero
%         nDofs
%         nRig
%         nFlex
%         % ── cache / data paths
%         cacheFile
%         dataPath
%         % ── modal data
%         phi0    % [Nstruct×Nm] mass‑normalised displacement modes
%         phi1
%         phi2
%         phi_sA  % Clamped Node Rigid Modes
%         Omega   % [Nm×Nm] diag(ω_i)
%         Khat    % [Nm×Nm] diag(ω_i²)
%         Mhat    % [Nm×Nm] identity after mass‑norm
%         Mglobal
%         Kglobal
%         eigVals
%         Gamma1  % cell tensors (optional)
%         Gamma2  % cell tensors (optional)
%         eta_e   % strain offset (optional)
%         cfg
%         Sigma
%         red
%         Pz, Pr
%         Zva, Zvr
%         Pr_phi2, Pz_phi2
%     end
% 
%     methods
%         function obj = BeamModel(cfg)
%             % Constructor: store filenames and load/calc mesh
%             obj.Nm       = cfg.Nm;
%             % Temporary if statement because M = 20 case is default
%             if obj.Nm ==20
%                 obj.femFile  = "pazy_ROM.fem.h5";
%                 obj.aeroFile = "pazy_ROM.aero.h5";
%                 obj.cacheFile= fullfile(cfg.paths.cache,"modal.mat");
%                 obj.dataPath  = cfg.paths.data;
%             else
%                 cfg.case_name = append('pazy_ROM_',string(cfg.Nm),'M');
%                 obj.case_name = cfg.case_name;
%                 obj.femFile = append(obj.case_name,'.fem.h5');
%                 obj.aeroFile =append(obj.case_name,'.aero.h5');
% 
%                 if not(isfolder(fullfile(cfg.paths.cache,obj.case_name)))
%                     mkdir(fullfile(cfg.paths.cache,obj.case_name))
%                 end            
%                 if not(isfolder(fullfile(cfg.paths.data,obj.case_name)))
%                     mkdir(fullfile(cfg.paths.data,obj.case_name))
%                 end 
%                 obj.cacheFile= fullfile(cfg.paths.cache,obj.case_name,"modal.mat");
%                 obj.dataPath  = fullfile(cfg.paths.data,obj.case_name);
% 
%             end
% 
%             obj.cfg = cfg;
% 
%             obj = obj.loadMesh();
%             obj = obj.loadOrImportModal(cfg);
% 
%         end
% 
%         %% ----------------------- mesh loading -------------------------
%         function obj = loadMesh(obj)
% 
%             h5_file = fullfile(obj.dataPath, "pazy_h5");
%             aeroLoc = fullfile(h5_file, obj.aeroFile);
%             femLoc = fullfile(h5_file, obj.femFile);
% 
% 
%             % If cache exists, load from it; otherwise read H5 and cache.
% 
%             if exist(obj.cacheFile,"file")
%                 data = load(obj.cacheFile,"fem", "aero");
%                 obj.fem = data.fem;
%                 obj.aero = data.aero;
%             elseif exist(aeroLoc, "file") && exist(femLoc, "file")
%                 data.aero = AeroFlex.aero.readAEROFromH5(aeroLoc);
%                 data.fem = AeroFlex.beam.readFEMFromH5(femLoc);
%                 obj.fem = data.fem; %fem =obj.fem;
%                 obj.aero = data.aero; %aero = obj.aero;
%                 % save(obj.cacheFile,"fem", "aero", "-v7.3");
%             else
%                 obj.fem = AeroFlex.beam.readFEMFromH5(obj.femFile);
%                 obj.aero = AeroFlex.aero.readAEROFromH5(obj.aeroFile);
%                 % Ensure cache directory exists (robust for fresh clones)
%                 cacheDir = fileparts(obj.cacheFile);
%                 if ~exist(cacheDir, "dir"), mkdir(cacheDir); end
%                 % Save ONLY the fem struct for now (step 3‑A)
%                 % fem = obj.fem;
%                 % aero = obj.aero;
%                 % save(obj.cacheFile,"fem", "aero", "-v7.3");
%             end
%             dofPerNode = 6;
%             allDofs    = 1:(obj.fem.num_node*dofPerNode);
% 
%             % uniqueBCs = unique(fem.boundary_conditions);
%             % fprintf('Unique boundary condition values: ');
%             % disp(uniqueBCs);
%             % 
%             bcs = obj.fem.boundary_conditions;  
%             clampedNodes = find(bcs==1);
%             removeDofs = [];
%             for cN = clampedNodes(:)'
%                 these = (cN-1)*6 + (1:6);
%                 removeDofs = [removeDofs, these];
%             end
%             keepDofs = setdiff(allDofs, removeDofs);
%             obj.fem.keepDofs = keepDofs;
%             fem = obj.fem;
%             aero = obj.aero;
%             save(obj.cacheFile,"fem", "aero", "-v7.3");
%         end
% 
%         %% ---------- Step 3‑B: modal import or build ---------------------
%         function obj = loadOrImportModal(obj, cfg)
%            % 0) cache hit with modal fields
%             if isfile(obj.cacheFile)
%                 data = load(obj.cacheFile);
%                 if isfield(data,"phi0") && isfield(data,"Omega")
%                     obj = obj.assign(data); return; end
%             end
% 
%             % 1) Try fast path: K_global.dat / M_global.dat -------------
%             if cfg.struct.rigid_body_modes 
%                 Kfile = fullfile(obj.dataPath,'Kglobal_rigid.dat');
%                 Mfile = fullfile(obj.dataPath,'Mglobal_rigid.dat');
%                 Mrigid = load(Mfile);
%                 obj.nDofs = size(Mrigid,1);
%                 obj.nFlex = round(round((obj.nDofs - 4)/6) -1)*6;
%                 obj.nRig = obj.nDofs - obj.nFlex;
%             else
%                 Kfile = fullfile(obj.dataPath,'Kglobal.dat');
%                 Mfile = fullfile(obj.dataPath,'Mglobal.dat');
%                 Mflex = load(Mfile);
%                 obj.nDofs = size(Mflex,1);
%                 obj.nFlex = obj.nDofs;
%                 obj.nRig = 0;
%             end
%             Mflex = load(fullfile(obj.dataPath,'Mglobal.dat'));
%             rigidFile = fullfile(cfg.paths.cache,'phi_sA.mat');   % 96×6
%             % if isfile(rigidFile)
%             %     tmp = load(rigidFile);          % variable name: phi_rigid
%             %     obj.phiRigid = tmp.phi_rigid;   
%             % else
%             %     obj.phiRigid = buildRigidModes(obj.fem);          % analytical fallback
%             %     save(rigidFile,'phi_rigid','-v7');
%             % end
%             % 
% 
% 
%             if isfile(Kfile) && isfile(Mfile)
%                 K = load(Kfile);  M = load(Mfile);
%                 sizeM = size(M,2);
%                 [V,lam] = eig(K,M);             % full eig
%                 [lam_sort,~]  = sort(diag(lam));
%                 lam_sort = lam_sort(1:obj.Nm);
% 
%                 [w,ix]  = sort(sqrt(diag(lam)));
%                 % V_unsort = V;
%                 % V_full = V(:,ix);
%                 % nQuat = obj.nRig -6;
%                 % V_Nm = V_full(:, nQuat+1:nQuat+obj.Nm);
%                 % % V_Nm = V_full(:, obj.nRig+1:obj.nRig+obj.Nm);
%                 % % V_Nm = V_full(:, 1:obj.Nm);
%                 % % V_Nm = V_full;
%                 % 
%                 % if isfile(rigidFile)
%                 %     tmp = load(rigidFile);          % variable name: phi_rigid
%                 %     obj.phi_sA = tmp.phi_sA;   
%                 % else
%                 %     obj.phi_sA = V_Nm(obj.nFlex +1: obj.nFlex + 6, :);     
%                 %     phi_sA = obj.phi_sA;
%                 %     save(rigidFile,'phi_sA','-v7');
%                 % end
%                 % 
%                 % rern = null(obj.phi_sA);
%                 % rern_r = null(obj.phi_sA, 'r');
%                 % 
%                 % rrthA = orth(obj.phi_sA);
%                 % rA_rref = rref(obj.phi_sA);
%                 % Zr = rA_rref.';
%                 % 
%                 % % Pz = rern/(rern.'*rern)*rern.';
%                 % Pz = rern_r/(rern_r.'*rern_r)*rern_r.';
%                 % Pr = Zr/(Zr.'*Zr)*Zr.';
%                 % 
%                 % keep96   = 1:obj.nFlex;      % or 1:96
%                 % Vr_red   = obj.rbmProject(V_full, keep96, Mflex);
%                 % phiRigid = Vr_red;
%                 % 
%                 % 
%                 % 
%                 % if cfg.struct.rigid_body_modes 
%                 %     Vrigid = V_full(:, 1:obj.nRig);
%                 % else
%                 %     Vrigif = [];
%                 % end
%                 % Vflex = V_full( :, obj.nRig+1 :end);
% 
% 
%                 V       = V(:,ix(1:obj.Nm));    % pick Nm lowest
%                 Vr = V;
% 
%                 w       = w (1:obj.Nm);
%                 % mass normalise φ₀ -----------------------------------
%                 Mred = V.'*M*V;
%                 scale = 1./sqrt(diag(Mred));
%                 V = V.*scale.';
%                 % lam_norm
%                 % for i=1:obj.Nm
%                 %    sc = sqrt(Vr(:,i)'*(M*Vr(:,i)));
%                 %    Vr(:,i) = Vr(:,i)/sc;
%                 %    lam_sort(i) = lam_sort(i)/(sc^2);
%                 % end
%                 % omegaRad = sqrt(lam_sort);
% 
%                 % split or keep
%                 if ~cfg.struct.retain_rigid
%                     V = V(1:obj.nFlex,:);          % 96×Nm
%                 else
%                     V = Vr;          % 96×Nm
%                     % obj.phi1_flex = [obj.phi1_local , obj.phiRigid];  % 96×(Nm+6)
%                 end
% 
%                 obj.phi0  = V;
%                 % obj.phi0  = Vr;
%                 obj.Omega = diag(w);
%                 obj.Mhat  = eye(obj.Nm);
%                 obj.Khat  = diag(w.^2);
%                 obj.eigVals = lam_sort;
%                 obj.Kglobal = K;
%                 obj.Mglobal = M;
%             else
%                 error('No Kfile or Mfile found');
%             end
% 
%             if isempty(obj.Gamma1) || isempty(obj.Gamma2)
%                 red = AeroFlex.beam.runRedtestWrapper(obj.fem, obj.aero, obj.phi0, obj.Nm, obj.cfg, obj);     % <- pass same cfg
%                 obj.phi0   = red.ModeVars_discrete.phi0_local;
%                 obj.phi1   = red.ModeVars_discrete.phi1_local;
%                 obj.phi2   = red.ModeVars_discrete.phi2_local;
%                 obj.Gamma1 = red.ModeVars_discrete.Gamma1;
%                 obj.Gamma2 = red.ModeVars_discrete.Gamma2;
%                 obj.eta_e  = red.ModeVars_discrete.etaVecExt;
%                 obj.red = red;
%                 % obj.Omega  = diag(red.Omega_vec);
%                 % obj.Khat   = red.ModeVars_discrete.Khat;
%                 % obj.Mhat   = red.ModeVars_discrete.Mhat;
%             end
%             phi0 = obj.phi0; phi1 = obj.phi1; phi2 = obj.phi2;
%             Gamma1 = obj.Gamma1; Gamma2 = obj.Gamma2; eta_e = obj.eta_e;
%             if cfg.struct.damping 
%                 % need to implement KV damping or newmark
%                 % Then modal project
%             else
%                 obj.Sigma = zeros(obj.Nm,1);
%                 Sigma = obj.Sigma;
%             end
%             %% Phi1 Projectors
%             obj.Zva = null(obj.red.phi1_sA);
%             % obj.Zva = null(obj.red.phi1_sA, 'r');
%             % obj.Zva = null(obj.red.phi1_sA, 'rational');
%             obj.Pz = obj.Zva/(obj.Zva.'*obj.Zva)*obj.Zva.';
%             obj.Zvr = orth(obj.red.phi1_sA.');
%             obj.Pr = obj.Zvr/(obj.Zvr.'*obj.Zvr)*obj.Zvr.';
% 
% 
%             %% Phi2 Projectors
%             Zva_phi2= null(obj.red.phi2_sA);
%             % obj.Zva = null(obj.red.phi1_sA, 'r');
%             obj.Pz_phi2 = Zva_phi2/(Zva_phi2.'*Zva_phi2)*Zva_phi2.';
%             Zvr_phi2 = orth(obj.red.phi2_sA.');
%             obj.Pr_phi2 = Zvr_phi2/(Zvr_phi2.'*Zvr_phi2)*Zvr_phi2.';
%             % 3) save full modal set back to cache ---------------------
%             save(obj.cacheFile,'-append', ...
%                  'phi0','phi1','phi2','Gamma1','Gamma2','eta_e', 'Sigma');
%         end
% 
%         %% Utility: copy fields if present -----------------------------
%         function obj = assign(obj,S)
%             f = fieldnames(S);
%             for k = 1:numel(f)
%                 if isprop(obj,f{k}) && ~isempty(S.(f{k}))
%                     obj.(f{k}) = S.(f{k}); end
%             end
%         end
% 
%         function Vr_red = rbmProject(obj, Vr_full, Ired, M96)
%             % Vr_full : 106×6 eigenvectors (full model)
%             % Ired    : 96‑row index or logical mask
%             % M96     : 96×96 truncated mass matrix
%             %
%             % returns 96×6 mass‑normal rigid‑body modes
%             Vr_full = Vr_full(:,obj.nRig - 6+1:obj.nRig);
%             Vr_raw = Vr_full(Ired,:);          % row purge
%             % mass‑orthonormalise -------------------------------------------------
%             Mr = Vr_raw.' * M96 * Vr_raw;      % 6×6
%             [U,S] = eig((Mr+Mr.')/2);          % enforce symmetry
%             Vr_red = Vr_raw / U * diag(1./sqrt(diag(S)));
% 
%             % optional sign convention -------------------------------------------
%             for k = 1:6
%                 col = Vr_red(:,k);
%                 idx = find(abs(col)>1e-8,1,'first');
%                 if col(idx) < 0, Vr_red(:,k) = -col; end
%             end
%         end
% 
%     end
% end
