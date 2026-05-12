function out = postProcess(t, X, cfg, beam, aero, base, log)
% POSTPROCESS  –  decode SimRunner state history, return physical signals
%   Minimal edits: same plots as before, but now saved under a 'plots' folder
%   with improved publication styling (fonts/line widths/white background).

arguments
    t
    X
    cfg
    beam
    aero
    base
    log
end
%======================================================================
% Coupled-state pre-processing
%======================================================================
% Legacy flexible ROM state length:
%   xFlex = [q1(Nm); q2(Nm); qxi(Nm+1); qGamma(Na); chi(3)]
%
% Therefore:
%   nxFlex = 3*Nm + Na + 4
%
% If X includes appended rigid-body states, split them here and continue
% all legacy flexible-wing reconstruction using Xflex only.

Nm = beam.Nm;

if isfield(cfg,'Na') && ~isempty(cfg.Na)
    Na2 = cfg.Na;
elseif isfield(aero,'Na') && ~isempty(aero.Na)
    Na2 = aero.Na;
else
    % Robust inference if X is flexible-only:
    %   size(X,1) = 3*Nm + Na + 4
    Na2 = size(X,1) - (3*Nm + 4);

    % If X includes appended rigid state, this inference is contaminated.
    % In that case, cfg.Na or aero.Na should be provided.
    if Na2 < 0
        error('postProcess:NaInference', ...
              'Cannot infer Na from X. Provide cfg.Na or aero.Na.');
    end
end

nxFlex = 3*Nm + Na2 + 4;

if size(X,1) < nxFlex
    error('postProcess:StateDimension', ...
          'X has %d rows, but expected at least nxFlex = %d.', ...
          size(X,1), nxFlex);
end

% Preserve raw input.
out_raw_X = X;

% Split flexible and rigid parts if present.
Xrb_from_X = [];

if size(X,1) > nxFlex
    Xrb_from_X = X(nxFlex+1:end,:);
    X = X(1:nxFlex,:);
end

% If log.rb exists, prefer it. Otherwise use appended rigid part from X.
rb = struct();
hasRB = false;

if isfield(log,'rb') && isstruct(log.rb) && isfield(log.rb,'state') && ~isempty(log.rb.state)
    rb.state = log.rb.state;
    hasRB = true;
elseif ~isempty(Xrb_from_X)
    rb.state = Xrb_from_X;
    hasRB = true;
end

if hasRB
    % Ensure rb.state is 12 x Nt_rb.
    if size(rb.state,1) ~= 12 && size(rb.state,2) == 12
        rb.state = rb.state.';
    end

    if size(rb.state,1) ~= 12
        warning('postProcess:RigidState', ...
                'Rigid-body state should be 12 x Nt. Ignoring rb state.');
        hasRB = false;
        rb = struct();
    else
        rb.r_I     = rb.state(1:3,:);
        rb.v_B     = rb.state(4:6,:);
        rb.euler   = rb.state(7:9,:);
        rb.omega_B = rb.state(10:12,:);
    end
end

% Reconcile t and X lengths.
Nt = size(X,2);

if numel(t) ~= Nt
    warning('postProcess:TimeLength', ...
            'numel(t)=%d but size(X,2)=%d. Truncating to common length.', ...
            numel(t), Nt);

    Ncommon = min(numel(t), Nt);
    t = t(1:Ncommon);
    X = X(:,1:Ncommon);
    Nt = Ncommon;

    if hasRB
        rb.state = rb.state(:,1:min(size(rb.state,2),Ncommon));
        rb.r_I = rb.r_I(:,1:min(size(rb.r_I,2),Ncommon));
        rb.v_B = rb.v_B(:,1:min(size(rb.v_B,2),Ncommon));
        rb.euler = rb.euler(:,1:min(size(rb.euler,2),Ncommon));
        rb.omega_B = rb.omega_B(:,1:min(size(rb.omega_B,2),Ncommon));
    end
end

% ------------------------ publish defaults -------------------------------
% (higher than default readability; safe for papers/presentations)
set(0,'DefaultAxesFontSize',16, ...
      'DefaultLegendFontSize',15, ...
      'DefaultLineLineWidth',2.2, ...
      'DefaultAxesLineWidth',1.4);

% Where to save plots (prefer the run's plots folder if present)
plots_dir = local_get_plots_dir(cfg);
if ~exist(plots_dir, 'dir'); mkdir(plots_dir); end
% small helper to export figures in PNG+PDF+FIG with 300 dpi
save_plot = @(fh, base) local_export_fig(fh, fullfile(plots_dir, base));

% ------------------------ estimator sampling period ----------------------
if isfield(cfg,'case') && strcmpi(cfg.case,'openLoop')
    log.xhat = zeros(size(X));
    Ts = cfg.sim.dt;
else
    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'Ts'), Ts = cfg.ctrl.Ts;
    else, Ts = cfg.sim.dt; end
end

% Nt = numel(t);
% Nm = beam.Nm;
% if isfield(cfg,'Na') && ~isempty(cfg.Na)
%     Na2 = cfg.Na;
% else
%     Na2 = size(X,1) - (4*Nm + 1 + 3);
%     if Na2 < 0
%         error('postProcess: cannot infer Na from X size. Provide cfg.Na.');
%     end
% end

% % Keep legacy behavior: override X from cached log in non-openLoop
% if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
%     try
%         openL_log = load(fullfile(fullfile(cfg.paths.cache, cfg.case_name), 'log'),'log');
%         if isfield(openL_log,'log') && isfield(openL_log.log,'x') && ~isempty(openL_log.log.x)
%             X  = openL_log.log.x;
%             Nt = size(X,2);
%         end
%     catch
%     end
% end

% ------------------------ modal partitions (legacy indexing) -------------
q1   =  X( 1:Nm              , :);
q2   =  X( Nm+(1:Nm)         , :);
qxi  =  X( 2*Nm+(1:Nm+1)     , :);
qGam =  X( 3*Nm+1+(1:Na2)    , :);
chi  =  X( end-2:end         , :);

if ~isfield(log,'xhat') || isempty(log.xhat), log.xhat = zeros(size(X)); end
q1_hat   =  log.xhat( 1:Nm              , :);
q2_hat   =  log.xhat( Nm+(1:Nm)         , :);
qxi_hat  =  log.xhat( 2*Nm+(1:Nm+1)     , :);
qGam_hat =  log.xhat( 3*Nm+1+(1:Na2)    , :);
chi_hat  =  log.xhat( end-2:end         , :);

% ------------------------ physical reconstruction -----------------------
q0      = -beam.Omega \ q2;
q0_hat  = -beam.Omega \ q2_hat;

x0      = beam.phi0 * q0;     % [6N x Nt]
x1      = beam.phi1 * q1;
x2      = beam.phi2 * q2;
% disp(size(x0))
% disp(size(beam.phi0 ))
x0_hat  = beam.phi0 * q0_hat;
x1_hat  = beam.phi1 * q1_hat;
x2_hat  = beam.phi2 * q2_hat;

% ------------------------ ξ-field & A-frame rotation ---------------------
Xi_full      = base.phi_xi_modes * qxi;
Xi_hat_full  = base.phi_xi_modes * qxi_hat;

Nnode       = beam.fem.num_node - 1;        % legacy
NdofPerNode = beam.nFlex / Nnode;
% selTip legacy choice (unchanged); true tip kept as comment
% selTip = double((Nnode-1)*NdofPerNode) + (1:3)
% selTip = double((Nnode/2-1)*(NdofPerNode+1)) + (1:3)
% selTip = double((Nnode/2)*(NdofPerNode+1)) + (1:3);
selTip = double((Nnode/2+1)*(NdofPerNode+1)) + (1:3);

cellT   = repmat({eye(7,6)},1,Nnode);
% cellT   = repmat({eye(6,6)},1,Nnode);
Tcache  = blkdiag(cellT{:});
T_hat   = Tcache;

TipA         = zeros(Nt,3);
vel_tip      = zeros(Nt,3);

ctrN = max(1, round(Ts/max(eps,cfg.sim.dt)));
nHat = max(1, floor(Nt/ctrN));
TipA_hat      = zeros(nHat,3);
vel_tip_hat   = zeros(nHat,3);

kSam = 0;
for k = 1:Nt
    quatAll = reshape( Xi_full(:,k) , 4,[] ).';
    for n = 1:Nnode
        idx6  = double((n-1)*6) + (1:6);
        idx7  = double((n-1)*7) + (1:7);
        Tcache(idx7,idx6) = blkdiag( ...
            AeroFlex.core.T_phi_quat(quatAll(n,:)), ...
            AeroFlex.core.halfTangentialOperator(quatAll(n,:)) );
        % Tcache(idx6,idx6) = blkdiag( ...
        %     AeroFlex.core.T_phi_quat(quatAll(n,:)), ...
        %     zeros(3,3) );
    end
    % if k<2
    % % disp((Tcache));
    % end
    x0A = Tcache * x0(:,k);
    % x0A = x0(:,k);
    x1A = Tcache * x1(:,k);
    % x0A =  x0(:,k);
    % x1A =  x1(:,k);
    TipA(k,:)    = x0A(selTip).';
    vel_tip(k,:) = x1A(selTip).';
    quatTipT(k,:) = x0A(selTip(end)+(1:4));
    % quatTipT(k,:) = x0(selTip(end)+(1:4));
    eulTipSel = quat2eul(quatTipT(k,:));
    for b= 1:length(eulTipSel)
        if eulTipSel(b) < -pi/2
            eulTipSel(b) = eulTipSel(b)+ 2*pi;
        end
    end
    eulTipT(k,:) = eulTipSel;

    x0A_sim(:,k) = x0A;
    if mod(k, ctrN) == 0
        kSam = kSam + 1;
        if kSam <= size(Xi_hat_full,2)
            quat_hat_All = reshape( Xi_hat_full(:,kSam) , 4,[] ).';
            for n = 1:Nnode
                idx6  = double((n-1)*6) + (1:6);
                idx7  = double((n-1)*7) + (1:7);
                T_hat(idx7,idx6) = blkdiag( ...
                    AeroFlex.core.T_phi_quat(quat_hat_All(n,:)), ...
                    AeroFlex.core.halfTangentialOperator(quat_hat_All(n,:)) );
            end
            x0_hatA = T_hat * x0_hat(:,kSam);
            x1_hatA = T_hat * x1_hat(:,kSam);
            TipA_hat(kSam,:)  = x0_hatA(selTip).';
            vel_tip_hat(kSam,:) = x1_hatA(selTip).';
        end
    end
end
t_ctr = linspace(0, cfg.sim.t_end, size(TipA_hat, 1));

% ------------------------ aero force proxy (legacy) ----------------------
try
    Faero = (aero.forceMap.C_Gamma * qGam).';  % Nt x Nm
catch
    Faero = zeros(Nt, Nm);
end

% ------------------------ package legacy outputs -------------------------
out.t        = t(:).';      % keep your row vector
out.t_col    = t(:);        % convenience col vector
out.disp     = TipA;        % Nx3 tip displacement (A-frame)
out.vel      = vel_tip;     % Nx3 tip velocity (A-frame)
out.TipA_hat = TipA_hat;    % estimator sample track
out.t_ctr    = t_ctr;
out.q0 = q0; out.q1 = q1; out.q2 = q2; out.qxi = qxi; 
% out.chi = chi;
out.chi = chi+ base.FM.chi_bar;
out.Faero    = Faero;


%======================================================================
% Rigid-body outputs, if coupled state is available
%======================================================================
out.raw_X_input = out_raw_X;
out.nxFlex = nxFlex;
out.hasRigidBody = hasRB;

if hasRB
    % Match rigid logs to output time length.
    NtRB = min(size(rb.state,2), Nt);

    rb_t = t(1:NtRB);

    out.rb = struct();
    out.rb.t        = rb_t(:).';
    out.rb.state    = rb.state(:,1:NtRB);
    out.rb.r_I      = rb.r_I(:,1:NtRB);
    out.rb.v_B      = rb.v_B(:,1:NtRB);
    out.rb.euler    = rb.euler(:,1:NtRB);
    out.rb.omega_B  = rb.omega_B(:,1:NtRB);

    out.rb.euler_deg   = rad2deg(out.rb.euler);
    out.rb.omega_deg_s = rad2deg(out.rb.omega_B);

    uB = out.rb.v_B(1,:);
    vB = out.rb.v_B(2,:);
    wB = out.rb.v_B(3,:);

    out.rb.U     = sqrt(uB.^2 + vB.^2 + wB.^2);
    out.rb.alpha = atan2(wB,uB);
    out.rb.beta  = asin(max(-1,min(1,vB ./ max(out.rb.U,eps))));

    out.rb.alpha_deg = rad2deg(out.rb.alpha);
    out.rb.beta_deg  = rad2deg(out.rb.beta);
else
    out.rb = struct();
end

%======================================================================
% Coupled load outputs
%======================================================================
out.loads = struct();

if isfield(log,'loads') && isstruct(log.loads)

    loadFields = fieldnames(log.loads);

    for ii = 1:numel(loadFields)
        f = loadFields{ii};
        val = log.loads.(f);

        if isnumeric(val) && ~isempty(val)
            % Force/moment logs are usually Nload x NtStep, one less than state time.
            out.loads.(f) = val;
        end
    end

    if isfield(out.loads,'Clamp6') && ~isempty(out.loads.Clamp6)

        % Clamp6 = [Fx; Fy; Fz; Mx; My; Mz].
        if isfield(cfg,'metrics') && isfield(cfg.metrics,'rootMomentIndex')
            iMR = cfg.metrics.rootMomentIndex;
        else
            iMR = 5;  % default: My as vertical bending root moment
        end

        iMR = max(1,min(6,iMR));

        MR = out.loads.Clamp6(iMR,:);

        % Loads are usually logged for each plant/control step, length Nt-1.
        if numel(MR) == numel(t)-1
            tLoad = t(1:end-1);
        else
            tLoad = linspace(t(1),t(end),numel(MR));
        end

        out.loads.t = tLoad(:).';
        out.loads.M_R = MR(:).';
        out.loads.M_R_abs = abs(MR(:).');
    end
end

%======================================================================
% Control / estimator / solver outputs
%======================================================================
out.control = struct();
out.estimator = struct();
out.solver = struct();

% -------------------- controls ---------------------------------------
if isfield(log,'U') && isnumeric(log.U) && ~isempty(log.U)
    Ulog = log.U;

    % Expected Ulog = nu x Nctrl.
    if size(Ulog,1) > size(Ulog,2) || size(Ulog,1) <= 20
        out.control.U = Ulog;
    else
        out.control.U = Ulog.';
    end

    nU = size(out.control.U,1);
    nK = size(out.control.U,2);

    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'Ts')
        tU = (0:nK-1)*cfg.ctrl.Ts;
    else
        tU = linspace(t(1),t(end),nK);
    end

    out.control.t = tU;

    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'n_surf')
        nSurf = cfg.ctrl.n_surf;
    else
        nSurf = min(2,nU);
    end

    if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'var_per')
        varPer = cfg.ctrl.var_per;
    else
        varPer = 1;
    end

    if varPer == 2 && nU >= 2*nSurf
        out.control.delta = out.control.U(1:nSurf,:);
        out.control.rate_cmd = out.control.U(nSurf+1:2*nSurf,:);
    else
        out.control.delta = out.control.U(1:min(nSurf,nU),:);
        out.control.rate_cmd = [];
    end

    % Finite-difference realized command rate from deflection command.
    if size(out.control.delta,2) >= 2
        dtU = mean(diff(tU));
        if isfinite(dtU) && dtU > 0
            out.control.rate_fd = [zeros(size(out.control.delta,1),1), ...
                                   diff(out.control.delta,1,2)/dtU];
        else
            out.control.rate_fd = zeros(size(out.control.delta));
        end
    else
        out.control.rate_fd = zeros(size(out.control.delta));
    end

    out.control.delta_deg = rad2deg(out.control.delta);
    out.control.rate_fd_deg_s = rad2deg(out.control.rate_fd);

    if ~isempty(out.control.rate_cmd)
        out.control.rate_cmd_deg_s = rad2deg(out.control.rate_cmd);
    end
end

% -------------------- gusts / MHE disturbance -------------------------
if isfield(log,'wTrue') && isnumeric(log.wTrue)
    out.estimator.wTrue = log.wTrue(:).';
end

if isfield(log,'wHorizon') && isnumeric(log.wHorizon)
    out.estimator.wHorizon = log.wHorizon;

    % If wHorizon is Ne x Nt, the newest slot is the last row.
    if size(log.wHorizon,1) <= size(log.wHorizon,2)
        out.estimator.wHatCurrent = log.wHorizon(end,:);
    else
        out.estimator.wHatCurrent = log.wHorizon(:,end).';
    end
end

% -------------------- solver diagnostics ------------------------------
if isfield(log,'diag') && isstruct(log.diag)
    out.solver = log.diag;
end

% % existing lines
% out.t    = t.';               out.disp = TipA;
% out.q0   = q0;                out.q1  = q1;  out.q2 = q2;
% out.qxi  = qxi;               out.chi = chi; 

% >>> ADD: span-based normalization that matches SHARPy (b_ref in SHARPy is semi-span)
span_half = cfg.struct.L/2;         % your span is cfg.struct.L (1.1 m for Pazy), half-span = 0.55 m
out.tip_z = out.disp(:,3);          % A-frame vertical tip [m]
out.tip_z_norm_pct = 100 * (out.tip_z ./ span_half);   % %. of (b/2)

% ------------------------ plots (unchanged content) + saving -------------

% disp(size(x0A_sim));
% 1) Tip estimation overlay (non-openLoop only)
if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
    fh = figure('Color','w','Position',[100 100 1120 520]);
    plot(t_ctr,TipA_hat(:,1)*100/cfg.debug.plt_scale); hold on
    plot(t_ctr,TipA_hat(:,2)*100/cfg.debug.plt_scale);
    plot(t_ctr,TipA_hat(:,3)*100/cfg.debug.plt_scale);
    xlabel('t [s]'); ylabel('Tip [% b/2]'); grid on
    legend({'x','y','z'},'Location','best')
    title('Wing-tip estimation (A-frame, % span/2)');
    hold off
    save_plot(fh, 'tip_estimation');
end

% 2) Modal q0/q1/q2 (debugPlots)
if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'debugPlots') && cfg.plotOpts.debugPlots
    fh = figure('Color','w','Position',[100 100 1120 840]);
    tiledlayout(3,1,'TileSpacing','tight');
    nexttile; plot(t,q0'); title('q_0'); grid on
    nexttile; plot(t,q1'); title('q_1'); grid on
    nexttile; plot(t,q2'); title('q_2'); grid on
    save_plot(fh, 'modal_qs');
end

% 3) Wing-tip heave (true vs estimate if available)
if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'wingTip') && cfg.plotOpts.wingTip
    fh = figure('Color','w','Position',[100 100 1120 520]);
    % plot(t,  out.disp(:,3)*100/cfg.debug.plt_scale); hold on
    plot(t,  out.tip_z_norm_pct); hold on
    if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
        plot(t_ctr, out.TipA_hat(:,3)*100/cfg.debug.plt_scale);
        legend({'True','Estimate'},'Location','best');
    end
    xlabel('t [s]'); ylabel('Tip z [% b/2]'); grid on
    title('Wing-tip heave (A-frame, % span/2)');
    hold off
    save_plot(fh, 'wingtip_heave');
end

% figure; plot(t, eulTipT); title('Test');
% ------------------------ energies (legacy math + plots) -----------------
Tkin = 0.5*sum(q1.^2, 1);       % 1 x Nt
Vpot = 0.5*sum(q2.^2, 1);       % 1 x Nt
if     cfg.struct.damping

    Edot = sum((q1.') * ((beam.Sigma)*q1 + cfg.flight.a*cfg.flight.t_inf*Faero.'), 1);
else
    Edot = sum((q1.') * (diag(beam.Sigma)*q1 + cfg.flight.a*cfg.flight.t_inf*Faero.'), 1);

end
Etot = Tkin + Vpot;

out.Tkin = Tkin.';  out.Vpot = Vpot.';  out.Etot = Etot.';

if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'debugPlots') && cfg.plotOpts.debugPlots
    fh = figure('Name','Energy check','Color','w','Position',[100 100 1120 520]);
    plot(t, Tkin, 'b-', t, Vpot, 'r--', t, Etot, 'k');
    grid on; xlabel('t  [s]'); ylabel('Energy  [J]');
    legend('Kinetic','Potential','Total','Location','best');
    title('Structural modal energy (flexible part)');
    save_plot(fh, 'energy_check');

    fh = figure('Name','Energy Rate check','Color','w','Position',[100 100 1120 520]);
    plot(t, Edot, 'b-');
    grid on; xlabel('t  [s]'); ylabel('Energy rate  [J/s]');
    title('Structural modal energy rate');
    save_plot(fh, 'energy_rate');
end


%%
node_xyz_hist = cell(1,Nnode);
node_xyzQuat_hist = cell(1,Nnode);
coords_hist = cell(1, length(t));
% if ~isfolder("node_xyz_hist"), mkdir("node_xyz_hist"), end
% if ~isfolder("node_xyzQuat_hist"), mkdir("node_xyzQuat_hist"), end
coords = beam.fem.coordinates;

% for k = 1:Nnode
%     idx_xyz = (k -1)*NdofPerNode +(1:NdofPerNode );
%     node_xyz_hist{k} = x0(idx_xyz,:).';
% 
%     idx_xyzQuat = (k -1)*(NdofPerNode+1) +(1:NdofPerNode +1);
%     node_xyzQuat_hist{k}=x0A(idx_xyzQuat,:).';
% 
% 
% 
% 
%     % writematrix(node_xyz_hist{k},               "node_xyz_hist/"+ int2str(k)+".xls");
%     % writematrix(node_xyzQuat_hist{k},               "node_xyzQuat_hist/"+ int2str(k)+".xls");
%     % % 
%     % writematrix(node_xyz_hist{k},               "node_xyz_hist/"+ int2str(k)+".csv");
%     % writematrix(node_xyzQuat_hist{k},               "node_xyzQuat_hist/"+ int2str(k)+".csv");
% 
% 
% end
%%
% if ~isfolder("timeseries"), mkdir("timeseries"), end
% 
% node_xyz_H = zeros(size(coords));
% for l = 1:length(t)
% 
%     for k = 1:Nnode
%     idx_xyz = double((k -1)*NdofPerNode )+(1:3);
%     node_xyz_H(k+1,:) = x0(idx_xyz,l).';
% 
%     % idx_xyzQuat = (k -1)*(NdofPerNode+1) +(1:NdofPerNode +1);
%     % node_xyzQuat_hist{k}=x0A(idx_xyzQuat,:).';
% 
%     end
%     coords_hist{l} = coords +node_xyz_H;
%     if l<20
%     T = table(t(l)*ones(Nnode+1,1),coords_hist{l}(:,1), coords_hist{l}(:,2), coords_hist{l}(:,3), ...
%               'VariableNames', {'time','x','y','z'});
%     writetable(T, "timeseries/" + int2str(l)+".csv");
%     end
% end
%%
if ~isfolder("timeseries"), mkdir("timeseries"), end

node_xyz_H = zeros(size(length(t),3));
coords_hist_time = cell(1, Nnode+1);
% disp(coords);
for k = 1:Nnode
    idx_xyz = double((k -1)*NdofPerNode )+(1:3);
    node_xyz_H = x0(idx_xyz,:).';
    % disp(node_xyz_H(1:3,1:3))
    % idx_xyzQuat = (k -1)*(NdofPerNode+1) +(1:NdofPerNode +1);
    % node_xyzQuat_hist{k}=x0A(idx_xyzQuat,:).';
    % for l = 1:length(t)
    % 
    % end
    coords_hist_time{k} = coords(k+1,:)+node_xyz_H;
    % if l<20
        % disp(coords_hist_time{k}(1:3,1:3))

    T = table(coords_hist_time{k}(:,1), coords_hist_time{k}(:,2), coords_hist_time{k}(:,3), t(:),...
              'VariableNames', {'x','y','z','time'});
    writetable(T, "timeseries/" + int2str(k)+".csv");
    % end
end

T = table(zeros(length(t),1), zeros(length(t),1),zeros(length(t),1), t(:), ...
              'VariableNames', {'x','y','z','time'});
    writetable(T, "timeseries/0.csv");
%%
coord_hist1 = coords_hist_time{1};
  % writematrix(coords_hist{1},   "coords_hist1.csv");

% writecell(coords_hist, "coords_hist.csv");
%%
% T = table(coord_hist1(:,1), coord_hist1(:,2), coord_hist1(:,3), ...
%               'VariableNames', {'x1','y1','z1'});
%     % writetable(T, 'timeseries.csv', 'Sheet','timeseries');
%     writetable(T, 'timeseries/.csv');

    % T = table(t(:), out.tip_z(:), out.tip_z_norm_pct(:), ...
    %           angles_deg(:,1), angles_deg(:,2), angles_deg(:,3), ...
    %           'VariableNames', {'time','tip_z','tip_z_norm_pct','roll_deg','pitch_deg','yaw_deg'});
    % writetable(T, fullfile(xcell_dir,'timeseries.xlsx'), 'Sheet','timeseries');

% writematrix(x0A,               'x0A.csv');
% writematrix(x0,               'x0.csv');
% writematrix(t,               'time.csv');


% ------------------------ convenience for sim_run exports ----------------
b2 = max(eps, cfg.debug.plt_scale);
out.tip.xyz        = out.disp;                  % alias
out.tip.z_norm_pct = 100*out.disp(:,3)/b2;      % Nt x 1
out.angles.deg     = rad2deg(out.chi.') + rad2deg(ones(size(out.chi.')).*base.chi0.');        % Nt x 3 [roll pitch yaw]
out.energy.Tkin    = out.Tkin;
out.energy.Vpot    = out.Vpot;
out.energy.Etot    = out.Etot;
out.q.q0 = out.q0; out.q.q1 = out.q1; out.q.q2 = out.q2;
out.plots_dir = plots_dir;                      % tell caller where plots live


%======================================================================
% Additional coupled/control/performance plots
%======================================================================

% -------------------- rigid-body plots --------------------------------
if out.hasRigidBody
    fh = figure('Color','w','Position',[100 100 1120 760]);
    tiledlayout(3,1,'TileSpacing','tight');

    nexttile;
    plot(out.rb.t,out.rb.euler_deg(1,:), ...
         out.rb.t,out.rb.euler_deg(2,:), ...
         out.rb.t,out.rb.euler_deg(3,:));
    grid on; ylabel('Euler [deg]');
    legend('\phi','\theta','\psi','Location','best');
    title('Rigid-body attitude');

    nexttile;
    plot(out.rb.t,out.rb.omega_deg_s(1,:), ...
         out.rb.t,out.rb.omega_deg_s(2,:), ...
         out.rb.t,out.rb.omega_deg_s(3,:));
    grid on; ylabel('Rates [deg/s]');
    legend('p','q','r','Location','best');

    nexttile;
    plot(out.rb.t,out.rb.U, ...
         out.rb.t,out.rb.alpha_deg, ...
         out.rb.t,out.rb.beta_deg);
    grid on; xlabel('t [s]');
    ylabel('Flight condition');
    legend('U [m/s]','\alpha [deg]','\beta [deg]','Location','best');

    save_plot(fh,'rigid_body_states');

    fh = figure('Color','w','Position',[100 100 1120 520]);
    plot3(out.rb.r_I(1,:),out.rb.r_I(2,:),out.rb.r_I(3,:));
    grid on; axis equal;
    xlabel('x_I [m]'); ylabel('y_I [m]'); zlabel('z_I [m]');
    title('Rigid-body inertial trajectory');
    save_plot(fh,'rigid_body_trajectory');
end

% -------------------- control plots -----------------------------------
if isfield(out,'control') && isfield(out.control,'U') && ~isempty(out.control.U)

    fh = figure('Color','w','Position',[100 100 1120 760]);

    if isfield(out.control,'rate_cmd') && ~isempty(out.control.rate_cmd)
        tiledlayout(3,1,'TileSpacing','tight');

        nexttile;
        plot(out.control.t,out.control.delta_deg.');
        grid on; ylabel('\delta [deg]');
        title('Control surface deflections');
        legend(local_surface_labels(size(out.control.delta,1),'\delta'),'Location','best');

        nexttile;
        plot(out.control.t,out.control.rate_cmd_deg_s.');
        grid on; ylabel('\dot{\delta}_{cmd} [deg/s]');
        title('Commanded actuator-rate channels');
        legend(local_surface_labels(size(out.control.rate_cmd,1),'\dot{\delta}_{cmd}'),'Location','best');

        nexttile;
        plot(out.control.t,out.control.rate_fd_deg_s.');
        grid on; xlabel('t [s]'); ylabel('\Delta\delta/T_s [deg/s]');
        title('Finite-difference realized deflection rate');
        legend(local_surface_labels(size(out.control.rate_fd,1),'\dot{\delta}_{FD}'),'Location','best');
    else
        tiledlayout(2,1,'TileSpacing','tight');

        nexttile;
        plot(out.control.t,out.control.delta_deg.');
        grid on; ylabel('\delta [deg]');
        title('Control surface deflections');
        legend(local_surface_labels(size(out.control.delta,1),'\delta'),'Location','best');

        nexttile;
        plot(out.control.t,out.control.rate_fd_deg_s.');
        grid on; xlabel('t [s]'); ylabel('\Delta\delta/T_s [deg/s]');
        title('Finite-difference realized deflection rate');
        legend(local_surface_labels(size(out.control.rate_fd,1),'\dot{\delta}_{FD}'),'Location','best');
    end

    save_plot(fh,'control_commands');
end

% -------------------- gust estimate plot -------------------------------
if isfield(out.estimator,'wTrue') || isfield(out.estimator,'wHatCurrent')
    fh = figure('Color','w','Position',[100 100 1120 520]);
    hold on;

    if isfield(out.estimator,'wTrue')
        tw = linspace(t(1),t(end),numel(out.estimator.wTrue));
        plot(tw,out.estimator.wTrue,'k--');
    end

    if isfield(out.estimator,'wHatCurrent')
        th = linspace(t(1),t(end),numel(out.estimator.wHatCurrent));
        plot(th,out.estimator.wHatCurrent,'b-');
    end

    hold off; grid on;
    xlabel('t [s]'); ylabel('w_g');
    legend({'true gust','MHE estimate'},'Location','best');
    title('Gust / disturbance estimate');
    save_plot(fh,'gust_estimate');
end

% -------------------- loads plots --------------------------------------
if isfield(out,'loads') && isfield(out.loads,'M_R')
    fh = figure('Color','w','Position',[100 100 1120 520]);
    plot(out.loads.t,out.loads.M_R);
    grid on; xlabel('t [s]'); ylabel('M_R [N m]');
    title('Wing-root bending moment proxy');
    save_plot(fh,'root_bending_moment');
end

if isfield(out,'loads') && isfield(out.loads,'Ftot_B') && isfield(out.loads,'Mtot_B')
    tLoad = linspace(t(1),t(end),size(out.loads.Ftot_B,2));

    fh = figure('Color','w','Position',[100 100 1120 760]);
    tiledlayout(2,1,'TileSpacing','tight');

    nexttile;
    plot(tLoad,out.loads.Ftot_B.');
    grid on; ylabel('F_B [N]');
    legend('F_x','F_y','F_z','Location','best');
    title('Total rigid-body force');

    nexttile;
    plot(tLoad,out.loads.Mtot_B.');
    grid on; xlabel('t [s]'); ylabel('M_B [N m]');
    legend('M_x','M_y','M_z','Location','best');
    title('Total rigid-body moment');

    save_plot(fh,'total_forces_moments');
end

% -------------------- solver diagnostics -------------------------------
if isfield(out,'solver') && isfield(out.solver,'tCtrl') && ~isempty(out.solver.tCtrl)
    valid = isfinite(out.solver.tCtrl);

    if any(valid)
        tc = out.solver.tCtrl(valid);

        fh = figure('Color','w','Position',[100 100 1120 760]);
        tiledlayout(3,1,'TileSpacing','tight');

        nexttile;
        if isfield(out.solver,'mheTime') && isfield(out.solver,'mpcTime')
            plot(tc,1e3*out.solver.mheTime(valid), ...
                 tc,1e3*out.solver.mpcTime(valid));
            legend('MHE','MPC','Location','best');
        end
        grid on; ylabel('solve [ms]');
        title('Estimator/controller solve times');

        nexttile;
        if isfield(out.solver,'mheCont') && isfield(out.solver,'mpcCont')
            semilogy(tc,out.solver.mheCont(valid), ...
                     tc,out.solver.mpcCont(valid));
            legend('MHE','MPC','Location','best');
        end
        grid on; ylabel('||c_{eq}||');

        nexttile;
        if isfield(out.solver,'mheFlag') && isfield(out.solver,'mpcFlag')
            stairs(tc,out.solver.mheFlag(valid)); hold on;
            stairs(tc,out.solver.mpcFlag(valid));
            hold off;
            legend('MHE','MPC','Location','best');
        end
        grid on; xlabel('t [s]'); ylabel('exit flag');

        save_plot(fh,'solver_diagnostics');
    end
end

%======================================================================
% Performance metrics and summary tables
%======================================================================
metrics = struct();

% -------------------- loads and deflections ---------------------------
metrics.tip_z_peak_abs = max(abs(out.tip_z),[],'omitnan');
metrics.tip_z_rms      = sqrt(mean(out.tip_z.^2,'omitnan'));

metrics.tip_z_norm_pct_peak_abs = max(abs(out.tip_z_norm_pct),[],'omitnan');
metrics.tip_z_norm_pct_rms      = sqrt(mean(out.tip_z_norm_pct.^2,'omitnan'));

if isfield(out,'loads') && isfield(out.loads,'M_R')
    MR = out.loads.M_R(:);
    metrics.M_R_peak_abs = max(abs(MR),[],'omitnan');
    metrics.M_R_rms      = sqrt(mean(MR.^2,'omitnan'));

    if isfield(out.loads,'t') && numel(out.loads.t) == numel(MR)
        Tload = out.loads.t(end) - out.loads.t(1);
        if Tload > 0
            metrics.M_R_rms_trapz = sqrt((1/Tload)*trapz(out.loads.t,MR.'.^2));
        else
            metrics.M_R_rms_trapz = metrics.M_R_rms;
        end
    else
        metrics.M_R_rms_trapz = metrics.M_R_rms;
    end
else
    metrics.M_R_peak_abs = nan;
    metrics.M_R_rms = nan;
    metrics.M_R_rms_trapz = nan;
end

% -------------------- rigid-body tracking -----------------------------
if out.hasRigidBody
    eul = out.rb.euler;

    metrics.roll_peak_deg  = max(abs(rad2deg(eul(1,:))),[],'omitnan');
    metrics.pitch_peak_deg = max(abs(rad2deg(eul(2,:))),[],'omitnan');
    metrics.yaw_peak_deg   = max(abs(rad2deg(eul(3,:))),[],'omitnan');

    metrics.roll_rms_deg  = sqrt(mean(rad2deg(eul(1,:)).^2,'omitnan'));
    metrics.pitch_rms_deg = sqrt(mean(rad2deg(eul(2,:)).^2,'omitnan'));
    metrics.yaw_rms_deg   = sqrt(mean(rad2deg(eul(3,:)).^2,'omitnan'));

    metrics.U_mean = mean(out.rb.U,'omitnan');
    metrics.U_peak_dev = max(abs(out.rb.U - metrics.U_mean),[],'omitnan');
else
    metrics.roll_peak_deg = nan;
    metrics.pitch_peak_deg = nan;
    metrics.yaw_peak_deg = nan;
    metrics.roll_rms_deg = nan;
    metrics.pitch_rms_deg = nan;
    metrics.yaw_rms_deg = nan;
    metrics.U_mean = nan;
    metrics.U_peak_dev = nan;
end

% -------------------- actuator usage ----------------------------------
if isfield(out,'control') && isfield(out.control,'delta')
    delta = out.control.delta;
    rateFD = out.control.rate_fd;

    metrics.u_peak_rad = max(abs(delta(:)),[],'omitnan');
    metrics.u_rms_rad  = sqrt(mean(delta(:).^2,'omitnan'));
    metrics.u_peak_deg = rad2deg(metrics.u_peak_rad);
    metrics.u_rms_deg  = rad2deg(metrics.u_rms_rad);

    metrics.udot_peak_rad_s = max(abs(rateFD(:)),[],'omitnan');
    metrics.udot_rms_rad_s  = sqrt(mean(rateFD(:).^2,'omitnan'));
    metrics.udot_peak_deg_s = rad2deg(metrics.udot_peak_rad_s);
    metrics.udot_rms_deg_s  = rad2deg(metrics.udot_rms_rad_s);

    if isfield(cfg,'uL') && isfield(cfg,'uU')
        uLim = max(abs([cfg.uL(:); cfg.uU(:)]));
        if isscalar(uLim) && isfinite(uLim) && uLim > 0
            metrics.u_saturation_fraction = mean(abs(delta(:)) >= 0.999*uLim,'omitnan');
        else
            metrics.u_saturation_fraction = nan;
        end
    else
        metrics.u_saturation_fraction = nan;
    end
else
    metrics.u_peak_rad = nan;
    metrics.u_rms_rad = nan;
    metrics.u_peak_deg = nan;
    metrics.u_rms_deg = nan;
    metrics.udot_peak_rad_s = nan;
    metrics.udot_rms_rad_s = nan;
    metrics.udot_peak_deg_s = nan;
    metrics.udot_rms_deg_s = nan;
    metrics.u_saturation_fraction = nan;
end

% -------------------- solver metrics ----------------------------------
if isfield(out,'solver') && isfield(out.solver,'mheTime')
    valid = isfinite(out.solver.mheTime);
    metrics.mhe_solve_mean_ms = 1e3*mean(out.solver.mheTime(valid),'omitnan');
    metrics.mhe_solve_max_ms  = 1e3*max(out.solver.mheTime(valid),[],'omitnan');
else
    metrics.mhe_solve_mean_ms = nan;
    metrics.mhe_solve_max_ms = nan;
end

if isfield(out,'solver') && isfield(out.solver,'mpcTime')
    valid = isfinite(out.solver.mpcTime);
    metrics.mpc_solve_mean_ms = 1e3*mean(out.solver.mpcTime(valid),'omitnan');
    metrics.mpc_solve_max_ms  = 1e3*max(out.solver.mpcTime(valid),[],'omitnan');
else
    metrics.mpc_solve_mean_ms = nan;
    metrics.mpc_solve_max_ms = nan;
end

if isfield(out,'solver') && isfield(out.solver,'mheFlag')
    valid = isfinite(out.solver.mheFlag);
    metrics.mhe_success_pct = 100*mean(out.solver.mheFlag(valid) > 0,'omitnan');
else
    metrics.mhe_success_pct = nan;
end

if isfield(out,'solver') && isfield(out.solver,'mpcFlag')
    valid = isfinite(out.solver.mpcFlag);
    metrics.mpc_success_pct = 100*mean(out.solver.mpcFlag(valid) > 0,'omitnan');
else
    metrics.mpc_success_pct = nan;
end

out.metrics = metrics;

% -------------------- summary table -----------------------------------
metricNames = fieldnames(metrics);
metricVals  = nan(numel(metricNames),1);

for ii = 1:numel(metricNames)
    v = metrics.(metricNames{ii});
    if isnumeric(v) && isscalar(v)
        metricVals(ii) = v;
    end
end

out.metricsTable = table(metricNames,metricVals, ...
    'VariableNames',{'Metric','Value'});

% Export summary table.
try
    writetable(out.metricsTable,fullfile(plots_dir,'performance_metrics.csv'));
catch ME
    warning('postProcess:MetricsExport', ...
            'Could not write performance_metrics.csv: %s', ME.message);
end

% -------------------- print compact summary ---------------------------
fprintf('\n');
fprintf('======================================================================\n');
fprintf(' POST-PROCESS PERFORMANCE SUMMARY\n');
fprintf('======================================================================\n');
fprintf('  Tip z peak abs        : %10.4e m\n', metrics.tip_z_peak_abs);
fprintf('  Tip z RMS             : %10.4e m\n', metrics.tip_z_rms);
fprintf('  Tip z peak abs        : %10.4f %% b/2\n', metrics.tip_z_norm_pct_peak_abs);
fprintf('  Tip z RMS             : %10.4f %% b/2\n', metrics.tip_z_norm_pct_rms);
fprintf('  M_R peak abs          : %10.4e N m\n', metrics.M_R_peak_abs);
fprintf('  M_R RMS               : %10.4e N m\n', metrics.M_R_rms_trapz);
fprintf('  Control peak          : %10.4f deg\n', metrics.u_peak_deg);
fprintf('  Control RMS           : %10.4f deg\n', metrics.u_rms_deg);
fprintf('  Control rate peak     : %10.4f deg/s\n', metrics.udot_peak_deg_s);
fprintf('  Control rate RMS      : %10.4f deg/s\n', metrics.udot_rms_deg_s);
if out.hasRigidBody
fprintf('  Roll/Pitch/Yaw RMS    : %10.4f / %10.4f / %10.4f deg\n', ...
        metrics.roll_rms_deg, metrics.pitch_rms_deg, metrics.yaw_rms_deg);
end
fprintf('  MHE success           : %10.2f %%\n', metrics.mhe_success_pct);
fprintf('  MPC success           : %10.2f %%\n', metrics.mpc_success_pct);
fprintf('  MHE solve mean/max    : %10.3f / %10.3f ms\n', ...
        metrics.mhe_solve_mean_ms, metrics.mhe_solve_max_ms);
fprintf('  MPC solve mean/max    : %10.3f / %10.3f ms\n', ...
        metrics.mpc_solve_mean_ms, metrics.mpc_solve_max_ms);
fprintf('======================================================================\n\n');

if isfield(log,'baseline') && isfield(log.baseline,'M_R_open') && ...
        isfield(out,'loads') && isfield(out.loads,'M_R')

    MRopen = log.baseline.M_R_open(:);
    MRclosed = out.loads.M_R(:);

    maxOpen = max(abs(MRopen),[],'omitnan');
    maxClosed = max(abs(MRclosed),[],'omitnan');

    rmsOpen = sqrt(mean(MRopen.^2,'omitnan'));
    rmsClosed = sqrt(mean(MRclosed.^2,'omitnan'));

    metrics.GLA_peak_pct = 100*(maxOpen - maxClosed)/max(maxOpen,eps);
    metrics.GLA_rms_pct  = 100*(rmsOpen - rmsClosed)/max(rmsOpen,eps);

else
    metrics.GLA_peak_pct = nan;
    metrics.GLA_rms_pct = nan;
end
fprintf('  GLA peak reduction    : %10.3f %%\n', metrics.GLA_peak_pct);
fprintf('  GLA RMS reduction     : %10.3f %%\n', metrics.GLA_rms_pct);
end % postProcess


% ============================== helpers ===================================
function labels = local_surface_labels(n,prefix)
%LOCAL_SURFACE_LABELS Generate legend labels for surface channels.
    labels = cell(1,n);
    for ii = 1:n
        labels{ii} = sprintf('%s_%d',prefix,ii);
    end
end
function plots_dir = local_get_plots_dir(cfg)
% Prefer cfg.paths.plots (set by sim_init/sim_run). If unavailable,
% fall back to a local ./plots next to the current working directory.
plots_dir = '';
if isfield(cfg,'paths')
    if isfield(cfg.paths,'plots') && ~isempty(cfg.paths.plots)
        plots_dir = cfg.paths.plots;
        return
    elseif isfield(cfg.paths,'reports') && ~isempty(cfg.paths.reports)
        plots_dir = fullfile(cfg.paths.reports,'plots');
        return
    end
end
plots_dir = fullfile(pwd,'plots');
end

function local_export_fig(fh, fullbase)
% Save PNG (300 dpi), PDF (vector), and FIG. Robust to MATLAB version.
png = [fullbase '.png']; pdf = [fullbase '.pdf']; fig = [fullbase '.fig'];
try
    exportgraphics(fh, png, 'Resolution', 300);
catch
    try, print(fh, '-dpng', '-r300', png); catch, saveas(fh, png); end
end
try
    exportgraphics(fh, pdf);
catch
    try, print(fh, '-dpdf', pdf); catch, saveas(fh, [fullbase '.pdf']); end
end
try, savefig(fh, fig); catch, end
% if ishghandle(fh), close(fh); end

end


% function out = postProcess(t, X, cfg, beam, aero, base, log)
% % POSTPROCESS  –  decode SimRunner state history, return physical signals
% %
% %   • Keeps legacy behavior/plots (so your "old" graphs match)
% %   • q0 = -Ω^{-1} q2  (Artola App. (46))  — no double integration
% %   • Quaternion-to-A frame rotation per node (same tip selection)
% %   • Enriches `out` with fields used by sim_run (tip %, angles, energies)
% %
% % --------------------------------------------------------------------
% arguments
%     t                        % (Nt x 1) double or row
%     X                        % (Nx x Nt) double
%     cfg                      % struct
%     beam                     % BeamModel
%     aero                     % AeroROM
%     base                     % struct  (phi_xi_modes, etc.)
%     log                      % struct (may omit fields; robust defaults)
% end
% 
% % -------------------- sampler step for estimator alignment ---------------
% if isfield(cfg,'case') && strcmpi(cfg.case,'openLoop')
%     % mimic legacy behavior: no estimator overlay in openLoop
%     log.xhat = zeros(size(X));
%     Ts = cfg.sim.dt;
% else
%     if isfield(cfg,'ctrl') && isfield(cfg.ctrl,'Ts')
%         Ts = cfg.ctrl.Ts;
%     else
%         Ts = cfg.sim.dt; % safe fallback
%     end
% end
% 
% % -------------------- visual defaults ------------------------------------
% set(0,'DefaultAxesFontSize',14,...
%       'DefaultLegendFontSize',14,...
%       'DefaultLineLineWidth',2,...
%       'DefaultAxesLineWidth',1.2)
% 
% Nt = numel(t);
% Nm = beam.Nm;
% 
% % Robust Na inference if cfg.Na not present
% if isfield(cfg,'Na') && ~isempty(cfg.Na)
%     Na2 = cfg.Na;
% else
%     % X stacking: q1(Nm) + q2(Nm) + qxi(Nm+1) + qGam(Na) + chi(3)
%     Na2 = size(X,1) - (4*Nm + 1 + 3);
%     if Na2 < 0
%         error('postProcess: cannot infer Na from X size. Provide cfg.Na.');
%     end
% end
% 
% % Keep legacy cache override for non-openLoop (unchanged behavior)
% if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
%     try
%         openL_log = load(fullfile(fullfile(cfg.paths.cache, cfg.case_name), 'log'),'log');
%         if isfield(openL_log,'log') && isfield(openL_log.log,'x') && ~isempty(openL_log.log.x)
%             X = openL_log.log.x;  % legacy override
%             Nt = size(X,2);
%         end
%     catch
%         % ignore if not present
%     end
% end
% 
% % -------------------- modal partitions (legacy indexing) ------------------
% q1   =  X( 1:Nm              , :);
% q2   =  X( Nm+(1:Nm)         , :);
% qxi  =  X( 2*Nm+(1:Nm+1)     , :);
% qGam =  X( 3*Nm+1+(1:Na2)    , :);
% chi  =  X( end-2:end         , :);    % roll,pitch,yaw (rad), legacy
% 
% % Estimator overlays (robust if absent)
% if ~isfield(log,'xhat') || isempty(log.xhat)
%     log.xhat = zeros(size(X));
% end
% q1_hat   =  log.xhat( 1:Nm              , :);
% q2_hat   =  log.xhat( Nm+(1:Nm)         , :);
% qxi_hat  =  log.xhat( 2*Nm+(1:Nm+1)     , :);
% qGam_hat =  log.xhat( 3*Nm+1+(1:Na2)    , :);
% chi_hat  =  log.xhat( end-2:end         , :);
% 
% % -------------------- displacement/velocity reconstruction ----------------
% % q0 via Artola identity (mass-normalised modal basis => Ω diagonal)
% q0      = -beam.Omega \ q2;
% q0_hat  = -beam.Omega \ q2_hat;
% 
% x0      = beam.phi0 * q0;     % [6N x Nt]
% x1      = beam.phi1 * q1;
% x2      = beam.phi2 * q2;
% 
% x0_hat  = beam.phi0 * q0_hat;
% x1_hat  = beam.phi1 * q1_hat;
% x2_hat  = beam.phi2 * q2_hat;
% 
% % -------------------- ξ-field (quaternions) and A-frame rotation ----------
% Xi_full      = base.phi_xi_modes * qxi;      % (4*Nnode) x Nt
% Xi_hat_full  = base.phi_xi_modes * qxi_hat;
% 
% % Legacy node/tip selection (unchanged)
% Nnode       = beam.fem.num_node - 1;         % last is tip (legacy convention)
% NdofPerNode = beam.nFlex / Nnode;            % should be 6
% 
% % Legacy "selTip" (unchanged: uses mid-span style from your old script)
% % Note: The commented "true tip" kept for reference; we keep current one.
% % selTip = double((Nnode-1)*NdofPerNode) + (1:3);              % true tip
% selTip = double((Nnode/2-1)*(NdofPerNode+1)) + (1:3);          % legacy selection
% 
% % Block-diagonal T(q) shape: each block is 7x6 => total 7N x 6N
% cellT   = repmat({eye(7,6)},1,Nnode);
% Tcache  = blkdiag(cellT{:});                  % 7N x 6N
% T_hat   = Tcache;                             % preallocate for overlay
% 
% TipA         = zeros(Nt,3);
% out.vel      = zeros(Nt,3);
% 
% ctrN = max(1, round(Ts/max(eps,cfg.sim.dt)));
% nHat = max(1, floor(Nt/ctrN));
% TipA_hat      = zeros(nHat,3);
% out.vel_hat   = zeros(nHat,3);
% 
% kSam = 0;
% for k = 1:Nt
%     % Current quaternions (Nnode x 4)
%     quatAll = reshape( Xi_full(:,k) , 4,[] ).';
%     % Build Tcache by overwriting blocks
%     for n = 1:Nnode
%         idx6  = double((n-1)*6) + (1:6);
%         idx7  = double((n-1)*7) + (1:7);
%         Tcache(idx7,idx6) = blkdiag( ...
%             AeroFlex.core.T_phi_quat(quatAll(n,:)), ...
%             AeroFlex.core.halfTangentialOperator(quatAll(n,:)) );
%     end
%     % Rotate displacements/velocities to A-frame
%     x0A = Tcache * x0(:,k);
%     x1A = Tcache * x1(:,k);
% 
%     TipA(k,:)   = x0A(selTip).';
%     out.vel(k,:)= x1A(selTip).';
% 
%     % Estimator sampling (only every ctrN steps)
%     if mod(k, ctrN) == 0
%         kSam = kSam + 1;
%         if kSam <= size(Xi_hat_full,2)
%             quat_hat_All = reshape( Xi_hat_full(:,kSam) , 4,[] ).';
%             % Build T_hat for this sample
%             for n = 1:Nnode
%                 idx6  = double((n-1)*6) + (1:6);
%                 idx7  = double((n-1)*7) + (1:7);
%                 T_hat(idx7,idx6) = blkdiag( ...
%                     AeroFlex.core.T_phi_quat(quat_hat_All(n,:)), ...
%                     AeroFlex.core.halfTangentialOperator(quat_hat_All(n,:)) );
%             end
%             x0_hatA = T_hat * x0_hat(:,kSam);
%             x1_hatA = T_hat * x1_hat(:,kSam);
%             TipA_hat(kSam,:)    = x0_hatA(selTip).';
%             out.vel_hat(kSam,:) = x1_hatA(selTip).';
%         end
%     end
% end
% 
% t_ctr = linspace(0, cfg.sim.t_end, size(TipA_hat, 1));
% 
% % -------------------- aero force proxy (legacy) ---------------------------
% try
%     out.Faero = (aero.forceMap.C_Gamma * qGam).';  % Nt x Nm (rough)
% catch
%     out.Faero = zeros(Nt, Nm); % keep field present
% end
% 
% % -------------------- legacy package + plots (UNCHANGED) ------------------
% out.t      = t.';                 % legacy row vector
% out.t_col  = t(:);                % extra column vector (for sim_run convenience)
% 
% out.disp   = TipA;                % legacy name (A-frame translations at selTip)
% out.q0     = q0;   out.q1 = q1;   out.q2 = q2;
% out.qxi    = qxi;  out.chi = chi;
% out.TipA_hat = TipA_hat; 
% out.t_ctr    = t_ctr;
% 
% % ---- Legacy debug plots (exactly as before) ------------------------------
% if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
%     figure; plot(t_ctr,TipA_hat(:,1)*100/cfg.debug.plt_scale);
%     hold on
%     plot(t_ctr,TipA_hat(:,2)*100/cfg.debug.plt_scale);
%     plot(t_ctr,TipA_hat(:,3)*100/cfg.debug.plt_scale);
%     xlabel('t [s]'); ylabel('Tip [% b/2]'); grid on
%     legend({'x','y','z'})
%     title('Wing-tip estimation (A-frame, % span/2)');
%     hold off
% end
% 
% if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'debugPlots') && cfg.plotOpts.debugPlots
%     figure; tiledlayout(3,1,'TileSpacing','tight');
%     nexttile; plot(t,q0'); title('q_0'); grid on
%     nexttile; plot(t,q1'); title('q_1'); grid on
%     nexttile; plot(t,q2'); title('q_2'); grid on
% end
% 
% if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'wingTip') && cfg.plotOpts.wingTip
%     figure; plot(t,TipA(:,3)*100/cfg.debug.plt_scale);
%     if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
%         hold on; plot(t_ctr,TipA_hat(:,3)*100/cfg.debug.plt_scale); hold off
%         legend({'True','Estimate'})
%     end
%     xlabel('t [s]'); ylabel('Tip z [% b/2]'); grid on
%     title('Wing-tip heave (A-frame, % span/2)');
% end
% 
% % -------------------- energies (legacy math kept) -------------------------
% Tkin = 0.5*sum(q1.^2, 1);           % 1 x Nt
% Vpot = 0.5*sum(q2.^2, 1);           % 1 x Nt
% Edot = sum((q1.') * (diag(beam.Sigma)*q1 + cfg.flight.a*cfg.flight.t_inf*out.Faero.'), 1);
% Etot = Tkin + Vpot;
% 
% out.Tkin = Tkin.';   out.Vpot = Vpot.';   out.Etot = Etot.';  % Nt x 1
% 
% if isfield(cfg,'plotOpts') && isfield(cfg.plotOpts,'debugPlots') && cfg.plotOpts.debugPlots
%     figure('Name','Energy check');
%     plot(t, Tkin, 'b-', t, Vpot, 'r--', t, Etot, 'k');
%     grid on; xlabel('t  [s]');
%     ylabel('Energy  [J]');
%     legend('Kinetic','Potential','Total','Location','best');
%     title('Structural modal energy (flexible part)');
% 
%     figure('Name','Energy Rate check');
%     plot(t, Edot, 'b-');
%     grid on; xlabel('t  [s]');
%     ylabel('Energy rate  [J/s]');
%     title('Structural modal energy rate');
% end
% 
% % ====================== EXTRA FIELDS FOR sim_run EXPORTS ==================
% % (These do not change any legacy behavior; they just add convenience.)
% b2 = max(eps, cfg.debug.plt_scale);        % span/2 scaling you use
% out.tip.xyz        = out.disp;              % alias: Nt x 3 (x,y,z in A-frame)
% out.tip.z_norm_pct = 100*out.disp(:,3)/b2;  % Nt x 1
% out.angles.deg     = rad2deg(out.chi.');    % Nt x 3 (roll,pitch,yaw)
% out.energy.Tkin    = out.Tkin;              % aliases as struct
% out.energy.Vpot    = out.Vpot;
% out.energy.Etot    = out.Etot;
% out.q.q0 = out.q0; out.q.q1 = out.q1; out.q.q2 = out.q2;  % nested copy
% 
% end
