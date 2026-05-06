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

Nt = numel(t);
Nm = beam.Nm;
if isfield(cfg,'Na') && ~isempty(cfg.Na)
    Na2 = cfg.Na;
else
    Na2 = size(X,1) - (4*Nm + 1 + 3);
    if Na2 < 0
        error('postProcess: cannot infer Na from X size. Provide cfg.Na.');
    end
end

% Keep legacy behavior: override X from cached log in non-openLoop
if ~(isfield(cfg,'case') && strcmpi(cfg.case,'openLoop'))
    try
        openL_log = load(fullfile(fullfile(cfg.paths.cache, cfg.case_name), 'log'),'log');
        if isfield(openL_log,'log') && isfield(openL_log.log,'x') && ~isempty(openL_log.log.x)
            X  = openL_log.log.x;
            Nt = size(X,2);
        end
    catch
    end
end

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

end % postProcess


% ============================== helpers ===================================
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
if ishghandle(fh), close(fh); end

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
