clear all; clc;
% close all;
% Will add comments later
cfg = nominalConfig;
% delete(cfg.paths.cache)
% cfg.sim.storeSens = false;           % true if also want ∂x/∂x0
cfg.sim.storeSens = true;           % true if also want ∂x/∂x0
cfg.NetworkPath = 1;
% if NetworkPath 
%     % addpath("\\malexeev\users")
% end
aero = AeroFlex.aero.AeroROM(cfg);
cfg.NetworkPath = 0;

beam = AeroFlex.beam.BeamModel(cfg);

% beam.Sigma = -10*ones(size(beam.Sigma)); % Temp Dampings
% beam.Sigma = -ones(size(beam.Sigma)); % Temp Dampings
base = AeroFlex.core.build_xi_coupling(beam, aero, cfg);

if cfg.trim.do
    % after BeamModel & AeroROM exist
    trim = AeroFlex.sim.trimAircraft(cfg,beam,aero,base);
end

idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
% [L, blk] = AeroFlex.core.assemble_L_matrix(cfg, beam, aero, base);
cfg.sim.storeSens = false;           % true if also want ∂x/∂x0
% cfg.fJac = AeroFlex.sim.buildNonlinearJacobian(beam.Gamma1,beam.Gamma2,base.Gamma_g,base.Gamma_xi, idx);

sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);


% x0 = zeros(size(sim.L,1),1);         % trim initialization
% % x0(2*beam.Nm+1) = 1;
% % x0(1:beam.Nm) = -1e-1;
% % x0(beam.Nm+1:2*beam.Nm) = 1e-1;
% % x0(end-2:end) = [0; .017; 0];
% 
% % state0 = load("state0.mat");
% % x0 = state0.state0;
% 
% q1start0 = aero.DataMatrix.q1start;
% q1_nodal = q1start0(1:size(beam.phi0,1));
% q10 = beam.phi0.'*beam.Mglobal*q1_nodal;
% 
% q0start0 = aero.DataMatrix.q0start;
% q0_nodal = q0start0(1:size(beam.phi0,1));
% q00 = beam.phi0.'*beam.Mglobal*q0_nodal;
% % q10 = beam.phi0.'*q1_nodal;
% x0(1:beam.Nm) = q10;
% 
% q_inf = 0.5*cfg.flight.rho*cfg.flight.U_inf^2;
% a = q_inf*cfg.flight.b_ref^2;
% 
% 
% 
% q02 = -beam.Omega*q00;
% % q02 = aero.DataMatrix.forces_aero_beam_dof;
% % q02 = -beam.phi2.'/beam.Kglobal*aero.DataMatrix.forces_structure(1:size(beam.phi2,1));
% % q02 = -beam.phi0.'*aero.DataMatrix.forces_structure(1:size(beam.phi2,1));
% 
% % x0(beam.Nm+1:2*beam.Nm) = q02;
% 
% % x0(beam.Nm+1:2*beam.Nm) = -q02;
% 
% x0(end-2:end) = q0start0(end-2:end);
% 
% perm = @(g, shep) permute(g, shep);
% shep = [2 1 3];
% % shep = [2 3 1];
% % shep = [3 2 1];
% % Config Gamma Initial States
% gamma0 = aero.DataMatrix.aero_gamma_state.gamma0;
% gamma0 = perm(gamma0, shep);
% gamma0_list = reshape(gamma0, [], 1, 1);
% 
% gamma_n = aero.DataMatrix.aero_gamma_state.gamma;
% gamma_n = perm(gamma_n, shep);
% gamma_list = reshape(gamma_n, [], 1, 1);
% 
% gamma_dot = aero.DataMatrix.aero_gamma_state.gamma_dot;
% gamma_dot = perm(gamma_dot, shep);
% gamma_dot_list = reshape(gamma_dot, [], 1, 1);
% 
% gamma_star = aero.DataMatrix.aero_gamma_state.gamma_star;
% gamma_star = perm(gamma_star, shep);
% gamma_star_list = reshape(gamma_star, [], 1, 1);
% 
% gamma_star0  = aero.DataMatrix.aero_gamma_state.gamma_star0;
% gamma_star0 = perm(gamma_star0, shep);
% gamma_star_list0 = reshape(gamma_star0, [], 1, 1);
% 
% xa_gust_tmp = zeros(5,1); % Assume 0 Gust for now?
% % x_a = [xa_gust_tmp; gamma_list; gamma_star_list; gamma_dot_list; gamma0_list];
% x_a = [xa_gust_tmp; gamma_list; gamma_star_list; cfg.sim.dt*gamma_dot_list; gamma0_list];
% % x_a = [xa_gust_tmp; gamma0_list; gamma_star_list; gamma_dot_list; zeros(size(gamma0_list))];
% % x_a = [xa_gust_tmp; gamma_list; gamma_star_list; gamma_dot_list; zeros(size(gamma0_list))];
% 
% project_V = aero.DataMatrix.aero_gamma_state.project_V;
% if norm(project_V.i) ==0 % no imaginary projectors
%     project_V = project_V.r;
% else
%     project_V = project_V.r.' + project_V.i.';
% end
% 
% q_Gamma = project_V*x_a;
% % x0(sim.idx.qGam)= q_Gamma;
% 
% xi_vec = reshape(base.xi_bar.', [], 1);
% % q_xi = phi_xi.'*xi_vec;
% q_xi = base.phi_xi_modes\xi_vec;
% x0(sim.idx.qxi)= q_xi;


% x0 = zeros(size(sim.L,1),1);         % trim initialization
x0 = trim.states;
u_ctrl = zeros(4,1); % Openloop
[t, x, log, sensEq] = sim.run(x0, cfg.sim.t_end, u_ctrl,cfg.sim.storeSens);   % no gust, no control
%%
% state0 = x(:,end);
% save('state0', "state0")
% plotOpts = cfg.plotOpts;
% EstPts = pts at 0L 0.25L 0.5L 0.75L and L
out = AeroFlex.sim.postProcess(t, x, cfg, beam, aero, base, log);
time = 0:cfg.sim.dt:cfg.sim.t_end;

%%
% % quick plot – modal displacement q₁
% 
q2_hist = x(sim.idx.q2,:);
% q0_hist = -beam.Omega\q2_hist; % Modal displacements
% q1_hist = x(sim.idx.q1,:);
% q_xi_hist = x(sim.idx.qxi,:);
% q_Gamma_hist = x(sim.idx.qGam,:);
% hat_chi_hist = x(sim.idx.chi,:);
% 
% %% Recover Displacemetns (need to use real formulas, currently incorrect)
phi0 = beam.phi0;
% phi1 = beam.phi1;
phi2 = beam.phi2;
% phi_xi = base.phi_xi_modes;
% 
% idx_subxi = size(phi_xi,1);
idx_sub0 = size(phi0, 1);
% 
% Xi_hist = zeros(idx_subxi,length(time));
% X0_hist = zeros(idx_sub0, length(time));
% X1_hist = zeros(idx_sub0, length(time));
X2_hist = zeros(idx_sub0, length(time));
% 
% % for k = 1:size(q_xi_hist,1)
% %     Xi_t = phi_xi(:,k)*(q_xi_hist(k,:));
% %     Xi_hist = Xi_hist+Xi_t;
% % end
% xi_bar = base.xi_bar;
% T_AB = AeroFlex.core.T_phi_quat(xi_bar(1,:));
% UT = AeroFlex.core.halfTangentialOperator( xi_bar(1,:) );
% %%
% T_op = blkdiag(T_AB, UT);
% T_op = repmat({T_op} ,1, beam.nFlex/6 );
% T_op = blkdiag(T_op{:});
% d_X0_hist = zeros(size(T_op, 1), length(time));
% X0rot_hist = zeros(size(T_op, 1), length(time));
% X1rot_hist = zeros(size(T_op, 1), length(time));
% %%
% Xi_hist = phi_xi*q_xi_hist;
for j = 1:beam.Nm
% % Temp to see disp estim.
%     d_x0_t = T_op * phi0(:,j)*q0_hist(j,:);
%     x0_t = phi0(:,j)*q0_hist(j,:);
%     X0_hist = X0_hist+x0_t;
%     d_X0_hist = d_X0_hist +d_x0_t;
% 
% 
% 
% 
%     x1_t = phi1(:,j)*q1_hist(j,:);
%     X1_hist = X1_hist+x1_t;
% 
    x2_t = phi2(:,j)*q2_hist(j,:);
    X2_hist = X2_hist+x2_t;
% 
%     % disp_transl_unrot = phi_omega(:,j)*q1_hist(:,j).';
end
% %%
% for i = 1:length(time)
%     for k = 1: beam.fem.num_node-1
%         idx_quat = double(k-1)*4 +(1:4);
%         % idx_T1 = double(k-1)*6 + (1:6);
%         % idx_T2 = double(k-1)
%         quat_ti = Xi_hist(idx_quat, i);
%         T_AB_t = AeroFlex.core.T_phi_quat(quat_ti);
%         UT_t = AeroFlex.core.halfTangentialOperator( quat_ti );
%         T_opt = blkdiag(T_AB_t, UT_t);
%         T_opt = repmat({T_opt} ,1, beam.nFlex/6 );
%         T_opt = blkdiag(T_opt{:});
% 
%         X0rot_hist(:, i) = T_opt*X0_hist(:, i);
%         X1rot_hist(:, i) = T_opt*X1_hist(:, i);
%     end
% end
% %%
% idx_wing_rot =  52;
% 
% 
scl = .55/2;
% % scl = 1;
% %%
% figure;
% % load(append(save_file, '.mat'))
% % plot(time, X0_hist_openLoop(45,:)/scl*100);
% hold on;
% plot(time, X2_hist(45,:)/scl*100);
% % plot(time, X0_hist(45,:)); % Unscaled
% % plot(time, equib0/scl*100);
% title('Z-Wing Rotated Disp Normalized With Span Wing 9'); xlabel('Time'); ylabel('Amplitude% (Disp/span/2)');
% legend({'open-loop'})
% grid on
% hold off
% %%

%%
figure;
% load(append(save_file, '.mat'))
% plot(time, X0_hist_openLoop(45,:)/scl*100);
hold on;
% plot(time, X2_hist(52,:));
plot(time, X2_hist(52,:));
plot(time, .6*(.24+10*X2_hist(46,:))); % Unscaled
plot(time, .09+X2_hist(53,:)); % Unscaled
% plot(time, X2_hist(59,:)); % Unscaled
% plot(time, equib0/scl*100);
title('Root Bending Moment'); xlabel('Time'); ylabel('N*m');
legend('open-loop', 'PID P = 10', 'nMHE/nMPC')
grid on
hold off
%%
% figure;
% % load(append(save_file, '.mat'))
% % plot(time, X0_hist_openLoop(45,:)/scl*100);
% hold on;
% plot(time, X0_hist(45,:)/scl*100);
% % plot(time, X0_hist(45,:)); % Unscaled
% % plot(time, equib0/scl*100);
% title('Z-Wing Disp Normalized With Span Wing 9'); xlabel('Time'); ylabel('Amplitude% (Disp/span/2)');
% legend({'open-loop'})
% grid on
% hold off
% 
% figure;
% % load(append(save_file, '.mat'))
% % plot(time, X0_hist_openLoop(45,:)/scl*100);
% hold on;
% plot(time, d_X0_hist(idx_wing_rot,:)/scl*100);
% % plot(time, X0_hist(45,:)); % Unscaled
% % plot(time, equib0/scl*100);
% title('dZ-Wing Disp Normalized With Span Wing 9'); xlabel('Time'); ylabel('Amplitude% (Disp/span/2)');
% legend({'open-loop'})
% grid on
% hold off


%% Plot orientation angles
figure;
plot(time, rad2deg(x(end-2:end,:)));

title('Euler Angles');
xlabel('Time [s]');
ylabel('Angle [deg]');
legend('Roll', 'Pitch', 'Yaw');
% plot(t, x(1,:)), grid on
% xlabel('time [s]'), ylabel('q_1')
% title('Open‑loop free decay – mode 1')



% %% 1) flutter analysis ---------------------------------------------------
% flutCfg = struct('Umin',10,'Umax',90,'dU',2, ...   % m/s sweep
%                  'rho',cfg.flight.rho, ...
%                  'savePlots',true);
% flutterAnalysis(beam,aero,base,flutCfg);
% 
% %% 2) frequency‑response (open‑loop FRF) ---------------------------------
% frCfg = struct('U',50, ...              % operating airspeed
%                'wMin',1,'wMax',120,'nPts',1000, ...
%                'outNode','tip', ...     % tip vertical displacement
%                'inType','gust');        % or 'flap'
% frequencyResponsePlot(beam,aero,base,frCfg);
%%%%

% testPos8 = load("struct_pos_node8.dat");
% testPos8(:,4) = testPos8(:,4)/(.55)*100;
% figure(2), hold on, plot(t, testPos8(1:1601,4))
% legend('Built ROM', 'SHARPy Reference')