clear all; close all; clc;

% Will add comments later
% cfg = nominalConfig;
% delete(cfg.paths.cache)
cfg = nominalConfig;
%%
% cfg.NetworkPath = 1;
% if NetworkPath 
%     % addpath("\\malexeev\users")
% end
%%
cfg.sim.storeSens = false;           % true if also want ∂x/∂x0

aero = AeroFlex.aero.AeroROM(cfg);
beam = AeroFlex.beam.BeamModel(cfg);

% beam.Sigma = -10*ones(size(beam.Sigma)); % Temp Dampings
% beam.Sigma = -ones(size(beam.Sigma)); % Temp Dampings
base = AeroFlex.core.build_xi_coupling(beam, aero, cfg);

if cfg.trim.do
    % after BeamModel & AeroROM exist
    trim = AeroFlex.sim.trimAircraft(cfg,beam,aero,base);
end
%% First Doing OpenLoop Case to get sensor dat for tuning of Estimator

idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
% [L, blk] = AeroFlex.core.assemble_L_matrix(cfg, beam, aero, base);
cfg.sim.storeSens = false;           % true if also want ∂x/∂x0
cfg.trim.do = false;
% cfg.fJac = AeroFlex.sim.buildNonlinearJacobian(beam.Gamma1,beam.Gamma2,base.Gamma_g,base.Gamma_xi, idx);

sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);

x0 = trim.states;   % trim initialization

[t, x, log, sensEq] = sim.run(x0, cfg.sim.t_end, cfg.sim.storeSens);   % no gust, no control
%% CASE CONTROLLER 
cfg = GLAConfig(cfg);

% [L, blk] = AeroFlex.core.assemble_L_matrix(cfg, beam, aero, base);

% cfg.SwTest = 1e-3;   % For testing the gradient of gust profile
sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);
sim.debug = true;       % debug plots and trackers


x0 = trim.states;
% %% Sw test to be ran only once
% cfg.SwTest = 1e-3;   % For testing the gradient of gust profile
% sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);
% sim.debug = false;       % debug plots and trackers
% [t1, x1, log1, sensEq1] = sim.run(x0, cfg.sim.dt, cfg.sim.storeSens);   % no gust, no control
% 
% 
% cfg.SwTest = -1e-3;   % For testing the gradient of gust profile
% sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);
% sim.debug = false;       % debug plots and trackers
% [t2, x2, log2, sensEq2] = sim.run(x0, cfg.sim.dt, cfg.sim.storeSens);   % no gust, no control
% 
% Sw_num = (x1(:,end) - x2(:,end))/(2*1e-3);
% 
% 
% cfg.SwTest = 0; 
% sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base, trim);
% sim.debug = false;       % debug plots and trackers
% [~, x_test, log_test, sensEq_test] = sim.run(x0, cfg.sim.dt, cfg.sim.storeSens);   % no gust, no control
% Sw_buf = sensEq_test(:,length(x0)+1,end);
% fprintf('‖Sw_true‖ %.2e   ‖Sw_num‖ %.2e   rel.err %.2e\n',...
%             norm(Sw_buf), norm(Sw_num), ...
%             norm(Sw_buf-Sw_num)/norm(Sw_num));
% %%

[t, x, log, sensEq] = sim.run(x0, cfg.sim.t_end, cfg.sim.storeSens);   % no gust, no control
% [t, x, log, sensEq] = sim.run(x0, cfg.ctrl.Ts, cfg.sim.storeSens);   % no gust, no control
%%
out = AeroFlex.sim.postProcess(t, x, cfg, beam, aero, base, log);
time = 0:cfg.sim.dt:cfg.sim.t_end;

%% Control Plots
time_ctrl = linspace(0, cfg.sim.t_end, size(log.xhat,2));

delta_init = 0;
figure; 
subplot(2,1,1);
plot(time_ctrl, (log.U(1,1:length(time_ctrl)))*180/pi);
hold on;
plot(time_ctrl, (log.U(2,1:length(time_ctrl)))*180/pi);
title('Control Deflection Log'); xlabel('Time'); ylabel('deg ');
legend({'left wing','right wing'})
hold off

subplot(2,1,2);
plot(time_ctrl, (log.U(3,1:length(time_ctrl)))*180/pi);
hold on;
plot(time_ctrl, (log.U(4,1:length(time_ctrl)))*180/pi);
title('Control Rate Log'); 
xlabel('Time'); ylabel('deg/s ');
legend({'left wing','right wing'})
hold off

%% wEst vs wReal
time_ctrl = linspace(0, cfg.sim.t_end, size(log.xhat,2));

figure;
plot(time_ctrl, log.wEst); 
hold on;
plot(t, [log.wTrue, 0]); 
title('Gust Estimated Vs Real'); 
xlabel('Time'); ylabel('gust ');
legend({'Estimated','True'})
hold off;


%% Plot orientation angles
% figure;
% plot(time, rad2deg(x(end-2:end,:)));
% 
% title('Euler Angles');
% xlabel('Time [s]');
% ylabel('Angle [deg]');
% legend('Roll', 'Pitch', 'Yaw');


% plot(t, x(1,:)), grid on
% xlabel('time [s]'), ylabel('q_1')
% % title('Open‑loop free decay – mode 1')
% 
% figure
% plot(time_ctrl , log.wTrue      , 'k--' , 'DisplayName','true gust')
% hold on
% plot(time_ctrl , log.wEst       , 'r-'  , 'DisplayName','estimated')
% xlabel('time [s]'), ylabel('w_g [m/s]')
% legend, grid on

% figure
% plot(log.t , log.x      , 'k--' , 'DisplayName','true state')
% hold on
% plot(log.t , log.xhat       , 'r-'  , 'DisplayName','estimated')
% xlabel('time [s]'), ylabel('w_g [m/s]')
% legend, grid on