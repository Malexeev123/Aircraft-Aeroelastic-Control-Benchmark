clear; clc;
rgParam = RigidBody.methods.paramsRigid_PazyUAV();
cfg = nominalConfig;

rigidPlant = RigidBody.methods.EOM(cfg, rgParam);
R6EoM = RigidBody.methods.EOM(cfg, rgParam);
uExt.F = zeros(3,1); 
% uExt.F = [0;0;-10]; % negative Z in up direction
uExt.M = zeros(3,1); 
uExt.uB = zeros(3,1);
% uExt.uB = [0;0;-15];
x = [1;0;0;0; ...
    10;0;0; ... 
    0;0;0; ...
    0;0;0];
dt =cfg.sim.dt ;
Ts = 10;
time = 0:dt:Ts;

xhist = zeros(length(x),length(time));
eulHist = zeros(3,length(time));
rHist = zeros(3,length(time));
for i = 1:length(time)
    x_n = rigidPlant.step(x,uExt,dt);
    xhist(:,i) = x_n;
    x = x_n;


    q= x(1:4).';
    eul = quat2eul(q);
    eulHist(:, i) = eul.';
    rHist(:,i) = x(11:13);
end
%% Plots
% ------------------ Euler -------------------------------
figure(1); clf
subplot(3,1,1); plot(time,rad2deg(eulHist(1,:)));  grid on
ylabel('\theta  [deg]'); title('Rigid–body pitch response')

subplot(3,1,2); plot(time,rad2deg(eulHist(2,:)));  grid on
ylabel('\phi  [deg]');

subplot(3,1,3); plot(time,rad2deg(eulHist(3,:)));  grid on
ylabel('\psi  [deg]'); xlabel('t  [s]');
% --------------------- Velocities ------------------------------

figure(2); clf
subplot(2,1,1); plot(time,xhist(5,:));  hold on;
plot(time,xhist(6,:)); plot(time,xhist(7,:)); hold off; grid on
ylabel('Vel.  [m/s]'); title('Rigid–body velocity response'); legend('Vx', ...
    'Vy', 'Vz')

subplot(2,1,2); plot(time,xhist(8,:));  hold on;
plot(time,xhist(9,:)); plot(time,xhist(10,:)); hold off; grid on
ylabel('Angular Vel  [rad/s]'); title('Rigid–body angular velocity response')
legend('Wx', 'Wy', 'Wz')

% --------------------- Position (r) ------------------------------
figure(3); clf
plot(time,xhist(11,:)); hold on;
plot(time,xhist(12,:)); plot(time,xhist(13,:)); hold off; grid on;
ylabel('Position  [m]'); title('Rigid–body position response'); legend('rx', ...
    'ry', 'rz')
