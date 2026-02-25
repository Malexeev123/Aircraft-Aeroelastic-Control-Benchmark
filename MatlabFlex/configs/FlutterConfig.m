function cfg = FlutterConfig
cfg.case   = 'Flutter';
cfg.dt     = 0.01;
cfg.Nc     = 48;  cfg.Ne = 12;
cfg.controller = nMPC_Artola(cfg);
cfg.estimator  = nMHE_Artola(cfg);
cfg.sensor     = WingVelSensor(loadBeamModel);

Nm = 9;  cfg.Qc = blkdiag(eye(Nm),zeros(Nm));  cfg.Pc = 5*cfg.Qc;
cfg.Rc = 1e-2*eye(2);      % flap + rate penalty
cfg.uL = [deg2rad(-15); -deg2rad(100)];
cfg.uU = [deg2rad( 15);  deg2rad(100)];

cfg.Qe = 0.1*eye(5);   cfg.Pe = blkdiag(zeros(Nm),eye(Nm));
cfg.Re = eye(2)*0.01;
end
