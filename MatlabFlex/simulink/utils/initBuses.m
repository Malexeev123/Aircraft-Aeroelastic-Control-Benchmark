% src/matlab/initBuses.m
function initBuses()
% -------- state bus: (rigid-body + flex-wing gen. coords) --------------
elems(1) = Simulink.BusElement; elems(1).Name = 'quatBody';   elems(1).Dimensions = 4;
elems(2) = Simulink.BusElement; elems(2).Name = 'omegaBody';  elems(2).Dimensions = 3;
elems(3) = Simulink.BusElement; elems(3).Name = 'uvwBody';    elems(3).Dimensions = 3;
elems(4) = Simulink.BusElement; elems(4).Name = 'etaWing';    elems(4).Dimensions = AeroFlex.cfg.Nm;   % modal coords
STATE_BUS = Simulink.Bus; STATE_BUS.Elements = elems; save('data/STATE_BUS.mat','STATE_BUS');

% -------- control bus: four aileron commands ---------------------------
cu(1) = Simulink.BusElement; cu(1).Name = 'deltaFlap';  cu(1).Dimensions = 4;
CTRL_BUS = Simulink.Bus; CTRL_BUS.Elements = cu;        save('data/CTRL_BUS.mat','CTRL_BUS');
end
