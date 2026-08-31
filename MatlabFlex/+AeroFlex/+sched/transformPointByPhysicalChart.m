function point = transformPointByPhysicalChart(source,transform)
%TRANSFORMPOINTBYPHYSICALCHART Transform one package into a physical chart.

point = source;
T = transform.T;
Ti = transform.Tinv;
T1 = transform.q1.T; S1 = transform.q1.Tinv;
T2 = transform.q2.T; S2 = transform.q2.Tinv;
Tx = transform.qxi.T; Sx = transform.qxi.Tinv;
Tg = transform.qGam.T; Sg = transform.qGam.Tinv;
point.L = T*source.L*Ti;
point.x_eq = T*source.x_eq(:);
if isfield(source,'physicalChartEquilibrium')
    point.physicalChartEquilibrium = source.physicalChartEquilibrium;
    point.physicalChartEquilibrium.state = ...
        T*source.physicalChartEquilibrium.state(:);
end
point.parConst = AeroFlex.sched.transformParConstToReference( ...
    source.parConst,transform,source.idx);
if isfield(source.parConst,'affineOffset') && ...
        ~isempty(source.parConst.affineOffset)
    point.parConst.affineOffset = T*source.parConst.affineOffset(:);
end

beam = source.beam;
beam.Pz = T1*beam.Pz*S1;
beam.Pr = T1*beam.Pr*S1;
beam.eta_e = T1*beam.eta_e(:);
beam.Gamma1 = AeroFlex.sched.transformTensor3( ...
    beam.Gamma1,T1,S1,S1);
beam.Gamma2 = AeroFlex.sched.transformTensor3( ...
    beam.Gamma2,T1,S2,S2);
red = beam.red;
red.phi1_sA = red.phi1_sA*S1;
red.phi2_sA = red.phi2_sA*S2;
disc = red.ModeVars_discrete;
disc.phi0_local = disc.phi0_local*S1;
disc.phi1_local = disc.phi1_local*S1;
disc.phi2_local = disc.phi2_local*S2;
if isfield(disc,'psi1_local'), disc.psi1_local = disc.psi1_local*T1.'; end
if isfield(disc,'psi2_local'), disc.psi2_local = disc.psi2_local*T2.'; end
disc.etaVecExt = T1*disc.etaVecExt(:);
disc.Gamma1 = beam.Gamma1;
disc.Gamma2 = beam.Gamma2;
red.ModeVars_discrete = disc;
beam.red = red;
point.beam = beam;

base = source.base;
base.Gamma_g = AeroFlex.sched.transformTensor3( ...
    base.Gamma_g,T1,Sx,Sx);
base.Gamma_xi = AeroFlex.sched.transformTensor3( ...
    base.Gamma_xi,Tx,Sx,S1);
base.phi_xi_modes = base.phi_xi_modes*Sx;
base.phiXi_sA = base.phiXi_sA*Sx;
if isfield(base,'FM')
    if isfield(base.FM,'Dchi'), base.FM.Dchi = T1*base.FM.Dchi; end
    if isfield(base.FM,'Bchi'), base.FM.Bchi = Tg*base.FM.Bchi; end
    if isfield(base.FM,'X'), base.FM.X = base.FM.X*S1; end
end
point.base = base;

fm = source.aero.forceMap;
fm.A_Gamma = Tg*fm.A_Gamma*Sg;
fm.B0 = Tg*fm.B0*S1;
fm.B1 = Tg*fm.B1*S1;
fm.C_Gamma = T1*fm.C_Gamma*Sg;
fm.D0 = T1*fm.D0*S1;
fm.D1 = T1*fm.D1*S1;
fm.Bw = Tg*fm.Bw;
fm.Dw = T1*fm.Dw;
fm.B_delta = Tg*fm.B_delta;
fm.D_delta = T1*fm.D_delta;
fm.B_ddelta = Tg*fm.B_ddelta;
fm.D_ddelta = T1*fm.D_ddelta;
point.aero = source.aero;
point.aero.forceMap = fm;
point.aero.forces_aero_beam_dof = ...
    T1*source.aero.forces_aero_beam_dof(:);

if isfield(point,'p5') && isfield(point.p5,'r2')
    r2 = point.p5.r2;
    r2.rootC = r2.rootC*Sg;
    directTransform = blkdiag(S1,S1,1,eye(4));
    r2.rootD = r2.rootD*directTransform;
    r2.Omega = T2*r2.Omega*S1;
    point.p5.r2 = r2;
end
point.compatibleCoordinates = true;
point.physicalChartTransform = transform;
end
