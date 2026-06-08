function parStar = transformParConstBlocks(par, Tblocks)
%TRANSFORMPARCONSTBLOCKS Transform nonlinear tensors under state coordinates.
%
% This is an optional Stage-2 tool.  The first production implementation
% should keep q1, q2 and qxi coordinates fixed and only schedule the Gamma
% tensors directly.  Use this only if q1/q2/qxi bases are transformed.
%
% Required fields in Tblocks:
%   Tq1, Tq2, Txi   maps local coordinates to common coordinates.

parStar = par;
Tq1 = Tblocks.Tq1; Tq1i = inv(Tq1);
Tq2 = Tblocks.Tq2; Tq2i = inv(Tq2);
Txi = Tblocks.Txi; Txii = inv(Txi);

parStar.Gamma1   = AeroFlex.sched.transformTensor3(par.Gamma1,   Tq1, Tq1i, Tq1i);
parStar.Gamma2   = AeroFlex.sched.transformTensor3(par.Gamma2,   Tq1, Tq2i, Tq2i);
parStar.Gamma_g  = AeroFlex.sched.transformTensor3(par.Gamma_g,  Tq1, Txii, Txii);
parStar.Gamma_xi = AeroFlex.sched.transformTensor3(par.Gamma_xi, Txi, Txii, Tq1i);

parStar.forces_0 = Tq1 * par.forces_0;
parStar.N_Thrust = Tq1 * par.N_Thrust;
parStar.Dw       = Tq1 * par.Dw;
parStar.Ddel     = Tq1 * par.Ddel;
parStar.Dddel    = Tq1 * par.Dddel;
end
