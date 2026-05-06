function result = smoketest_hinged_damped_beam()
% Simply supported (hinged–hinged) beam, strong flexibility, modest KV.
L=1.0; Nn=41;
E=3.5e9; A=1.0e-4; Iy=1.0e-8; Iz=1.0e-8; GJ=0.5; rho=1600;
eta_f=2e3; eta_m=2e-3;

[fem, keepDofs, chain, elemLen, Eall] = AeroFlex.beam.make_uniform_beam(L,Nn,'hinged',E,A,Iy,Iz,GJ,rho);
[M,K,Vr,w] = AeroFlex.beam.linear_modal_analysis(fem, keepDofs);
Mred = Vr.'*M*Vr; scl = 1./sqrt(diag(Mred)); phi0_mn = Vr.*scl.'; omega = w(:);

opts = struct('eta_f',eta_f,'eta_m',eta_m);
Dlin = AeroFlex.beam.build_linearized_kv_damping(fem, keepDofs, chain, elemLen, Eall, opts);
Dr   = AeroFlex.beam.project_kv_damping_to_modes(phi0_mn, Dlin);

% report first 3 modes’ ζ
zeta = diag(Dr)./(2*omega);
result = table((1:min(3,length(omega))).', omega(1:min(3,end)), zeta(1:min(3,end)), ...
               'VariableNames', {'mode','omega','zeta'});
disp(result)
end
