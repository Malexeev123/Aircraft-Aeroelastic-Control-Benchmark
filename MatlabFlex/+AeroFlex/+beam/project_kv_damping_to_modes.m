function Dr = project_kv_damping_to_modes(phi0_mn, Dlin)
% PROJECT_KV_DAMPING_TO_MODES
% Build a node-space matrix Dglob = D0 + D1_sub*S1 + D2_sub*S2,
% then project with mass-normalised φ0 (so Mhat = I).
Dglob = Dlin.D0_node + Dlin.gather( Dlin.D1_sub * Dlin.S1 + Dlin.D2_sub * Dlin.S2 );
Dr = phi0_mn.' * Dglob * phi0_mn;   % reduced (Nm x Nm)
% enforce symmetry for numerical hygiene
Dr = 0.5*(Dr+Dr.');
end
