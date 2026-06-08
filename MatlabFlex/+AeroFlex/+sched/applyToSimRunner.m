function sim = applyToSimRunner(sim, sched, cfg)
%APPLYTOSIMRUNNER Install a scheduled ROM into an existing SimRunner.
sim.L = sched.L;
sim.idx = sched.idx;
sim.parConst = sched.parConst;
if isprop(sim,'Pz') && isfield(sched.beam,'Pz')
    sim.Pz = sched.beam.Pz;
end
fac = AeroFlex.sched.factorIMEX(sim.L, cfg.sim.dt);
sim.Lfac = fac.Lfac;
sim.Ufac = fac.Ufac;
sim.piv = fac.piv;
end
