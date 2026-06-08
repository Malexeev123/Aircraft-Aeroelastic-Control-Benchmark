function model = applyToROMIntegrator(model, sched, cfg)
%APPLYTOROMINTEGRATOR Install a scheduled ROM into an existing integrator.
model.L = sched.L;
model.idx = sched.idx;
model.parConst = sched.parConst;
if isfield(model.parConst,'dt')
    model.dt = model.parConst.dt;
elseif isfield(cfg,'sim') && isfield(cfg.sim,'dt')
    model.dt = cfg.sim.dt;
end
fac = AeroFlex.sched.factorIMEX(model.L, model.dt);
model.gamma = fac.gamma;
model.delta = fac.delta;
model.Lfac = fac.Lfac;
model.Ufac = fac.Ufac;
model.piv = fac.piv;
end
