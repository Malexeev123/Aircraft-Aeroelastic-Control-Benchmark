function plant = applyToPlant(plant, sched)
%APPLYTOPLANT Apply a schedule to PlantRunTime without rebuilding from SHARPy.
plant.model = AeroFlex.sched.applyToROMIntegrator(plant.model, sched, plant.cfg);
plant.idx = sched.idx;
plant.base.Gamma_xi = sched.base.Gamma_xi;
plant.base.Gamma_g  = sched.base.Gamma_g;
if isfield(sched.base,'xi_bar')
    plant.base.xi_bar = sched.base.xi_bar;
end
plant.last.sched_mu = sched.mu;
plant.last.sched_weights = sched.weights;
plant.last.sched_pointIds = sched.pointIds;
plant.last.sched_info = sched.info;
end
