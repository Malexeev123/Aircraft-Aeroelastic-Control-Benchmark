function plant = applyToPlant(plant, sched)
%APPLYTOPLANT Apply a Stage-2 scheduled ROM to PlantRunTime.
%
% Stage-2 compatible libraries use common reduced coordinates.  The plant
% state is therefore not transformed here.  This function only commits an
% already-accepted scheduled ROM package.

plant.sched = sched;

plant.model = AeroFlex.sched.applyToROMIntegrator( ...
    plant.model, sched, plant.cfg);
plant.dt = plant.model.dt;
plant.cfg.sim.dt = plant.model.dt;

% Load recovery and rate projection must use the active scheduled/reference
% maps.  The beam maps should be identical across a finalized Stage-2
% library, but we still commit the active maps so diagnostics and recovery
% use the same object.
if isfield(sched,'beam')
    if isfield(sched.beam,'Pz'),  plant.beam.Pz  = sched.beam.Pz;  end
    if isfield(sched.beam,'Pr'),  plant.beam.Pr  = sched.beam.Pr;  end
    if isfield(sched.beam,'red'), plant.beam.red = sched.beam.red; end
end

plant.idx = sched.idx;

if isfield(sched,'base')
    plant.base = sched.base;
end

% Keep diagnostic force-map fields synchronized. Propagation uses
% plant.model.parConst, but several diagnostics inspect plant.aero.forceMap.
if isfield(sched,'parConst') && isfield(plant,'aero') && isfield(plant.aero,'forceMap')
    plant.aero.forceMap.Bw       = sched.parConst.Bw;
    plant.aero.forceMap.Dw       = sched.parConst.Dw;
    plant.aero.forceMap.B_delta  = sched.parConst.Bdel;
    plant.aero.forceMap.D_delta  = sched.parConst.Ddel;
    plant.aero.forceMap.B_ddelta = sched.parConst.Bddel;
    plant.aero.forceMap.D_ddelta = sched.parConst.Dddel;
end

plant.nx = size(plant.model.L,1);

if isfield(plant.model.parConst,'Bdel') && isfield(plant.model.parConst,'Bddel')
    plant.nu = size(plant.model.parConst.Bdel,2) + ...
               size(plant.model.parConst.Bddel,2);
end

if isfield(plant.model.parConst,'Bw')
    plant.nw = size(plant.model.parConst.Bw,2);
end

if numel(plant.xFlex) ~= plant.nx
    error('applyToPlant:StateSize', ...
        'xFlex length is %d, but scheduled ROM has nx=%d.', ...
        numel(plant.xFlex), plant.nx);
end

plant.last.sched_mu       = sched.mu;
plant.last.sched_weights  = sched.weights;
plant.last.sched_pointIds = sched.pointIds;
plant.last.sched_info     = sched.info;
plant.last.qRatio         = 1;

if isfield(sched,'x_eq')
    plant.last.sched_x_eq = sched.x_eq;
end

if isfield(sched,'u_eq')
    plant.last.sched_u_eq = sched.u_eq;
end
end

% function plant = applyToPlant(plant, sched)
% %APPLYTOPLANT Apply a Stage-2 scheduled ROM to PlantRunTime.
% %
% % The flexible state remains in the common Stage-2 coordinates.  Runtime
% % scheduling only replaces L/parConst/recovery maps and rebuilds the IMEX
% % factors through applyToROMIntegrator.
% 
% plant.sched = sched;
% 
% plant.model = AeroFlex.sched.applyToROMIntegrator( ...
%     plant.model, sched, plant.cfg);
% 
% % Load recovery and rate projection must use the active scheduled/reference
% % maps.  Keeping the pre-schedule maps here makes trim replay look correct
% % but lets propagation use the wrong root constraints.
% if isfield(sched,'beam')
%     if isfield(sched.beam,'Pz'),  plant.beam.Pz  = sched.beam.Pz;  end
%     if isfield(sched.beam,'Pr'),  plant.beam.Pr  = sched.beam.Pr;  end
%     if isfield(sched.beam,'red'), plant.beam.red = sched.beam.red; end
% end
% 
% plant.idx = sched.idx;
% 
% if isfield(sched,'base')
%     plant.base = sched.base;
% end
% 
% % Keep diagnostic force-map fields synchronized.  Propagation itself uses
% % plant.model.parConst.
% if isfield(sched,'parConst') && isfield(plant,'aero') && isfield(plant.aero,'forceMap')
%     plant.aero.forceMap.Bw       = sched.parConst.Bw;
%     plant.aero.forceMap.Dw       = sched.parConst.Dw;
%     plant.aero.forceMap.B_delta  = sched.parConst.Bdel;
%     plant.aero.forceMap.D_delta  = sched.parConst.Ddel;
%     plant.aero.forceMap.B_ddelta = sched.parConst.Bddel;
%     plant.aero.forceMap.D_ddelta = sched.parConst.Dddel;
% end
% 
% plant.nx = size(plant.model.L,1);
% if isfield(plant.model.parConst,'Bdel') && isfield(plant.model.parConst,'Bddel')
%     plant.nu = size(plant.model.parConst.Bdel,2) + ...
%                size(plant.model.parConst.Bddel,2);
% end
% if isfield(plant.model.parConst,'Bw')
%     plant.nw = size(plant.model.parConst.Bw,2);
% end
% 
% if numel(plant.xFlex) ~= plant.nx
%     error('applyToPlant:StateSize', ...
%         'xFlex length is %d, but scheduled ROM has nx=%d.', ...
%         numel(plant.xFlex), plant.nx);
% end
% 
% plant.last.sched_mu       = sched.mu;
% plant.last.sched_weights  = sched.weights;
% plant.last.sched_pointIds = sched.pointIds;
% plant.last.sched_info     = sched.info;
% plant.last.qRatio         = 1;
% 
% if isfield(plant.cfg,'debug') && isfield(plant.cfg.debug,'schedule') && plant.cfg.debug.schedule
%     plant.debugScheduleConsistency('applyToPlant', sched);
% end
% end
% 
% 
% 
% % 
% % function plant = applyToPlant(plant, sched)
% % %APPLYTOPLANT Apply a Stage-2 scheduled ROM to PlantRunTime.
% % %
% % % Stage-2 libraries are already in common reduced coordinates.  The plant
% % % state is not transformed here.  Only the active ROM data and recovery maps
% % % are replaced.
% % 
% %     plant.sched = sched;
% % 
% %     plant.model = AeroFlex.sched.applyToROMIntegrator( ...
% %         plant.model, sched, plant.cfg);
% % 
% %     % Load recovery and projection must use the scheduled/reference maps.
% %     if isfield(sched,'beam')
% %         if isfield(sched.beam,'Pz')
% %             plant.beam.Pz = sched.beam.Pz;
% %         end
% %         if isfield(sched.beam,'Pr')
% %             plant.beam.Pr = sched.beam.Pr;
% %         end
% %         if isfield(sched.beam,'red')
% %             plant.beam.red = sched.beam.red;
% %         end
% %     end
% % 
% %     plant.idx = sched.idx;
% % 
% %     if isfield(sched,'base')
% %         plant.base = sched.base;
% %     end
% % 
% %     % Keep diagnostic force-map fields synchronized. Propagation uses
% %     % plant.model.parConst, but several checks inspect plant.aero.forceMap.
% %     if isfield(sched,'parConst') && isfield(plant,'aero') && isfield(plant.aero,'forceMap')
% %         plant.aero.forceMap.Bw       = sched.parConst.Bw;
% %         plant.aero.forceMap.Dw       = sched.parConst.Dw;
% %         plant.aero.forceMap.B_delta  = sched.parConst.Bdel;
% %         plant.aero.forceMap.D_delta  = sched.parConst.Ddel;
% %         plant.aero.forceMap.B_ddelta = sched.parConst.Bddel;
% %         plant.aero.forceMap.D_ddelta = sched.parConst.Dddel;
% %     end
% % 
% %     plant.nx = size(plant.model.L,1);
% %     plant.nu = size(plant.model.parConst.Bdel,2) + ...
% %                size(plant.model.parConst.Bddel,2);
% %     plant.nw = size(plant.model.parConst.Bw,2);
% % 
% %     if numel(plant.xFlex) ~= plant.nx
% %         error('applyToPlant:StateSize', ...
% %             'xFlex length is %d, but scheduled ROM has nx=%d.', ...
% %             numel(plant.xFlex), plant.nx);
% %     end
% % 
% %     plant.last.sched_mu       = sched.mu;
% %     plant.last.sched_weights  = sched.weights;
% %     plant.last.sched_pointIds = sched.pointIds;
% %     plant.last.sched_info     = sched.info;
% %     plant.last.qRatio         = 1;
% % 
% %     if isfield(plant.cfg,'debug') && isfield(plant.cfg.debug,'schedule') && plant.cfg.debug.schedule
% %         plant.debugScheduleConsistency('applyToPlant',sched);
% %     end
% % end
% % % function plant = applyToPlant(plant, sched)
% % % %APPLYTOPLANT Apply a scheduled Stage-2 ROM to PlantRunTime.
% % % %
% % % % The flexible state is not transformed here.  Stage-2 compatibility already
% % % % puts all library vertices into one reduced coordinate system.
% % % 
% % % plant.model = AeroFlex.sched.applyToROMIntegrator(plant.model, sched, plant.cfg);
% % % 
% % % % Runtime load recovery and projection must use the scheduled/reference beam
% % % % maps.  Do not keep the pre-scheduling beam maps here.
% % % if isfield(sched,'beam')
% % %     if isfield(sched.beam,'Pz')
% % %         plant.beam.Pz = sched.beam.Pz;
% % %     end
% % %     if isfield(sched.beam,'Pr')
% % %         plant.beam.Pr = sched.beam.Pr;
% % %     end
% % %     if isfield(sched.beam,'red')
% % %         plant.beam.red = sched.beam.red;
% % %     end
% % % end
% % % 
% % % plant.idx = sched.idx;
% % % 
% % % if isfield(sched,'base')
% % %     plant.base = sched.base;
% % % end
% % % 
% % % % Keep force-map fields available for debugging/dimension checks.  The
% % % % actual propagation uses plant.model.parConst.
% % % if isfield(plant,'aero') && isfield(sched,'parConst')
% % %     plant.aero.forceMap.Bw       = sched.parConst.Bw;
% % %     plant.aero.forceMap.Dw       = sched.parConst.Dw;
% % %     plant.aero.forceMap.B_delta  = sched.parConst.Bdel;
% % %     plant.aero.forceMap.D_delta  = sched.parConst.Ddel;
% % %     plant.aero.forceMap.B_ddelta = sched.parConst.Bddel;
% % %     plant.aero.forceMap.D_ddelta = sched.parConst.Dddel;
% % % end
% % % 
% % % plant.last.sched_mu       = sched.mu;
% % % plant.last.sched_weights  = sched.weights;
% % % plant.last.sched_pointIds = sched.pointIds;
% % % plant.last.sched_info     = sched.info;
% % % 
% % % if isfield(plant.cfg,'debug') && isfield(plant.cfg.debug,'schedule') && plant.cfg.debug.schedule
% % %     plant.debugScheduleConsistency('applyToPlant',sched);
% % % end
% % % end
% % % % function plant = applyToPlant(plant, sched)
% % % % %APPLYTOPLANT Apply a scheduled ROM to PlantRunTime.
% % % % 
% % % %     plant.model = AeroFlex.sched.applyToROMIntegrator(plant.model, sched, plant.cfg);
% % % %     % plant.beam.Pz = sched.beam.Pz;
% % % %     % plant.beam.Pr = sched.beam.Pr;
% % % %     % plant.beam.red = sched.beam.red;
% % % %     % Pz/Pr/red stay in the plant beam basis.
% % % %     plant.model.parConst.RateProject = struct( ...
% % % %         'projSet', true, ...
% % % %         'Pz', plant.beam.Pz);
% % % % 
% % % %     if ismethod(plant.model,'rebuildConstrainedOperator')
% % % %         plant.model = plant.model.rebuildConstrainedOperator();
% % % %     end
% % % % 
% % % %     plant.idx = sched.idx;
% % % % 
% % % %     plant.base.Gamma_xi = sched.base.Gamma_xi;
% % % %     plant.base.Gamma_g  = sched.base.Gamma_g;
% % % % 
% % % %     if isfield(sched.base,'xi_bar')
% % % %         plant.base.xi_bar = sched.base.xi_bar;
% % % %     end
% % % % 
% % % %     % Note: Do not overwrite plant.beam.Pz/Pr/red in Stage 1.
% % % %     % These maps are basis-sensitive.
% % % % 
% % % %     plant.last.sched_mu       = sched.mu;
% % % %     plant.last.sched_weights  = sched.weights;
% % % %     plant.last.sched_pointIds = sched.pointIds;
% % % %     plant.last.sched_info     = sched.info;
% % % % 
% % % %     if isfield(plant.cfg,'debug') && isfield(plant.cfg.debug,'schedule') && plant.cfg.debug.schedule
% % % %         plant.debugScheduleConsistency('applyToPlant',sched);
% % % %     end
% % % % end
