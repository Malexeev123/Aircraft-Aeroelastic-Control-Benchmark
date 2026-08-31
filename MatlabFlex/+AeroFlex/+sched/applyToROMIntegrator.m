function model = applyToROMIntegrator(model, sched, cfg)
%APPLYTOROMINTEGRATOR Install a scheduled ROM into an existing ROMIntegrator.
%
% The propagated state is already in the Stage-2 common coordinate system.
% This function replaces the active ROM operators, restores runtime
% bookkeeping fields, applies the projection policy, and rebuilds the IMEX
% factorization.

assert(~isempty(model), ...
    'applyToROMIntegrator received an empty ROMIntegrator.');

if nargin < 3 || isempty(cfg)
    cfg = struct();
end

req = {'L','idx','parConst'};
for k = 1:numel(req)
    if ~isfield(sched, req{k})
        error('applyToROMIntegrator:MissingField', ...
            'Scheduled ROM is missing field "%s".', req{k});
    end
end

oldPc = struct();
if isprop(model,'parConst') && isstruct(model.parConst)
    oldPc = model.parConst;
end

model.sched    = sched;
model.L        = sched.L;
model.idx      = sched.idx;
model.parConst = sched.parConst;
if isfield(sched,'internalCoupledCoordinate')
    model.internalCoupledCoordinate=sched.internalCoupledCoordinate;
else
    model.internalCoupledCoordinate=struct('enabled',false);
end

% Preserve runtime-only bookkeeping fields if they existed.
if isfield(oldPc,'u_ctrl')
    model.parConst.u_ctrl = oldPc.u_ctrl;
end

if isfield(oldPc,'gust')
    model.parConst.gust = oldPc.gust;
end

% The atomic scheduled package owns the physical runtime step. Source sample
% times were consumed by discrete-to-continuous conversion and are provenance.
if isfield(model.parConst,'dt') && ~isempty(model.parConst.dt)
    model.dt = model.parConst.dt;
elseif isfield(cfg,'sim') && isfield(cfg.sim,'dt') && ~isempty(cfg.sim.dt)
    model.dt = cfg.sim.dt;
    model.parConst.dt = cfg.sim.dt;
end

% Runtime projection is owned by the receiving model/body case. A trim
% construction option must not silently change wing-only propagation when a
% scheduled package is installed.
useProjection = false;
if isfield(oldPc,'RateProject') && isstruct(oldPc.RateProject) && ...
        isfield(oldPc.RateProject,'projSet') && ...
        ~isempty(oldPc.RateProject.projSet)
    useProjection = logical(oldPc.RateProject.projSet);
end

useProjection = AeroFlex.sched.runtimeRateProjectionPolicy( ...
    cfg,useProjection);

Pz = [];
if isfield(sched,'beam') && isfield(sched.beam,'Pz') && ~isempty(sched.beam.Pz)
    Pz = sched.beam.Pz;
elseif isfield(model.parConst,'RateProject') && ...
        isfield(model.parConst.RateProject,'Pz') && ...
        ~isempty(model.parConst.RateProject.Pz)
    Pz = model.parConst.RateProject.Pz;
end

if useProjection
    if isempty(Pz)
        error('applyToROMIntegrator:MissingPz', ...
            'Rate projection is enabled, but scheduled beam.Pz is missing.');
    end

    model.parConst.RateProject = struct( ...
        'projSet', true, ...
        'Pz', Pz);
else
    model.parConst.RateProject = struct( ...
        'projSet', false, ...
        'Pz', []);
end

% Fuselage/CG thrust is handled by the rigid force balance unless explicitly
% enabled as a modal wing force.
if isfield(cfg,'trim') && isfield(cfg.trim,'thrustActsOnWing') && ...
        ~logical(cfg.trim.thrustActsOnWing)
    model.parConst.N_Thrust = zeros(numel(model.idx.q1),1);
end

model = model.rebuildConstrainedOperator();
end

% function model = applyToROMIntegrator(model, sched, cfg)
% %APPLYTOROMINTEGRATOR Install a scheduled ROM into an existing ROMIntegrator.
% %
% % Stage-2 scheduled libraries are already in common reduced coordinates.
% % Therefore the propagated state is not transformed here.  This function
% % only replaces the active ROM operators, applies the trim/runtime projection
% % policy, and rebuilds the IMEX factorization.
% 
% assert(~isempty(model), ...
%     'applyToROMIntegrator received an empty ROMIntegrator.');
% 
% if nargin < 3 || isempty(cfg)
%     cfg = struct();
% end
% 
% req = {'L','idx','parConst'};
% for k = 1:numel(req)
%     if ~isfield(sched, req{k})
%         error('applyToROMIntegrator:MissingField', ...
%             'Scheduled ROM is missing field "%s".', req{k});
%     end
% end
% 
% model.sched    = sched;
% model.L        = sched.L;
% model.idx      = sched.idx;
% model.parConst = sched.parConst;
% 
% % Use the MATLAB simulation step.  The SHARPy saved Ts is only the source
% % discretization used before continuous-time conversion.
% if isfield(cfg,'sim') && isfield(cfg.sim,'dt') && ~isempty(cfg.sim.dt)
%     model.dt = cfg.sim.dt;
%     model.parConst.dt = cfg.sim.dt;
% elseif isfield(model.parConst,'dt') && ~isempty(model.parConst.dt)
%     model.dt = model.parConst.dt;
% end
% 
% useProjection = true;
% if isfield(cfg,'trim') && isfield(cfg.trim,'useRateProjection') && ~isempty(cfg.trim.useRateProjection)
%     useProjection = logical(cfg.trim.useRateProjection);
% end
% 
% Pz = [];
% if isfield(sched,'beam') && isfield(sched.beam,'Pz') && ~isempty(sched.beam.Pz)
%     Pz = sched.beam.Pz;
% elseif isfield(model.parConst,'RateProject') && isfield(model.parConst.RateProject,'Pz')
%     Pz = model.parConst.RateProject.Pz;
% end
% 
% if useProjection
%     if isempty(Pz)
%         error('applyToROMIntegrator:MissingPz', ...
%             'Rate projection is enabled, but scheduled beam.Pz is missing.');
%     end
% 
%     model.parConst.RateProject = struct( ...
%         'projSet', true, ...
%         'Pz', Pz);
% else
%     model.parConst.RateProject = struct( ...
%         'projSet', false, ...
%         'Pz', []);
% end
% 
% % Fuselage/CG thrust is handled by the rigid force balance unless explicitly
% % enabled as a modal wing force.
% if isfield(cfg,'trim') && isfield(cfg.trim,'thrustActsOnWing') && ...
%         ~logical(cfg.trim.thrustActsOnWing)
%     model.parConst.N_Thrust = zeros(numel(model.idx.q1),1);
% end
% 
% model = model.rebuildConstrainedOperator();
% end
