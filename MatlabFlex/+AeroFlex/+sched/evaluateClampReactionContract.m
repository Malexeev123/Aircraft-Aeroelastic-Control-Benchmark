function contract = evaluateClampReactionContract(sched,x,u,loadLedger)
%EVALUATECLAMPREACTIONCONTRACT Diagnose independent root-reaction constructions.
%   C1 evaluates the stored intrinsic force-state interpolation. C2 recovers
%   the constraint reaction from the raw, unprojected q1 equation and splits
%   its owned terms. C3 uses the independent rigid-body load balance when a
%   trim/runtime load ledger is supplied. All wrenches use the active body
%   frame, wing-root reference, full-wing multiplicity, and N/N m units.

arguments
    sched (1,1) struct
    x double
    u double
    loadLedger struct = struct()
end

x = x(:);
u = u(:);
idx = sched.idx;
required = {'L','beam','parConst','idx'};
for k = 1:numel(required)
    if ~isfield(sched,required{k})
        error('AeroFlex:sched:ReactionContractField', ...
            'Scheduled package is missing field "%s".',required{k});
    end
end
if numel(x) ~= size(sched.L,2) || numel(u) ~= 4
    error('AeroFlex:sched:ReactionContractDimensions', ...
        'Expected %d states and four wing-control inputs.',size(sched.L,2));
end

pc = sched.parConst;
pc.u_ctrl = u;
pc.gust = [];
if ~isfield(pc,'N_Thrust') || isempty(pc.N_Thrust)
    pc.N_Thrust = zeros(numel(idx.q1),1);
end

linearQ1 = sched.L(idx.q1,:)*x;

pcStructural = pc;
pcStructural.forces_0 = zeros(size(pc.forces_0));
pcStructural.N_Thrust = zeros(size(pc.N_Thrust));
pcStructural.u_ctrl = [];
pcStructural.gust = [];
if isfield(pcStructural,'affineOffset')
    pcStructural.affineOffset = [];
end
structuralN = AeroFlex.sim.nonlinear_terms(x,pcStructural,idx);
structuralQ1 = structuralN(idx.q1);

preloadQ1 = pc.Fscale*sched.beam.eta_e(:);
steadyAeroQ1 = pc.Fscale*(pc.forces_0(:)-sched.beam.eta_e(:));
thrustQ1 = pc.N_Thrust(:);

pcNoControl = pc;
pcNoControl.u_ctrl = [];
noControlN = AeroFlex.sim.nonlinear_terms(x,pcNoControl,idx);
fullN = AeroFlex.sim.nonlinear_terms(x,pc,idx);
controlQ1 = fullN(idx.q1)-noControlN(idx.q1);

affineQ1 = zeros(numel(idx.q1),1);
if isfield(pc,'affineOffset') && ~isempty(pc.affineOffset)
    affineQ1 = pc.affineOffset(idx.q1);
end

ownedQ1 = linearQ1+structuralQ1+preloadQ1+steadyAeroQ1+ ...
    thrustQ1+controlQ1+affineQ1;
raw = sched.L*x+fullN;
rawQ1 = raw(idx.q1);

components = struct( ...
    'linearOperator',linearQ1, ...
    'structuralNonlinear',structuralQ1, ...
    'preload',preloadQ1, ...
    'steadyAerodynamic',steadyAeroQ1, ...
    'thrust',thrustQ1, ...
    'control',controlQ1, ...
    'affineOffset',affineQ1);
names = fieldnames(components);
wrenchComponents = struct();
for k = 1:numel(names)
    wrenchComponents.(names{k}) = recover(components.(names{k}));
end

c2 = recover(rawQ1);
c1 = sched.beam.red.phi2_sA*x(idx.q2);

c3 = nan(6,1);
if ~isempty(fieldnames(loadLedger))
    requiredLoads = {'Ftot_B','Mtot_B','Ftail_B','Mtail_B','Ffin_B', ...
        'Mfin_B','Fgrav_B','Fthrust_B','Mthrust_B'};
    for k = 1:numel(requiredLoads)
        if ~isfield(loadLedger,requiredLoads{k})
            error('AeroFlex:sched:ReactionLedgerField', ...
                'Load ledger is missing field "%s".',requiredLoads{k});
        end
    end
    c3(1:3) = loadLedger.Ftot_B(:)-loadLedger.Ftail_B(:)- ...
        loadLedger.Ffin_B(:)-loadLedger.Fgrav_B(:)- ...
        loadLedger.Fthrust_B(:);
    c3(4:6) = loadLedger.Mtot_B(:)-loadLedger.Mtail_B(:)- ...
        loadLedger.Mfin_B(:)-loadLedger.Mthrust_B(:);
end

sumWrench = zeros(6,1);
for k = 1:numel(names)
    sumWrench = sumWrench+wrenchComponents.(names{k});
end

contract = struct();
contract.schemaVersion = 3;
contract.frame = 'body';
contract.referencePoint = 'wing root';
contract.units = {'N','N','N','N m','N m','N m'};
contract.multiplicity = 'full wing; no additional mirroring';
contract.componentOrder = {'Fx','Fy','Fz','Mx','My','Mz'};
contract.C1 = struct('method','phi2_sA*q2 intrinsic candidate', ...
    'wrench',c1,'differenceFromC2',c1-c2);
contract.C2 = struct('method','raw constraint-bearing q1 equation', ...
    'wrench',c2,'generalizedReaction',sched.beam.Pr*rawQ1, ...
    'components',wrenchComponents,'componentClosure',sumWrench-c2, ...
    'rawOwnershipClosure',ownedQ1-rawQ1);
contract.C3 = struct('method','rigid-body global load ledger', ...
    'wrench',c3,'differenceFromC2',c3-c2);
contract.ownership = struct( ...
    'wingOnly',true,'tailInClamp',false,'gravityInClamp',false, ...
    'rigidThrustInClamp',false,'directAeroAddedSeparately',false, ...
    'affineOffsetAddedOnce',true,'outputFeedsDynamics',false);

    function wrench = recover(q1Contribution)
        wrench = AeroFlex.beam.recoverRootWrench( ...
            sched.beam.red.phi1_sA,sched.beam.Pr, ...
            sched.beam.Pr*q1Contribution(:));
    end
end
