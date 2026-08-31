function sched = installP5R2T2Anchor(sched,trim,category)
%INSTALLP5R2T2ANCHOR Install a model-owned T2 anchor after final trim.

arguments
    sched (1,1) struct
    trim (1,1) struct
    category (1,1) string {mustBeMember(category, ...
        ["SOURCE_TRAINING_T2_ANCHOR","MODEL_PREDICTED_QUERY_T2_ANCHOR"])}
end

if ~trim.converged || numel(trim.states) ~= size(sched.L,1)
    error('AeroFlex:sched:R2AnchorFinalTrim', ...
        'A final converged trim with the active state dimension is required.');
end
u = [trim.deltaWing(1);trim.deltaWing(1);0;0];
stateHash = vectorHash(trim.states(:));
commandHash = vectorHash(u);
structuralHash = contractHash(struct('Pz',sched.beam.Pz, ...
    'Pr',sched.beam.Pr,'eta_e',sched.beam.eta_e, ...
    'phi1_sA',sched.beam.red.phi1_sA));
rawLHash = contractHash(sched.L);
recoveryHash = contractHash(struct('Pr',sched.beam.Pr, ...
    'phi1_sA',sched.beam.red.phi1_sA));

sched.x_eq = trim.states(:);
sched.u_eq = u;
sched = recenterAtFinalTrim(sched);
sched.p5.t2AnchorSource = struct('query',sched.mu, ...
    'trimStateHash',stateHash,'trimCommandHash',commandHash, ...
    'structuralPackageHash',structuralHash,'rawLHash',rawLHash, ...
    'recoveryMapHash',recoveryHash);
anchor = sched.p5.r2.anchor;
anchor.provenanceCategory = char(category);
anchor.construction = 'raw_L_T2_reaction_at_final_matching_equilibrium';
anchor.interpolated = false;
anchor.wrench = trim.debug.loads.Clamp6(:);
anchor.query = sched.mu;
anchor.trimStateHash = stateHash;
anchor.trimCommandHash = commandHash;
anchor.structuralPackageHash = structuralHash;
anchor.rawLHash = rawLHash;
anchor.recoveryMapHash = recoveryHash;
sched.p5.r2.anchor = anchor;
if isfield(sched,'equilibriumCentered') && ...
        isfield(sched.equilibriumCentered,'recoveryVertices')
    sched.equilibriumCentered.recoveryAnchorWrench = anchor.wrench(:);
end
sched.p5.packageHash = packageHash(sched);
AeroFlex.sched.assertR2AnchorProvenance(sched);
end

function sched = recenterAtFinalTrim(sched)
% The fixed-basis package is initially centered before the final coupled
% trim. Installing a different final state/command must translate the affine
% term to that declared equilibrium without changing the state or operators.
if ~isfield(sched,'equilibriumCentered') || ...
        ~isfield(sched.equilibriumCentered,'enabled') || ...
        ~sched.equilibriumCentered.enabled || ...
        ~isfield(sched.parConst,'affineOffset') || ...
        isempty(sched.parConst.affineOffset)
    return
end

pc = sched.parConst;
pc.u_ctrl = sched.u_eq(:);
pc.gust = [];
affineBefore = pc.affineOffset(:);
residualBefore = sched.L*sched.x_eq(:)+ ...
    AeroFlex.sim.nonlinear_terms(sched.x_eq(:),pc,sched.idx);
pc.affineOffset = affineBefore-residualBefore;
residualAfter = sched.L*sched.x_eq(:)+ ...
    AeroFlex.sim.nonlinear_terms(sched.x_eq(:),pc,sched.idx);
if norm(residualAfter) > 1e-10
    error('AeroFlex:sched:FinalTrimRecenteringResidual', ...
        'Final trim affine re-centering residual is %.3e.', ...
        norm(residualAfter));
end

sched.parConst = pc;
sched.equilibriumCentered.rEq = residualAfter;
sched.equilibriumCentered.finalTrimRecentering = struct( ...
    'applied',true,'equationChanged',false,'stateChanged',false, ...
    'operatorChanged',false,'residualBefore2',norm(residualBefore), ...
    'residualBeforeQGamma2',norm(residualBefore(sched.idx.qGam)), ...
    'residualAfter2',norm(residualAfter), ...
    'affineBeforeHash',contractHash(affineBefore), ...
    'affineAfterHash',contractHash(pc.affineOffset));
if isfield(sched,'p5')
    sched.p5.aerodynamicResidual = norm(residualAfter(sched.idx.qGam));
    sched.p5.centeredResidual = norm(residualAfter);
    sched.p5.finalTrimRecentering = ...
        sched.equilibriumCentered.finalTrimRecentering;
end
end

function digest = vectorHash(value)
engine=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
engine.update(typecast(double(value(:)),'uint8'));
digest=hexDigest(engine.digest());
end

function digest = contractHash(value)
engine=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
engine.update(uint8(unicode2native(jsonencode(value),'UTF-8')));
digest=hexDigest(engine.digest());
end

function digest = packageHash(sched)
payload=struct('mu',sched.mu,'L',full(sched.L),'xEq',sched.x_eq, ...
    'uEq',sched.u_eq,'forces0',sched.parConst.forces_0, ...
    'affineOffset',sched.parConst.affineOffset, ...
    'sourceHash',sched.p5.sourceContractHash, ...
    'referenceHash',sched.p5.referenceContractHash);
digest=contractHash(payload);
end

function digest = hexDigest(bytes)
bytes=typecast(bytes,'uint8');
digest=lower(string(reshape(dec2hex(bytes,2).',1,[])));
end
