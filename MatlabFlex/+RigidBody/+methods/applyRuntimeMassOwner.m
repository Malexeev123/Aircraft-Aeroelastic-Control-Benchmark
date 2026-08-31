function [cfg,info] = applyRuntimeMassOwner(cfg,owner)
%APPLYRUNTIMEMASSOWNER Install one rigid-body mass/reference contract.
%   The non-wing owner is used when the retained wing source supplies the
%   wing inertial reaction.  The legacy owner remains available for exact
%   historical reproduction.

arguments
    cfg struct
    owner (1,1) string
end

parameters = RigidBody.methods.paramsRigid_PazyUAV();
owner = lower(strtrim(owner));
switch owner
    case "nonwing_source_owned_wing_reaction"
        selected = parameters.nonwing;
        sourceWingOwned = true;
    case "legacy_combined_project_tips"
        selected = parameters.legacyIncludingProjectTipMass;
        sourceWingOwned = false;
    otherwise
        error("RigidBody:RuntimeMassOwner", ...
            "Unsupported runtime mass owner %s.",owner);
end

assert(isfinite(selected.mass) && selected.mass>0 && ...
    isequal(size(selected.J),[3,3]) && ...
    all(isfinite(selected.J),"all"), ...
    "RigidBody:RuntimeMassDefinition", ...
    "The selected runtime mass definition is incomplete or nonfinite.");
selected.J = 0.5*(selected.J+selected.J.');
assert(min(eig(selected.J))>0, ...
    "RigidBody:RuntimeMassInertia", ...
    "The selected runtime inertia must be positive definite.");

if ~isfield(cfg,"rigidEOMset") || ~isstruct(cfg.rigidEOMset)
    cfg.rigidEOMset = struct();
end
cfg.rigidEOMset.mass = selected.mass;
cfg.rigidEOMset.I_B = selected.J;
if ~isfield(cfg.rigidEOMset,"rThrust_B") || ...
        isempty(cfg.rigidEOMset.rThrust_B)
    cfg.rigidEOMset.rThrust_B = zeros(3,1);
end

if sourceWingOwned
    cfg.rigidEOMset.rWingRoot_B = selected.wingRootFromCM;
    cfg.rigidEOMset.rTail_B = selected.tailArmFromCM;
    if ~isfield(cfg,"tail") || ~isstruct(cfg.tail)
        cfg.tail = struct();
    end
    cfg.tail.r_B = selected.tailArmFromCM;
end

changeId = ...
    "phase18c-v17a-nonlinear-runtime-nonwing-mass-owner-promotion-v1";
cfg.rigidEOMset.massOwnership = struct( ...
    "owner",owner, ...
    "sourceWingOwned",sourceWingOwned, ...
    "changeControlId",changeId, ...
    "massKilograms",selected.mass, ...
    "centerOfMassFromWingRootMeters",selected.rCM(:), ...
    "inertiaAboutSelectedCenterKgM2",selected.J);

info = cfg.rigidEOMset.massOwnership;
info.projectTipMassExcludedKilograms = 0;
if sourceWingOwned
    info.projectTipMassExcludedKilograms = ...
        parameters.legacyIncludingProjectTipMass.projectTipMass;
end
end
