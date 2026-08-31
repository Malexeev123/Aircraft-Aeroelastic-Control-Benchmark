function result = realizeVirtualDemandAudit(model, zDeviation, uNominal, ...
        wingRealized, bounds, wingRate)
%REALIZEVIRTUALDEMANDAUDIT Audit-only tail/thrust realization of a virtual demand.
%   The nominal three-channel command defines the requested next rigid state.
%   The shared wing is fixed to the already authorized realized value; this
%   routine chooses only elevator and thrust residuals within existing bounds.
%   The two-variable convex problem is solved by exhaustive bound active sets,
%   not by a general-purpose optimizer.

arguments
    model struct
    zDeviation double
    uNominal double
    wingRealized (1,1) double
    bounds struct
    wingRate (1,1) struct = struct()
end

requiredModel = {"A","B","stateScale","inputScale"};
for name = requiredModel
    assert(isfield(model,name{1}),"AeroFlex:VirtualDemand:Model", ...
        "Missing model.%s.",name{1});
end
assert(isequal(size(model.A),[6,6]) && isequal(size(model.B),[6,3]), ...
    "AeroFlex:VirtualDemand:Dimensions", ...
    "The realization model must be six-state and three-input.");
zDeviation = zDeviation(:);
uNominal = uNominal(:);
assert(numel(zDeviation)==6 && numel(uNominal)==3 && ...
    all(isfinite([zDeviation;uNominal;wingRealized]),"all"), ...
    "AeroFlex:VirtualDemand:Finite", "State and command inputs must be finite.");
assert(isfield(bounds,"elevatorLower") && isfield(bounds,"elevatorUpper") && ...
    isfield(bounds,"thrustLower") && isfield(bounds,"thrustUpper"), ...
    "AeroFlex:VirtualDemand:Bounds", "All tail/thrust bounds are required.");
lower = [double(bounds.elevatorLower);double(bounds.thrustLower)];
upper = [double(bounds.elevatorUpper);double(bounds.thrustUpper)];
lower = lower(:);
upper = upper(:);
assert(numel(lower)==2 && numel(upper)==2 && all(isfinite([lower;upper])) && ...
    all(lower<=upper), ...
    "AeroFlex:VirtualDemand:Bounds", "The realization bounds are invalid.");

stateScale = double(model.stateScale(:));
inputScale = double(model.inputScale(:));
assert(numel(stateScale)==4 && numel(inputScale)==3 && ...
    all(stateScale>0) && all(inputScale>0), ...
    "AeroFlex:VirtualDemand:Scale", "Positive four-state and three-input scales are required.");

% Only physical rigid outputs [u,w,theta,q] are matched.  Fusion and servo
% memory remain in the source-derived prediction but are not fictitious
% tracking targets.
rigidRows = 1:4;
W = diag(1./stateScale);
G = W*double(model.B(rigidRows,2:3));
desired = double(model.A(rigidRows,:))*zDeviation + ...
    double(model.B(rigidRows,:))*uNominal;
fixed = double(model.A(rigidRows,:))*zDeviation + ...
    double(model.B(rigidRows,1))*wingRealized;
rateCorrection = 0;
if ~isempty(fieldnames(wingRate))
    assert(isfield(model,"wingRateB") && ...
        isequal(size(model.wingRateB),[size(model.B,1),1]) && ...
        all(isfinite(model.wingRateB),'all') && ...
        isfield(wingRate,"nominal") && isfield(wingRate,"realized") && ...
        isscalar(wingRate.nominal) && isscalar(wingRate.realized) && ...
        all(isfinite([wingRate.nominal;wingRate.realized])), ...
        "AeroFlex:VirtualDemand:WingRate", ...
        "The physical wing-rate realization contract is invalid.");
    desired = desired+double(model.wingRateB(rigidRows))* ...
        double(wingRate.nominal);
    fixed = fixed+double(model.wingRateB(rigidRows))* ...
        double(wingRate.realized);
    rateCorrection = double(wingRate.realized-wingRate.nominal);
end
target = W*(desired-fixed);

% The secondary nominal-command term is normalized by declared source input
% scales.  It resolves rank deficiency without changing the primary rigid
% demand or introducing an empirical tail boost.
R = diag(1./inputScale(2:3));
H = G.'*G + R.'*R;
g = G.'*target + R.'*R*uNominal(2:3);
[tailThrust,objective,activeSet] = localBoundedLeastSquares(H,g,lower,upper);
assert(all(isfinite(tailThrust)),"AeroFlex:VirtualDemand:Algebraic", ...
    "The bounded algebraic realization returned a nonfinite point.");
realized = [wingRealized;tailThrust];
predictedRigid = double(model.A(rigidRows,:))*zDeviation + ...
    double(model.B(rigidRows,:))*realized;
residual = predictedRigid-desired;
result = struct("accepted",true,"solver","explicit_two_variable_active_set", ...
    "objective",objective,"activeSet",activeSet,"nominalCommand",uNominal, ...
    "realizedCommand",realized,"wingCorrection",wingRealized-uNominal(1), ...
    "wingRateCorrection",rateCorrection, ...
    "desiredRigidNext",desired,"predictedRigidNext",predictedRigid, ...
    "rigidResidual",residual,"normalizedRigidResidual",W*residual, ...
    "normalizedRigidResidualInf",norm(W*residual,inf), ...
    "totalThrustLowerSatisfied",isfield(bounds,"trimThrust") && ...
        double(bounds.trimThrust)+tailThrust(2)>=-10*eps(max(1,double(bounds.trimThrust))));
end

function [solution,bestObjective,activeSet] = localBoundedLeastSquares(H,g,lower,upper)
% Enumerate the nine possible free/lower/upper statuses for two variables.
status = [0 0;0 -1;0 1;-1 0;-1 -1;-1 1;1 0;1 -1;1 1];
bestObjective = inf;
solution = nan(2,1);
activeSet = strings(2,1);
for row = 1:size(status,1)
    candidate = zeros(2,1);
    fixed = status(row,:)~=0;
    if any(fixed)
        candidate(fixed) = lower(fixed).*(status(row,fixed)'<0) + ...
            upper(fixed).*(status(row,fixed)'>0);
    end
    free = ~fixed;
    if any(free)
        candidate(free) = H(free,free)\(g(free)-H(free,fixed)*candidate(fixed));
    end
    if any(candidate<lower-1e-12) || any(candidate>upper+1e-12) || ...
            any(~isfinite(candidate))
        continue
    end
    objective = candidate.'*H*candidate-2*g.'*candidate;
    if objective<bestObjective
        bestObjective = objective;
        solution = min(max(candidate,lower),upper);
        activeSet = strings(2,1);
        activeSet(status(row,:)'<0) = "lower";
        activeSet(status(row,:)'==0) = "free";
        activeSet(status(row,:)'>0) = "upper";
    end
end
assert(isfinite(bestObjective),"AeroFlex:VirtualDemand:ActiveSet", ...
    "No bounded two-variable active-set realization was feasible.");
end
