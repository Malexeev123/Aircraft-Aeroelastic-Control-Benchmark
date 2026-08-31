function [stateNext,sensitivity,rootWrench,stateHistory] = ...
        scheduledReciprocalIntervalKernelAudit(state,disturbance,wingTotal, ...
        wingIncrement,elevatorIncrement,thrustIncrement,rigidState, ...
        needJacobian,packet)
%SCHEDULEDRECIPROCALINTERVALKERNELAUDIT Scheduled aggregate interval map.
%   This default-inactive numerical kernel evaluates the existing 620-state
%   common-coordinate aggregate over fourteen 0.01/14 s source steps.

state = state(:);
elevatorIncrement = elevatorIncrement(:).';
thrustIncrement = thrustIncrement(:).';
assert(numel(state)==620 && isscalar(disturbance) && ...
    isequal(size(wingTotal),[4 14]) && ...
    isequal(size(wingIncrement),[4 14]) && ...
    isequal(size(elevatorIncrement),[1 14]) && ...
    isequal(size(thrustIncrement),[1 14]) && ...
    isequal(size(rigidState),[9 14]), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalDimensions', ...
    'The scheduled reciprocal interval inputs have invalid dimensions.');
assert(all(isfinite(state)) && isfinite(disturbance) && ...
    all(isfinite(wingTotal),'all') && ...
    all(isfinite(wingIncrement),'all') && ...
    all(isfinite(elevatorIncrement)) && ...
    all(isfinite(thrustIncrement)) && ...
    all(isfinite(rigidState),'all'), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalFinite', ...
    'The scheduled reciprocal interval inputs must be finite.');
assert(all(packet.queryTrimThrust+thrustIncrement>=0), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalThrust', ...
    'The scheduled reciprocal interval violates T >= 0.');
assert(norm(wingTotal-(packet.queryTrimWing+wingIncrement),inf)<= ...
    100*eps(max(1,norm(wingTotal,inf))), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalWingOwner', ...
    'Total and trim-relative scheduled wing commands disagree.');
assert(isscalar(needJacobian), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalJacobianFlag', ...
    'The scheduled interval Jacobian flag must be scalar.');
needJacobian = logical(needJacobian);

sourceMask = logical(packet.sensitivitySourceMask(:));
sourceCounts = double(packet.hiddenStateCounts(:));
assert(numel(sourceMask)==14 && numel(sourceCounts)==14 && ...
    all(sourceCounts==39) && sum(sourceCounts)==546 && ...
    numel(packet.interpolationWeights)==14 && ...
    all(sourceMask(packet.interpolationWeights(:)>0)), ...
    'AeroFlex:ctrl:ScheduledReciprocalIntervalStencil', ...
    'The scheduled sensitivity stencil omits an active source member.');
useReducedSensitivity = needJacobian && ~all(sourceMask);
if useReducedSensitivity
    activeHidden = localActiveHiddenIndices(sourceMask,sourceCounts);
    activeState = [reshape(1:74,[],1);74+activeHidden];
    inactiveHidden = localInactiveHiddenIndices(activeHidden);
    assert(norm(packet.Avh(:,inactiveHidden),inf)<=1e-14 && ...
        norm(packet.Ahh(activeHidden,inactiveHidden),inf)<=1e-14 && ...
        norm(packet.Ahh(inactiveHidden,activeHidden),inf)<=1e-14 && ...
        norm(packet.rootWrenchStateMap(:,83+inactiveHidden),inf)<=1e-14, ...
        'AeroFlex:ctrl:ScheduledReciprocalIntervalStencilStructure', ...
        'An omitted source member can influence an active scheduled output.');
    stateSensitivity = eye(numel(activeState));
    disturbanceSensitivity = zeros(numel(activeState),1);
    wingSensitivity = zeros(numel(activeState),4);
else
    activeHidden = reshape(1:546,[],1);
    activeState = reshape(1:620,[],1);
    stateSensitivity = zeros(620,620);
    disturbanceSensitivity = zeros(620,1);
    wingSensitivity = zeros(620,4);
    if needJacobian
        stateSensitivity = eye(620);
    end
end
rootWrench = zeros(6,14);
stateHistory = zeros(620,14);

for substep = 1:14
    nativeState = state(1:74);
    hiddenState = state(75:620);
    sourceInput = [disturbance;wingIncrement(:,substep); ...
        elevatorIncrement(substep);thrustIncrement(substep)];
    runtimeDeparture = packet.runtimeFromNative*nativeState+ ...
        packet.runtimeFromRigid*rigidState(:,substep)- ...
        packet.correctionStateEquilibrium;
    inputDeparture = sourceInput-packet.inputEquilibrium;
    visibleCorrection = packet.deltaAFullState*runtimeDeparture+ ...
        packet.Avh*hiddenState+packet.deltaBv*inputDeparture;
    hiddenNext = packet.hiddenFromFullState*runtimeDeparture+ ...
        packet.Ahh*hiddenState+packet.Bh*inputDeparture;

    if needJacobian
        [nativeNext,nativeSensitivity] = localNativeStep( ...
            nativeState,wingTotal(:,substep),disturbance,packet);
    else
        nativeNext = localNativeValueStep( ...
            nativeState,wingTotal(:,substep),disturbance,packet);
        nativeSensitivity = zeros(74,79);
    end
    nativeNext = nativeNext+packet.correctionToNative*visibleCorrection;
    stateNext = [nativeNext;hiddenNext];
    assert(all(isfinite(stateNext)), ...
        'AeroFlex:ctrl:ScheduledReciprocalIntervalOutputFinite', ...
        'The scheduled reciprocal interval produced a nonfinite state.');
    rootWrench(:,substep) = packet.rootWrenchEquilibrium+ ...
        packet.rootWrenchStateMap*[runtimeDeparture;hiddenState]+ ...
        packet.rootWrenchInputMap*inputDeparture;
    assert(all(isfinite(rootWrench(:,substep))), ...
        'AeroFlex:ctrl:ScheduledReciprocalIntervalWrenchFinite', ...
        'The scheduled reciprocal interval produced a nonfinite root wrench.');

    if needJacobian
        nativeA = nativeSensitivity(:,1:74);
        nativeInput = nativeSensitivity(:,75:79);
        transitionA = [ ...
            nativeA+packet.correctionToNative*packet.deltaAFullState* ...
                packet.runtimeFromNative, ...
            packet.correctionToNative*packet.Avh; ...
            packet.hiddenFromFullState*packet.runtimeFromNative, ...
            packet.Ahh];
        sourceInputNative = zeros(74,7);
        sourceInputNative(:,1:5) = nativeInput;
        transitionInput = [ ...
            sourceInputNative+packet.correctionToNative*packet.deltaBv; ...
            packet.Bh];
        if useReducedSensitivity
            transitionReduced = transitionA(activeState,activeState);
            inputReduced = transitionInput(activeState,:);
            disturbanceSensitivity = transitionReduced* ...
                disturbanceSensitivity+inputReduced(:,1);
            wingSensitivity = transitionReduced*wingSensitivity+ ...
                inputReduced(:,2:5);
            stateSensitivity = transitionReduced*stateSensitivity;
        else
            disturbanceSensitivity = transitionA*disturbanceSensitivity+ ...
                transitionInput(:,1);
            wingSensitivity = transitionA*wingSensitivity+ ...
                transitionInput(:,2:5);
            stateSensitivity = transitionA*stateSensitivity;
        end
    end
    state = stateNext;
    stateHistory(:,substep) = stateNext;
end
if useReducedSensitivity
    sensitivity = zeros(620,625);
    sensitivity(activeState,activeState) = stateSensitivity;
    sensitivity(activeState,621) = disturbanceSensitivity;
    sensitivity(activeState,622:625) = wingSensitivity;
else
    sensitivity = [stateSensitivity,disturbanceSensitivity,wingSensitivity];
end
end

function indices = localInactiveHiddenIndices(activeIndices)
% Preserve ascending source-memory order without a generated set operation.
active = false(546,1);
active(activeIndices) = true;
indices = zeros(546-numel(activeIndices),1);
cursor = 0;
for index = 1:546
    if ~active(index)
        cursor = cursor+1;
        indices(cursor) = index;
    end
end
end

function indices = localActiveHiddenIndices(sourceMask,sourceCounts)
offsets = [0;cumsum(sourceCounts(:))];
indices = zeros(sum(sourceCounts(sourceMask)),1);
cursor = 0;
for source = 1:numel(sourceMask)
    if sourceMask(source)
        count = sourceCounts(source);
        indices(cursor+(1:count)) = offsets(source)+(1:count);
        cursor = cursor+count;
    end
end
end

function stateNext = localNativeValueStep(state,wingTotal,disturbance,packet)
par = packet.par;
par.gust = disturbance;
par.u_ctrl = wingTotal;
dt = packet.dt;
gamma = packet.gamma;
delta = packet.delta;
k1Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(state,par));
k1Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+gamma*k1Nonlinear));
stage1 = state+gamma*(k1Linear+k1Nonlinear);
k2Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(stage1,par));
k2Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+(1-gamma)*k1Linear+delta*k1Nonlinear+ ...
    (1-delta)*k2Nonlinear));
stage2 = state+(1-gamma)*k1Linear+gamma*k2Linear+ ...
    delta*k1Nonlinear+(1-delta)*k2Nonlinear;
k3Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(stage2,par));
stateNext = state+(1-gamma)*(k1Linear+k2Nonlinear)+ ...
    gamma*(k2Linear+k3Nonlinear);
end

function [stateNext,sensitivityNext] = ...
        localNativeStep(state,wingTotal,disturbance,packet)
par = packet.par;
par.gust = disturbance;
par.u_ctrl = wingTotal;
dt = packet.dt;
gamma = packet.gamma;
delta = packet.delta;
sensitivity = [eye(74),zeros(74,5)];

k1Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(state,par));
[jacobian1,inputJacobian1] = localNonlinearJacobian(state,par);
j1Nonlinear = packet.rateProjection* ...
    (dt*(jacobian1*sensitivity+inputJacobian1));
j1Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (sensitivity+gamma*j1Nonlinear));
sensitivity1 = sensitivity+gamma*(j1Linear+j1Nonlinear);
k1Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+gamma*k1Nonlinear));
stage1 = state+gamma*(k1Linear+k1Nonlinear);

k2Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(stage1,par));
[jacobian2,inputJacobian2] = localNonlinearJacobian(stage1,par);
j2Nonlinear = packet.rateProjection* ...
    (dt*(jacobian2*sensitivity1+inputJacobian2));
j2Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (sensitivity+(1-gamma)*j1Linear+delta*j1Nonlinear+ ...
    (1-delta)*j2Nonlinear));
sensitivity2 = sensitivity+(1-gamma)*j1Linear+gamma*j2Linear+ ...
    delta*j1Nonlinear+(1-delta)*j2Nonlinear;
k2Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+(1-gamma)*k1Linear+delta*k1Nonlinear+ ...
    (1-delta)*k2Nonlinear));
stage2 = state+(1-gamma)*k1Linear+gamma*k2Linear+ ...
    delta*k1Nonlinear+(1-delta)*k2Nonlinear;

k3Nonlinear = packet.rateProjection*(dt*localNonlinearTerms(stage2,par));
[jacobian3,inputJacobian3] = localNonlinearJacobian(stage2,par);
j3Nonlinear = packet.rateProjection* ...
    (dt*(jacobian3*sensitivity2+inputJacobian3));
sensitivityNext = sensitivity+(1-gamma)*(j1Linear+j2Nonlinear)+ ...
    gamma*(j2Linear+j3Nonlinear);
stateNext = state+(1-gamma)*(k1Linear+k2Nonlinear)+ ...
    gamma*(k2Linear+k3Nonlinear);
end

function value = localImplicitSolve(packet,rightHandSide)
permuted = rightHandSide(packet.piv,:);
value = packet.Ufac\(packet.Lfac\permuted);
end

function nonlinear = localNonlinearTerms(state,par)
q1 = state(1:10);
q2 = state(11:20);
qxi = state(21:31);
nonlinear = zeros(74,1);
gamma1q1 = squeeze(pagemtimes(par.Gamma1,'none',q1,'none'));
gamma2q2 = squeeze(pagemtimes(par.Gamma2,'none',q2,'none'));
gammaGqxi = squeeze(pagemtimes(par.Gamma_g,'none',qxi,'none'));
gammaXiqxi = squeeze(pagemtimes(par.Gamma_xi,'none',qxi,'none'));
gamma2TransposeQ1 = squeeze(pagemtimes( ...
    par.Gamma2,'transpose',q1,'none'));
nonlinear(1:10) = -(gamma1q1*q1)-(gamma2q2*q2)+ ...
    gammaGqxi*qxi+par.Fscale*par.forces_0+par.N_Thrust;
nonlinear(11:20) = gamma2TransposeQ1*q2;
nonlinear(21:31) = gammaXiqxi*q1;
nonlinear(32:71) = par.Bw*par.gust+ ...
    (1/par.t_inf)*par.Bdel*par.u_ctrl(1:2)+ ...
    par.Bddel*par.u_ctrl(3:4);
nonlinear(1:10) = nonlinear(1:10)+ ...
    par.scaleAero*par.Dw*par.gust+ ...
    par.Ddel*par.u_ctrl(1:2)+ ...
    par.t_inf*par.Dddel*par.u_ctrl(3:4);
nonlinear = nonlinear+par.affineOffset;
end

function [stateJacobian,inputJacobian] = localNonlinearJacobian(state,par)
q1 = state(1:10);
q2 = state(11:20);
qxi = state(21:31);
gamma1 = par.Gamma1;
gamma2 = par.Gamma2;
gammaG = par.Gamma_g;
gammaXi = par.Gamma_xi;
termB = localGammaContract(permute(gamma1,[1 3 2]),q1);
termA = localGammaContract(gamma1,q1);
termD = localGammaContract(permute(gamma2,[1 3 2]),q2);
termC = localGammaContract(gamma2,q2);
termF = localGammaContract(gammaG,qxi);
termE = localGammaContract(permute(gammaG,[1 3 2]),qxi);
termJ = localGammaContract(gammaXi,qxi);
termI = localGammaContract(permute(gammaXi,[1 3 2]),q1);
termG = zeros(10,10);
termH = zeros(10,10);
for mode = 1:10
    termG = termG+q2(mode)*gamma2(:,:,mode).';
    termH(:,mode) = gamma2(:,:,mode).'*q1;
end
stateJacobian = zeros(74,74);
stateJacobian(1:10,1:10) = -termA-termB;
stateJacobian(1:10,11:20) = -termC-termD;
stateJacobian(1:10,21:31) = termE+termF;
stateJacobian(11:20,1:10) = termG;
stateJacobian(11:20,11:20) = termH;
stateJacobian(21:31,1:10) = termJ;
stateJacobian(21:31,21:31) = termI;
inputJacobian = zeros(74,79);
inputJacobian(1:10,75) = par.scaleAero*par.Dw;
inputJacobian(32:71,75) = par.Bw;
inputJacobian(1:10,76:77) = par.Ddel;
inputJacobian(1:10,78:79) = par.t_inf*par.Dddel;
inputJacobian(32:71,76:77) = (1/par.t_inf)*par.Bdel;
inputJacobian(32:71,78:79) = par.Bddel;
end

function contracted = localGammaContract(tensor,coordinates)
pageCoordinates = repmat(reshape(coordinates,[],1,1), ...
    1,1,size(tensor,3));
contractedPages = pagemtimes(tensor,pageCoordinates);
contracted = reshape(contractedPages,size(tensor,1),size(tensor,3));
end
