function [stateNext,sensitivity,rootWrench,stateHistory] = ...
        fixedReciprocalIntervalKernelAudit(state,disturbance,wingTotal, ...
        wingIncrement,elevatorIncrement,thrustIncrement,rigidState, ...
        needJacobian,packet)
%FIXEDRECIPROCALINTERVALKERNELAUDIT Fixed-size exact Case-A interval map.
%   This default-inactive numerical kernel reproduces the approved
%   ReciprocalControllerModelProvider interval transition. It retains all
%   sixteen 0.000625 s source steps and returns the complete endpoint,
%   state/disturbance sensitivity, and source-root wrench history.

state = state(:);
elevatorIncrement = elevatorIncrement(:).';
thrustIncrement = thrustIncrement(:).';
substepCount = size(wingTotal,2);
assert(numel(state)==121 && isscalar(disturbance) && ...
    isequal(size(wingTotal),[4 16]) && ...
    isequal(size(wingIncrement),[4 16]) && ...
    isequal(size(elevatorIncrement),[1 16]) && ...
    isequal(size(thrustIncrement),[1 16]) && ...
    isequal(size(rigidState),[9 16]) && substepCount==16, ...
    'AeroFlex:ctrl:FixedReciprocalIntervalDimensions', ...
    'The fixed reciprocal interval inputs have invalid dimensions.');
assert(all(isfinite(state)) && isfinite(disturbance) && ...
    all(isfinite(wingTotal),'all') && ...
    all(isfinite(wingIncrement),'all') && ...
    all(isfinite(elevatorIncrement)) && ...
    all(isfinite(thrustIncrement)) && ...
    all(isfinite(rigidState),'all'), ...
    'AeroFlex:ctrl:FixedReciprocalIntervalFinite', ...
    'The fixed reciprocal interval inputs must be finite.');

assert(isscalar(needJacobian), ...
    'AeroFlex:ctrl:FixedReciprocalIntervalJacobianFlag', ...
    'The fixed reciprocal interval Jacobian flag must be scalar.');
needJacobian = logical(needJacobian);
stateSensitivity = zeros(121,121);
disturbanceSensitivity = zeros(121,1);
wingSensitivity = zeros(121,4);
if needJacobian
    stateSensitivity = eye(121);
end
rootWrench = zeros(6,16);
stateHistory = zeros(121,16);

for substep = 1:16
    nativeState = state(1:74);
    hiddenState = state(75:121);
    sourceInput = [disturbance;wingIncrement(:,substep); ...
        elevatorIncrement(substep);thrustIncrement(substep)];
    runtimeState = packet.runtimeFromNative*nativeState+ ...
        packet.runtimeFromRigid*rigidState(:,substep);
    runtimeDeparture = runtimeState-packet.correctionStateEquilibrium;
    inputDeparture = sourceInput-packet.inputEquilibrium;

    visibleCorrection = ...
        packet.deltaAFullState*runtimeDeparture+ ...
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

    rootWrench(:,substep) = packet.rootWrenchEquilibrium+ ...
        packet.rootWrenchStateMap*[runtimeDeparture;hiddenState]+ ...
        packet.rootWrenchInputMap*inputDeparture;
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
        disturbanceSensitivity = transitionA*disturbanceSensitivity+ ...
            transitionInput(:,1);
        wingSensitivity = transitionA*wingSensitivity+ ...
            transitionInput(:,2:5);
        stateSensitivity = transitionA*stateSensitivity;
    end
    state = stateNext;
    stateHistory(:,substep) = stateNext;
end

sensitivity = [stateSensitivity,disturbanceSensitivity,wingSensitivity];
end

function stateNext = localNativeValueStep( ...
        state,wingTotal,disturbance,packet)
%LOCALNATIVEVALUESTEP Three-stage IMEX transition without sensitivities.
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
%LOCALNATIVESTEP Three-stage IMEX transition and analytical sensitivity.
par = packet.par;
par.gust = disturbance;
par.u_ctrl = wingTotal;
dt = packet.dt;
gamma = packet.gamma;
delta = packet.delta;
sensitivity = [eye(74),zeros(74,5)];

k1Nonlinear = dt*localNonlinearTerms(state,par);
k1Nonlinear = packet.rateProjection*k1Nonlinear;
[jacobian1,inputJacobian1] = localNonlinearJacobian(state,par);
j1Nonlinear = dt*(jacobian1*sensitivity+inputJacobian1);
j1Nonlinear = packet.rateProjection*j1Nonlinear;
j1Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (sensitivity+gamma*j1Nonlinear));
sensitivity1 = sensitivity+gamma*(j1Linear+j1Nonlinear);
k1Linear = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+gamma*k1Nonlinear));
stage1 = state+gamma*(k1Linear+k1Nonlinear);

k2Nonlinear = dt*localNonlinearTerms(stage1,par);
k2Nonlinear = packet.rateProjection*k2Nonlinear;
[jacobian2,inputJacobian2] = localNonlinearJacobian(stage1,par);
j2Nonlinear = dt*(jacobian2*sensitivity1+inputJacobian2);
j2Nonlinear = packet.rateProjection*j2Nonlinear;
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
k3Nonlinear = dt*localNonlinearTerms(stage2,par);
k3Nonlinear = packet.rateProjection*k3Nonlinear;
[jacobian3,inputJacobian3] = localNonlinearJacobian(stage2,par);
j3Nonlinear = dt*(jacobian3*sensitivity2+inputJacobian3);
j3Nonlinear = packet.rateProjection*j3Nonlinear;

sensitivityNext = sensitivity+(1-gamma)*(j1Linear+j2Nonlinear)+ ...
    gamma*(j2Linear+j3Nonlinear);
stateNext = state+(1-gamma)*(k1Linear+k2Nonlinear)+ ...
    gamma*(k2Linear+k3Nonlinear);
end

function value = localImplicitSolve(packet,rightHandSide)
%LOCALIMPLICITSOLVE Match ROMIntegrator's pivoted triangular solve order.
permuted = rightHandSide(packet.piv,:);
value = packet.Ufac\(packet.Lfac\permuted);
end

function nonlinear = localNonlinearTerms(state,par)
%LOCALNONLINEARTERMS Exact fixed-size Artola nonlinear residual.
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

function [stateJacobian,inputJacobian] = ...
        localNonlinearJacobian(state,par)
%LOCALNONLINEARJACOBIAN Exact fixed-size analytical residual derivatives.
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
%LOCALGAMMACONTRACT Contract the second tensor index with coordinates.
pageCoordinates = repmat(reshape(coordinates,[],1,1), ...
    1,1,size(tensor,3));
contractedPages = pagemtimes(tensor,pageCoordinates);
contracted = reshape(contractedPages,size(tensor,1),size(tensor,3));
end
