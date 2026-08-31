function output = scheduledReciprocalPlantIntervalCoreAudit( ...
        xFlexInitial,rigidInitial,hiddenInitial,command,gust, ...
        rigidIncrement,packet)
%SCHEDULEDRECIPROCALPLANTINTERVALCOREAUDIT Exact packed 14-substep P1 core.

xFlex = double(xFlexInitial(:));
rigid = double(rigidInitial(:));
hidden = double(hiddenInitial(:));
command = double(command(:));
gust = double(gust);
rigidIncrement = double(rigidIncrement(:));
assert(numel(xFlex)==74 && numel(rigid)==12 && ...
    numel(hidden)==546 && numel(command)==4 && isscalar(gust) && ...
    numel(rigidIncrement)==2 && ...
    all(isfinite([xFlex;rigid;hidden;command;gust;rigidIncrement])), ...
    'AeroFlex:sim:PlantIntervalInput', ...
    'The P1 plant interval input has invalid dimensions or values.');
substepCount = double(packet.substepCount);
assert(isscalar(substepCount) && isfinite(substepCount) && ...
    substepCount==round(substepCount) && substepCount>=1 && ...
    substepCount<=14, ...
    'AeroFlex:sim:PlantIntervalSubstepCount', ...
    'The P1 audit core requires between one and fourteen substeps.');
refinementUpdateLimit = 2;
if isfield(packet,'rigidWrenchRefinementUpdateLimit')
    refinementUpdateLimit = ...
        double(packet.rigidWrenchRefinementUpdateLimit);
end
assert(ismember(refinementUpdateLimit,[2,3]), ...
    'AeroFlex:sim:PlantIntervalRigidWrenchRefinementLimit', ...
    'The P1 rigid-wrench refinement update limit is invalid.');

wingTotal = packet.trimWing+command;
elevatorTotal = packet.trimElevator+rigidIncrement(1);
thrustTotal = packet.trimThrust+rigidIncrement(2);
assert(thrustTotal>=0, ...
    'AeroFlex:sim:PlantIntervalThrust', ...
    'The P1 plant interval requires nonnegative total thrust.');
sourceInput = [gust;wingTotal-packet.queryTrimWing; ...
    elevatorTotal-packet.queryTrimElevator; ...
    thrustTotal-packet.queryTrimThrust];
inputDeparture = sourceInput-packet.inputEquilibrium;
symmetryResidual = norm([sourceInput(2)-sourceInput(3); ...
    sourceInput(4)-sourceInput(5)],inf);
symmetryTolerance = 5e-13+1e-10*max(1,norm(sourceInput(2:5),inf));
assert(symmetryResidual<=symmetryTolerance, ...
    'AeroFlex:sim:PlantIntervalSymmetry', ...
    'The P1 symmetric-longitudinal input contains differential wing motion.');

stateHistory = zeros(632,substepCount+1);
clampHistory = zeros(6,substepCount);
wingWrenchHistory = zeros(6,substepCount);
tailWrenchHistory = zeros(6,substepCount);
finWrenchHistory = zeros(6,substepCount);
gravityWrenchHistory = zeros(6,substepCount);
thrustWrenchHistory = zeros(6,substepCount);
reciprocalWrenchHistory = zeros(6,substepCount);
totalWrenchHistory = zeros(6,substepCount);
sourceRootHistory = zeros(6,substepCount);
wrenchCorrectionHistory = zeros(6,substepCount);
visibleCorrectionHistory = zeros(69,substepCount);
fixedRootIncrementHistory = zeros(6,substepCount);
bodyAppliedWingWrenchHistory = zeros(6,substepCount);
rigidChartCorrectionHistory = zeros(3,substepCount);
rigidStateBeforeHistory = zeros(12,substepCount);
flightConditionHistory = zeros(9,substepCount);
sourceRatios = zeros(14,4,substepCount);
sourceAccepted = true(14,substepCount);
stepAccepted = false(1,substepCount);
domainRejected = false(1,substepCount);
domainReason = strings(1,substepCount);
domainDetails = cell(1,substepCount);
rigidResidual = zeros(1,substepCount);
refinementCount = zeros(1,substepCount);
stateHistory(:,1) = [xFlex;rigid;hidden];
for substep = 1:substepCount
    rigidStateBeforeHistory(:,substep) = rigid;
    flightConditionHistory(:,substep) = localFlightCondition(rigid, ...
        packet.fallbackSpeed,packet.fallbackAlpha);
    runtimeState = localRuntimeState(xFlex,rigid);
    clamp = localCenteredRootWrench(xFlex,wingTotal,gust, ...
        packet.recovery);
    [tailForce,tailMoment] = localTailLoads(rigid,elevatorTotal, ...
        packet.tail,packet.fallbackSpeed,packet.fallbackAlpha);
    gravityForce = localGravity(rigid,packet.mass,packet.gravity);
    thrustForce = [thrustTotal;0;0];
    thrustMoment = cross(packet.thrustArm,thrustForce);
    totalForce = clamp(1:3)+tailForce+gravityForce+thrustForce;
    totalMoment = clamp(4:6)+tailMoment+thrustMoment;

    visibleDeparture = runtimeState-packet.correctionStateEquilibrium;
    visibleCorrection = packet.visibleCorrectionMap*visibleDeparture+ ...
        packet.Avh*hidden+packet.deltaBv*inputDeparture;
    hiddenNext = packet.hiddenStateMap*visibleDeparture+ ...
        packet.Ahh*hidden+packet.Bh*inputDeparture;
    fixedRootIncrement = clamp-packet.fixedRootWingEquilibrium;
    totalForce = totalForce-fixedRootIncrement(1:3);
    totalMoment = totalMoment-fixedRootIncrement(4:6);
    sourceRoot = packet.sourceRootEquilibrium+ ...
        packet.sourceRootStateMap*[visibleDeparture;hidden]+ ...
        packet.sourceRootInputMap*inputDeparture;

    finiteCandidate = all(isfinite([visibleCorrection;hiddenNext;sourceRoot]));
    if finiteCandidate
        [domainOk,ratios,memberAccepted,currentDomainReason, ...
            currentDomainDetails] = localSourceDomain( ...
            runtimeState,inputDeparture,hidden,packet);
    else
        domainOk = false;
        ratios = zeros(14,4);
        memberAccepted = true(14,1);
        currentDomainReason = "nonfinite reciprocal prediction";
        currentDomainDetails = struct([]);
    end
    sourceRatios(:,:,substep) = ratios;
    sourceAccepted(:,substep) = memberAccepted;
    domainRejected(substep) = ~domainOk;
    domainReason(substep) = currentDomainReason;
    domainDetails{substep} = currentDomainDetails;
    accepted = false;
    wrenchCorrection = zeros(6,1);
    chartCorrection = zeros(3,1);
    realizedResidual = 0;
    refinements = 0;
    if domainOk
        operatorInput = [visibleDeparture;hidden;inputDeparture];
        wrenchCorrection = ...
            packet.rigidWrenchCommandOperator*operatorInput;
        desiredRigidCorrection = visibleCorrection(61:69);
        baselineNext = localRigidStep( ...
            rigid,totalForce,totalMoment,packet);
        for refinement = 0:refinementUpdateLimit
            correctedNext = localRigidStep(rigid, ...
                totalForce+wrenchCorrection(1:3), ...
                totalMoment+wrenchCorrection(4:6),packet);
            realizedRigidCorrection = ...
                localReciprocalRigidState(correctedNext)- ...
                localReciprocalRigidState(baselineNext);
            dynamicResidual = desiredRigidCorrection(1:6)- ...
                realizedRigidCorrection(1:6);
            realizedResidual = norm(dynamicResidual,inf);
            realizationTolerance = 5e-13+1e-8*max( ...
                norm(desiredRigidCorrection(1:6),inf),1e-10);
            if realizedResidual<=realizationTolerance || ...
                    refinement==refinementUpdateLimit
                break
            end
            wrenchCorrection = wrenchCorrection+ ...
                packet.dynamicWrenchMap\dynamicResidual;
            refinements = refinements+1;
        end
        accepted = realizedResidual<=realizationTolerance && ...
            all(isfinite([visibleCorrection;hiddenNext; ...
                wrenchCorrection;realizedRigidCorrection]));
        chartCorrection = desiredRigidCorrection(7:9)- ...
            realizedRigidCorrection(7:9);
    end
    if accepted
        totalForce = totalForce+wrenchCorrection(1:3);
        totalMoment = totalMoment+wrenchCorrection(4:6);
    else
        totalForce = totalForce+fixedRootIncrement(1:3);
        totalMoment = totalMoment+fixedRootIncrement(4:6);
        visibleCorrection = zeros(69,1);
        hiddenNext = hidden;
        wrenchCorrection = zeros(6,1);
        chartCorrection = zeros(3,1);
        if ~domainOk
            sourceRoot = zeros(6,1);
        end
    end
    if accepted
        bodyAppliedWingWrench = ...
            packet.fixedRootWingEquilibrium+wrenchCorrection;
    else
        bodyAppliedWingWrench = clamp+wrenchCorrection;
    end

    xFlex = localNativeStep(xFlex,wingTotal,gust,packet.native)+ ...
        packet.native.correctionToNative*visibleCorrection;
    rigid = localRigidStep(rigid,totalForce,totalMoment,packet);
    rigid(7:9) = rigid(7:9)+chartCorrection;
    hidden = hiddenNext;
    assert(all(isfinite([xFlex;rigid;hidden])), ...
        'AeroFlex:sim:PlantIntervalOutputFinite', ...
        'The P1 plant interval produced a nonfinite state.');

    stateHistory(:,substep+1) = [xFlex;rigid;hidden];
    clampHistory(:,substep) = clamp;
    wingWrenchHistory(:,substep) = clamp;
    tailWrenchHistory(:,substep) = [tailForce;tailMoment];
    gravityWrenchHistory(:,substep) = [gravityForce;zeros(3,1)];
    thrustWrenchHistory(:,substep) = [thrustForce;thrustMoment];
    reciprocalWrenchHistory(:,substep) = wrenchCorrection;
    totalWrenchHistory(:,substep) = [totalForce;totalMoment];
    sourceRootHistory(:,substep) = sourceRoot;
    wrenchCorrectionHistory(:,substep) = wrenchCorrection;
    visibleCorrectionHistory(:,substep) = visibleCorrection;
    fixedRootIncrementHistory(:,substep) = fixedRootIncrement;
    bodyAppliedWingWrenchHistory(:,substep) = bodyAppliedWingWrench;
    rigidChartCorrectionHistory(:,substep) = chartCorrection;
    stepAccepted(substep) = accepted;
    rigidResidual(substep) = realizedResidual;
    refinementCount(substep) = refinements;
end
output = struct( ...
    'stateHistory',stateHistory, ...
    'clampHistory',clampHistory, ...
    'wingWrenchHistory',wingWrenchHistory, ...
    'tailWrenchHistory',tailWrenchHistory, ...
    'finWrenchHistory',finWrenchHistory, ...
    'gravityWrenchHistory',gravityWrenchHistory, ...
    'thrustWrenchHistory',thrustWrenchHistory, ...
    'reciprocalWrenchHistory',reciprocalWrenchHistory, ...
    'totalWrenchHistory',totalWrenchHistory, ...
    'sourceRootHistory',sourceRootHistory, ...
    'wrenchCorrectionHistory',wrenchCorrectionHistory, ...
    'visibleCorrectionHistory',visibleCorrectionHistory, ...
    'fixedRootIncrementHistory',fixedRootIncrementHistory, ...
    'bodyAppliedWingWrenchHistory',bodyAppliedWingWrenchHistory, ...
    'rigidChartCorrectionHistory',rigidChartCorrectionHistory, ...
    'rigidStateBeforeHistory',rigidStateBeforeHistory, ...
    'flightConditionHistory',flightConditionHistory, ...
    'sourceRatios',sourceRatios, ...
    'sourceAccepted',sourceAccepted, ...
    'stepAccepted',stepAccepted, ...
    'domainRejected',domainRejected, ...
    'domainReason',domainReason, ...
    'domainDetails',{domainDetails}, ...
    'rigidResidual',rigidResidual, ...
    'refinementCount',refinementCount, ...
    'inputDeparture',inputDeparture, ...
    'wingTotal',wingTotal, ...
    'elevatorTotal',elevatorTotal, ...
    'thrustTotal',thrustTotal);
end

function condition = localFlightCondition(rigid,fallbackSpeed,fallbackAlpha)
velocity = rigid(4:6);
speed = norm(velocity);
if speed<1e-9
    speed = fallbackSpeed;
    alpha = fallbackAlpha;
    beta = 0;
else
    alpha = atan2(velocity(3),velocity(1));
    beta = asin(max(-1,min(1,velocity(2)/speed)));
end
condition = [speed;alpha;beta;rigid(7:9);rigid(10:12)];
end

function state = localRuntimeState(xFlex,rigid)
state = [xFlex(1:10);xFlex(11:20);xFlex(32:71); ...
    rigid(4:6);rigid(10:12);rigid(7:9);xFlex(21:31);xFlex(72:74)];
end

function state = localReciprocalRigidState(rigid)
state = [rigid(4:6);rigid(10:12);rigid(7:9)];
end

function [accepted,ratios,memberAccepted,reason,details] = localSourceDomain( ...
        runtimeState,inputDeparture,hidden,packet)
ratios = zeros(14,4);
memberAccepted = true(14,1);
accepted = true;
reason = "";
details = repmat(struct('sourceId',"",'inputRatio',0, ...
    'visibleRatio',0,'hiddenRatio',0,'pulseHiddenRatio',0, ...
    'memoryOwner',"one_step_pulse",'rootWrenchRatio',0, ...
    'accepted',true),14,1);
for memberIndex = 1:14
    if packet.weights(memberIndex)<=0
        continue
    end
    member = packet.members(memberIndex);
    localHidden = (memberIndex-1)*39+(1:39);
    departure = member.runtimeFromQueryState*runtimeState- ...
        member.correctionStateEquilibrium;
    runtimeDeviation = departure-member.queryRuntimeDeparture;
    queryInput = inputDeparture+packet.inputEquilibrium;
    actualInput = [queryInput(1); ...
        packet.queryTrimWing+queryInput(2:5); ...
        packet.queryTrimElevator+queryInput(6); ...
        packet.queryTrimThrust+queryInput(7)];
    sourceInputDeparture = [actualInput(1); ...
        actualInput(2:5)-member.trimWing; ...
        actualInput(6)-member.trimElevator; ...
        actualInput(7)-member.trimThrust]-member.inputEquilibrium;
    inputDeviation = sourceInputDeparture-member.queryInputDeparture;
    hiddenDeviation = hidden(localHidden)-member.queryHiddenEquilibrium;
    visible = runtimeDeviation(member.sourceVisible);
    rootWrench = member.rootEquilibrium+member.rootStateMap* ...
        [departure;hiddenDeviation+member.queryHiddenEquilibrium]+ ...
        member.rootInputMap*sourceInputDeparture;
    rootWrenchDeviation = rootWrench-member.queryRootWrench;
    symmetricInput = [inputDeviation(1);mean(inputDeviation(2:3)); ...
        mean(inputDeviation(4:5));inputDeviation(6:7)];
    inputSymmetry = norm([inputDeviation(2)-inputDeviation(3); ...
        inputDeviation(4)-inputDeviation(5)],inf);
    symmetryTolerance = 100*eps(max(1,norm(inputDeviation,inf)));
    inputRatio = norm(symmetricInput./member.inputBounds,inf);
    visibleRatio = norm(visible,inf)/member.visibleBound;
    pulseHiddenRatio = norm(hiddenDeviation,inf)/member.hiddenBound;
    hiddenRatio = pulseHiddenRatio;
    memoryOwner = "one_step_pulse";
    if member.memoryCertificateEnabled
        assert(isequal(size(member.memoryLyapunovFactor),[39,39]) && ...
            isfinite(member.memoryInvariantRadius) && ...
            member.memoryInvariantRadius>0, ...
            'AeroFlex:sim:PlantIntervalMemoryCertificate', ...
            'The interval recurrent-memory certificate is invalid.');
        hiddenRatio = norm(member.memoryLyapunovFactor'*hiddenDeviation,2)/ ...
            member.memoryInvariantRadius;
        memoryOwner = "discrete_lyapunov_invariant";
    end
    rootRatio = norm(rootWrenchDeviation,inf)/member.rootWrenchBound;
    ratios(memberIndex,:) = ...
        [inputRatio,visibleRatio,hiddenRatio,rootRatio];
    tolerance = 100*eps(max(1,max(ratios(memberIndex,:))));
    memberAccepted(memberIndex) = inputSymmetry<=symmetryTolerance && ...
        max(ratios(memberIndex,:))<=1+tolerance;
    details(memberIndex) = struct('sourceId',string(member.sourceId), ...
        'inputRatio',inputRatio,'visibleRatio',visibleRatio, ...
        'hiddenRatio',hiddenRatio,'pulseHiddenRatio',pulseHiddenRatio, ...
        'memoryOwner',memoryOwner,'rootWrenchRatio',rootRatio, ...
        'accepted',memberAccepted(memberIndex));
    if ~memberAccepted(memberIndex)
        accepted = false;
        reason = "source_domain_"+string(member.sourceId);
        break
    end
end
end

function wrench = localCenteredRootWrench(state,control,gust,recovery)
dx = state-recovery.referenceState;
du = control-recovery.referenceControl;
wrench = recovery.anchor;
for index = 1:numel(recovery.members)
    member = recovery.members(index);
    vertexState = member.xEq+dx;
    vertexControl = member.uEq+du;
    current = localVertexReaction(vertexState,vertexControl,gust,member);
    wrench = wrench+member.weight*(current-member.referenceWrench);
end
end

function wrench = localVertexReaction(state,control,gust,member)
nonlinear = localNonlinear(state,control,gust,member.par);
q1raw = member.Lq1*state+nonlinear(1:10);
wrench = member.rootWrenchOperator*q1raw;
end

function [force,moment] = localTailLoads( ...
        rigid,elevator,tail,fallbackSpeed,fallbackAlpha)
velocity = rigid(4:6);
omega = rigid(10:12);
U = norm(velocity);
if U<1e-9
    U = fallbackSpeed;
    alpha = fallbackAlpha;
else
    alpha = atan2(velocity(3),velocity(1));
end
Ueff = max(abs(U),1e-9);
qbar = 0.5*tail.rho*Ueff^2;
vRef = [Ueff*cos(alpha);0;Ueff*sin(alpha)];
vTail = vRef+cross([0;omega(2);0],tail.arm);
alphaTail = atan2(vTail(3),vTail(1));
alphaEffective = alphaTail-tail.incidence-tail.etaDelta*elevator;
cl = tail.clAlpha*alphaEffective;
cd = tail.cd0+cl^2/(pi*tail.efficiency*tail.aspectRatio);
lift = qbar*tail.area*cl;
drag = qbar*tail.area*cd;
force = [-drag;0;-lift];
momentAc = [0;qbar*tail.area*tail.chord* ...
    (tail.cm0+tail.cmDelta*elevator);0];
moment = momentAc+cross(tail.arm,force);
end

function force = localGravity(rigid,mass,gravity)
dcm = localDcm(rigid(7:9));
force = dcm.'*[0;0;mass*gravity];
end

function next = localRigidStep(rigid,force,moment,packet)
dt = packet.dt;
k1 = localRigidDerivative(rigid,force,moment,packet);
k2 = localRigidDerivative(rigid+0.5*dt*k1,force,moment,packet);
k3 = localRigidDerivative(rigid+0.5*dt*k2,force,moment,packet);
k4 = localRigidDerivative(rigid+dt*k3,force,moment,packet);
next = rigid+(dt/6)*(k1+2*k2+2*k3+k4);
end

function derivative = localRigidDerivative(rigid,force,moment,packet)
velocity = rigid(4:6);
euler = rigid(7:9);
omega = rigid(10:12);
positionRate = localDcm(euler)*velocity;
velocityRate = force/packet.mass-cross(omega,velocity);
omegaRate = packet.inertia\ ...
    (moment-cross(omega,packet.inertia*omega));
eulerRate = localEulerRates(euler,omega);
derivative = [positionRate;velocityRate;eulerRate;omegaRate];
end

function dcm = localDcm(euler)
phi=euler(1); theta=euler(2); psi=euler(3);
cphi=cos(phi); sphi=sin(phi); cth=cos(theta); sth=sin(theta);
cpsi=cos(psi); spsi=sin(psi);
dcm = [cth*cpsi, sphi*sth*cpsi-cphi*spsi, ...
    cphi*sth*cpsi+sphi*spsi; ...
    cth*spsi, sphi*sth*spsi+cphi*cpsi, ...
    cphi*sth*spsi-sphi*cpsi; ...
    -sth, sphi*cth, cphi*cth];
end

function rate = localEulerRates(euler,omega)
phi=euler(1); theta=euler(2); ctheta=cos(theta);
assert(abs(ctheta)>=1e-8, ...
    'AeroFlex:sim:PlantIntervalEulerSingularity', ...
    'P1 does not admit a near-singular Euler chart.');
transform = [1,sin(phi)*tan(theta),cos(phi)*tan(theta); ...
    0,cos(phi),-sin(phi);0,sin(phi)/ctheta,cos(phi)/ctheta];
rate = transform*omega;
end

function stateNext = localNativeStep(state,control,gust,packet)
par = packet.par;
dt=packet.dt; gamma=packet.gamma; delta=packet.delta;
k1N = packet.rateProjection*(dt*localNonlinear(state,control,gust,par));
k1L = localImplicitSolve(packet,dt*packet.Ldyn*(state+gamma*k1N));
stage1 = state+gamma*(k1L+k1N);
k2N = packet.rateProjection*(dt*localNonlinear(stage1,control,gust,par));
k2L = localImplicitSolve(packet,dt*packet.Ldyn* ...
    (state+(1-gamma)*k1L+delta*k1N+(1-delta)*k2N));
stage2 = state+(1-gamma)*k1L+gamma*k2L+delta*k1N+(1-delta)*k2N;
k3N = packet.rateProjection*(dt*localNonlinear(stage2,control,gust,par));
stateNext = state+(1-gamma)*(k1L+k2N)+gamma*(k2L+k3N);
end

function value = localImplicitSolve(packet,rightHandSide)
value = packet.Ufac\(packet.Lfac\rightHandSide(packet.piv,:));
end

function nonlinear = localNonlinear(state,control,gust,par)
q1=state(1:10); q2=state(11:20); qxi=state(21:31);
nonlinear=zeros(74,1);
gamma1q1=squeeze(pagemtimes(par.Gamma1,'none',q1,'none'));
gamma2q2=squeeze(pagemtimes(par.Gamma2,'none',q2,'none'));
gammaGqxi=squeeze(pagemtimes(par.Gamma_g,'none',qxi,'none'));
gammaXiqxi=squeeze(pagemtimes(par.Gamma_xi,'none',qxi,'none'));
gamma2TransposeQ1=squeeze(pagemtimes( ...
    par.Gamma2,'transpose',q1,'none'));
nonlinear(1:10)=-(gamma1q1*q1)-(gamma2q2*q2)+ ...
    gammaGqxi*qxi+par.Fscale*par.forces_0+par.N_Thrust;
nonlinear(11:20)=gamma2TransposeQ1*q2;
nonlinear(21:31)=gammaXiqxi*q1;
nonlinear(32:71)=par.Bw*gust+(1/par.t_inf)* ...
    par.Bdel*control(1:2)+par.Bddel*control(3:4);
nonlinear(1:10)=nonlinear(1:10)+par.scaleAero*par.Dw*gust+ ...
    par.Ddel*control(1:2)+par.t_inf*par.Dddel*control(3:4);
nonlinear=nonlinear+par.affineOffset;
end
