function trim = TrimRBwFlex(cfg, beam, aero, base)
%TRIMRBWFLEX Coupled longitudinal trim for the Pazy flexible-aircraft ROM.
%
% Trim variables:
%   z(1) = alpha      [rad]  flight/rigid-body angle of attack
%   z(2) = deltaWing  [rad]  symmetric wing control: [d; d; 0; 0]
%   z(3) = deltaElev  [rad]  elevator/tail command used by the local tail model
%   z(4) = thrust     [N]    rigid-body +X thrust
%
% Residual:
%   R = [Fx; Fz; My] in body axes.
%
% The flexible inner solve is the clamped-frame solve:
%   Pz*q1dot -> 0
% Root loads are recovered afterward from the raw RHS with Pr*q1dot.  This
% separation is intentional; do not use the projected RHS for load recovery.

% -------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------
cfg = localFillTrimDefaults(cfg);

idx0 = AeroFlex.core.buildIndexStruct(beam.Nm, aero.Na);
Pz0  = beam.Pz;
Pr0  = beam.Pr;

useSched = isfield(cfg,'library') && isstruct(cfg.library) && ...
           isfield(cfg.library,'enable') && logical(cfg.library.enable) && ...
           isfield(cfg.library,'path') && exist(cfg.library.path,'file') == 2;

if useSched
    ROMlib = AeroFlex.sched.loadLibrary(cfg.library.path);
    if ~isfield(ROMlib,'compatibleCoordinates') || ~logical(ROMlib.compatibleCoordinates)
        error('TrimRBwFlex:RawLibrary', ...
            'Scheduled coupled trim requires the Stage-2 compatible ROM library.');
    end
    cfg.library.noExtrapolate = true;

    % Finite differences call evalLibrary repeatedly.  Keep interpolation
    % prints off unless the schedule itself is being debugged.
    if ~localGetNested(cfg, {'debug','trimSchedule'}, false)
        cfg.library.debug = false;
    end
else
    ROMlib = [];
end

fixedExactOwner = localFixedExactPackageTrimOwner(cfg, ROMlib, useSched);

rigid = RigidBody.methods.paramsRigid_PazyUAV();
if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'mass') && ~isempty(cfg.rigidEOMset.mass)
    mTot = cfg.rigidEOMset.mass;
elseif isfield(rigid,'mass') && ~isempty(rigid.mass)
    mTot = rigid.mass;
else
    error('TrimRBwFlex:MissingMass', ...
        'Rigid-body mass unavailable. Set cfg.rigidEOMset.mass.');
end

if isfield(cfg,'rigidEOMset') && ...
        isfield(cfg.rigidEOMset,'massOwnership') && ...
        string(cfg.rigidEOMset.massOwnership.owner)== ...
        "nonwing_source_owned_wing_reaction"
    expected = rigid.nonwing;
    assert(abs(mTot-expected.mass)<=1e-14 && ...
        isfield(cfg.rigidEOMset,'I_B') && ...
        norm(cfg.rigidEOMset.I_B-expected.J,'fro')<=1e-14 && ...
        isfield(cfg.rigidEOMset,'rWingRoot_B') && ...
        norm(cfg.rigidEOMset.rWingRoot_B(:)- ...
            expected.wingRootFromCM,inf)<=1e-14 && ...
        isfield(cfg,'tail') && isfield(cfg.tail,'r_B') && ...
        norm(cfg.tail.r_B(:)-expected.tailArmFromCM,inf)<=1e-14, ...
        'TrimRBwFlex:MixedRigidMassOwnership', ...
        ['The V17A non-wing trim mass, inertia, wing-root arm, and tail ', ...
         'arm must come from one owner.']);
end

g0 = localGetNested(cfg, {'flight','g'}, 9.807);
W  = mTot*g0;

bounds = localTrimBounds(cfg, ROMlib, useSched, W, fixedExactOwner);
z      = localInitialTrimVariables(cfg, bounds);
opts   = localTrimOptions(cfg, W, bounds, z);

[xGuess, ~] = localInitialState(cfg, beam, aero, base, idx0, ROMlib, useSched, z(1));

history = localInitHistory(opts.maxIt);
lambda = opts.lambda0;
converged = false;
stallCount = 0;
log_trim = [];

if opts.debug
    fprintf('\n[TrimRBwFlex] Scheduled ROM: %d\n', useSched);
    fprintf('[TrimRBwFlex] Bounds: alpha=[%.3f %.3f] deg, wing=[%.3f %.3f] deg, elev=[%.3f %.3f] deg, T=[%.3f %.3f] N\n', ...
        rad2deg(bounds.lb(1)), rad2deg(bounds.ub(1)), ...
        rad2deg(bounds.lb(2)), rad2deg(bounds.ub(2)), ...
        rad2deg(bounds.lb(3)), rad2deg(bounds.ub(3)), ...
        bounds.lb(4), bounds.ub(4));
    fprintf('[TrimRBwFlex] Targets: |Fx|<=%.1e N, |Fz|<=%.1e N, |My|<=%.1e N*m, |Pz*q1dot|<=%.1e\n', ...
        opts.tolPhys(1), opts.tolPhys(2), opts.tolPhys(3), opts.tolSS);
end

% -------------------------------------------------------------------------
% Trust-region Levenberg iteration.
%
% The main difference from earlier attempts is that a trial point is accepted
% only if the flexible inner solve remains valid.  Force-only descent is not
% allowed to accept a point with a bad inner residual.
% -------------------------------------------------------------------------
for it = 1:opts.maxIt

    [R, xSS, log_trim, diag0] = localResidual(z, xGuess, cfg, beam, aero, base, ...
        ROMlib, useSched, idx0, Pz0, Pr0, W, opts);

    if ~all(isfinite(R))
        error('TrimRBwFlex:NonFiniteResidual', ...
            'Residual became non-finite at alpha=%.6g deg.', rad2deg(z(1)));
    end

    if ~localInnerAcceptable(diag0, opts)
        warning('TrimRBwFlex:CurrentInnerNotConverged', ...
            'Current iterate has inner residual %.3e. Outer residual is not reliable.', ...
            localGet(diag0.inner,'residual',NaN));
    end

    % Final axial cleanup, but only when Fz/My are already good.
    [didPolish, zPolish, RPolish, dTpolish] = localTryThrustPolish(z, R, bounds, cfg, opts);
    if didPolish
        z = zPolish;
        R = RPolish;
    end

    [resScale, phase, forceNorm, momentScale] = localTrimScale(R, opts);
    [phi0, rScaled] = localMerit(z, R, resScale, opts);
    rNorm = norm(rScaled,2);
    
    % Endgame block polish: once the force/moment balance is close, solve the
    % [Fz; My] block using wing/elevator, then clean Fx with thrust.  This avoids
    % wasting many LM iterations on a nearly singular 3-by-4 system.
    didBlockPolish = false;
    if opts.endgamePolish && phase == "full" && rNorm <= opts.endgamePolishNorm
        [didBlockPolish, zBlk, RBlk, xBlk, logBlk, diagBlk, infoBlk] = ...
            localTryEndgameBlockPolish(z, R, xSS, log_trim, diag0, bounds, ...
            cfg, beam, aero, base, ROMlib, useSched, idx0, Pz0, Pr0, W, opts);
    
        if didBlockPolish
            z = zBlk;
            R = RBlk;
            xSS = xBlk;
            log_trim = logBlk;
            diag0 = diagBlk;
    
            [didAxialPolish, zAx, RAx, dTax] = localTryThrustPolish(z, R, bounds, cfg, opts);
            if didAxialPolish
                z = zAx;
                R = RAx;
                dTpolish = dTax;
                didPolish = true;
            end
    
            [resScale, phase, forceNorm, momentScale] = localTrimScale(R, opts);
            [phi0, rScaled] = localMerit(z, R, resScale, opts);
            rNorm = norm(rScaled,2);
    
            if opts.debug
                fprintf('       endgame polish: %s | dz=[%.3e %.3e %.3e %.3e]\n', ...
                    string(infoBlk.reason), infoBlk.dz(1), infoBlk.dz(2), ...
                    infoBlk.dz(3), infoBlk.dz(4));
            end
        end
    end
    
    history.z(:,it) = z;
    history.R(:,it) = R;
    history.RscaledNorm(it) = rNorm;
    history.cost(it) = phi0;
    history.lambda(it) = lambda;
    history.phase(it) = phase;
    history.innerResidual(it) = localGet(diag0.inner,'residual',NaN);
    history.innerStop(it) = string(localGet(diag0.inner,'stopReason','unknown'));
    history.thrustPolish(it) = didPolish;

    if opts.debug
        fprintf(['----------------------------------------------------------\n', ...
                 'It %2d | alpha=%8.4f deg | wing=%8.4f deg | elev=%8.4f deg | T=%9.4f N\n', ...
                 '       R=[%+.3e %+.3e %+.3e] | |R_s|=%.3e | phase=%s | Mscale=%.3e\n'], ...
            it, rad2deg(z(1)), rad2deg(z(2)), rad2deg(z(3)), z(4), ...
            R(1), R(2), R(3), rNorm, string(phase), momentScale);
        fprintf('       inner: %s | iter=%d | res=%.3e | initial=%.3e | best=%.3e | lambda=%.3e\n', ...
            string(localGet(diag0.inner,'stopReason','unknown')), ...
            localGet(diag0.inner,'iterations',-1), ...
            localGet(diag0.inner,'residual',NaN), ...
            localGet(diag0.inner,'initialResidual',NaN), ...
            localGet(diag0.inner,'bestResidual',NaN), lambda);
        if didPolish
            fprintf('       thrust polish: dT=%+.6e N -> T=%.6e N\n', dTpolish, z(4));
        end
    end

    physOK = all(abs(R(:)) <= opts.tolPhys(:));
    if localInnerAcceptable(diag0, opts) && physOK
        converged = true;
        xGuess = xSS;
        break
    end

    [J, jacInfo] = localFiniteDifferenceJacobian(z, xSS, R, bounds, cfg, beam, aero, base, ...
        ROMlib, useSched, idx0, Pz0, Pr0, W, opts);

    if ~all(isfinite(J(:)))
        error('TrimRBwFlex:NonFiniteJacobian', ...
            'Trim Jacobian contains Inf/NaN.');
    end

    JsDebug = diag(1./resScale) * J * diag(opts.varScale);
    history.condJ(it) = cond(JsDebug);

    if opts.debug && (it == 1 || mod(it,opts.debugEvery) == 0 || rNorm < opts.fullDebugNorm)
        fprintf('       J scaled rank=%d cond=%.3e colNorms=[%s]\n', ...
            rank(JsDebug), cond(JsDebug), num2str(vecnorm(JsDebug), ' %.3e'));
        fprintf('       dR/dalpha=[%+.3e %+.3e %+.3e], dR/dwing=[%+.3e %+.3e %+.3e]\n', ...
            J(1,1), J(2,1), J(3,1), J(1,2), J(2,2), J(3,2));
        fprintf('       dR/delev =[%+.3e %+.3e %+.3e], dR/dT   =[%+.3e %+.3e %+.3e]\n', ...
            J(1,3), J(2,3), J(3,3), J(1,4), J(2,4), J(3,4));
    end

    % Residual-dependent trust region.  This is not a fixed step size:
    % the allowed step shrinks as the residual approaches the target.
    maxStepIt = localActiveMaxStep(opts, R, rNorm, phase);
    
    dzGN = localLevenbergStep(J, R, z, resScale, lambda, opts);
    dzGN = localLimitStep(dzGN, maxStepIt);
    dzGN = localRespectBounds(z, dzGN, bounds);
    
    [accepted, best, stepFactor] = localLineSearch(z, dzGN, phi0, resScale, ...
        xSS, cfg, beam, aero, base, ROMlib, useSched, idx0, Pz0, Pr0, W, opts, bounds);
    
    usedFallback = false;
    
    if ~accepted
        dzSD = localSteepestDescentStep(J, R, z, resScale, opts);
        dzSD = localLimitStep(dzSD, 0.50*maxStepIt);
        dzSD = localRespectBounds(z, dzSD, bounds);
    
        [accepted, best, stepFactor] = localLineSearch(z, dzSD, phi0, resScale, ...
            xSS, cfg, beam, aero, base, ROMlib, useSched, idx0, Pz0, Pr0, W, opts, bounds);
    
        usedFallback = accepted;
    end
    
    if accepted
        z = best.z;
        xGuess = best.x;
        log_trim = best.log;
        history.lineStep(it) = stepFactor;
    
        actualReduction = (phi0 - best.phi)/max(phi0, eps);
    
        % Do not blindly drive lambda to 1e-8.  If line-search only accepts a
        % small step, keep or increase damping.  If a full step gives meaningful
        % reduction, reduce damping.
        if stepFactor >= 0.75 && actualReduction > 0.25
            lambda = max(lambda/opts.lambdaDown, opts.lambdaMin);
        elseif stepFactor < 0.25 || actualReduction < 0.05
            lambda = min(lambda*opts.lambdaUpSoft, opts.lambdaMax);
        end
    
        stallCount = 0;
    
        if opts.debug && usedFallback
            fprintf('       accepted steepest-descent fallback, stepFactor=%.3f\n', stepFactor);
        end
    else
        lambda = min(lambda*opts.lambdaUp, opts.lambdaMax);
        history.lineStep(it) = 0;
        stallCount = stallCount + 1;
    
        if opts.debug
            fprintf('[TrimRBwFlex] No feasible descent step. lambda -> %.3e\n', lambda);
            fprintf('       rejected GN dz=[%+.3e %+.3e %+.3e %+.3e] = [%.4f %.4f %.4f %.4f] deg/N\n', ...
                dzGN(1), dzGN(2), dzGN(3), dzGN(4), ...
                rad2deg(dzGN(1)), rad2deg(dzGN(2)), rad2deg(dzGN(3)), dzGN(4));
        end
    end

    if lambda >= opts.lambdaMax || stallCount >= opts.maxStall
        if opts.debug
            fprintf('[TrimRBwFlex] Stopping outer loop: lambda=%.3e, stallCount=%d.\n', ...
                lambda, stallCount);
        end
        xGuess = xSS;
        break
    end
end

% -------------------------------------------------------------------------
% Final residual and output
% -------------------------------------------------------------------------
[Rfinal, xFinal, log_trim, diagFinal] = localResidual(z, xGuess, cfg, beam, aero, base, ...
    ROMlib, useSched, idx0, Pz0, Pr0, W, opts);

[didBlockFinal, zBlock, RBlock, xBlock, logBlock, diagBlock, infoBlock] = ...
    localTryEndgameBlockPolish(z, Rfinal, xFinal, log_trim, diagFinal, bounds, ...
    cfg, beam, aero, base, ROMlib, useSched, idx0, Pz0, Pr0, W, opts);

if didBlockFinal
    z = zBlock;
    Rfinal = RBlock;
    xFinal = xBlock;
    log_trim = logBlock;
    diagFinal = diagBlock;

    if opts.debug
        fprintf('[TrimRBwFlex] final endgame polish: %s | dz=[%.3e %.3e %.3e %.3e]\n', ...
            string(infoBlock.reason), infoBlock.dz(1), infoBlock.dz(2), ...
            infoBlock.dz(3), infoBlock.dz(4));
    end
end

[didPolishFinal, zPolish, RPolish, dTpolish] = localTryThrustPolish(z, Rfinal, bounds, cfg, opts);
if didPolishFinal
    z = zPolish;
    Rfinal = RPolish;
    if opts.debug
        fprintf('[TrimRBwFlex] final thrust polish: dT=%+.6e N | T=%.6e N | R=[%+.3e %+.3e %+.3e]\n', ...
            dTpolish, z(4), Rfinal(1), Rfinal(2), Rfinal(3));
    end
end

[resScaleFinal, ~, ~, ~] = localTrimScale(Rfinal, opts);
rFinalNorm = norm(Rfinal./resScaleFinal, 2);
physFinalOK = all(abs(Rfinal(:)) <= opts.tolPhys(:));
innerFinalOK = localInnerAcceptable(diagFinal, opts);

if ~converged && physFinalOK && innerFinalOK
    converged = true;
end

if ~converged
    warning('TrimRBwFlex:NotConverged', ...
        'Trim did not converge. Final residual = [%+.3e %+.3e %+.3e], inner residual = %.3e.', ...
        Rfinal(1), Rfinal(2), Rfinal(3), localGet(diagFinal.inner,'residual',NaN));
end

trim = localBuildOutput(z, xFinal, log_trim, Rfinal, rFinalNorm, ...
    converged, bounds, history, useSched, opts, diagFinal, resScaleFinal);

localPrintFinalDiagnostics(trim, diagFinal, opts);

end

% =============================================================================
% Residual and model construction
% =============================================================================
function [R, xSS, log_trim, diagOut] = localResidual(z, xGuess, cfg, beam, aero, base, ...
    ROMlib, useSched, idxIn, PzIn, PrIn, W, opts) %#ok<INUSD>

alpha = z(1);
deltaWingScalar = z(2);
deltaElev = z(3);
T = z(4);
deltaWing = [deltaWingScalar; deltaWingScalar; 0; 0];

try
    [model, baseUse, schedUse] = localMakeTrimModel(cfg, beam, aero, base, ROMlib, useSched, alpha);
catch ME
    if opts.debug
        fprintf('[TrimRBwFlex] model build failed at alpha=%.4f deg: %s\n', ...
            rad2deg(alpha), ME.message);
    end
    R = [NaN; NaN; NaN];
    xSS = xGuess;
    log_trim = [];
    diagOut = struct('error', ME.message);
    return
end

idx = model.idx;
[PzUse, PrUse, phi1_sA_Use] = localActiveBeamMaps(beam, schedUse, useSched, PzIn, PrIn);
x0 = localPrepareInitialState(xGuess, model, baseUse, schedUse, idx, cfg, useSched, alpha);

[xSS, model, log_trim, innerInfo] = localInnerClampedSteady(x0, deltaWing, model, idx, PzUse, opts);

if ~all(isfinite(xSS))
    R = [NaN; NaN; NaN];
    diagOut = struct('inner', innerInfo, 'sched', schedUse, 'FlexPlant', model);
    return
end

pc = model.parConst;
pc.u_ctrl = deltaWing;
pc.gust = localZeroGust(pc);
if ~localGetNested(cfg, {'trim','thrustActsOnWing'}, false)
    pc.N_Thrust = zeros(numel(idx.q1),1);
else
    pc.N_Thrust = phi1_sA_Use.' * [T;0;0;0;0;0];
end

% Raw RHS for root reaction.  Do not use model.Ldyn here.
xdot = AeroFlex.sim.nonlinear_terms(xSS, pc, idx) + model.L*xSS;
q1dot = xdot(idx.q1);

if isstruct(schedUse) && isfield(schedUse,'equilibriumCentered') && ...
        isfield(schedUse.equilibriumCentered,'enabled') && ...
        schedUse.equilibriumCentered.enabled
    Clamp6 = AeroFlex.sched.recoverCenteredRootWrench( ...
        schedUse,xSS,deltaWing);
else
    Clamp6 = AeroFlex.beam.recoverRootWrench( ...
        phi1_sA_Use,PrUse,PrUse*q1dot);
end
Clamp6 = localApplyWingSymmetry(Clamp6, cfg);

Fwing_B = Clamp6(1:3);
Mwing_B = Clamp6(4:6);

rA_CG_B = localGetNested(cfg, {'rigidEOMset','rWingRoot_B'}, ...
          localGetNested(cfg, {'geom','rWingRoot_B'}, [0;0;0]));
Mwing_B = Mwing_B + cross(rA_CG_B(:), Fwing_B);

U = localGetNested(cfg, {'flight','U_inf'}, 0);
[Ftail_B, Mtail_B] = localTailAeroForceMoment(U, alpha, deltaElev, cfg);
[Ffin_B,  Mfin_B ] = localFinAeroForceMoment(U, 0.0, 0.0, cfg);
Fgrav_B = localGravityBody(W, alpha);

rT = localGetNested(cfg, {'rigidEOMset','rThrust_B'}, [0;0;0]);
Fthrust_B = [T;0;0];
Mthrust_B = cross(rT(:), Fthrust_B);

Ftot_B = Fwing_B + Ftail_B + Ffin_B + Fgrav_B + Fthrust_B;
Mtot_B = Mwing_B + Mtail_B + Mfin_B + Mthrust_B;

R = [Ftot_B(1); Ftot_B(3); Mtot_B(2)];

diagOut = struct();
diagOut.Clamp6 = Clamp6;
diagOut.Fwing_B = Fwing_B;
diagOut.Mwing_B = Mwing_B;
diagOut.Ftail_B = Ftail_B;
diagOut.Mtail_B = Mtail_B;
diagOut.Ffin_B = Ffin_B;
diagOut.Mfin_B = Mfin_B;
diagOut.Fgrav_B = Fgrav_B;
diagOut.Fthrust_B = Fthrust_B;
diagOut.Mthrust_B = Mthrust_B;
diagOut.Ftot_B = Ftot_B;
diagOut.Mtot_B = Mtot_B;
diagOut.FpreT_B = Fwing_B + Ftail_B + Ffin_B + Fgrav_B;
diagOut.MpreT_B = Mwing_B + Mtail_B + Mfin_B;
diagOut.inner = innerInfo;
diagOut.sched = schedUse;
diagOut.FlexPlant = model;
diagOut.model = model;
diagOut.idx = idx;
diagOut.Pz = PzUse;
diagOut.Pr = PrUse;
diagOut.phi1_sA = phi1_sA_Use;
diagOut.raw_q1dot = q1dot;
diagOut.raw_xdot = xdot;
diagOut.Tused = T;
end

function [model, baseUse, schedUse] = localMakeTrimModel(cfg, beam, aero, base, ROMlib, useSched, alpha)
baseUse = base;
schedUse = [];

if useSched
    mu = localTrimSchedulePoint(cfg, alpha);
    schedUse = localEvalSchedule(ROMlib,mu,cfg.library);
    localRequireExactSchedule(cfg, schedUse, mu);
    baseUse = schedUse.base;

    model = AeroFlex.sim.ROMIntegrator(cfg, beam, aero, baseUse);
    model = AeroFlex.sched.applyToROMIntegrator(model, schedUse, cfg);

    PzUse = schedUse.beam.Pz;
else
    [~, Gamma_xi, Gamma_g, ~, xi_bar_trim] = ...
        AeroFlex.core.solve_baseline_quaternion_xi( ...
            beam.fem, aero.aeroMesh, beam.phi1, beam.Nm, ...
            beam.Mglobal, rad2deg(alpha), beam.phi0);

    baseUse.Gamma_xi = Gamma_xi;
    baseUse.Gamma_g  = Gamma_g;
    if exist('xi_bar_trim','var') && ~isempty(xi_bar_trim)
        baseUse.xi_bar = xi_bar_trim;
    end

    model = AeroFlex.sim.ROMIntegrator(cfg, beam, aero, baseUse);
    PzUse = beam.Pz;
end

model.parConst.RateProject = struct('projSet', true, 'Pz', PzUse);
model.parConst.N_Thrust = zeros(numel(model.idx.q1),1);
model.parConst.gust = localZeroGust(model.parConst);
if isfield(model.parConst,'gust_input')
    model.parConst.gust_input = zeros(size(model.parConst.gust_input));
end
model = model.rebuildConstrainedOperator();
end

function [PzUse, PrUse, phi1_sA_Use] = localActiveBeamMaps(beam, schedUse, useSched, PzIn, PrIn)
PzUse = PzIn;
PrUse = PrIn;
phi1_sA_Use = beam.red.phi1_sA;

if useSched && isstruct(schedUse) && isfield(schedUse,'beam')
    if isfield(schedUse.beam,'Pz'), PzUse = schedUse.beam.Pz; end
    if isfield(schedUse.beam,'Pr'), PrUse = schedUse.beam.Pr; end
    if isfield(schedUse.beam,'red') && isfield(schedUse.beam.red,'phi1_sA')
        phi1_sA_Use = schedUse.beam.red.phi1_sA;
    end
end
end

function x0 = localPrepareInitialState(xGuess, model, baseUse, schedUse, idx, cfg, useSched, alpha)
x0 = xGuess(:);

if isempty(x0) || numel(x0) ~= size(model.L,1)
    x0 = zeros(size(model.L,1),1);
    if useSched && isstruct(schedUse) && isfield(schedUse,'x_eq') && numel(schedUse.x_eq) == numel(x0)
        x0 = schedUse.x_eq(:);
    elseif isfield(baseUse,'xi_bar') && ~isempty(baseUse.xi_bar)
        x0(idx.qxi(1)) = baseUse.xi_bar(1,1);
    end
end

if isfield(idx,'chi') && ~isempty(idx.chi)
    if useSched
        x0(idx.chi) = zeros(3,1);
    else
        aoaRef = deg2rad(localGetNested(cfg, {'flight','aoa_deg'}, 0));
        x0(idx.chi) = [0; alpha - aoaRef; 0];
    end
end
end

function [x, model, log_trim, info] = localInnerClampedSteady(x0, deltaWing, model, idx, Pz, opts)
x = x0(:);
log_trim = [];

info = struct();
info.iterations = 0;
info.residual = inf;
info.initialResidual = inf;
info.bestResidual = inf;
info.converged = false;
info.stopReason = "maxIter";

gust0 = localZeroGust(model.parConst);
if isfield(idx,'chi') && ~isempty(idx.chi)
    chiHold = x(idx.chi);
else
    chiHold = [];
end

model.parConst.N_Thrust = zeros(numel(model.idx.q1),1);
model.parConst.u_ctrl = deltaWing;
model.parConst.gust = gust0;

[res0, raw0] = evalFlexResidual(x);
info.initialResidual = res0;

bestX = x;
bestRes = res0;
bestIter = 0;

if res0 < opts.tolSS
    info.iterations = 0;
    info.residual = res0;
    info.bestResidual = res0;
    info.converged = true;
    info.stopReason = "initial";
    return
end

% During the outer optimization, an already-small projected residual is
% adequate.  This prevents thousands of inner integration steps when the
% warm start is already inside the admissible clamped manifold.
if res0 <= opts.innerAcceptTol
    info.iterations = 0;
    info.residual = res0;
    info.bestResidual = res0;
    info.converged = false;
    info.stopReason = "acceptable";
    return
end

for k = 1:opts.innerMaxIt
    [x,~] = model.step(x, deltaWing, gust0, [], false);

    if ~isempty(chiHold)
        x(idx.chi) = chiHold;
    end

    if ~all(isfinite(x))
        x = bestX;
        info.iterations = k;
        info.residual = bestRes;
        info.bestResidual = bestRes;
        info.stopReason = "nonFiniteState";
        return
    end

    if k == 1 || mod(k, opts.innerCheckEvery) == 0 || k == opts.innerMaxIt
        [resNorm, raw] = evalFlexResidual(x);

        if resNorm < bestRes
            bestRes = resNorm;
            bestX = x;
            bestIter = k;
        end

        if resNorm < opts.tolSS
            info.iterations = k;
            info.residual = resNorm;
            info.bestResidual = bestRes;
            info.converged = true;
            info.stopReason = "tol";
            return
        end

        if k - bestIter >= opts.innerPatience
            x = bestX;
            info.iterations = k;
            info.residual = bestRes;
            info.bestResidual = bestRes;
            info.stopReason = "plateau";
            return
        end

        if ~isfinite(resNorm) || norm(raw) > opts.innerDivergeRawNorm
            x = bestX;
            info.iterations = k;
            info.residual = bestRes;
            info.bestResidual = bestRes;
            info.stopReason = "diverged";
            return
        end
    end
end

x = bestX;
info.iterations = opts.innerMaxIt;
info.residual = bestRes;
info.bestResidual = bestRes;
info.stopReason = "maxIter";

    function [resNorm, raw] = evalFlexResidual(xEval)
        pc = model.parConst;
        pc.u_ctrl = deltaWing;
        pc.gust = gust0;
        pc.N_Thrust = zeros(numel(model.idx.q1),1);
        raw = AeroFlex.sim.nonlinear_terms(xEval, pc, model.idx) + model.L*xEval;
        resNorm = norm(Pz*raw(model.idx.q1));
    end
end

% =============================================================================
% Jacobian and optimizer
% =============================================================================
function [J, info] = localFiniteDifferenceJacobian(z, xSS, R0, bounds, cfg, beam, aero, base, ...
    ROMlib, useSched, idx, Pz, Pr, W, opts)

J = zeros(3,4);
info = struct();

for j = 1:2
    baseStep = opts.fdStep(j);
    gotColumn = false;

    for nTry = 1:opts.fdMaxHalves
        h = baseStep / 2^(nTry-1);
        zp = z; zm = z;
        zp(j) = min(bounds.ub(j), z(j) + h);
        zm(j) = max(bounds.lb(j), z(j) - h);

        if abs(zp(j)-zm(j)) < 1e-14
            break
        end

        [Rp, ~, ~, dp] = localResidual(zp, xSS, cfg, beam, aero, base, ROMlib, useSched, idx, Pz, Pr, W, opts);
        [Rm, ~, ~, dm] = localResidual(zm, xSS, cfg, beam, aero, base, ROMlib, useSched, idx, Pz, Pr, W, opts);

        okP = all(isfinite(Rp)) && localInnerAcceptableFD(dp, opts);
        okM = all(isfinite(Rm)) && localInnerAcceptableFD(dm, opts);

        if okP && okM
            J(:,j) = (Rp - Rm)/(zp(j)-zm(j));
            gotColumn = true;
            break
        elseif okP
            J(:,j) = (Rp - R0)/(zp(j)-z(j));
            gotColumn = true;
            break
        elseif okM
            J(:,j) = (R0 - Rm)/(z(j)-zm(j));
            gotColumn = true;
            break
        end
    end

    if ~gotColumn
        warning('TrimRBwFlex:BadFDColumn', ...
            'Could not compute a reliable finite-difference column for variable %d. Using zero column.', j);
        J(:,j) = zeros(3,1);
    end
end

% Elevator only affects the local rigid tail model used in this trim file.
hE = opts.fdStep(3);
U = localGetNested(cfg, {'flight','U_inf'}, 0);
[Ftp, Mtp] = localTailAeroForceMoment(U, z(1), z(3) + hE, cfg);
[Ftm, Mtm] = localTailAeroForceMoment(U, z(1), z(3) - hE, cfg);
dFtail = (Ftp(:) - Ftm(:))/(2*hE);
dMtail = (Mtp(:) - Mtm(:))/(2*hE);
J(:,3) = [dFtail(1); dFtail(3); dMtail(2)];

% Thrust is a rigid-body +X force, plus any moment from rThrust_B.
rT = localGetNested(cfg, {'rigidEOMset','rThrust_B'}, [0;0;0]);
dM_dT = cross(rT(:), [1;0;0]);
J(:,4) = [1; 0; dM_dT(2)];

info.rank = rank(J);
end

function dz = localLevenbergStep(J, R, z, resScale, lambda, opts)
%LOCALLEVENBERGSTEP Scaled damped Gauss-Newton step with weak trim selection.
%
% Residual equations alone give 3 equations in 4 variables.  The regularizer
% selects the low-alpha / low-control solution without overriding the force
% and moment balance.

Wr = diag(1./resScale);
Xs = diag(opts.varScale);

Js = Wr*J*Xs;
rs = Wr*R;

zErr = (z(:) - opts.zRef(:))./opts.varScale(:);
wReg = opts.regWeight(:);

H = Js.'*Js + diag(wReg) + lambda*eye(4);
g = Js.'*rs + wReg.*zErr;

if rcond(H) < 1e-14
    dy = -pinv(H)*g;
else
    dy = -H\g;
end

dz = Xs*dy;
end

function dz = localSteepestDescentStep(J, R, z, resScale, opts)
%LOCALSTEEPESTDESCENTSTEP Fallback direction in scaled variables.

Wr = diag(1./resScale);
Xs = diag(opts.varScale);

Js = Wr*J*Xs;
rs = Wr*R;

zErr = (z(:) - opts.zRef(:))./opts.varScale(:);
wReg = opts.regWeight(:);

g = Js.'*rs + wReg.*zErr;
ng = norm(g);

if ng < eps
    dz = zeros(4,1);
    return
end

dy = -g/ng;
dz = Xs*dy;
end

function [accepted, best, stepFactor] = localLineSearch(z, dz, phi0, resScale, ...
    xSS, cfg, beam, aero, base, ROMlib, useSched, idx, Pz, Pr, W, opts, bounds)

accepted = false;
stepFactor = 1.0;
best = struct('phi',inf,'z',z,'x',xSS,'R',nan(3,1),'log',[],'diag',struct());

for ls = 1:opts.lineMax
    zTry = localClamp(z + stepFactor*dz, bounds.lb, bounds.ub);

    if norm(zTry - z, inf) < 1e-14
        stepFactor = 0.5*stepFactor;
        continue
    end

    [RTry, xTry, logTry, diagTry] = localResidual(zTry, xSS, cfg, beam, aero, base, ...
        ROMlib, useSched, idx, Pz, Pr, W, opts);

    if all(isfinite(RTry)) && localInnerAcceptable(diagTry, opts)
        [phiTry, ~] = localMerit(zTry, RTry, resScale, opts);

        if phiTry < best.phi
            best = struct('phi',phiTry,'z',zTry,'x',xTry,'R',RTry,'log',logTry,'diag',diagTry);
        end

        sufficientDecrease = phiTry < phi0 - opts.armijoC*stepFactor*max(phi0, 1e-12);

        if sufficientDecrease
            accepted = true;
            return
        end
    end

    stepFactor = 0.5*stepFactor;
end
end

function [phi, rScaled] = localMerit(z, R, resScale, opts)
%LOCALMERIT Residual merit plus weak trim-selection penalty.

rScaled = R(:)./resScale(:);

zErr = (z(:) - opts.zRef(:))./opts.varScale(:);
reg = sqrt(opts.regWeight(:)).*zErr;

phi = 0.5*(rScaled.'*rScaled + reg.'*reg);
end

function [resScale, phase, forceNorm, momentScale] = localTrimScale(R, opts)
%LOCALTRIMSCALE Smoothly transition from force-dominated to full trim merit.
%
% The previous switch used forceNorm/forceSwitchNorm.  With forceSwitchNorm
% = 0.05, Mscale stayed frozen for almost the whole run.  This blend starts
% increasing My weight earlier but does not over-weight My at iteration 1.

forceNorm = norm(R(1:2)./opts.forceScale);

hi = opts.forceBlendHi;
lo = opts.forceBlendLo;

if hi <= lo
    error('TrimRBwFlex:BadBlend', ...
        'forceBlendHi must be greater than forceBlendLo.');
end

tau = (forceNorm - lo)/(hi - lo);
tau = min(1.0, max(0.0, tau));

momentScale = opts.momentScaleFull + tau*(opts.momentScaleForce - opts.momentScaleFull);

resScale = [opts.forceScale(1); opts.forceScale(2); momentScale];

if tau > 0.95
    phase = "force";
elseif tau < 0.05
    phase = "full";
else
    phase = "blend";
end
end

function ok = localInnerAcceptable(diagOut, opts)
ok = false;
if ~isstruct(diagOut) || ~isfield(diagOut,'inner') || ~isstruct(diagOut.inner)
    return
end
r = localGet(diagOut.inner,'residual',inf);
conv = localGet(diagOut.inner,'converged',false);
ok = logical(conv) || (isfinite(r) && r <= opts.innerAcceptTol);
end

function maxStepUse = localActiveMaxStep(opts, R, rNorm, phase) %#ok<INUSD>
%LOCALACTIVEMAXSTEP Residual-dependent componentwise trust region.
%
% Alpha is deliberately capped more tightly than wing/elevator.  This keeps
% the trim from using angle of attack as the cheap actuator when control
% surfaces can remove the residual.

forceNorm = norm(R(1:2)./opts.forceScale);

hi = opts.stepBlendHi;
lo = opts.stepBlendLo;

sigma = (forceNorm - lo)/(hi - lo);
sigma = min(1.0, max(0.0, sigma));

maxStepUse = opts.maxStepFine + sigma*(opts.maxStepCoarse - opts.maxStepFine);
end

function dz = localLimitStep(dz, maxStep)
for k = 1:numel(dz)
    dz(k) = min(max(dz(k), -maxStep(k)), maxStep(k));
end
end

function dz = localRespectBounds(z, dz, bounds)
dz = localClamp(z + dz, bounds.lb, bounds.ub) - z;
end

function y = localClamp(x, lb, ub)
y = min(max(x, lb), ub);
end

function ok = localInnerAcceptableFD(diagOut, opts)
%LOCALINNERACCEPTABLEFD Slightly looser gate for finite-difference columns.
%
% Finite-difference points are off the accepted trajectory.  Rejecting every
% perturbation with a residual just above innerAcceptTol makes the Jacobian
% column collapse to zero, which was visible for the wing column.

ok = false;

if ~isstruct(diagOut) || ~isfield(diagOut,'inner') || ~isstruct(diagOut.inner)
    return
end

r = localGet(diagOut.inner,'residual',inf);
conv = localGet(diagOut.inner,'converged',false);

ok = logical(conv) || (isfinite(r) && r <= opts.fdInnerAcceptTol);
end

function [didPolish,zNew,RNew,dT] = localTryThrustPolish(z, R, bounds, cfg, opts)
didPolish = false;
zNew = z;
RNew = R;
dT = 0;

if ~opts.thrustPolish || localGetNested(cfg, {'trim','thrustActsOnWing'}, false)
    return
end

% ready = abs(R(2)) <= opts.thrustPolishTol(1) && ...
%         abs(R(3)) <= opts.thrustPolishTol(2);
% if ~ready
%     return
% end

% Do not let thrust cleanup corrupt pitch trim if the thrust line has a
% vertical offset.  Fz is not part of this gate because thrust does not act in
% body z in this trim model.
myTol = localGet(opts, 'thrustPolishMyTol', opts.thrustPolishTol(2));
if abs(R(3)) > myTol
    return
end


Told = z(4);
Tnew = min(max(Told - R(1), bounds.lb(4)), bounds.ub(4));
dT = Tnew - Told;

if abs(dT) < 1e-12
    return
end

rT = localGetNested(cfg, {'rigidEOMset','rThrust_B'}, [0;0;0]);
dM_dT = cross(rT(:), [1;0;0]);

Rtry = R;
Rtry(1) = R(1) + dT;
Rtry(3) = R(3) + dM_dT(2)*dT;

if abs(Rtry(1)) < abs(R(1)) && abs(Rtry(3)) <= max(abs(R(3)), opts.tolPhys(3))
    didPolish = true;
    zNew = z;
    zNew(4) = Tnew;
    RNew = Rtry;
else
    dT = 0;
end
end

function [didPolish,zNew,RNew,xNew,logNew,diagNew,info] = ...
    localTryEndgameBlockPolish(z, R, xSS, log_trim, diag0, bounds, ...
    cfg, beam, aero, base, ROMlib, useSched, idx, Pz, Pr, W, opts)
%LOCALTRYENDGAMEBLOCKPOLISH Small block Newton step for [Fz; My].
%
% The full trim problem has three residuals and four variables.  Near
% convergence, Fx is best handled by thrust and [Fz; My] by wing/elevator.
% This routine performs only a small accepted-iterate correction; it is not
% used inside localResidual and does not alter the main trim path.

didPolish = false;
zNew = z;
RNew = R;
xNew = xSS;
logNew = log_trim;
diagNew = diag0;

info = struct();
info.reason = "not_attempted";
info.dz = zeros(4,1);

if ~localInnerAcceptable(diag0, opts)
    info.reason = "bad_inner";
    return
end

if norm(R([2 3])./[opts.tolPhys(2); opts.tolPhys(3)]) <= 1.0
    info.reason = "already_within_tolerance";
    return
end

[J, ~] = localFiniteDifferenceJacobian(z, xSS, R, bounds, cfg, beam, aero, base, ...
    ROMlib, useSched, idx, Pz, Pr, W, opts);

cols = opts.endgamePolishVars(:).';
rows = [2 3];

A = J(rows, cols);
b = -R(rows);

if any(~isfinite(A(:))) || any(~isfinite(b(:))) || rank(A) < numel(cols)
    info.reason = "singular_block";
    return
end

if rcond(A) < opts.endgameBlockRcond
    info.reason = "ill_conditioned_block";
    return
end

dz = zeros(4,1);
dz(cols) = A\b;

% Keep this as a polish, not a new outer step.
dz = localLimitStep(dz, opts.endgamePolishMaxStep);
dz = localRespectBounds(z, dz, bounds);

if norm(dz,inf) < 1e-14
    info.reason = "zero_step";
    return
end

e0 = norm(R(rows)./[opts.tolPhys(2); opts.tolPhys(3)]);

for ls = 1:opts.endgamePolishLineMax
    sf = 0.5^(ls-1);
    zTry = localClamp(z + sf*dz, bounds.lb, bounds.ub);

    [RTry, xTry, logTry, diagTry] = localResidual(zTry, xSS, cfg, beam, aero, base, ...
        ROMlib, useSched, idx, Pz, Pr, W, opts);

    if ~all(isfinite(RTry)) || ~localInnerAcceptable(diagTry, opts)
        continue
    end

    eTry = norm(RTry(rows)./[opts.tolPhys(2); opts.tolPhys(3)]);

    % Fx may change slightly through wing/elevator; that is acceptable
    % because thrust polish handles Fx exactly after this block.
    fxOK = abs(RTry(1)) <= abs(R(1)) + opts.endgameFxRelax;

    if eTry < opts.endgameAcceptFrac*e0 && fxOK
        didPolish = true;
        zNew = zTry;
        RNew = RTry;
        xNew = xTry;
        logNew = logTry;
        diagNew = diagTry;
        info.reason = "accepted";
        info.dz = sf*dz;
        return
    end
end

info.reason = "no_accepted_line_step";
info.dz = dz;
end

% =============================================================================
% Initialization and options
% =============================================================================
function history = localInitHistory(maxIt)
history = struct();
history.z = nan(4,maxIt);
history.R = nan(3,maxIt);
history.RscaledNorm = nan(1,maxIt);
history.cost = nan(1,maxIt);
history.lambda = nan(1,maxIt);
history.lineStep = nan(1,maxIt);
history.condJ = nan(1,maxIt);
history.phase = strings(1,maxIt);
history.innerResidual = nan(1,maxIt);
history.innerStop = strings(1,maxIt);
history.thrustPolish = false(1,maxIt);
end

function z = localInitialTrimVariables(cfg, bounds)
alpha0Deg = localGetNested(cfg, {'flight','aoa_deg'}, 0);
if isfield(cfg,'trim') && isfield(cfg.trim,'alphaGuessDeg') && ~isempty(cfg.trim.alphaGuessDeg)
    alpha0Deg = cfg.trim.alphaGuessDeg;
end

deltaWing0 = 0;
if isfield(cfg,'trim') && isfield(cfg.trim,'deltaWing') && ~isempty(cfg.trim.deltaWing)
    deltaWing0 = cfg.trim.deltaWing(1);
elseif isfield(cfg,'trim') && isfield(cfg.trim,'deltaDeg') && ~isempty(cfg.trim.deltaDeg)
    d0 = cfg.trim.deltaDeg;
    if max(abs(d0(:))) > 2*pi
        deltaWing0 = deg2rad(d0(1));
    else
        deltaWing0 = d0(1);
    end
end

deltaElev0 = 0;
if isfield(cfg,'trim') && isfield(cfg.trim,'deltaElev') && ~isempty(cfg.trim.deltaElev)
    deltaElev0 = cfg.trim.deltaElev(1);
elseif isfield(cfg,'trim') && isfield(cfg.trim,'deltaElevDeg') && ~isempty(cfg.trim.deltaElevDeg)
    deltaElev0 = deg2rad(cfg.trim.deltaElevDeg(1));
end

T0 = 0;
if isfield(cfg,'trim') && isfield(cfg.trim,'thrustScalar') && ~isempty(cfg.trim.thrustScalar)
    T0 = cfg.trim.thrustScalar;
elseif isfield(cfg,'trim') && isfield(cfg.trim,'thrust') && isscalar(cfg.trim.thrust)
    T0 = cfg.trim.thrust;
end

z = localClamp([deg2rad(alpha0Deg); deltaWing0; deltaElev0; T0], bounds.lb, bounds.ub);
end

function opts = localTrimOptions(cfg, W, bounds, z0) %#ok<INUSD>
Fref = max(W, 1.0);
Lref = localGetNested(cfg, {'struct','L'}, 1.1);
cRef = localGetNested(cfg, {'geom','c_ref'}, ...
       localGetNested(cfg, {'struct','c_ref'}, ...
       localGetNested(cfg, {'flight','c_ref'}, 0.10)));

momentArmScale = localGetNested(cfg, {'trim','momentArmScale'}, cRef);

opts = struct();
opts.tolSS      = localGetNested(cfg, {'trim','tolSteady'}, 1e-6);
opts.tolNewton  = localGetNested(cfg, {'trim','tolNewtonScaled'}, ...
                  localGetNested(cfg, {'trim','tolNewton'}, 1e-4));
opts.tolPhys    = localGetNested(cfg, {'trim','tolPhysical'}, [1e-4; 1e-4; 1e-4]);

opts.maxIt      = localGetNested(cfg, {'trim','maxNewtonIter'}, 80);
opts.innerMaxIt = localGetNested(cfg, {'trim','itersInt'}, 20000);
opts.debug      = localGetNested(cfg, {'debug','trim'}, true);

opts.forceScale = localGetNested(cfg, {'trim','forceScale'}, [Fref; Fref]);

opts.momentScaleForce = localGetNested(cfg, {'trim','momentScaleForce'}, ...
    max(W*max(Lref,0.1), 1.0));

opts.momentScaleFull = localGetNested(cfg, {'trim','momentScaleFull'}, ...
    max(W*momentArmScale, 0.05));

% Smooth merit transition.  My begins to matter before the force residual is
% tiny, but it is not over-weighted at the first iteration.
opts.forceBlendHi = localGetNested(cfg, {'trim','forceBlendHi'}, 0.30);
opts.forceBlendLo = localGetNested(cfg, {'trim','forceBlendLo'}, 0.03);

% Residual-dependent trust-region blend.
opts.stepBlendHi = localGetNested(cfg, {'trim','stepBlendHi'}, 0.30);
opts.stepBlendLo = localGetNested(cfg, {'trim','stepBlendLo'}, 0.03);

% Variable metric.  Alpha has a smaller scale than control surfaces, so a
% unit scaled alpha step is physically smaller than a unit scaled delta step.
opts.varScale = localGetNested(cfg, {'trim','varScale'}, ...
    [deg2rad(0.25); deg2rad(2.00); deg2rad(3.00); max(0.20,0.05*W)]);

opts.maxStepCoarse = localGetNested(cfg, {'trim','maxStepCoarse'}, ...
    [deg2rad(0.08); deg2rad(0.35); deg2rad(0.75); max(0.08,0.02*W)]);

opts.maxStepFine = localGetNested(cfg, {'trim','maxStepFine'}, ...
    [deg2rad(0.015); deg2rad(0.08); deg2rad(0.20); max(0.02,0.005*W)]);

% Retain legacy fields as aliases for older config files.
opts.maxStep = localGetNested(cfg, {'trim','maxStep'}, opts.maxStepCoarse);
opts.fineStepNorm = localGetNested(cfg, {'trim','fineStepNorm'}, 5e-2);

opts.lambda0   = localGetNested(cfg, {'trim','lambda0'}, 1e-2);
opts.lambdaMin = localGetNested(cfg, {'trim','lambdaMin'}, 1e-4);
opts.lambdaMax = localGetNested(cfg, {'trim','lambdaMax'}, 1e6);
opts.lambdaDown = localGetNested(cfg, {'trim','lambdaDown'}, 2);
opts.lambdaUp   = localGetNested(cfg, {'trim','lambdaUp'}, 5);
opts.lambdaUpSoft = localGetNested(cfg, {'trim','lambdaUpSoft'}, 2);

opts.lineMax = localGetNested(cfg, {'trim','lineSearchMax'}, 10);
opts.armijoC = localGetNested(cfg, {'trim','armijoC'}, 1e-4);
opts.maxStall = localGetNested(cfg, {'trim','maxStall'}, 10);

opts.fdStep = localGetNested(cfg, {'trim','fdStep'}, ...
    [deg2rad(0.04); deg2rad(0.04); deg2rad(0.03); max(0.01,0.002*W)]);

opts.fdMaxHalves = localGetNested(cfg, {'trim','fdMaxHalves'}, 3);

opts.innerCheckEvery = localGetNested(cfg, {'trim','innerCheckEvery'}, 10);
opts.innerPatience = localGetNested(cfg, {'trim','innerPatience'}, 2500);
opts.innerDivergeRawNorm = localGetNested(cfg, {'trim','innerDivergeRawNorm'}, 1e12);

opts.innerAcceptTol = localGetNested(cfg, {'trim','innerAcceptTol'}, ...
    max(5*opts.tolSS, 5e-6));

opts.fdInnerAcceptTol = localGetNested(cfg, {'trim','fdInnerAcceptTol'}, ...
    max(20*opts.tolSS, 2e-5));

opts.thrustPolish = localGetNested(cfg, {'trim','thrustPolish'}, true);
opts.thrustPolishTol = localGetNested(cfg, {'trim','thrustPolishTol'}, [1e-4; 1e-4]);

% Axial thrust cleanup only needs pitch moment to be small.  Fz is not
% affected by a body-X thrust in the current model.
opts.thrustPolishMyTol = localGetNested(cfg, {'trim','thrustPolishMyTol'}, ...
    max(opts.tolPhys(3), 2e-3));

% Endgame block polish for [Fz; My] using wing/elevator.
opts.endgamePolish = localGetNested(cfg, {'trim','endgamePolish'}, true);
opts.endgamePolishNorm = localGetNested(cfg, {'trim','endgamePolishNorm'}, 1.0e-2);
opts.endgamePolishVars = localGetNested(cfg, {'trim','endgamePolishVars'}, [2 3]);
opts.endgamePolishMaxStep = localGetNested(cfg, {'trim','endgamePolishMaxStep'}, ...
    [0; deg2rad(0.03); deg2rad(0.08); 0]);
opts.endgamePolishLineMax = localGetNested(cfg, {'trim','endgamePolishLineMax'}, 5);
opts.endgameAcceptFrac = localGetNested(cfg, {'trim','endgameAcceptFrac'}, 0.75);
opts.endgameBlockRcond = localGetNested(cfg, {'trim','endgameBlockRcond'}, 1e-10);
opts.endgameFxRelax = localGetNested(cfg, {'trim','endgameFxRelax'}, 5e-2);

% Weak Tikhonov term for non-unique trim.  Alpha is penalized most strongly;
% wing/elevator regularization is only for nullspace selection.
opts.regWeight = localGetNested(cfg, {'trim','regWeight'}, ...
    [1e-3; 1e-6; 1e-7; 0]);

alphaRef = deg2rad(localGetNested(cfg, {'flight','aoa_deg'}, 0));
opts.zRef = localGetNested(cfg, {'trim','zRef'}, ...
    [alphaRef; 0; 0; z0(4)]);

opts.debugEvery = localGetNested(cfg, {'trim','debugEvery'}, 1);
opts.fullDebugNorm = localGetNested(cfg, {'trim','fullDebugNorm'}, 5e-2);
end

function bounds = localTrimBounds(cfg, ROMlib, useSched, W, fixedExactOwner)
if useSched && isstruct(ROMlib) && isfield(ROMlib,'mu') && size(ROMlib.mu,2) >= 2
    aMin = min(ROMlib.mu(:,2));
    aMax = max(ROMlib.mu(:,2));
else
    aMin = -10;
    aMax = 15;
end

alphaBoundsDeg = localGetNested(cfg, {'trim','alphaBoundsDeg'}, [aMin, aMax]);
if fixedExactOwner.enabled
    alphaBoundsDeg = fixedExactOwner.physicalAlphaBoundsDeg;
elseif useSched
    alphaBoundsDeg(1) = max(alphaBoundsDeg(1), aMin);
    alphaBoundsDeg(2) = min(alphaBoundsDeg(2), aMax);
end

wingBoundsDeg = localGetNested(cfg, {'trim','wingDeltaBoundsDeg'}, [-20,20]);
elevBoundsDeg = localGetNested(cfg, {'trim','elevatorBoundsDeg'}, [-25,25]);
thrustBounds  = localGetNested(cfg, {'trim','thrustBoundsN'}, [0, max(5,2*W)]);

bounds.lb = [deg2rad(alphaBoundsDeg(1)); deg2rad(wingBoundsDeg(1)); deg2rad(elevBoundsDeg(1)); thrustBounds(1)];
bounds.ub = [deg2rad(alphaBoundsDeg(2)); deg2rad(wingBoundsDeg(2)); deg2rad(elevBoundsDeg(2)); thrustBounds(2)];

assert(all(bounds.lb < bounds.ub), 'TrimRBwFlex:BadBounds', 'Invalid trim bounds.');
end

function owner = localFixedExactPackageTrimOwner(cfg, ROMlib, useSched)
owner = struct('enabled',false,'physicalAlphaBoundsDeg',[NaN,NaN]);
if ~isfield(cfg,'trim') || ~isstruct(cfg.trim) || ...
        ~isfield(cfg.trim,'fixedExactPackageOwner')
    return
end

request = cfg.trim.fixedExactPackageOwner;
assert(isstruct(request) && isscalar(request), ...
    'TrimRBwFlex:FixedExactOwnerType', ...
    'fixedExactPackageOwner must be a scalar structure.');
if ~isfield(request,'enabled') || ~logical(request.enabled)
    return
end

expectedChange = "phase18c-v17a-fixed-alpha-exact-package-coupled-trim-owner-v1";
expectedSource = "CASEA_U040_ALPHA_P01";
expectedHash = "7c501bb9782332395b7be7d9f14aa310c17398b7640f565fcf51ce256f9dec02";
expectedQuery = [40,1];
required = {'auditOnly','changeControlId','sourceId','packageSha256', ...
    'query','physicalAlphaBoundsDeg'};
assert(all(isfield(request,required)), ...
    'TrimRBwFlex:FixedExactOwnerFields', ...
    'The fixed exact-package trim request is incomplete.');
assert(logical(request.auditOnly) && ...
    string(request.changeControlId)==expectedChange && ...
    string(request.sourceId)==expectedSource && ...
    lower(string(request.packageSha256))==expectedHash, ...
    'TrimRBwFlex:FixedExactOwnerProvenance', ...
    'The fixed exact-package trim request is not bound to the approved source.');

query = double(request.query(:).');
alphaBounds = double(request.physicalAlphaBoundsDeg(:).');
assert(useSched && isstruct(ROMlib) && size(ROMlib.mu,1)==1 && ...
    numel(query)==2 && all(isfinite(query)) && ...
    norm(query-expectedQuery,inf)<=1e-12 && ...
    norm(double(ROMlib.mu(1,:))-query,inf)<=1e-12, ...
    'TrimRBwFlex:FixedExactOwnerQuery', ...
    'The audit selector requires the singleton exact package at [40,1].');
assert(numel(alphaBounds)==2 && all(isfinite(alphaBounds)) && ...
    alphaBounds(1)<alphaBounds(2), ...
    'TrimRBwFlex:FixedExactOwnerBounds', ...
    'Physical alpha bounds must be finite and strictly increasing.');
assert(isfield(cfg,'library') && isfield(cfg.library,'freezeInOptimizer') && ...
    logical(cfg.library.freezeInOptimizer) && ...
    isfield(cfg.library,'schedulerQueryPoint') && ...
    norm(double(cfg.library.schedulerQueryPoint(:).')-query,inf)<=1e-12, ...
    'TrimRBwFlex:FixedExactOwnerFreeze', ...
    'The exact package query must remain frozen throughout trim.');

owner.enabled = true;
owner.physicalAlphaBoundsDeg = alphaBounds;
end

function [x0, sched0] = localInitialState(cfg, beam, aero, base, idx, ROMlib, useSched, alpha)
sched0 = [];
Nx = 3*beam.Nm + 1 + aero.Na + 3;

x0 = zeros(Nx,1);
if isfield(cfg,'trim') && isfield(cfg.trim,'stateGuess') && ...
        numel(cfg.trim.stateGuess) == Nx && all(isfinite(cfg.trim.stateGuess(:)))
    x0 = cfg.trim.stateGuess(:);
    if isfield(idx,'chi'), x0(idx.chi) = zeros(3,1); end
elseif useSched
    mu = localTrimSchedulePoint(cfg, alpha);
    sched0 = localEvalSchedule(ROMlib,mu,cfg.library);
    localRequireExactSchedule(cfg, sched0, mu);
    if isfield(sched0,'x_eq') && numel(sched0.x_eq) == Nx
        x0 = sched0.x_eq(:);
    elseif isfield(sched0,'trim') && isfield(sched0.trim,'states') && numel(sched0.trim.states) == Nx
        x0 = sched0.trim.states(:);
    elseif isfield(base,'xi_bar') && ~isempty(base.xi_bar)
        x0(idx.qxi(1)) = base.xi_bar(1,1);
    end
    if isfield(idx,'chi'), x0(idx.chi) = zeros(3,1); end
else
    if isfield(base,'xi_bar') && ~isempty(base.xi_bar)
        x0(idx.qxi(1)) = base.xi_bar(1,1);
    end
    if isfield(idx,'chi')
        aoaRef = deg2rad(localGetNested(cfg, {'flight','aoa_deg'}, 0));
        x0(idx.chi) = [0; alpha - aoaRef; 0];
    end
end
end

function sched = localEvalSchedule(ROMlib,mu,libraryCfg)
if isfield(libraryCfg,'frozenPackage') && ...
        isstruct(libraryCfg.frozenPackage) && ...
        ~isempty(fieldnames(libraryCfg.frozenPackage))
    sched = libraryCfg.frozenPackage;
    if ~isfield(sched,'mu') || norm(sched.mu(:)-mu(:),inf) > 1e-12
        error('TrimRBwFlex:FrozenPackageQueryMismatch', ...
            'Frozen package query differs from the trim scheduler query.');
    end
else
    sched = AeroFlex.sched.evalLibrary(ROMlib,mu,libraryCfg);
end
end

function mu = localTrimSchedulePoint(cfg, alpha)
freezePackage = isfield(cfg,'library') && ...
    isfield(cfg.library,'freezeInOptimizer') && ...
    logical(cfg.library.freezeInOptimizer);
if freezePackage && isfield(cfg.library,'schedulerQueryPoint') && ...
        numel(cfg.library.schedulerQueryPoint) == 2
    mu = double(cfg.library.schedulerQueryPoint(:).');
else
    mu = [cfg.flight.U_inf, rad2deg(alpha)];
end
end

function localRequireExactSchedule(cfg, sched, mu)
requireExact = isfield(cfg,'library') && ...
    isfield(cfg.library,'requireExactNode') && logical(cfg.library.requireExactNode);
if ~requireExact
    return
end
tol = 1e-10;
if isfield(cfg.library,'interpTol'), tol = cfg.library.interpTol; end
isExact = numel(sched.pointIds) == 1 && numel(sched.weights) == 1 && ...
    abs(sched.weights(1) - 1) <= tol && ...
    norm(sched.pointMu(1,:) - mu, inf) <= tol;
if ~isExact
    error('TrimRBwFlex:ExactNodeRequired', ...
        'Trim requires one exact source node at U=%.12g m/s, alpha=%.12g deg.', ...
        mu(1), mu(2));
end
end

% =============================================================================
% Output and diagnostics
% =============================================================================
function trim = localBuildOutput(z, xFinal, log_trim, Rfinal, rFinalNorm, ...
    converged, bounds, history, useSched, opts, diagFinal, resScaleFinal)

trim = struct();
trim.converged = converged;
trim.alphaDeg = rad2deg(z(1));
trim.alphaRad = z(1);

% Legacy aliases used by some setup/runtime print paths.
trim.aoaDeg = trim.alphaDeg;
trim.aoaRad = trim.alphaRad;
trim.alpha  = trim.alphaRad;

trim.deltaWing = [z(2); z(2); 0; 0];
trim.deltaWingDeg = rad2deg(trim.deltaWing);
trim.deltaDeg = trim.deltaWingDeg;

trim.deltaElev = z(3);
trim.deltaElevDeg = rad2deg(z(3));
trim.thrust = z(4);
trim.thrustScalar = z(4);

trim.states = xFinal(:);
trim.sensor = log_trim;
trim.residual = Rfinal(:);
trim.residualScaled = Rfinal(:)./resScaleFinal(:);
trim.residualScaledNorm = rFinalNorm;
trim.history = history;
trim.bounds = bounds;
trim.useSched = useSched;
trim.opts = opts;

if useSched && isfield(diagFinal,'sched')
    trim.sched = diagFinal.sched;
end

idx = diagFinal.FlexPlant.idx;
if isfield(idx,'chi')
    trim.chi0 = trim.states(idx.chi);
else
    trim.chi0 = [0;0;0];
end

trim.uWingTrim = zeros(4,1);
trim.uWingTrim(1:2) = trim.deltaWing(1:2);
trim.uWingTrim(3:4) = 0;
trim.u_ctrl = trim.deltaWing;

trim.debug = struct();
trim.debug.loads = struct();
trim.debug.loads.Clamp6    = diagFinal.Clamp6(:);
trim.debug.loads.Fwing_B   = diagFinal.Fwing_B(:);
trim.debug.loads.Mwing_B   = diagFinal.Mwing_B(:);
trim.debug.loads.Ftail_B   = diagFinal.Ftail_B(:);
trim.debug.loads.Mtail_B   = diagFinal.Mtail_B(:);
trim.debug.loads.Ffin_B    = diagFinal.Ffin_B(:);
trim.debug.loads.Mfin_B    = diagFinal.Mfin_B(:);
trim.debug.loads.Fgrav_B   = diagFinal.Fgrav_B(:);
trim.debug.loads.Fthrust_B = diagFinal.Fthrust_B(:);
trim.debug.loads.Mthrust_B = diagFinal.Mthrust_B(:);
trim.debug.loads.Ftot_B    = diagFinal.Ftot_B(:);
trim.debug.loads.Mtot_B    = diagFinal.Mtot_B(:);
trim.debug.loads.R         = Rfinal(:);
trim.debug.inner           = diagFinal.inner;
trim.debug.Pz_q1dot_norm   = norm(diagFinal.Pz*diagFinal.raw_q1dot(:));
trim.debug.Pr_q1dot_norm   = norm(diagFinal.Pr*diagFinal.raw_q1dot(:));
end

function localPrintFinalDiagnostics(trim, diagFinal, opts)
idx = diagFinal.FlexPlant.idx;
pc = diagFinal.FlexPlant.parConst;
pc.u_ctrl = trim.u_ctrl;
pc.gust = localZeroGust(pc);
pc.N_Thrust = zeros(numel(idx.q1),1);

raw = AeroFlex.sim.nonlinear_terms(trim.states, pc, idx) + ...
      diagFinal.FlexPlant.L*trim.states;

projectedQ1 = norm(diagFinal.Pz*raw(idx.q1));
reactionQ1  = norm(diagFinal.Pr*raw(idx.q1));

[xNextProjected,~] = diagFinal.FlexPlant.step(trim.states, trim.u_ctrl, pc.gust, [], false);
dxProjected = xNextProjected - trim.states;

fprintf('\n[TrimRBwFlex verification]\n');
fprintf('  continuous projected |Pz*q1dot| = %.6e\n', projectedQ1);
fprintf('  continuous reaction  |Pr*q1dot| = %.6e\n', reactionQ1);
fprintf('  projected one-step   |dx|/dt    = %.6e\n', norm(dxProjected)/diagFinal.FlexPlant.dt);
fprintf('  loads residual       Fx=%+.3e N, Fz=%+.3e N, My=%+.3e N*m\n', ...
    trim.residual(1), trim.residual(2), trim.residual(3));

if opts.debug
    fprintf('[TrimRBwFlex] Final: converged=%d | alpha=%.5f deg | wing=%.5f deg | elev=%.5f deg | T=%.5f N | |R_s|=%.3e\n\n', ...
        trim.converged, trim.alphaDeg, trim.deltaWingDeg(1), ...
        trim.deltaElevDeg, trim.thrust, trim.residualScaledNorm);
end
end

% =============================================================================
% Force models and utilities
% =============================================================================
function Clamp6 = localApplyWingSymmetry(Clamp6, cfg)
if localGetNested(cfg, {'trim','mirrorWingClamp'}, true)
    S = diag([1,-1,1,-1,1,-1]);
    Clamp6 = Clamp6 + S*Clamp6;
elseif isfield(cfg,'sim') && isfield(cfg.sim,'wingMultiplicity')
    Clamp6 = cfg.sim.wingMultiplicity*Clamp6;
end
end

function Fgrav_B = localGravityBody(W, alpha)
Fgrav_B = [-W*sin(alpha); 0; W*cos(alpha)];
end

function [F_tail, M_tail] = localTailAeroForceMoment(U, alpha, delta_e, cfg)
rho = localGetNested(cfg, {'flight','rho'}, 1.225);
qbar = 0.5*rho*U^2;

S_t  = localGetNested(cfg, {'tail','S'}, localGetNested(cfg, {'tail','area'}, 0.035));
CL_a = localGetNested(cfg, {'tail','CL_alpha'}, 4.4);
eta  = localGetNested(cfg, {'tail','eta_delta'}, 1.0);
i_t  = deg2rad(localGetNested(cfg, {'tail','incidence_deg'}, 0.01));
CD0  = localGetNested(cfg, {'tail','CD0'}, 0.01);
e    = localGetNested(cfg, {'tail','e'}, 0.8);
AR   = localGetNested(cfg, {'tail','AR'}, 4.0);
c_t  = localGetNested(cfg, {'tail','c'}, localGetNested(cfg, {'tail','chord'}, 0.10));
Cm0  = localGetNested(cfg, {'tail','Cm0'}, 0.0);
Cm_delta = localGetNested(cfg, {'tail','Cm_delta'}, 0.0);
l_t  = localGetNested(cfg, {'tail','arm'}, 0.60);

if isfield(cfg,'tail') && isfield(cfg.tail,'r_B') && ~isempty(cfg.tail.r_B)
    rTail_B = cfg.tail.r_B(:);
elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'rTail_B') && ~isempty(cfg.rigidEOMset.rTail_B)
    rTail_B = cfg.rigidEOMset.rTail_B(:);
elseif isfield(cfg,'geom') && isfield(cfg.geom,'rTail_B') && ~isempty(cfg.geom.rTail_B)
    rTail_B = cfg.geom.rTail_B(:);
else
    rTail_B = [-l_t; 0; 0];
end

CLt = CL_a*(alpha - i_t - eta*delta_e);
CDt = CD0 + CLt^2/(pi*e*AR);

Lt = qbar*S_t*CLt;
Dt = qbar*S_t*CDt;

F_tail = [-Dt; 0; -Lt];
M_ac_B = [0; qbar*S_t*c_t*(Cm0 + Cm_delta*delta_e); 0];
M_tail = M_ac_B + cross(rTail_B, F_tail);
end

function [F_fin, M_fin] = localFinAeroForceMoment(~, ~, ~, ~)
F_fin = zeros(3,1);
M_fin = zeros(3,1);
end

function g = localZeroGust(pc)
if isfield(pc,'Bw') && ~isempty(pc.Bw)
    g = zeros(size(pc.Bw,2),1);
else
    g = 0;
end
end

function cfg = localFillTrimDefaults(cfg)
if ~isfield(cfg,'trim') || ~isstruct(cfg.trim), cfg.trim = struct(); end
if ~isfield(cfg,'debug') || ~isstruct(cfg.debug), cfg.debug = struct(); end
if ~isfield(cfg.debug,'trim'), cfg.debug.trim = true; end
if ~isfield(cfg.trim,'useRateProjection'), cfg.trim.useRateProjection = true; end
if ~isfield(cfg.trim,'thrustActsOnWing'), cfg.trim.thrustActsOnWing = false; end
if isfield(cfg,'library') && isstruct(cfg.library) && ~isfield(cfg.library,'requireCompatible')
    cfg.library.requireCompatible = true;
end
end

function val = localGet(S, field, defaultVal)
if isstruct(S) && isfield(S,field) && ~isempty(S.(field))
    val = S.(field);
else
    val = defaultVal;
end
end

function val = localGetNested(S, fields, defaultVal)
val = defaultVal;
try
    tmp = S;
    for i = 1:numel(fields)
        if ~isstruct(tmp) || ~isfield(tmp, fields{i}) || isempty(tmp.(fields{i}))
            return
        end
        tmp = tmp.(fields{i});
    end
    val = tmp;
catch
    val = defaultVal;
end
end
