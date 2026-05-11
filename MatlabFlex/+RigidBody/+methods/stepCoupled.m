function xNext = stepCoupled(obj, u_cmd, gust)
%STEPCOUPLED One coupled rigid-flexible time step.

    % -------------------------------------------------------------
    % 1. Extract current rigid-body state
    % -------------------------------------------------------------
    rb = obj.rbState;

    % Compute local flight condition at wing root.
    Uinf  = rb.Uinf;
    alpha = rb.alpha;
    beta  = rb.beta;
    chi   = rb.euler;
    omega = rb.bodyRates;

    % -------------------------------------------------------------
    % 2. Update wing ROM parameters from rigid-body motion
    % -------------------------------------------------------------
    obj.updateWingFlightCondition(Uinf, alpha, beta, chi, omega);

    % -------------------------------------------------------------
    % 3. Propagate flexible wing ROM
    % -------------------------------------------------------------
    xFlexNext = obj.stepWingROM(u_cmd, gust);

    % -------------------------------------------------------------
    % 4. Compute wing reaction/resultant at the fuselage/root
    % -------------------------------------------------------------
    Clamp6 = obj.computeWingClampReaction();
    Fwing_B = Clamp6(1:3);
    Mwing_B = Clamp6(4:6);

    % -------------------------------------------------------------
    % 5. Compute rigid-body-only contributions
    % -------------------------------------------------------------
    [Ftail_B, Mtail_B] = obj.computeTailLoads(u_cmd);
    [Ffin_B,  Mfin_B]  = obj.computeFinLoads(u_cmd);
    [Fth_B,   Mth_B]   = obj.computeThrustLoads();
    Fgrav_B            = obj.computeGravityBody();

    % -------------------------------------------------------------
    % 6. Sum total aircraft forces/moments
    % -------------------------------------------------------------
    Ftot_B = Fwing_B + Ftail_B + Ffin_B + Fth_B + Fgrav_B;
    Mtot_B = Mwing_B + Mtail_B + Mfin_B + Mth_B;

    % -------------------------------------------------------------
    % 7. Propagate rigid-body equations
    % -------------------------------------------------------------
    rbNext = obj.stepRigidBody(Ftot_B, Mtot_B);

    % -------------------------------------------------------------
    % 8. Commit states
    % -------------------------------------------------------------
    obj.xFlex  = xFlexNext;
    obj.rbState = rbNext;

    xNext = obj.packCoupledState();
end