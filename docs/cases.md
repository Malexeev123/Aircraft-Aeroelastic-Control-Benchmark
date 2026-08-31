# Benchmark cases and model workflows

## Formal cases

| ID | Command and disturbance | Status |
| --- | --- | --- |
| A1 | commanded attitude, no gust | qualified |
| A2 | attitude hold, frozen gust | qualified |
| A3 | commanded attitude, frozen gust | qualified |
| B1 | thrust-led scheduled speed transition, no gust | experimental scheduled case |
| B2 | matched scheduled speed transition, frozen gust | experimental scheduled case |
| C | longitudinal trajectory and attitude tracking under gust | deferred; definition not frozen |

The Case-A execution owner is unchanged. The Case-B plan binds the executed
V90 experimental profile. Its `auditDisposition` field is the machine-readable
record of retained speedups and intentionally disabled experiments.

V90 retains compiled estimator/controller kernels, active-stencil
sensitivities, condensed RTI, packet reuse, native reduced horizons, the P1
plant interval, rigid-wrench refinement, safe endpoint realization, and the
corrected reciprocal `qGam` bound owner. It deliberately keeps
prepared-horizon reuse, accepted-replay reuse, full online fmincon correction,
future scheduled-package sequencing, and the later standalone condensation
experiments disabled. These selections are part of the retained executable
profile; the facade does not infer or reactivate discarded speed experiments.

## General model workflows

The lower-level setup and execution path remains available for:

- wing-only open-loop propagation;
- wing-only nMHE/nMPC gust-load-alleviation control;
- coupled-full open-loop no-gust trim/replay verification;
- coupled-full controlled development runs.

These paths are prepared with `sim_init` and executed with `sim_run`. The
wing-only configuration keeps rate projection off; the coupled configuration
uses its qualified coupled policy. The general workflows are not aliases for
formal free-flight Case A.

## Custom and combined maneuvers

`AeroFlex.benchmark.customCaseDefinition` describes custom speed, attitude,
gust, and initial-condition histories without modifying a formal case. A
combined attitude-and-speed maneuver uses one scheduled runtime across
the complete timeline; it does not splice the exact-source Case-A runner into
the Case-B runner or reset estimator/controller memory at an arbitrary label.

The formal Case-A and Case-B initial conditions are not interchangeable. If a
requested phase boundary cannot be reached continuously inside the scheduled
domain, it must be represented as two separately manifested runs or supplied
with a separately qualified transition. The scheduled coupled-full speed/pitch
adapter executes only with `AllowUnqualified=true`; altitude/lateral guidance
and nonzero initial perturbations remain fail-closed.

## Qualification rule

A trajectory that appears useful or stable is not automatically a benchmark
pass. Qualification also requires the declared estimator/controller success,
constraint, actuator, thrust, loads, source-domain, and fallback gates. An
experimental Case-B run preserves and plots failed evidence instead of
discarding it or relabeling it.
