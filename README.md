# Pazy Aeroelastic Control Benchmark

This repository implements a coupled nonlinear aeroelastic estimation and
control benchmark for a Pazy-wing aircraft. The MATLAB model combines
intrinsic/modal structural dynamics, unsteady aerodynamic states, rigid-body
motion, trim, ROM scheduling, outer-loop control, nonlinear moving-horizon
estimation, nonlinear model predictive control, command fusion, and actuator
realization.

SHARPy and XBeam are external source-generation dependencies. The benchmark
uses their published interfaces unchanged; conversion, scheduling, control,
validation, and post-processing are project-owned.

## Current benchmark status

| Scenario | Description | Status |
| --- | --- | --- |
| `A1` | Commanded attitude, no gust | Qualified |
| `A2` | Attitude hold under gust | Qualified |
| `A3` | Commanded attitude under gust | Qualified |
| `B1` | Thrust-led speed transition, no gust | Experimental scheduled case |
| `B2` | Matched speed transition under gust | Experimental scheduled case |
| `C` | Longitudinal trajectory tracking under gust | Deferred; maneuver not frozen |

The Case-B interface retains the complete V90 scheduled runtime, including
scheduled package/history ownership, current-package forecast context, state transport,
scheduled guidance, compiled interval and horizon kernels, condensed RTI,
packet reuse, the accelerated plant interval, rigid-wrench refinement, and
safe endpoint behavior. B1 and B2 fail closed unless the caller explicitly
requests an unqualified research reproduction.

The V90 binding deliberately leaves the future scheduled-package sequence,
online full-fmincon correction, R1 one-correction experiment, and the later
standalone condensation experiments disabled. Those audits were not part of
the executed V90 owner. The resolved `plan.protectedRuntime.auditDisposition`
lists every retained and intentionally disabled owner so later integration
cannot silently lose a speedup or reactivate a superseded experiment.

The clamped wing-only gust-load-alleviation prerequisite remains separate from
formal free-flight Case A. Wing-only uses rate projection off; coupled cases
use their qualified coupled projection policy.

## Dependencies and qualified environment

| Dependency | Needed for | Qualified configuration |
| --- | --- | --- |
| Git and WSL Ubuntu | Clone, source generation, and the Linux-side model workspace | Repository stored in WSL and opened from MATLAB through `\\wsl.localhost` |
| SHARPy with XBeam | Generate or regenerate structural/aerodynamic source models | Existing project SHARPy environment in WSL; upstream sources remain unchanged |
| Windows MATLAB | Setup, simulation, control, post-processing, and tests | R2025b Update 5, 64-bit Windows |
| Control System Toolbox | State-space conversion, LQR design/application, and response analysis | Matching the MATLAB release |
| Optimization Toolbox | `fmincon`, `quadprog`, `lsqnonlin`, and constrained trim/control references | Matching the MATLAB release |
| MATLAB Coder | Build the project-owned native interval and horizon kernels | Matching the MATLAB release |
| Supported C/C++ compiler | Compile MEX acceleration binaries | Microsoft Visual C++ 2022 is the qualified compiler |
| HDF5 support | Exchange unchanged SHARPy model products and physical-output fields | MATLAB HDF5 functions and SHARPy's established Python environment |

SHARPy/XBeam are required when generating model sources, but a packaged
benchmark release can run from its supplied hash-locked MATLAB/HDF5 assets.
MATLAB Coder and a supported C/C++ compiler are required to build the native
tools. Exact MATLAB implementations remain available when compatible binaries
are absent, although scheduled controlled cases can be substantially slower.

Confirm MATLAB toolbox and compiler availability with:

```matlab
ver
mex.getCompilerConfigurations("C++","Selected")
```

If no compiler is selected, run `mex -setup C++` after installing a compiler
supported by the installed MATLAB release.

Native kernels are strongly recommended for controlled scheduled cases. Exact
MATLAB implementations remain available when a compatible binary is absent,
but runtime can be substantially longer.

## Clone and configure MATLAB

Clone the repository in WSL, then open the same directory from Windows MATLAB
through its WSL network path:

```matlab
cd('\\wsl.localhost\Ubuntu\home\<user>\Aircraft-Aeroelastic-Control-Benchmark')
project = setupProject(ChangeCurrentFolder=true);
```

`setupProject` adds only the correct MATLAB package roots. Do not use `genpath`
on `+AeroFlex` or `+RigidBody` directories.

Verify the locked runtime-model payload supplied with a benchmark release:

```matlab
assets = prepareBenchmarkReleaseAssets(Action="check");
assert(assets.passed)
```

Release maintainers can materialize the same hash-locked payload from a
qualified source workspace into a clean packaging root:

```matlab
assets = prepareBenchmarkReleaseAssets( ...
    Action="stage",DestinationRoot="C:\pazy-release-stage");
assert(assets.passed)
```

Staging preserves the runtime-relative layout expected by the qualified
owners, copies only the selected numerical and provenance files, reuses an
existing byte-identical target, and refuses to overwrite a mismatch.

To inspect or stage the complete curated source-and-data package, use
`prepareBenchmarkReleasePackage`. It derives the live MATLAB dependency
closure instead of excluding files solely from historical names. See
[docs/release-package.md](docs/release-package.md) for the supplied-library,
optional-regeneration, history, cache, and validation-example policy.

Verify the Beam, Aero, Core/Base, plant, scheduling, control, trim, registry,
and native-tool integration without running a model:

```matlab
status = verifyBenchmarkInstallation;
assert(status.passed)
```

To require all accelerated components:

```matlab
status = verifyBenchmarkInstallation(RequireNativeKernels=true);
```

## Build the C/C++ MEX components

Select a supported compiler once:

```matlab
mex -setup C++
```

Build every project-owned native component in dependency order:

```matlab
report = buildBenchmarkTools;
assert(report.passed)
```

The build covers:

1. fixed-source reciprocal interval propagation;
2. scheduled reciprocal interval propagation;
3. scheduled estimator/controller horizon and reduced-tangent kernels;
4. scheduled estimator/controller value-horizon kernels;
5. scheduled estimator/controller causal-rollout kernels.

Each builder records the MATLAB release, architecture, compiler, source
signature, binary hash, build time, numerical parity, and timing result.
Caches are separated by MATLAB release, architecture, and source hash. Stale
or incompatible binaries are rejected. A forced clean-source rebuild is:

```matlab
report = buildBenchmarkTools(Force=true,RunParity=true);
```

A read-only cache and toolchain check is:

```matlab
report = buildBenchmarkTools(Action="check");
```

Run the standalone native-tool unit and integration checks with:

```matlab
results = runtests("tests/test_native_tools.m");
assertSuccess(results)
```

This suite builds missing kernels, reuses compatible caches, repeats numerical
parity for all five native families, verifies cache manifests and binary
hashes, and requires the complete accelerated installation to resolve.

Compilation uses two compact, hash-locked fixtures under
`MatlabFlex/configs/benchmark/native-build-fixtures`. They preserve the
accepted code-generation dimensions and parity inputs without depending on a
private simulation checkpoint. The builders verify the fixture hash, schema,
and source-checkpoint provenance before compilation. These fixtures are build
inputs only; benchmark execution continues to use the supplied V17/V17A
runtime library and the selected case configuration.

SHARPy and XBeam are not built, patched, or copied by these commands.

## Run from the polished script

Open [Run_Pazy_Benchmark.m](Run_Pazy_Benchmark.m), edit the short user-settings
section, and run the script. It defaults to plan-only mode so a long simulation
cannot start accidentally. Set `settings.entryMode="benchmark"` for formal
A/B cases or `settings.entryMode="model_workflow"` for wing-only and general
coupled setup/execution.

The programmatic equivalent is:

```matlab
[summary,plan] = runBenchmarkCase("A1",Execute=false);
disp(plan)

summary = runBenchmarkCase("A1", ...
    FiguresVisible=true, ...
    SavePlots=true, ...
    NativeKernelPolicy="required");
```

Experimental Case-B reproduction is deliberate and visibly marked:

```matlab
[~,plan] = runBenchmarkCase("B1", ...
    Execute=false,AllowUnqualified=true);

summary = runBenchmarkCase("B1", ...
    AllowUnqualified=true, ...
    NativeKernelPolicy="required");
```

This override does not qualify B1 or B2 and does not relax a physical,
solver, actuator, or source-domain threshold.

## Custom commands and combined maneuvers

Set `settings.entryMode="custom"` in `Run_Pazy_Benchmark.m`, or construct a
portable custom definition programmatically:

```matlab
reference = struct( ...
    "timeSeconds",[0,5,10], ...
    "speedMps",[15,17,20], ...
    "pitchRad",deg2rad([10,10.5,10]), ...
    "phaseLabel",["attitude","combined","speed"]);

definition = AeroFlex.benchmark.customCaseDefinition( ...
    Name="combined_maneuver", ...
    DurationSeconds=10, ...
    Reference=reference);
[~,plan] = runCustomBenchmarkCase(definition,Execute=false);
```

After reviewing the resolved plan, execute the currently supported scheduled
speed/pitch scope explicitly:

```matlab
summary = runCustomBenchmarkCase(definition, ...
    Execute=true,AllowUnqualified=true);
```

A varying-speed request selects the scheduled runtime for the entire
maneuver. Custom phase labels do not switch between formal runners or reset
plant, estimator, controller, scheduler, fusion, or actuator state. When two
operating points cannot be joined continuously inside the supported domain,
they remain separate manifested runs until a transition is qualified. Formal
A1--A3 and B1/B2 definitions remain unchanged by custom configurations.
Altitude/lateral guidance and nonzero initial perturbations remain fail-closed
until their dedicated runtime owners are qualified.

## Wing-only and coupled open-loop workflows

The retained general setup/execution layer supports these model families:

| Body model | Simulation mode | Purpose |
| --- | --- | --- |
| `wingOnly` | `openloop` | Clamped flexible/aerodynamic propagation |
| `wingOnly` | `nmhe_nmpc` | Clamped gust-load-alleviation control |
| `coupledFull` | `openloop` | No-control coupled trim/replay verification |
| `coupledFull` | `nmhe_nmpc` | Full rigid-flexible estimation and control |

Generate a setup from the unchanged SHARPy outputs:

```matlab
[setup,ok] = sim_init(fullfile(pwd,"TestBenchPazy"), ...
    'case_name',"pazy_krylov_ROM", ...
    'body_case',"wingOnly", ...
    'sim_case',"openloop", ...
    'runner',"PlantROM", ...
    'gustOn',true);
assert(ok)
```

The same operation is available from the centralized interface:

```matlab
[~,plan] = runPazyModelWorkflow( ...
    BodyCase="coupledFull",SimulationMode="openloop", ...
    GustEnabled=false,Execute=false);

result = runPazyModelWorkflow( ...
    BodyCase="wingOnly",SimulationMode="nmhe_nmpc", ...
    GustEnabled=true,Execute=true);
```

Then execute the serialized setup:

```matlab
[ok,history] = sim_run("pazy_krylov_ROM","wingOnly", ...
    'setup_dir',setup.paths.run_dir);
assert(ok)
```

Use `sim_case="nmhe_nmpc"` for the closed wing-only path. For a coupled
no-gust open-loop verification use `body_case="coupledFull"`,
`sim_case="openloop"`, and `gustOn=false`. Setup and execution are separate so
ROM loading, trim, hashes, and native initialization are not charged to online
per-step timing.

## SHARPy-to-MATLAB model workflow

The reproducible source workflow is:

1. generate declared source points through project wrappers in the established
   WSL SHARPy environment;
2. retain unchanged SHARPy/XBeam source trees;
3. export structural, aerodynamic, and premodal HDF5 products;
4. assemble fixed-coordinate and physical-output sidecars in MATLAB;
5. validate each source node and the production registry hashes;
6. construct the wing-only or coupled runtime package;
7. build or verify compatible native kernels;
8. resolve a case plan and execute it.

Source generation is explicit. MATLAB setup never downloads, regenerates, or
silently substitutes a model. Extrapolation is rejected unless a separate
validated study enables it.

## Results, metrics, and plots

Formal facade runs use versioned result directories:

```text
results/<case>/<run-id>/
  manifest.json
  configuration/
  data/
  metrics/
  logs/
  checkpoints/
  plots/diagnostic/
  plots/publication/
  captions/
```

Computational metrics include total wall time, online time, real-time factor,
mean/p95/maximum component times, scheduling, sensing, nMHE, nMPC, allocation,
fusion, actuator, plant, logging, plotting, solver acceptance, fallback, and
checkpoint progress. Initialization, compilation, first-call/JIT, hashing,
serialization, and plotting are reported separately from online timing.

Physical metrics include speed, altitude, attitude and trajectory tracking;
estimator truth/error and gust reconstruction; symmetric/differential wing,
elevator and thrust commands; actuator position/rate/saturation; wingtip
motion; root forces and moments; flexible-state response; load alleviation;
constraint margins; and minimum thrust. Held commands and estimates use stair
plots; propagated truth uses continuous lines.

Plots are toggleable, saved before closure, use deterministic names and units,
and do not call `close all` from reusable production code. Every run writes a
combined benchmark overview plus separate PNG files for airspeed, pitch,
altitude, estimation, wingtip, actuator, thrust, root-load, constraint,
solver, and timing results.

For journal figures, enable publication mode:

```matlab
summary = runBenchmarkCase("A1", ...
    PublicationMode=true,FiguresVisible=false,SavePlots=true);
```

Publication mode removes plot titles, applies consistent article-scale fonts,
and saves each standalone figure as both a 300-dpi PNG and a vector PDF under
`plots/publication/`. Axes, physical units, reference curves, bounds, and
legends remain present so the exported figure is interpretable in a captioned
paper layout.

## Repository organization

- `MatlabFlex/+AeroFlex/+beam`: structural model and recovery operators;
- `MatlabFlex/+AeroFlex/+aero`: aerodynamic ROM and force maps;
- `MatlabFlex/+AeroFlex/+core`: shared indexing, coupling, and assembly;
- `MatlabFlex/+AeroFlex/+sim`: trim, integration, plant, and post-processing;
- `MatlabFlex/+AeroFlex/+sched`: source-library scheduling and transport;
- `MatlabFlex/+AeroFlex/+ctrl`: LQR, nMHE, nMPC, fusion, and actuators;
- `MatlabFlex/+AeroFlex/+benchmark`: case definitions, production runners,
  metrics, plots, and reproducibility manifests;
- `MatlabFlex/+RigidBody`: rigid-body equations and coupled trim;
- `tools/matlab`: reproducible native builds and MATLAB utilities;
- `TestBenchPazy`: local SHARPy/MATLAB exchange and generated run products.

## Troubleshooting

If MATLAB cannot resolve a package, rerun `setupProject` and `rehash path`.
If native verification fails, run `buildBenchmarkTools(Force=true)` and inspect
the component report. If a scheduled source or package is missing, inspect the
release-asset check, registry, and resolved plan; do not substitute a nearby
source or enable extrapolation. Resume long cases only from a hash-matched
checkpoint.

## Standalone validation examples

The [`tests`](tests/README.md) directory includes a quick longitudinal linear
validation with pole, frequency-response, and step-response plots, plus a
checkpointed extended source-trim and flexible-linearization suite. These analyses use
the accepted model products and do not alter the runtime configuration.

## Reproducibility and limitations

Run manifests record the repository revision, MATLAB release, case/profile,
source registry, resolved plan, selected kernels, hashes, timing, outputs, and
qualification state. A visually stable trajectory alone is not a validation
pass.

Case B is distributed as an experimental scheduled workflow with the same
physical, solver, actuator, thrust, and source-domain acceptance thresholds as
the qualified cases. Case C remains deferred until its maneuver definition is
frozen.

## License and attribution

Project-owned source code is released under the
[BSD 3-Clause License](LICENSE). Cite the benchmark paper and archived software
release as described in [CITATION.cff](CITATION.cff). The contributor list is
maintained in [AUTHORS.md](AUTHORS.md), and independently licensed upstream
generation dependencies are identified in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
