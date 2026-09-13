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

## Installation and setup

Use MATLAB R2023b or newer with Control System Toolbox and Optimization
Toolbox. Building C++ MEX acceleration additionally requires MATLAB Coder and
a C++ compiler supported by your MATLAB release.

The build accepts 64-bit Windows, Linux, and macOS (Intel or Apple silicon).
Each machine builds its own binaries and compares them with MATLAB results.
Full benchmark validation was performed on Windows MATLAB R2025b Update 5;
native Linux still needs full-case validation. macOS is currently untested
and remains a release limitation.

Clone the repository from a terminal:

```sh
git clone https://github.com/Malexeev123/Aircraft-Aeroelastic-Control-Benchmark.git
cd Aircraft-Aeroelastic-Control-Benchmark
```

In MATLAB, open that directory as the Current Folder. It must contain
`setupProject.m` and `Run_Pazy_Benchmark.m`. Then run:

```matlab
project = setupProject(ChangeCurrentFolder=true);
assets = prepareBenchmarkReleaseAssets(Action="check");
mex -setup C++
report = buildBenchmarkTools;
installation = verifyBenchmarkInstallation(RequireNativeKernels=true);
```

Both build and installation summaries should say `PASS` before running
a controlled case. If a build fails, `report.components` contains its errors.
Run `setupProject` after restarting MATLAB. It adds the package and tool
paths; do not use `genpath` on package folders.

The supplied model data let you run MATLAB without installing SHARPy.
Without MATLAB Coder or a compiler, omit the build and call
`verifyBenchmarkInstallation` without requiring native kernels.
`NativeKernelPolicy="auto"` uses the MATLAB implementations when needed,
although controlled scheduled cases can take substantially longer.

### Operating systems and compilers

| System | MATLAB native build | SHARPy regeneration |
| --- | --- | --- |
| Windows | Supported MinGW or Microsoft Visual C++; prefer a local-drive checkout | Linux through WSL or the upstream container |
| Linux | Supported GCC/G++; WSL is unnecessary | Native Linux installation |
| macOS | Supported Xcode/Clang, matching MATLAB's architecture | Upstream macOS environment, including Fortran and numerical libraries |

Choose a compiler from the [MathWorks requirements for your MATLAB release](https://www.mathworks.com/support/requirements/previous-releases.html).
Inspect the selection with `mex.getCompilerConfigurations("C++","Selected")`.

Windows MATLAB can also open a WSL checkout through its network path.
A local-drive checkout avoids WSL file-notification warnings and network-file
access problems for MATLAB-only runs. Linux and macOS use ordinary local
paths; do not copy a Windows or WSL path from another machine.

### Native tools

`buildBenchmarkTools` builds fixed and scheduled interval kernels, full
and reduced-tangent horizons, value-only horizons, and causal rollouts.
The two fixtures in `MatlabFlex/configs/benchmark/native-build-fixtures`
supply the code-generation dimensions and comparison inputs. Generated C++
and MEX files are local build products.

```matlab
report = buildBenchmarkTools(Action="check"); % inspect existing binaries
report = buildBenchmarkTools(Force=true);     % rebuild and compare with MATLAB
```

Caches are separated by MATLAB release, architecture, and source content.
Rebuild after changing kernel code or moving to another platform.
The builders use MATLAB Coder's compiler setup; a separate project CMake
build is unnecessary. Numerical comparisons remain required before new
binaries are accepted.

## SHARPy setup and library regeneration

This is needed only to regenerate or extend model data. Keep an existing
working SHARPy environment. For a new installation, follow the
[upstream installation guide](https://ic-sharpy.readthedocs.io/en/latest/content/installation.html)
for Linux, macOS, Apple silicon, WSL, or a container. SHARPy builds XBeam
and UVLM with CMake and needs C++, Fortran, Eigen, BLAS, and LAPACK.

On Debian/Ubuntu, the system dependencies can be installed with:

```sh
sudo apt install cmake g++ gfortran libblas-dev liblapack-dev libeigen3-dev
```

Follow the upstream environment and installation steps to build those
libraries. Record the upstream revision and settings with a new source
library; different source versions or settings need fresh validation.

From the benchmark root, in the Python environment containing SHARPy:

```sh
python TestBenchPazy/sweep_pazy_rom_library.py --speed 40 --alpha 1
python TestBenchPazy/sweep_pazy_rom_library.py --speed 40 --alpha 1 --execute --open-loop-reference
```

The first command previews the grid; the second generates it and its
open-loop reference. Repeat `--speed` and `--alpha` for a grid (all combinations).
Outputs go to `TestBenchPazy/generated`, with each source preserved under
`library_source/pazy_krylov_ROM/`. Use `--root` for another output directory.
Keep these files separate from the supplied benchmark models.

Convert each source to a MATLAB setup, then assemble the library:

```matlab
sourceRoot = fullfile(pwd,"TestBenchPazy","generated","library_source", ...
    "pazy_krylov_ROM","pt_U040_alpha_p01");
[setup,ok] = sim_init(sourceRoot, ...
    'case_name',"pazy_krylov_ROM",'body_case',"wingOnly", ...
    'sim_case',"openloop",'runner',"PlantROM",'debug',false);
% Continue only when ok=true. Add other successful setup folders here.
setupDirs = {setup.paths.run_dir};
library = AeroFlex.sched.buildLibraryFromSetups(setupDirs, ...
    'library_name',"my_library", ...
    'save_path',fullfile(pwd,"TestBenchPazy","generated","my_library.mat"));
```

A one-point library permits only that condition. A scheduled library also
needs coordinate compatibility, source-node checks, trim replay, and tests
at intermediate conditions. Assembling files does not qualify interpolation.
Do not overwrite individual production members or edit their checksums to
insert a new model. Rebuilding source data and rebuilding MEX acceleration
are separate steps.

The [notebook](Benchmark.ipynb) and [packaging notes](docs/release-package.md)
describe the retained generation tools and inputs.

## Benchmark cases

| Scenario | Description |
| --- | --- |
| `A1` | Commanded attitude, no gust |
| `A2` | Attitude hold under gust |
| `A3` | Commanded attitude under gust |
| `B1` | Thrust-led speed transition, no gust |
| `B2` | Matched speed transition under gust |
| `C` | Longitudinal trajectory tracking under gust |


Open [Run_Pazy_Benchmark.m](Run_Pazy_Benchmark.m) and edit its user settings.
It starts in plan-only mode. Select `settings.entryMode="benchmark"`
for formal cases, `"custom"` for custom maneuvers, or `"model_workflow"`
for wing-only and general coupled runs. Enable execution after checking the plan.

```matlab
[~,plan] = runBenchmarkCase("A1",Execute=false);
summary = runBenchmarkCase("A1", ...
    FiguresVisible=true,SavePlots=true,NativeKernelPolicy="required");
```

Use A2, A3, B1, or B2 in the same call. Case B retains scheduling, state
transport, guidance, and compiled acceleration. Its full-duration validation
is pending; a successful installation does not close that work. See the
[case definitions](docs/cases.md) for commands and flight conditions.

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
summary = runCustomBenchmarkCase(definition,Execute=true);
```

A varying-speed request selects the scheduled runtime for the entire
maneuver. Custom phase labels do not switch between formal runners or reset
plant, estimator, controller, scheduler, fusion, or actuator state. When two
operating points cannot be joined continuously inside the supported domain,
they remain separate manifested runs until a transition is qualified. Formal
A1--A3 and B1/B2 definitions remain unchanged by custom configurations.
Altitude/lateral guidance and nonzero initial perturbations are not yet
supported by this custom interface.

## Wing-only and coupled open-loop workflows

The retained general setup/execution layer supports these model families:

| Body model | Simulation mode | Purpose |
| --- | --- | --- |
| `wingOnly` | `openloop` | Clamped flexible/aerodynamic propagation |
| `wingOnly` | `nmhe_nmpc` | Clamped gust-load-alleviation control |
| `coupledFull` | `openloop` | No-control coupled trim/replay verification |
| `coupledFull` | `nmhe_nmpc` | Full rigid-flexible estimation and control |

A clean release includes the compact, model inputs needed by
these workflows under `TestBenchPazy/`. `verifyBenchmarkInstallation` checks
all nine files before reporting the general workflow ready. SHARPy/XBeam are
needed only to regenerate or extend those supplied source products.

Generate a setup from the unchanged SHARPy outputs:

```matlab
[setup,ok] = sim_init(fullfile(pwd,"TestBenchPazy"), ...
    'case_name',"pazy_krylov_ROM", ...
    'body_case',"wingOnly", ...
    'sim_case',"openloop", ...
    'runner',"PlantROM", ...
    'gustOn',true,'debug',false);
```

The same operation is available from the centralized interface:

```matlab
[~,plan] = runPazyModelWorkflow( ...
    BodyCase="coupledFull",SimulationMode="openloop", ...
    GustEnabled=false,Execute=false);

result = runPazyModelWorkflow( ...
    BodyCase="wingOnly",SimulationMode="nmhe_nmpc", ...
    GustEnabled=true,DurationSeconds=0.10,Execute=true);
```

By default, the general workflow uses the operating point recorded in the
supplied source package. Set both `TargetSpeedMps` and
`TargetAngleOfAttackDeg` to request an explicit condition; an unscheduled
run fails closed when that request does not match the loaded source data.

Then execute the serialized setup:

```matlab
[ok,history] = sim_run("pazy_krylov_ROM","wingOnly", ...
    'setup_dir',setup.paths.run_dir);
```

Use `sim_case="nmhe_nmpc"` for the closed wing-only path. For a coupled
no-gust open-loop verification use `body_case="coupledFull"`,
`sim_case="openloop"`, and `gustOn=false`. Setup and execution are separate so
ROM loading, trim, hashes, and native initialization are not charged to online
per-step timing.

## Results, metrics, and plots

Formal runs use a separate directory for each execution:

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

## Repository organization and further reading

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


The main installation sequence is above. Detailed references remain in
[case definitions](docs/cases.md), [outputs and plots](docs/outputs-and-plots.md),
[standalone tests](tests/README.md), and [release packaging](docs/release-package.md).

## Troubleshooting

| Symptom | What to check |
| --- | --- |
| `setupProject` unrecognized | Open the checkout containing that file, not its parent or an older checkout. |
| Build tools or packages unrecognized | Rerun `setupProject`; check `which setupProject -all` for another checkout. |
| All native families pass but Linux/macOS overall fails | Older releases had a Windows-only reporting check. Update the build and installation code. |
| No C++ compiler selected | Install a compiler supported by your MATLAB release, then run `mex -setup C++`. |
| Native cache stale or incompatible | Run `buildBenchmarkTools(Force=true)` and inspect `report.components`. |
| Gain path contains backslashes on Linux | Update the path-handling fix and verify the release assets; do not rename data directories. |
| `project_V.r` empty in wing-only setup | Retain `save_pmor_data/pazy_krylov_ROM_krylov_aerorob.h5` alongside `savedata` under the source output. |
| WSL change-notification warning | Use a local-drive checkout if Windows MATLAB access to the WSL share is unreliable. |
| No figure window | Use `FiguresVisible=true` for formal runs, or inspect saved plots when figures are hidden. |
| Error after structural-test figures | These are setup diagnostics, not a completed run. Inspect the following error; 'debug',false omits the diagnostic. |
| Source-domain rejection | Check flight conditions and run diagnostics. Disabling the check does not repair a model inconsistency. |

`prepareBenchmarkReleaseAssets(Action="check")` checks runtime data.
`AeroFlex.benchmark.verifyGeneralModelAssets(pwd,PrintSummary=true)`
checks the nine shared inputs and returns individual file records.
For a machine without a display, use hidden figures and MATLAB `-batch`.
Avoid `-nojvm`: verification uses Java and figure generation needs graphics.

See [troubleshooting details](docs/troubleshooting.md) for cache,
data, and interrupted-run checks.

## Standalone validation examples

The [`tests`](tests/README.md) directory includes a quick longitudinal linear
validation with pole, frequency-response, and step-response plots, plus a
checkpointed extended source-trim and flexible-linearization suite. These analyses use
the accepted model products and do not alter the runtime configuration.

An accepted 40 m/s, 1 degree SHARPy/MATLAB wingtip dataset and standalone
reproduction script are provided under
[`results/validation/sharpy-wingtip-comparison`](results/validation/sharpy-wingtip-comparison)
and [`tests/Run_SHARPy_Wingtip_Comparison.m`](tests/Run_SHARPy_Wingtip_Comparison.m).
The comparison is trim-relative and retains rate projection disabled for the
wing-only replay.

## Reproducibility and limitations

Run manifests record the repository revision, MATLAB release, case/profile,
source registry, resolved plan, selected kernels, hashes, timing, outputs, and
qualification state. A visually stable trajectory alone is not a validation
pass.

Case B is distributed with its production scheduled workflow. Execution is
available from the shared interface, while formal B1/B2 qualification remains
a separate recorded state governed by the same physical, solver, actuator,
thrust, and source-domain acceptance thresholds as the qualified cases. Case C
remains deferred until its maneuver definition is frozen.

## License and attribution

Project-owned source code is released under the
[BSD 3-Clause License](LICENSE). Cite the benchmark paper and archived software
release as described in [CITATION.cff](CITATION.cff). The contributor list is
maintained in [AUTHORS.md](AUTHORS.md), and independently licensed upstream
generation dependencies are identified in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
