# Getting started

## MATLAB environment

Use Windows MATLAB R2025b Update 5 against the repository's WSL network
path. The established SHARPy environment remains in WSL and is not replaced
by a Windows Python installation.

Required MATLAB products are Control System Toolbox and Optimization Toolbox.
Building the strongly recommended native kernels additionally requires MATLAB
Coder and a supported C/C++ compiler; Microsoft Visual C++ 2022 is the
qualified Windows compiler. SHARPy/XBeam and their established WSL Python
environment are required for source-model generation, but not for replaying a
release that already contains the verified runtime assets.

```matlab
cd('\\wsl.localhost\Ubuntu\home\<user>\Aircraft-Aeroelastic-Control-Benchmark')
project = setupProject(ChangeCurrentFolder=true);
installation = verifyBenchmarkInstallation(RequireNativeKernels=true, ...
    ProjectInfo=project);
assert(installation.passed)
```

The verification resolves the Beam, Aero, Core/Base, plant, scheduling,
control, actuator, trim, general `sim_init`/`sim_run`, and benchmark-facade
entry points. It also checks the production Case-A/Case-B runners, locked V17
registry, and native binaries.

## Runtime model assets

Check the numerical payload before building or launching a case:

```matlab
assets = prepareBenchmarkReleaseAssets(Action="check");
assert(assets.passed)
```

The release manifest binds every selected source package, certificate,
full-coordinate field, scheduled reciprocal member, controller product, and
case-input contract to its accepted SHA-256. Packaging from a qualified source
workspace uses `Action="stage"` with a clean `DestinationRoot`; changed target
files are never overwritten.

Use `prepareBenchmarkReleasePackage(Action="check")` to inspect the complete
curated source-and-data boundary. The package supplies the verified runtime
library for immediate case reproduction and retains the unchanged external
generation path as an optional workflow. See `docs/release-package.md` for the
history, cache, validation-example, and exclusion policy.

## Native components

Select a compiler and build all project-owned acceleration kernels:

```matlab
mex -setup C++
report = buildBenchmarkTools(Force=false,RunParity=true);
assert(report.passed)
```

Each active cache is specific to the MATLAB release, architecture, and source
signature. `buildBenchmarkTools(Action="check")` is read-only and verifies
the current source hashes as well as the binaries. Exact MATLAB paths remain
available as correctness fallbacks.

## First run

Open `Run_Pazy_Benchmark.m`. The script defaults to plan-only mode. Inspect
the resolved plan, then opt in to execution. Programmatic use is:

```matlab
[~,plan] = runBenchmarkCase("A1",Execute=false);
summary = runBenchmarkCase("A1",NativeKernelPolicy="required");
```

Case B remains experimental and requires `AllowUnqualified=true`. That flag
does not relax any physical or numerical acceptance threshold.

For wing-only or general coupled workflows, select `model_workflow` in the
top-level script or call `runPazyModelWorkflow`. Plan-only mode verifies the
request without loading a model or creating output.
