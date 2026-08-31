# SHARPy model data and source generation

`TestBenchPazy` is the exchange boundary between the unchanged SHARPy/XBeam
generation environment and the MATLAB benchmark. A release contains two
deliberately separate parts:

- a compact, hash-locked runtime payload under `cases/` and `output/`; and
- the project-owned scripts needed to regenerate or extend raw SHARPy ROM
  source points.

Normal MATLAB execution uses the supplied runtime payload. It does not launch
Python and does not require SHARPy to be run first.

## Runtime payload

The required FEM, aerodynamic, saved-state, Krylov ROM, operating-point, and
modal matrix files are declared with SHA-256 hashes in
`MatlabFlex/configs/benchmark/pazy_general_model_assets_v1.json`. Verify them
before a run:

```matlab
setupProject;
status = verifyBenchmarkInstallation;
assert(status.generalWorkflows.modelAssets.passed)
```

Generated MATLAB products are written beneath `sim_setup/` and `sim_run/`.
Each run owns its own `plots/`, `timeseries/`, logs, MAT/HDF5 exchange files,
and compact summary. These directories are local results and are not source
dependencies.

## Generate a SHARPy source point

Source generation runs in WSL using the established SHARPy environment. It
requires Python, NumPy, SciPy, h5py, Matplotlib, SHARPy, and XBeam. Install
SHARPy/XBeam using their official upstream instructions; this repository does
not patch or vendor either dependency.

The sweep driver defaults to plan-only mode:

```bash
python TestBenchPazy/sweep_pazy_rom_library.py --speed 40 --alpha 1
```

After activating the SHARPy environment, generate the ROM source and matched
nonlinear open-loop reference:

```bash
python TestBenchPazy/sweep_pazy_rom_library.py \
  --speed 40 --alpha 1 \
  --open-loop-reference --extract-premodal --execute
```

Repeat `--speed` and `--alpha` to form a Cartesian grid. Use `--fixed-dt` when
a study requires one common source timestep. Existing snapshots are protected;
replacement requires the explicit `--overwrite` option.

The driver performs the following operations for each point:

1. creates the Pazy FEM and aerodynamic inputs through SHARPy's public case
   template;
2. executes static coupling, one dynamic step, modal analysis, linear
   assembly, and Krylov reduction;
3. optionally exports the project-owned premodal physical-coordinate contract;
4. writes MATLAB operating-point metadata; and
5. snapshots the complete point beneath
   `library_source/pazy_krylov_ROM/pt_U...`.

The four `get_settings_*.py` modules retain the Krylov, open-loop/UDP, modal,
and structural SHARPy configurations. `phase18_premodal_extractor.py` is a
project-owned postprocessor and is registered only when `--extract-premodal`
is requested.

## Assemble a MATLAB research library

For every generated point, create a MATLAB setup from that point's own root:

```matlab
pointRoot = fullfile(pwd,"TestBenchPazy","library_source", ...
    "pazy_krylov_ROM","pt_U040_alpha_p01");
[setup,ok] = sim_init(pointRoot, ...
    "case_name","pazy_krylov_ROM", ...
    "body_case","wingOnly", ...
    "sim_case","openloop", ...
    "runner","PlantROM", ...
    "date_only_runs",false, ...
    "run_id","matlab_source_setup");
assert(ok)
```

Collect the resulting setup directories and assemble the library:

```matlab
setupDirectories = [string(setup.paths.run_dir)]; % add every source point
libraryPath = fullfile(pwd,"TestBenchPazy","rom_library", ...
    "custom_pazy_library.mat");
ROMlib = AeroFlex.sched.buildLibraryFromSetups(setupDirectories, ...
    "library_name","custom_pazy_library", ...
    "save_path",libraryPath, ...
    "make_compatible",true);
```

The file is written in MATLAB v7 format. A newly generated library is a
research candidate until exact-node identity, coordinate compatibility,
derivatives, trim replay, open-loop propagation, and applicable scheduling
tests pass. It does not replace the supplied production registry merely
because assembly completed.

After changing source dimensions or rebuilding a scheduled library, rebuild
and recheck the project-owned C/C++ acceleration tools:

```matlab
report = buildBenchmarkTools(Force=true,RunParity=true);
assert(report.passed)
```

## Notebook comparison

[`Benchmark.ipynb`](../Benchmark.ipynb) provides a compact 40 m/s, 1 degree
walkthrough and a SHARPy/MATLAB wingtip comparison. The comparison is
trim-relative: SHARPy's clamped-wing equilibrium and the MATLAB free-flight
source equilibrium are different absolute-coordinate contracts. Wing-only
propagation retains `RateProject` off; coupled-full propagation uses its
qualified projection policy.

The accepted numerical dataset and a standalone MATLAB reproduction script
are also supplied under
`results/validation/sharpy-wingtip-comparison/` and
`tests/Run_SHARPy_Wingtip_Comparison.m`, respectively. This path does not
require rerunning SHARPy.

## Material intentionally not distributed

Local run histories, setup caches, Python bytecode, notebook checkpoints,
Phase-18 campaign diagnostics, rejected interpolation candidates, profiler
captures, and private audit artifacts are not part of the public generation
pipeline. They document development history but are neither runtime inputs nor
necessary to reproduce a new raw SHARPy source point.
