# Validation examples

This directory contains standalone MATLAB entry scripts for inspecting the
benchmark independently of the formal case runners.

## Clean-clone workflow checks

Run:

```matlab
results = runtests("tests/test_general_model_workflows.m");
assertSuccess(results)
```

The test verifies the eight supplied model files against their SHA-256 manifest,
proves that a missing payload is rejected, and resolves the shared
`wingOnly/openloop`, `wingOnly/nmhe_nmpc`, and `coupledFull/openloop` plans.

## Quick linear validation

Run:

```matlab
run("tests/Run_Linear_Validation.m")
```

The script verifies the accepted compact longitudinal model, compares its
open- and closed-loop poles, evaluates the frequency response, propagates a
small symmetric-wing step, and writes both numerical data and a PNG summary.

## SHARPy wingtip comparison

Run:

```matlab
run("tests/Run_SHARPy_Wingtip_Comparison.m")
```

The script verifies the supplied 40 m/s, 1 degree trim-relative wingtip
dataset, reproduces the accepted SHARPy/MATLAB comparison metrics, and saves a
separate PNG, JSON summary, and MATLAB v7 data file. Set
`comparisonSettings.publicationMode = true` to omit the plot title.

## Extended validation

Run the preflight first:

```matlab
run("tests/Run_Extended_Validation.m")
```

Then select `validationSettings.mode = "full"` in the user-settings section
to verify the locked source hashes, trim ownership, nonlinear residuals, and
flexible/aerodynamic poles. The extended suite writes a checkpoint after each
source so an interrupted run retains its completed evidence.

Neither script changes the benchmark configuration or production source.
Figures are created after the numerical work, saved under `results/validation`,
and closed by handle without affecting unrelated MATLAB figures.

## Native-tool tests

After selecting a C++ compiler, run the public MEX build and parity suite with:

```matlab
mex -setup C++
results = runtests("tests/test_native_tools.m");
assertSuccess(results)
```

The suite builds any missing native families from the supplied hash-locked
fixtures, reuses compatible caches, checks MATLAB release and ABI ownership,
compares every C++ MEX family with its MATLAB implementation, verifies cache
manifests and binary hashes, and finishes with a native-required installation
check. A clean first build can take several minutes; cached checks are faster.
