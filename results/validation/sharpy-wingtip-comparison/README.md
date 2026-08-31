# SHARPy wingtip comparison

This directory contains the accepted 40 m/s, 1 degree wing-only gust
comparison between unchanged SHARPy output and the MATLAB reduced-order model.
The comparison uses symmetric displacement relative to each model's initial
equilibrium. This avoids mixing the absolute coordinates of the clamped
SHARPy equilibrium with those of the coupled MATLAB source package.

The supplied CSV and JSON summary are immutable reference data. Reproduce the
metrics, a MATLAB v7 data file, and a separate PNG plot with:

```matlab
run("tests/Run_SHARPy_Wingtip_Comparison.m")
```

Set `comparisonSettings.publicationMode = true` before running to omit the
plot title. New products are written to a timestamped directory beneath
`runs/`, so the supplied reference data are not overwritten.

The wing-only replay retains rate projection disabled. This is the qualified
configuration for comparison with the clamped SHARPy response.
