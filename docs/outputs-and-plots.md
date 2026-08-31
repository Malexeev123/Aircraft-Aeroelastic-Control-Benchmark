# Outputs, metrics, and plots

Each facade run writes a new versioned directory:

```text
results/<case>/<run-id>/
  manifest.json
  configuration/
  data/standardized_results.mat
  metrics/metrics.json
  metrics/metrics.csv
  logs/
  checkpoints/
  plots/diagnostic/
  plots/publication/<case>_benchmark_trends.png
  plots/publication/<case>_airspeed.png
  plots/publication/<case>_pitch.png
  plots/publication/<case>_<signal>.png
  captions/
```

The manifest records the resolved case plan, qualification state, repository
revision, MATLAB environment, production-registry hash, and native-kernel
policy.
MAT artifacts use MATLAB v7 format.

The standardized overview contains airspeed, pitch, altitude, gust truth and
estimate, state-estimation error, wingtip response, actuator positions and
rates, thrust, root force and bending moment, active-source gate,
solver acceptance, and component timings. Propagated truth uses continuous
lines; held commands and estimates use stairs.

Each overview panel is also saved as a separate figure for inspection and
article assembly. Set `settings.publicationMode=true` in
`Run_Pazy_Benchmark.m`, or pass `PublicationMode=true` to `runBenchmarkCase`,
to remove panel titles and add vector PDF copies while retaining axis labels,
units, legends, bounds, and reference curves.

Metrics include RMS, peak, and terminal tracking errors where available;
root loads and wingtip response; actuator use; minimum thrust; source-domain
ratio; nMHE/nMPC success and fallback; and mean, p95, and maximum online
component time. Preparation, plot generation, metric assembly, and
serialization are reported separately from online component timing.
`onlineWallSeconds` is assembled from recurring scheduling, sensing,
estimation, control, allocation, fusion, actuator, and plant work.
`runnerExecutionWallSeconds` retains the complete numerical-owner call for
diagnosis, including any runner-side loading or final checks. Neither setup nor
plot generation is silently charged to the online real-time factor.

If a Case-B run reaches its final qualification
assertion, the facade recovers the runner's saved MAT artifact and generates
the same standardized products while retaining a `VALIDATION_PENDING_FAIL...`
status. Other execution errors are recorded and rethrown.

General wing-only and coupled workflows retain their established `sim_run`
directory and physical post-processing products. The centralized
`runPazyModelWorkflow` entry adds `model_workflow_summary.mat` beside those
products. It records setup time, simulation-plus-post-processing wall time, a
conservative real-time factor, available per-component timing distributions,
the physical metric structure, and exact plot/output locations. Because the
legacy runner performs post-processing before returning, this conservative
wall time is not presented as setup-excluded online time.
