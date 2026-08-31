# Troubleshooting

## MATLAB path or package resolution

Run `setupProject`, then `rehash path`. Do not recursively add package folders
with `genpath`. `verifyBenchmarkInstallation` reports every unresolved or
shadowed production entry point.

## Native cache rejected

Run:

```matlab
check = buildBenchmarkTools(Action="check");
build = buildBenchmarkTools(Force=true,RunParity=true);
```

A cache is rejected when its MATLAB release, architecture, source hashes,
binary hashes, or parity manifest do not match. The exact MATLAB fallback is
retained, but scheduled controlled cases may become much slower.

## Scheduled source-domain rejection

Do not change interpolation weights, envelopes, or gate limits. Inspect the
saved maximum source-domain ratio, active stencil, scheduler transition,
state-transport history, and plant/provider parity. Resume only after the
reported source and query coordinates agree with the frozen case plan.

## Runtime asset missing or changed

Run:

```matlab
assets = prepareBenchmarkReleaseAssets(Action="check");
```

Inspect the first record whose `passed` field is false. Restore the exact
release asset whose expected SHA-256 is recorded in the manifest; do not use a
nearby operating point, regenerate a single member with different settings,
or disable its source check.

## Long or interrupted run

Use a new output root and resume only from a checkpoint whose configuration,
source registry, runtime profile, and checkpoint hashes match. Initialization,
JIT, loading, hashing, and serialization belong to total wall accounting but
not to the per-control-interval online timing distribution.
