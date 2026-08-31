# Release package boundary

The public benchmark uses a hybrid distribution. Verified, hash-locked ROM
and controller products are supplied so the predefined cases can be run
without rebuilding the aerodynamic library. The unchanged SHARPy/XBeam
generation path and project-owned MATLAB wrappers remain documented for users
who want to add flight conditions, modes, or scheduled source points.

Ordinary simulation histories are outputs and are not source dependencies.
The package may include compact reference metrics or selected comparison data
when they are required for a regression or a published figure, but it does not
include local run directories, checkpoints, profiler sessions, or campaign
logs.

The staged package includes `AUTHORS.md`, `CITATION.cff`, and
`THIRD_PARTY_NOTICES.md`. SHARPy and XBeam are obtained separately under their
upstream BSD 3-Clause licenses; neither source tree is copied into the public
MATLAB package.

## Curated source staging

From a qualified development workspace, inspect the package boundary with:

```matlab
status = prepareBenchmarkReleasePackage(Action="check");
```

The command computes the MATLAB dependency closure of the documented entry
points and combines it with the locked runtime-asset manifest and the eight
explicitly declared general-workflow model files. This is
important because a small number of production numerical owners retain
historical filenames; they are included when the executable call graph
requires them and are not discarded by a filename heuristic. Dynamic MAT/HDF5
loads are manifest-owned because static dependency analysis cannot discover
their filenames.

After reviewing the package inventory, stage a clean candidate into an empty
directory outside the source repository:

```matlab
status = prepareBenchmarkReleasePackage( ...
    Action="stage",DestinationRoot="C:\work\pazy-benchmark-release");
assert(status.publicationReady)
```

The staging command never overwrites a file with different bytes. It writes a
machine-readable inventory containing every selected source and runtime asset
hash.

The staged `.gitattributes` preserves those byte-locked hashes across Git
clones, including accepted MATLAB sources whose established line endings are
intentionally unchanged by release packaging.
Required MAT/HDF5 libraries are tracked explicitly; `results/` and the
general-workflow run directories contain generated outputs and remain visible
through tracked README files. The compact SHARPy/MATLAB wingtip dataset,
summary, and figure under `results/validation/sharpy-wingtip-comparison/` are
the sole selected reference-result exception. Timestamped reproduction runs,
caches, and code-generation directories remain excluded by the public
`.gitignore`.

## Clean native build

The package contains every project-owned builder and dynamically selected
kernel source required by `buildBenchmarkTools`. It also contains two compact
native-build fixtures. The fixtures are extracted from accepted runtime
checkpoints and preserve only the type, dimension, and parity inputs needed by
MATLAB Coder; ordinary checkpoints and simulation histories are not release
dependencies.

From a newly staged package, select a supported C++ compiler and run:

```matlab
setupProject;
report = buildBenchmarkTools(Force=true,RunParity=true);
assert(report.passed)

status = verifyBenchmarkInstallation(RequireNativeKernels=true);
assert(status.passed)
```

Every builder verifies the fixture hash, schema, and source-checkpoint
provenance before generating code. The resulting MEX files are local cache
products and are not part of the curated source inventory. Exact MATLAB
implementations remain the fallback when compatible binaries are unavailable.

## Excluded material

The public package omits private development context and conversations,
campaign-only audit scripts, checkpoints and logs, ordinary simulation
histories, rejected numerical candidates, generated MEX/codegen caches,
Python bytecode and environments, editor state, credentials, machine-specific
configuration, and publisher PDFs without verified redistribution rights.

The selected standalone linear-response and extended-source validation
scripts remain included because they document the pole, frequency-response,
step-response, source identity, and physical-linearization checks used by the
benchmark. Long validation runs write their results to versioned output
directories rather than the source tree.

The release also retains the plan-safe SHARPy sweep, its four settings
modules, the optional premodal extractor, and the notebook walkthrough.
Historical source-placement campaigns, error-budget studies, profiler data,
and their generated source grids are excluded because they are development
evidence rather than inputs to the supported generation or runtime workflow.
