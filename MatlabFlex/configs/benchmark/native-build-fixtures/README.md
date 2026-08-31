# Native build fixtures

These compact MAT files provide the type, dimension, and parity inputs needed
to compile the optional project-owned MEX kernels from a clean release. They
are extracted from hash-locked accepted runtime checkpoints; ordinary runtime
histories are not distributed.

The fixture loader verifies each file hash, schema, and source-checkpoint hash
before a builder can use it. Both files use MATLAB v7 MAT format for reliable
Windows MATLAB access through the WSL repository path.
