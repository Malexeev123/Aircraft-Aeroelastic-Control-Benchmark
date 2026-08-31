# Third-party software notices

The benchmark uses external software to generate aeroelastic source data. The
public MATLAB package does not redistribute or modify the source trees listed
below; users install them separately and remain responsible for complying with
their licenses.

## SHARPy

SHARPy (Simulation of High Aspect Ratio aeroplanes in Python) is developed by
the Aeroelastics Group at Imperial College London and distributed under the
BSD 3-Clause License.

- Project: <https://github.com/ImperialCollegeLondon/sharpy>
- Documentation: <https://ic-sharpy.readthedocs.io/>

## XBeam

XBeam is the geometrically exact beam library used by SHARPy. The upstream
project is maintained by Imperial College London and distributed under the
BSD 3-Clause License.

- Project: <https://github.com/ImperialCollegeLondon/xbeam>

The benchmark distributes a compact set of generated SHARPy/XBeam model data,
but not either upstream source tree. The generated runtime libraries, MATLAB
interfaces, and native acceleration kernels are project-owned release
products. Their hashes and provenance are recorded in the release inventory
and model-library manifests.
