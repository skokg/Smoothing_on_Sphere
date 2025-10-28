# The Smoothing on Sphere Python Software Package

#### Description:

The package can be used for smoothing fields defined on a sphere. It includes two approaches, one based on k-d trees and the second on overlap detection. While each approach has its strengths and weaknesses, both are potentially fast enough to make the smoothing of high-resolution meteorological global fields feasible. A detailed description of the methodology is available in the paper listed in the References section. 

The underlying code is written in C++, and a Python ctypes-based wrapper is provided for easy use in the Python environment. A precompiled shared library file is available for easy use with Python on Linux systems (the C++ source code is available in the `source_for_Cxx_shared_library` folder - it can be used to compile the shared library for other types of systems). 

#### Usage:

To see how the package can be used in practice, please refer to the examples in the files:
1) PY_smoothing_on_sphere_example_A_kdtree.py - smoothing a single field using kdtree-based methodology without saving the smoothing data.
2) PY_smoothing_on_sphere_example_A_kdtree_multiple_fields_simultaneously.py  - smoothing multiple fields simultaneously using kdtree-based methodology without saving the smoothing data. This is usually significantly faster than smoothing each field separately.
3) PY_smoothing_on_sphere_example_A_kdtree_using_saved_smoothing_data.py - smoothing a single field using kdtree-based methodology by first generating the smoothing data that is later used for the smoothing. This is significantly faster than smoothing without using the smoothing data; however, it can require a substantial amount of memory and disk space for the smoothing data.
4) PY_smoothing_on_sphere_example_B_overlap_detection.py - smoothing a single field using overlap-based methodology by first generating the smoothing data that is later used for the smoothing.
5) PY_smoothing_on_sphere_example_C_missing_data.py - smoothing a single field with missing data using kdtree-based methodology.

#### Author:

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia

Email: Gregor.Skok@fmf.uni-lj.si

#### References:

Skok, G. and Kosovelj, K., 2025: Smoothing and spatial verification of global fields, Geosci. Model Dev., 18, 7417–7433, https://doi.org/10.5194/gmd-18-7417-2025.
