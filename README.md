# The Smoothing on Sphere Python Software Package

#### Description:

The package can be used for smoothing of fields defined on a sphere. It includes two appraches with one based on k-d trees and the second on overlap detection. While each approach has its strengths and weaknesses, both are potentially fast enough to make the smoothing of high-resolution meteorological global fields feasible. Detailed description of the methodology is available in the Skok (2024) preprint paper listed in the References section. 

The underlying code is written in C++ and a Python ctypes-based wrapper is provided for easy used in the Python enviroment. A precompiled shared library file is available for easy use with python on Linux systems (the C++ source code is available in the `source_for_Cxx_shared_library` folder - it can be used to compile the shared library for other types of systems). 

#### Usage:

Please refer to the two examples in the files PY_smoothing_on_sphere_example_A_kdtree.py and PY_smoothing_on_sphere_example_B_overlap_detection.py to see how the package can be used in practice. 

#### Author:

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia

Email: Gregor.Skok@fmf.uni-lj.si

#### References:

Arxiv paper
