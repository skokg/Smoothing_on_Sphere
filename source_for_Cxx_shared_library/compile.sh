#!/usr/bin/env bash

# After compilation, copy the library beside the Python wrapper in the parent directory.
g++ -std=c++11 -fopenmp -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas -shared -fPIC CC_smoothing_on_sphere_python_lib.cc -o smoothing_on_sphere_Cxx_shared_library.so &&
    printf '%s\n' 'smoothing_on_sphere_Cxx_shared_library.so was created successfully. Copy it to the parent folder before using the Python library.'
