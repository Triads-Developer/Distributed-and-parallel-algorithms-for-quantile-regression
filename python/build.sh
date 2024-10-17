#!/bin/bash

 g++ -O3 -shared -std=c++17 -fPIC $(python3.8 -m pybind11 --includes) qpadmslack.cpp -o qpadmslack$(python3-config --extension-suffix) -Wl,-undefined,dynamic_lookup -I /usr/local/lib/python3.8/dist-packages/numpy/core/include/  -I libraries/carma -I libraries/RcppArmadillo -llapack_atlas

