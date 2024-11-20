#!/bin/bash


rm -rf build
mkdir build
cd build
cmake -G Ninja ..  \
    -DCMAKE_BUILD_TYPE=Release  \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
ninja
