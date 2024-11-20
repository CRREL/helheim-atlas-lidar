
# Building on Linux

## Create your environment

    conda create -n atlas
    conda activate atlas

## Install PDAL and a compiler

    conda install -c conda-forge pdal compilers ninja cmake -y

## Build

    ./build.sh

## Run

    ./build/atlas

