#!/bin/sh

echo "Clean project"
sh ./clean.sh

echo "Download external projects"
git submodule init
git submodule update --init --recursive

echo "Create folder build"
mkdir -p build
cd build

echo "Run CMake"
cmake -DCMAKE_BUILD_TYPE=Release ..

NUM_CORES=$(nproc)

echo "Run make with $NUM_CORES parallel jobs"
make -j"$NUM_CORES"

echo "DONE!!!"