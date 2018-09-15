#!/usr/bin/env bash

mkdir build
cd build
cmake \
  -LAH \
  -DCMAKE_PREFIX_PATH=${PREFIX} \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  ..
VERBOSE=1 make install
#VERBOSE=1 make
