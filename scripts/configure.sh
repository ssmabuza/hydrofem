#!/bin/bash

HYDROFEM_NG="${CODES_DIR}/hydrofem_ng/hydrofem_ng"

rm -rf CMakeCache.txt CMakeFiles
\
cmake \
  -D EIGEN_DIR="${SOFTWARE}/eigen/include/eigen3/"\
${HYDROFEM_NG}
