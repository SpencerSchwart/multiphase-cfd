#!/bin/bash

CPPFLAGS="-std=c++20 -Wall -O2" # performance
# CPPFLAGS="-std=c++20 -Wall -pg -g -O2" # for debugging + performance
OUTPUT="lid"

UPWIND=${UPWIND:-0} # default is central

# normal methods, default is MYC
YOUNGS=${YOUNGS:-0}
CCM=${CCM:-0}

# VOF advection methods, default is the eulerian implicit method
LE_VOF=${LE_VOF:-0}
BASILISK_VOF=${BASILISK_VOF:-0}

CONSTANTS="-DUPWIND=$UPWIND -DYOUNGS=$YOUNGS -DCCM=$CCM -DLE_VOF=$LE_VOF -DBASILISK_VOF=$BASILISK_VOF"

SOURCES="src/2D-grid.cpp src/events.cpp src/timestep.cpp src/utils.cpp src/navier-stokes.cpp src/geometry.cpp src/normals.cpp src/vof.cpp lid.cpp"

g++ $CPPFLAGS $SOURCES $CONSTANTS -o $OUTPUT
