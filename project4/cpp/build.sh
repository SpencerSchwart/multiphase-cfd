#!/bin/bash

# CPPFLAGS="-std=c++20 -Wall -O2"
CPPFLAGS="-std=c++20 -Wall -pg -g -O2" # for debugging + performance
OUTPUT="lid"

UPWIND=${UPWIND:-0} # default is central

RE=${RE:-400}

CONSTANTS="-DUPWIND=$UPWIND -DRE=$RE"

SOURCES="src/2D-grid.cpp src/events.cpp src/timestep.cpp src/utils.cpp src/navier-stokes.cpp src/vof.cpp lid.cpp"

g++ $CPPFLAGS $SOURCES $CONSTANTS -o $OUTPUT
