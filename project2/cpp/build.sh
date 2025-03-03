#!/bin/bash

CPPFLAGS="-std=c++20 -Wall -O2"
# CPPFLAGS="-std=c++17 -Wall -pg -g -O2" # for debugging + performance
OUTPUT="lid"

UPWIND=${UPWIND:-0}

RE=${RE:-0}

CONSTANTS="-DUPWIND=$UPWIND -DRE=$RE"

SOURCES="src/2D-grid.cpp src/events.cpp src/timestep.cpp src/navier-stokes.cpp src/tracer.cpp lid.cpp"

g++ $CPPFLAGS $SOURCES $CONSTANTS -o $OUTPUT
