#!/bin/bash

CPPFLAGS="-std=c++17 -Wall -O2"
# CPPFLAGS="-std=c++17 -Wall -pg -g -O2" # for debugging + performance
OUTPUT="lid"

UPWIND=${UPWIND:-0}

RE=${RE:-0}

CONSTANTS="-DUPWIND=$UPWIND -DMU=$MU"

SOURCES="src/timestep.cpp src/events.cpp src/2D-grid.cpp lid.cpp"

g++ $CPPFLAGS $SOURCES $CONSTANTS -o $OUTPUT
