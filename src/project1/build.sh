#!/bin/bash

CPPFLAGS="-std=c++17 -Wall -O2"
OUTPUT="vortex"

SOURCES="vortex.cpp 2D-grid.cpp"

g++ $CPPFLAGS $SOURCES -o $OUTPUT
