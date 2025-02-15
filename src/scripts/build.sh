#!/bin/bash

CPPFLAGS="-stc=c++17 -Wall -Wetra -02"
OUTPUT="vortex"

SOURCES="vortex.cpp 2D-grid.cpp"

g++ $CPPFLAGS $SOURCES -o $OUTPUT
