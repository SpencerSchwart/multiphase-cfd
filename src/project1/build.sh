#!/bin/bash

# CPPFLAGS="-std=c++17 -Wall -O2 -fsanitize=address -g"
CPPFLAGS="-std=c++17 -Wall -O2"
OUTPUT="vortex"

SOURCES="events.cpp 2D-grid.cpp tracer.cpp vortex.cpp"

g++ $CPPFLAGS $SOURCES -o $OUTPUT
