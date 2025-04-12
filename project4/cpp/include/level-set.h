#ifndef TRACER_2D_H
#define TRACER_2D_H

#include "2D-grid.h"
#include "methods.h"
#include "utils.h"

/*
This header file solves the 2D advection equation for a passive
scalar field represented by a level-set function.
*/

extern ScalarField tracer; // passive tracer field
extern ScalarField heavy;

extern FaceVectorField flux;
extern FaceVectorField tracerf;

extern double epsilon;

double get_heavyside (ScalarField& phi, int i, int j);

#endif // TRACER_2D_H
