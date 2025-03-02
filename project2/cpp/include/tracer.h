#ifndef TRACER_2D_H
#define TRACER_2D_H

#include "2D-grid.h"
#include "methods.h"
/*
This header file solves the 2D advection-diffusion equation for a passive
scalar field.
*/

extern ScalarField tracer; // passive tracer field

extern FaceVectorField flux;
extern FaceVectorField tracerf;

//void solve();


#endif // TRACER_2D_H
