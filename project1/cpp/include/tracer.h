#ifndef TRACER_2D_H
#define TRACER_2D_H

#include "2D-grid.h"
#include "methods.h"

/*
This header file solves the 2D advection-diffusion equation for a passive
scalar field.
*/

extern ScalarField tracer; // passive tracer field

extern VectorField u;      // cell center velocity field
extern FaceVectorField uf;     // cell face velocity field

extern FaceVectorField flux;
extern FaceVectorField tracerf;

extern const double mu;

//void solve();


#endif // TRACER_2D_H
