#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H

#include "2D-grid.h"
#include "methods.h"

extern VectorField u;
extern VectorField us;

extern ScalarField error;

extern FaceVectorField uf;
extern FaceVectorField mu;

extern ScalarField p;

struct GSsolver
{
    int i;       // # of iterations
    double maxe; // maximum error
};

extern GSsolver project;

#endif // NAVIER_STOKES_H