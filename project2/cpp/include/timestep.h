#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "2D-grid.h"

extern const double cfl;

extern const double timeEnd;

void set_timestep (FaceVectorField& uf, double& dt, double dx, FaceVectorField mu);


#endif // TIMESTEP_H
