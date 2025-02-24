#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "2D-grid.h"

extern const double cfl;
extern FaceVectorField uf;
extern FaceVectorField mu;
extern const double timeEnd;
extern double dtmax;

void set_timestep (FaceVectorField& uf, double& dt, double dx, FaceVectorField mu);

void stability (int istep, double t, double& dt);

#endif // TIMESTEP_H
