#ifndef VOF_H
#define VOF_H

#include "2D-grid.h"
#include "methods.h"
#include "utils.h"

extern ScalarField f;

void levelset_to_vof (VertexField& phi, ScalarField& c);

Vector youngs_normal (Point pc, ScalarField& c);

double line_alpha (double c, Vector n);

void reconstruction (ScalarField& c, ScalarField& alpha, VectorField& nf);

int interface_points (Vector n, double alpha, Coord intPoints[2]);

void output_interface (ScalarField& c, char* file);

double rectangle_fraction (Vector nf, double alpha, Coord lp, Coord rp);

void vof_advection (int istep, double t, double dt);

#endif // VOF_H