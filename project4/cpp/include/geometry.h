#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "2D-grid.h"

double line_intersect (Vector n, double alpha, double x = NODATA, double y = NODATA);

int interface_points (Vector n, double alpha, Coord points[2]);

#endif // GEOMETRY_H