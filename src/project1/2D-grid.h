#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "timstep.h"


typedef struct Vector 
{
    double x;
    double y;
}


typedef struct Grid 
{
    double width, height;
    int level;

    double delta;
    int row, col;

    double *x, *y;

    Vector *u;
    double *p;
} Grid;


const int ghost = 1;    // layers of ghost cells around a domain boundary
double delta = 1;
Grid grid;


inline double* x(int i, int j)
{
    return &grid.x[i * grid.col + j];
}


inline double* y(int i, int j)
{
    return &grid.y[i * grid.col + j];
}


inline Vector* u(int i, int j)
{
    return &grid.u[i * grid.col + j];
}


inline double* p(int i, int j)
{
    return &grid.p[i * grid.col + j];
}


#define FOREACH_CELL() \
    for (int i = 0; i < grid.col + 2 * ghost; ++i) \
        for (int j = 0; j < grid.row + 2 * ghost; ++j)

#define FOREACH() \
    for (int i = ghost; i < grid.col + ghost; ++i) \
        for (int j = ghost; j < grid.row + ghost; ++j)


void init_grid (double width, double height, int level)
{
    grid.width = width;
    grid.height = height;
    grid.level = level;

    delta = grid.delta = max(width, height) / ((double) pow(2, level));
    grid.col = grid.row = pow(2, level);

    grid.x = (double*) malloc((grid.col + 2*ghost) * (grid.row + 2*ghost) * sizeof(double));
    grid.y = (double*) malloc((grid.col + 2*ghost) * (grid.row + 2*ghost) * sizeof(double));

    grid.u = (Vector*) malloc((grid.col + 2*ghost) * (grid.row + 2*ghost) * sizeof(Vector));
    grid.p = (double*) malloc((grid.col + 2*ghost) * (grid.row + 2*ghost) * sizeof(double));
    
    // initalize x & y coordinates and the velocity and pressure fields
    FOREACH_CELL() {        
        *x(i,j) =  (i * grid.delta - ghost) + grid.delta / 2;
        *y(i,j) =  (j * grid.delta - ghost) + grid.delta / 2;
        u(i,j)->x = 0.;
        u(i,j)->y = 0.;
        *p(i,j) = 0.;
    }
}


bool is_ghost_cell (int i, int j)
{
    return (i < ghost || j < ghost || i > grid.col + 2*ghost || j > grid.row + 2*ghost);
}

void destroy_grid ()
{
    if (grid.x) { free(grid.x); grid.x = NULL; }
    if (grid.y) { free(grid.y); grid.y = NULL; }
    if (grid.u) { free(grid.u); grid.u = NULL; }
    if (grid.p) { free(grid.p); grid.p = NULL; }
}




