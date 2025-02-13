#include "2D-grid.h"


void init_grid (double width, double height, int level)
{
    grid(width, height, level);

    // initalize x & y coordinates
    FOREACH_CELL() 
    {        
        grid.x(i,j) =  (i * grid.delta - 2*ghost) + grid.delta / 2;
        grid.y(i,j) =  (j * grid.delta - 2*ghost) + grid.delta / 2;
    }
}


bool is_ghost_cell (int i, int j)
{
    return (i < ghost || j < ghost || i > grid.col + 2*ghost || j > grid.row + 2*ghost);

}


void update_boundary_impl (VectorField& vf)
{

    // left boundary
    for (int i = 0; i < ghost; ++i) 
    {
        for (int j = 0; j < grid.row; ++j) 
        {
            if (vf.boundary.dirichhlet_x)
                vf.x(i,j) = 2*vf.boundary.x[LEFT]() - vf.x(i+1,j);
            else
                vf.x(i,j) = vf.boundary.x[LEFT]();

            if (vf.boundary.dirichhlet_y)
                vf.y(i,j) = 2*vf.boundary.y[LEFT]() - vf.y(i+1,j);
            else
                vf.y(i,j) = vf.boundary.y[LEFT]();
        }
    }
}


void update_boundary_impl (ScalarField& sf) {
    // modify sf
}

VectorField f();
VectorField u();
VectorField uf();
ScalarField p();

update_boundary (f, u, uf, p);

VectorField f();
f.boundary.x[LEFT] = [&]() -> double {
    return x + y;
}

if (f.boundary.dirichlet_x)
f.x(i,j) = f.boundary.x[LEFT](); 
