#include "../../include/2D-grid.h"


double delta = 1.;
Grid grid;
Array2D& x = grid.x, y = grid.y;

extern VectorField pf;

void init_grid (double width, double height, int level)
{
    grid.resize(width, height, level);

    std::cout << "### Initalized Gride ###\n" << 
                 "width | height | lvl | row | col | delta\n" <<
                 grid.width << "\t" << grid.height << "\t" << grid.level << "\t" <<
                 grid.row << "\t" << grid.col << "\t" << grid.delta << "\n";

    delta = grid.delta;
    // initalize x & y coordinates
    FOREACH_CELL() 
    {
        grid.x(i,j) =  i*delta + delta/2 - ghost*delta;
        grid.y(i,j) =  j*delta + delta/2 - ghost*delta;
    }

    pf.resize(grid.width, grid.height, 0.0);
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
            if (vf.boundary.dirichlet_x)
                vf.x(i,j) = 2*vf.boundary.x[LEFT]() - vf.x(i+1,j);
            else
                vf.x(i,j) = vf.boundary.x[LEFT]();

            if (vf.boundary.dirichlet_y)
                vf.y(i,j) = 2*vf.boundary.y[LEFT]() - vf.y(i+1,j);
            else
                vf.y(i,j) = vf.boundary.y[LEFT]();
        }
    }
}


void update_boundary_impl (ScalarField& sf) {
    // modify sf
}

#if 0
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
#endif
