#include "../include/2D-grid.h"


double delta = 1.;
Grid grid;
Array2D& x = grid.x, y = grid.y;

std::vector<VectorField*> VectorField::allVectorFields;
std::vector<ScalarField*> ScalarField::allScalarFields;
std::vector<FaceVectorField*> FaceVectorField::allFaceVectorFields;


void init_grid (double width, double height, int level)
{
    grid.resize(width, height, level);

    std::cout << "### Initalized Grid ###\n" << 
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

    for (int k = 0; k < (int)VectorField::allVectorFields.size(); ++k)
    {
        VectorField::allVectorFields[k]->resize(grid.row, grid.col, 0.0);
        //std::cout << VectorField::allVectorFields[k]->nx << " " << VectorField::allVectorFields[k]->ny << "\n";
    }

    for (int k = 0; k < (int)ScalarField::allScalarFields.size(); ++k)
    {
        ScalarField::allScalarFields[k]->resize(grid.row, grid.col, 0.0);
    }

    for (int k = 0; k < (int)FaceVectorField::allFaceVectorFields.size(); ++k)
    {
        FaceVectorField::allFaceVectorFields[k]->resize(grid.row, grid.col, 0.0);
    }
}


bool is_ghost_cell (int i, int j)
{
    return (i < ghost || j < ghost || i > grid.col - ghost || j > grid.row - ghost);

}


void no_slip_boundary (VectorField& vf, Boundaries boundary, double val)
{
    vf.boundary.dirichlet_x[boundary] = true;
    vf.boundary.dirichlet_y[boundary] = true;
    vf.boundary.x[boundary] = [&]() -> double {return val;};
    vf.boundary.y[boundary] = [&]() -> double {return val;};
}


void no_slip_boundary (FaceVectorField& vf, Boundaries boundary, double val)
{
    vf.boundary.dirichlet_x[boundary] = true;
    vf.boundary.dirichlet_y[boundary] = true;
    vf.boundary.x[boundary] = [&]() -> double {return val;};
    vf.boundary.y[boundary] = [&]() -> double {return val;};
}


void outlet_boundary (ScalarField& sf, Boundaries boundary)
{
    sf.dirichlet[boundary] = true;
    sf.boundary[boundary] = [&]() -> double {return 0;};
}


void update_boundary_impl (VectorField& vf)
{
    // left boundary
    for (int i = 0; i < ghost; ++i) 
    {
        for (int j = 0; j < grid.row; ++j) 
        {
            if (vf.boundary.dirichlet_x[LEFT])
                vf.x(i,j) = 2*vf.boundary.x[LEFT]() - vf.x(i+1,j);
            else
                vf.x(i,j) = vf.x(i+1,j);

            if (vf.boundary.dirichlet_y[LEFT])
                vf.y(i,j) = 2*vf.boundary.y[LEFT]() - vf.y(i+1,j);
            else
                vf.y(i,j) = vf.y(i+1,j);
        }
    }
    // right boundary
    for (int i = grid.col - ghost; i < grid.col; ++i) 
    {
        for (int j = ghost; j < grid.row; ++j) 
        {
            
            if (vf.boundary.dirichlet_x[RIGHT])
                vf.x(i,j) = 2*vf.boundary.x[RIGHT]() - vf.x(i-1,j);
            else
            vf.x(i,j) = vf.x(i-1,j);

            if (vf.boundary.dirichlet_y[RIGHT])
                vf.y(i,j) = 2*vf.boundary.y[RIGHT]() - vf.y(i-1,j);
            else
                vf.y(i,j) = vf.y(i-1,j);
        }
    }
    // top boundary
    for (int i = 0; i < grid.col; ++i) 
    {
        for (int j = grid.row - ghost; j < grid.row; ++j) 
        {
            if (vf.boundary.dirichlet_x[TOP])
                vf.x(i,j) = 2*vf.boundary.x[TOP]() - vf.x(i,j-1);
            else
                vf.x(i,j) = vf.x(i,j-1);

            if (vf.boundary.dirichlet_y[TOP])
                vf.y(i,j) = 2*vf.boundary.y[TOP]() - vf.y(i,j-1);
            else
                vf.y(i,j) = vf.y(i,j-1);
        }
    }
    // bottom boundary
    for (int i = 0; i < grid.col; ++i) 
    {
        for (int j = 0; j < ghost; ++j) 
        {
            if (vf.boundary.dirichlet_x[BOTTOM])
                vf.x(i,j) = 2*vf.boundary.x[BOTTOM]() - vf.x(i,j+1);
            else
                vf.x(i,j) = vf.x(i,j+1);

            if (vf.boundary.dirichlet_y[BOTTOM])
                vf.y(i,j) = 2*vf.boundary.y[BOTTOM]() - vf.y(i,j+1);
            else
                vf.y(i,j) = vf.y(i,j+1);
        }
    }
}


void update_boundary_impl (ScalarField& sf) 
{
    // left boundary
    for (int i = 0; i < ghost; ++i) 
    {
        for (int j = 0; j < grid.row; ++j) 
        {
            if (sf.dirichlet[LEFT])
                sf(i,j) = 2*sf.boundary[LEFT]() - sf(i+1,j);
            else
                sf(i,j) = sf(i+1,j);
        }
    }
    // right boundary
    for (int i = grid.col - ghost; i < grid.col; ++i) 
    {
        for (int j = ghost; j < grid.row; ++j) 
        {
            if (sf.dirichlet[RIGHT])
                sf(i,j) = 2*sf.boundary[RIGHT]() - sf(i-1,j);
            else
                sf(i,j) = sf(i-1,j);
        }
    }
    // top boundary
    for (int i = 0; i < grid.col; ++i) 
    {
        for (int j = grid.row - ghost; j < grid.row; ++j) 
        {
            if (sf.dirichlet[TOP])
                sf(i,j) = 2*sf.boundary[TOP]() - sf(i,j-1);
            else
                sf(i,j) = sf(i,j-1);
        }
    }
    // bottom boundary
    for (int i = 0; i < grid.col; ++i) 
    {
        for (int j = 0; j < ghost; ++j) 
        {
            if (sf.dirichlet[BOTTOM])
                sf(i,j) = 2*sf.boundary[BOTTOM]() - sf(i,j+1);
            else
                sf(i,j) = sf(i,j+1);
        }
    }
}


void update_boundary_impl (FaceVectorField& vf) 
{
    // left boundary
    for (int i = 0; i <= ghost; ++i)
    {
        for (int j = 0; j < grid.row; ++j)
        {
            if (vf.boundary.dirichlet_x[LEFT])
                vf.x(i,j) = vf.boundary.x[LEFT]();
            else
                vf.x(i,j) = vf.x(i+1,j);

            if (i < ghost)
            {
                if (vf.boundary.dirichlet_y[LEFT])
                    vf.y(i,j) = 2*vf.boundary.y[LEFT]() - vf.y(i+1,j);
                else
                    vf.y(i,j) = vf.y(i+1,j);
            }
        }
    }

    // right boundary
    for (int i = grid.col - ghost; i < grid.col; ++i)
    {
        for (int j = 0; j < grid.row; ++j)
        {
            if (vf.boundary.dirichlet_x[RIGHT])
                vf.x(i,j) = vf.boundary.x[RIGHT]();
            else
                vf.x(i,j) = vf.x(i-1,j);

            if (vf.boundary.dirichlet_y[RIGHT])
                    vf.y(i,j) = 2*vf.boundary.y[RIGHT]() - vf.y(i-1,j);
            else
                    vf.y(i,j) = vf.y(i-1,j);
        }
    }

    // top boundary
    for (int i = 0; i < grid.col; ++i)
    {
        for (int j = grid.row - ghost; j < grid.row; ++j)
        {
            if (vf.boundary.dirichlet_x[TOP])
                vf.x(i,j) = 2*vf.boundary.x[TOP]() - vf.x(i,j-1);
            else
                vf.x(i,j) = vf.x(i,j-1);

            if (vf.boundary.dirichlet_y[TOP])
                    vf.y(i,j) = vf.boundary.y[TOP]();
            else
                    vf.y(i,j) = vf.y(i,j-1);
        }
    }

    // bottom boundary
    for (int i = 0; i < grid.col; ++i)
    {
        for (int j = 0; j <= ghost; ++j)
        {
            if (j < ghost) {
                if (vf.boundary.dirichlet_x[BOTTOM])
                    vf.x(i,j) = 2*vf.boundary.x[BOTTOM]() - vf.x(i,j+1);
                else
                    vf.x(i,j) = vf.x(i,j+1);
            } 
            if (vf.boundary.dirichlet_y[BOTTOM])
                    vf.y(i,j) = vf.boundary.y[BOTTOM]();
            else
                    vf.y(i,j) = vf.y(i,j+1);
        }
    }
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
