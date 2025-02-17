#include "../include/tracer.h"


ScalarField tracer; // passive tracer field

VectorField u;      // cell center velocity field
FaceVectorField uf;     // cell face velocity field

FaceVectorField flux;
FaceVectorField tracerf;


static void initial_field (FaceVectorField& uf, ScalarField& tracer, 
                           VectorField& u, FaceVectorField& tracerf)
{
    // set boundary conditions 
    no_slip_boundary (uf, LEFT, 0);
    no_slip_boundary (uf, RIGHT, 0);
    no_slip_boundary (uf, TOP, 0);
    no_slip_boundary (uf, BOTTOM, 0);
   
    FOREACH() 
    {
        // initalize tracer field
        double term1 = -pow(((grid.x(i,j)) - 0.75) / (0.025 * M_PI), 2);
        double term2 =  pow(((grid.y(i,j)) - 0.75) / (0.025 * M_PI), 2);
        tracer(i,j) = exp (term1 - term2);

        // initalize velocity field
        u.x(i, j) = -sin(M_PI * (grid.x(i, j))) * cos(M_PI * (grid.y(i, j)));
        u.y(i, j) =  cos(M_PI * (grid.x(i, j))) * sin(M_PI * (grid.y(i, j)));
    }

    FOREACH_FACE_X()
    {
#if UPWIND
        uf.x(i,j) = (u.x(i,j) > 0)? u.x(i-1,j): u.x(i,j);
#else // central
        uf.x(i,j) = face_value_x (u.x, i, j);
#endif    
    }

    FOREACH_FACE_Y()
    {
#if UPWIND
        uf.y(i,j) = (u.y(i,j) > 0)? u.y(i,j-1): u.y(i,j);
#else // central
        uf.y(i,j) = face_value_y (u.y, i, j);
#endif
    }

    update_boundary(uf);
}



static void tracer_flux (ScalarField& tracer, FaceVectorField& tracerf, FaceVectorField& flux,
                         FaceVectorField& uf, double mu)
{
    FOREACH_FACE_X()
    {
#if UPWIND
        tracerf.x(i,j) = (u.x(i,j) > 0)? tracer(i-1,j): tracer(i,j);
#else        
        tracerf.x(i,j) = face_value_x (tracer, i, j);
#endif        

        flux.x(i,j) = tracerf.x(i,j) * uf.x(i,j) - mu * face_flux_x(tracer, i, j, grid.delta);
    }

    FOREACH_FACE_Y()
    {
#if UPWIND
        tracerf.y(i,j) = (u.y(i,j) > 0)? tracer(i,j-1): tracer(i,j);
#else
        tracerf.y(i,j) = face_value_y (tracer, i, j);
#endif        

        flux.y(i,j) = tracerf.y(i,j) * uf.y(i,j) - mu * face_flux_y(tracer, i, j, grid.delta);
    }
}


static void solve (int istep, double time, double dt) 
{
    // 0. set B.C for uf and flux fields
    if (istep == 0)
        initial_field (uf, tracer, u, tracerf);

    // 1. Calculate flux for each cell
    tracer_flux (tracer, tracerf, flux, uf, mu);

    // 2. Use forward Euler to step in time
    FOREACH()
    {
        tracer(i,j) -= dt*(flux.x(i+1,j) - flux.x(i,j) +
                           flux.y(i,j+1) - flux.y(i,j))/grid.delta;
    }

    update_boundary(tracer);
}


static SetEvent event(solve);


