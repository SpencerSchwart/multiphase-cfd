#include "../include/tracer.h"


ScalarField tracer; // passive tracer field

FaceVectorField tracerf;

extern VectorField u;

static void tracer_flux (ScalarField& tracer, FaceVectorField& tracerf, FaceVectorField& flux)
{
    FOREACH_FACE_X()
    {
#if UPWIND
        tracerf.x(i,j) = (u.x(i,j) > 0)? tracer(i-1,j): tracer(i,j);
#else        
        tracerf.x(i,j) = face_value_x (tracer, i, j);
#endif        
        flux.x(i,j) = tracerf.x(i,j) * uf.x(i,j) - mu.x(i,j) * face_gradient_x(tracer, i, j);
    }

    FOREACH_FACE_Y()
    {
#if UPWIND
        tracerf.y(i,j) = (u.y(i,j) > 0)? tracer(i,j-1): tracer(i,j);
#else
        tracerf.y(i,j) = face_value_y (tracer, i, j);
#endif        
        flux.y(i,j) = tracerf.y(i,j) * uf.y(i,j) - mu.y(i,j) * face_gradient_y(tracer, i, j);
    }
}


static void solve (int istep, double time, double dt) 
{
    FaceVectorField flux;
    // 1. Calculate flux for each cell
    tracer_flux (tracer, tracerf, flux);

    // 2. Use forward Euler to step in time
    FOREACH()
    {
        tracer(i,j) -= dt*(flux.x(i+1,j) - flux.x(i,j) +
                           flux.y(i,j+1) - flux.y(i,j))/delta;
    }

    update_boundary(tracer);
}


static Event event(solve, "tracer");


