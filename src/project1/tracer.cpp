#include "../../include/tracer.h"


static void initial_field (ScalarField& tracer, FaceVectorField& uf)
{ 
    uf.boundary.x[LEFT] = [&]() -> double {return 0;};
    uf.boundary.x[RIGHT] = [&]() -> double {return 0;};
    uf.boundary.x[BOTTOM] = [&]() -> double {return 0;};

    uf.boundary.y[LEFT] = [&]() -> double {return 0;};
    uf.boundary.y[RIGHT] = [&]() -> double {return 0;};
    uf.boundary.y[BOTTOM] = [&]() -> double {return 0;};
}


static void get_face_velocity (VectorField& u, FaceVectorField& uf)
{
    FOREACH()
    {
        // upwind
        uf.x(i,j) = (u.x(i,j) > 0)? u.x(i-1,j): u.x(i,j);
        uf.y(i,j) = (u.y(i,j) > 0)? u.y(i,j-1): u.y(i,j);

        // central
        // uf.x(i,j) = face_value_x (u.x, i, j);
        // uf.y(i,j) = face_value_y (u.y, i, j);
    }
}


static void tracer_flux (ScalarField& tracer, FaceVectorField& flux, FaceVectorField& uf, double mu)
{
    FOREACH() 
    {
        flux.x(i,j) = tracer(i,j) * uf.x(i,j) - mu * face_flux_x(tracer, i, j, grid.delta);
        flux.y(i,j) = tracer(i,j) * uf.y(i,j) - mu * face_flux_y(tracer, i, j, grid.delta);
    }
}


static void solve (int i, double time, double dt) 
{
    initial_field (tracer, uf);

    // 1. Calculate the face velocity
    get_face_velocity (u, uf);

    //update_boundary(uf);

    // 2. Calculate flux for each cell
    FaceVectorField flux;
    tracer_flux (tracer, flux, uf, mu);
    
    // 3. Use forward Euler to step in time
    FOREACH()
    {
        tracer(i,j) -= dt*((flux.x(i+1,j) - flux.x(i,j))/grid.delta +
                           (flux.y(i,j+1) - flux.y(i,j))/grid.delta);
    }

    std::cout << "In solver! i = " << i << " time = " << time << " dt = " << dt << "\n";
}


static SetEvent event(solve);


#if 0
void solve()
{
    tracer = (double*) malloc(grid.col * grid.row * sizeof(double));
    init_tracer();

    Vecotr *flux = (Vector*) malloc((grid.col + 2*ghost) * (grid.row + 2*ghost) * sizeof(Vector));

    for (int _i = 0; t <= tEnd; ++_i) {
        set_timestep(grid.u, dx, mu);
        
        tracer_flux (tracer, flux);

        FOREACH() {
            
        }

        t += dt;
    }

    free(tracer); tracer = NULL;
}
#endif
