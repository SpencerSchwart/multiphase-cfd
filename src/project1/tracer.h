#ifndef TRACER
#define TRACER
/*
This header file solves the 2D advection-diffusion equation for a passive
scalar field.
*/

double *tracer; // passive tracer field

extern const double tEnd;
extern const double mu;

inline double* tracer(int i, int j)
{
    return &tracer[i * grid.col + j];
}


void init_tracer();


double tracer_flux (double *tracer, Vector *flux)
{
    FOREACH_CELL() {
        // central difference
        flux[i * grid.col + j].x = *tracer(i,j) * u(i,j)->x - mu * central_x(tracer, i, j, delta);
        flux[i * grid.col + j].y = *tracer(i,j) * u(i,j)->y - mu * central_y(tracer, i, j, delta);

        // forward difference
        //flux[i * grid.col + j].x = *tracer(i,j) * u(i,j)->x - mu * upwind_x(tracer, i, j, delta);
        //flux[i * grid.col + j].y = *tracer(i,j) * u(i,j)->y - mu * upwind_y(tracer, i, j, delta);
    }
}

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
