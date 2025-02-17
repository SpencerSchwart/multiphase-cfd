#ifndef BASILISK_HEADER_1
#define BASILISK_HEADER_1
#line 1 "./tracer-ad.h"
#include "run.h"

scalar tracer[];

vector u[];

face vector uf[];
face vector mu[];

// is this necessary?
tracer[left] = neumann (0);
tracer[right] = neumann(0);
tracer[top] = neumann(0);
tracer[bottom] = neumann(0);

extern const double viscosity;

double timestep (const face vector u, double dtmax)
{
    static double previous = 0.;
    double delta = 0;
    if (t == 0.) previous = 0.;
    dtmax /= CFL;
    foreach_face() {
        if (u.x[] != 0.) {
            double dt = Delta/fabs(u.x[]);
            if (dt < dtmax) dtmax = dt;
            delta = Delta;
        }
    }
    dtmax *= CFL;

    if (viscosity) {
        double dtvisc = (sq(delta) / (4*viscosity)) / 1.2;
        dtmax = min(dtmax, dtvisc);
    }

    if (dtmax > previous)
        dtmax = (previous + 0.1*dtmax)/1.1;
    previous = dtmax;
    return dtmax;
}


face vector tracerf[];
face vector flux[];

event tracer (i++)
{

    double dt = timestep(uf, DT);
    dt = dtnext(dt);
    
    // 1. calculate the flux for each face
    foreach_face() {
#if UPWIND
        tracerf.x[] = u.x[] > 0? tracer[-1]: tracer[];
#else // CENTRAL
        tracerf.x[] = face_value (tracer, 0);
#endif
        flux.x[] = tracerf.x[] * uf.x[] - mu.x[] * face_gradient_x(tracer, 0);
    }

    // 2. solve the advection-diffusion equation and use forward euler to march in time
    foreach() {
        tracer[] -= dt * (flux.x[1] - flux.x[] + flux.y[0,1] - flux.y[])/Delta;
    }

    boundary ({tracer});
}


#endif
