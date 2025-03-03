scalar tracer[];

#if 0
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
#endif

face vector tracerf[];
face vector tflux[];

event tracer (i++)
{

    //double dt = timestep(uf, DT);
    //dt = dtnext(dt);
    
    // 1. calculate the flux for each face
    foreach_face() {
#if UPWIND
        tracerf.x[] = u.x[] > 0? tracer[-1]: tracer[];
#else // CENTRAL
        tracerf.x[] = face_value (tracer, 0);
#endif
        tflux.x[] = tracerf.x[] * uf.x[] - mu.x[] * face_gradient_x(tracer, 0);
    }

    // 2. solve the advection-diffusion equation and use forward euler to march in time
    foreach()
        tracer[] -= dt * (tflux.x[1] - tflux.x[] + tflux.y[0,1] - tflux.y[])/Delta;

    boundary ({tracer});
}

