scalar tracer[];

face vector tracerf[];
face vector tflux[];

event tracer (i++)
{
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

