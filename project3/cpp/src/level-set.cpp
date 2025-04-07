#include "../include/level-set.h"


ScalarField tracer; // level-set function
ScalarField heavy;  // heavyside function of level-set funciton
                    // heavy = 0 (inside object)
                    // heavy = 1 (outside object)

FaceVectorField tracerf; // tracer flux

const double dtau = 0.00005; // 0.00005 for 128x128 and 0.000025 for 256x256
double epsilon = 2*delta; // default interface thickness is 2 cell widths
                          // doesn't work right now

extern VectorField u;
extern FaceVectorField uf;

static void tracer_flux (ScalarField& tracer, FaceVectorField& tracerf, FaceVectorField& flux)
{
    // upwind only due to no viscous term
    FOREACH_FACE_X()
    {
        tracerf.x(i,j) = (u.x(i,j) > 0)? tracer(i-1,j): tracer(i,j);
        flux.x(i,j) = tracerf.x(i,j) * uf.x(i,j);
    }

    FOREACH_FACE_Y()
    {
        tracerf.y(i,j) = (u.y(i,j) > 0)? tracer(i,j-1): tracer(i,j);
        flux.y(i,j) = tracerf.y(i,j) * uf.y(i,j);
    }
}


double get_heavyside (ScalarField& phi, int i, int j)
{
    if (phi(i,j) < -epsilon)
        return 0;
    else if (phi(i,j) >= -epsilon && phi(i,j) <= epsilon) // define diffused interface
        return 0.5 + phi(i,j)/(2*epsilon) + sin(M_PI*phi(i,j)/epsilon)/(2*M_PI);
    else
        return 1;
}


static double phi_sign (ScalarField& s, int i, int j, double eps = epsilon)
{
    return s(i,j) / sqrt( sq(s(i,j)) + sq(eps) );
}


static void advect_tracer (int istep, double time, double dt) 
{
    FaceVectorField flux;
    // 1. Calculate flux for each cell
    tracer_flux (tracer, tracerf, flux);

    // 2. Use forward Euler to step in time
    FOREACH()
    {
        tracer(i,j) -= dt*(flux.x(i+1,j) - flux.x(i,j) +
                           flux.y(i,j+1) - flux.y(i,j))/delta;
        heavy(i,j) = get_heavyside(tracer, i, j);
    }

    update_boundary(tracer);
}
static Event event1(advect_tracer, "tracer advection");


static void reinitialize (ScalarField& phi)
{
    ScalarField phi0;
    FOREACH()
        phi0(i,j) = phi(i,j);

    int iter = 0, itermax = 7;
    while (iter < itermax)
    {
        ScalarField phin;
        FOREACH()
            phin(i,j) = phi(i,j);

        FOREACH()
        {
            //double dpdxb = face_gradient_x(phin, i, j);
            //double dpdxf = face_gradient_x(phin, i+1, j);
            //double dpdyb = face_gradient_y(phin, i, j);
            //double dpdyf = face_gradient_y(phin, i, j+1);

            //double dpdx = uf.x(i,j) > 0? dpdxb: dpdxf;
            //double dpdy = uf.y(i,j) > 0? dpdyb: dpdyf;
            double dpdx = central_x (phin, i, j, delta);
            double dpdy = central_y (phin, i, j, delta);
            double gradphi = sqrt(sq(dpdx) + sq(dpdy));

#if 0
            double phixf = max(max(dpdxb, 0), min(dpdxf, 0));
            double phixb = min(min(dpdxb, 0), max(dpdxf, 0));
            double phiyf = max(max(dpdyb, 0), min(dpdyf, 0));
            double phiyb = min(min(dpdyb, 0), max(dpdyf, 0));

            double gradphi = sqrt( max(sq(phixf), sq(phixb)) +
                                   max(sq(phiyf), sq(phiyb)));
#elif 0
            double gradphi = 0;
            if (phi0(i,j) > 0.) 
                gradphi = sqrt( max(sq(max(dpdxb, 0)), sq(min(dpdxf, 0))) +
                                max(sq(max(dpdyb, 0)), sq(min(dpdyf, 0))) );
            else if (phi0(i,j) < 0.)
                gradphi = sqrt( min(sq(min(dpdxb, 0)), sq(max(dpdxf, 0))) +
                                max(sq(min(dpdyb, 0)), sq(max(dpdyf, 0))) );
#endif
            phi(i,j) += dtau*phi_sign(phi0, i, j, epsilon)*(1 - gradphi);
        }
        update_boundary(phi);
        iter++;
    }
}

static void reinitialization (int istep, double time, double dt)
{
    //if (istep % 200 == 0) {
    if (istep % 100 == 0) {
        reinitialize (tracer);
        update_boundary(tracer);
        FOREACH()
        {
            heavy(i,j) = get_heavyside(tracer, i, j);
        }
    }
}
static Event event2(reinitialization, "tracer reinitialization");
