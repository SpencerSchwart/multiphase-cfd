#include "run.h"

#define P_TOLERANCE 1e-3
#define N_MAX 500

vector u[];
vector us[];
face vector uf[];

scalar p[];

face vector mu[];

typedef struct Solvert
{
    int i;
    double maxe;
} Solvert;

Solvert project;


event stability (i++)
{
    double dt = DT;
    dt = dtnext(dt);
}


event advection_term (i++)
{
    vector dp[];
    foreach()
        foreach_dimension()
            dp.x[] = center_gradient(p);

    foreach_face() {
    #if UPWIND
        double ufavg = u.x[] > 0? u.x[-1]: u.x[];
        double dpavg = u.x[] > 0? dp.x[-1]: dp.x[];
    #else
        double ufavg = face_value(u.x,0);
        double dpavg = face_value(dp.x,0);        
    #endif
        double dpdxf = face_gradient_x(p,0);
        uf.x[] = ufavg - dt*(dpdxf - dpavg);  // Rhie-Chow correction
    }
    boundary({uf});

    vector flux[];
    foreach() {
        double utf = (u.x[0,1] + u.x[])/2.; // should also include Rhie-Chow correction
        double ubf = (u.x[] + u.x[0,-1])/2.;
        flux.x[] = sq(uf.x[1]) - sq(uf.x[]) + uf.y[0,1]*utf - uf.y[]*ubf;
        
        us.x[] = u.x[] - dt*flux.x[]/Delta;

        double vrf = (u.y[1] + u.y[])/2.;
        double vlf = (u.y[] + u.y[-1])/2.;
        flux.y[] = sq(uf.y[0,1]) - sq(uf.y[]) + uf.x[1]*vrf - uf.x[]*vlf;

        us.y[] = u.y[] - dt*flux.y[]/Delta;
    }
}


event viscous_term (i++)
{
    face vector taux[], tauy[];
    foreach_face(x) {
        taux.x[] = mu.x[]*2*face_gradient_x(u.x,0);
        tauy.x[] = mu.x[]*((u.y[] - u.y[-1])/Delta +
                           (u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(4*Delta));
    }
    foreach_face(y) {
        taux.y[] = mu.y[]*((u.x[] - u.x[0,-1])/Delta +
                           (u.y[1] - u.y[-1] + u.y[1,-1] - u.y[-1,-1])/(4*Delta));
        tauy.y[] = mu.y[]*2*face_gradient_y(u.y,0);
    }

    foreach() {
        double fluxx = taux.x[1] - taux.x[] + taux.y[0,1] - taux.y[];
        us.x[] += dt * fluxx / Delta;

        double fluxy = tauy.x[1] - tauy.x[] + tauy.y[0,1] - tauy.y[];
        us.y[] += dt * fluxy / Delta;
    } 
    boundary({us});
}


void correction (vector dp, double dt)
{
    foreach() {
        foreach_dimension()
            u.x[] = us.x[] + dt * dp.x[];
    }
}


void centered_gradient (scalar p, vector dp)
{
    face vector dpf[];
    foreach_face()
        dpf.x[] = -face_gradient_x(p,0);

    foreach()
        foreach_dimension()
            dp.x[] = (dpf.x[] + dpf.x[1])/2;
}

scalar error[];

event pressure_correction (i++)
{
    foreach_face()
        uf.x[] = face_value(u.x,0);

    scalar div[];
    foreach()
        div[] = (uf.x[1] - uf.x[] + uf.y[0,1] - uf.y[])/(dt*Delta);
     
    int iter = 0;
    for (iter = 1; iter < N_MAX; ++iter) {
        foreach()
            p[] = (p[1] + p[-1] + p[0,1] + p[0,-1] - dv()*div[])/4.;
        boundary({p});

        double maxerror = -1e30;
        foreach() {
            error[] = (p[1] + p[-1] + p[0,1] + p[0,-1] - dv()*div[])/4. - p[];
            if (fabs(error[]) > maxerror)
                maxerror = error[];
        }
        project.maxe = maxerror;
        if (fabs(maxerror) < P_TOLERANCE)
            break;
    }
    project.i = iter;

    foreach_face()
        uf.x[] -= dt*face_gradient_x(p,0);

    vector dp[];
    centered_gradient(p, dp);
    correction(dp, dt);

    boundary({u, uf});
}

