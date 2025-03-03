#include "../include/navier-stokes.h"

#define P_TOLERANCE 1e-3
#define N_MAX 500 // maximum iterations of pressure solver

VectorField u;
VectorField us; // intermediate velocity field

FaceVectorField uf;
FaceVectorField mu;

ScalarField p;

GSsolver project; // to hold stats about pressure solver

static void advection_term (int istep, double t, double dt)
{
    VectorField dp;

    FOREACH()
    {
        dp.x(i,j) = central_x (p, i, j, delta);
        dp.y(i,j) = central_y (p, i, j, delta);
    }
    FOREACH_FACE_X()
    {
    #if UPWIND
        double ufAvg = (u.x(i,j) > 0)? u.x(i-1,j): u.x(i,j);
        double dpAvg = (u.x(i,j) > 0)? dp.x(i-1,j): dp.x(i,j);
    #else
        double ufAvg = face_value_x(u.x,i,j); // central
        double dpAvg = face_value_x(dp.x,i,j);
    #endif
        double dpdxf = face_gradient_x(p,i,j,delta);
        uf.x(i,j) = ufAvg - dt*(dpdxf - dpAvg);  // Rhie-Chow correction
    }
    FOREACH_FACE_Y()
    {
    #if UPWIND
        double ufAvg = (u.y(i,j) > 0)? u.y(i,j-1): u.y(i,j); // upwind
        double dpAvg = (u.y(i,j) > 0)? dp.y(i-1,j): dp.y(i,j);
    #else
        double ufAvg = face_value_y(u.y,i,j); // central
        double dpAvg = face_value_y(dp.y,i,j);
    #endif
        double dpdyf = face_gradient_y(p,i,j,delta);
        uf.y(i,j) = ufAvg - dt*(dpdyf - dpAvg);  // Rhie-Chow correction
    }

    update_boundary(uf);

    FaceVectorField flux;
    FOREACH()
    {
        flux.x(i,j) = sq(uf.x(i+1,j)) - sq(uf.x(i,j)) + 
                      uf.y(i,j+1)*face_value_y(u.x,i,j+1) - 
                      uf.y(i,j)*face_value_y(u.x,i,j);
        flux.y(i,j) = sq(uf.y(i,j+1)) - sq(uf.y(i,j)) + 
                      uf.x(i+1,j)*face_value_x(u.y,i+1,j) -
                      uf.x(i,j)*face_value_x(u.y,i,j);

        us.x(i,j) = u.x(i,j) - dt * flux.x(i,j) / delta;
        us.y(i,j) = u.y(i,j) - dt * flux.y(i,j) / delta;
    }
}
Event advection_event (advection_term, "advection term");


void viscous_term (int istep, double t, double dt)
{
    FaceVectorField taux;
    FaceVectorField tauy;

    FOREACH_FACE_X()
    {
        taux.x(i,j) = mu.x(i,j) * 2 * face_gradient_x(u.x,i,j,delta);
        tauy.x(i,j) = mu.x(i,j) * ((u.y(i,j) - u.y(i-1,j))/delta +
                      (u.x(i,j+1) - u.x(i,j-1) + u.x(i-1,j+1) - u.x(i-1,j-1))/(4*delta));
    }
    FOREACH_FACE_Y()
    {
        taux.y(i,j) = mu.y(i,j) * ((u.x(i,j) - u.x(i,j-1))/delta + 
                      (u.y(i+1,j) - u.y(i-1,j) + u.y(i+1,j-1) - u.y(i-1,j-1))/(4*delta));
        tauy.y(i,j) = mu.y(i,j) * 2 * face_gradient_y(u.y,i,j,delta);
    }
    FOREACH()
    {
        double fluxx = taux.x(i+1,j) - taux.x(i,j) + taux.y(i,j+1) - taux.y(i,j);
        us.x(i,j) += dt * fluxx / delta;

        double fluxy = tauy.x(i+1,j) - tauy.x(i,j) + tauy.y(i,j+1) - tauy.y(i,j);
        us.y(i,j) += dt * fluxy / delta;
    }
    update_boundary(us);
}
Event viscous_event (viscous_term, "viscous term");


void correction (VectorField& dp, double dt)
{
    FOREACH()
    {
        u.x(i,j) = us.x(i,j) + dt*dp.x(i,j);
        u.y(i,j) = us.y(i,j) + dt*dp.y(i,j);
    }
}


void centered_gradient (ScalarField& p, VectorField& dp)
{
    FaceVectorField dpf;
    FOREACH_FACE_X()
        dpf.x(i,j) = -face_gradient_x(p,i,j,delta);
    FOREACH_FACE_Y()
        dpf.y(i,j) = -face_gradient_y(p,i,j,delta);
    FOREACH() 
    {
        dp.x(i,j) = (dpf.x(i,j) + dpf.x(i+1,j))/2;
        dp.y(i,j) = (dpf.y(i,j) + dpf.y(i,j+1))/2;
    }
}

ScalarField error;

static void pressure_correction (int istep, double t, double dt)
{
    fill_facefield(uf, u);
    update_boundary(uf);

    ScalarField div;
    FOREACH()
        div(i,j) = (uf.x(i+1,j)-uf.x(i,j)+uf.y(i,j+1)-uf.y(i,j))/(dt*delta);
    
    int iter = 0;
    for(iter = 1; iter < N_MAX; ++iter)
    {
        FOREACH()
            p(i,j) = (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - dv()*div(i,j))/4.;
        update_boundary(p);
        double max_error = -1e30;
        FOREACH()
        {
            error(i,j) = (p(i+1,j) + p(i-1,j) + p(i,j+1) + p(i,j-1) - dv()*div(i,j))/4. - p(i,j);
            if (fabs(error(i,j)) > max_error)
                max_error = error(i,j);
        }
        project.maxe = max_error;
        if (fabs(max_error) < P_TOLERANCE)
            break;
    }
    project.i = iter;

    FOREACH_FACE_X()
        uf.x(i,j) -= dt*face_gradient_x(p,i,j,delta);
    FOREACH_FACE_Y()
        uf.y(i,j) -= dt*face_gradient_y(p,i,j,delta);

    VectorField dp;
    centered_gradient(p, dp);
    correction(dp, dt);

    update_boundary(u,uf);
}
Event pressure_event (pressure_correction, "pressure term");



