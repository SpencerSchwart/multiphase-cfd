#include "../include/navier-stokes.h"

VectorField u;
VectorField us; // intermediate velocity field

FaceVectorField uf;
FaceVectorField mu;

ScalarField p;

static void advection_term (int istep, double t, double dt)
{
    VectorField dp;

    // 1.
    FOREACH()
    {
        dp.x(i,j) = central_x (p, i, j, delta);
        dp.y(i,j) = central_y (p, i, j, delta);
    }

    // 2.
    FOREACH_FACE_X()
    {
        //double ufAvg = (u.x(i,j) > 0)? u.x(i-1,j): u.x(i,j); // upwind
        double ufAvg = face_value_x(u.x,i,j); // central
        double dpAvg = face_value_x(dp.x,i,j);
        double dpdxf = face_gradient_x(p,i,j,delta);
        uf.x(i,j) = ufAvg - dt*(dpdxf - dpAvg);  // Rhie-Chow correction
    }
    FOREACH_FACE_Y()
    {
        //double ufAvg = (u.y(i,j) > 0)? u.y(i,j-1): u.y(i,j); // upwind
        double ufAvg = face_value_y(u.y,i,j); // central
        double dpAvg = face_value_y(dp.y,i,j);
        double dpdyf = face_gradient_y(p,i,j,delta);
        uf.y(i,j) = ufAvg - dt*(dpdyf - dpAvg);  // Rhie-Chow correction
    }

    update_boundary(uf);

    FaceVectorField flux;
    // 3.
    FOREACH()
    {
        flux.x(i,j) = (uf.x(i,j)*-uf.x(i,j) + uf.x(i+1,j)*uf.x(i+1,j))*delta;
        flux.y(i,j) = (uf.y(i,j)*-uf.y(i,j) + uf.y(i,j+1)*uf.y(i,j+1))*delta;
        us.x(i,j) = u.x(i,j) -(dt/dv()) * (flux.x(i,j));
        us.y(i,j) = u.y(i,j) -(dt/dv()) * (flux.y(i,j));
    }

    update_boundary(u);
}
Event advection_event (advection_term, "advection term");


static void viscous_term (int istep, double t, double dt)
{
    FaceVectorField flux;
    FaceVectorField du, dv;
    FOREACH_FACE_X()
    {
        du.x(i,j) = face_gradient_x(u.x,i,j,delta);
        dv.x(i,j) = face_gradient_x(u.y,i,j,delta);
    }
    FOREACH_FACE_Y()
    {
        du.y(i,j) = face_gradient_y(u.x,i,j,delta);
        dv.y(i,j) = face_gradient_y(u.y,i,j,delta);
    }
    FOREACH()
    {
        flux.x(i,j) = delta*(mu.x(i+1,j)*2*du.x(i+1,j) -
                             mu.x(i,j)*2*du.x(i,j) +
                             mu.y(i,j+1)*du.y(i,j+1)+dv.x(i,j+1) -
                             mu.y(i,j)*du.y(i,j)+dv.x(i,j));
        us.x(i,j) += (dt/dv()) * flux.x(i,j);

        flux.y(i,j) = delta*(mu.x(i+1,j)*du.y(i+1,j)*dv.x(i+1,j) -
                             mu.x(i,j)*du.y(i,j)*dv.x(i,j) +
                             mu.y(i,j+1)*2*dv.y(i,j+1) -
                             mu.y(i,j)*2*dv.y(i,j));
        us.y(i,j) += (dt/dv()) * flux.y(i,j);                             
    }
}
Event viscous_event (viscous_term, "viscous term");
 

static void pressure_correction (int istep, double t, double dt)
{
    FOREACH_FACE_X()
        uf.x(i,j) = face_value_x(us.x,i,j);
    FOREACH_FACE_Y()
        uf.y(i,j) = face_value_y(us.y,i,j);
    
}
Event pressure_event (pressure_correction, "pressure term");



