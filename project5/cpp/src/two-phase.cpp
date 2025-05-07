#include "../include/two-phase.h"

double rho1 = 1, rho2 = 1;
double mu1 = 1, mu2 = 1;

static void properties (int istep, double t, double dt)
{
    FOREACH()
        rho(i,j) = rho1*(f(i,j)) + rho2*(1 - f(i,j));
    FOREACH_FACE_X()
    {
        double ff = (f(i,j) + f(i-1,j))/2.;
        alpha.x(i,j) = 1/(rho1*ff + rho2*(1 - ff));
        mu.x(i,j) = mu1*ff + mu2*(1 - ff);
    }
    FOREACH_FACE_X()
    {
        double ff = (f(i,j) + f(i,j-1))/2.;
        alpha.y(i,j) = 1/(rho1*ff + rho2*(1 - ff));
        mu.y(i,j) = mu1*ff + mu2*(1 - ff);
    }
}
Event properties_event (properties, "two_phase_properties", 1);