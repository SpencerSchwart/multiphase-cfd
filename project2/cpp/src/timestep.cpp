#include "../include/timestep.h"

const double cfl = 0.5;
double dtmax = 0.01;

/*
TODO: add feature to search for max value of mu
*/

double set_timestep (FaceVectorField& uf, FaceVectorField& mu, double dx = delta)
{
    double maxu = -1e30;
    FOREACH_CELL()
    {
        if (fabs(uf.x(i,j)) > maxu && fabs(uf.x(i,j)) > 0) 
            maxu = fabs(uf.x(i,j));

        if (fabs(uf.y(i,j)) > maxu && fabs(uf.y(i,j)) > 0) 
            maxu = fabs(uf.y(i,j));
    }
    
    assert (maxu > 0 && "ERROR: All velocity in the domain is zero");
    double dt = cfl * dx / maxu;
    
    if (mu.x(2,2)) // IMPROVE
    {
        double dtVisc = (sq(dx) / (4 * mu.x(2,2))) / 1.2; // viscous timestep constraint
        dt = fmin(dt, dtVisc);
    }
    return fmin(dt, dtmax);
}

void stability (int istep, double t, double& dt)
{
    if (uf.empty())
        dt = dtmax;
    else {
        dt = set_timestep (uf, mu);
    }
}

// why can't I set dt in an event? scope issue?
//Event timestep_event (stability, "stability");
