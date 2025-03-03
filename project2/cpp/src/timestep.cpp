#include "../include/timestep.h"

const double cfl = 0.5;
double dtmax = 0.0001;
/*
TODO: add feature to search for max value of mu
*/

void set_timestep (FaceVectorField& uf, double& dt, double dx, FaceVectorField& mu)
{
    double maxu = -1e30;
    FOREACH()
    {
        if (fabs(uf.x(i,j)) > maxu && fabs(uf.x(i,j)) > 0) 
            maxu = fabs(uf.x(i,j));

        if (fabs(uf.y(i,j)) > maxu && fabs(uf.y(i,j)) > 0) 
            maxu = fabs(uf.y(i,j));
    }
    
    assert (maxu > 0 && "ERROR: All velocity in the domain is zero");
    dt = cfl * dx / maxu;
    
    if (mu.x(2,2)) // IMPROVE
    {
        double dtVisc = (sq(dx) / (4 * mu.x(2,2))) / 1.2; // viscous timestep constraint
        //std::cout << dt << " " << dtVisc << " " << dx << " " << mu.x(2,2) << "\n";
        dt = fmin(dt, dtVisc);
    }
    dt = fmin(dt, dtmax);
}

void stability (int istep, double t, double& dt)
{
    if (uf.empty())
        dt = dtmax;
    else
        set_timestep (uf, dt, grid.delta, mu);
}

//Event timestep_event (stability, "stability");
