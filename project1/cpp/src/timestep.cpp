#include "../include/timestep.h"


const double cfl = 0.5;


void set_timestep (FaceVectorField& uf, double& dt, double dx, double mu)
{
    double maxu = 0;
    FOREACH()
    {
        if (fabs(uf.x(i,j)) > maxu && uf.x(i,j) > 0) 
            maxu = fabs(uf.x(i,j));

        if (fabs(uf.y(i,j)) > maxu && uf.y(i,j) > 0) 
            maxu = fabs(uf.y(i,j));
    }
    
    assert (maxu > 0 && "ERROR: All velocity in the domain is zero");
    dt = cfl * dx / maxu;
    
    if (mu) 
    {
        double dtVisc = (dx*dx / (4 * mu)) / 1.2; // viscous timestep constraint
        dt = fmin(dt, dtVisc);
    }
}
