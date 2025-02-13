#ifndef TIMESTEP
#define TIMESTEP


double cfl = 0.5;
double t = 0;
double dt = 0.1;

extern const double tEnd;
extern const double mu;

void set_timestep (Vector *u, double dx, double mu)
{
    double max = 0;
    for (int i = 0; i < grid.col; ++i) {
        for (int j = 0; j < grid.row; ++j) {
            if (fabs(u(i,j)->x) > maxu && u(i,j)->x > 0) {
                maxu = u(i,j)->x;
            }
            if (fabs(u(i,j)->y) > maxu && u(i,j)->y > 0) {
                maxu = u(i,j)->y;
            }
        }
    }
    
    assert (fabs(maxu) > 0 && "ERROR: All velocity in the domain is zero");
    dt = cfl * dx / fabs(maxu);
    
    if (mu) {
        double dtDiff = dx*dx / (4 * mu);
        dt = min(dt, dtDiff);
    }
}

void solve();

void destroy_grid();

void run()
{
    solve();
    destroy_grid();
}

#endif
