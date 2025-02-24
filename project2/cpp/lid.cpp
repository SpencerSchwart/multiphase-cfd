//#include "include/2D-grid.h"
#include "include/navier-stokes.h"

const double endTime = 1;
const double size = 1;
const double u0 = 1;
const int level = 7;

static void init (int istep, double t, double dt)
{
    if (istep == 0) {
        // Boundary Conditions
        dirichlet_boundary (uf, LEFT, 0);
        dirichlet_boundary (uf, RIGHT, 0);
        dirichlet_boundary (uf, BOTTOM, 0);
        dirichlet_boundary_x (uf, TOP, u0);
        dirichlet_boundary_y (uf, TOP, 0);

        dirichlet_boundary (u, LEFT, 0);
        dirichlet_boundary (u, RIGHT, 0);
        dirichlet_boundary (u, BOTTOM, 0);
        dirichlet_boundary_x (u, TOP, u0);
        dirichlet_boundary_y (u, TOP, 0);
    }
}
Event init_event (init, "init event", 0);


int main()
{
    events.printEvents();
    init_grid(size, size, level);

    dtmax = 0.01;
    run_events();
    events.printEvents();

    return 0;
}
