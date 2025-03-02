#include "include/2D-grid.h"
#include "include/navier-stokes.h"
#include "include/tracer.h"

#ifndef RE
#define RE 1000.
#endif

//const double endTime = 0.00382;
const double endTime = 10;
const double size = 1;
const double u0 = 1;
const int level = 7;

Vector line = {0.5,0};

void init (int istep, double t, double dt)
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

        dirichlet_boundary (us, LEFT, 0);
        dirichlet_boundary (us, RIGHT, 0);
        dirichlet_boundary (us, BOTTOM, 0);
        dirichlet_boundary_x (us, TOP, u0);
        dirichlet_boundary_y (us, TOP, 0);

        update_boundary(u, uf);

        FOREACH_FACE_X()
            mu.x(i,j) = (u0*size)/RE;
        FOREACH_FACE_Y()
            mu.y(i,j) = (u0*size)/RE;
#if 1
        FOREACH() 
        {
            // initalize tracer field
            double term1 = -pow(((grid.x(i,j)) - 0.75) / (0.025 * M_PI), 2);
            double term2 =  pow(((grid.y(i,j)) - 0.75) / (0.025 * M_PI), 2);
            tracer(i,j) = exp (term1 - term2);
        }
#endif        
        update_boundary(mu);
    }
}
Event init_event (init, "init event", 0);


void snapshot()
{
    std::ofstream snapshot;
    snapshot.open("snapshot");
    Point cut = locate (line.x, line.y);
    std::cout << cut.i << " " << cut.j << "\n";
    for (int j = 0; j < grid.row; ++j)
    {
        Point temp = {cut.i, j};
        //double uxint = interpolate_bi (temp, u.x, line.x, grid.y(cut.i,j));
        snapshot << grid.x(cut.i,j) << " " << grid.y(cut.i,j) << " " 
                 << u.x(cut.i,j) << "\n";
    }

}



void output_grid()
{
    std::ofstream log;
    log.open("log");
    for (int i = 0; i < grid.col; ++i) {
        for (int j = 0; j < grid.row; ++j) {
            log << grid.x(i,j) << " " << grid.y(i,j) << " " 
                << u.x(i,j) << " " << u.y(i,j) << " " << p(i,j) 
                << " " << error(i,j) << " " << tracer(i,j) << "\n";
        }
        log << "\n";
    }
    log.close();
}


int main()
{
    events.printEvents();
    init_grid(size, size, level);

    dtmax = 0.01;
    run_events();
    events.printEvents();

    output_grid();
    snapshot();

    return 0;
}
