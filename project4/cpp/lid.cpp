#include "include/2D-grid.h"
#include "include/navier-stokes.h"
#include "include/vof.h"

#define U_TOLERANCE 1e-10
double uerror = 0;

double tend = 10;
double RE = 400;
double R = 0.025*M_PI;
const double size = 1;
const double u0 = 1;
int level = 8;

VertexField phi;
VectorField nfg;
ScalarField alphag;

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
        update_boundary(mu);

        FOREACH_VERTEX() 
            phi(i,j) = R - sqrt(sq(grid.xv(i,j) - 0.75) + sq(grid.yv(i,j) - 0.75));
        levelset_to_vof(phi, f);
        reconstruction(f, alphag, nfg);
    }
}
Event init_event (init, "init event", 0);


void output_grid(char* name)
{
    std::ofstream log;
    log.open(name);
    log << "x | y | u.x | u.y | p | e | f\n";
    for (int i = 0; i < grid.col; ++i) {
        for (int j = 0; j < grid.row; ++j) {
            log << grid.x(i,j) << " " << grid.y(i,j) << " " 
                << u.x(i,j) << " " << u.y(i,j) << " " << p(i,j) 
                << " " << error(i,j) << " " << f(i,j) << " " 
                << phi(i,j) << " " << nfg.x(i,j) << " " << nfg.y(i,j)
                << " " << alphag(i,j) << "\n";
        }
        log << "\n";
    }
    log.close();
}


void output (int istep, double t, double dt)
{
    if (istep % 2000 == 0)
    {
        char name[80];
        sprintf(name, "%d-snapshot-%g", istep, t);
        output_grid(name);
        sprintf(name, "%d-interface-%g", istep, t);
        output_interface(f, name);
    }
}
Event output_event (output, "output");


void logfile (int istep, double t, double dt)
{
    double dropv = get_volume (f);
    std::cout << istep << " " << t << " " << dt << " " 
              << Event::allEvents.size() << " " << project.i << " " << project.maxe 
              << " " << uerror << " " << dropv << "\n";
    std::ofstream log;
    if (istep == 0)
    {
        log.open("log");
        log << "1:i | 2:t | 3:dt | 4:events | 5:p.itr | 6:p.maxe | 7:umaxe | 8:volume\n";
    }
    else 
        log.open("log",std::ios_base::app);
    log << istep << " " << t << " " << dt << " " 
        << Event::allEvents.size() << " " << project.i << " " << project.maxe 
        << " " << uerror << " " << dropv << "\n";
}
Event log_event (logfile, "log");

int main(int argc, char** argv)
{
    dtmax = 0.0005;  // for lvl 8
    //dtmax = 0.001; // for lvl 7

    if (argc > 1)
        level = atoi (argv[1]);
    if (argc > 2)
        tend = atof (argv[2]);
    if (argc > 3)
        RE = atof (argv[3]);
    if (argc > 4)
        dtmax = atof (argv[4]);
    if (argc > 5)
        R = atof (argv[5]);

    events.printEvents();
    init_grid(size, size, level);

    run_events();
    events.printEvents();

    char name[80];
    sprintf(name, "final-snapshot");
    output_grid(name);
    sprintf(name, "final-interface");
    output_interface(f, name);

    return 0;
}
