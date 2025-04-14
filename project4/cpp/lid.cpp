#include "include/2D-grid.h"
#include "include/navier-stokes.h"
#include "include/vof.h"

#ifndef RE
#define RE 400.
#endif

#define R 0.025*M_PI

#define U_TOLERANCE 1e-10
double uerror = 0;

const double tend = 0.002;
const double size = 1;
const double u0 = 1;
const int level = 8;

Vector line = {0.5,0.5};

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


void output_xcut(char* name)
{
    std::ofstream snapshot;
    snapshot.open(name);
    snapshot << "x | y | u.x | u.y | uint | vint | p\n";
    Point cut = locate (line.x, 0);
    for (int j = ghost; j < grid.row - ghost; ++j)
    {
        Point temp = {cut.i, j};
        double uxint = interpolate_bi (temp, u.x, line.x, grid.y(cut.i,j));
        double uyint = interpolate_bi (temp, u.y, line.x, grid.y(cut.i,j));
        snapshot << grid.x(cut.i,j) << " " << grid.y(cut.i,j) << " " 
                 << u.x(cut.i,j) << " " << u.y(cut.i,j) << " " 
                 << uxint << " " << uyint << " " << p(cut.i,j) << "\n";
    }
}

void output_ycut(char* name)
{
    std::ofstream snapshot;
    snapshot.open(name);
    snapshot << "x | y | u.x | u.y | uint | vint | p\n";
    Point cut = locate (0, line.y);
    for (int i = ghost; i < grid.col - ghost; ++i)
    {
        Point temp = {i, cut.j};
        double uxint = interpolate_bi (temp, u.x, grid.x(i,cut.j), line.y);
        double uyint = interpolate_bi (temp, u.y, grid.x(i,cut.j), line.y);
        snapshot << grid.x(i,cut.j) << " " << grid.y(i,cut.j) << " " 
                 << u.x(i,cut.j) << " " << u.y(i,cut.j) << " " 
                 << uxint << " " << uyint << " " << p(i,cut.j) << "\n";
    }
}


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


VectorField un;
void converge_check (int istep, double t, double dt)
{
    double maxdiff = -1e30;
    FOREACH()
    {
        double diffx = un.x(i,j) - u.x(i,j);
        double diffy = un.y(i,j) - u.y(i,j);
        if (fabs(diffx) > maxdiff)
            maxdiff = diffx;
        if (fabs(diffy) > maxdiff)
            maxdiff = diffy;
    }
    uerror = maxdiff;
    if (fabs(maxdiff) < U_TOLERANCE && istep > 0) 
        stop = true;

    FOREACH()
    {
        un.x(i,j) = u.x(i,j);
        un.y(i,j) = u.y(i,j);
    }
}
Event converge (converge_check, "u_converge");

int main()
{
    events.printEvents();
    init_grid(size, size, level);

    //dtmax = 0.0005;
    dtmax = 0.001;
    run_events();
    events.printEvents();

    char name[80];
    sprintf(name, "final-snapshot");
    output_grid(name);
    sprintf(name, "final-xcut");
    output_xcut(name);
    sprintf(name, "final-ycut");
    output_ycut(name);
    sprintf(name, "final-interface");
    output_interface(f, name);

    return 0;
}
