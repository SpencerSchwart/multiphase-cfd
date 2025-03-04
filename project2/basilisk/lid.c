#include "grid/multigrid.h"
#include "./navier-stokes.h"
#include "./tracer.h"

#define LEVEL 7
#define L0 1.
#define U0 1.

#ifndef RE
#define RE 100.
#endif

const double tend = 15;
coord line = {0.5,0.5};

u.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);
u.n[top] = dirichlet(0);
u.n[bottom] = dirichlet(0);

u.t[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[top] = dirichlet(U0);
u.t[bottom] = dirichlet(0);

us.n[left] = dirichlet(0);
us.n[right] = dirichlet(0);
us.n[top] = dirichlet(0);
us.n[bottom] = dirichlet(0);

us.t[left] = dirichlet(0);
us.t[right] = dirichlet(0);
us.t[top] = dirichlet(U0);
us.t[bottom] = dirichlet(0);

uf.n[left] = dirichlet(0);
uf.n[right] = dirichlet(0);
uf.n[top] = dirichlet(0);
uf.n[bottom] = dirichlet(0);

uf.t[left] = dirichlet(0);
uf.t[right] = dirichlet(0);
uf.t[top] = dirichlet(U0);
uf.t[bottom] = dirichlet(0);

int main()
{
    N = 1 << LEVEL;
    DT = 0.0005;

    run();
    return 0;
}

event init (i = 0)
{
    foreach() {
        double term1 = -pow((x - 0.75) / (0.025 * M_PI), 2);
        double term2 =  pow((y - 0.75) / (0.025 * M_PI), 2);
        tracer[] = exp (term1 - term2);

        foreach_dimension()
            u.x[] = 0;
    }

    foreach_face()
        mu.x[] = L0 * U0 / RE;
}

event logfile (i++; t <= tend)
{
    fprintf(stderr, "%d %g %d %g\n", i, t, project.i, project.maxe);
}

event output (i += 500)
{
    char name[80];
    sprintf(name, "%d-snapshot-%g", i, t);
    FILE * fp = fopen(name, "w");
    output_field ({u.x, u.y, p, error, tracer}, fp = fp);   
}

event final_output (t = end)
{
    double delta = L0 / (pow(2,LEVEL));

    char name[80];
    sprintf(name, "final-xcut");
    FILE * fp = fopen(name, "w");
    for (double j = 0; j <= L0; j+=delta) {
        foreach_point(0.5, j) {
            double uint = interpolate_linear (point, u.x, line.x, y, 0);
            double vint = interpolate_linear (point, u.y, line.x, y, 0);
            fprintf(fp, "%g %g %g %g %g %g %g\n", 
                         x, y, u.x[], u.y[], uint, vint, p[]);
        }
    }
    sprintf(name, "final-ycut");
    FILE * fp1 = fopen(name, "w");
    for (double i = 0; i <= L0; i+=delta) {
        foreach_point(i, line.y) {
            double uint = interpolate_linear (point, u.x, x, line.y, 0);
            double vint = interpolate_linear (point, u.y, x, line.y, 0);
            fprintf(fp1, "%g %g %g %g %g %g %g\n", 
                         x, y, u.x[], u.y[], uint, vint, p[]);
        }
    }
    sprintf(name, "final-snapshot");
    FILE * fp2 = fopen(name, "w");
    output_field ({u.x, u.y, p, error, tracer}, fp = fp2);

}
