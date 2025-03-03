#include "grid/multigrid.h"
#include "./navier-stokes.h"
//#include "./tracer.h"

#define LEVEL 7
#define L0 1.
#define U0 1.
#define RE 1000.

const double tend = 10;

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
    DT = 0.0015;
    run();
    return 0;
}

scalar tracer[];
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
