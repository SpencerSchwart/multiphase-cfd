#include "grid/multigrid.h"
#define UPWIND 1
#include "tracer-ad.h"

#define L0 1.
#define LEVEL 7

const double end_time = 10;
const double viscosity = 0.;

uf.n[left] = dirichlet(0);
uf.n[right] = dirichlet(0);
uf.n[top] = dirichlet(0);
uf.n[bottom] = dirichlet(0);

uf.t[left] = dirichlet(0);
uf.t[right] = dirichlet(0);
uf.t[top] = dirichlet(0);
uf.t[bottom] = dirichlet(0);

int main()
{
    N = 1 << LEVEL;
    DT = 0.01;
    run();
}

event init (i = 0)
{
    foreach() {
        u.x[] = -sin(M_PI * x) * cos(M_PI * y);
        u.y[] =  cos(M_PI * x) * sin(M_PI * y);
        double term1 = -pow((x - 0.75) / (0.025 * M_PI), 2);
        double term2 =  pow((y - 0.75) / (0.025 * M_PI), 2);
        tracer[] = exp (term1 - term2);
    }

    foreach_face() {
#if UPWIND
        uf.x[] = u.x[] > 0? u.x[-1]: u.x[];
#else
        uf.x[] = face_value (u.x, 0);
#endif
        mu.x[] = viscosity;
    }
    boundary ({uf, tracer});
}

event log (i++; t < end_time)
{
    fprintf (stderr,"%d %g\n", i, t);
}

event output (i += 100)
{
    char name[80];
    sprintf (name, "tracer/%d-tracer-%g.png", i, t);
    output_ppm (tracer, file = name, box = {{0,0},{L0,L0}},
                n = 500, min = 0, max = 0.5);
}
