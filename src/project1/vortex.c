#include "2D-grid.h"
#include "tracer.h"


const double width = 1;
const int level = 6;
const double tEnd = 5;
const double mu = 1;


void init_u()
{
    // initalize velocity field for a non-decaying Taylor-Green vortex
    FOREACH_CELL() {
        u(i, j)->x = -sin(M_PI * (*x(i, j))) * cos(M_PI * (*y(i, j)));
        u(i, j)->y =  cos(M_PI * (*x(i, j))) * sin(M_PI * (*y(i, j)));
    }
}


void init_tracer()
{
    FOREACH_CELL() {
        double term1 = -pow(((*x(i,j)) - 0.75) / (0.025 * M_PI), 2);
        double term2 =  pow(((*y(i,j)) - 0.75) / (0.025 * M_PI), 2);
        *tracer(i,j) = exp (term1 - term2);
    }
}


int main() {
    init_grid (width, width, level);
    init_u();
    run();

    return 0;
}

