#include "../include/utils.h"

double get_volume (ScalarField& s)
{
    double volume = 0.;
    FOREACH()
    {
        volume += s(i,j) * sq(grid.delta);
    }
    return volume;
}

double l2_scalar (ScalarField& s0, ScalarField& s)
{
    double l2norm = 0.;
    FOREACH()
    {
        l2norm += sq(s(i,j) - s0(i,j));
    }
    return l2norm;
}


