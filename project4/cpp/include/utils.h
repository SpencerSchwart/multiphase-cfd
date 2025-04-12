#ifndef UTILS_H
#define UTILS_H

#include "2D-grid.h"

template <typename T, typename R>
inline int min (T a, R b)
{
    return a < b? a: b;
}

template <typename T, typename R>
inline int max (T a, R b)
{
    return a > b? a: b;
}

double get_volume (ScalarField& s);

double l2_scalar (ScalarField& s0, ScalarField& s);

#endif // UTILS_H
