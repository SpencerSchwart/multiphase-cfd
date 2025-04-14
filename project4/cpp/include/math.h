#ifndef MATH_H
#define MATH_H
#include <cmath>

#define SEPS 1e-30  // small epsilon

template <typename T>
inline T sq(T x)
{
    return x * x;
}


template <typename T>
inline int sign(T x)
{
    return x < 0? -1: 1;
}


template <typename T>
inline int sign2(T x)
{
    
    return x < 0? -1: (x > 0? 1: 0);
}


struct Coord 
{
    double x;
    double y;
};


struct Vector
{
    double x, y;
    double mag, abs_mag;

    Vector() = default;

    Vector(double xf, double yf)
        : x(xf), y(yf), mag(sqrt(sq(xf) + sq(yf)) + SEPS), abs_mag(fabs(xf) + fabs(yf) + SEPS)
    {
    }

    void normalize()
    {
        mag = sqrt(sq(x) + sq(y)) + SEPS;
        x /= mag; y /= mag;
    }

    void abs_normalize()
    {
        abs_mag = fabs(x) + fabs(y) + SEPS;
        x /= abs_mag; y /= abs_mag;
    }
};

#endif // MATH_H
