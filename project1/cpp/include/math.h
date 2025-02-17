#ifndef MATH_H
#define MATH_H
#include <cmath>


template <typename T>
inline T sq(T x)
{
    return x * x;
};


struct Coord 
{
    double x;
    double y;
};


struct Vector
{
    double x;
    double y;
    double mag;

    Vector() = default;

    Vector(double xf, double yf)
        : x(xf), y(yf), mag(sqrt(sq(xf) + sq(yf)))
    {
    }

    void normalize()
    {
        mag = sqrt(sq(x) + sq(y));
        x /= mag; y /= mag;
    }
};

#endif // MATH_H
