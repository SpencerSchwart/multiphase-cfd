#ifndef METHODS_H
#define METHODS_H


template <typename T>
inline double upwind_x (T& field, int i, int j, double dx)
{
    int k = -sign(field(i,j));

    return k * (field(i,j) - field(i - k,j)) / dx;
}


template <typename T>
inline double upwind_y (T& field, int i, int j, double dy)
{
    int k = -sign(field(i,j));

    return k * ((field(i,j) ) - (field(i,j - k))) / dy;
}


template <typename T>
inline double central_x (T& field, int i, int j, double dx)
{
    return (field(i + 1,j) - field(i - 1, j)) / (2 * dx);
}


template <typename T>
inline double central_y (T& field, int i, int j, double dy)
{
    return (field(i,j + 1) - field(i, j - 1)) / (2 * dy);
}


template <typename T>
inline double face_flux_x (T& field, int i, int j, double dx)
{
    return (field(i,j) - field(i-1,j)) / dx;
}


template <typename T>
inline double face_flux_y (T& field, int i, int j, double dy)
{
    return (field(i,j) - field(i,j-1)) / dy;
}


template <typename T>
inline double face_value_x (T& field, int i, int j)
{
    return (field(i,j) + field(i-1,j)) / 2.;
}


template <typename T>
inline double face_value_y (T& field, int i, int j)
{
    return (field(i,j) + field(i,j-1)) / 2.;
}

#endif // METHODS_H
