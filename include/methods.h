

inline double upwind_x (double *s, int i, int j, double dx)
{
    int k = -sign(*s(i,j));

    return k * ((*s(i,j) ) - (*s(i - k,j))) / dx;
}


inline double upwind_y (double *s, int i, int j, double dy)
{
    int k = -sign(*s(i,j));

    return k * ((*s(i,j) ) - (*s(i,j - k))) / dy;
}


inline double central_x (double *s, int i, int j, double dx)
{
    return (*s(i + 1,j) - *s(i - 1, j)) / (2 * dx);
}


inline double central_y (double *s, int i, int j, double dy)
{
    return (*s(i,j + 1) - *s(i, j + 1)) / (2 * dy);
}

