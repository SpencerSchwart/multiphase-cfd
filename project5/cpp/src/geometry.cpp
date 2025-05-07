#include "../include/geometry.h"


/*
line_intersect returns the x (resp. y) coordinate along a line given y (resp. x).
The return value is in the cell's local coordinate system, i.e. [-0.5,-0.5]x[0.5,0.5].
*/

double line_intersect (Vector n, double alpha, double x, double y)
{
    double inter = 0;
    if (x == NODATA && y != NODATA)
        inter = (alpha - n.y*y) / (n.x+NZERO);
    else
        inter = (alpha - n.x*x) / (n.y+NZERO);
    return inter;
}

// fills points[2] with the coordinates of where the interface intersects the cell's faces
// returns the # of intersecting points
int interface_points (Vector n, double alpha, Coord points[2])
{
    int count = 0;

    // check x faces first
    for (double xint = -0.5; xint <= 0.5; xint += 1) 
    {
        double yint = line_intersect (n, alpha, xint, NODATA);
        if (fabs(yint) <= 0.5) 
        {
            points[count].fill(xint, yint); 
            ++count;
        }
    }

    // then check y faces
    for (double yint = -0.5; yint <= 0.5; yint += 1) 
    {
        double xint = line_intersect (n, alpha, NODATA, yint);
        if (fabs(xint) <= 0.5 && count < 2) 
        {
            points[count].fill(xint, yint); 
            ++count;
        }
    }

    return count;
}
