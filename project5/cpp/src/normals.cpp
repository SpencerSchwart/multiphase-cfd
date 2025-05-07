#include "../include/normals.h"

// Youngs finit difference method
Vector youngs_normal (Point pc, ScalarField& c)
{
    int i = pc.i, j = pc.j;

    Vector m;
    m.x = -(c(i+1,j+1) + 2*c(i+1,j) - c(i-1,j+1) - 
            2*c(i-1,j) - c(i-1,j-1) + c(i+1,j-1))/(8*delta);
    m.y = -(c(i+1,j+1) + 2*c(i,j+1) - c(i+1,j-1) -
            2*c(i,j-1) - c(i-1,j-1) + c(i-1,j+1))/(8*delta);
    m.abs_normalize();

    return m;
}


double col_height (Point pc, ScalarField&c, double height)
{
    int stencil = (height - 1)/2; // get index of stencil given i = 0 corresponds to the center cell
    double yi = 0;

    for (int j = -stencil; j <= stencil; ++j)
        yi += c(pc.i, pc.j+j);

    return yi*grid.delta;
}

double row_height (Point pc, ScalarField&c, double height)
{
    int stencil = (height - 1)/2; // get index of stencil given i = 0 corresponds to the center cell
    double xi = 0;

    for (int i = -stencil; i <= stencil; ++i)
        xi += c(pc.i+i, pc.j);

    return xi*grid.delta;
}


double height_function (Point pc, ScalarField&c, double height, HeightType type)
{
    if (type == COLUMN)
        return col_height(pc, c, height);
    else if (type == ROW)
        return row_height(pc, c, height);
    else
        return 0;
}


// centered-column method (height-function method)
Vector height_normal (Point pc, ScalarField& c, double height)
{
    if (c(pc.i, pc.j) == 0 || c(pc.i, pc.j) == 1) // cell is full or empty
        return {0,0};

    // column wise
    double yi1 = height_function({pc.i-1,pc.j}, c, height, COLUMN);
    double yi2 = height_function({pc.i+1,pc.j}, c, height, COLUMN);

    double mxc = -(yi2 - yi1)/(2*grid.delta);

    // row wise
    double xi1 = height_function({pc.i,pc.j-1}, c, height, ROW);
    double xi2 = height_function({pc.i,pc.j+1}, c, height, ROW);

    double myc = -(xi2 - xi1)/(2*grid.delta);

    Vector n;
    if (fabs(mxc) < fabs(myc))
        n.fill(mxc, sign(myc));
    else
        n.fill(sign(mxc), myc);
    
    n.abs_normalize();

    return n;
}


// mixed youngs and central scheme (MYCS)
Vector mycs_normal (Point pc, ScalarField& c, double height)
{
    if (c(pc.i, pc.j) == 0 || c(pc.i, pc.j) == 1) // cell is full or empty
        return {0,0};

    int i = pc.i, j = pc.j;

    // Youngs methods
    Vector ny;
    ny.x = -(c(i+1,j+1) + 2*c(i+1,j) - c(i-1,j+1) - 
            2*c(i-1,j) - c(i-1,j-1) + c(i+1,j-1))/(8*delta);
    ny.y = -(c(i+1,j+1) + 2*c(i,j+1) - c(i+1,j-1) -
            2*c(i,j-1) - c(i-1,j-1) + c(i-1,j+1))/(8*delta);
    
    // Centered-Column method
    // column wise
    double yi1 = height_function({pc.i-1,pc.j}, c, height, COLUMN);
    double yi2 = height_function({pc.i+1,pc.j}, c, height, COLUMN);

    double mxc = -(yi2 - yi1)/(2*grid.delta);

    // row wise
    double xi1 = height_function({pc.i,pc.j-1}, c, height, ROW);
    double xi2 = height_function({pc.i,pc.j+1}, c, height, ROW);

    double myc = -(xi2 - xi1)/(2*grid.delta);

    Vector n;
    if (fabs(mxc) < fabs(myc))
    {
        double mxy = fabs(ny.x)/fabs(ny.y);
        if (fabs(mxy) > fabs(mxc))
            n.fill(ny.x,ny.y);
        else
            n.fill(mxc, sign(myc));
    }
    else
    {
        double myy = fabs(ny.y)/fabs(ny.x);
        if (fabs(myy) > fabs(myc))
            n.fill(ny.x,ny.y);
        else
            n.fill(sign(mxc), myc);
    }

    n.abs_normalize();
    return n;
}
