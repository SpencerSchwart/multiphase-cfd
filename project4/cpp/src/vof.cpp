#include "../include/vof.h"

ScalarField f;
/*
levelset (phi) notation:
    + = inside
    0 = on interface
    - = outisde
*/

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

// Taken from Basilisk "geometry.h"
// http://basilisk.fr/src/geometry.h#line_alpha
double line_alpha (double c, Vector n)
{
    double alpha = 0;
    n.x = fabs(n.x), n.y = fabs(n.y);

    if (n.x > n.y)
        swap(n.x, n.y);

    c = clamp(c, 0, 1);

    double c1 = n.x/(2*n.y);
    if (c >= 0 && c < c1)
        alpha = sqrt(2*n.x*n.y*c);
    else if (c >= c1 && c <= 0.5)
        alpha = n.y*c + n.x/2.;
    else
        alpha = n.x + n.y - sqrt(2.*n.x*n.y*(1. - c));

    if (n.x < 0.)
        alpha += n.x;
    if (n.y < 0.)
        alpha += n.y;

    return alpha - (n.x + n.y)/2.;
}

// Taken from Basilisk "geometry.h"
// http://basilisk.fr/src/geometry.h#line_area
double line_area (Vector n, double alpha)
{
    double a, v, area;

    alpha += (n.x + n.y)/2.;
    if (n.x < 0.)
    {
        alpha -= n.x;
        n.x = -n.x;
    }
    if (n.y < 0.)
    {
        alpha -= n.y;
        n.y = -n.y;
    }

    if (alpha <= 0.)
        return 0.;

    if (alpha >= n.x + n.y)
        return 1.;
    
    if (n.x < 1e-10)
        area = alpha/n.y;
    else if (n.y < 1e-10)
        area = alpha/n.x;
    else
    {
        v = sq(alpha);
        a = alpha - n.x;
        if (a > 0.)
            v -= a*a;
        a = alpha - n.y;
        if (a > 0.)
            v -= a*a;
        area = v/(2.*n.x*n.y);
    }

    return clamp(area, 0., 1.);
}

void levelset_to_vof (VertexField& phi, ScalarField& c)
{
    FaceVectorField cf;
    FOREACH_FACE_X()
    {
        if (sign2(phi(i,j)) == sign2(phi(i,j+1)))
            cf.x(i,j) = sign2(phi(i,j)) >= 0? 1: 0;
        else if (phi(i,j) == 0) // interface perfectly cuts the vertex
            cf.x(i,j) = sign(phi(i,j+1)) == 1? 1: 0;
        else  // interface cuts face
        {
            cf.x(i,j) = -phi(i,j)/(phi(i,j+1)-phi(i,j));
            if (phi(i,j) < 0) // orient the interface correctly
                cf.x(i,j) = 1 - cf.x(i,j);
        }
    }
    FOREACH_FACE_Y()
    {
        if (sign2(phi(i,j)) == sign2(phi(i+1,j)))
            cf.y(i,j) = sign2(phi(i,j)) >= 0? 1: 0;
        else if (phi(i,j) == 0) // interface perfectly cuts the vertex
            cf.y(i,j) = sign(phi(i+1,j)) == 1? 1: 0;
        else  // interface cuts face
        {
            cf.y(i,j) = -phi(i,j)/(phi(i+1,j)-phi(i,j));
            if (phi(i,j) < 0) // orient the interface correctly
                cf.y(i,j) = 1 - cf.y(i,j);
        }
    }

    FOREACH()
    {
        Vector n;
        n.x = cf.x(i,j) - cf.x(i+1,j);
        n.y = cf.y(i,j) - cf.y(i,j+1);
        n.abs_normalize();

        if (n.abs_mag == 0.)
            c(i,j) = cf.x(i,j);
        else
        {
            double alpha = 0., ni = 0.;
            for (int k = 0; k <= 1; k++)
            {
                if (cf.y(i,j+k) > 0. && cf.y(i,j+k) < 1.)
                {
                    double a = sign(phi(i,j+k))*(cf.y(i,j+k) - 0.5);
                    alpha += n.x*a + n.y*(k - 0.5);
                    ni++;
                }
                if (cf.x(i+k,j) > 0. && cf.x(i+k,j) < 1.)
                {
                    double a = sign(phi(i+k,j))*(cf.x(i+k,j) - 0.5);
                    alpha += n.x*(k - 0.5) + n.y*a;
                    ni++;
                }
            }
            if (ni == 0)
                c(i,j) = max(cf.x(i,j), cf.y(i,j));
            else if (ni != 4)
                c(i,j) = line_area (n, alpha/ni);
            else
                c(i,j) = 0;
        }
    }

}


void reconstruction (ScalarField& c, ScalarField& alpha, VectorField& nf)
{
    FOREACH()
    {
        if (c(i,j) > 0. && c(i,j) < 1.)
        {
            Vector n = youngs_normal({i,j}, c);
            nf.x(i,j) = n.x, nf.y(i,j) = n.y;

            alpha(i,j) = line_alpha (c(i,j), n);
        }
        else
        {
            nf.x(i,j) = 0., nf.y(i,j) = 0.;
            alpha(i,j) = 0.;
        }
    }
}


int interface_points (Vector n, double alpha, Coord intPoints[2])
{
    int count = 0;

    // check on x faces first
    for (double xp = -0.5; xp <= 0.5; xp += 1)
    {
        double yint = (alpha - n.x*xp)/(n.y + SEPS);
        if (fabs(yint) <= 0.5)
        {
            assert(count < 2);
            intPoints[count] = {xp, yint};
            count++;
        }
    }

    // then check y faces
    for (double yp = -0.5; yp <= 0.5; yp += 1)
    {
        double xint = (alpha - n.y*yp)/(n.x + SEPS);
        if (fabs(xint) <= 0.5)
        {
            assert(count < 2);
            intPoints[count] = {xint, yp};
            count++;
        }
    }

    return count;
}


// Taken from Basilisk "fractions.h"
// http://basilisk.fr/src/fractions.h#interface-output
void output_interface (ScalarField& c, char* file)
{
    std::ofstream log;
    log.open(file);

    FOREACH()
    {
        if (c(i,j) > 0. && c(i,j) < 1.)
        {
            double x = grid.x(i,j), y = grid.y(i,j);
            Vector n = youngs_normal({i,j}, c);
            double alpha = line_alpha(c(i,j), n);

            Coord intPoints[2];
            if (interface_points (n, alpha, intPoints) == 2)
            {
                log << x + intPoints[0].x*delta << " " << y + intPoints[0].y*delta << "\n"
                    << x + intPoints[1].x*delta << " " << y + intPoints[1].y*delta << "\n\n";
            }
        }
    }
    log.close();
}


double rectangle_fraction (double c, Coord lp, Coord rp, Vector nf, double alpha)
{

}


static void sweep_x(ScalarField& c, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.x(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(uf.x(i,j));
        int k = s == 1? 0: -1;

        if (c(i+s,j) == 1.)
            cf = 1.;
        else
        {
            Coord lp = {-0.5, -0.5}, rp = {un - 0.5,0.5};
            Vector n(s*nf.x(i,j), nf.y(i,j));
            cf = rectangle_fraction(c(i+k,j), lp, rp, n, alpha(i,j));
        }

        flux(i,j) = cf*uf.x(i,j);
    }
}