#include "../include/vof.h"

ScalarField f;
/*
levelset (phi) notation:
    + = inside
    0 = on interface
    - = outisde
*/

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

Vector interface_normal (Point pc, ScalarField& c)
{
#if YOUNGS
    return youngs_normal(pc, c);
#elif CCM
    return height_normal(pc, c);
#else // MYC
    return mycs_normal(pc, c);
#endif
}


void reconstruction (ScalarField& c, ScalarField& alpha, VectorField& nf)
{
    FOREACH()
    {
        if (c(i,j) > 0. && c(i,j) < 1.)
        {
            Vector n = interface_normal ({i,j}, c);
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

// Taken from Basilisk "fractions.h"
// http://basilisk.fr/src/fractions.h#interface-output
void output_interface (ScalarField& c, char* file)
{
    std::ofstream log;
    log.open(file);

    FOREACH()
    {
        if (c(i,j) > 0.+INT_TOL && c(i,j) < 1.-INT_TOL)
        {
            double x = grid.x(i,j), y = grid.y(i,j);
            Vector n = interface_normal ({i,j}, c);
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


double rectangle_fraction (Vector nf, double alpha, Coord lp, Coord rp)
{
    Vector nf1 = {nf.x*(rp.x - lp.x), nf.y*(rp.y - lp.y)};
    alpha -= nf.x*(lp.x + rp.x)/2.;
    alpha -= nf.y*(lp.y + rp.y)/2.;

    return line_area (nf1, alpha);
}

#if BASILISK_VOF
static void sweep_x(ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.x(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i+k,j) >= 1. || c(i+k,j) <= 0.)
            cf = c(i+k,j);
        else
        {
            Coord lp = {-0.5, -0.5}, rp = {s*un - 0.5,0.5};
            Vector n(-s*nf.x(i+k,j), nf.y(i+k,j));
            cf = rectangle_fraction(n, alpha(i+k,j), lp, rp);
        }

        flux(i,j) = cf*uf.x(i,j);
    }

    FOREACH()
        c(i,j) += dt*(flux(i,j) - flux(i+1,j) + cc(i,j)*(uf.x(i+1,j) - uf.x(i,j)))/(delta);
}


static void sweep_y (ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.y(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i,j+k) >= 1. || c(i,j+k) <= 0.)
            cf = c(i,j+k);
        else
        {
            Coord lp = {-0.5, -0.5}, rp = {s*un - 0.5,0.5};
            Vector n(-s*nf.y(i,j+k), nf.x(i,j+k));
            cf = rectangle_fraction(n, alpha(i,j+k), lp, rp);
        }

        flux(i,j) = cf*uf.y(i,j);
    }

    FOREACH()
        c(i,j) += dt*(flux(i,j) - flux(i,j+1) + cc(i,j)*(uf.y(i,j+1) - uf.y(i,j)))/(delta);
}
#endif

static void ei_sweep_x(ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.x(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i+k,j) >= 1. || c(i+k,j) <= 0.)
            cf = c(i+k,j);
        else
        {
            Coord lp = {-0.5, -0.5}, rp = {s*un - 0.5,0.5};
            Vector n(-s*nf.x(i+k,j), nf.y(i+k,j));
            cf = rectangle_fraction(n, alpha(i+k,j), lp, rp);
        }

        flux(i,j) = dt*cf*uf.x(i,j)/delta;
    }

    FOREACH()
        c(i,j) = (c(i,j) + flux(i,j) - flux(i+1,j)) /
                 (1 - dt*(uf.x(i+1,j) - uf.x(i,j))/delta);
}


static void ei_sweep_y (ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.y(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i,j+k) >= 1. || c(i,j+k) <= 0.)
            cf = c(i,j+k);
        else
        {
            Coord lp = {-0.5, -0.5}, rp = {s*un - 0.5,0.5};
            Vector n(-s*nf.y(i,j+k), nf.x(i,j+k));
            cf = rectangle_fraction(n, alpha(i,j+k), lp, rp);
        }

        flux(i,j) = dt*cf*uf.y(i,j)/delta;
    }

    FOREACH()
        c(i,j) = (c(i,j) + flux(i,j) - flux(i,j+1)) /
                 (1 - dt*(uf.y(i,j+1) - uf.y(i,j))/delta);

}


// rectangle_fraction but with two mapping procedures
// 1. map interface line (n*x = alpha) to cell defined by lp and rp
// 2. compute volume fraction of area enclosed by lp1 and rp1
double map_rectangle_fraction (Vector nf, double alpha, Coord lp, Coord rp, Coord lp1, Coord rp1)
{
    Vector nf1 = {nf.x*(rp.x - lp.x), nf.y*(rp.y - lp.y)};
    alpha -= nf.x*(lp.x + rp.x)/2.;
    alpha -= nf.y*(lp.y + rp.y)/2.;

    return rectangle_fraction (nf1, alpha, lp1, rp1);
}

#if LE_VOF
static void le_sweep_x(ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.x(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i+k,j) >= 1. || c(i+k,j) <= 0.)
            cf = c(i+k,j);
        else
        {
            if (un == 0)
                continue;
            double unl = uf.x(i-s,j)*dt/grid.delta;
            Coord lp, rp, lp1, rp1;
            if (s == 1) 
            {
                lp.fill (-0.5 + unl, -0.5), rp.fill(0.5 + un, 0.5);
                lp1.fill(0.5, -0.5), rp1.fill(0.5 + un,0.5);
            }
            else
            {
                lp.fill (-0.5 + un, -0.5), rp.fill (0.5 + unl, 0.5);
                lp1.fill(-0.5 + un, -0.5), rp1.fill(-0.5, 0.5);
            }
            Vector n(-s*nf.x(i+k,j), nf.y(i+k,j));
            cf = map_rectangle_fraction(n, alpha(i+k,j), lp, rp, lp1, rp1);
            //std::cout << c(i+k,j) << " " << cf << std::endl;
        }

        flux(i,j) = dt*cf*uf.x(i,j)/delta;
    }

    FOREACH()
        c(i,j) = c(i,j)*(1 + dt*(uf.x(i+1,j) - uf.x(i,j))/delta) + flux(i,j) - flux(i+1,j);
}


static void le_sweep_y(ScalarField& c, ScalarField& cc, double dt)
{
    ScalarField alpha, flux;
    VectorField nf;

    reconstruction(c, alpha, nf);

    FOREACH()
    {
        double un = uf.y(i,j)*dt/grid.delta, cf;
        assert (un < cfl && "ERROR: CFL condition is violated\n");

        int s = sign(un);
        int k = s == 1? -1: 0;

        if (c(i,j+k) >= 1. || c(i,j+k) <= 0.)
            cf = c(i,j+k);
        else
        {
            double unl = uf.y(i,j-s)*dt/grid.delta;
            Coord lp, rp, lp1, rp1;
            if (s == 1) 
            {
                lp.fill (-0.5 + unl, -0.5), rp.fill(0.5 + un, 0.5);
                lp1.fill(0.5, -0.5), rp1.fill(0.5 + un,0.5);
            }
            else
            {
                lp.fill (-0.5 + un, -0.5), rp.fill (0.5 + unl, 0.5);
                lp1.fill(-0.5 + un, -0.5), rp1.fill(-0.5, 0.5);
            }
            Vector n(-s*nf.y(i,j+k), nf.x(i,j+k));
            cf = map_rectangle_fraction(n, alpha(i,j+k), lp, rp, lp1, rp1);
        }

        flux(i,j) = dt*cf*uf.y(i,j)/delta;
    }

    FOREACH()
        c(i,j) = c(i,j)*(1 + dt*(uf.y(i,j+1) - uf.y(i,j))/delta) + flux(i,j) - flux(i,j+1);
}
#endif

void vof_advection (int istep, double t, double dt)
{
    ScalarField cc;
    FOREACH()
        cc(i,j) = f(i,j) > 0.5;

    #if LE_VOF
    if (istep % 2 == 0) {
        le_sweep_x(f, cc, dt);
        le_sweep_y(f, cc, dt);
    }
    else {
        //le_sweep_y(f, cc, dt);
        le_sweep_x(f, cc, dt);
    }
    #elif BASILISK_VOF
    if (istep % 2 == 0) {
        sweep_x(f, cc, dt);
        sweep_y(f, cc, dt);
    }
    else {
        sweep_y(f, cc, dt);
        sweep_x(f, cc, dt);
    }
    #else // EI_VOF
    if (istep % 2 == 0) {
        ei_sweep_x(f, cc, dt);
        ei_sweep_y(f, cc, dt);
    }
    else {
        ei_sweep_y(f, cc, dt);
        ei_sweep_x(f, cc, dt);
    }
    #endif

    update_boundary(f);

    FOREACH()
        f(i,j) = clamp(f(i,j), 0., 1.);
}
Event vof_advection_event (vof_advection, "vof advection", 0);
