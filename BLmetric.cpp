#include "BLmetric.hpp"

namespace BL
{
    double delta(double M, double a, double r, double theta)
    {
        double _delta = SQR(r) - 2 * M * r + SQR(a);
        return _delta;
    }

    double sigma(double M, double a, double r, double theta)
    {
        double _sigma = SQR(r) + SQR(a * cos(theta));
        return _sigma;
    }

    // Covariant components of g_ij: only defined for j >= i
    double g00l(double M, double a, double r, double theta)
    {
        double _g00 = - (1 - 2 * M * r / sigma(M, a, r, theta));
        return _g00;
    }

    double g01l(double M, double a, double r, double theta)
    {
        return 0;
    }

    double g02l(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g03l(double M, double a, double r, double theta)
    {
        double _g03 = - (2 * M * a * r / sigma(M, a, r, theta)) * SQR(sin(theta));
        return _g03;
    }

    double g11l(double M, double a, double r, double theta)
    {
        double _g11 = sigma(M, a, r, theta) / delta(M, a, r, theta);
        return _g11;
    }

    double g12l(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g13l(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g22l(double M, double a, double r, double theta)
    {
        double _g22 = sigma(M, a, r, theta);
        return _g22;
    }

    double g23l(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g33l(double M, double a, double r, double theta)
    {
        double _g33 = (SQR(r) + SQR(a) + (2 * M * r * SQR(a) / sigma(M, a, r, theta)) * SQR(sin(theta))) * SQR(sin(theta));
        return _g33;
    }

    // Contravariant componetns of g^ij with j>=i
    double g00u(double M, double a, double r, double theta)
    {
        double _g00 = -1. -2 * M * r * (SQR(r) + SQR(a)) / (delta(M, a, r, theta) * sigma(M, a, r, theta));
        return _g00;
    }

    double g01u(double M, double a, double r, double theta)
    {
        double _g01 = 0;
        return _g01;
    }

    double g02u(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g03u(double M, double a, double r, double theta)
    {
        double _g03 = - (2. * M * a * r) / (delta(M, a, r, theta) * sigma(M, a, r, theta));
        return _g03;
    }

    double g11u(double M, double a, double r, double theta)
    {
        double _g11 = delta(M, a, r, theta) / sigma(M, a, r, theta);
        return _g11;
    }

    double g12u(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g13u(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g22u(double M, double a, double r, double theta)
    {
        double _g22 = 1 / sigma(M, a, r, theta);
        return _g22;
    }

    double g23u(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g33u(double M, double a, double r, double theta)
    {
        double _g33 = (sigma(M, a, r, theta) - 2 * M * r) / (delta(M, a, r, theta) * sigma(M, a, r, theta) * SQR(sin(theta)));
        return _g33;
    }
}

void GetBLCovariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double a = a_BH;
    double M = M_BH;
    g[0][0] = BL::g00l(M, a, r, theta);
    g[0][1] = BL::g01l(M, a, r, theta);
    g[0][2] = BL::g02l(M, a, r, theta);
    g[0][3] = BL::g03l(M, a, r, theta);
    g[1][0] = g[0][1];
    g[1][1] = BL::g11l(M, a, r, theta);
    g[1][2] = BL::g12l(M, a, r, theta);
    g[1][3] = BL::g13l(M, a, r, theta);
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    g[2][2] = BL::g22l(M, a, r, theta);
    g[2][3] = BL::g23l(M, a, r, theta);
    g[3][0] = g[0][3];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
    g[3][3] = BL::g33l(M, a, r, theta);
    return;
}
void GetBLContravariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double a = a_BH;
    double M = M_BH;
    g[0][0] = BL::g00u(M, a, r, theta);
    g[0][1] = BL::g01u(M, a, r, theta);
    g[0][2] = BL::g02u(M, a, r, theta);
    g[0][3] = BL::g03u(M, a, r, theta);
    g[1][0] = g[0][1];
    g[1][1] = BL::g11u(M, a, r, theta);
    g[1][2] = BL::g12u(M, a, r, theta);
    g[1][3] = BL::g13u(M, a, r, theta);
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    g[2][2] = BL::g22u(M, a, r, theta);
    g[2][3] = BL::g23u(M, a, r, theta);
    g[3][0] = g[0][3];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
    g[3][3] = BL::g33u(M, a, r, theta);
    return;
}

void GetBLConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM])
{
    double r       = X[1];
    double rfactor = 1;
    double theta   = X[2];
    double sint    = sin(theta);
    double cost    = cos(theta);
    double sigma   = r * r + a_BH* a_BH* cost * cost;
    double delta   = r * r + a_BH* a_BH- 2. * r;
    double A       = (r * r + a_BH* a_BH) * (r * r + a_BH* a_BH) - delta * a_BH* a_BH*
                    sint * sint;
    double sigma3  = sigma * sigma * sigma;

    // Unique, non-zero connection elements
    gamma[1][0][0] =  delta / sigma3 * (2. * r * r - sigma) / rfactor;
    gamma[2][0][0] = -2. * a_BH* a_BH* r * sint * cost / sigma3;
    gamma[1][1][1] =  (1. - r) / delta + r / sigma;
    gamma[2][1][1] =  a_BH* a_BH* sint * cost / (sigma * delta) * rfactor *
                        rfactor;
    gamma[1][2][2] = -r * delta / sigma / rfactor;
    gamma[2][2][2] = -a_BH * a_BH* sint * cost / sigma;
    gamma[1][3][3] = -delta * sint * sint / sigma * (r - a_BH* a_BH* sint *
                                                        sint / (sigma * sigma) * (2. * r * r - sigma)) /
                        rfactor;
    gamma[2][3][3] = -sint * cost / sigma3 * ((r * r + a_BH* a_BH) *
                                                A - sigma * delta * a_BH* a_BH* sint * sint);
    gamma[0][0][1] =  (r * r + a_BH* a_BH) / (sigma * sigma * delta) *
                        (2. * r * r - sigma) * rfactor;
    gamma[3][0][1] =  a_BH/ (sigma * sigma * delta) * (2. * r * r - sigma)
                        * rfactor;
    gamma[0][0][2] = -2. * a_BH* a_BH* r * sint * cost / (sigma * sigma);
    gamma[3][0][2] = -2. * a_BH* r * cost / (sigma * sigma * sint);
    gamma[1][0][3] = -a_BH * delta * sint * sint / sigma3 *
                        (2. * r * r - sigma) / rfactor;
    gamma[2][0][3] =  2. * a_BH* r * (r * r + a_BH* a_BH) * sint * cost /
                        sigma3;
    gamma[1][1][2] = -a_BH * a_BH* sint * cost / sigma;
    gamma[2][1][2] =  r / sigma * rfactor;
    gamma[0][1][3] = -a_BH * sint * sint / (sigma * delta) * (2. * r *
                                                            r / sigma * (r * r + a_BH* a_BH) + r * r - a_BH* a_BH) *
                        rfactor;
    gamma[3][1][3] =  (r / sigma - a_BH* a_BH* sint * sint / (sigma * delta) *
                        (r - 1. + 2. * r * r / sigma)) * rfactor;
    gamma[0][2][3] =  2. * a_BH* a_BH* a_BH* r * sint * sint * sint * cost /
                        (sigma * sigma);
    gamma[3][2][3] =  cost / sint * (1. + 2. * a_BH* a_BH* r * sint * sint /
                                        (sigma * sigma));

    // Take symmetries into account
    gamma[0][1][0] = gamma[0][0][1];
    gamma[3][1][0] = gamma[3][0][1];
    gamma[0][2][0] = gamma[0][0][2];
    gamma[3][2][0] = gamma[3][0][2];
    gamma[1][3][0] = gamma[1][0][3];
    gamma[2][3][0] = gamma[2][0][3];
    gamma[1][2][1] = gamma[1][1][2];
    gamma[2][2][1] = gamma[2][1][2];
    gamma[0][3][1] = gamma[0][1][3];
    gamma[3][3][1] = gamma[3][1][3];
    gamma[0][3][2] = gamma[0][2][3];
    gamma[3][3][2] = gamma[3][2][3];

        return;
}