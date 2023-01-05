#include "SCmetric.hpp"

namespace SC
{
    double g00l(double M, double r, double theta)
    {
        double _g00 =  -(1. - 2. * M / r);
        return _g00;
    }

    double g11l(double M, double r, double theta)
    {
        double _g11 = 1. / (1. - 2 * M / r);
        return _g11;
    }

    double g22l(double M, double r, double theta)
    {
        double _g22 = SQR(r);
        return _g22;
    }

    double g33l(double M, double r, double theta)
    {
        double _g33 = SQR(r * sin(theta));
        return _g33;
    }

    double g00u(double M, double r, double theta)
    {
        double _g00 =  -1. / (1. - 2. * M / r);
        return _g00;
    }

    double g11u(double M, double r, double theta)
    {
        double _g11 =  (1. - 2. * M / r);
        return _g11;
    }

    double g22u(double M, double r, double theta)
    {
        double _g22 = 1. / SQR(r);
        return _g22;
    }

    double g33u(double M, double r, double theta)
    {
        double _g33 = 1. / SQR(r * sin(theta));
        return _g33;
    }
}

void GetSCCovariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double M = M_BH;
    g[0][0] = SC::g00l(M, r, theta);
    g[0][1] = 0.;
    g[0][2] = 0.;
    g[0][3] = 0.;
    g[1][0] = 0.;
    g[1][1] = SC::g11l(M, r, theta);
    g[1][2] = 0.;
    g[1][3] = 0.;
    g[2][0] = 0.;
    g[2][1] = 0.;
    g[2][2] = SC::g22l(M, r, theta);
    g[2][3] = 0.;
    g[3][0] = 0.;
    g[3][1] = 0.;
    g[3][2] = 0.;
    g[3][3] = SC::g33l(M, r, theta);
    return;
}
void GetSCContravariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double M = M_BH;
    g[0][0] = SC::g00u(M, r, theta);
    g[0][1] = 0.;
    g[0][2] = 0.;
    g[0][3] = 0.;
    g[1][0] = 0.;
    g[1][1] = SC::g11u(M, r, theta);
    g[1][2] = 0.;
    g[1][3] = 0.;
    g[2][0] = 0.;
    g[2][1] = 0.;
    g[2][2] = SC::g22u(M, r, theta);
    g[2][3] = 0.;
    g[3][0] = 0.;
    g[3][1] = 0.;
    g[3][2] = 0.;
    g[3][3] = SC::g33u(M, r, theta);
    return;
}

void GetSCConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM])
{
    std::cout<<"WARNING: SC connection coefficients not supported yet"<<std::endl;
    std::cout<<"Call the numerical version";
    return;
}