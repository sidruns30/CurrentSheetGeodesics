#include "KSmetric.hpp"

namespace KS
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
        double _g01 = 2 * M * r / sigma(M, a, r, theta);
        return _g01;
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
        double _g11 = 1 + 2 * M * r / sigma(M, a, r, theta);
        return _g11;
    }

    double g12l(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g13l(double M, double a, double r, double theta)
    {
        double _g13 = - (1 + 2 * M * r / sigma(M, a, r, theta)) * a * SQR(sin(theta));
        return _g13;
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
        double _g00 = - (1 + 2 * M * r / sigma(M, a, r, theta));
        return _g00;
    }

    double g01u(double M, double a, double r, double theta)
    {
        double _g01 = 2 * M * r / sigma(M, a, r, theta);
        return _g01;
    }

    double g02u(double M, double a, double r, double theta)
    {
        return 0.;
    }

    double g03u(double M, double a, double r, double theta)
    {
        return 0.;
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
        double _g13 = a / sigma(M, a, r, theta);
        return _g13;
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
        double _g33 = 1 / (sigma(M, a, r, theta) * SQR(sin(theta)));
        return _g33;
    }
}

void GetKSCovariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double a = a_BH;
    double M = M_BH;
    g[0][0] = KS::g00l(M, a, r, theta);
    g[0][1] = KS::g01l(M, a, r, theta);
    g[0][2] = KS::g02l(M, a, r, theta);
    g[0][3] = KS::g03l(M, a, r, theta);
    g[1][0] = g[0][1];
    g[1][1] = KS::g11l(M, a, r, theta);
    g[1][2] = KS::g12l(M, a, r, theta);
    g[1][3] = KS::g13l(M, a, r, theta);
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    g[2][2] = KS::g22l(M, a, r, theta);
    g[2][3] = KS::g23l(M, a, r, theta);
    g[3][0] = g[0][3];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
    g[3][3] = KS::g33l(M, a, r, theta);
    return;
}

void GetKSContravariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double r = X[1];
    double theta = X[2];
    double a = a_BH;
    double M = M_BH;
    g[0][0] = KS::g00u(M, a, r, theta);
    g[0][1] = KS::g01u(M, a, r, theta);
    g[0][2] = KS::g02u(M, a, r, theta);
    g[0][3] = KS::g03u(M, a, r, theta);
    g[1][0] = g[0][1];
    g[1][1] = KS::g11u(M, a, r, theta);
    g[1][2] = KS::g12u(M, a, r, theta);
    g[1][3] = KS::g13u(M, a, r, theta);
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    g[2][2] = KS::g22u(M, a, r, theta);
    g[2][3] = KS::g23u(M, a, r, theta);
    g[3][0] = g[0][3];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
    g[3][3] = KS::g33u(M, a, r, theta);
    return;
}

// Kerr Schild connection coefficients
void GetKSConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM])
{
  double r1 = X[1];
    double r2 = r1*r1;
    double r3 = r2*r1;
    double r4 = r3*r1;

    double th = X[2];
    double dthdx2 = 1.0;
    double d2thdx22 = 0.0;
    double dthdx22 = dthdx2*dthdx2;
    double sth=sin(th);
    double cth=cos(th);
    double sth2 = sth*sth;
    double r1sth2 = r1*sth2;
    double sth4 = sth2*sth2;
    double cth2 = cth*cth;
    double cth4 = cth2*cth2;
    double s2th = 2.*sth*cth;
    double c2th = 2.*cth2 - 1.;

    double a = a_BH;
    double a2 = a*a;
    double a2sth2 = a2*sth2;
    double a2cth2 = a2*cth2;
    double a3 = a2*a;
    double a4 = a3*a;

    double rho2 = r2 + a2cth2;
    double rho22 = rho2*rho2;
    double rho23 = rho22*rho2;
    double irho2 = 1./rho2;
    double irho22 = irho2*irho2;
    double irho23 = irho22*irho2;
    double irho23_dthdx2 = irho23/dthdx2;

    double fac1 = r2 - a2cth2;
    double fac1_rho23 = fac1*irho23;
    double fac2 = a2 + 2.*r2 + a2*c2th;
    double fac3 = a2 + r1*(-2. + r1);

    double fac4 = r2 + a2 * cth2;
    double fac5 = r1 * (r1 + 2.);

    gamma[0][0][0] = 2.*r1*fac1_rho23;
    gamma[0][0][1] = (r2 - a2 * cth2) * (r1 * (2. + r1) + a2 * cth2) / pow(r2 + a2 * cth2, 3.); //
    gamma[0][0][2] = -a2*r1*s2th*dthdx2*irho22;
    gamma[0][0][3] = -2.*a*r1sth2*fac1_rho23;

    gamma[0][1][0] = gamma[0][0][1];
    gamma[0][1][1] = 2. * (r2 - a2 * cth2) * (r1 + r2 + a2 * cth2) / pow(r2 + a2 * cth2, 3.); //
    gamma[0][1][2] = -(2. * a2 * r1 * cth * sth) / pow(r2 + a2 * cth2, 2.); //
    gamma[0][1][3] = a * (-r2 + a2 * cth2) * (r1 * (2. + r1) + a2 * cth2) *sth2 / pow(r2 + a2 * cth2, 3.); //

    gamma[0][2][0] = gamma[0][0][2];
    gamma[0][2][1] = gamma[0][1][2];
    gamma[0][2][2] = -2.*r2*dthdx22*irho2;
    gamma[0][2][3] = a3*r1sth2*s2th*dthdx2*irho22;

    gamma[0][3][0] = gamma[0][0][3];
    gamma[0][3][1] = gamma[0][1][3];
    gamma[0][3][2] = gamma[0][2][3];
    gamma[0][3][3] = 2.*r1sth2*(-r1*rho22 + a2sth2*fac1)*irho23;

    gamma[1][0][0] = (r2 - a2 * cth2) * (0.5 * (a2 + r1 * (r1 - 2.)) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + r1 * (2. + r1)) * cos(2. * th)) - a2 * ((r1 - 2.) * r1 + a2 * cth2) * fac3 * sth2) / (pow(r2 + a2 * cth2, 3.) * (pow(a2 + r2, 2.) - a2 * (r1 * (2. + r1) + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
    gamma[1][0][1] = 2. * (-a2 + 4. * r1 + a2 * cos(2. * th)) * (a2 - 2. * r2 + a2 * cos(2. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //
    gamma[1][0][2] = 0.; //
    gamma[1][0][3] = a * (a2 - 2. * r2 + a2 * cos(2. * th)) * sth2 * ((a2 + (r1 - 2.) * r1) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + fac5) * cos(2. * th)) - a2 * (a2 + 2. * (r1 - 2.) * r1 + a2 * cos(2. * th)) * fac3 * sth2) /
                        (4. * pow(r2 + a2 * cth2, 3.) * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //

    gamma[1][1][0] = gamma[1][0][1];
    gamma[1][1][1] = 4. * (a2 - 2. * r2 + a2 * cos(2. * th)) * (fac5 + a2 * cos(2. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

    gamma[1][1][2] = -a2 * sin(2. * th) / (a2 + 2. * r2 + a2 * cos(2. * th)); //
    gamma[1][1][3] = a * (a4 - 8. * a2 * r1 + 3. * a4 * r1 - 4. * a2 * r2 + 16. * r3 + 8. * a2 * r3 + 8. * r4 * r1 + 4. * a2 * r1 * (-2. + a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) * sth2 / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

    gamma[1][2][0] = gamma[1][0][2];
    gamma[1][2][1] = gamma[1][1][2];
    gamma[1][2][2] = -(r1 * ((a2 + r1 * (r1 - 2.)) * (a4 + 2. * r4 + a2 * r1 * (3. * r1 - 2.) + a2 * (a2 + fac5) * cos(2. * th)) - 2. * a2 * (r1 * (r1 - 2.) + a2 * cth2) * fac3 * sth2)) / (2. * fac4 * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
    gamma[1][2][3] = 0.; //

    gamma[1][3][0] = gamma[1][0][3];
    gamma[1][3][1] = gamma[1][1][3];
    gamma[1][3][2] = gamma[1][2][3];
    gamma[1][3][3] = -(a2 + r1 * (r1 - 2.)) * (8. * r4 * r1 + 4. * a2 * r2 * (2. * r1 - 1.) + a4 * (1. + 3. * r1) + 4. * a2 * r1 * (a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) * sth2 / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

    gamma[2][0][0] = -a2*r1*s2th*irho23_dthdx2;
    gamma[2][0][1] = -2. * a2 * r1 * cth * sth / pow(r2 + a2 * cth2, 3.); //
    gamma[2][0][2] = 0.0;
    gamma[2][0][3] = a*r1*(a2+r2)*s2th*irho23_dthdx2;

    gamma[2][1][0] = gamma[2][0][1];
    gamma[2][1][1] = -2. * a2 * r1 * cth * sth / pow(r2 + a2 * cth2, 3.); //
    gamma[2][1][2] = r1 / (r2 + a2 * cth2); //
    gamma[2][1][3] = a * cth * sth * (r3 * (2. + r1) + 2. * a2 * r1 * (1. + r1) * cth2 + a4 * cth4 + 2. * a2 * r1 * sth2) / pow(r2 + a2 * cth2, 3.); //

    gamma[2][2][0] = gamma[2][0][2];
    gamma[2][2][1] = gamma[2][1][2];
    gamma[2][2][2] = -a2*cth*sth*dthdx2*irho2 + d2thdx22/dthdx2;
    gamma[2][2][3] = 0.0;

    gamma[2][3][0] = gamma[2][0][3];
    gamma[2][3][1] = gamma[2][1][3];
    gamma[2][3][2] = gamma[2][2][3];
    gamma[2][3][3] = -cth * sth * (rho23 + a2sth2 * rho2 * (r1 * (4. + r1) +
                                                            a2cth2) + 2. * r1 * a4 * sth4) * irho23_dthdx2;

    gamma[3][0][0] = a * fac1_rho23;
    gamma[3][0][1] = (a * r2 - a3 * cth2) / ((r2 + a2 * cth2) * (pow(a2 + r2, 2.) - a2 * (fac5 + a2 * cth2) * sth2 - a2 * fac3 * sth2)); //
    gamma[3][0][2] = -2. * a * r1 * cth * dthdx2 / (sth * rho22);
    gamma[3][0][3] = -a2sth2 * fac1_rho23;

    gamma[3][1][0] = gamma[3][0][1];
    gamma[3][1][1] = 8. * (a * r2 - a3 * cth2) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //
    gamma[3][1][2] = (a * (fac5 + a2 * cth2) * (1. / tan(th))) / (-pow(a2 + r2, 2.) + a2 * (fac5 + a2 * cth2) * sth2 + a2 * fac3 * sth2); //
    gamma[3][1][3] = (8. * r4 * r1 + 4. * a2 * r2 * (2. * r1 - 1.) + a4 * (1. + 3. * r1) + 4. * a2 * r1 * (a2 + r1 + 2. * r2) * cos(2. * th) + a4 * (r1 - 1.) * cos(4. * th)) / pow(a2 + 2. * r2 + a2 * cos(2. * th), 3.); //

    gamma[3][2][0] = gamma[3][0][2];
    gamma[3][2][1] = gamma[3][1][2];
    gamma[3][2][2] = -a * r1 * dthdx22 * irho2;
    gamma[3][2][3] = dthdx2 * (0.25 * fac2 * fac2 * cth / sth + a2 * r1 *
                                s2th) * irho22;

    gamma[3][3][0] = gamma[3][0][3];
    gamma[3][3][1] = gamma[3][1][3];
    gamma[3][3][2] = gamma[3][2][3];
    gamma[3][3][3] = (-a*r1sth2*rho22 + a3*sth4*fac1)*irho23;

    return;
}
/****** Kerr Schild to Boyer Lindquist ******/
// Helpful functions taken from eq 55 of https://arxiv.org/pdf/0706.0622.pdf

double t_func(double r)
{
    return 2 * M_BH * r / (SQR(r) - 2 * M_BH * r + SQR(a_BH));
}

double phi_func(double r)
{
    return  a_BH / (SQR(r) - 2 * M_BH * r + SQR(a_BH));
}


namespace KS_int
{
    double rmin = 1.e-2;
    double dr = 0.0001;
    double dr2 = dr/2;
    double sum = 0.;
    double _r = rmin;
}

double t_integral(double r)
{
    double dr = KS_int::dr;
    double dr2 = KS_int::dr2;
    double sum = KS_int::sum;
    double rmin = KS_int::rmin;
    double _r = KS_int::_r;

    while(_r <= r)
    {
        sum += 0.5 * (t_func(_r + dr2) + t_func(_r - dr2)) * dr;
        _r += dr;
    }

    return sum;
}

double phi_integral(double r)
{
    double dr = KS_int::dr;
    double dr2 = KS_int::dr2;
    double sum = KS_int::sum;
    double rmin = KS_int::rmin;
    double _r = KS_int::_r;
    while(_r <= r)
    {
        sum += 0.5 * (phi_func(_r + dr2) + phi_func(_r - dr2)) * dr;
        _r += dr;
    }

    return sum;
}

// Coordinate transformations - don't transform like the 4 vector
void X_KSToBL(double XKS[NDIM], double XBL[NDIM])
{
    XBL[0] = XKS[0] - t_integral(XKS[1]);
    XBL[1] = XKS[1];
    XBL[2] = XKS[2];
    XBL[3] = XKS[3] - phi_integral(XKS[1]);
    return;
}

void X_BLToKS(double XBL[NDIM], double XKS[NDIM])
{
    XKS[0] = XBL[0] + t_integral(XBL[1]);
    XKS[1] = XBL[1];
    XKS[2] = XBL[2];
    XKS[3] = XBL[3] + phi_integral(XBL[1]);
    return;
}

// Now transform an arbitrary tensor from KS to BL
void T_KSToBL(double TKS[NDIM], double TBL[NDIM], double XKS[NDIM])
{
    TBL[0] = TKS[0] - t_func(XKS[1]) * TKS[1];
    TBL[1] = TKS[1];
    TBL[2] = TKS[2];
    TBL[3] = TKS[3] - phi_func(XKS[1]) * TKS[1];
    return;
}

void T_BLToKS(double TBL[NDIM], double TKS[NDIM], double XBL[NDIM])
{
    TKS[0] = TBL[0] + t_func(XBL[1]) * TBL[1];
    TKS[1] = TBL[1];
    TKS[2] = TBL[2];
    TKS[3] = TBL[3] + phi_func(XBL[1]) * TBL[1];
    return;
}

// Convert from Cartesianized spatial coordinates to Kerr Schild spatial coordinates
void CartToKS(double x, double y, double z, double XKS[NDIM])
{
    // Time coordinate is meaningless since we do not have the information
    XKS[0] = 0.;
    XKS[1] = sqrt(SQR(x) + SQR(y) + SQR(z));
    XKS[2] = acos(z / XKS[1]);
    XKS[3] = atan2(y, x);
    return;
}

// Get the vector fields in Kerr Schild from cartesianized components (taken from regrid.py)
void GetKSField(double x, double y, double z, double b1, double b2, double b3, double &br, 
                double &bth, double &bph)
{
    double r = sqrt(SQR(x) + SQR(y) + SQR(z));
    double a = sqrt(SQR(x) + SQR(y));

    br = (x/r)*b1 + (y/r)*b2 + (z/r)*b3;
    bth = (x*z)/(SQR(r)*a)*b1 + (y*z)/(SQR(r)*a)*b2 - a/SQR(r) * b3;
    bph = (-y/SQR(a))*b1 + (x/SQR(a))*b2;
    return;
}

// Convert a tensor in Cartesian coordinates to KS coordinates (basically cartesian to spherical)
void T_3CartTo3_KS(double T3Cart[NDIM-1], double T3KS[NDIM-1], double XKS[NDIM])
{
    double J[NDIM-1][NDIM-1];
    double r = XKS[1];
    double theta = XKS[2];
    double phi = XKS[3];

    J[0][0] = cos(phi) * sin(theta);
    J[1][0] = sin(phi) * sin(theta);
    J[2][0] = cos(theta);

    J[0][1] = r * cos(phi) * cos(theta);
    J[1][1] = r * sin(phi) * cos(theta);
    J[2][1] = -r * sin(theta);
    
    J[0][2] = -r * sin(phi) * sin(theta);
    J[1][2] = r * cos(phi) * sin(theta);
    J[2][2] = 0;

    // Now invert J
    double Jinv[NDIM-1][NDIM-1];
    bool success = Invert3Matrix(J, Jinv);

    if (!success)
    {
        std::cout<<"Failed to invert matrix \n";
        return;
    }

    T3KS[0] = Jinv[0][0] * T3Cart[0] + Jinv[0][1] * T3Cart[1] + Jinv[0][2] * T3Cart[2];
    T3KS[1] = Jinv[1][0] * T3Cart[0] + Jinv[1][1] * T3Cart[1] + Jinv[1][2] * T3Cart[2];
    T3KS[2] = Jinv[2][0] * T3Cart[0] + Jinv[2][1] * T3Cart[1] + Jinv[2][2] * T3Cart[2];

    return;
}


// Convert a tensor in KS coordinates to Cartesian coordinates (basically spherical to cartesian) for HARM
void T_3KSTo_3Cart(double T3KS[NDIM-1], double T3Cart[NDIM-1], double XKS[NDIM])
{
    double J[NDIM-1][NDIM-1];
    double r = XKS[1];
    double theta = XKS[2];
    double phi = XKS[3];

    J[0][0] = cos(phi) * sin(theta);
    J[1][0] = sin(phi) * sin(theta);
    J[2][0] = cos(theta);

    J[0][1] = r * cos(phi) * cos(theta);
    J[1][1] = r * sin(phi) * cos(theta);
    J[2][1] = -r * sin(theta);
    
    J[0][2] = -r * sin(phi) * sin(theta);
    J[1][2] = r * cos(phi) * sin(theta);
    J[2][2] = 0;


    T3Cart[0] = J[0][0] * T3KS[0] + J[0][1] * T3KS[1] + J[0][2] * T3KS[2];
    T3Cart[1] = J[1][0] * T3KS[0] + J[1][1] * T3KS[1] + J[1][2] * T3KS[2];
    T3Cart[2] = J[2][0] * T3KS[0] + J[2][1] * T3KS[1] + J[2][2] * T3KS[2];

    return;
}

// Convert a tensor in log KS coordinates to regular KS
void T_3logKSTo3_KS(double T3logKS[NDIM-1], double T3KS[NDIM-1], double XKS[NDIM])
{
    double r = XKS[1];
    T3KS[0] = r * T3logKS[0];
    T3KS[1] = T3logKS[1];
    T3KS[2] = T3logKS[2];
    return;
}


// Convert a tensor in regular KS coordinates to log KS
void T_3KSTo3_logKS(double T3KS[NDIM-1], double T3logKS[NDIM-1], double XKS[NDIM])
{
    double r = XKS[1];
    T3logKS[0] = T3KS[0] / r;
    T3logKS[1] = T3KS[1];
    T3logKS[2] = T3KS[2];
    return;
}