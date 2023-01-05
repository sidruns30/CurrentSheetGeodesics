#include "MKSmetric.hpp"

/******** Modified Kerr Schild to Kerr Schild *****************/
double _theta(double vartheta)
{
    return vartheta + 0.5 * h_ks * sin(2 * vartheta);
}

void GetMKSCovariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double XKS[NDIM];
    TransformCoordinates(3, 2, X, XKS);

    GetKSCovariantMetric(g, XKS);
    double J[NDIM][NDIM], temp[NDIM][NDIM];
    int i, j;
    for (i=0;i<NDIM;i++)
    {
        for (j=0;j<NDIM;j++)
        {
            J[i][j] = 0.;
            temp[i][j] = 0.;
        }
    }

    J[0][0] = 1.;
    J[1][1] = exp(X[1]);
    J[2][2] = 1. + 2*h_ks + 12*h_ks*(SQR(X[2]/PI) - X[2]/PI);
    J[3][3] = 1.;

    Multiply4Matrices(J, g, temp);
    Multiply4Matrices(temp, J, g);
    return;
}

void GetMKSContravariantMetric(double g[NDIM][NDIM], double X[NDIM])
{
    double XKS[NDIM];
    TransformCoordinates(3, 2, X, XKS);

    GetKSContravariantMetric(g, X);
    double J[NDIM][NDIM], temp[NDIM][NDIM];
    int i, j;
    for (i=0;i<NDIM;i++)
    {
        for (j=0;j<NDIM;j++)
        {
            J[i][j] = 0.;
            temp[i][j] = 0.;
        }
    }

    J[0][0] = 1.;
    J[1][1] = 1./exp(X[1]);
    J[2][2] = 1./(1. + 2*h_ks + 12*h_ks*(SQR(X[2]/PI) - X[2]/PI));
    J[3][3] = 1.;

    Multiply4Matrices(J, g, temp);
    Multiply4Matrices(temp, J, g);
    return;
}

// Reference https://arxiv.org/pdf/1611.09720.pdf
void X_MKSToKS(double XMKS[NDIM], double XKS[NDIM])
{
    XKS[0] = XMKS[0];
    XKS[1] = R0_ks + exp(XMKS[1]);
    XKS[2] = XMKS[2] + (2 * h_ks * XMKS[2] / SQR(PI)) * (PI - 2*XMKS[2]) * (PI - XMKS[2]);
    XKS[3] = XMKS[3];
    return;
}

// Second version of the coordinate trasnforms
void X_MKSToKS_v2(double XMKS[NDIM], double XKS[NDIM])
{
    XKS[0] = XMKS[0];
    XKS[1] = R0_ks + exp(XMKS[1]);
    XKS[2] = XMKS[2] + (h_ks/2.) * sin(2 * XMKS[2]);
    XKS[3] = XMKS[3];
    return;
}


// Taken from https://link.springer.com/content/pdf/10.1186/s40668-017-0020-2.pdf
void X_KSToMKS(double XKS[NDIM], double XMKS[NDIM])
{
    XMKS[0] = XKS[0];
    XMKS[1] = log(XKS[1] - R0_ks);
    // This is a mess...
    double theta = XKS[2];
    double R_theta = pow(h_ks*pow(-3*h_ks*(-108*h_ks*SQR(theta) + 108*PI*h_ks*theta + 
                    (h_ks - 4)*SQR(2*PI*h_ks + PI)), 0.5) + 9*(PI - 2*theta)*SQR(h_ks), (1./3));
    
    XMKS[2] = (1./12)*pow(PI,(2./3))*((-2. * pow(2., (1./3))*pow(3*PI, (2./3))*(h_ks-1))/(R_theta) - 
    (pow(2., (2./3)) * pow(3., (1./3))* R_theta)/h_ks + 6*pow(PI, (1./3)));
    XMKS[3] = XKS[3];
    return;
}

void X_KSToMKS_v2(double XKS[NDIM], double XMKS[NDIM])
{
    XMKS[0] = XKS[0];
    XMKS[1] = log(XKS[1] - R0_ks);
    XMKS[3] = XKS[3];
    // We can do a crappy root find
    // Want to get varphi such that theta = varphi + h_ks/2 sin(var_phi)
    // Limits of varphi: 0 - pi
    // Good thing is that the functions are roughly monotonic except for varphi=pi/2
    double low = 0., high=PI;
    double tol = 1.e-6;
    double mid = (low + high)/2;
    double theta = XKS[2], eval;
    eval = _theta(mid);
    int niter = 0;
    while(fabs(eval - theta) > tol)
    {
        if (eval < theta)
        {
            low = mid;
            mid = 0.5 * (low + high);
        }
        else if (eval > theta)
        {
            high = mid;
            mid = 0.5 * (low + high);
        }
        eval = _theta(mid);
        niter ++;
    }
    XMKS[2] = mid;
    return;
}

void T_MKSToKS(double TMKS[NDIM], double TKS[NDIM], double XMKS[NDIM])
{
    TKS[0] = TMKS[0];
    TKS[1] = TMKS[1] * exp(XMKS[1]);
    TKS[2] = TMKS[2] * (1 + 2 * h_ks + 12 * h_ks * (SQR(XMKS[2]/PI) - XMKS[2]/PI));
    TKS[3] = TMKS[3];
    return;
}


void T_MKSToKS_v2(double TMKS[NDIM], double TKS[NDIM], double XMKS[NDIM])
{
    TKS[0] = TMKS[0];
    TKS[1] = TMKS[1] * exp(XMKS[1]);
    TKS[2] = TMKS[2] * (1 + h_ks * cos(2 * XMKS[2]));
    TKS[3] = TMKS[3];
    return;
}


void T_KSToMKS(double TKS[NDIM], double TMKS[NDIM], double XKS[NDIM])
{
    double XMKS[NDIM];
    TransformCoordinates(2, 3, XKS, XMKS);
    TMKS[0] = TKS[0];
    TMKS[1] = TKS[1] / exp(XMKS[1]);
    TMKS[2] = TKS[2] / (1 + 2 * h_ks + 12 * h_ks * (SQR(XMKS[2]/PI) - XMKS[2]/PI));
    TMKS[3] = TKS[3];
    return;
}

void T_KSToMKS_v2(double TKS[NDIM], double TMKS[NDIM], double XKS[NDIM])
{
    double XMKS[NDIM];
    TransformCoordinates(2, 3, XKS, XMKS);
    TMKS[0] = TKS[0];
    TMKS[1] = TKS[1] / exp(XMKS[1]);
    TMKS[2] = TKS[2] / (1 + h_ks * cos(2 * XMKS[2]));
    TMKS[3] = TKS[3];
    return;
}


// Same as above but hard coded 
void GetMKSCovariantMetric_v2(double g[NDIM][NDIM], double X[NDIM])
{
    double XKS[NDIM];
    TransformCoordinates(3, 2, X, XKS);

    GetKSCovariantMetric(g, XKS);
    double ds_dsks = exp(X[1]);
    double dth_dthks = 1 + h_ks * cos(2 * X[2]);

    g[0][1] *= ds_dsks;
    g[0][2] *= dth_dthks;
    g[1][0] *= ds_dsks;
    g[1][1] *= SQR(ds_dsks);
    g[1][2] *= ds_dsks * dth_dthks;
    g[1][3] *= ds_dsks;
    g[2][0] *= dth_dthks;
    g[2][1] *= ds_dsks * dth_dthks;
    g[2][2] *= SQR(dth_dthks);
    g[2][3] *= dth_dthks;
    g[3][1] *= ds_dsks;
    g[3][2] *= dth_dthks;
    return;
}


// Same as above but hard coded 
void GetMKSContravariantMetric_v2(double g[NDIM][NDIM], double X[NDIM])
{
    double XKS[NDIM];
    TransformCoordinates(3, 2, X, XKS);
    
    GetKSContravariantMetric(g, XKS);
    double ds_dsks = 1/exp(X[1]);
    double dth_dthks = 1/(1 + h_ks * cos(2 * X[2]));

    g[0][1] *= ds_dsks;
    g[0][2] *= dth_dthks;
    g[1][0] *= ds_dsks;
    g[1][1] *= SQR(ds_dsks);
    g[1][2] *= ds_dsks * dth_dthks;
    g[1][3] *= ds_dsks;
    g[2][0] *= dth_dthks;
    g[2][1] *= ds_dsks * dth_dthks;
    g[2][2] *= SQR(dth_dthks);
    g[2][3] *= dth_dthks;
    g[3][1] *= ds_dsks;
    g[3][2] *= dth_dthks;
    return;
}

void GetMKSConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM])
{
    double r           = R0_ks + exp(X[1]);
    double r2          = r * r;
    double rprime      = exp(X[1]);
    double rprime2     = rprime * rprime;
    double rprimeprime = rprime;
    double theta       = M_PI * X[2] + 0.5 * (1. - h_ks) * sin(2. * M_PI * X[2]);
    double thetaprime  = M_PI * (1. + (1. - h_ks) * cos(2. * M_PI * X[2]));
    double thetaprime2 = thetaprime * thetaprime;
    double thetaprimeprime = -2. * M_PI * M_PI * (1. - h_ks) * sin(2. * M_PI * X[2]);
    double costh       = cos(theta);
    double cos2th      = costh * costh;
    double sinth       = sin(theta);

    double sin2th      = sinth * sinth;
    double sintwoth    = sin(2. * theta);
    double cotth       = 1. / tan(theta);
    double a           = a_BH;
    double a2          = a * a;

    double Sigma       = r2 + a2 * cos2th;
    double Sigma2      = Sigma * Sigma;
    double Delta       = r2 - 2. * r + a2;
    double A           = Sigma * Delta + 2. * r * (r2 + a2);
    double Sigmabar    = 1. / Sigma;
    double Sigmabar2   = Sigmabar * Sigmabar;
    double Sigmabar3   = Sigmabar * Sigmabar * Sigmabar;
    double B           = (2. * r2 - Sigma) * Sigmabar3;
    double C           = r * Sigmabar - a2 * B * sin2th;

    // Gamma[t][mu][nu]
    gamma[0][0][0] = 2. * r * B;
    gamma[0][0][1] = (Sigma + 2. * r) * B * rprime;
    gamma[0][1][0] = gamma[0][0][1];
    gamma[0][0][2] = -a2 * r * sintwoth * Sigmabar2 * thetaprime;
    gamma[0][2][0] = gamma[0][0][2];
    gamma[0][0][3] = -2. * a * r * B * sin2th;
    gamma[0][3][0] = gamma[0][0][3];
    gamma[0][1][1] = 2. * (Sigma + r) * B * rprime2;
    gamma[0][1][2] = gamma[0][0][2] * rprime;
    gamma[0][2][1] = gamma[0][1][2];
    gamma[0][1][3] = -a * sin2th * gamma[0][0][1];
    gamma[0][3][1] = gamma[0][1][3];
    gamma[0][2][2] = -2. * r2 * Sigmabar * thetaprime2;
    gamma[0][2][3] = -a * sin2th * gamma[0][0][2];
    gamma[0][3][2] = gamma[0][2][3];
    gamma[0][3][3] = -2. * r * C * sin2th;

    // Gamma[r][mu][nu]
    gamma[1][0][0] = Delta * B / rprime;
    gamma[1][0][1] = (Delta - Sigma) * B;
    gamma[1][1][0] = gamma[1][0][1];
    gamma[1][2][0] = 0;
    gamma[1][0][2] = gamma[1][2][0];
    gamma[1][0][3] = -a * Delta * B * sin2th / rprime;
    gamma[1][3][0] = gamma[1][0][3];
    gamma[1][1][1] = rprimeprime / rprime - (2. * Sigma - Delta) * B * rprime;
    gamma[1][1][2] = -a2 * sinth * costh * Sigmabar * thetaprime;
    gamma[1][3][2] = 0;
    gamma[1][2][1] = gamma[1][1][2];
    gamma[1][2][3] = gamma[1][3][2];
    gamma[1][1][3] = a * (r * Sigmabar + (Sigma - Delta) * B) * sin2th;
    gamma[1][3][1] = gamma[1][1][3];
    gamma[1][2][2] = -r * Delta * Sigmabar /rprime * thetaprime2;
    gamma[1][3][3] = -Delta * C * sin2th / rprime;

    // Gamma[theta][mu][nu]
    gamma[2][0][0] = Sigmabar * gamma[0][0][2] / thetaprime2;
    gamma[2][0][1] = gamma[2][0][0] * rprime;
    gamma[2][1][0] = gamma[2][0][1];
    gamma[2][2][0] = 0;
    gamma[2][0][2] = gamma[2][2][0];
    gamma[2][0][3] = a * r * (r2 + a2) * sintwoth * Sigmabar3 / thetaprime;
    gamma[2][3][0] = gamma[2][0][3];
    gamma[2][1][1] = gamma[2][0][0] * rprime2;
    gamma[2][1][2] = r * Sigmabar * rprime;
    gamma[2][2][1] = gamma[2][1][2];
    gamma[2][1][3] = a * sinth * costh * (A + Sigma * (Sigma - Delta)) * Sigmabar3 * rprime / thetaprime;
    gamma[2][3][1] = gamma[2][1][3];
    gamma[2][3][2] = 0;
    gamma[2][2][3] = gamma[2][3][2];
    gamma[2][2][2] = thetaprimeprime / thetaprime + gamma[1][1][2];
    gamma[2][3][3] = -sinth * costh * (Delta * Sigma2 + 2. * r * (r2 + a2) * (r2 + a2)) * Sigmabar3 / thetaprime;

    // Gamma[phi][mu][nu]
    gamma[3][0][0] = a * B;
    gamma[3][0][1] = gamma[3][0][0] * rprime;
    gamma[3][1][0] = gamma[3][0][1];
    gamma[3][0][2] = -2. * a * r * cotth * Sigmabar2 * thetaprime;
    gamma[3][2][0] = gamma[3][0][2];
    gamma[3][0][3] = -a2 * B * sin2th;
    gamma[3][3][0] = gamma[3][0][3];
    gamma[3][1][1] = gamma[3][0][0] * rprime2;
    gamma[3][1][2] = -a * (Sigma + 2. * r) * cotth * Sigmabar2 * rprime * thetaprime;
    gamma[3][2][1] = gamma[3][1][2];
    gamma[3][1][3] = C * rprime;
    gamma[3][3][1] = gamma[3][1][3];
    gamma[3][2][2] = -a * r * Sigmabar * thetaprime2;
    gamma[3][2][3] = (cotth + a2 * r * sintwoth * Sigmabar2) * thetaprime;
    gamma[3][3][2] = gamma[3][2][3];
    gamma[3][3][3] = -a * C * sin2th;
    return;   
}

// Convert 4 vector in MKS to 4 vector in Cart
void T_3MKSTo3Cart(double T3MKS[NDIM-1], double T3Cart[NDIM-1], double XMKS[NDIM])
{
    double J[NDIM-1][NDIM-1];
    double es = exp(XMKS[1]);
    double vtheta = XMKS[2];
    double phi = XMKS[3];
    
    // This is basically theta KS
    double alpha = vtheta + 0.5 * h_ks * sin(2 * vtheta);
    // Similarly this is R_ks
    double r = R0_ks + es;

    J[0][0] = es * cos(phi) * sin(alpha);
    J[1][0] = es * sin(phi) * sin(alpha);
    J[2][0] = es * cos(alpha);

    J[0][1] = r * cos(phi) * cos(alpha) * (1. + h_ks * cos(2 * vtheta));
    J[1][1] = r * sin(phi) * cos(alpha) * (1. + h_ks * cos(2 * vtheta));
    J[2][1] = -r * sin(alpha) * (1 + h_ks * cos(2 * vtheta));
    
    J[0][2] = -r * sin(phi) * sin(alpha);
    J[1][2] = r * cos(phi) * sin(alpha);
    J[2][2] = 0;

    T3Cart[0] = J[0][0] * T3MKS[0] + J[0][1] * T3MKS[1] + J[0][2] * T3MKS[2];
    T3Cart[1] = J[1][0] * T3MKS[0] + J[1][1] * T3MKS[1] + J[1][2] * T3MKS[2];
    T3Cart[2] = J[2][0] * T3MKS[0] + J[2][1] * T3MKS[1] + J[2][2] * T3MKS[2];

    return;
}

void T_3CartTo3MKS(double T3Cart[NDIM-1], double T3MKS[NDIM-1], double XMKS[NDIM])
{
    double J[NDIM-1][NDIM-1];
    double es = exp(XMKS[1]);
    double vtheta = XMKS[2];
    double phi = XMKS[3];
    
    // This is basically theta KS
    double alpha = vtheta + 0.5 * h_ks * sin(2 * vtheta);
    // Similarly this is R_ks
    double r = R0_ks + es;

    J[0][0] = es * cos(phi) * sin(alpha);
    J[1][0] = es * sin(phi) * sin(alpha);
    J[2][0] = es * cos(alpha);

    J[0][1] = r * cos(phi) * cos(alpha) * (1. + h_ks * cos(2 * vtheta));
    J[1][1] = r * sin(phi) * cos(alpha) * (1. + h_ks * cos(2 * vtheta));
    J[2][1] = -r * sin(alpha) * (1 + h_ks * cos(2 * vtheta));
    
    J[0][2] = -r * sin(phi) * sin(alpha);
    J[1][2] = r * cos(phi) * sin(alpha);
    J[2][2] = 0;

    // Now invert J
    double Jinv[NDIM-1][NDIM-1];
    bool success = Invert3Matrix(J, Jinv);

    if (!success)
    {
        std::cout<<"Failed to invert matrix \n";
        return;
    }

    T3MKS[0] = Jinv[0][0] * T3Cart[0] + Jinv[0][1] * T3Cart[1] + Jinv[0][2] * T3Cart[2];
    T3MKS[1] = Jinv[1][0] * T3Cart[0] + Jinv[1][1] * T3Cart[1] + Jinv[1][2] * T3Cart[2];
    T3MKS[2] = Jinv[2][0] * T3Cart[0] + Jinv[2][1] * T3Cart[1] + Jinv[2][2] * T3Cart[2];

    return;
}

// Try to invert using the second way (non transcedental way)
void T_3CartTo3MKS_v2(double T3Cart[NDIM-1], double T3MKS[NDIM-1], double XMKS[NDIM])
{
    double J[NDIM-1][NDIM-1];
    double es = exp(XMKS[1]);
    double vtheta = XMKS[2];
    double phi = XMKS[3];
    
    // This is basically theta KS
    double alpha = vtheta + 0.5 * h_ks * sin(2 * vtheta);
    double dvthdth = (1 + 2 * h_ks + 12 * h_ks * (SQR(XMKS[2]/PI) - XMKS[2]/PI));

    // Similarly this is R_ks
    double r = R0_ks + es;

    J[0][0] = es * cos(phi) * sin(alpha);
    J[1][0] = es * sin(phi) * sin(alpha);
    J[2][0] = es * cos(alpha);

    J[0][1] = r * cos(phi) * cos(alpha) * dvthdth;
    J[1][1] = r * sin(phi) * cos(alpha) * dvthdth;
    J[2][1] = -r * sin(alpha) * dvthdth;
    
    J[0][2] = -r * sin(phi) * sin(alpha);
    J[1][2] = r * cos(phi) * sin(alpha);
    J[2][2] = 0;

    // Now invert J
    double Jinv[NDIM-1][NDIM-1];
    bool success = Invert3Matrix(J, Jinv);

    if (!success)
    {
        std::cout<<"Failed to invert matrix \n";
        return;
    }

    T3MKS[0] = Jinv[0][0] * T3Cart[0] + Jinv[0][1] * T3Cart[1] + Jinv[0][2] * T3Cart[2];
    T3MKS[1] = Jinv[1][0] * T3Cart[0] + Jinv[1][1] * T3Cart[1] + Jinv[1][2] * T3Cart[2];
    T3MKS[2] = Jinv[2][0] * T3Cart[0] + Jinv[2][1] * T3Cart[1] + Jinv[2][2] * T3Cart[2];

    return;
}