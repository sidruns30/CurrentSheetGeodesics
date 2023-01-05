#include "geodesic.hpp"

// In KS coordinates
double startx[NDIM] = {0., 1., 0., 0.};
double stopx[NDIM] = {0., 10., 0., 0.};

double stepsize(double X[NDIM], double K[NDIM])
{
    double dl, dlx1,dlxx1, dlx2, dlx3;
    double idlx1,idlxx1, idlx2, idlx3;

    dlx1 = EPS  / (fabs(K[1]) + SMALL);
    dlxx1 = EPS *X[1]/ (fabs(K[1]) + SMALL);
    dlx2 = EPS * fmin(X[2], stopx[2] - X[2]) / (fabs(K[2]) + SMALL);
    dlx3 = EPS / (fabs(K[3]) + SMALL);
    idlx1 = 1. / (fabs(dlx1) + SMALL);
    idlxx1 = 1. / (fabs(dlxx1) + SMALL);
    idlx2 = 1. / (fabs(dlx2) + SMALL);
    idlx3 = 1. / (fabs(dlx3) + SMALL);
    if (X[1] > stopx[1] * 1.05)
    {
        dl = 1/ (idlxx1 + idlx2 + idlx3);
    }
    else
    {
        dl = 1. / (idlx1 + idlx2 + idlx3);
    }
    return (dl);
}

void init_dKdlam(int metric_type, double X[], double Kcon[], double dK[])
{
  int k;
  double lconn[NDIM][NDIM][NDIM];
  double lconn_v2[NDIM][NDIM][NDIM];


  GetConnectionAnalytic(metric_type, X, lconn);

  for (k = 0; k < 4; k++) {

    dK[k] =
        -2. * (Kcon[0] * (lconn[k][0][1] * Kcon[1] + lconn[k][0][2] * Kcon[2] +
                          lconn[k][0][3] * Kcon[3]) +
               Kcon[1] * (lconn[k][1][2] * Kcon[2] + lconn[k][1][3] * Kcon[3]) +
               lconn[k][2][3] * Kcon[2] * Kcon[3]);

    dK[k] -= (lconn[k][0][0] * Kcon[0] * Kcon[0] +
              lconn[k][1][1] * Kcon[1] * Kcon[1] +
              lconn[k][2][2] * Kcon[2] * Kcon[2] +
              lconn[k][3][3] * Kcon[3] * Kcon[3]);
  }

    return;
}

void push_photon(int metric_type, double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                 double dl, double E0, int n)
{
   //  double lconn[NDIM][NDIM][NDIM];
    double lconn[NDIM][NDIM][NDIM];
    double Kcont[NDIM], K[NDIM], dK;
    double Xcpy[NDIM], Kcpy[NDIM], dKcpy[NDIM];
    double Gcov[NDIM][NDIM], E1;
    double dl_2, err, errE;
    int i, k, iter;
    
    if (X[1] < startx[1])
        return;

    FAST_CPY(X, Xcpy);
    FAST_CPY(Kcon, Kcpy);
    FAST_CPY(dKcon, dKcpy);

    dl_2 = 0.5 * dl;
    
    /* Step the position and estimate new wave vector */
    for (i = 0; i < NDIM; i++) 
    {
        dK = dKcon[i] * dl_2;
        Kcon[i] += dK;
        K[i] = Kcon[i] + dK;
        X[i] += Kcon[i] * dl;
    }

    GetConnectionAnalytic(metric_type, X, lconn);

    /* We're in a coordinate basis so take advantage of symmetry in the connection
    */
    iter = 0;
    do 
    {
        iter++;
        FAST_CPY(K, Kcont);

        err = 0.;
        for (k = 0; k < 4; k++) 
        {
            dKcon[k] =
                -2. *
                (Kcont[0] * (lconn[k][0][1] * Kcont[1] + lconn[k][0][2] * Kcont[2] +
                            lconn[k][0][3] * Kcont[3]) +
                Kcont[1] * (lconn[k][1][2] * Kcont[2] + lconn[k][1][3] * Kcont[3]) +
                lconn[k][2][3] * Kcont[2] * Kcont[3]);

            dKcon[k] -= (lconn[k][0][0] * Kcont[0] * Kcont[0] +
                        lconn[k][1][1] * Kcont[1] * Kcont[1] +
                        lconn[k][2][2] * Kcont[2] * Kcont[2] +
                        lconn[k][3][3] * Kcont[3] * Kcont[3]);

            K[k] = Kcon[k] + dl_2 * dKcon[k];
            err += fabs((Kcont[k] - K[k]) / (K[k] + SMALL));
        }

    } while (err > ETOL && iter < MAX_ITER);

    FAST_CPY(K, Kcon);

    GcovFunc(metric_type, X, Gcov);
    
    E1 = -(Kcon[0] * Gcov[0][0] + Kcon[1] * Gcov[0][1] + Kcon[2] * Gcov[0][2] +
         Kcon[3] * Gcov[0][3]);

    errE = fabs((E1 - (E0)) / (E0));

    if (n < 7 && (errE > 1.e-4 || err > ETOL || std::isnan(err) || std::isinf(err))) 
    {
        FAST_CPY(Xcpy, X);
        FAST_CPY(Kcpy, Kcon);
        FAST_CPY(dKcpy, dKcon);
        push_photon(metric_type, X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
        push_photon(metric_type, X, Kcon, dKcon, 0.5 * dl, E0, n + 1);
    }

    /* done! */
    return;
}

// The main geodesic integrator - saves the orbit of the photon every 0.01 Rg
void GetGeodesic(ARRAY2D &gdsc, int metric_type, double X[NDIM], double K[NDIM])
{
    gdsc.clear();

    // Initalize the derivative
    double dKdlam[NDIM];
    init_dKdlam(metric_type, X, K, dKdlam);

    ARRAY temp;
    temp.push_back(X[0]);
    temp.push_back(X[1]);
    temp.push_back(X[2]);
    temp.push_back(X[3]);

    gdsc.push_back(temp);
    double Kcov[NDIM];

    int nsteps = 0;
    while (temp[1] < stopx[1])
    {
        double step = stepsize(X, K);
        UpperToLower(metric_type, X, K, Kcov);

        double E0 = Kcov[0];
        push_photon(metric_type, X, K, dKdlam, step, E0, 0);

        // return if the timestep is negative
        if (nsteps > 10000 || temp[1]<startx[1])
        {
            break;
        }

        temp.clear();
        temp.push_back(X[0]);
        temp.push_back(X[1]);
        temp.push_back(X[2]);
        temp.push_back(X[3]);

        if (fabs(gdsc.back()[1] - X[1]) > 0.05)
        {
            gdsc.push_back(temp);
        }

        nsteps ++;
    }
    return;
}

// Output Geodesic DataFile
void GetGeodesics(ARRAY3D &gdscs, const ARRAY2D &X_K)
{
    size_t i,j, N=X_K.size();
    gdscs.clear();

    // Now start integrating geodesics
    // _x and _k contain the starting position and momentum vector of each geodesic
    double _x[NDIM], _k[NDIM], _xks[NDIM], _kks[NDIM];
    for (i=0; i<N; i++)
    {
        ARRAY2D gdsc;
        _x[0] = X_K[i][0];
        _x[1] = X_K[i][1];
        _x[2] = X_K[i][2];
        _x[3] = X_K[i][3];        
        _k[0] = X_K[i][4];
        _k[1] = X_K[i][5];
        _k[2] = X_K[i][6];
        _k[3] = X_K[i][8];

        // Geodesic integration is currently done in Kerr Schild
        TransformCoordinates(3, 2, _x, _xks);
        TransformFourVectors(3, 2, _x, _k, _kks);

        // now integrate the geodesic
        GetGeodesic(gdsc, 2, _xks, _kks);
        gdscs.push_back(gdsc);
    }
    return;
}