#include "geodesic.hpp"

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
  /*GetConnectionNumerical(metric_type, X, lconn_v2);

  int a,b,c;
  double diff[NDIM];
  for (a=0;a<NDIM;a++)
  {
    for (b=0;b<NDIM;b++)
    {
        for (c=0;c<NDIM;c++)
        {
            diff[c] = lconn[a][b][c] - lconn_v2[a][b][c];
        }
        print4("numerical", lconn_v2[a][b]);
        print4("analytical", lconn[a][b]);
        print4("diff", diff);
    }
  }
*/

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
    
    // NOTE TO SELF: Ask Jordy about commenting out this section of the code
    //#if MKS
    if (X[1] < startx[1])
        return;
    //#endif

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
std::vector <std::vector <double>> GetGeodesic(int metric_type, double X[NDIM], double K[NDIM])
{
    // Initalize the derivative
    double dKdlam[NDIM];
    init_dKdlam(metric_type, X, K, dKdlam);

    std::vector <std::vector <double>> X_traj;
    std::vector <double> temp;
    temp.push_back(X[0]);
    temp.push_back(X[1]);
    temp.push_back(X[2]);
    temp.push_back(X[3]);

    X_traj.push_back(temp);
    printvar("r init", temp[1]);
    print4("K", K);
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

        if (fabs(X_traj.back()[1] - X[1]) > 0.05)
        {
            X_traj.push_back(temp);
        }

        nsteps ++;
    }
    std::cout<<"\n";
    std::cout<<"nsteps: "<<nsteps<<'\n';
    printvar("r", temp[1]);
    return X_traj;
}

// Output Geodesic DataFile
void GetGeodesics(std::string Infile, std::string Outfile)
{
    std::ifstream infile(Infile);
    std::string line;
    int i,j;
    double row_data[2*NDIM];
    std::vector <double> x[NDIM], k[NDIM];
    while (std::getline(infile, line))
    {
        std::istringstream row(line);
        std::string row_s;
        for (i=0; i<2*NDIM; i++)
        {
            row >> row_s;
            row_data[i] = std::stod(row_s);
        }
        x[0].push_back(row_data[0]);
        x[1].push_back(row_data[1]);
        x[2].push_back(row_data[2]);
        x[3].push_back(row_data[3]);
        k[0].push_back(row_data[4]);
        k[1].push_back(row_data[5]);
        k[2].push_back(row_data[6]);
        k[3].push_back(row_data[7]);
    }

    std::cout<<"Loaded " << x[0].size() << " geodesics" << std::endl;
    // Now start integrating geodesics
    std::vector <std::vector <double>> x_geodesic;
    // _x and _k contain the starting vector of each geodesic
    double _x[NDIM], _k[NDIM], _xks[NDIM], _kks[NDIM];
    
    std::ofstream file(Outfile);
    if (file.is_open())
    {
        for (i=0;i<x[0].size();i++)
        {
            std::cout<<"Geodesic number: "<<i<<"\n";
            _x[0] = x[0][i];
            _x[1] = x[1][i];
            _x[2] = x[2][i];
            _x[3] = x[3][i];        
            _k[0] = k[0][i];
            _k[1] = k[1][i];
            _k[2] = k[2][i];
            _k[3] = k[3][i];

            TransformCoordinates(3, 2, _x, _xks);
            TransformFourVectors(3, 2, _x, _k, _kks);
            //X_MKSToKS_v2(_x, _xks);
            //T_MKSToKS_v2(_k, _kks, _x);


            // now integrate the geodesic
            x_geodesic = GetGeodesic(2, _xks, _kks);
            
            
            // write to file
            for (j=0; j<x_geodesic.size(); j++)
            {
                    file << x_geodesic[j][0] << "\t";
            }
            file<<std::endl;
            for (j=0; j<x_geodesic.size(); j++)
            {
                    file << x_geodesic[j][1] << "\t";
            }
            file<<std::endl;
            for (j=0; j<x_geodesic.size(); j++)
            {
                    file << x_geodesic[j][2] << "\t";
            }
            file<<std::endl;
            for (j=0; j<x_geodesic.size(); j++)
            {
                    file << x_geodesic[j][3] << "\t";
            }
            file<<std::endl;
        }
            file.close();
        
    
    }
    return;
}