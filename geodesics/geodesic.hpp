#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif
#include "../fluid/getk.hpp"

#define ETOL 1.e-3
#define MAX_ITER 10
#define EPS (0.01)
#define SMALL (1.e-40)

#define FAST_CPY(in, out)                                                      \
  {                                                                            \
    out[0] = in[0];                                                            \
    out[1] = in[1];                                                            \
    out[2] = in[2];                                                            \
    out[3] = in[3];                                                            \
  }

double stepsize(double X[NDIM], double K[NDIM]);
void push_photon(int metric_type, double X[NDIM], double Kcon[NDIM], double dKcon[NDIM],
                 double dl, double E0, int n);
void init_dKdlam(int metric_type, double X[], double Kcon[], double dK[]);
std::vector <std::vector <double>> GetGeodesic(int metric_type, double X[NDIM], double K[NDIM]);
void GetGeodesics(std::string Infile, std::string Outfile);
extern double startx[NDIM];
extern double stopx[NDIM];