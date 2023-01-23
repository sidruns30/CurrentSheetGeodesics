/*
    All the problem defitions and imports go in here!
*/
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <bits/stdc++.h>
#include <omp.h>
#include "cnpy/cnpy.h"

// OpenMP support
#if defined(_OPENMP)
    const int maxthreads = omp_get_max_threads();
#else 
    const int maxthreads = 1;
#endif

#define NDIM (4)

// Definitions appropriate for the metric
#define M_BH (1)
#define a_BH (0.9375)
#define SQR(x) ((x)*(x))
#define PI (3.141592653589793238)

// Definitions that may be useful for the MKS coordinatees
#define h_ks (0.9)
#define R0_ks (0.)

// Definitions for the geodesics
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


// Variables useful to find the current sheet
const double SIGMA_THRESHOLD = 3.;
const double BOX_THRESHOLD = 0.01;
const double BSQR_THRESHOLD = 10.;

// Necessary for field calculations
const bool idealMHD = false;

// Define types
typedef std::vector <double> ARRAY;
typedef std::vector <ARRAY> ARRAY2D;
typedef std::vector <ARRAY2D> ARRAY3D;

/*
extern ARRAY Bsqr_sim,bfluid0_sim,bfluid1_sim,bfluid2_sim,bfluid3_sim;
extern ARRAY efluid0_sim,efluid1_sim,efluid2_sim,efluid3_sim;
extern ARRAY b2_sim, e2_sim, rMKS_sim, thetaMKS_sim;;
extern ARRAY2D COORDS;
extern ARRAY2D PRIMS;
*/

extern ARRAY3D COORDS_BLOCKS;
extern ARRAY3D PRIMS_BLOCKS;


extern std::unordered_map <std::string, int> iprim;

struct photon
{
    double X[NDIM];
    double K[NDIM];
    double dKdlam[NDIM];
    double E0;
};

typedef struct photon photon;
extern double startx[4], stopx[4];

bool Invert4Matrix(const double m[16], double invOut[16]);
bool Invert3Matrix(const double m[3][3], double minv[3][3]);
void Multiply4Matrices(double a[NDIM][NDIM], double b[NDIM][NDIM], double c[NDIM][NDIM]);
void Multiply3Matrices(double a[NDIM-1][NDIM-1], double b[NDIM-1][NDIM-1], double c[NDIM-1][NDIM-1]);
void TransposeMatrix(double m[NDIM][NDIM], double minv[NDIM][NDIM]);
std::vector <size_t> sort_indices(const std::vector <double> &v);


void WriteVectorToFile(std::string fname, std::vector <std::vector <double>> &data);
void WriteVectorToNumpyArray(std::string fname, ARRAY &data);

// For debugging
template <typename T>
void printvar(std::string out, T var)
{
    std::cout<<out<<": "<<var<<"\n";
}

void print(std::string out);
void print3(std::string name, double vec[NDIM]);
void print4(std::string name, double vec[NDIM]);
void print3M(std::string name, double mat[NDIM-1][NDIM-1]);
void print4M(std::string name, double mat[NDIM][NDIM]);