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



#define NDIM (4)

// Definitions appropriate for the metric
#define M_BH (1)
#define a_BH (0.9375)
#define SQR(x) ((x)*(x))
#define PI (3.141592653589793238)

// Definitions that may be useful for the MKS coordinatees
#define h_ks (0.9)
#define R0_ks (0.)

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


//extern void WriteVectorToFile(std::string fname, std::vector <std::vector <double>> data);

void WriteVectorToFile(std::string fname, std::vector <std::vector <double>> data);
void print(std::string out);
void printvar(std::string out, double var);
void print3(std::string name, double vec[NDIM]);
void print4(std::string name, double vec[NDIM]);
void print3M(std::string name, double mat[NDIM-1][NDIM-1]);
void print4M(std::string name, double mat[NDIM][NDIM]);