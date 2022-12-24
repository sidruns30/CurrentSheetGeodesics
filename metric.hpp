#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

#include "BLmetric.hpp"
#include "KSmetric.hpp"
#include "SCmetric.hpp"
#include "MKSmetric.hpp"

#define SMALL_VECTOR	1.e-30


void GcovFunc(int metric_type, double X[NDIM], double gcov[NDIM][NDIM]);
void GconFunc(int metric_type, double X[NDIM], double gcon[NDIM][NDIM]);
void GetConnectionNumerical(int metric_type, double X[NDIM], double conn[NDIM][NDIM][NDIM]);
void GetConnectionAnalytic(int metric_type, double X[NDIM], double gamma[NDIM][NDIM][NDIM]);
void TransformCoordinates(int metric_type_1, int metric_type_2, double X1[NDIM], double X2[NDIM]);
void TransformFourVectors(int metric_type_1, int metric_type_2, double X1[NDIM], double T1[NDIM], double T2[NDIM]);
double Detgcov(int metric_type, double X[NDIM]);
void GetEta(double ep[NDIM][NDIM][NDIM][NDIM]);

void UpperToLower(int metric_type, double X[NDIM], double Ucon[NDIM], double Ucov[NDIM]);
double delta(int i, int j);
void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov);
bool Get4Velocity(int metric_type, double X[NDIM], double ui_con[NDIM-1], double Ui_con[NDIM], 
					double b_con[NDIM]);
void Test4Velocity( int id, int metric_type, double X[NDIM], double U_con[NDIM], 
                    double b_con[NDIM], double _lfac, double alpha);
void make_tetrad(double Ucon[NDIM],
		 double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
		 double Ecov[NDIM][NDIM]);
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM]);
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM]);
double Norm_con(int metric_type, double X[NDIM], double V_con[NDIM]);
double Norm_cov(int metric_type, double X[NDIM], double V_cov[NDIM]);
double ADotB(int metric_type, int vec_type, double X[NDIM], double a[NDIM], double b[NDIM]);
#define MUNULOOP for(int mu=0; mu < NDIM; mu++)                                \
                 for(int nu=0; nu < NDIM; nu++)