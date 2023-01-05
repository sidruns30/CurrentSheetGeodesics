#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif

#include "KSmetric.hpp"

extern void TransformCoordinates(int metric_type_1, int metric_type_2, double X1[NDIM], double X2[NDIM]);

void GetMKSCovariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetMKSContravariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetMKSConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM]);
void X_MKSToKS(double XMKS[NDIM], double XKS[NDIM]);
void T_MKSToKS(double TMKS[NDIM], double TKS[NDIM], double XMKS[NDIM]);
void T_KSToMKS(double TKS[NDIM], double TMKS[NDIM], double XMKS[NDIM]);
void X_KSToMKS(double XKS[NDIM], double XMKS[NDIM]);

// The transcedental versions
void GetMKSCovariantMetric_v2(double g[NDIM][NDIM], double X[NDIM]);
void GetMKSContravariantMetric_v2(double g[NDIM][NDIM], double X[NDIM]);
void X_KSToMKS_v2(double XKS[NDIM], double XMKS[NDIM]);
void X_MKSToKS_v2(double XMKS[NDIM], double XKS[NDIM]);
void T_KSToMKS_v2(double TKS[NDIM], double TMKS[NDIM], double XKS[NDIM]);
void T_MKSToKS_v2(double TMKS[NDIM], double TKS[NDIM], double XMKS[NDIM]);


void T_3MKSTo3Cart(double T3MKS[NDIM-1], double T3Cart[NDIM-1], double XMKS[NDIM]);
void T_3CartTo3MKS(double T3Cart[NDIM-1], double T3MKS[NDIM-1], double XMKS[NDIM]);
void T_3CartTo3MKS_v2(double T3Cart[NDIM-1], double T3MKS[NDIM-1], double XMKS[NDIM]);


extern bool Invert3Matrix(const double m[3][3], double minv[3][3]);
extern void Multiply3Matrices(double a[NDIM-1][NDIM-1], double b[NDIM-1][NDIM-1], double c[NDIM-1][NDIM-1]);
extern void Multiply4Matrices(double a[NDIM][NDIM], double b[NDIM][NDIM], double c[NDIM][NDIM]);
