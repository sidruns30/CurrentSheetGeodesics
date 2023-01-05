#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif

void GetKSCovariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetKSContravariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetKSConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM]);
void X_KSToBL(double XKS[NDIM], double XBL[NDIM]);
void X_BLToKS(double XBL[NDIM], double XKS[NDIM]);
void T_KSToBL(double TKS[NDIM], double TBL[NDIM], double XKS[NDIM]);
void T_BLToKS(double TBL[NDIM], double TKS[NDIM], double XBL[NDIM]);
void CartToKS(double x, double y, double z, double XKS[NDIM]);
void GetKSField(double x, double y, double z, double b1, double b2, double b3, double &br, 
                double &bth, double &bph);