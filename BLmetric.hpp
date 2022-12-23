#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

void GetBLCovariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetBLContravariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetBLConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM]);
