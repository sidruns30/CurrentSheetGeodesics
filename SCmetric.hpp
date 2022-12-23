#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

void GetSCCovariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetSCContravariantMetric(double g[NDIM][NDIM], double X[NDIM]);
void GetSCConnection(double X[NDIM], double gamma[NDIM][NDIM][NDIM]);