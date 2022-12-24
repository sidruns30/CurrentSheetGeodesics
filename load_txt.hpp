#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

const int NVARS = 20;

// Define types
typedef std::vector <double> ARRAY;
typedef std::vector <std::vector <double>> ARRAY2D;


ARRAY x,y,z;
ARRAY rho,P,V1,V2,V3,B1,B2,B3,lfac,E1,E2,E3;
ARRAY Bsqr_sim,bfluid0_sim,bfluid1_sim,bfluid2_sim,bfluid3_sim;
ARRAY b2_sim, rMKS_sim, thetaMKS_sim;

// Necessary for field calculations
const bool idealMHD = true;

void InitializeArrays(std::string FILE_NAME);
extern const int NVARS;

extern void WriteVectorToFile(std::string fname, ARRAY2D data);
