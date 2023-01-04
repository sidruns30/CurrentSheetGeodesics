#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

// Coordinates also imports metric
#include "metric.hpp"
#include "BHAC_MHD.hpp"



const int NVARS_SHEET = 19;

/*
extern ARRAY x,y,z;
extern ARRAY V1, V2, V3, B1, B2, B3, E1, E2, E3, lfac, P, rho;
extern ARRAY Bsqr_sim,bfluid0_sim,bfluid1_sim,bfluid2_sim,bfluid3_sim;
extern ARRAY b2_sim, rMKS_sim, thetaMKS_sim;
*/

// Main functions to be called in main
void FindCurrentSheet(std::vector <std::vector <size_t>> &indices, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
void ConstructWavevectors(ARRAY2D &X_K, int mode,
                        const std::vector <std::vector <size_t>> &indices, 
                        ARRAY2D &COORDS_BLOCK, ARRAY2D &PRIMS_BLOCK);