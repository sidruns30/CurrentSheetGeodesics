#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif

// Coordinates also imports metric
#include "../metric/metric.hpp"
#include "BHAC_MHD.hpp"

// Main functions to be called in main
void FindCurrentSheet(std::vector <std::vector <size_t>> &indices, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
void ConstructWavevectors(ARRAY2D &X_K, int mode,
                        const std::vector <std::vector <size_t>> &indices, 
                        ARRAY2D &COORDS_BLOCK, ARRAY2D &PRIMS_BLOCK);