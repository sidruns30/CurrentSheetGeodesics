#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "../defs.hpp"
#endif

#include "../metric/metric.hpp"


// The coords and primitives are given by the reader. The file only requires to be linked to the
// metric file

// Check if MHD is ideal
extern const bool idealMHD;


// Return global fluid quantities
namespace HAMR_MHD
{
    void Getbfluid(ARRAY &bfluid0, ARRAY &bfluid1, ARRAY &bfluid2, ARRAY &bfluid3, ARRAY &b2,
                    ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);

    void Getefluid(ARRAY &efluid0, ARRAY &efluid1, ARRAY &efluid2, ARRAY &efluid3, 
                    ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
    void GetU(ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);

    void GetBsqr(ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
    void GetEsqr(ARRAY &Esqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
    
    void GetTemp(ARRAY &temp, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);    
    void GetSigma(ARRAY &sigma, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
    void GetBeta(ARRAY &beta, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK);
};