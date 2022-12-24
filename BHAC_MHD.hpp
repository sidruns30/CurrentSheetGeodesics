#include "metric.hpp"
#include "load_txt.hpp"
#include <unordered_map>

// The coords and primitives are given by the reader. The file only requires to be linked to the
// metric file

// Define a typedef
typedef std::vector <std::vector <double>> ARRAY2D;
typedef std::vector <double> ARRAY;

// Check if MHD is ideal
extern const int idealMHD;

// Set the coordinates:
extern ARRAY x, y, z;
ARRAY2D COORDS = {x, y, z};

// Set the primitives:
extern ARRAY V1, V2, V3, B1, B2, B3, E1, E2, E3, lfac, P, rho;
ARRAY2D PRIMS= {V1, V2, V3, B1, B2, B3, E1, E2, E3, lfac, P, rho};

// Create a map to access PRIMS by array names
std::unordered_map <std::string, int> iprim =   {{"V1",0}, {"V2",1},{"V3",2}, {"B1",3},{"B2",4}, {"B3",5},
                                            {"E1",6}, {"E2",7},{"E3",8}, {"lfac",9},{"P",10}, {"rho",11}
                                            }; 

// Set the other fluid quantities to be computed
ARRAY   blfuid0, bfluid1, bfluid2, bfluid3,
        Bsqr, b2, u0, u1, u2, u3, beta, sigma,
        efluid0, efluid1, efluid2, efluid3, Esqr,
        e2, temp;

// Return global fluid quantities
namespace BHAC_MHD
{
    void Getbfluid(ARRAY bfluid0, ARRAY bfluid1, ARRAY bfluid2, ARRAY bfluid3, ARRAY2D coords, ARRAY2D PRIMS);
    void Getefluid(ARRAY efluid0, ARRAY efluid1, ARRAY efluid2, ARRAY efluid3, ARRAY2D coords, ARRAY2D PRIMS);
    void GetU(ARRAY u0, ARRAY u1, ARRAY u2, ARRAY u3, ARRAY2D coords, ARRAY2D PRIMS);

    void GetBsqr(ARRAY Bsqr, ARRAY2D coords, ARRAY2D PRIMS);
    void GetEsqr(ARRAY Esqr, ARRAY2D coords, ARRAY2D PRIMS);
    
    void Getb2(ARRAY b2, ARRAY2D coords, ARRAY2D PRIMS);
    void Gete2(ARRAY e2, ARRAY2D coords, ARRAY2D PRIMS);
    
    void GetTemp(ARRAY temp, ARRAY2D coords, ARRAY2D PRIMS);    
    void GetSigma(ARRAY sigma, ARRAY2D coords, ARRAY2D PRIMS);
    void GetBeta(ARRAY beta, ARRAY2D coords, ARRAY2D PRIMS);
};