#include "metric.hpp"
#include <unordered_map>

// The coords and primitives are given by the reader. The file only requires to be linked to the
// metric file

// Define a typedef
typedef std::vector <std::vector <double>> ARRAY;

// Set the coordinates:
extern std::vector <double> x, y, z;
ARRAY COORDS = {x, y, z};

// Set the primitives:
extern std::vector <double> u1, u2, u3, B1, B2, B3, E1, E2, E3, lfac, P, rho;
ARRAY PRIMS= {u1, u2, u3, B1, B2, B3, E1, E2, E3, lfac, P, rho};

// Create a map to access PRIMS by array names
std::unordered_map <std::string, int> iprim =   {{"u1",0}, {"u2",1},{"u3",2}, {"B1",3},{"B2",4}, {"B3",5},
                                            {"E1",6}, {"E2",7},{"E3",8}, {"lfac",9},{"P",10}, {"rho",11}
                                            }; 

// Set the other fluid quantities to be computed
std::vector <double>    blfuid0, bfluid1, bfluid2, bfluid3,
                        Bsqr, b2, u0, u1, u2, u3, beta, sigma,
                        efluid0, efluid1, efluid2, efluid3, Esqr,
                        e2;

// Return global fluid quantities
namespace BHAC_MHD
{
    std::vector <double> Getbfluid(ARRAY coords, ARRAY PRIMS);
    std::vector <double> Getefluid(ARRAY coords, ARRAY PRIMS);
    std::vector <double> GetBsqr(ARRAY coords, ARRAY PRIMS);
    std::vector <double> GetEsqr(ARRAY coords, ARRAY PRIMS);
    std::vector <double> GetU(ARRAY coords, ARRAY PRIMS);
    std::vector <double> Getb2(ARRAY coords, ARRAY PRIMS);
    std::vector <double> Gete2(ARRAY coords, ARRAY PRIMS);
    std::vector <double> GetTemp(ARRAY coords, ARRAY PRIMS);    
    std::vector <double> GetSigma(ARRAY coords, ARRAY PRIMS);
    std::vector <double> GetBeta(ARRAY coords, ARRAY PRIMS);
};