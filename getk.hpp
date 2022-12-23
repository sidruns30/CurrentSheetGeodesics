#ifndef DEF_HEAD
#define DEF_HEAD (1)
#include "defs.hpp"
#endif

// Coordinates also imports metric
#include "metric.hpp"

// The sigma threshold for the current sheet
const double SIGMA_THRESHOLD = 3.;
// The size of the box around each cell to to check 
// for the upstream magnetic field
const double BOX_THRESHOLD = 0.01;
// Minimum upstream magnetic field strength
const double BSQR_THRESHOLD = 10.;

const int NVARS_SHEET = 19;

extern std::vector<double> x,y,z;
extern std::vector<double> rho,p,u1,u2,u3,b1,b2,b3;
extern std::vector<double> Bsqr,bfluid0,bfluid1,bfluid2,bfluid3;
extern std::vector<double> lfac, B2, thetaMKS, rMKS;

void ReadSheetDataFromFile (std::string FILE_NAME, std::vector <double> data[11]);
void WriteCurrentSheetData(std::string Outfile);
void WriteX_kToFile(std::string Infile, std::string OutName);