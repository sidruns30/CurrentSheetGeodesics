#include "load_txt.hpp"
#include "MKSmetric.hpp"
#include "geodesic.hpp"

int main()
{

    int i, j;
    std::string FNAME("../input/BHAC_data.txt");
    std::string INFILE("currentsheetdata.txt");
    std::string OUTNAME("wavevectors.txt");

    std::vector <size_t> i_sheet;
    std::vector <double> data[2*NDIM];
    InitializeArrays(FNAME);
    WriteCurrentSheetData(INFILE);
    
    /*
    // Step 1: Get Current Sheet
    std::vector <size_t> i_sheet;
    std::vector <double> data[2*NDIM];
    InitializeArrays(FNAME);
    WriteCurrentSheetData(INFILE);
    */
    // Step 2: Get Wavevectors of photons
    // mode sets the direction of the wave vector
    /* 
    WriteX_kToFile(INFILE,OUTNAME);
    INFILE.clear();
    OUTNAME.clear();
    */
    // Step 3: Integrate the geodesics
    /*
    INFILE = "wavevectors.txt";
    OUTNAME = "geodesics.txt";
    GetGeodesics(INFILE, OUTNAME);
    */
    return 0;
}