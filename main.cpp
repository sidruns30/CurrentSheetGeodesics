#include "load_txt.hpp"
#include "BHAC_MHD.hpp"
#include "partition.hpp"
#include "geodesic.hpp"

// To do: read up on openmp

int main()
{

    size_t i, j;
    std::string FNAME("../input/BHAC_data.txt");
    
    ARRAY2D COORDS, PRIMS;
    InitializeArrays(COORDS, PRIMS, FNAME);

    // For openmp COORDS, PRIMS -> COORDS_BLOCKS, PRIMS_BLOCKS
    PartitionGrid(COORDS, PRIMS);

    // Determine how to make the wavevectors (parallel, antiparallel, perp)
    int mode = 0;

    #pragma omp parallel for schedule(static)
    for (i=0; i<maxthreads; i++)
    {
        ARRAY2D COORDS_BLOCK = COORDS_BLOCKS[i];
        ARRAY2D PRIMS_BLOCK = PRIMS_BLOCKS[i];

        // Find the current sheet indices for each block and initialize geodesics
        std::vector <std::vector <size_t>> indices;
        FindCurrentSheet(indices, COORDS_BLOCK, PRIMS_BLOCK);
        
        // Now get the wavevectors
        ARRAY2D X_K;
        ConstructWavevectors(X_K, mode, indices, COORDS_BLOCK, PRIMS_BLOCK);
        //MakeGeodesics(COORDS_BLOCK, PRIMS_BLOCK);
    }


    // Step 3: Integrate the geodesics
    /*
    INFILE = "wavevectors.txt";
    OUTNAME = "geodesics.txt";
    GetGeodesics(INFILE, OUTNAME);
    */
    return 0;
}