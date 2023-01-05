#include "input/load_txt.hpp"
#include "fluid/BHAC_MHD.hpp"
#include "partition.hpp"
#include "geodesics/geodesic.hpp"

int main()
{

    size_t i, j;
    std::string FNAME("../input/BHAC_data.txt");
    std::string OUTNAME("geodesics.txt");
    
    ARRAY2D COORDS, PRIMS;
    InitializeArrays(COORDS, PRIMS, FNAME);

    // For openmp COORDS, PRIMS -> COORDS_BLOCKS, PRIMS_BLOCKS
    PartitionGrid(COORDS, PRIMS);

    // Determine how to make the wavevectors (parallel, antiparallel, perp)
    ARRAY3D gdscs;
    int mode = 0;

    #pragma omp parallel for
    for (i=0; i<maxthreads; i++)
    {
        ARRAY2D COORDS_BLOCK = COORDS_BLOCKS[i];
        ARRAY2D PRIMS_BLOCK = PRIMS_BLOCKS[i];

        // Find the current sheet indices for each block
        std::vector <std::vector <size_t>> indices;
        FindCurrentSheet(indices, COORDS_BLOCK, PRIMS_BLOCK);
        
        // Now get the positions and wavevectors of geodesics
        ARRAY2D X_K;
        ConstructWavevectors(X_K, mode, indices, COORDS_BLOCK, PRIMS_BLOCK);

        // Finally integrate the geodesics
        ARRAY3D gdscs_local;
        GetGeodesics(gdscs_local, X_K);

        #pragma omp critical
        gdscs.insert(gdscs.end(), gdscs_local.begin(), gdscs_local.end());
    }

    WriteGeodesics(gdscs, OUTNAME);

    return 0;
}