/*
    Code to partition the grid into various phi slices for parallelization 
    General structure:
    - Take input from load_txt to partition the grid
    - compute phi over the grid
    - sort phi
    - get a dphi overlap between branches
    - flag all cells that are within the dphi (on the right side)
    - store indices of each block
    - make new code block vectors
*/
#include "partition.hpp"

ARRAY3D COORDS_BLOCKS;
ARRAY3D PRIMS_BLOCKS;

// Return an array of arrays of indices over which to partition the grid
// Also clears the input COORDS and PRIMS arrays
void PartitionGrid(ARRAY2D &COORDS, ARRAY2D &PRIMS)
{
    std::cout<<"Partitioning grid"<<std::endl;
    ARRAY x, y, z, rho, P, B1, B2, B3, E1, E2, E3, V1, V2, V3;

    x = COORDS[0];
    y = COORDS[1];
    ARRAY phi;

    int i, j, N = x.size();

    // Populate phi
    for (i=0; i<N; i++)
    {
        double _x = x[i];
        double _y = y[i];
        phi.push_back(atan2(_y, _x));
    }

    // Get the sorted indices of phi (useful for finding cells in global arrays)
    // Only sort in phi when openmp is enabled
    std::vector <size_t> phi_indices(phi.size());
    #if defined(_OPENMP)
    #else
        COORDS_BLOCKS.push_back(COORDS);
        PRIMS_BLOCKS.push_back(PRIMS);
        COORDS.clear();
        PRIMS.clear();
        iota(phi_indices.begin(), phi_indices.end(), 0);
        std::cout << "Program not called with OpenMP, skipping partition " << std::endl;
        return;
    #endif

    phi_indices = sort_indices(phi);
    // Partition based on resources
    const int block_size = N / maxthreads;
    const int extra_cells = N % maxthreads;

    if (extra_cells > 0)
    {
        std::cout<< "Warning: ignoring "<< extra_cells <<" cells from incomplete division \n";
    }

    size_t ind;
    x.clear();
    y.clear();
    
    for (i=0; i<maxthreads; i++)
    {
        ARRAY2D COORDS_BLOCKS_local;
        ARRAY2D PRIMS_BLOCKS_local;
        for (j=0; j<block_size; j++)
        {
            ind = phi_indices[i*block_size + j];
            x.push_back(COORDS[0][ind]);
            y.push_back(COORDS[1][ind]);
            z.push_back(COORDS[2][ind]);
            rho.push_back(PRIMS[iRHO][ind]);
            P.push_back(PRIMS[iP][ind]);
            B1.push_back(PRIMS[iBX][ind]);
            B2.push_back(PRIMS[iBY][ind]);
            B3.push_back(PRIMS[iBZ][ind]);
            V1.push_back(PRIMS[iVX][ind]);
            V2.push_back(PRIMS[iVY][ind]);
            V3.push_back(PRIMS[iVZ][ind]);
            if (COORDS.size() == 11)
            {
                E1.push_back(PRIMS[iEX][ind]);
                E2.push_back(PRIMS[iEY][ind]);
                E3.push_back(PRIMS[iEZ][ind]);
            }
        }

        COORDS_BLOCKS_local.push_back(x);
        COORDS_BLOCKS_local.push_back(y);
        COORDS_BLOCKS_local.push_back(z);

        PRIMS_BLOCKS_local.push_back(rho);
        PRIMS_BLOCKS_local.push_back(P);
        PRIMS_BLOCKS_local.push_back(V1);
        PRIMS_BLOCKS_local.push_back(V2);
        PRIMS_BLOCKS_local.push_back(V3);
        PRIMS_BLOCKS_local.push_back(B1);
        PRIMS_BLOCKS_local.push_back(B2);
        PRIMS_BLOCKS_local.push_back(B3);
        if (COORDS.size() == 1)
        {
            PRIMS_BLOCKS_local.push_back(E1);
            PRIMS_BLOCKS_local.push_back(E2);
            PRIMS_BLOCKS_local.push_back(E3);
        }


        COORDS_BLOCKS.push_back(COORDS_BLOCKS_local);
        PRIMS_BLOCKS.push_back(PRIMS_BLOCKS_local);

        x.clear();
        y.clear();
        z.clear();
        rho.clear();
        P.clear();
        B1.clear();
        B2.clear();
        B3.clear();
        E1.clear();
        E2.clear();
        E3.clear();
        V1.clear();
        V2.clear();
        V3.clear();
    }

    COORDS.clear();
    PRIMS.clear();
    return;
}