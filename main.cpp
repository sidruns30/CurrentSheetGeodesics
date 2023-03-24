#include "input/load_txt.hpp"
#include "partition.hpp"
#include "geodesics/geodesic.hpp"
#include "regrid.hpp"

int main()
{

    int i;
    //std::string FNAME("../input/BHAC_data.txt");
    std::string FNAME("../input/npy_data");
    std::string OUTNAME("../geodesics.txt");
    
    ARRAY2D COORDS, PRIMS;

    std::cout<<"Loading data \n";
    InitializeNumpyArrays(COORDS, PRIMS, FNAME);
    
    // For openmp COORDS, PRIMS -> COORDS_BLOCKS, PRIMS_BLOCKS
    PartitionGrid(COORDS, PRIMS);

    // Determine how to make the wavevectors (parallel, antiparallel, perp)
    ARRAY3D gdscs;
    //int mode = 0;

    //ARRAY u0, u1, u2, u3, Bsqr, b2, sigma, beta, rho, x, y, z;
    ARRAY b_mag_mean, mean_b_mag;
    #if defined(_OPENMP)
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(8);
    #endif
    #pragma omp parallel for ordered schedule(static,1)
    for (i=0; i<maxthreads; i++)
    {
        ARRAY2D COORDS_BLOCK = COORDS_BLOCKS[i];
        ARRAY2D PRIMS_BLOCK = PRIMS_BLOCKS[i];
        ARRAY b_mag_mean_loc, mean_b_mag_loc;
        double dl = 1.e-1;
        Grid *p_Grid = new Grid(COORDS_BLOCK, dl);
        p_Grid->ComputeQuantity("vec_mag_mean", PRIMS_BLOCK[iBX], PRIMS_BLOCK[iBY], PRIMS_BLOCK[iBZ]);
        b_mag_mean_loc = p_Grid->GetGridData();
        p_Grid->ComputeQuantity("mean_vec_mag", PRIMS_BLOCK[iBX], PRIMS_BLOCK[iBY], PRIMS_BLOCK[iBZ]);
        mean_b_mag_loc = p_Grid->GetGridData();
        #pragma omp ordered
        {
            AppendArray(b_mag_mean, b_mag_mean_loc);
            AppendArray(mean_b_mag, mean_b_mag_loc);
        }
        /*
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
       #pragma omp ordered
       {
        //GetArrayInfo("rho_loc", rho_loc);
        //GetArrayInfo("rho original", PRIMS_BLOCK[iRHO]);
        AppendArray(u0, u0_loc);
        AppendArray(u1, u1_loc);
        AppendArray(u2, u2_loc);
        AppendArray(u3, u3_loc);
        AppendArray(sigma, sigma_loc);
        AppendArray(Bsqr, Bsqr_loc);
        AppendArray(beta, beta_loc);
        AppendArray(rho, rho_loc);
        AppendArray(x, x_loc);
        AppendArray(y, y_loc);
        AppendArray(z, z_loc);
        AppendArray(b2, b2_loc);
       }
       */
      p_Grid->GridToVTK("out.vtk", false);
      p_Grid->ClearCells();
    }
    //WriteVectorToNumpyArray("../input/mean_b_mag_dl_0.03.npy", mean_b_mag);
    //WriteVectorToNumpyArray("../input/b_mag_mean_dl_0.03.npy", b_mag_mean);
    //WriteGeodesics(gdscs, OUTNAME);
    return 0;
}