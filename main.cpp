#include "input/load_txt.hpp"
#include "partition.hpp"
#include "geodesics/geodesic.hpp"

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

    ARRAY u0, u1, u2, u3, Bsqr, sigma, beta, rho, x, y, z;
    #pragma omp parallel for ordered schedule(static,1)
    for (i=0; i<maxthreads; i++)
    {
        ARRAY2D COORDS_BLOCK = COORDS_BLOCKS[i];
        ARRAY2D PRIMS_BLOCK = PRIMS_BLOCKS[i];

        ARRAY u0_loc, u1_loc, u2_loc, u3_loc, x_loc, y_loc, z_loc;
        ARRAY Bsqr_loc, sigma_loc, beta_loc, rho_loc;
        ARRAY bfluid0_loc, bfluid1_loc, bfluid2_loc, bfluid3_loc, b2_loc;

        HAMR_MHD::GetU(u0_loc, u1_loc, u2_loc, u3_loc, COORDS_BLOCK, PRIMS_BLOCK);
        HAMR_MHD::GetBsqr(Bsqr_loc, COORDS_BLOCK, PRIMS_BLOCK);
        HAMR_MHD::GetSigma(sigma_loc, Bsqr_loc, COORDS_BLOCK, PRIMS_BLOCK);
        HAMR_MHD::GetBeta(beta_loc, Bsqr_loc, COORDS_BLOCK, PRIMS_BLOCK);
        HAMR_MHD::Getbfluid(bfluid0_loc, bfluid1_loc, bfluid2_loc, bfluid3_loc, 
                            b2_loc, u0_loc, u1_loc, u2_loc, u3_loc, COORDS_BLOCK, PRIMS_BLOCK);
        CopyArray(rho_loc, PRIMS_BLOCK[iRHO]);
        CopyArray(x_loc, COORDS_BLOCK[0]);
        CopyArray(y_loc, COORDS_BLOCK[1]);
        CopyArray(z_loc, COORDS_BLOCK[2]);

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
        */
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
       }
    }

    WriteGeodesics(gdscs, OUTNAME);
    WriteVectorToNumpyArray("../input/u0.npy", u0);
    WriteVectorToNumpyArray("../input/u1.npy", u1);
    WriteVectorToNumpyArray("../input/u2.npy", u2);
    WriteVectorToNumpyArray("../input/u3.npy", u3);
    WriteVectorToNumpyArray("../input/bsqr.npy", Bsqr);
    WriteVectorToNumpyArray("../input/sigma.npy", sigma);
    WriteVectorToNumpyArray("../input/beta.npy", beta);
    WriteVectorToNumpyArray("../input/rho.npy", rho);
    WriteVectorToNumpyArray("../input/x.npy", x);
    WriteVectorToNumpyArray("../input/y.npy", y);
    WriteVectorToNumpyArray("../input/z.npy", z);

    ARRAY2D data = {x, y, z, u0, rho, Bsqr, sigma, beta};
    WriteVectorToFile("../input/text_output.txt", data);
    return 0;
}