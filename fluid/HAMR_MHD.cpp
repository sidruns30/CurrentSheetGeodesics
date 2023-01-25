#include "HAMR_MHD.hpp"

/*
    Similar to BHAC MHD but compute fluid quantities for HAMR data
    Analysis is done in regular KS coordinates
*/

void HAMR_MHD::GetU(ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3,
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    print("Computing four velocity");
    // Input Vis are v^i in Cartesian coordinates
    size_t i, N = COORDS_BLOCK[0].size();
    u0.clear();
    u1.clear();
    u2.clear();
    u3.clear();

    size_t bad = 0;
    for (i=0; i<N; i++)
    {
        // First get the coordinates in KS
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];
        double XKS[NDIM];
        CartToKS(x, y, z, XKS);

        // Now get the arrays in KS
        double Vinp[NDIM-1] = {PRIMS_BLOCK[iVX][i], PRIMS_BLOCK[iVY][i], PRIMS_BLOCK[iVZ][i]};
        double vi_con[NDIM-1], vi_cov[NDIM-1];

        T_3CartTo3_KS(Vinp, vi_con, XKS);
        UpperToLower3(2, XKS, vi_con, vi_cov);

        // Now solve for the lfac
        double sum = (vi_con[0]*vi_cov[0] + vi_con[1]*vi_cov[1] + vi_con[2]*vi_cov[2]);
        double _lfac = 1 / sqrt(1 - sum);

        // Metric stuff
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(2, XKS, gcov);
        GconFunc(2, XKS, gcon);

        // Lapse, shift (contravariant) and normal vector (contravariant)
        double alpha = sqrt(-1/gcon[0][0]);
        double beta_con[NDIM-1] = {SQR(alpha)*gcon[0][1], SQR(alpha)*gcon[0][2], SQR(alpha)*gcon[0][3]};

        // Finally the four velocities
        double _u0 = _lfac / alpha;
        double _u1 = _lfac * (vi_con[0] - beta_con[0] / alpha);     
        double _u2 = _lfac * (vi_con[1] - beta_con[1] / alpha);     
        double _u3 = _lfac * (vi_con[2] - beta_con[2] / alpha);     

        // Sanity check for the dot product
        double u_con[NDIM] = {_u0, _u1, _u2, _u3};
        sum = Norm_con(2, XKS, u_con);
        
        if (std::isnan(sum) || abs(sum + 1) > 1.e-2)
        {
            bad++;
        }

        u0.push_back(_u0);
        u1.push_back(_u1);
        u2.push_back(_u2);
        u3.push_back(_u3);

    }
    std::cout<<(double)100*bad/N<<" percent bad cells"<<std::endl;
    return;
}

void HAMR_MHD::GetTemp(ARRAY &temp, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    print("Computing temp");
    size_t i, N = COORDS_BLOCK[0].size();
    temp.clear();

    for (i=0; i<N; i++)
    {
        double _rho = PRIMS_BLOCK[iRHO][i], _P = PRIMS_BLOCK[iP][i];
        temp.push_back(_rho / _P);
    }
    return;
}

void HAMR_MHD::GetBsqr(ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    print("Computing bsqr");
    size_t i, N = COORDS_BLOCK[0].size();
    Bsqr.clear();

    for (i=0; i<N; i++)
    {
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];

        // Get KS coordinates
        double XKS[NDIM];
        CartToKS(x, y, z, XKS);
        
        // Now get the arrays in KS
        double Binp[NDIM-1] = {PRIMS_BLOCK[iBX][i], PRIMS_BLOCK[iBY][i], PRIMS_BLOCK[iBZ][i]};
        double Bi_con[NDIM-1], Bi_cov[NDIM-1];

        T_3CartTo3_KS(Binp, Bi_con, XKS);
        UpperToLower3(2, XKS, Bi_con, Bi_cov);
        
        double sum = (Bi_con[0]*Bi_cov[0] + Bi_con[1]*Bi_cov[1] + Bi_con[2]*Bi_cov[2]);
        Bsqr.push_back(sum);
    }
    return;
}

void HAMR_MHD::GetSigma(ARRAY &sigma, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    print("Computing sigma");
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    sigma.clear();

    // Check if Bsqr is computed
    if (Bsqr.size() != N)
    {
        HAMR_MHD::GetBsqr(Bsqr, COORDS_BLOCK, PRIMS_BLOCK);
    }

    for (i=0; i<N; i++)
    {
        sigma.push_back(Bsqr[i] / PRIMS_BLOCK[iRHO][i]);
    }

    return;
}

void HAMR_MHD::GetBeta(ARRAY &beta, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    print("Computing beta");
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    beta.clear();
    
    // Check if Bsqr is computed
    if (Bsqr.size() != N)
    {   
        HAMR_MHD::GetBsqr(Bsqr, COORDS_BLOCK, PRIMS_BLOCK);
    }

    for (i=0; i<N; i++)
    {
        beta.push_back(PRIMS_BLOCK[iP][i] / Bsqr[i]);
    }    
    return;
}