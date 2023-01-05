#include "BHAC_MHD.hpp"

/*
    Compute various (ideal) MHD quantities from the primitive inputs.
    The inputs include:
    u1, u2, u3, B1, B2, B3, E1, E2, E3, lfac, P, rho (vectors are cartesianized)
    The metric and the coordinates must be provided
    The following fluid variables can be calculated:
    bfluid[0-3], Bsqr, b^2, u[0-3], T, sigma, beta
    Based on 
*/

/*ARRAY   blfuid0, bfluid1, bfluid2, bfluid3,
        Bsqr, b2, u0, u1, u2, u3, beta, sigma,
        efluid0, efluid1, efluid2, efluid3, Esqr,
        e2, temp;
*/


void BHAC_MHD::Getbfluid(ARRAY &bfluid0, ARRAY &bfluid1, ARRAY &bfluid2, ARRAY &bfluid3, 
                    ARRAY &u0,  ARRAY &u1,  ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    int j, k, l, m, n;
    N = COORDS_BLOCK[0].size();

    bfluid0.clear();
    bfluid1.clear();
    bfluid2.clear();
    bfluid3.clear();

    // Need to call this after 4 velocity (sets MKS u)
    if (u0.size() < N)
    {
        BHAC_MHD::GetU(u0, u1, u2, u3, COORDS_BLOCK, PRIMS_BLOCK);
    }

    for (i=0; i<N; i++)
    {
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];

        // Get MKS coordinates
        double XKS[NDIM], XMKS[NDIM];
        CartToKS(x, y, z, XKS);
        TransformCoordinates(2, 3, XKS, XMKS);
        
        // Get MKS E and B fields
        double EMKS[NDIM-1], ECart[NDIM-1] =    {PRIMS_BLOCK[iprim["E1"]][i], 
                                                PRIMS_BLOCK[iprim["E2"]][i], PRIMS_BLOCK[iprim["E3"]][i]};
        T_3CartTo3MKS(ECart, EMKS, XMKS);
        double BMKS[NDIM-1], BCart[NDIM-1] =    {PRIMS_BLOCK[iprim["B1"]][i], 
                                    PRIMS_BLOCK[iprim["B2"]][i], PRIMS_BLOCK[iprim["B3"]][i]};
        T_3CartTo3MKS(BCart, BMKS, XMKS);


        // Get the metric
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(3, XMKS, gcov);
        GconFunc(3, XMKS, gcon);

        // Lapse, shift (contravariant) and normal vector (contravariant)
        double alpha = sqrt(-1/gcon[0][0]);
        double beta_con[NDIM-1] = {SQR(alpha)*gcon[0][1], SQR(alpha)*gcon[0][2], SQR(alpha)*gcon[0][3]};
        double n_con[NDIM] = {1/alpha, -beta_con[0]/alpha, -beta_con[1]/alpha, -beta_con[2]/alpha};
        double n_cov[NDIM] = {-alpha, 0, 0, 0};

        double lfac = PRIMS_BLOCK[iprim["lfac"]][i];

        double u_con[NDIM] = {u0[i], u1[i], u2[i], u3[i]}, u_cov[NDIM];
        UpperToLower(3, XMKS, u_con, u_cov);
        double v_cov[NDIM-1] = {u_cov[1]/lfac, u_cov[2]/lfac, u_cov[3]/lfac};
        
        double Bivi = BMKS[0]*v_cov[0] + BMKS[1]*v_cov[1] + BMKS[2]*v_cov[2];

        // different approaches for ideal and non ideal MHD
        // taken from eq 19 of 
        // https://link.springer.com/content/pdf/10.1186/s40668-017-0020-2.pdf?pdf=button
        if (idealMHD)
        {   
            double b0 = lfac * Bivi / alpha;
            bfluid0.push_back(b0);
            bfluid1.push_back((BMKS[0] + alpha * b0 * u_con[1]) / lfac);
            bfluid2.push_back((BMKS[1] + alpha * b0 * u_con[2]) / lfac);
            bfluid3.push_back((BMKS[2] + alpha * b0 * u_con[3]) / lfac);
        }

        else
        {
            // Taken from equations 35 and 36 of https://iopscience.iop.org/article/10.3847/1538-4365/ab3922/pdf

            // Get MKS E field (covariant)
            double E4MKS[NDIM] = {0., EMKS[0], EMKS[1], EMKS[2]};
            double E4MKS_cov[NDIM];
            UpperToLower(3, XMKS, E4MKS, E4MKS_cov);

            // Make the 4 vectors
            double v4_cov[NDIM], v4_con[NDIM];
            double B4MKS[NDIM] = {0, BMKS[0], BMKS[1], BMKS[2]};
            
            // eq 24 of N. Bucciantini and L. Del Zanna (2013)
            for (j = 0; j<NDIM; j++)
            {
                v4_con[j] = u_con[j] / lfac - n_con[j];
            }
            
            UpperToLower(3, XMKS, v4_con, v4_cov);
            double b_calc[NDIM];


            for (j=0; j<NDIM; j++)
            {

                // Below eq 27 of N. Bucciantini and L. Del Zanna (2013)
                // Determinant factor of the gamma matric
                double detgam = Detgammacov(3, XMKS);

                // First the contraction with the leveicivita symbol
                double vcrossE = 0.;

                for (k=1; k<NDIM; k++)
                {
                    for (l=1; l<NDIM; l++)
                    {
                        vcrossE += n_cov[0] * GetEta(j,k,l,0) * v4_cov[k] * E4MKS_cov[l] / sqrt(detgam);
                    }
                }

                // Set the b fluid
                    b_calc[j] = lfac * (Bivi * n_con[j] + B4MKS[j] - vcrossE);
            }
            bfluid0.push_back(b_calc[0]);
            bfluid1.push_back(b_calc[1]);
            bfluid2.push_back(b_calc[2]);
            bfluid3.push_back(b_calc[3]);

        }
    }
    return;
}

void BHAC_MHD::Getefluid(ARRAY &efluid0, ARRAY &efluid1, ARRAY &efluid2, ARRAY &efluid3, 
                    ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    int j, k, l, m, n;
    N = COORDS_BLOCK[0].size();

    efluid0.clear();
    efluid1.clear();
    efluid2.clear();
    efluid3.clear();

    // Need to call this after 4 velocity (sets MKS u)
    if (u0.size() < N)
    {
        BHAC_MHD::GetU(u0, u1, u2, u3, COORDS_BLOCK, PRIMS_BLOCK);
    }


    for (i=0; i<N; i++)
    {
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];

        // Get MKS coordinates
        double XKS[NDIM], XMKS[NDIM];
        CartToKS(x, y, z, XKS);
        TransformCoordinates(2, 3, XKS, XMKS);

        // Get MKS E and B fields
        double EMKS[NDIM-1], ECart[NDIM-1] =    {PRIMS_BLOCK[iprim["E1"]][i], 
                                                PRIMS_BLOCK[iprim["E2"]][i], PRIMS_BLOCK[iprim["E3"]][i]};
        T_3CartTo3MKS(ECart, EMKS, XMKS);
        double BMKS[NDIM-1], BCart[NDIM-1] =    {PRIMS_BLOCK[iprim["B1"]][i], 
                                    PRIMS_BLOCK[iprim["B2"]][i], PRIMS_BLOCK[iprim["B3"]][i]};
        T_3CartTo3MKS(BCart, BMKS, XMKS);

        // Get the metric
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(3, XMKS, gcov);
        GconFunc(3, XMKS, gcon);

        // Lapse, shift (contravariant) and normal vector (contravariant)
        double alpha = sqrt(-1/gcon[0][0]);
        double beta_con[NDIM-1] = {SQR(alpha)*gcon[0][1], SQR(alpha)*gcon[0][2], SQR(alpha)*gcon[0][3]};
        double n_con[NDIM] = {1/alpha, -beta_con[0]/alpha, -beta_con[1]/alpha, -beta_con[2]/alpha};
        double n_cov[NDIM] = {-alpha, 0, 0, 0};

        double lfac = PRIMS_BLOCK[iprim["lfac"]][i];

        double u_con[NDIM] = {u0[i], u1[i], u2[i], u3[i]}, u_cov[NDIM];
        UpperToLower(3, XMKS, u_con, u_cov);
        double v_cov[NDIM-1] = {u_cov[1]/lfac, u_cov[2]/lfac, u_cov[3]/lfac};
        
        double Eivi = EMKS[0]*v_cov[0] + EMKS[1]*v_cov[1] + EMKS[2]*v_cov[2];

        if (idealMHD)
        {
            efluid0.push_back(0.);
            efluid1.push_back(0.);
            efluid2.push_back(0.);
            efluid3.push_back(0.);
        }

        else
        {
         // Taken from equations 35 and 36 of https://iopscience.iop.org/article/10.3847/1538-4365/ab3922/pdf


            // Get MKS E field (contravariant)

            // Get MKS E field (covariant)
            double B4MKS[NDIM] = {0., BMKS[0], BMKS[1], BMKS[2]};
            double B4MKS_cov[NDIM];
            UpperToLower(3, XMKS, B4MKS, B4MKS_cov);

            // Make the 4 vectors
            double v4_cov[NDIM], v4_con[NDIM];
            double E4MKS[NDIM] = {0, EMKS[0], EMKS[1], EMKS[2]};
            
            // eq 24 of N. Bucciantini and L. Del Zanna (2013)
            for (j = 0; j<NDIM; j++)
            {
                v4_con[j] = u_con[j] / lfac - n_con[j];
            }
            
            UpperToLower(3, XMKS, v4_con, v4_cov);
            double e_calc[NDIM];

            for (j=0; j<NDIM; j++)
            {

                // Below eq 27 of N. Bucciantini and L. Del Zanna (2013)
                // Determinant factor of the gamma matric
                double detgam = Detgammacov(3, XMKS);

                // First the contraction with the leveicivita symbol
                double vcrossB = 0.;

                for (k=1; k<NDIM; k++)
                {
                    for (l=1; l<NDIM; l++)
                    {
                        vcrossB += n_cov[0] * GetEta(j,k,l,0) * v4_cov[k] * B4MKS_cov[l] / sqrt(detgam);
                    }
                }

                // Set the b fluid
                    e_calc[j] = lfac * (Eivi * n_con[j] + E4MKS[j] + vcrossB);
            }
            
            efluid0.push_back(e_calc[0]);
            efluid1.push_back(e_calc[1]);
            efluid2.push_back(e_calc[2]);
            efluid3.push_back(e_calc[3]);
        }
    }
    return;
}

// Get the 4 velocity arrays
void BHAC_MHD::GetU(ARRAY &u0, ARRAY &u1, ARRAY &u2, ARRAY &u3, 
                    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{   
    size_t i,j,k,N;
    N = COORDS_BLOCK[0].size();
    // Count bad cells
    int Nbad = 0;

    u0.clear();
    u1.clear();
    u2.clear();
    u3.clear();


    for (i=0; i<N; i++)
    {
        double x, y, z, V1, V2, V3, lfac;
        V1 = PRIMS_BLOCK[iprim["V1"]][i];
        V2 = PRIMS_BLOCK[iprim["V2"]][i];
        V3 = PRIMS_BLOCK[iprim["V3"]][i];
        lfac = PRIMS_BLOCK[iprim["lfac"]][i];
        x = COORDS_BLOCK[0][i];
        y = COORDS_BLOCK[1][i];
        z = COORDS_BLOCK[2][i];

        // Transform coordinates to MKS
        double XKS[NDIM], XMKS[NDIM], U[NDIM], VMKS[NDIM-1], V[NDIM-1] = {V1, V2, V3};
        CartToKS(x, y, z, XKS);
        TransformCoordinates(2, 3, XKS, XMKS);

        // Get the MKS matrices
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(3, XMKS, gcov);
        GconFunc(3, XMKS, gcon);

        // Get upper and lower lapse and shift vectors
        double  alpha, beta_con[NDIM-1], beta_cov[NDIM-1], 
                gamma_con[NDIM-1][NDIM-1], gamma_cov[NDIM-1][NDIM-1];
        alpha = sqrt(-1/gcon[0][0]);
        for (j=0;j<NDIM-1;j++)
        {
            beta_con[j] = gcon[0][j+1] * SQR(alpha);
            beta_cov[j] = gcov[0][j+1];
        }
        for (k=0; k<NDIM-1; k++)
        {
            for (j=0; j<NDIM-1; j++)
            {
                gamma_con[k][j] = gcon[k+1][j+1] + beta_con[k] * beta_con[j] / SQR(alpha);
                gamma_cov[k][j] = gcov[k+1][j+1];
            }
        }

        // Flux in the MKS units
        T_3CartTo3MKS(V, VMKS, XMKS);

        // u con (Input is lfac v^i)
        double ui_con[NDIM-1];
        for(j=0;j<NDIM-1;j++)
        {
            ui_con[j] = VMKS[j] - lfac * beta_con[j]/alpha;
        }

        // Solve the quadratic to find the 4 velocity
        double A=0., B=0., C=0.;
        A = gcov[0][0];
        for (j=0;j<NDIM-1;j++)
        {
            B += 2 * gcov[0][j+1] * ui_con[j];
            for (k=0; k<NDIM-1;k++)
            {
                C += gcov[j+1][k+1] * ui_con[j] * ui_con[k];
            }
            U[j+1] = ui_con[j];
        }
        C +=1.;
        // Obtain two roots for the answer
        double root1, root2;
        root1 = (-B + sqrt(SQR(B) - 4*A*C))/(2*A);
        root2 = (-B - sqrt(SQR(B) - 4*A*C))/(2*A);

        // Pick the root that matches the lorentz factor within a tolerance
        double tol = 1.e-3;
        U[0] = (fabs(1 - root1*alpha/lfac) < fabs(1 - root2*alpha/lfac)) ? root1 : root2;

        if (fabs(1 - U[0]*alpha/lfac) > tol)
        {
            Nbad ++;
            U[0] = lfac / alpha;
        }

        // Push to the 4 velocity arrays
        u0.push_back(U[0]);
        u1.push_back(U[1]);
        u2.push_back(U[2]);
        u3.push_back(U[3]);
        
    }
    if (Nbad > 0)
    {
        std::cout<<"WARNING: computed four velocity lfac does not match simulation lfac for "<<
                Nbad<<" out of "<< N << " cells to within 0.1 percent \n";

    }

    return;
}

void BHAC_MHD::GetBsqr(ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, j, N;
    N = COORDS_BLOCK[0].size();

    Bsqr.clear();

    for (i=0; i<N; i++)
    {
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];

        // Get MKS coordinates
        double XKS[NDIM], XMKS[NDIM];
        CartToKS(x, y, z, XKS);
        TransformCoordinates(2, 3, XKS, XMKS);
        
        // Get MKS B fields
        double BMKS[NDIM-1], BCart[NDIM-1] =    {PRIMS_BLOCK[iprim["B1"]][i], 
                                    PRIMS_BLOCK[iprim["B2"]][i], PRIMS_BLOCK[iprim["B3"]][i]};
        T_3CartTo3MKS(BCart, BMKS, XMKS);

        // Get the metric
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(3, XMKS, gcov);
        GconFunc(3, XMKS, gcon);
        
        // Make the covariant and contravariant vector fields
        double B4MKS_con[NDIM] = {0., BMKS[0], BMKS[1], BMKS[2]};
        double B4MKS_cov[NDIM];
        UpperToLower(3, XMKS, B4MKS_con, B4MKS_cov);

        double sum = 0.;
        for (j=0; j<NDIM; j++)
        {
            sum += B4MKS_con[j] * B4MKS_cov[j];
        }
        Bsqr.push_back(sum);
    }
    return;
}
    
void BHAC_MHD::GetEsqr(ARRAY &Esqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, j, N;
    N = COORDS_BLOCK[0].size();

    Esqr.clear();

    for (i=0; i<N; i++)
    {
        double x = COORDS_BLOCK[0][i], y = COORDS_BLOCK[1][i], z = COORDS_BLOCK[2][i];

        // Get MKS coordinates
        double XKS[NDIM], XMKS[NDIM];
        CartToKS(x, y, z, XKS);
        TransformCoordinates(2, 3, XKS, XMKS);
        
        // Get MKS E fields
        double EMKS[NDIM-1], ECart[NDIM-1] =    {PRIMS_BLOCK[iprim["E1"]][i], 
                                    PRIMS_BLOCK[iprim["E2"]][i], PRIMS_BLOCK[iprim["E3"]][i]};
        T_3CartTo3MKS(ECart, EMKS, XMKS);

        // Get the metric
        double gcov[NDIM][NDIM];
        double gcon[NDIM][NDIM];
        GcovFunc(3, XMKS, gcov);
        GconFunc(3, XMKS, gcon);
        
        // Make the covariant and contravariant vector fields
        double E4MKS_con[NDIM] = {0., EMKS[0], EMKS[1], EMKS[2]};
        double E4MKS_cov[NDIM];
        UpperToLower(3, XMKS, E4MKS_con, E4MKS_cov);

        double sum = 0.;
        for (j=0; j<NDIM; j++)
        {
            sum += E4MKS_con[j] * E4MKS_cov[j];
        }

        Esqr.push_back(sum);

    }
    return;
}
    
void BHAC_MHD::Getb2(ARRAY &b2, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}

void BHAC_MHD::Gete2(ARRAY &e2, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}
    
void BHAC_MHD::GetTemp(ARRAY &temp, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    temp.clear();

    for (i=0; i<N; i++)
    {
        temp.push_back(PRIMS_BLOCK[iprim["rho"]][i] / PRIMS_BLOCK[iprim["P"]][i]);
    }    
    return;
}

void BHAC_MHD::GetSigma(ARRAY &sigma, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    sigma.clear();

    // Check if Bsqr is computed
    if (Bsqr.size() != N)
    {
        BHAC_MHD::GetBsqr(Bsqr, COORDS_BLOCK, PRIMS_BLOCK);
    }

    for (i=0; i<N; i++)
    {
        sigma.push_back(Bsqr[i] / PRIMS_BLOCK[iprim["rho"]][i]);
    }

    return;
}

void BHAC_MHD::GetBeta(ARRAY &beta, ARRAY &Bsqr, const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    size_t i, N;
    N = COORDS_BLOCK[0].size();

    beta.clear();
    
    // Check if Bsqr is computed
    if (Bsqr.size() != N)
    {   
        BHAC_MHD::GetBsqr(Bsqr, COORDS_BLOCK, PRIMS_BLOCK);
    }

    for (i=0; i<N; i++)
    {
        beta.push_back(PRIMS_BLOCK[iprim["P"]][i] / Bsqr[i]);
    }    
    return;
}