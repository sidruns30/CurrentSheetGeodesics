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

ARRAY   blfuid0, bfluid1, bfluid2, bfluid3,
        Bsqr, b2, u0, u1, u2, u3, beta, sigma,
        efluid0, efluid1, efluid2, efluid3, Esqr,
        e2, temp;

void BHAC_MHD::Getbfluid(ARRAY bfluid0, ARRAY bfluid1, ARRAY bfluid2, ARRAY bfluid3, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {

    }
    return;
}

void BHAC_MHD::Getefluid(ARRAY efluid0, ARRAY efluid1, ARRAY efluid2, ARRAY efluid3, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}

// Get the 4 velocity arrays
void BHAC_MHD::GetU(ARRAY u0, ARRAY u1, ARRAY u2, ARRAY u3, ARRAY2D COORDS, ARRAY2D PRIMS)
{    
    size_t i,j,k,N;
    N = COORDS[0].size();

    u0.clear();
    u1.clear();
    u2.clear();
    u3.clear();

    for (i=0; i<N; i++)
    {
        double x, y, z, V1, V2, V3, lfac;
        V1 = PRIMS[iprim["V1"]][i];
        V2 = PRIMS[iprim["V2"]][i];
        V3 = PRIMS[iprim["V3"]][i];
        lfac = PRIMS[iprim["lfac"]][i];
        x = COORDS[0][i];
        y = COORDS[1][i];
        z = COORDS[2][i];

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
        double tol = 1.e-2;
        if (fabs(1 - root1*alpha/lfac) < tol)
        {
            U[0] = root1;
        }
        else if (fabs(1 - root2*alpha/lfac) < tol)
        {
            U[0] = root2;
        }
        else
        {
            print("WARNING: 4 velocity lfac does not match simulation lfac, assigning 0");
            U[0] = 0;
        }

        // Push to the 4 velocity arrays
        u0.push_back(U[0]);
        u1.push_back(U[1]);
        u2.push_back(U[2]);
        u3.push_back(U[3]);
        
    }
    return;
}

void BHAC_MHD::GetBsqr(ARRAY Bsqr, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}
    
void BHAC_MHD::GetEsqr(ARRAY Esqr, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}
    
void BHAC_MHD::Getb2(ARRAY b2, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}

void BHAC_MHD::Gete2(ARRAY e2, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }
    return;
}
    
void BHAC_MHD::GetTemp(ARRAY temp, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }    
    return;
}

void BHAC_MHD::GetSigma(ARRAY sigma, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }    
    return;
}

void BHAC_MHD::GetBeta(ARRAY beta, ARRAY2D coords, ARRAY2D PRIMS)
{
    size_t i, N;
    N = coords[0].size();

    for (i=0; i<N; i++)
    {
        
    }    
    return;
}