#include "getk.hpp"

// Find all cells which have lie within a box_threshold of the cell
// in the current sheet

void GetClosebyCells(size_t ind, std::vector<size_t> &goodindices,
                    const ARRAY2D &COORDS_BLOCK)
{
    double xcoord = COORDS_BLOCK[0][ind];
    double ycoord = COORDS_BLOCK[1][ind];
    double zcoord = COORDS_BLOCK[2][ind];
    goodindices.clear();
    
    std::vector <size_t> _temp1, _temp2;
    size_t i, j, k, N;
    N = COORDS_BLOCK[0].size();

    // search along x
    for (i=0; i<N; i++)
    {
        if (fabs(COORDS_BLOCK[0][i] - xcoord) < BOX_THRESHOLD)
        {
            _temp1.push_back(i);
        }
    }

    // search along y
    for (i=0; i<_temp1.size(); i++)
    {
        if (fabs(COORDS_BLOCK[1][_temp1[i]] - ycoord) < BOX_THRESHOLD)
        {
            _temp2.push_back(_temp1[i]);
        }
    }

    // search along z
    for (i=0; i<_temp2.size(); i++)
    {
        if (fabs(COORDS_BLOCK[2][_temp2[i]] - zcoord) < BOX_THRESHOLD)
        {
            goodindices.push_back(_temp2[i]);
        }
    }

    return;
}

// Return an array of tuples of indices where first index corresponds to
// cell in the current sheet and the second index corresponds to the index 
// of the nearest cell not in the current sheet, which has a high B field
void FindCurrentSheet(std::vector <std::vector <size_t>> &indices, 
    const ARRAY2D &COORDS_BLOCK, const ARRAY2D &PRIMS_BLOCK)
{
    // indices of cells containing the current sheet
    std::vector <size_t> i_sheet;
    size_t i, j, k, N;
    N = PRIMS_BLOCK[iprim["rho"]].size();

    ARRAY Bsqr;
    BHAC_MHD::GetBsqr(Bsqr, COORDS_BLOCK, PRIMS_BLOCK);
    
    // Get all of the current sheet
    for (i=0; i<N; i++)
    {
        double sigma = Bsqr[i]/PRIMS_BLOCK[iprim["rho"]][i];
        if (sigma < SIGMA_THRESHOLD)
        {
            i_sheet.push_back(i);
        }
    }

    indices.clear();
    // Now find upstream magnetic field cells close to the current sheet cells
    N = i_sheet.size();
    for (i=0; i<N; i++)
    {
        std::vector<size_t> closecells;
        size_t imax;
        double x, y, z, r, bsqr_max=0.;
        bool flag = false;

        x = COORDS_BLOCK[0][i_sheet[i]];
        y = COORDS_BLOCK[1][i_sheet[i]];
        z = COORDS_BLOCK[2][i_sheet[i]];

        r = sqrt(SQR(x) + SQR(y) + SQR(z));

        GetClosebyCells(i_sheet[i], closecells, COORDS_BLOCK);
        for (k=0; k<closecells.size(); k++)
        {
            if (Bsqr[closecells[k]] > BSQR_THRESHOLD/r)
            {
                if (Bsqr[closecells[k]] > bsqr_max)
                {
                    bsqr_max = Bsqr[closecells[k]];
                    imax = closecells[k];
                    flag = true;
                }
            }
        }

        if (flag)
        {
            std::vector <size_t> temp = {i_sheet[i], imax};
            indices.push_back(temp);
        }
    }
    
    printvar("Number of cells in the current sheet: ", indices.size());

    return;
}

// Returns the wavevector in MKS
void Initialize_XK(ARRAY &X_K, double _x, double _y, double _z, 
                                    double _u0, double _u1, double _u2, double _u3,
                                    double _bfluid0, double _bfluid1, double _bfluid2, double _bfluid3,
                                    int mode)
{
    // mode = {0, 1, 2} for parallel, antiparallel and perpendicular wavevectors
    size_t i,j;

    // MKS
    int metric_type = 3;

    // Get MKS coordinates
    double XKS[NDIM], XMKS[NDIM];
    CartToKS(_x, _y, _z, XKS);
    TransformCoordinates(2, 3, XKS, XMKS);

    // Get both upper and lower metrics
    double gcov[NDIM][NDIM];
    double gcon[NDIM][NDIM];

    GcovFunc(metric_type, XMKS, gcov);
    GconFunc(metric_type, XMKS, gcon);


    // Input fields and 4 velocity
    double U_con_MKS[NDIM] = {_u0, _u1, _u2, _u3},
            b_con[NDIM] = {_bfluid0, _bfluid1, _bfluid2, _bfluid3};

    // Make the tetrad
    double Econ[NDIM][NDIM];
    double Ecov[NDIM][NDIM];
    make_tetrad(U_con_MKS, gcov, Econ, Ecov);


    double b_con_tetrad[NDIM];
    // Get bfluid in the fluid drame
    coordinate_to_tetrad(Ecov, b_con, b_con_tetrad);

    // The photon wavevector in MKS
    double k_con_tetrad[NDIM], k_con_MKS[NDIM], k_cov_tetrad[NDIM], k_cov_MKS[NDIM];
    // parallel
    if (mode == 0)
    {
        for (i=1;i<NDIM;i++)
        {
            k_con_tetrad[i] = b_con_tetrad[i];
        }
    }
    // anti parallel
    else if (mode == 1)
    {
        for (i=1; i<NDIM;i++)
        {
            k_con_tetrad[i] = -b_con_tetrad[i];
        }
    }

    // perpendicular
    else if (mode == 2)
    {
        k_con_tetrad[1] = 1.;
        k_con_tetrad[2] = 0.;
        k_con_tetrad[3] = -b_con_tetrad[1]/b_con_tetrad[3];
    }

    // Normalize the wavevector
    k_con_tetrad[0] = sqrt(SQR(k_con_tetrad[1]) + SQR(k_con_tetrad[2]) + SQR(k_con_tetrad[3]));

    // Get the coordinate frame
    tetrad_to_coordinate(Econ, k_con_tetrad, k_con_MKS);

    // Write to full vector
    X_K.clear();
    for (i=0; i<NDIM; i++)
    {
        X_K.push_back(XMKS[i]);
    }
    for (i=0; i<NDIM; i++)
    {
        X_K.push_back(k_con_MKS[i]);
    }

    return;
}

// Return an array of arrays X_k that contains the initial position and 
// wavevector of the photon. Needs the indices of the cells in the current
// sheet as well as the indices of the upstream cells
void ConstructWavevectors(ARRAY2D &X_K, int mode, 
                        const std::vector <std::vector <size_t>> &indices, 
                        ARRAY2D &COORDS_BLOCK, ARRAY2D &PRIMS_BLOCK)
{

    // First get the bfluid of all the cells
    ARRAY u0, u1, u2, u3, bfluid0, bfluid1, bfluid2, bfluid3;
    BHAC_MHD::Getbfluid(bfluid0, bfluid1, bfluid2, bfluid3, 
                        u0, u1, u2, u3, COORDS_BLOCK, PRIMS_BLOCK);
    
    // Now iterate through the cells of the current sheet to get
    // wavevectors of all the cells
    size_t i, j, N;
    N = indices.size();

    // Clear input array
    X_K.clear();

    for (i=0; i<N; i++)
    {
        size_t i_sheet, i_upstream;
        i_sheet = indices[i][0];
        i_upstream = indices[i][1];

        // Use magnetic fields and 4 velocity of the upstream cells
        // but the position vector of the current sheet cell
        double _x, _y, _z;
        double _u0, _u1, _u2, _u3;
        double _bfluid0, _bfluid1, _bfluid2, _bfluid3;
        _x = COORDS_BLOCK[0][i_sheet];
        _y = COORDS_BLOCK[1][i_sheet];
        _z = COORDS_BLOCK[2][i_sheet];
        _u0 = u0[i_upstream];
        _u1 = u1[i_upstream];
        _u2 = u2[i_upstream];
        _u3 = u3[i_upstream];
        _bfluid0 = bfluid0[i_upstream];
        _bfluid1 = bfluid1[i_upstream];
        _bfluid2 = bfluid2[i_upstream];
        _bfluid3 = bfluid3[i_upstream];

        ARRAY X_K_local;
        Initialize_XK(X_K_local, _x, _y, _z, _u0, _u1, _u2, _u3, _bfluid0, 
                        _bfluid1, _bfluid2, _bfluid3, mode);
        X_K.push_back(X_K_local);
    }

    printvar("Number of geodesics initialized: ", X_K.size());
    return;
}