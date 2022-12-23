#include "getk.hpp"

// Find all cells which have lie within a box_threshold of the cell
// in the current sheet

void GetClosebyCells(size_t ind, std::vector<size_t> &goodindices)
{
    double xcoord = x[ind];
    double ycoord = y[ind];
    double zcoord = z[ind];
    goodindices.clear();
    
    std::vector <size_t> _temp1, _temp2;
    size_t i;

    // search along x
    for (i=0; i<x.size(); i++)
    {
        if (fabs(x[i] - xcoord) < BOX_THRESHOLD)
        {
            _temp1.push_back(i);
        }
    }

    // search along y
    for (i=0; i<_temp1.size(); i++)
    {
        if (fabs(y[_temp1[i]] - ycoord) < BOX_THRESHOLD)
        {
            _temp2.push_back(_temp1[i]);
        }
    }

    // search along z
    for (i=0; i<_temp2.size(); i++)
    {
        if (fabs(z[_temp2[i]] - zcoord) < BOX_THRESHOLD)
        {
            goodindices.push_back(_temp2[i]);
        }
    }

    return;
}

// Get indices of cells in the current sheet
void GetCurrentSheet(std::vector <size_t> &i_sheet, bool store_data)
{
    i_sheet.clear();
    size_t i;
    for (i=0; i<x.size(); i++)
    {
        double sigma = Bsqr[i]/rho[i];
        if (sigma < SIGMA_THRESHOLD)
        {
            i_sheet.push_back(i);
        }
    }

    std::cout<<"Number of cells in the current sheet = " << i_sheet.size();

    // Write current sheet
    if (store_data)
    {
        std::vector <int> indices; 
        std::vector <double> xcpy;
        std::vector <double> ycpy;
        std::vector <double> zcpy;

        for (i=0; i<i_sheet.size(); i++)
        {
            indices.push_back(i_sheet[i]);
            xcpy.push_back(x[i_sheet[i]]);
            ycpy.push_back(y[i_sheet[i]]);
            zcpy.push_back(z[i_sheet[i]]);
        }

        std::ofstream file("fullsheet.txt");
        
        if (file.is_open())
        {
            for (i=0; i<xcpy.size(); i++)
            {
                file << indices[i] << ", " << xcpy[i] << ", "<< ycpy[i] << ", "<< zcpy[i] << "\n";
            }
            file.close();
        }
    }
    return;
}


// Write current sheet data (current sheet coordinates and b field values) to file
void WriteCurrentSheetData(std::string Outfile)
{
    std::vector <size_t> i_sheet;
    GetCurrentSheet(i_sheet, false);

    std::vector <int> indices;
    std::vector <double> xcpy;
    std::vector <double> ycpy;
    std::vector <double> zcpy;
    std::vector <double> xfldcpy;
    std::vector <double> yfldcpy;
    std::vector <double> zfldcpy;
    std::vector <double> u1cpy;
    std::vector <double> u2cpy;
    std::vector <double> u3cpy;
    std::vector <double> bsqrcpy;
    std::vector <double> bfluid0cpy;
    std::vector <double> bfluid1cpy;
    std::vector <double> bfluid2cpy;
    std::vector <double> bfluid3cpy;
    std::vector <double> lfaccpy;
    std::vector <double> B2cpy;
    std::vector <double> thetaMKScpy;
    std::vector <double> rMKScpy;
    size_t i;

    for (i=0; i<i_sheet.size(); i++)
    {
        std::vector<size_t> closecells;
        size_t k, imax;
        double bsqr_max=0.;
        std::vector<double> x_k;
        
        GetClosebyCells(i_sheet[i], closecells);
        for (k=0; k<closecells.size(); k++)
        {
            if (Bsqr[closecells[k]] > BSQR_THRESHOLD)
            {
                if (Bsqr[closecells[k]] > bsqr_max)
                {
                    bsqr_max = Bsqr[closecells[k]];
                    imax = closecells[k];
                }
            }
        }
        // imax: index of the cell in vicinty with high Bsqr (not in current sheet)
        // i_sheet[i]: index of the corresponding cell in the current sheet
        if (bsqr_max != 0.)
        {
            indices.push_back(i_sheet[i]);
            xcpy.push_back(x[i_sheet[i]]);
            ycpy.push_back(y[i_sheet[i]]);
            zcpy.push_back(z[i_sheet[i]]);
            xfldcpy.push_back(x[imax]);
            yfldcpy.push_back(y[imax]);
            zfldcpy.push_back(z[imax]);
            u1cpy.push_back(u1[imax]);
            u2cpy.push_back(u2[imax]);
            u3cpy.push_back(u3[imax]);
            bsqrcpy.push_back(Bsqr[imax]);
            bfluid0cpy.push_back(bfluid0[imax]);
            bfluid1cpy.push_back(bfluid1[imax]);
            bfluid2cpy.push_back(bfluid2[imax]);
            bfluid3cpy.push_back(bfluid3[imax]);
            lfaccpy.push_back(lfac[imax]);
            B2cpy.push_back(B2[imax]);
            thetaMKScpy.push_back(thetaMKS[imax]);
            rMKScpy.push_back(rMKS[imax]);
        }
    }

    std::cout<<"Number of cells in the outer sheet: "<<xcpy.size()<<std::endl;
    std::ofstream file(Outfile);
    
    if (file.is_open())
    {
        for (i=0; i<xcpy.size(); i++)
        {
            file << indices[i] << ", " << xcpy[i] << ", "<< ycpy[i] << ", "<< zcpy[i] << ", "
                << xfldcpy[i] << ", " << yfldcpy[i] << ", " << zfldcpy[i] << ", " <<
                u1cpy[i] << ", " << u2cpy[i] << ", " << u3cpy[i] << ", " << bsqrcpy[i] << 
                ", " << bfluid0cpy[i] << ", " << bfluid1cpy[i] << ", " << bfluid2cpy[i] <<
                ", " << bfluid3cpy[i] << ", "<< lfaccpy[i]<< ", "<< B2cpy[i] << 
                ", "<<thetaMKScpy[i]<<", "<<rMKScpy[i]<< "\n";
        }
        file.close();
    }

    return;
}

void ReadSheetDataFromFile (std::string FILE_NAME, std::vector <double> data[NVARS_SHEET])
{
    std::ifstream infile(FILE_NAME);
    std::string line;
    int i, line_number = 0;
    // Stored arrays are: ind, x, y, z, xs, ys, zs, u1, u2, u3, bsq, bfld0, bfld1, bfld2, bfld3, lfac, B2
    double row_data[NVARS_SHEET];
    for (i=0; i<NVARS_SHEET; i++)
    {
        data[i].clear();
    }
    while (std::getline(infile, line))
    {
        std::istringstream row(line);
        std::string row_s;
        for (i=0; i<NVARS_SHEET; i++)
        {
            row >> row_s;
            row_data[i] = std::stod(row_s);
            data[i].push_back(row_data[i]);
        }
    }
    return;
}


// We are doing everything in MKS (i.e. metric_type = 3)
std::vector <double> Initialize_XK( int metric_type, double _x, double _y, double _z, double _u1, double _u2, double _u3,
                                    double _bfluid0, double _bfluid1, double _bfluid2, double _bfluid3,
                                    double _lfac, int mode)
{
    // mode = {0, 1, 2} for parallel, antiparallel and perpendicular wavevectors
    int i,j;

    // Get MKS coordinates
    double XKS[NDIM], XMKS[NDIM];
    CartToKS(_x, _y, _z, XKS);
    TransformCoordinates(2, 3, XKS, XMKS);

    // Get both upper and lower metrics
    double gcov[NDIM][NDIM];
    double gcon[NDIM][NDIM];

    GcovFunc(metric_type, XMKS, gcov);
    GconFunc(metric_type, XMKS, gcon);

    // Get upper and lower lapse and shift vectors
    double  alpha, beta_con[NDIM-1], beta_cov[NDIM-1], 
            gamma_con[NDIM-1][NDIM-1], gamma_cov[NDIM-1][NDIM-1];

    alpha = sqrt(-1/gcon[0][0]);
    for (i=0;i<NDIM-1;i++)
    {
        beta_con[i] = gcon[0][i+1] * SQR(alpha);
        beta_cov[i] = gcov[0][i+1];
    }
    for (i=0; i<NDIM-1; i++)
    {
        for (j=0; j<NDIM-1; j++)
        {
            gamma_con[i][j] = gcon[i+1][j+1] + beta_con[i] * beta_con[j] / SQR(alpha);
            gamma_cov[i][j] = gcov[i+1][j+1];
        }
    }

    // Get 3 velocity in MKS
    double ui_input_Cart[NDIM-1] = {_u1, _u2, _u3}, ui_input_MKS[NDIM-1];
    T_3CartTo3MKS(ui_input_Cart, ui_input_MKS, XMKS);

    // The contravariant magnetic field
    double b_con[NDIM] = {_bfluid0, _bfluid1, _bfluid2, _bfluid3};
    double ui_con[NDIM-1];
    double  U_con_MKS[NDIM];

    // Input is lfac v^i
    for(i=0;i<NDIM-1;i++)
    {
        ui_con[i] = ui_input_MKS[i] - _lfac * beta_con[i]/alpha;
    }

    bool success = false;
    success = Get4Velocity(metric_type, XMKS, ui_con, U_con_MKS, b_con);
    if (!success)
    {
        std::vector <double> X_k = {0.};
        return X_k;
    }
    
    //Test4Velocity(3, 3, XMKS, U_con_MKS, b_con, _lfac, alpha);

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

    // Write to full vector + tetrad stuff
    std::vector <double> X_k = {XMKS[0], XMKS[1], XMKS[2], XMKS[3], k_con_MKS[0], k_con_MKS[1], k_con_MKS[2], k_con_MKS[3]};
    return X_k;
}



// Write the wavevectors to file: give the current sheet arrays and the b field data
void WriteX_kToFile(std::string Infile, std::string OutName)
{
    // Read the existing data from file
    std::vector <double> all_data[NVARS_SHEET];
    ReadSheetDataFromFile (Infile, all_data);

    std::vector <double> x = all_data[1];
    std::vector <double> y = all_data[2];
    std::vector <double> z = all_data[3];    
    std::vector <double> x_f = all_data[4];
    std::vector <double> y_f = all_data[5];
    std::vector <double> z_f = all_data[6];
    std::vector <double> u1 = all_data[7];
    std::vector <double> u2 = all_data[8];
    std::vector <double> u3 = all_data[9];
    std::vector <double> Bsqr = all_data[10];
    std::vector <double> bfluid0 = all_data[11];
    std::vector <double> bfluid1 = all_data[12];
    std::vector <double> bfluid2 = all_data[13];
    std::vector <double> bfluid3 = all_data[14];
    std::vector <double> lfac = all_data[15];
    std::vector <double> B2 = all_data[16];
    std::vector <double> thetaMKS = all_data[17];
    std::vector <double> rMKS = all_data[18];

    std::cout<<"Data File Read: number of elements = "<<x.size()<<std::endl;

    size_t ind, i=0, j=0;
    std::vector <std::vector <double>> data;
    std::vector <double> X_k;
    
    for (ind=0; ind<x.size();ind++)
    {
        // Mks metric
        int metric_type = 3;
        // Parallel wavevectors
        int mode = 0;
        X_k = Initialize_XK(metric_type, x_f[ind], y_f[ind], z_f[ind], u1[ind], u2[ind], u3[ind],
                            bfluid0[ind], bfluid1[ind], bfluid2[ind], bfluid3[ind], lfac[ind], mode);

        if (X_k.size() == 2*NDIM)
        {
            // write the position and wave vector
            data.push_back(X_k);
            i++;
        }

        // Anti-parallel wave vectors
        mode = 1;
        X_k = Initialize_XK(metric_type, x_f[ind], y_f[ind], z_f[ind], u1[ind], u2[ind], u3[ind],
                            bfluid0[ind], bfluid1[ind], bfluid2[ind], bfluid3[ind], lfac[ind], mode);

        if (X_k.size() == 2*NDIM)
        {
            data.push_back(X_k);
            i++;
        }
    }

    std::cout<<"Number of Geodesics Initialized = "<<data.size()<<std::endl;
    WriteVectorToFile(OutName, data);

    return;
}