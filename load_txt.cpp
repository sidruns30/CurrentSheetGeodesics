#include "load_txt.hpp"

void InitializeArrays(std::string FILE_NAME)
{
    std::ifstream infile(FILE_NAME);
    std::string line;
    int i, line_number = 0;
    double row_data[NVARS];

    while (std::getline(infile, line))
    {
        // The first two lines form the user header
        if (line_number >=2)
        {
            std::istringstream row(line);
            for (i=0; i<NVARS; i++)
            {
                row >> row_data[i];
            }
            
            x.push_back(row_data[0]);
            y.push_back(row_data[1]);
            z.push_back(row_data[2]);
            rho.push_back(row_data[3]);
            P.push_back(row_data[4]);
            V1.push_back(row_data[5]);
            V2.push_back(row_data[6]);
            V3.push_back(row_data[7]);
            B1.push_back(row_data[8]);
            B2.push_back(row_data[9]);
            B3.push_back(row_data[10]);
            Bsqr_sim.push_back(row_data[11]);
            bfluid0_sim.push_back(row_data[12]);
            bfluid1_sim.push_back(row_data[13]);
            bfluid2_sim.push_back(row_data[14]);
            bfluid3_sim.push_back(row_data[15]);
            lfac.push_back(row_data[16]);
            b2_sim.push_back(row_data[17]);
            thetaMKS_sim.push_back(row_data[18]);
            rMKS_sim.push_back(row_data[19]);
            // Normally read fields, set to 0 for now
            if (!idealMHD)
            {
                E1.push_back(0.);
                E2.push_back(0.);
                E3.push_back(0.);
            }
        }
        line_number ++;
    }
    int _vector_size = x.size();
    std::cout<<"Loaded the fluid data; number of cells = " << _vector_size << std::endl;
    return;
}
