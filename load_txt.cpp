#include "load_txt.hpp"

// Global fields
ARRAY Bsqr_sim,bfluid0_sim,bfluid1_sim,bfluid2_sim,bfluid3_sim;
ARRAY efluid0_sim,efluid1_sim,efluid2_sim,efluid3_sim;
ARRAY b2_sim, e2_sim, rMKS_sim, thetaMKS_sim;

ARRAY2D COORDS;
ARRAY2D PRIMS;

void InitializeArrays(ARRAY2D &COORDS, ARRAY2D &PRIMS, 
                    std::string FILE_NAME)
{
    std::ifstream infile(FILE_NAME);
    std::string line;
    int i, line_number = 0;
    double row_data[NVARS];

    COORDS.clear();
    PRIMS.clear();

    ARRAY x,y,z;
    ARRAY rho,P,V1,V2,V3,B1,B2,B3,lfac,E1,E2,E3;

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
            
            // Coordinates and prims
            x.push_back(row_data[0]);
            y.push_back(row_data[1]);
            z.push_back(row_data[2]);
            V1.push_back(row_data[3]);
            V2.push_back(row_data[4]);
            V3.push_back(row_data[5]);
            B1.push_back(row_data[6]);
            B2.push_back(row_data[7]);
            B3.push_back(row_data[8]);
            E1.push_back(row_data[9]);
            E2.push_back(row_data[10]);
            E3.push_back(row_data[11]);
            lfac.push_back(row_data[12]);
            P.push_back(row_data[13]);
            rho.push_back(row_data[14]);
            // Extra arrays
            Bsqr_sim.push_back(row_data[15]);
            b2_sim.push_back(row_data[16]);
            e2_sim.push_back(row_data[17]);
            bfluid0_sim.push_back(row_data[18]);
            bfluid1_sim.push_back(row_data[19]);
            bfluid2_sim.push_back(row_data[20]);
            bfluid3_sim.push_back(row_data[21]);
            efluid0_sim.push_back(row_data[22]);
            efluid1_sim.push_back(row_data[23]);
            efluid2_sim.push_back(row_data[24]);
            efluid3_sim.push_back(row_data[25]);
            rMKS_sim.push_back(row_data[26]);
            thetaMKS_sim.push_back(row_data[27]);
            // Normally read fields, set to 0 for now
            if (idealMHD)
            {
                // Do something with the fields
            }
        }
        line_number ++;
    }

    // Push back data to global arrays
    COORDS.push_back(x);
    COORDS.push_back(y);
    COORDS.push_back(z);

    PRIMS.push_back(V1);
    PRIMS.push_back(V2);
    PRIMS.push_back(V3);
    PRIMS.push_back(B1);
    PRIMS.push_back(B2);
    PRIMS.push_back(B3);
    PRIMS.push_back(E1);
    PRIMS.push_back(E2);
    PRIMS.push_back(E3);
    PRIMS.push_back(lfac);
    PRIMS.push_back(P);
    PRIMS.push_back(rho);
       
    int _vector_size = x.size();
    std::cout<<"Loaded the fluid data; number of cells = " << _vector_size << std::endl;
    return;
}
