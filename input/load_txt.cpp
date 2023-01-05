#include "load_txt.hpp"

std::vector<double> x,y,z;
std::vector<double> rho,p,u1,u2,u3,b1,b2,b3;
std::vector<double> Bsqr,bfluid0,bfluid1,bfluid2,bfluid3;

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
            p.push_back(row_data[4]);
            u1.push_back(row_data[5]);
            u2.push_back(row_data[6]);
            u3.push_back(row_data[7]);
            b1.push_back(row_data[8]);
            b2.push_back(row_data[9]);
            b3.push_back(row_data[10]);
            Bsqr.push_back(row_data[11]);
            bfluid0.push_back(row_data[12]);
            bfluid1.push_back(row_data[13]);
            bfluid2.push_back(row_data[14]);
            bfluid3.push_back(row_data[15]);
        }
        line_number ++;
    }
    int _vector_size = x.size();
    std::cout<<"Loaded the fluid data; number of cells = " << _vector_size << std::endl;
    return;
}
