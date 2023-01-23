#include "load_txt.hpp"
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "../cnpy/cnpy.h"

// Global fields
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
    ARRAY rho,P,V1,V2,V3,B1,B2,B3,E1,E2,E3;

    // Count thr number of lines

    while (std::getline(infile, line))
    {
        if (line_number == 0)
        {
            std::istringstream row(line);
            size_t size;
            row >> size;
            std::cout << "Size of file is " << size << std::endl;
            x.reserve(size);
            y.reserve(size);
            z.reserve(size);
            rho.reserve(size);
            P.reserve(size);
            V1.reserve(size);
            V2.reserve(size);
            V3.reserve(size);
            B1.reserve(size);
            B2.reserve(size);
            B3.reserve(size);
            E1.reserve(size);
            E2.reserve(size);
            E3.reserve(size);
        }
    
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
            P.push_back(row_data[12]);
            rho.push_back(row_data[13]);
            // Extra arrays
        }
        line_number ++;
    }

    // Push back data to global arrays
    COORDS.push_back(x);
    COORDS.push_back(y);
    COORDS.push_back(z);

    PRIMS.push_back(rho);
    PRIMS.push_back(P);
    PRIMS.push_back(V1);
    PRIMS.push_back(V2);
    PRIMS.push_back(V3);
    PRIMS.push_back(B1);
    PRIMS.push_back(B2);
    PRIMS.push_back(B3);
    PRIMS.push_back(E1);
    PRIMS.push_back(E2);
    PRIMS.push_back(E3);


    int _vector_size = x.size();
    std::cout<<"Loaded the fluid data; number of cells = " << _vector_size << std::endl;
    return;
}

// Void initialize single numpy array and store it in a vector
void LoadNumpyArray(std::string fname, std::string arrname, ARRAY &vec)
{
    std::cout<<"Loading array "<<arrname<<std::endl;
    cnpy::NpyArray *parray = new cnpy::NpyArray;
    *parray = cnpy::npz_load(fname, arrname);
    size_t i, nelements = parray->shape[0];
    vec.reserve(nelements);
    double *data = parray->data<double>();
    for (i=0; i<nelements; i++)
    {
        vec.push_back(data[i]);
    }
    delete parray;
    return;
}


// Initialize numpy arrays
void InitializeNumpyArrays(ARRAY2D &COORDS, ARRAY2D &PRIMS, std::string fname)
{
    ARRAY x, y, z;
    ARRAY rho, p, u1, u2, u3, b1, b2, b3;

    LoadNumpyArray(fname, "x", x);
    LoadNumpyArray(fname, "y", y);
    LoadNumpyArray(fname, "z", z);
    LoadNumpyArray(fname, "u1", u1);
    LoadNumpyArray(fname, "u2", u2);
    LoadNumpyArray(fname, "u3", u3);
    LoadNumpyArray(fname, "b1", b1);
    LoadNumpyArray(fname, "b2", b2);
    LoadNumpyArray(fname, "b3", b3);
    LoadNumpyArray(fname, "rho", rho);
    LoadNumpyArray(fname, "p", p);

    // Push back data to global arrays    
    COORDS.push_back(x);
    COORDS.push_back(y);
    COORDS.push_back(z);
    PRIMS.push_back(rho);
    PRIMS.push_back(p);  
    PRIMS.push_back(u1);
    PRIMS.push_back(u2);
    PRIMS.push_back(u3);
    PRIMS.push_back(b1);
    PRIMS.push_back(b2);
    PRIMS.push_back(b3);

    
    x.clear();
    y.clear();
    z.clear();
    rho.clear();
    p.clear();
    u1.clear();
    u2.clear();
    u3.clear();
    b1.clear();
    b2.clear();
    b3.clear();
    
    return;
}

/*
// Memory mapped reader for faster io
void InitializeArraysMmap(std::string fname)
{
    // Initial vectors
    ARRAY x,y,z;
    ARRAY rho,P,V1,V2,V3,B1,B2,B3,lfac,E1,E2,E3;

    // Convert to char *
    char* fname_c = const_cast<char*>(fname.c_str());

    int fd = open(fname_c, O_RDONLY, S_IRUSR);
    struct stat sb;

    if (fstat(fd, &sb) == -1)
    {
        print("Could not read file");
    }

    const char* file_in_memory = static_cast<const char*>(mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));

    int i, varctr=0, linectr=0, skip_lines=NVARS+2, tructr=0;
    std::string number;

    for (i=0; i<2000; i++)
    {
        if (file_in_memory[i] != '\n' && file_in_memory[i] != '\t')
        {
            number.push_back(file_in_memory[i]);
        }
        else
        {
            linectr ++;

            if (linectr >= skip_lines && number.length() > 1)
            {
                    varctr = (linectr - (NVARS+2)) % NVARS;
                    printvar("num", number);
                    printvar("line number", linectr);
                    printvar("var number", varctr);
                
                number.clear();
            }

        }


        
         (number.length() >= 1)
        {
            if (linectr >= NVARS + 2)
            {
                varctr = ((linectr - (NVARS + 2)) % NVARS + 1);
                printvar("num", number);
                number.clear();
                printvar("line number", linectr);
                printvar("var number", varctr);
            }
            linectr++;

        }
        
    }
    return;
}
*/