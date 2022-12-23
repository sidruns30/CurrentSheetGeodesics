#include "defs.hpp"

// Generic print (not sure where to place for now)
void print(std::string out)
{
    std::cout<<out<<"\n";
}

void printvar(std::string out, double var)
{
    std::cout<<out<<": "<<var<<"\n";
}
// Print 3 vector
void print3(std::string name, double vec[NDIM-1])
{
    std::cout<<name<<": "<<vec[0]<<", "<<vec[1]<<", "<<vec[2]<< "\n";
}

// Print 4 vector
void print4(std::string name, double vec[NDIM])
{
    std::cout<<name<<": "<<vec[0]<<", "<<vec[1]<<", "<<vec[2]<< ", "<<vec[3]<<"\n";
}

// Print 3 matrix
void print3M(std::string name, double mat[NDIM-1][NDIM-1])
{
    std::cout<<name<<"\n";
    int i,j;
    for(i=0;i<NDIM-1;i++)
    {
        for (j=0;j<NDIM-1; j++)
        {
            std::cout << "("<<i<<","<<j<<"): "<<mat[i][j] << "\t";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

// Print 4 matrix
void print4M(std::string name, double mat[NDIM][NDIM])
{
    std::cout<<name<<"\n";
    int i,j;
    for(i=0;i<NDIM;i++)
    {
        for (j=0;j<NDIM; j++)
        {
            std::cout << "("<<i<<","<<j<<"): "<<mat[i][j] << "\t";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
}

// Input vector referenced as data[row_id, col_id]
void WriteVectorToFile(std::string fname, std::vector <std::vector <double>> data)
{
    size_t nrows, ncols;
    nrows = data.size();
    std::cout << "Number of rows in data: "<<nrows << std::endl;
    if (nrows == 0)
    {
        std::cout<<"Nothing to write"<<std::endl;
        return;
    }
    else
    {
        ncols = data[0].size();
        std::cout << "Number of columns in data: "<<ncols << std::endl;

        size_t i, j;
        
        std::ofstream file(fname);
    
        if (file.is_open())
        {
            for (i=0; i<nrows; i++)
            {
                for (j=0; j<ncols; j++)
                {
                    file << data[i][j] << "\t";
                }
                file << std::endl;
            }
            file.close();
            std::cout << "Successfully wrote to file "<< fname << std::endl;
        }

        else 
        {
            std::cout << "Could not open file " << std::endl;
        }
    }

}

// Matrix manipulations for 3 and 4 dimensional matrices

// Code to invert matrix, taken from https://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
bool Invert4Matrix(const double m[16], double invOut[16])
{
    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}

// Taken from https://stackoverflow.com/questions/983999/simple-3x3-matrix-inverse-code-c
bool Invert3Matrix(const double m[3][3], double minv[3][3])
{
    double determinant =    m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2])
                        -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                        +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
    
    if (determinant == 0.)
    {
        return false;
    }
    double invdet = 1/determinant;
    minv[0][0] =  (m[1][1]*m[2][2]-m[2][1]*m[1][2])*invdet;
    minv[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1])*invdet;
    minv[0][2] =  (m[0][1]*m[1][2]-m[0][2]*m[1][1])*invdet;
    minv[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0])*invdet;
    minv[1][1] =  (m[0][0]*m[2][2]-m[0][2]*m[2][0])*invdet;
    minv[1][2] = -(m[0][0]*m[1][2]-m[1][0]*m[0][2])*invdet;
    minv[2][0] =  (m[1][0]*m[2][1]-m[2][0]*m[1][1])*invdet;
    minv[2][1] = -(m[0][0]*m[2][1]-m[2][0]*m[0][1])*invdet;
    minv[2][2] =  (m[0][0]*m[1][1]-m[1][0]*m[0][1])*invdet;
    return true;
}

// Returns c = a * b
void Multiply4Matrices(double a[NDIM][NDIM], double b[NDIM][NDIM], double c[NDIM][NDIM])
{
	int i, j, k;
	for (i=0; i<NDIM; i++)
	{
		for (j=0; j<NDIM; j++)
		{
			c[i][j] = 0.;
			for (k=0; k<NDIM; k++)
			{
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return;
}

void Multiply3Matrices(double a[NDIM-1][NDIM-1], double b[NDIM-1][NDIM-1], double c[NDIM-1][NDIM-1])
{
	int i, j, k;
	for (i=0; i<NDIM-1; i++)
	{
		for (j=0; j<NDIM-1; j++)
		{
			c[i][j] = 0.;
			for (k=0; k<NDIM-1; k++)
			{
				c[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return;
}

// Get the matrix transpose
void TransposeMatrix(double m[NDIM][NDIM], double minv[NDIM][NDIM])
{
    int i,j;
    for (i=0;i<NDIM;i++)
    {
        for (j=0;j<NDIM;j++)
        {
            minv[i][j] = m[j][i];
        }
    }
    return;
}
