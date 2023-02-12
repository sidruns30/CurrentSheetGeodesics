/*
    File to copy the array data into a cartesianized grid for binning and correlation
    analysis.The grid cell size is bigger than the cell size of the simulation
*/

#include "regrid.hpp"

Cell::Cell()
{
    this->x = 0.;
    this->y = 0.;
    this->z = 0.;
    this->dx = 0.;
    this->dy = 0.;
    this->dz = 0.;
    this->data = 0.;
}

Cell::Cell(double x, double y, double z)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->dx = 0.;
    this->dy = 0.;
    this->dz = 0.;
    this->data = 0.;
}

Cell::Cell(double x, double y, double z, double dx, double dy, double dz)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->data = 0.;
}

void Cell::GetCoordinates(double &x, double &y, double &z)
{
    x = this->x;
    y = this->y;
    z = this->z;
    return;
}

void Cell::GetIndices(std::vector <size_t> &indices)
{
    indices.clear();
    for (size_t i=0; i < this->indices.size(); i++)
    {
        indices.push_back(this->indices[i]);
    }
    return;
}

void Cell::AddIndex(size_t ind)
{
    this->indices.push_back(ind);
    return;
}

void Cell::GetSize(size_t &size)
{
    size = this->indices.size();
    return;
}

double Cell::GetVolume()
{
    return this->dx * this->dy * this->dz;
}

void Cell::GetData(double &data)
{
    data = this->data;
    return;
}
// Average all the points of an array that fall inside that cell
void Cell::AverageArray(const ARRAY &arr)
{
    size_t arr_size = arr.size();
    double average = 0.;
    for (size_t i=0; i<this->indices.size(); i++)
    {
        if (this->indices[i] >= arr_size)
        {
            throw std::out_of_range("Attempted to access array element outside its size");
        }
        average += arr[this->indices[i]];
    }
    this->data = average / (this->indices.size());
    return;
}

// Return the root mean squared of an array
void Cell::RMSArray(const ARRAY &arr)
{
    size_t arr_size = arr.size();
    double rms = 0.;
    for (size_t i=0; i<this->indices.size(); i++)
    {
        if (this->indices[i] >= arr_size)
        {
            throw std::out_of_range("Attempted to access array element outside its size");
        }
        double element = arr[this->indices[i]];
        rms += element * element;
    }
    this->data = sqrt(rms / (this->indices.size()));
    return;
}

// Return the mean magnitude of the vectors
void Cell::AverageVecMagnitude(const ARRAY &arrx, const ARRAY &arry, const ARRAY &arrz)
{
    size_t arr_size = arrx.size();
    if ((arr_size != arry.size()) || (arry.size() != arrz.size()))
    {
        throw std::invalid_argument("arrays don't have equal sizes");
    }
    double averagex=0., averagey=0., averagez=0.;
    for (size_t i=0; i<this->indices.size(); i++)
    {
        if (this->indices[i] >= arr_size)
        {
            throw std::out_of_range("Attempted to access array element outside its size");
        }
        size_t index = this->indices[i];
        averagex += arrx[index];
        averagey += arry[index];
        averagez += arrz[index];
    }
    averagex /= this->indices.size();
    averagey /= this->indices.size();
    averagez /= this->indices.size();
    this->data = sqrt(SQR(averagex) + SQR(averagey) + SQR(averagez));
    return;
}

// Return the mean of the magnitude of vectors
void Cell::AverageMagnitudeVec(const ARRAY &arrx, const ARRAY &arry, const ARRAY &arrz)
{
    size_t arr_size = arrx.size();
    if ((arr_size != arry.size()) || (arry.size() != arrz.size()))
    {
        throw std::invalid_argument("arrays don't have equal sizes");
    }
    double average =0.;
    for (size_t i=0; i<this->indices.size(); i++)
    {
        if (this->indices[i] >= arr_size)
        {
            throw std::out_of_range("Attempted to access array element outside its size");
        }
        size_t index = this->indices[i];
        average += sqrt(SQR(arrx[index]) + SQR(arry[index]) + SQR(arrz[index]));
    }
    average /= this->indices.size();
    this->data = average;
    return;
}

// Get the location of the cell that contains (x, y, z)
Cell* Grid::GetCell(double x, double y, double z)
{
    size_t i, j, k;
    i = static_cast<size_t>((x - this->xmin) / this->dx);
    j = static_cast<size_t>((y - this->ymin) / this->dy);
    k = static_cast<size_t>((z - this->zmin) / this->dz);
    size_t i_cell = static_cast<size_t>((k * this->ny * this->nx) + (j * this->nx) + i);
    // Print cell coordiantes for sanity check
    double x_cell, y_cell, z_cell;
    if (i_cell > this->p_cell.size())
    {
        std::cout<<"Requested coordinates: " << x << ", " << y << ", " << z << "\n";
        std::cout<<"Requested indices: "<<i<<", "<<j<<", "<<k<<"\n";
        std::cout<<"Requested index: " << i_cell << "\n";
        std::cout<<"Grid info: \n";
        PrintGrid();
        throw std::out_of_range("cell index is outside the cell array");
    }
    return this->p_cell[i_cell];
}

Grid::Grid()
{
    this->ncells = 0;
    this->nx = 0;
    this->ny = 0;
    this->nz = 0;
    this->populated = false;
}

Grid::Grid(double xmin, double xmax, double ymin, double ymax, double zmin, 
            double zmax, double dx, double dy, double dz)
{
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->nx = static_cast<size_t>((xmax-xmin)/dx) + 1;
    this->ny = static_cast<size_t>((ymax-ymin)/dy) + 1;
    this->nz = static_cast<size_t>((zmax-zmin)/dz) + 1;
    this->ncells = this->nx * this->ny * this->nz;
    this->xarr.reserve(this->dx);
    this->yarr.reserve(this->dy);
    this->zarr.reserve(this->dz);
    size_t i;
    for (i=0; i < this->nx; i++)
    {
        this->xarr.push_back((xmin + i * dx));
    }
    for (i=0; i < this->ny; i++)
    {
        this->yarr.push_back((ymin + i * dy));
    }
    for (i=0; i < this->nz; i++)
    {
        this->zarr.push_back((zmin + i * dz));
    }
    this->MakeCells();
    this->populated = false;

}

// Make a grid from the coordinates (200 pts in each direction)
Grid::Grid(const ARRAY2D &COORDS_BLOCK, const double dl)
{
    double xmin = DBL_MAX, ymin = DBL_MAX, zmin = DBL_MAX;
    double xmax = DBL_MIN, ymax = DBL_MIN, zmax = DBL_MIN;
    for (size_t i=0; i<COORDS_BLOCK[0].size(); i++)
    {
        double x = COORDS_BLOCK[0][i];
        double y = COORDS_BLOCK[1][i];
        double z = COORDS_BLOCK[2][i];
        xmin = std::min(x, xmin);
        ymin = std::min(y, ymin);
        zmin = std::min(z, zmin);
        xmax = std::max(x, xmax);
        ymax = std::max(y, ymax);
        zmax = std::max(z, zmax);
    }
    // Make the grid slightly larger than the domain
    xmin -= dl;
    ymin -= dl;
    zmin -= dl;
    xmax += dl;
    ymax += dl;
    xmax += dl;

    //size_t cells = 200;
    double dx = dl;//(xmax - xmin) / (double)cells;
    double dy = dl;//(ymax - ymin) / (double)cells;
    double dz = dl;//(zmax - zmin) / (double)cells;

    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->nx = static_cast<size_t>((xmax-xmin)/dx) + 1;
    this->ny = static_cast<size_t>((ymax-ymin)/dy) + 1;
    this->nz = static_cast<size_t>((zmax-zmin)/dz) + 1;
    this->ncells = this->nx * this->ny * this->nz;
    this->xarr.reserve(this->dx);
    this->yarr.reserve(this->dy);
    this->zarr.reserve(this->dz);
    size_t i;
    for (i=0; i < this->nx; i++)
    {
        this->xarr.push_back((xmin + i * dx));
    }
    for (i=0; i < this->ny; i++)
    {
        this->yarr.push_back((ymin + i * dy));
    }
    for (i=0; i < this->nz; i++)
    {
        this->zarr.push_back((zmin + i * dz));
    }
    this->MakeCells();
    FillGrid(COORDS_BLOCK);
}

void Grid::MakeCells()
{
    std::cout<<"Making "<< this->ncells << " cells \n";
    this->p_cell.reserve(this->ncells);
    double dx = this->dx, dy = this->dy, dz = this->dz;
    double x, y, z;
    for (size_t k=0; k < this->nz; k++)
    {
        for (size_t j=0; j < this->ny; j++)
        {
            for (size_t i=0; i < this->nx; i++)
            {
                x = this->xarr[i];
                y = this->yarr[j];
                z = this->zarr[k];
                this->p_cell.push_back(new Cell(x, y, z, dx, dy, dz));
            }
        }
    }
    return;
}

void Grid::ClearCells()
{
    for (Cell* pcell : this->p_cell)
    {
        delete pcell;
    }
    return;
}

// Take the coordinate arrays of simulation and populate the grid
void Grid::FillGrid(const ARRAY2D &COORDS_BLOCK)
{
    for (size_t i=0; i < COORDS_BLOCK[0].size(); i++)
    {
        double x = COORDS_BLOCK[0][i];
        double y = COORDS_BLOCK[1][i];
        double z = COORDS_BLOCK[2][i];

        Cell *p_cell = this->GetCell(x, y, z);
        p_cell->AddIndex(i);
    }
    this->populated = true;
    return;
}

// Display status about the grid
void Grid::PrintGrid()
{
    std::cout<< "*** Printing Grid data ***\n";
    std::cout<< "x range (min, max, step): (" << this->xmin << ", " << this->xmax << ", " << this->dx << ")\n";
    std::cout<< "y range (min, max, step): (" << this->ymin << ", " << this->ymax << ", " << this->dy << ")\n";
    std::cout<< "z range (min, max, step): (" << this->zmin << ", " << this->zmax << ", " << this->dz << ")\n";
    std::cout<< "Dimensions: " << this->nx << " x " << this->ny << " x " << this->nz << "\n";
    std::cout <<"Size: " << this->nx * this->ny * this->nz << "; ncells: "<<  this->ncells << "\n";
    std::cout<< "Populated: " << this->populated << "\n";
    return;
}

// Compute cell values in the grid
void Grid::ComputeQuantity(std::string method, const ARRAY &arr)
{
    for (size_t i=0; i<this->ncells; i++)
    {
        Cell *pcell = this->p_cell[i];
        if (method.compare("mean"))
        {
            pcell->AverageArray(arr);
        }
        else if (method.compare("std"))
        {
            pcell->RMSArray(arr);
        }
        else
        {
            throw std::invalid_argument("Given method does not exist, try \"mean\" or \"std\"");
        }
    }
    return;
}

// For vectors
void Grid::ComputeQuantity(std::string method, const ARRAY &arrx, const ARRAY &arry, const ARRAY &arrz)
{
    for (size_t i=0; i<this->ncells; i++)
    {
        Cell *pcell = this->p_cell[i];
        if (method.compare("vec_mag_mean"))
        {
            pcell->AverageMagnitudeVec(arrx, arry, arrz);
        }
        else if (method.compare("mean_vec_mag"))
        {
            pcell->AverageVecMagnitude(arrx, arry, arrz);
        }
        else
        {
            throw std::invalid_argument("Given method does not exist, try \"vec_mag_mean\" or \"mean_vec_mag\"");
        }
    }
    return;
}

// Write the cell data from grid into numpy arrays
// Write both the average and std dev
ARRAY Grid::GetGridData()
{
    if (! this->populated)
    {
        throw std::runtime_error("Grid is not populated");
    }
    ARRAY cell_data;
    cell_data.reserve(this->ncells);
    for (size_t i=0; i<this->ncells; i++)
    {
        Cell *pcell = this->p_cell[i];
        double x, y, z, data;
        size_t size = 0;
        pcell->GetSize(size);
        //if (size >= 1)
        {
            pcell->GetCoordinates(x, y, z);
            pcell->GetData(data);
            cell_data.push_back(x);
            cell_data.push_back(y);
            cell_data.push_back(z);
            cell_data.push_back(data);
        }
    }
    return cell_data;
}

// Write the grid data to a vtk format
// Note that binary is not supported
void Grid::GridToVTK(std::string outname, bool ASCII)
{
    std::ofstream vtk_oufile;
    if (ASCII)
    {
        vtk_oufile.open(outname, std::ios::out);
    }
    else
    {
        vtk_oufile.open(outname, std::ios::out | std::ios::binary);
    }

    // FILE VERSION
    std::string VERSION("# vtk DataFile Version 4.2\n");

    // HEADER
    std::string TITLE("Grid data output\n");

    // FILE FORMAT
    std::string FILE_FORMAT;
    if (ASCII)
    {
        FILE_FORMAT = "ASCII\n";
    }
    else
    {
        FILE_FORMAT = "BINARY\n";
    }

    // DATASET STRUCTURE
    std::string GEOMETRY("DATASET RECTILINEAR_GRID\n");
    std::string _nx = std::to_string(this->nx); 
    std::string _ny = std::to_string(this->ny); 
    std::string _nz = std::to_string(this->nz); 
    std::stringstream DIMENSIONS_HEAD;
    DIMENSIONS_HEAD << "DIMENSIONS " << _nx << " " << _ny << " " << _nz << "\n";
    std::string DIMENSIONS = DIMENSIONS_HEAD.str();

    std::stringstream X_COORDINATES_HEAD;
    X_COORDINATES_HEAD << "X_COORDINATES " << std::to_string(this->nx) << " double\n";
    std::string X_COORDINATES = X_COORDINATES_HEAD.str();
    
    std::stringstream Y_COORDINATES_HEAD;
    Y_COORDINATES_HEAD << "Y_COORDINATES " << std::to_string(this->ny) << " double\n";
    std::string Y_COORDINATES = Y_COORDINATES_HEAD.str();
    
    std::stringstream Z_COORDINATES_HEAD;
    Z_COORDINATES_HEAD << "Z_COORDINATES " << std::to_string(this->nz) << " double\n";
    std::string Z_COORDINATES = Z_COORDINATES_HEAD.str();

    // DATASET ATTRIBUTES
    std::string POINT_DATA = "POINT_DATA ";
    POINT_DATA += std::to_string(this->ncells);
    POINT_DATA += "\n";
    std::string SCALARS = "SCALARS GridData double\n";
    std::string LOOKUP_TABLE = "LOOKUP_TABLE default\n";
    double temp;
    if (vtk_oufile.is_open())
    {
        vtk_oufile << VERSION;
        vtk_oufile << TITLE;
        vtk_oufile << FILE_FORMAT;
        vtk_oufile << GEOMETRY;
        vtk_oufile << DIMENSIONS;
        vtk_oufile << X_COORDINATES;
        if (ASCII) for (size_t i=0; i<this->nx; i++) vtk_oufile << this->xarr[i] << " ";
        else for (size_t i=0; i<this->nx; i++) 
        {
            vtk_oufile.write((char*)&(this->xarr[i]), sizeof(double));
            vtk_oufile << " ";
        }
        vtk_oufile << "\n";
        vtk_oufile << Y_COORDINATES;
        if (ASCII) for (size_t i=0; i<this->ny; i++) vtk_oufile << this->yarr[i] << " ";
        else for (size_t i=0; i<this->ny; i++) 
        {
            vtk_oufile.write((char*)&(this->yarr[i]), sizeof(double));
            vtk_oufile << " ";
        }
        vtk_oufile << "\n";
        vtk_oufile << Z_COORDINATES;
        if (ASCII) for (size_t i=0; i<this->nz; i++) vtk_oufile << this->zarr[i] << " ";
        else for (size_t i=0; i<this->nz; i++) 
        {
            vtk_oufile.write((char*)&(this->zarr[i]), sizeof(double));
            vtk_oufile << " ";
        }
        vtk_oufile << "\n\n";
        vtk_oufile << POINT_DATA;
        vtk_oufile << SCALARS;
        vtk_oufile << LOOKUP_TABLE;
        for (size_t i=0; i<this->ncells; i++) 
        {
            this->p_cell[i]->GetData(temp);
            if (std::isnan(temp)) temp = 0.;
            if (ASCII)  vtk_oufile << temp;
            else vtk_oufile.write((char*)&temp, sizeof(double));
            vtk_oufile << "\n";
        }
        vtk_oufile.close();
    }
}