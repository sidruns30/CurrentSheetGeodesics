import numpy as np
import h5py
from scipy import interpolate
import vtk as v
import numpy_support as ah
import time

"""
    Load a vtk file using the vtk reader and convert the arrays of
    interest into a hdf5 file
"""

# Load the datacube using the default vtk reader
def LoadVtk(fname:str, array_names:list):
    # Reader to use for the simulation data
    datareader = v.vtkXMLUnstructuredGridReader()
    datareader.SetFileName(fname)
    datareader.Update()
    data = datareader.GetOutput()
    ncells = datareader.GetNumberOfCells()
    # Allocate memory for a numpy array (+ 3 from the cell centers)
    data_np = np.zeros((ncells, len(array_names) + 3))
    # Add the coordinates first
    for icell in range(ncells):
        pts=ah.vtk2array(data.GetCell(icell).GetPoints().GetData())
        # x, y and z coordinates respectively
        data_np[icell, 0] = np.average([pts[i][0] for i in range(8)])
        data_np[icell, 1] = np.average([pts[i][1] for i in range(8)])
        data_np[icell, 2] = np.average([pts[i][2] for i in range(8)])
    # Now add all the requested arrays
    for j, array in enumerate(array_names):
        print(array)
        var_arr = ah.vtk2array(data.GetCellData().GetArray(array))[0:ncells]
        data_np[:, j+3] = var_arr
    del data
    print("Successfully loaded VTK data")
    array_names = ["x", "y", "z"] + array_names
    return data_np, array_names

# Convert numpy data to hdf5 format
def ConvertToHDF5(data:np.ndarray, array_names:list):
    savefile = "BHAC_data.hdf5"
    with h5py.File(savefile, "w") as _file:
        grp = _file.create_group("data")
        for key, value in zip(array_names, data.T):
            grp.create_dataset(key, data=value)
    print("Successfully created HDF5 dataset")
    return

# Convert to text
def ConvertToText(data:np.ndarray, array_names:list):
    savefile = "BHAC_data.txt"
    ncells = data.shape[0]
    ncols = data.shape[1]
    with open(savefile, "w") as _file:
        _file.write("%d\t%d\n" % (ncells, ncols))
        for l in range(ncols):
            _file.write("%s\t" % array_names[l])
        _file.write("\n")
        for k in range(ncells):
            for l in range(ncols):
                _file.write("%f\t" % data[k,l])
            _file.write("\n")
    print("Successfully created text data")
    return



def main():
    fname = "/Users/siddhant/research/bh-flare/bhflare/analysis/data_convert1335.vtu"
    array_names = ["u1", "u2", "u3", "b1", "b2", "b3", "e1", "e2", "e3", 
                    "lfac", "p", "rho", 
                    "Bsqr", "B2", "E2",
                    "bfluid0", "bfluid1", "bfluid2", "bfluid3", 
                    "efluid0", "efluid1", "efluid2", "efluid3", "rMKS", "thetaMKS"
                    ]
    data, array_names = LoadVtk(fname, array_names)
    #ConvertToHDF5(data, array_names)
    ConvertToText(data, array_names)
    del data
    return

if __name__ == "__main__":
    main()