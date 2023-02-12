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

"""
    Convert numpy data to hdf5 format
"""
def ConvertToHDF5(data:np.ndarray, array_names:list):
    savefile = "BHAC_data.hdf5"
    with h5py.File(savefile, "w") as _file:
        grp = _file.create_group("data")
        for key, value in zip(array_names, data.T):
            grp.create_dataset(key, data=value)
    print("Successfully created HDF5 dataset")
    return

"""
    Convert to text
"""
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

"""
    Convert vtk to vtu
    Takes multiple fnames for each array and a dir for where
    the files are located
    Option to append numpy arrays
"""
def vtkTovtu(fnames, dir, outname="test.vtu"):
    print("Converting vtk to vtu...")
    import meshio

    point_data = {}
    print("Writing vtk data")
    for file in fnames:
        loc = dir + file
        mesh = meshio.read(loc, file_format="vtk")
        
        # append all the arrays into one 
        for key, value in mesh.point_data.items():
            if value.shape[-1] == 3:
                if file == "Bcart_reduced212_9537n0.vtk":
                    point_data["b1"] = value[:,0]
                    point_data["b2"] = value[:,1]
                    point_data["b3"] = value[:,2]
                elif file == "vcart_reduced212_9537n0.vtk":
                    point_data["u1"] = value[:,0]
                    point_data["u2"] = value[:,1]
                    point_data["u3"] = value[:,2]
            else: 
                point_data[key] = value

     # create a mesh object
    points = mesh.points
    cells = mesh.cells # suspicious -- do all arrays have the same cells?
    del mesh
    mesh = meshio.Mesh(points=points, cells=cells, point_data=point_data)
    #mesh.write(outname)
    print("done")
    print("cells are", cells)
    return cells

"""
    Convert numpy arrays into vtu
"""
def npzTovtu(dir, arrnames, cells, outname="test_calc.vtu"):
    print("Converting npy data to vtu")
    # Firt need to rearrange cells based on sorted elements
    point_data = {}
    points = [0, 0, 0]
    # First append existing numpy arrays if present
    for arr in arrnames:
        loc = dir + arr + ".npy"
        data = np.load(loc)
        print("length of ", arr, " is ", len(data))
        if arr == "x" :
            points[0] = data
        elif arr == "y":
            points[1] = data
        elif arr == "z":
            points[2] = data
        else:
            point_data[arr] = data
    import meshio
    points = np.array(points, dtype=float).T
    mesh = meshio.Mesh(points=points, point_data=point_data, cells=cells)
    mesh.write(outname)
    del mesh
    return

"""
    Covert vtk to npz arrays
"""
def vtkTonpz(fnames, dir, outdir):
    print("Converting vtk to npz")
    import meshio

    x = np.ndarray([])
    y = np.ndarray([])
    z = np.ndarray([])
    p = np.ndarray([])
    rho = np.ndarray([])
    u1 = np.ndarray([])
    u2 = np.ndarray([])
    u3 = np.ndarray([])
    b1 = np.ndarray([])
    b2 = np.ndarray([])
    b3 = np.ndarray([])

    for file in fnames:
        loc = dir + file
        mesh = meshio.read(loc, file_format="vtk")
        # append arrays into npz arrays
        for key, value in mesh.point_data.items():
            if file == "Bcart_reduced212_9537n0.vtk":
                b1 = value[:,0].reshape(-1)
                b2 = value[:,1].reshape(-1)
                b3 = value[:,2].reshape(-1)
            elif file == "vcart_reduced212_9537n0.vtk":
                u1 = value[:,0].reshape(-1)
                u2 = value[:,1].reshape(-1)
                u3 = value[:,2].reshape(-1)
            elif file == "p_reduced212_9537n0.vtk":
                p = value[:].reshape(-1)
            elif file == "Rho_reduced212_91_9537n0.vtk":
                rho = value[:].reshape(-1)
        
    x = mesh.points[:,0]
    y = mesh.points[:,1]
    z = mesh.points[:,2]

    # Set the same data type before saving
    x.astype(np.float64)
    y.astype(np.float64)
    z.astype(np.float64)
    u1.astype(np.float64)
    u2.astype(np.float64)
    u3.astype(np.float64)
    b1.astype(np.float64)
    b2.astype(np.float64)
    b3.astype(np.float64)
    p.astype(np.float64)
    rho.astype(np.float64)


    np.save(outdir + "/x.npy", x)
    np.save(outdir + "/y.npy", y)
    np.save(outdir + "/z.npy", z)
    np.save(outdir + "/rho.npy", rho)
    np.save(outdir + "/p.npy", p)
    np.save(outdir + "/u1.npy", u1)
    np.save(outdir + "/u2.npy", u2)
    np.save(outdir + "/u3.npy", u3)
    np.save(outdir + "/b1.npy", b1)
    np.save(outdir + "/b2.npy", b2)
    np.save(outdir + "/b3.npy", b3)

    return

"""
    Make an array of zeros if a certain prim vtk
    array is not available
"""
def makezeros(goodarray, loc):
    # open good array to get the size of the array
    shape = np.load(goodarray).shape
    p = np.zeros(shape)
    np.save(loc, p)
    return

def main():
    # first load the vtk and make a vtu
    datadir = "/Users/siddhant/research/bh-flare/data/highresdata/"
    fnames =    ["Rho_reduced212_91_9537n0.vtk"]#, "vcart_reduced212_9537n0.vtk", 
                #"Bcart_reduced212_9537n0.vtk"]
            # "p_reduced212_9537n0.vtk", "Bcart_reduced212_9537n0.vtk", 
    #vtkTonpz(fnames=fnames, dir=datadir, outdir="npy_data")
    #makezeros("npy_data/rho.npy", "npy_data/p.npy")
    cells = vtkTovtu(fnames, datadir, "test.vtu")
    npzTovtu(dir="", arrnames=["b2", "x", "y", "z", "u0", "rho", "beta", "bsqr", "sigma"], cells=cells)
    #fname = "/Users/siddhant/research/bh-flare/data/data_convert0303.vtu"
    #outname="test.vtu"
    #array_names = ["u1", "u2", "u3", "b1", "b2", "b3", 
    #                "p", "rho"]
    #data, array_names = LoadVtk(outname, array_names)
    #ConvertToHDF5(data, array_names)
    #ConvertToText(data, array_names)
    #del data
    return

if __name__ == "__main__":
    main()