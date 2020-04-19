

def read_cube(filename: str):
    """
    This function reads the wavefunction from a cube file and stores it in
    a 1D numpy array

    Args:

        filename ( string ): the name of the .cube file to read

    Returns:

        numpy.array : the 1D array of the wavefunctions for all the points on the grid
    

    Note: 

        Previously, it was called `read_1`

    """
    
    f = open(filename,'r')
    lines = f.readlines()
    f.close()

    natoms = int(lines[2].split()[0])
    # We skip a few lines in the cube files, that go as follows:
    # 2 lines - comments
    # 1 line  - the number of atoms, etc.
    # 3 lines - the grid spacing and number of grid points in each dimensions

    nstart = natoms+2+1+3        # the index of the first line containing wfc data
    nlines = len(lines)          # the total number of lines

    x = []
    for i in range(nstart,nlines):
        tmp = lines[i].split()
        ncols = len(tmp)

        for j in range(ncols):
            x.append(float(tmp[j]))
    
    data = np.array(x)
#     data = np.loadtxt(filename,skiprows=n)
    
    return data
