import numpy as np
def get_fractal():
    """
    Computes the lattice positions represented by complex numbers of the sierpinski carpet with 64 lattice points.
    :return:
    """
    offset_values = [-1, 0, 1]
    z_single = np.array([x + 1j * y for y in offset_values for x in offset_values if not (x == 0 and y == 0)]) # Draws a single square of lattice points
    z = np.array([]) # Container for lattice positions
    for y in offset_values: # Loops over all offsets
        for x in offset_values: # --||--
            if not (x == 0 and y == 0): # We dont want the middle filled out
                z = np.concatenate([z, z_single + 3 * x + 3j * y]) # Add next lattice point
    return z

def get_fractal_general(matrix, latvec, ngen):
  fullmat=np.array(matrix)
  for igen in range(ngen-1):
    fullmat2=np.zeros((fullmat.shape[0]*matrix.shape[0],fullmat.shape[0]*matrix.shape[0]), dtype=int)
    for i in range(fullmat.shape[0]):
     for j in range(fullmat.shape[1]): 
      if(fullmat[i,j]==1):
       for k in range(matrix.shape[0]):
        for l in range(matrix.shape[1]): 
         fullmat2[matrix.shape[0]*i+k, matrix.shape[1]*j+l]=matrix[k,l] 
    fullmat=np.array(fullmat2)
  zs=[]
  for i in range(fullmat.shape[0]):
   for j in range(fullmat.shape[1]):
    if(fullmat[i,j]==1):
     x=latvec[0][0]*i+latvec[1][0]*j
     y=latvec[0][1]*i+latvec[1][1]*j
     zs.append(x+1j*y)
  return zs

def get_sierpinski(ngen):
  matrix=np.zeros((2,2),dtype=int)
  matrix[0,0]=1
  matrix[0,1]=1
  matrix[1,0]=1
  latvec=[[1,0], [0.5, np.sqrt(3)/2.]]
  zs=get_fractal_general(matrix, latvec, ngen)
  return zs

def get_plane(side_length):
    """
    Computes the lattice points of a 9x9 grid represented by complex numbers
    :param side_length: Side length of grid, i.e. the grid contains side_length^2 sites
    :return:
    """
    #import numpy as np
    total_size = side_length ** 2
    low = -(side_length - 1) / 2
    high = (side_length - 1) / 2 + 1
    z = (np.arange(low, high, 1).reshape((1, side_length)) + 1j*np.arange(low, high, 1).reshape((side_length, 1))).reshape((total_size, ))
    return z


def get_rectangle(Nx, Ny):
    """
    Computes the lattice points of a 9x9 grid represented by complex numbers
    :param side_length: Side length of grid, i.e. the grid contains side_length^2 sites
    :return:
    """
    #import numpy as np
    #total_size = side_length ** 2
    #low = -(side_length - 1) / 2
    #high = (side_length - 1) / 2 + 1
    z = (np.arange(Nx).reshape((1, Nx)) + 1j*np.arange(Ny).reshape((Ny, 1))).reshape((Nx*Ny, ))
    return z

def test_lattice():
    #from matplotlib import pyplot as plt
    import numpy as np

    # ----------- Plane ----------- #
    side_length = np.random.randint(2, 20)
    z = get_plane(side_length)
    fig, ax = plt.subplots()
    ax.set_title(f'Plane with randomly chosen side length \n (side length = {side_length})')
    ax.grid()
    ax.axis('equal')
    ax.plot(np.real(z), np.imag(z), 'bo')

    # ----------- Fractal ----------- #
    z = get_fractal()
    fig, ax = plt.subplots()
    ax.set_title('Sierpinski carpet with 64 sites')
    ax.grid()
    ax.axis('equal')
    ax.plot(np.real(z), np.imag(z), 'bo')

if __name__ == '__main__':
    test_lattice()

