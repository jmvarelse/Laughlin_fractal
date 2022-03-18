def main():
    import numpy as np
    import sys
    from get_lattice import get_plane, get_rectangle
    from get_psi import get_psi_naive
    from matplotlib import pyplot as plt
    from scipy.stats.mstats import trimmed_mean, trimmed_std
    import matplotlib.pyplot as plt
    print("---------------------------------------------------------------------------------------------------------------------------")    
    print("-----                                     Starting calculations                                                      ------")    
    print("---------------------------------------------------------------------------------------------------------------------------")    

    # ---------------------- Initialize system ---------------------- #
    Nx=2
    Ny=2
    #z = get_fractal()
    z= get_rectangle(Nx,Ny)
    print("z:", z)
    N = len(z)
    n_particles = 2
    q = 2
    psi, basis, basis_dict = get_psi_naive(N, n_particles, q, z)
    z_dict = dict()
    for i, x in enumerate(z):
        z_dict[x] = i
    plt.scatter(np.real(z), np.imag(z))    
    for i, x in enumerate(z):
      plt.text(x.real, x.imag,str(i))

    # ---------------------- Calculate operators ---------------------- #
    operators=all_operators(basis, basis_dict, z, z_dict, [1,1j, 1j+1, -1j+1])
    # ------------ Calculate the correct linear combination of operators ------------ #
    #d, p = find_combination(psi, operators)
    d, p = find_combination_2nd_method(psi, operators)
    imin=np.argmin(np.real(d))
    print("d=", d)
    print("min(d)", d[imin])
    print("Operator coefficients:", p[:,imin])
    o = sum([p[i, imin]*operators[i] for i in range(len(p))])
    print("Hermiticity check:", np.max(np.conj(np.transpose(o)) @ o))
    psi_T = np.transpose(np.conj(psi))
    print(psi_T@ o @ o @ psi - (psi_T @ o @ psi)**2)
    l = np.transpose(np.conj(psi)) @ o @ psi
    err = np.abs(o @ psi - l * psi)
    print(f"mean(err) = {np.mean(err)}")
    print(f"median(err) = {np.median(err)}")
    print(f"std(err) = {np.std(err)}")
    print(f"max(err) = {max(err)}")
    print(f"trimmed_mean(err) = {trimmed_mean(err)}")
    print(f"trimmed_std(err) = {trimmed_std(err)}")
    bins = np.linspace(0, trimmed_mean(err) + 3*trimmed_std(err), 100)
    
    fig, ax = plt.subplots()
    ax.grid()
    ax.hist(err, bins=bins, density=True)
    plt.show()
"""
    generator = operator_generator(basis, basis_dict, z, z_dict)
    for operators, site1, site2 in(generator):
    #generator = operator_generator2(basis, basis_dict, z, z_dict)
    #for operators,operator_descriptions, site1, site2 in(generator):
      if(len(operators)>1):
        # ------------ Calculate the correct linear combination of operators ------------ #
        print(f"site1 = {site1}, site2 = {site2}")
        d, p = find_combination(psi, operators)
        imin=np.argmin(np.real(d))
        print("d=", d)
        print("min(d)", d[imin])
        print("Operator coefficients:", p[:,imin])
        o = sum([p[i, imin]*operators[i] for i in range(len(p))])
        print("Hermiticity check:", np.max(np.conj(np.transpose(o)) @ o))
        psi_T = np.transpose(np.conj(psi))
        print(psi_T@ o @ o @ psi - (psi_T @ o @ psi)**2)
        l = np.transpose(np.conj(psi)) @ o @ psi
        err = np.abs(o @ psi - l * psi)
        print(f"mean(err) = {np.mean(err)}")
        print(f"median(err) = {np.median(err)}")
        print(f"std(err) = {np.std(err)}")
        print(f"max(err) = {max(err)}")
        print(f"trimmed_mean(err) = {trimmed_mean(err)}")
        print(f"trimmed_std(err) = {trimmed_std(err)}")
        bins = np.linspace(0, trimmed_mean(err) + 3*trimmed_std(err), 100)
        
        fig, ax = plt.subplots()
        ax.grid()
        ax.hist(err, bins=bins, density=True)
""" 
