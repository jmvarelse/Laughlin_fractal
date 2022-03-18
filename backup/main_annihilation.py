import numpy as np
import sys
from get_lattice import get_plane, get_rectangle, get_sierpinski
from get_psi import get_psi_naive
from matplotlib import pyplot as plt
from scipy.stats.mstats import trimmed_mean, trimmed_std
from scipy.linalg import orth
import matplotlib.pyplot as plt
#import operafperators
from operators import *
from itertools import combinations
import pickle as pkl
print("---------------------------------------------------------------------------------------------------------------------------")    
print("-----                                     Starting calculations                                                      ------")    
print("---------------------------------------------------------------------------------------------------------------------------")    

# ---------------------- Initialize system ---------------------- #
Nx=3
Ny=3
#z = get_fractal()
#z= get_rectangle(Nx,Ny)
z=get_sierpinski(3)
plt.scatter(np.real(z), np.imag(z))
plt.show()
print("z:", z)
N = len(z)
n_particles = 6
q = 2
psi, basis, basis_dict = get_psi_naive(N, n_particles, q, z)
for i, coeff in enumerate(psi):
  print("i=", i, "conf=", basis[i], "psi[", i,"]=",psi[i])
basis2 = list(combinations(range(N), n_particles-1)) # Compute all basis states with n_partices on N sites
basis_dict2 = dict() # Initialize a dictionary for translating a given basis state into the index in the reduced basis
for i, b in enumerate(basis2):
    basis_dict2[b] = i # Loop through all basis states and assign each basis state its index in the reduced basis.
    print("i=", i, "conf=", b)

z_dict = dict()
for i, x in enumerate(z):
    z_dict[x] = i
plt.scatter(np.real(z), np.imag(z))    
for i, x in enumerate(z):
  plt.text(x.real, x.imag,str(i))

thres=0.005

operators_all=[]
r=1.1
for isite in range(N):
  operator_descriptions=[]
  for i in range(N):
    if(abs(z[isite]-z[i])<r):
      operator_descriptions.append(["annihilation", i, 1])
  #for i in range(N):
  #  if(abs(z[isite]-z[i])<r):
  #    if(isite!=i):
  #      operator_descriptions.append(["an", isite, i, 1])

  operators=operators_from_descriptions(operator_descriptions, [basis, basis_dict, basis2, basis_dict2], 1)

  #operators=[]
  #operator_descriptions=[]
  #[operators, operator_descriptions]=add_operator2("annihilation", [0], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("annihilation", [1], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("annihilation", [2], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("annihilation", [3], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("an", [0,1], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("an", [0,2], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #[operators, operator_descriptions]=add_operator2("an", [0,3], basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions)
  #print(operators)
  #print(operators2)
  #ref=[-3., -1., -1., -1.,  2.,  2.,  2.]
  #ref=ref/np.sqrt(np.dot(np.conj(ref), ref))
  #print("ref",ref)
  d, p = find_combination_2nd_method(psi, operators)
  print("d=", d)
  nullspace=[]
  for di, d0 in enumerate(d):
    if(abs(d0)<thres):
      nullspace.append(p[:,di])
  print("Matrix size:",len(d), "nullspace size:", len(nullspace))
  for inull, null in enumerate(nullspace):
    operator_descriptions2=update_coeffs(operator_descriptions, null)
    operators_all.append(operator_descriptions2)     
    #test_operator=0
    #for icoeff,coeff in enumerate(null):
    #  test_operator=test_operator+operators[icoeff]*coeff    
    #print("Test:")
    #print(test_operator @ psi)


f=open("operators.pkl", "wb")
pkl.dump(operators_all, f)
f.close()
print(*operators_all, sep="\n")

