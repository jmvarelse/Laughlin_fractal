import numpy as np
import sys
from get_lattice import get_plane, get_rectangle
from get_psi import get_psi_naive
from matplotlib import pyplot as plt
from scipy.stats.mstats import trimmed_mean, trimmed_std
import matplotlib.pyplot as plt
from itertools import combinations
#import operafperators
from operators import *
print("---------------------------------------------------------------------------------------------------------------------------")    
print("-----                                     Starting calculations                                                      ------")    
print("---------------------------------------------------------------------------------------------------------------------------")    

# ---------------------- Initialize system ---------------------- #
Nx=2
Ny=3
#z = get_fractal()
z= get_rectangle(Nx,Ny)
print("z:", z)
N = len(z)
n_particles = 3
q = 2
psi, basis, basis_dict = get_psi_naive(N, n_particles, q, z)
for i, coeff in enumerate(psi):
  print("i=", i, "conf=", basis[i], "psi[", i,"]=",psi[i])

print("Second basis:")

basis2 = list(combinations(range(N), n_particles-1)) # Compute all basis states with n_partices on N sites
basis_dict2 = dict() # Initialize a dictionary for translating a given basis state into the index in the reduced basis
for i, b in enumerate(basis2):
    basis_dict2[b] = i # Loop through all basis states and assign each basis state its index in the reduced basis.
    print("i=", i, "conf=", b)


i=0

operator=(q-2)*annihilation_operator(basis, basis2, basis_dict2, i)
for j in range(N):
 if(j!=i):
   omega=(z[i]+z[j])/(z[i]-z[j])
   operator=operator+omega*(annihilation_operator(basis, basis2, basis_dict2, j)-(q*annihilation_operator(basis, basis2, basis_dict2, i) @ number_operator(basis, j))+annihilation_operator(basis, basis2, basis_dict2, i))

coeffs0=np.zeros(2*N, dtype=complex)
coeffs0[i]=q-2
for j in range(N):
 if(j!=i):
  omega=(z[i]+z[j])/(z[i]-z[j])
  #omega=1./(z[i]-z[j])
  coeffs0[i]=coeffs0[i]+omega 
  coeffs0[j]=omega
  coeffs0[j+N]=-omega*q
coeffs0=coeffs0/np.sqrt(np.dot(np.conj(coeffs0), coeffs0))
print("coeffs", coeffs0)


coeffs1=np.array([ 0.56224144+0.j,  0.38894745+0.j,  0.51186138+0.j,  0.45040442+0.j,
  0,0.0110769 +0.j, -0.23475096+0.j, -0.11183703+0.j])
coeffs1=coeffs1/np.sqrt(np.dot(np.conj(coeffs1), coeffs1))
coeffs2=np.array([-0.17285295+0.j,  0.14006809+0.j,  0.39825395+0.j,  0.26916102+0.j,
 0, 0.70019983+0.j,  0.18382811+0.j,  0.44201397+0.j])
coeffs2=coeffs2/np.sqrt(np.dot(np.conj(coeffs2), coeffs2))
coeffs3=np.array([-0.00316801+0.j,  0.44831697+0.j, -0.13453394+0.j,  0.15689151+0.j,
 0, -0.42279138+0.j,  0.74291044+0.j,  0.16005953+0.j])
coeffs3=coeffs3/np.sqrt(np.dot(np.conj(coeffs3), coeffs3))

"""
print(abs(np.dot(coeffs1, coeffs0))**2, abs(np.dot(coeffs2, coeffs0))**2, abs(np.dot(coeffs3, coeffs0))**2)
print(abs(np.dot(coeffs1, coeffs0))**2+abs(np.dot(coeffs2, coeffs0))**2+abs(np.dot(coeffs3, coeffs0))**2)
"""

operator2=coeffs0[i]*annihilation_operator(basis, basis2, basis_dict2, i)


for j in range(N):
 if(j!=i):
  operator2=operator2+coeffs0[j]*annihilation_operator(basis, basis2, basis_dict2, j)
  operator2=operator2+coeffs0[j+N]*annihilation_operator(basis, basis2, basis_dict2, i) @ number_operator(basis, j)


#print(annihilation_operator(basis, basis2, basis_dict2, i) @ number_operator(basis, 1) @ psi)
print(operator @ psi)
print(operator2 @ psi)

plt.show()
