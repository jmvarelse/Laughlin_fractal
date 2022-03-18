import numpy as np
import sys
from get_lattice import get_plane, get_rectangle
from get_psi import get_psi_naive
from matplotlib import pyplot as plt
from scipy.stats.mstats import trimmed_mean, trimmed_std
from scipy.linalg import orth
import matplotlib.pyplot as plt
#import operafperators
from operators import *
from itertools import combinations
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

operator_descriptions=[
   ["annihilation", 0, 1],
   ["annihilation", 1, 1],
   ["annihilation", 2, 1],
   ["annihilation", 3, 1],
   ["an", 0, 1, 1],
   ["an", 0, 2, 1],
   ["an", 0, 3, 1],
]



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

ref=[-3., -1., -1., -1.,  2.,  2.,  2.]
ref=ref/np.sqrt(np.dot(np.conj(ref), ref))
print("ref",ref)


d, p = find_combination_2nd_method(psi, operators)
print("d=", d)
arg=np.argsort(d)
nullspace=[]
for i, a in enumerate(arg[:3]):
  nullspace.append(p[:,a])
  print("Operator coefficients:", p[:,a])
overlap=0j
for vec in nullspace:
  overlap=overlap+np.abs(np.dot(np.conj(ref), vec))**2

results=[]
  
oref = np.array(sum([ref[j]*operators[j] for j in range(len(ref))]).todense())

oref2=oref.flatten()

print(oref2.shape)
#oref2=oref2/np.dot(np.transpose(np.conj(oref2)), oref2)
oref2=np.array(oref2)/np.sqrt(np.dot(np.conj(oref2), oref2))
print(oref2)
#print(np.dot(np.transpose(np.conj(oref2)), oref2) )


nullspace2=[]
for i, a in enumerate(arg[:3]):
  o = sum([p[j, a]*operators[j] for j in range(len(p[:,a]))])
  o2=np.array(o.todense()).flatten()
  o2=o2/np.sqrt(np.dot(np.conj(o2), o2))
  nullspace2.append(o2)

#print(np.array(nullspace2).shape())


nullspace_basis=orth(np.transpose(np.array(nullspace2)))
print(nullspace_basis)


overlap2=0
for vec in np.transpose(nullspace_basis):
  print("vec", vec)
  print("overlap", np.abs(np.dot(np.conj(oref2), vec))**2)
  overlap2=overlap2+np.abs(np.dot(np.conj(oref2), vec))**2 
   



print("overlap:", overlap, overlap2)

#imin=np.argmin(np.real(d))
#print("min(d)", d[imin])
#print("Operator coefficients:", p[:,imin])


Hamiltonian=[]
# construct the Hamiltonian from the annihilation operators
for i, o1 in enumerate(operator_descriptions):
 for j, o2 in enumerate(operator_descriptions):
  Hamiltonian.append(combine_operators(o1, o2, np.conj(p[i,arg[0]])*p[j,arg[0]]))
for h in Hamiltonian:
 print(h)


"""
Hamiltonian=[]
# construct the Hamiltonian from the annihilation operators
for i, o1 in enumerate(operator_descriptions):
 for j, o2 in enumerate(operator_descriptions):
  if(o1[0]=="annihilation" and o2[0]=="annihilation"):
    site1=o1[1]
    site2=o2[1]
    if(site1==site2):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["number", site1,np.conj(p[i,imin])*p[j,imin]])
    elif(site2<site1):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["hopping", site1, site2, np.conj(p[i,imin])*p[j,imin]])
  if(o1[0]=="annihilation" and o2[0]=="an"):
    site1=o1[1]
    site2=o2[1]
    site3=o2[2]
    if(site1==site2):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["interaction", site1, site3,np.conj(p[i,imin])*p[j,imin]])
    elif(site2<site1):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["hn", site1, site2, site3,np.conj(p[i,imin])*p[j,imin]])
  if(o1[0]=="an" and o2[0]=="an"):
    site1=o1[1]
    site2=o1[2]
    site3=o2[1]
    site4=o2[2]
    if(site2==site2):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["3n", site2, site1, site4,np.conj(p[i,imin])*p[j,imin]])
    elif(site2<site1):
      #if(np.conj(p[i,imin])*p[j,imin]>1E-10):
        Hamiltonian.append(["hn", site2, site1, site3,site4,np.conj(p[i,imin])*p[j,imin]])
  if(o1[0]=="number" and o2[0]=="number"):
    site1=o1[1]
    site2=o2[1]
    if(site1==site2):
      Hamiltonian.append(["number", site1,np.conj(p[i,imin])*p[j,imin]])
    else:
      Hamiltonian.append(["interaction", site1, site2, np.conj(p[i,imin])*p[j,imin]])
  if (o1[0]=="constant" and o2[0]=="constant"):
    Hamiltonian.append(["constant", site1,np.conj(p[i,imin])*p[j,imin]])
print("Hamiltonian:")
print(Hamiltonian) 
"""
plt.show()
