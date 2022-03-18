def add_operator(otype, sites, basis, basis_dict, operators, operator_descriptions):
 if(otype=="number"):
   operator_descriptions.append(["number", sites[0], 1])
   operators.append(number_operator(basis, sites[0]))  # Append the number operator on site 1
 if(otype=="interaction"): 
   operator_descriptions.append(["interaction", sites[0], sites[1], 1])
   operators.append(interaction(basis, sites[0], sites[1]))
 if(otype=="rhopping"):
   operator_descriptions.append(["rhopping", sites[0], sites[1], 1])
   operators.append(hopping(basis, basis_dict, sites[1], sites[0]) + hopping(basis, basis_dict, sites[0], sites[1]))
 if(otype=="ihopping"):
   operator_descriptions.append(["ihopping", sites[0], sites[1], 1])
   operators.append(-1j*hopping(basis, basis_dict, sites[1], sites[0]) + 1j*hopping(basis, basis_dict, sites[0], sites[1]))
 if(otype=="annihilation"):
   operator_descriptions.append(["annihilation", sites[0], 1])
   operators.append(-1j*hopping(basis, basis_dict, sites[1], sites[0]) + 1j*hopping(basis, basis_dict, sites[0], sites[1]))



 return operators, operator_descriptions


def add_operator2(otype, sites, basis, basis_dict,basis2, basis_dict2, operators, operator_descriptions):
 if(otype=="number"):
   operator_descriptions.append(["number", sites[0], 1])
   operators.append(number_operator(basis, sites[0]))  # Append the number operator on site 1
 if(otype=="interaction"): 
   operator_descriptions.append(["interaction", sites[0], sites[1], 1])
   operators.append(interaction(basis, sites[0], sites[1]))
 if(otype=="rhopping"):
   operator_descriptions.append(["rhopping", sites[0], sites[1], 1])
   operators.append(hopping(basis, basis_dict, sites[1], sites[0]) + hopping(basis, basis_dict, sites[0], sites[1]))
 if(otype=="ihopping"):
   operator_descriptions.append(["ihopping", sites[0], sites[1], 1])
   operators.append(-1j*hopping(basis, basis_dict, sites[1], sites[0]) + 1j*hopping(basis, basis_dict, sites[0], sites[1]))
 if(otype=="annihilation"):
   operator_descriptions.append(["annihilation", sites[0], 1])
   operators.append(annihilation_operator(basis, basis2, basis_dict2, sites[0]))
 if(otype=="an"):
   operator_descriptions.append(["an", sites[0],sites[1], 1])
   operators.append(annihilation_operator(basis, basis2, basis_dict2, sites[0]) @ number_operator(basis, sites[1]))
 return operators, operator_descriptions


def operators_from_descriptions(operator_descriptions, bases, dN):
 operators=[]
 if(dN==0):
   basis, basis_dict = bases
   for op in operator_descriptions:
    if(op[0]=="number"):
     operators.append(number_operator(basis, op[1]))  # Append the number operator on site 1
    if(op[0]=="interaction"):
     operators.append(number_operator(basis, op[1], op[2]))  # Append the number operator on site 1
    if(op[0]=="rhopping"):
     operators.append(hopping(basis, basis_dict, op[2], op[1]) + hopping(basis, basis_dict, op[1], op[2]))
    if(op[0]=="ihopping"):
     operators.append(-1j*hopping(basis, basis_dict, op[2], op[1]) + 1j*hopping(basis, basis_dict, op[1], op[2])) 
 if(dN==1):
   basis, basis_dict,basis2, basis_dict2= bases
   for op in operator_descriptions:
    if(op[0]=="annihilation"):
      operators.append(annihilation_operator(basis, basis2, basis_dict2, op[1]))
    if(op[0]=="an"):
      operators.append(annihilation_operator(basis, basis2, basis_dict2, op[1]) @ number_operator(basis, op[2]) )
 if(dN==2):
   basis, basis_dict,basis2, basis_dict2,basis3, basis_dict3 = bases
   for op in operator_descriptions:
    if(op[0]=="aa"):
      operators.append(annihilation_operator(basis2, basis3, basis_dict3, op[1]) @   annihilation_operator(basis, basis2, basis_dict2, op[2]) )
 return operators

def update_coeffs(operator_descriptions, coeffs):
  from copy import deepcopy
  #print(operator_descriptions)
  operator_descriptions2=deepcopy(operator_descriptions)
  for i, o in enumerate(operator_descriptions2):
    o[-1]=coeffs[i]
  return operator_descriptions2

   
def operator_generator(basis, basis_dict, z, z_dict, offsets):
    """
    Generates a list of operators
    :param basis: List of basis elements
    :param z: List of lattice positions
    :return:
    """
    for site1, z1 in enumerate(z): # Loop through all lattice sites
        for offset in offsets: # Site 1 can interact with the site above it and to the right of it. 
               
            operators = [number_operator(basis, site1)]  # Append the number operator on site 1

            z2 = z1 + offset # Find the position of a neighbouring site
            if z2 in z: # If this site is a part of the lattice
                site2 = z_dict[z2] # Find the index of this lattice site
                # Append the interaction and hopping between these sites
                print(site1, site2, offset)
                operators = operators + [ #number_operator(basis, site1),
                                        interaction(basis, site1, site2),
                                         hopping(basis, basis_dict, site2, site1) + hopping(basis, basis_dict, site1, site2),
                                         -1j*hopping(basis, basis_dict, site2, site1) +1j*hopping(basis, basis_dict, site1, site2)]
            yield operators, site1, site2


def operator_generator2(basis, basis_dict, z, z_dict):
    """
    Generates a list of operators
    :param basis: List of basis elements
    :param z: List of lattice positions
    :return:
    """
    for site1, z1 in enumerate(z): # Loop through all lattice sites
        for offset in [1, 1j]: # Site 1 can interact with the site above it and to the right of it. 
            operators=[]
            operator_descriptions=[]
            operators, operator_descriptions = add_operator("number", [site1],basis, basis_dict, operators, operator_descriptions)# Append the number operator on site 1
            z2 = z1 + offset # Find the position of a neighbouring site
            if z2 in z: # If this site is a part of the lattice
                site2 = z_dict[z2] # Find the index of this lattice site
                # Append the interaction and hopping between these sites
                print(site1, site2, offset)
                operators, operator_descriptions = add_operator("interaction", [site1, site2], basis, basis_dict, operators, operator_descriptions)# Append the interaction between sites 1 and 2
                operators, operator_descriptions = add_operator("rhopping", [site1, site2], basis, basis_dict, operators, operator_descriptions)# Append the real hopping between sites 1 and 2
                operators, operator_descriptions = add_operator("ihopping", [site1, site2], basis, basis_dict, operators, operator_descriptions)# Append the imaginary hopping between sites 1 and 2
            yield operators,operator_descriptions, site1, site2

def all_operators(basis, basis_dict, z, z_dict, offsets):
  generator = operator_generator(basis, basis_dict, z, z_dict, offsets)
  operators=[]
  for operators0, site1, site2 in(generator):
    operators.extend(operators0)
  return operators

def generate_M_matrix(operator_descriptions, psi, basis, basis_dict, z, z_dict):
  import numpy as np
  sizeM=len(operator_descriptions)
  psi_T = np.conj(np.transpose(psi)) # Compute the hermitian conjugate of psi
  M=np.zeros((sizeM, sizeM), dtype=complex)
  for i,o1 in enumerate(operator_descriptions):
    for j,o2 in enumerate(operator_descriptions):
      operator=None
      if(o1[0]=="annihilation" and o2[0]=="annihilation"):
        site1=o1[1]
        site2=o2[1]
        if(site1==site2):
          operator = number_operator(basis, site1)  # Append the number operator on site 1
        else:
          operator=hopping(basis, basis_dict, site2, site1)
      if(o1[0]=="annihilation" and o2[0]=="an"):
        site1=o1[1]
        site2=o2[1]
        site3=o2[2]
        if(site1==site2)and(site2==site3):
          operator = number_operator(basis, site1)  
        if(site1==site2)and(site2!=site3):
          operator = interaction(basis, site1, site3)  
        if(site1!=site2):
          operator1 = hopping(basis, basis_dict, site1, site2)  
          operator2 = number_operator(basis, site3)
          operator = operator1 @ operator2 
      if(o1[0]=="an" and o2[0]=="annihilation"):
        site1=o1[1]
        site2=o1[2]
        site3=o2[1]
        if(site1==site3)and(site2==site3):
          operator = number_operator(basis, site1)  
        if(site1==site3)and(site2!=site3):
          operator = interaction(basis, site1, site2)  
        if(site1!=site3):
          operator1 = number_operator(basis, site2)
          operator2 = hopping(basis, basis_dict, site1, site3)  
          operator = operator1 @ operator2 
      if(o1[0]=="an" and o2[0]=="an"):
        site1=o1[1]
        site2=o1[2]
        site3=o2[1]
        site4=o2[2]        
        operator1 = number_operator(basis, site2)
        operator3 = number_operator(basis, site4)
        if(site1==site3):
          operator2 = number_operator(basis, site1)
        else:  
          operator2 = hopping(basis, basis_dict, site1, site3)  
        operator = operator1 @ operator2 @ operator3 
      if(o1[0]=="number" and o2[0]=="number"):
        site1=o1[1]
        site2=o2[1]
        if(site1==site2):
          operator = number_operator(basis, site1)  # n**2=n
        else:
          operator=interaction(basis, site1, site2)
      if operator is not None:
        print("Operator:")
        print(operator)
        M[i,j]=psi_T @ operator @ psi
      elif (o1[0]=="constant" and o2[0]=="constant"):
        M[i,j]=1
  return M    
def find_combination_annihilation(psi, operator_descriptions, basis, basis_dict, z, z_dict):
    """
    Uses the method from reference https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031029 to find
    the linear combinations of operators which has psi as an eigenvector
    :param psi: Wave function
    :param operators: List of operators
    :return:
    """
    import numpy as np
    psi_T = np.conj(np.transpose(psi)) # Compute the hermitian conjugate of psi
    M=generate_M_matrix(operator_descriptions, psi, basis, basis_dict, z, z_dict)
    print("M:", M)
    D, P = np.linalg.eig(M) # Diagonalize C to find the zero eigenvalues
    return D, P # Return eigenvalues and eigenvectors.
        

def find_combination(psi, operators):
    """
    Uses the method from reference https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031029 to find
    the linear combinations of operators which has psi as an eigenvector
    :param psi: Wave function
    :param operators: List of operators
    :return:
    """
    import numpy as np
    psi_T = np.conj(np.transpose(psi)) # Compute the hermitian conjugate of psi
    C = np.zeros((len(operators), len(operators)), dtype=np.complex128) # Initialize correlation matrix
    for i, o1 in enumerate(operators): # Loop through all operators
        print(f"i = {i}")
        for j, o2 in enumerate(operators): # --||--
            C[i, j] = (psi_T @ o1 @ o2 @ psi) - (psi_T @ o1 @ psi) * (psi_T @ o2 @ psi) # Compute the current matrix element

    D, P = np.linalg.eig(C) # Diagonalize C to find the zero eigenvalues
    return D, P # Return eigenvalues and eigenvectors.

def find_combination_2nd_method(psi, operators):
    """
    Uses the method from reference https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031029 to find
    the linear combinations of operators which has psi as an eigenvector
    :param psi: Wave function
    :param operators: List of operators
    :return:
    """
    import numpy as np
    psi_T = np.conj(np.transpose(psi)) # Compute the hermitian conjugate of psi
    C = np.zeros((len(operators), len(operators)), dtype=np.complex128) # Initialize correlation matrix
    for i, o1 in enumerate(operators): # Loop through all operators
        print(f"i = {i}")
        for j, o2 in enumerate(operators): # --||--
            C[i, j] = (psi_T @ np.conj(np.transpose(o1)) @ o2 @ psi) # Compute the current matrix element

    D, P = np.linalg.eig(C) # Diagonalize C to find the zero eigenvalues
    return D, P # Return eigenvalues and eigenvectors.

def number_operator(basis, site):
    """
    Computes the number operator on a site in a basis
    These matrices are very large, so we use sparse matrices.
    Note that lil_matrix efficiently builds a matrix while csc_matrix efficiently does calculations.
    :param basis: List of basis elements
    :param site: The site we compute the number operator for.
    :return:
    """
    from scipy.sparse import lil_matrix
    n = lil_matrix((len(basis), len(basis))) # Initialization of operator
    for i, b in enumerate(basis): # Loop through all basis elements
        if site in b: # If the site is in the basis
            n[i, i] = 1 # Add 1
    return n.tocsc() # Convert to csc format so arithmetic is efficient

def interaction(basis, site1, site2):
    """
    Computes the interaction operator between site 1 and site 2
    :param basis: List of basis elements
    :param site1: Site 1
    :param site2: Site 2
    :return:
    """
    from scipy.sparse import lil_matrix
    inter = lil_matrix((len(basis), len(basis))) # Initialize the interaction operator
    for i, b in enumerate(basis): # Loop though all basis elemnts
        if site1 in b and site2 in b: # If both site 1 and site 2 are in the basis element they interact
            inter[i, i] = 1 # Add 1
    return inter.tocsc()

def hopping(basis, basis_dict, site1, site2):
    """
    Computes the hopping operator between site 1 and site 2.
    Note that this operator is NOT hermitian. Add the conjugate to make it hermittian
    ie. hopping(basis, basis_dict, site1, site2) + hopping(basis, basis_dict, site2, site1) IS hermittian.
    :param basis: List of basis elements
    :param basis_dict: Dictionary which returns the index in the reduced basis corresponding to a basis element
    :param site1: Site 1
    :param site2: Site 2
    :return:
    """
    from scipy.sparse import lil_matrix
    hop = lil_matrix((len(basis), len(basis))) # Initialization of hopping operator
    for i, b in enumerate(basis): # Loop through all basis elements
        if site1 in b and site2 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site1)] + b[b.index(site1) + 1:] + (site2, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = 1 # Add 1 to this entry in the matrix.
    return hop.tocsc()


def annihilation_operator(basis, basis2, basis_dict2, site1):
    """
    Computes the hopping operator between site 1 and site 2.
    Note that this operator is NOT hermitian. Add the conjugate to make it hermittian
    ie. hopping(basis, basis_dict, site1, site2) + hopping(basis, basis_dict, site2, site1) IS hermittian.
    :param basis: List of basis elements
    :param basis_dict: Dictionary which returns the index in the reduced basis corresponding to a basis element
    :param site1: Site 1
    :param site2: Site 2
    :return:
    """
    from scipy.sparse import lil_matrix
    ann = lil_matrix((len(basis2), len(basis))) # Initialization of hopping operator
    for i, b in enumerate(basis): # Loop through all basis elements
        if site1 in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site1)] + b[b.index(site1) + 1:])) # Remove site 1 from the basis element, add site 2.
            j = basis_dict2[b_new] # Find the index of this new basis element
            #print("i:", i,"b:",  b, "b_new:", b_new, "j:", j)
            ann[j, i] = 1 # Add 1 to this entry in the matrix.
    return ann.tocsc()

def nhn(basis, basis_dict, site1, site2, site3,site4):
    """
    n_{site1}c^{dagger}_{site2}c_{site3}n_{site4}
    """
    from scipy.sparse import lil_matrix
    hop = lil_matrix((len(basis), len(basis))) # Initialization of hopping operator
    for i, b in enumerate(basis): # Loop through all basis elements
        mult=1
        if site1 not in b:
          mult=0
        if site4 not in b:
          mult=0
        if site2 in b and site3 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site2)] + b[b.index(site2) + 1:] + (site3, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = 1 # Add 1 to this entry in the matrix.
    return hop.tocsc()

def hopping_hermitian_real(basis, basis_dict, site1, site2):
    """
    Computes the hopping operator between site 1 and site 2 - hermitian version
    :param basis: List of basis elements
    :param basis_dict: Dictionary which returns the index in the reduced basis corresponding to a basis element
    :param site1: Site 1
    :param site2: Site 2
    :return:
    """
    from scipy.sparse import lil_matrix
    hop = lil_matrix((len(basis), len(basis))) # Initialization of hopping operator
    for i, b in enumerate(basis): # Loop through all basis elements
        if site1 in b and site2 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site1)] + b[b.index(site1) + 1:] + (site2, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = 1 # Add 1 to this entry in the matrix.
        if site2 in b and site1 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site2)] + b[b.index(site2) + 1:] + (site1, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = 1 # Add 1 to this entry in the matrix.
    return hop.tocsc()

def hopping_hermitian_imag(basis, basis_dict, site1, site2):
    """
    Computes the hopping operator between site 1 and site 2 - hermitian version with imaginary coefficient
    :param basis: List of basis elements
    :param basis_dict: Dictionary which returns the index in the reduced basis corresponding to a basis element
    :param site1: Site 1
    :param site2: Site 2
    :return:
    """
    from scipy.sparse import lil_matrix
    hop = lil_matrix((len(basis), len(basis))) # Initialization of hopping operator
    for i, b in enumerate(basis): # Loop through all basis elements
        if site1 in b and site2 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site1)] + b[b.index(site1) + 1:] + (site2, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = 1j # Add 1 to this entry in the matrix.
        if site2 in b and site1 not in b: # If site 1 is in the basis element but site 2 is not, then the particle on site 1 can hop to site 2
            b_new = tuple(sorted(b[:b.index(site2)] + b[b.index(site2) + 1:] + (site1, ))) # Remove site 1 from the basis element, add site 2.
            j = basis_dict[b_new] # Find the index of this new basis element
            hop[j, i] = -1j # Add 1 to this entry in the matrix.
    return hop.tocsc()


def combine_operators(operator1, operator2, value):
 if(operator1[0]=="an"):
  if(operator1[1]==operator1[2]):
   operator1=["annihilation", operator1[1], 1]
 if(operator1[0]=="an"):
  if(operator1[1]==operator1[2]):
   operator1=["annihilation", operator1[1], 1]
 if(operator1[0]=="annihilation")and(operator2[0]=="annihilation"):
  if(operator1[1]==operator2[1]):
    return ["number", operator1[1], value] 
  else:
    return ["hopping", operator1[1], operator2[1], value]
 if(operator1[0]=="an")and(operator2[0]=="annihilation"):
   site1=operator1[2] 
   site2=operator1[1] 
   site3=operator2[1] 
   print(operator1[0], operator2[0], site1, site2, site3)
   if(site1==site2)and(site2==site3):
     return ["number", site1, value] 
   if(site1!=site2)and(site2==site3):
     return ["interaction", site1, site2, value] 
   if(site1!=site2)and(site1!=site3):
     return ["nh", site1, site2, site2, value] 
 if(operator1[0]=="annihilation")and(operator2[0]=="an"):
   site1=operator1[1] 
   site2=operator2[1] 
   site3=operator2[2] 
   if(site1==site2)and(site2==site3):
     return ["number", site1, value] 
   if(site1==site2)and(site2!=site3):
     return ["interaction", site1, site2, value] 
   if(site1!=site2)and(site1!=site3):
     return ["hn", site1, site2, site2, value] 
 if(operator1[0]=="an")and(operator2[0]=="an"):
   site1=operator1[2] 
   site2=operator1[1] 
   site3=operator2[1] 
   site4=operator2[2]
   if(site2==site3):
     sites=[site1,site2,site4] 
     sites = list(dict.fromkeys(sites)) #remove repetitions
     print(site1,site2, site3, site4, sites)
     if(len(sites)==1):
       return ["number", sites[0], value] 
     if(len(sites)==2):
       return ["interaction", sites[0], sites[1], value] 
     if(len(sites)==3):
       return ["3n", sites[0], sites[1], sites[2], value] 
   else:
     if(site1!=site3)and(site2!=site4):
      return ["nhn", site1, site2, site3, site4, value] 
 

if __name__ == '__main__':
    main()
