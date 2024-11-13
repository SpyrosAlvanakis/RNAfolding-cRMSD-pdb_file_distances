import pandas as pd
from Bio.PDB import *
import numpy as np
import warnings
warnings.simplefilter('ignore')
import matplotlib.pyplot as plt
import random
parser = PDBParser(PERMISSIVE=1)


def get_ca_atoms(pdb_file, start_resnum=102, num_residues=50):
    """
    Parses a PDB file and returns the coordinates of alpha carbons (CA) for a
    specified range of residues.
    
    Args:
        pdb_file (str): The name of the PDB file to parse.
        start_resnum (int): The residue number to start extracting alpha carbons from.
        num_residues (int): The number of residues to extract alpha carbons from.
        
    Returns:
        A tuple containing two lists: the coordinates of alpha carbons (CA) and the
        residue names for each alpha carbon.
    """
    structure = parser.get_structure('6lu7', pdb_file)
    ca_atoms = []
    residue_count = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resnum = residue.get_id()[1]
                if resnum >= start_resnum and resnum < start_resnum + num_residues:
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            ca_atoms.append(atom.get_coord())
                            residue_count += 1
                            break
                if residue_count == num_residues:
                    break
            if residue_count == num_residues:
                break
    
    return np.array(ca_atoms)

def CM_mtrx(atoms):
    """
    Computes the Cayley-Menger matrix for a given set of 3D coordinates of atoms.

    The Cayley-Menger matrix is a matrix of squared distances between all pairs of atoms, 
    with the diagonal being all zeros. This matrix is used in the computation of the 
    rank of the Cayley-Menger matrix, which is a measure of the number of degrees of 
    freedom of the system.

    Args:
        atoms (numpy.ndarray): The 3D coordinates of the atoms, shape (n_atoms, 3).

    Returns:
        The Cayley-Menger matrix, shape (n_atoms+1, n_atoms+1).
    """
    n = atoms.shape[0]
    cm_matrix = np.zeros(shape=(n+1, n+1))
    cm_matrix[0,:]=1
    cm_matrix[:,0]=1
    cm_matrix[0,0]=0
    for i in range(n):
        for j in range(n):
            if i <= j:
                val_1 = atoms[i]
                val_2 = atoms[j] 
                dist = np.linalg.norm(val_1-val_2)
                cm_matrix[i+1][j+1] = (dist**2)/2
                cm_matrix[j+1][i+1] = (dist**2)/2
            elif i == j:
                cm_matrix[i,j] = 0

    return cm_matrix


def perturb(mtrx,perc):
    """
    Perturb the Cayley-Menger matrix by a given percentage.

    The Cayley-Menger matrix is a matrix of squared distances between all pairs of atoms, 
    with the diagonal being all zeros. This function perturbs the matrix by a given percentage, 
    by multiplying each off-diagonal element by a random number between 1-perc/100 and 1+perc/100.

    Args:
        mtrx (numpy.ndarray): The Cayley-Menger matrix, shape (n_atoms+1, n_atoms+1).
        perc (float): The percentage of perturbation.

    Returns:
        The perturbed Cayley-Menger matrix, shape (n_atoms+1, n_atoms+1).
    """
    pert_mtrx=np.zeros(shape=mtrx.shape)
    pert_mtrx[0,:]=1
    pert_mtrx[:,0]=1
    pert_mtrx[0,0]=0
    for i in range(1,mtrx.shape[0]):
        for j in range(1,mtrx.shape[1]):
            if i==0 or j==0: 
                continue
            elif i<j:  # since the diagonal is zero
                num = np.random.uniform(1-perc/100, 1+perc/100)
                pert_val = abs(num*mtrx[i,j])
                pert_mtrx[i,j] = pert_val
                pert_mtrx[j,i] = pert_val
    return pert_mtrx            
         
         
def Gram_mtrx(mtrx):
    """
    Computes the Gram matrix for a given Cayley-Menger matrix.

    The Gram matrix is a matrix of dot products of the columns of the Cayley-Menger
    matrix, and is used in the computation of the rank of the Cayley-Menger matrix.

    Args:
        mtrx (numpy.ndarray): The Cayley-Menger matrix, shape (n_atoms+1, n_atoms+1).

    Returns:
        The Gram matrix, shape (n_atoms, n_atoms).
    """
    gram_mtrx=np.zeros(shape=(mtrx.shape[0]-1,mtrx.shape[1]-1))
    for i in range(gram_mtrx.shape[0]):
        for j in range(gram_mtrx.shape[1]):
            gram_mtrx[i,j]=mtrx[i+1,1]+mtrx[j+1,1]-mtrx[i+1,j+1]
    return gram_mtrx        

def S_arr(mtrx):
    """
    Computes the diagonal matrix of singular values of a given matrix.

    The matrix is decomposed using the singular value decomposition (SVD) into
    the product of three matrices: U, S, and V. The diagonal matrix of singular
    values is then constructed by taking the first three singular values and
    arranging them in a diagonal matrix.

    Args:
        mtrx (numpy.ndarray): The matrix to be decomposed, shape (n_atoms, n_atoms).

    Returns:
        A tuple containing the diagonal matrix of singular values and the matrix V
        from the SVD decomposition.
    """
    _ , s, vT = np.linalg.svd(mtrx, full_matrices=True)
    sorted_values = np.sort(s.flatten(),kind='mergesort')
    rev_sort=sorted_values[::-1]
    arr = np.zeros((3,3))
    for i in range(mtrx.shape[0]):
        if i < 3:
            arr[i, i] = rev_sort[i]
        else:
            break
    return arr, vT 

def centered(array):
    """
    Centers the input array by subtracting the mean along each column.

    Args:
        array (numpy.ndarray): Input array to be centered, shape (n_samples, n_features).

    Returns:
        numpy.ndarray: Centered array with the mean of each column subtracted, same shape as input.
    """
    mean = np.mean(array, axis=0)
    arr_centered = array - mean
    return arr_centered


def crmsd(mtrx_1,mtrx_2):
    """
    Computes the centered root mean square deviation (cRMSD) between two conformation matrices.

    Args:
        mtrx_1 (numpy.ndarray): The first conformation matrix, shape (n_atoms, 3).
        mtrx_2 (numpy.ndarray): The second conformation matrix, shape (n_atoms, 3).

    Returns:
        float: The cRMSD between the two matrices.
    """
    temp_mtrx = mtrx_1.T@mtrx_2
    u, s, vT = np.linalg.svd(temp_mtrx, full_matrices=True)
    Q_mtrx = u@vT
    if np.linalg.det(Q_mtrx)<0:
        u[:,-1]*= -1
        Q_mtrx = u@vT
    A_mtrx = mtrx_1@Q_mtrx.T-mtrx_2
    value = np.sqrt(np.sum(pow((np.linalg.norm(A_mtrx)),2)/mtrx_1.shape[0]))
    return value