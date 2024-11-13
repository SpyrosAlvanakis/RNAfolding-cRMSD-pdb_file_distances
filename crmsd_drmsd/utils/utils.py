import numpy as np
import pandas as pd
from scipy.linalg import svd

def read_conformations(path):
    """
    Reads a conformation data file and returns a DataFrame with molecular data.

    Args:
        path (str): The file path to the conformation data, expected to be a tab-separated text file.

    Returns:
        pandas.DataFrame: A DataFrame containing columns 'X', 'Y', 'Z' for coordinates and 'mol' for molecule identifiers
        ranging from Mol_1 to Mol_10, with each molecule containing 369 records.
    """
    col_names=['X', 'Y', 'Z']
    data= pd.read_csv(path, sep='\t', skiprows=2, names=col_names)
    molecule_names= ['mol_{}'.format(i) for i in range(1, 11)]
    for i in range(1, 11):
        start_index= (i-1)*369
        end_index= i*369
        data.loc[start_index:end_index,'mol']= f'Mol_{i}'
    return data 


def crmsd(mtrx_1,mtrx_2):
    """
    Calculates the centered root mean square deviation (cRMSD) between two conformation matrices.

    Args:
        mtrx_1 (numpy.ndarray): The first conformation matrix, shape (n_atoms, 3).
        mtrx_2 (numpy.ndarray): The second conformation matrix, shape (n_atoms, 3).

    Returns:
        float: The cRMSD between the two matrices.
    """
    mtrx_1 -= np.mean(mtrx_1, axis = 0)
    mtrx_2 -= np.mean(mtrx_2, axis = 0)
    temp_mtrx = mtrx_1.T@mtrx_2
    u, s, vT = np.linalg.svd(temp_mtrx, full_matrices=True)
    Q_mtrx = u@vT
    if np.linalg.det(Q_mtrx)<0:
        u[:,-1]*= -1
        Q_mtrx = u@vT
    A_mtrx = mtrx_1@Q_mtrx.T-mtrx_2
    value = np.sqrt(np.sum(pow((np.linalg.norm(A_mtrx)),2)/mtrx_1.shape[0]))
    return value

def df_format(mtrx,names):
    """
    Formats a matrix into a DataFrame, and prints out the column name of the column with the minimum mean value.

    Args:
        mtrx (numpy.ndarray): The matrix to be formatted, shape (n_samples, n_features).
        names (list): A list of column names.

    Returns:
        pandas.DataFrame: The formatted DataFrame.

    """
    col_means = np.mean(mtrx, axis=0)
    min_col_index = np.argmin(col_means)
    data = pd.DataFrame(mtrx, columns=names)  
    min_col_name = data.columns[min_col_index]
    return data