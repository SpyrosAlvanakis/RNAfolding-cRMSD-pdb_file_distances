from utils import utils
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Get the current working directory
current_path = os.path.dirname(os.path.abspath(__file__))  # Get the current working directory
results_path = os.path.join(current_path, "results") # Add the results folder

# Read conformations
conformations = utils.read_conformations(f'{current_path}/data/10_conformations.txt')

# Calculate c-RMSD
start_time = time.time()
cRMSD_mtrx=np.zeros(shape=(10,10))
row=0
for base_mol in conformations.mol.unique():
    subset_base = conformations.loc[conformations.mol == base_mol, ['X', 'Y', 'Z']].values
    col=0
    for comp_mol in conformations.mol.unique():
        if row<col:
            subset_comp = conformations.loc[conformations.mol == comp_mol, ['X', 'Y', 'Z']].values
            mtrx_value = utils.crmsd(subset_base,subset_comp)
            cRMSD_mtrx[row, col] = mtrx_value
            cRMSD_mtrx[col, row] = mtrx_value
        elif row == col: cRMSD_mtrx[row,col] = 0
        col+=1
    row+=1              
end_time = time.time()

cRMSD_time= end_time - start_time
print(f'c-RMSD calculation time: {cRMSD_time}')

cRMSD_df = utils.df_format(cRMSD_mtrx,conformations.mol.unique())

# Plot the heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(cRMSD_df, annot=True, cmap='Blues', fmt='.0f')
plt.title('cRMSD Matrix')

# Save the plot to the results folder
output_file = os.path.join(results_path, "cRMSD_matrix.png")
plt.savefig(output_file)
print(f"Plot saved at {output_file}")

