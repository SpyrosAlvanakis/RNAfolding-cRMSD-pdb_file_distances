import matplotlib.pyplot as plt
from utils import utils
import os 
import numpy as np

# Get the directory where the script is located
current_path = os.path.dirname(os.path.abspath(__file__))  # Get the current working directory
results_path = os.path.join(current_path, "results") # Add the results folder

# Extract CA atoms
ca_atoms = utils.get_ca_atoms(f'{current_path}/data/6lu7.pdb',start_resnum=102, num_residues=50)

# Compute CM matrix
cm_matrix = utils.CM_mtrx(ca_atoms)

# Plot
fig, ax = plt.subplots(figsize=(15, 10)) 
im = ax.imshow(cm_matrix, cmap='plasma')
plt.colorbar(im)
plt.title('Cayley-Menger matrix')

# Save the plot to the results directory
output_file_cm = os.path.join(results_path, "cayley_menger_matrix.png")
plt.savefig(output_file_cm)

# Compute rank
rank=np.linalg.matrix_rank(cm_matrix, tol = 0.01)
print('############## Answer to question 1 ##############')
print(f'The rank of the Cayley-Menger matrix is {rank}')
print(f"Cayley-Menger matrix plot saved at {output_file_cm}")
print('##################################################')

# Compute perturbed matrix
pert_mtrx = utils.perturb(cm_matrix,5)

# Calculate rank
rank=np.linalg.matrix_rank(pert_mtrx, tol = 0.01)

# Compute Gram matrix
gram_mtrx=utils.Gram_mtrx(pert_mtrx)

# Compute SVD
Sigma_mtrx, vT = utils.S_arr(gram_mtrx)

# Compute coordinates
coordinates = (np.sqrt(Sigma_mtrx) @ vT[:3,:]).T

# Center
cent_ca_atoms = utils.centered(ca_atoms)
cen_coordinates = utils.centered(coordinates)

# Compute cRMSD
cRMSD = utils.crmsd(cent_ca_atoms,cen_coordinates)

# Plot of Gram matrix
fig, ax = plt.subplots(figsize=(15,10)) 
im = ax.imshow(gram_mtrx, cmap='inferno')
plt.colorbar(im)
plt.title('Gram matrix')

# Save the plot to the results directory
output_file_grm = os.path.join(results_path, "gram_matrix.png")
plt.savefig(output_file_grm)

print('############## Answer to question 2 ##############')
print(f'The cRMSD between the two matrices is {cRMSD}')
print(f'The rank of the perturbated matrix is {rank}')
print(f"Gram matrix plot saved at {output_file_grm}")
print('##################################################')