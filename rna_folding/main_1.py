import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os 

# Defining the sequence
seq='AAUACUCCGUUGCAGCAU'

# Defining the matrix
size= len(seq)
matrix= np.ones(shape=(size,size))*100

# Defining the score
def score(let1,let2):
    """
    This function takes two nucleotides and returns the score of the base pairing according to the following rules:
    A-U or U-A: -4
    G-U or U-G: 0
    All other pairs: 4
    """
    if (let1=='A' and let2=='U') or (let1=='U' and let2=='A') or (let1=='G' and let2=='C') or (let1== 'C' and let2==' G') :
        scor= -4
    elif (let1=='G' and let2=='U') or (let2== 'G' and let1== 'U'):
        scor=0
    else :
        scor=4
    return scor    

# Filling the matrix
for row in range(size-1,-1,-1):
    for col in range(row+5, size):
        let1= seq[col] 
        let2= seq[row] 
        damage= score(let1,let2) 
        temp1= matrix[row,col-1]
        temp2= matrix[row+1,col]
        temp3= matrix[row+1,col-1]+ damage
        tempk=[]
        for k in range(row+2,col):
            tempk.append(matrix[k,col]+matrix[row,k-1])    
        matrix[row,col]=min(temp1,temp2,temp3,min(tempk))
        
# Get the directory where the script is located
script_path = os.path.dirname(os.path.abspath(__file__)) 
results_path = os.path.join(script_path, "results") 

# Plotting
plt.figure(figsize=(15,10))
sns.heatmap(matrix, annot=True, cmap='Blues', fmt='.0f')
plt.title('Base Pairing Matrix')
output_file = os.path.join(results_path, "base_pairing_matrix.png")
plt.savefig(output_file)
print(f"Plot saved at {output_file}")