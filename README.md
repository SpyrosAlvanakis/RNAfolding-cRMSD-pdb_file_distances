# Algorithms in Structural Bioinformatics

This repository contains the solutions to **Assignment 1** for the course *Algorithms in Structural Bioinformatics*. Each part of the assignment is organized into a separate folder, containing the code, plots, and a PDF report detailing the solution.

## Repository Structure

The repository is structured as follows:

```
.
├── rna_folding/
│   ├── results/         # Visualizations of RNA secondary structures
│   └── main_1.py        # Script for part rna folding
│
├── crmsd_drmsd/
│   ├── data/            # Input data 
│   ├── utils/           # Functions of c-RMSD and d-RMSD calculations
│   ├── results/         # Relevant visualizations
│   └── main_2.py        # Script for part of c-RMSD and d-RMSD
│
├── distances/
│   ├── data/            # Input data of a .pdb file
│   ├── utils/           # Computation functions involving Cayley-Menger matrix and SVD
│   ├── results/         # Relevant visualizations
│   └── main_3.py        # Script for part of distances
│
├── requirements.txt     # Python 3.11 environment dependencies
├── Reports .pdf
└── README.md            # This file
```

## Description of Each Part

### Part 1: RNA Folding
- **Objective**: Find all optimal secondary structures of a given RNA sequence using a crude energy minimization algorithm.
- **Tasks**:
  - Implement the algorithm to fill the energy table.
  - Draw all optimal folds, showing bonds and backtrack paths.
- **Tools**: Implemented in Python (or a similar system).

### Part 2: c-RMSD and d-RMSD
- **Objective**: Calculate c-RMSD and d-RMSD for 10 molecular conformations and analyze centroid identification and performance.
- **Tasks**:
  1. Compute c-RMSD distances for all pairs and find the L1-centroid.
  2. Repeat using d-RMSD with two approaches: all distances and a random subset.
  3. Compare results and performance of the three methods.
- **Tools**: Utilizes linear algebra functions like SVD.

### Part 3: Distances
- **Objective**: Analyze the Cayley-Menger matrix for 50 α-carbon atoms from the SARS-CoV-2 protease structure.
- **Tasks**:
  1. Compute the rank of the matrix and justify results.
  2. Perturb the matrix and analyze changes in rank.
  3. Recover 3D coordinates and compute c-RMSD against the original structure.
- **Tools**: Uses Gram matrix computation and SVD decomposition.

## Installation and Setup

1. Clone this repository:
   ```bash
   git clone <repository_url>
   cd <repository_name>
   ```

2. Create a virtual environment with Python 3.11:
   ```bash
   python3.11 -m venv .venv
   source .venv/bin/activate 
   ```

3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```
