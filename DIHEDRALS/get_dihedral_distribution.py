########################### Code Overview #####################################
# This script computes the backbone dihedral angle distribution for a polymer
# system from LAMMPS trajectory data.
#
# The code is written for polypropylene (PP) and polystyrene (PS) chains with
# fixed monomer connectivity, and extracts only backbone dihedral angles from
# the atomic coordinates.
#
# Workflow:
#   1. Read the number of atoms and dihedrals, and the dihedral connectivity,
#      from `structure_300.lammps-data`.
#   2. Read the list of simulation seeds from `input.dat`.
#   3. Define the backbone atom IDs based on the polymer type (PP or PS).
#   4. Loop over all seeds and read atomic coordinates from
#      `concatenated_atom_info.dat`.
#   5. Extract the last 10 trajectory blocks from each seed.
#   6. Compute backbone dihedral angles for all valid dihedrals.
#   7. Combine angles from all seeds and frames.
#   8. Estimate the dihedral angle distribution using Gaussian kernel density
#      estimation (KDE).
#   9. Save the smooth dihedral distribution to `dihedrals.dat`.
#
# -------------------------------------------------------------------------------
# Expected input files and directory structure
#
# Current/
# │
# ├── input.dat                  # Contains simulation parameters
# ├── structure_300.lammps-data  # LAMMPS data file with topology
# │
# ├── Temp_300/
# │   ├── seed_1/
# │   │   └── concatenated_atom_info.dat (Trajectory file)
# │   ├── seed_2/
# │   │   └── concatenated_atom_info.dat (Trajectory File)
# │   └── ...
#
#
# -------------------------------------------------------------------------------
#
# input.dat should contain entries like:
#
# units = real
# T = [300]
# seeds = [1,2,3,4] # Random numbers
# timestep = 0.5
#
# -------------------------------------------------------------------------------
# Polymer types supported
#
#   PP : polypropylene
#   PS : polystyrene
#
# The backbone atom IDs are constructed assuming a 40-monomer chain and fixed
# monomer indexing:
#   - PP uses 9 atoms per repeat unit
#   - PS uses 16 atoms per repeat unit
#
# -------------------------------------------------------------------------------
# Output
#
# dihedrals.dat
#   Column 1 : Dihedral angle (degrees)
#   Column 2 : KDE-based probability density
#
# This file can be used later for plotting or comparing dihedral angle
# distributions between different polymer systems.
#
# -------------------------------------------------------------------------------
# Notes
#
# - Only dihedrals formed entirely by backbone atoms are included.
# - The last 10 saved trajectory blocks from each seed are used for analysis.
# - The distribution is smoothed using Gaussian KDE instead of a histogram.
#
###############################################################################


import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

pol_type = input('Polymer type (PP/PS):')

# Read the full LAMMPS data file to extract topology information
filename = 'structure_300.lammps-data'

natoms = ndihedrals = 0
line_dih = None

with open(filename) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if 'atoms' in line:
        natoms = int(line.strip().split()[0])
    if 'dihedrals' in line:
        ndihedrals = int(line.strip().split()[0])
    if line.startswith('Dihedrals'):
        line_dih = i

Dihedrals = []
for line in lines[line_dih+2:line_dih+ndihedrals+2]:
    data = line.strip().split()
    Dihedrals.append(data)

Dihedrals = np.array(Dihedrals).astype(int)
Dihedrals = Dihedrals[:, 2:]


with open('input.dat', 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith("T ="):
            temps_str = line.split('=')[1].strip().strip('[]')
            Temp = list(map(int, temps_str.split(',')))
        elif line.startswith("seeds ="):
            seeds_str = line.split('=')[1].strip().strip('[]')
            seeds = list(map(int, seeds_str.split(',')))

# Define dihedral angle function


if pol_type == 'PP':
    nunit = 9
if pol_type == 'PS':
    nunit = 16


bb_ids = []
id1, id2 = 1, 2
for _ in range(40):
    bb_ids.extend([id1, id2])
    id1 += nunit
    id2 += nunit



def compute_dihedral(p):
    p = np.array(p)
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    angle = np.degrees(np.arctan2( y, x ))
    if angle <0:
        angle+=360
    return angle

def get_dihedrals(dih_ids, dictionary):
    angles = []
    for a, b, c, d in dih_ids:
        if any(i not in bb_ids for i in [a, b, c, d]):
            continue  # Skip this dihedral

        coords = [dictionary[i][-1][1:] for i in [a, b, c, d]]
        angle = compute_dihedral(coords)
        angles.append(angle)
    return angles

# Extract dihedral angles
all_angles = []
for s in seeds:
    filepath = f'./Temp_300/seed_{s}/concatenated_atom_info.dat'
    with open(filepath, 'r') as f:
        lines = f.readlines()

    b_lines = 9 + natoms
    blocks = len(lines) // b_lines

    for i in range(blocks-10, blocks):
        
        start = b_lines * i + 9
        end = b_lines * (i + 1)

        my_dict = {i: [] for i in range(1, natoms + 1)}
        for line in lines[start:end]:
            line1 = line.strip().split()
            id, mol, type = map(int, line1[0:3])
            x, y, z = map(float, line1[9:12])
            my_dict[id].append([type, x, y, z])

        angles = get_dihedrals(Dihedrals, my_dict)
        all_angles.extend(angles)


all_angles = np.array(all_angles)

kde = gaussian_kde(all_angles, bw_method=0.2)  # Lower bandwidth = sharper
xvals = np.linspace(0, 360, 500)
yvals = kde(xvals)

with open('dihedrals.dat', 'w') as f:
    for x, y in zip(xvals, yvals):
        f.write(f"{x:.4f} {y:.6f}\n")
