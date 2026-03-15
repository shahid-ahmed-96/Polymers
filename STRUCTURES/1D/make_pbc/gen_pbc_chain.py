import numpy as np
import sklearn
import os, glob
import ase
import ase.io
from ase.io import read, write
from ase.visualize import view
from sklearn.decomposition import PCA
from scipy.spatial.transform import Rotation as R
import shutil

shutil.copy("eq1.data", "eq1.lammps-data")

chain = input('Name of polymer (PE, PP, PS): ').strip()

if chain == 'PP':
    monomer = 9
else:
    monomer = 16

print('No of atoms in monomer: ',monomer)

monomer = monomer*2


def periodic_params(data, natoms):
    # bond does not mean only bond info, can have sngles, dihedrals, impropers
    bond_join = [bond for bond in data if np.min(bond[2:]) <= monomer and np.max(bond[2:]) > monomer]
    new_bonds = []
    current_bond_id = len(data) + 1
    for bond in bond_join:
        bond_type = bond[1]  # Keep bond type unchanged
        new_bond = [current_bond_id, bond_type]
        for atom in bond[2:]:
            if atom <= monomer:
                new_atom = natoms - (monomer - atom)
            else:
                new_atom = atom - monomer
            new_bond.append(new_atom)

        new_bonds.append(new_bond)
        current_bond_id += 1  
    new_bonds = np.array(new_bonds, dtype=int)
    if len(new_bonds) > 0:
        data = np.vstack((data, new_bonds))
    return data


atoms = read('eq1.lammps-data', atom_style = 'full')
atom_nos = atoms.get_atomic_numbers()
atom_info = np.loadtxt('position_info.dat', skiprows=2)
ids = atom_info[:, 0].astype(int)
keys = atom_info[:, 2].astype(int)
atom_dict = dict(zip(ids, keys))
c_ids = [atom.index + 1 for atom in atoms if atom.symbol == 'C']
c_ids = [c for c in c_ids if atom_dict[c] == 1]

positions = atoms.get_positions()
n = len(positions)
del_id = [0, 1, 2, 3, n-1, n-2, n-3, n-4]
del atoms[[atom.index for atom in atoms if atom.index in del_id]]
positions = atoms.get_positions()

centroid = np.mean(positions, axis=0)
centered_positions = positions - centroid
pca = PCA(n_components=3)
pca.fit(centered_positions)
principal_axes = pca.components_

target_axis = np.array([1, 0, 0])  # x-axis
rotation_axis = np.cross(principal_axes[0], target_axis)

norm = np.linalg.norm(rotation_axis)
if norm > 1e-6:
    rotation_axis /= norm
else:
    rotation_axis = np.array([1, 0, 0])



angle = np.arccos(np.dot(principal_axes[0], target_axis))

def rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
        [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
        [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]
    ])

rotation_mat = rotation_matrix(rotation_axis, angle)
rotated_positions = np.dot(centered_positions, rotation_mat.T)
atoms.set_positions(rotated_positions)
positions = atoms.get_positions()
x = positions[:,0]
y = positions[:,1]
z = positions[:,2]
span_x = np.max(x)-np.min(x)
span_y = np.max(y)-np.min(y)
span_z = np.max(z)-np.min(z)
atoms.set_cell([span_x+1,span_y+20,span_z+20])
atoms.center()
positions = atoms.get_positions()
natoms = len(positions)

#view(atoms)

del_id = [1+d for d in del_id]
c_ids = [c for c in c_ids if c not in del_id]

# BOND_DATA
bond_data = np.loadtxt('bond_info.dat',skiprows=2).astype(int)
linenos_1 = [i for i in range(len(bond_data)) if bond_data[i, -2] in del_id]
linenos_2 = [i for i in range(len(bond_data)) if bond_data[i, -1] in del_id]
unique_linenos = sorted(set(linenos_1) | set(linenos_2))
bond_data = np.delete(bond_data, unique_linenos, axis=0)
bond_data = bond_data.astype(int)

bond_data[:,0] = np.arange(1, len(bond_data) + 1)
bond_data[:,2] = bond_data[:,2]-4
bond_data[:,3] = bond_data[:,3]-4

# ANGLE_DATA
angle_data = np.loadtxt('angle_info.dat',skiprows=2).astype(int)
linenos_1 = [i for i in range(len(angle_data)) if angle_data[i, -3] in del_id]
linenos_2 = [i for i in range(len(angle_data)) if angle_data[i, -2] in del_id]
linenos_3 = [i for i in range(len(angle_data)) if angle_data[i, -1] in del_id]
unique_linenos = sorted(set(linenos_1) | set(linenos_2) | set(linenos_3))
angle_data = np.delete(angle_data, unique_linenos, axis=0)
angle_data = angle_data.astype(int)

angle_data[:,0] = np.arange(1, len(angle_data) + 1)
angle_data[:,2] = angle_data[:,2]-4
angle_data[:,3] = angle_data[:,3]-4
angle_data[:,4] = angle_data[:,4]-4

angle_data = periodic_params(angle_data, natoms)

np.savetxt('new_angles.dat', angle_data, fmt="%s", delimiter=" ", header="Angles\n", comments='')

# DIHEDRAL_DATA
dihedral_data = np.loadtxt('dihedral_info.dat',skiprows=2).astype(int)
linenos_1 = [i for i in range(len(dihedral_data)) if dihedral_data[i, -4] in del_id]
linenos_2 = [i for i in range(len(dihedral_data)) if dihedral_data[i, -3] in del_id]
linenos_3 = [i for i in range(len(dihedral_data)) if dihedral_data[i, -2] in del_id]
linenos_4 = [i for i in range(len(dihedral_data)) if dihedral_data[i, -1] in del_id]
unique_linenos = sorted(set(linenos_1) | set(linenos_2) | set(linenos_3) | set(linenos_4))
dihedral_data = np.delete(dihedral_data, unique_linenos, axis=0)
dihedral_data = dihedral_data.astype(int)

dihedral_data[:,0] = np.arange(1, len(dihedral_data) + 1)
dihedral_data[:,2] = dihedral_data[:,2]-4
dihedral_data[:,3] = dihedral_data[:,3]-4
dihedral_data[:,4] = dihedral_data[:,4]-4
dihedral_data[:,5] = dihedral_data[:,5]-4

dihedral_data = periodic_params(dihedral_data, natoms)

np.savetxt('new_dihedrals.dat', dihedral_data, fmt="%s", delimiter=" ", header="Dihedrals\n", comments='')


# IMPROPER_DATA
if os.path.exists('improper_info.dat'):
    improper_data = np.loadtxt('improper_info.dat', skiprows=2).astype(int)
    linenos_1 = [i for i in range(len(improper_data)) if improper_data[i, -4] in del_id]
    linenos_2 = [i for i in range(len(improper_data)) if improper_data[i, -3] in del_id]
    linenos_3 = [i for i in range(len(improper_data)) if improper_data[i, -2] in del_id]
    linenos_4 = [i for i in range(len(improper_data)) if improper_data[i, -1] in del_id]
    unique_linenos = sorted(set(linenos_1) | set(linenos_2) | set(linenos_3) | set(linenos_4))
    improper_data = np.delete(improper_data, unique_linenos, axis=0)
    improper_data = improper_data.astype(int)
    
    improper_data[:,0] = np.arange(1, len(improper_data) + 1)
    improper_data[:,2] = improper_data[:,2]-4
    improper_data[:,3] = improper_data[:,3]-4
    improper_data[:,4] = improper_data[:,4]-4
    improper_data[:,5] = improper_data[:,5]-4
    improper_data = periodic_params(improper_data, natoms)
    
    np.savetxt('new_impropers.dat', improper_data, fmt="%s", delimiter=" ", header="Impropers\n", comments='')

def count(index):
    return np.sum((bond_data[:, -2] == index) | (bond_data[:, -1] == index))

C=[]
for c in c_ids:
    v = count(c-4)    
    if v<4:
        C.append(c-5) #To ensure atom indices

C0 = positions[C[0]]
C1 = positions[C[1]]
vec = C1 - C0
rotation_axis_x = np.array([1, 0, 0])
angle_x = np.arctan2(vec[2], vec[1])
rot_x = R.from_rotvec(-angle_x * rotation_axis_x)
positions = rot_x.apply(positions - C0) + C0

C0 = positions[C[0]]
C1 = positions[C[1]]
vec = C1 - C0
rotation_axis_z = np.array([0, 0, 1])
angle_z = np.arctan2(vec[1], vec[0])
rot_z = R.from_rotvec(-angle_z * rotation_axis_z)
positions = rot_z.apply(positions - C0) + C0



span_x = np.max(positions[:,0])-np.min(positions[:,0])
res = positions[C[0]][0]+span_x-positions[C[1]][0]
lx = span_x - res + 1.5

atoms.set_positions(positions)
atoms.set_cell([lx, 50, 50])
atoms.center()
positions = atoms.get_positions()
cell = atoms.get_cell()
#view(atoms)

bond_data = periodic_params(bond_data, natoms)
np.savetxt('new_bonds.dat', bond_data, fmt="%s", delimiter=" ", header="Bonds\n", comments='')

# INFO_DATA
input_file = "basic_info.dat"
output_file = "new_basics.dat"
with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        if 'xlo' in line:            
            line = f"0.000000 {lx:.6f} xlo xhi\n"
        if 'ylo' in line:
            line = "0.000000 50.000000 ylo yhi \n"
        if 'zlo' in line:
            line = "0.000000 50.000000 zlo zhi \n"
        if 'atoms' in line:            
            line = f"{int(len(positions))} atoms\n"
        if 'bonds' in line:            
            line = f"{int(len(bond_data))} bonds\n"
        if 'angles' in line:            
            line = f"{int(len(angle_data))} angles\n"
        if 'impropers' in line:            
            line = f"{int(len(improper_data))} impropers\n"
        if 'dihedrals' in line:            
            line = f"{int(len(dihedral_data))} dihedrals\n"    
        f_out.write(line)


ase.io.write('trial.lammps-data', atoms, format='lammps-data', atom_style='full')

natoms = 0
line_atom = None
filename = 'trial.lammps-data'
with open(filename) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if 'atoms' in line:
        natoms = int(line.strip().split()[0])  # Convert to integer
    if line.startswith('Atoms'):
        line_atom = i
        
with open('new_positions.dat', 'w') as f1:
    for line in lines[line_atom :line_atom + natoms+2]:
        line = line.strip().split()
        if len(line)>3: #just check that its not heading line
            line[2] = str(atom_dict[int(line[0])+4])
        f1.write(' '.join(line) + '\n')


files = [
    "new_basics.dat",
    "new_positions.dat",
    "new_bonds.dat",
    "new_angles.dat",
    "new_dihedrals.dat"
]

if os.path.exists("new_impropers.dat"):
    files.append("new_impropers.dat")

output_file = "single_chain.lammps-data"

with open(output_file, "w") as f_out:
    for file in files:
        if os.path.exists(file):
            with open(file, "r") as f_in:
                f_out.write(f_in.read())  
                f_out.write("\n")

[os.remove(f) for f in glob.glob('*.dat')]