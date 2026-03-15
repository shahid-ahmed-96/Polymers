import numpy as np
import os, glob

polymer_type = input('Enter polymer type (PP, PS): ')
tacticity = int(input('0 or 1 for syndio / iso: '))

################################## sPP START #########################################

def update_basic_info_sPP(nc, ymin, a):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * nc} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'ylo':
            ymax = ymin+nc*(0.5*a)
            new_lines.append(f"{ymin} {ymax} ylo yhi\n")            
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)


def syndiotactic_PP(a, b):
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    chain = atom_lines
        
    for i in range(nc):
        y_shift = i * 0.5 * a
        z_shift = 0.5 * b if i % 2 != 0 else 0.0

        for atom in chain:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * natoms)
            new_atom[1] = str(i + 1)
            x, y, z = map(float, atom[4:7])
            new_atom[4:7] = f"{x:.15f}", f"{y + y_shift:.15f}", f"{z + z_shift:.15f}"
            all_atoms.append(new_atom)

    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in all_atoms:
            f.write(" ".join(atom) + "\n")

################################## sPP END ###########################################

################################## iPP START #########################################

def update_basic_info_iPP(nc, ymin, a, b):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * nc} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'ylo':
            ymax = ymin+nc*(0.25*b)
            new_lines.append(f"{ymin} {ymax} ylo yhi\n")
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)

def mirrored(coords):
    Y = coords[:, 5].astype(float)
    y1 = np.mean(Y)
    Y = -Y
    y2 = np.mean(Y)
    shift = y1-y2
    Y = Y+shift
    return Y

def isotactic_PP(Lx, nc, a, b):
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    L_chain = atom_lines
    L_chain1 = np.array(L_chain)
    r_chain_y = mirrored(L_chain1)
    R_chain = []
    shift = 0.25*b

    for idx in range(len(L_chain)):
        atom = L_chain[idx]
        new_atom = atom.copy()
        new_atom[0] = str(int(atom[0]) + natoms)
        new_atom[1] = str(int(atom[1]) + 1)

        new_atom[5] = str(r_chain_y[idx] + shift)
        R_chain.append(new_atom)

    for i in range(nc//2):
        y_shift = i * 0.5 * b
        z_shift = 0.5 * a if i % 2 != 0 else 0.0

        for atom in L_chain:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * 2 * natoms)
            new_atom[1] = str(2*i + 1)
            x, y, z = map(float, atom[4:7])
            #print(new_atom)
            new_atom[4:7] = f"{x:.15f}", f"{y + y_shift:.15f}", f"{z - z_shift:.15f}"
            #print(new_atom)
            all_atoms.append(new_atom)

        for atom in R_chain:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * 2 * natoms)
            new_atom[1] = str(2*i + 2)
            x, y, z = map(float, atom[4:7])
            #print(new_atom)
            new_atom[4:7] = f"{x:.15f}", f"{y - y_shift:.15f}", f"{z - z_shift:.15f}"
            #print(new_atom)
            all_atoms.append(new_atom)

    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in all_atoms:
            f.write(" ".join(atom) + "\n")

################################## iPP END ###########################################

################################## sPS START #########################################

def update_basic_info_sPS(nc):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * nc} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'ylo':
            new_lines.append(f"{-25} {25} ylo yhi\n")
        elif len(parts) > 2 and parts[2] == 'zlo':
            new_lines.append(f"{-25} {25} zlo zhi\n")
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)


def modify_chain(chain, translate):
    chain_ = np.array(chain)
    X = chain_[:,4].astype(float)
    Y = chain_[:,5].astype(float)
    Z = chain_[:,6].astype(float)

    x_0 = X[0]
    y_0 = Y[0]
    z_0 = Z[0]

    X = X - x_0
    Y = Y - y_0
    Z = Z - z_0

    positions = np.column_stack((X, Y, Z))

    direction_y = positions[1][1] - positions[0][1]
    direction_z = positions[1][2] - positions[0][2]
    theta = np.arctan2(direction_z, direction_y)

    cos_theta = np.cos(-theta)
    sin_theta = np.sin(-theta)
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, cos_theta, sin_theta],
        [0, -sin_theta, cos_theta]
    ])

    rotated_positions = np.dot(positions, rotation_matrix.T)

    new_chain = []
    for i in range(len(chain)):
        new_atom = chain[i].copy()
        new_atom[4] = str(rotated_positions[i][0])
        new_atom[5] = str(rotated_positions[i][1] + translate[1])
        new_atom[6] = str(rotated_positions[i][2])
        new_chain.append(new_atom)
    return new_chain


def change_chain(chain, rotation):
    chain_ = np.array(chain)

    X = chain_[:,4].astype(float)
    Y = chain_[:,5].astype(float)
    Z = chain_[:,6].astype(float)

    positions = np.column_stack((X,Y,Z))

    theta = np.radians(rotation)  # Convert 120° to radians
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    rotation_matrix = np.array([
        [1, 0, 0],
        [0, cos_theta, -sin_theta],
        [0, sin_theta, cos_theta]
    ])

    rotated_positions = np.dot(positions, rotation_matrix.T)
    natoms = len(rotated_positions)

    new_chain = []
    for i in range(len(chain)):
        chain[i][0] = str( int(chain[i][0]) + natoms)
        chain[i][1] = str( int(chain[i][1]) + 1)
        chain[i][4] = str( rotated_positions[i][0])
        chain[i][5] = str( rotated_positions[i][1])
        chain[i][6] = str( rotated_positions[i][2])
        new_chain.append(chain[i])
    return new_chain

def syndiotactic_PS():
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    chain = atom_lines
    new_chain = modify_chain(chain, [0,5,0]) # Adjust y translate

    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in new_chain:
            f.write(" ".join(atom) + "\n")

    chain_2 = change_chain(new_chain, 120)
    with open('combined_position_info.dat', 'a') as f:
        for atom in chain_2:
            f.write(" ".join(atom) + "\n")

    chain_3 = change_chain(chain_2, 120)
    with open('combined_position_info.dat', 'a') as f:
        for atom in chain_3:
            f.write(" ".join(atom) + "\n")

################################## sPS END ###########################################

################################## iPS START ###########################################

def update_basic_info_iPS(nc, ymin, zmin, a, b):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * nc} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'ylo':
            ymax = ymin+a
            new_lines.append(f"{ymin} {ymax} ylo yhi\n")
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)


def isotactic_PS(a, b):
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    L_chain = atom_lines
    R_chain = []
    Z_L = []

    for atom in L_chain:
        z = float(atom[6])
        Z_L.append(z)
    shift = ((2*b)/3) + 2*np.mean(np.array(Z_L))

    for atom in L_chain:
        x, y, z = map(float, atom[4:7])
        z_mirror = - z+shift
        new_atom = atom.copy()
        new_atom[0] = str(int(atom[0]) + 6*natoms)  # Atom ID offset
        new_atom[1] = str(int(atom[1]) + 6)       # Molecule ID +1 (R-chain)
        new_atom[4:7] = f"{x:.15f}", f"{y:.15f}", f"{z_mirror:.15f}"
        R_chain.append(new_atom)

    for i in range(6):
        y_shift = i * (a/6)
        z_shift = b/2 if i % 2 != 0 else 0.0

        for atom in L_chain:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * natoms)
            new_atom[1] = str(i + 1)
            x, y, z = map(float, atom[4:7])
            new_atom[4:7] = f"{x:.15f}", f"{y + y_shift:.15f}", f"{z + z_shift:.15f}"
            all_atoms.append(new_atom)

        for atom in R_chain:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + (i* natoms))
            new_atom[1] = str(int(atom[1])+i)
            x, y, z = map(float, atom[4:7])
            new_atom[4:7] = f"{x:.15f}", f"{y + y_shift:.15f}", f"{z - z_shift:.15f}"
            all_atoms.append(new_atom)

    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in all_atoms:
            f.write(" ".join(atom) + "\n")

################################## iPS END ###########################################

def replicate_section(file_name, section_name, nc, natoms):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    header = lines[:2]
    section_lines = [line.strip().split() for line in lines[2:]]
    n_entries = len(section_lines)
    all_entries = []
    for entry in section_lines:
        all_entries.append(entry)
    for i in range(1, nc):
        for entry in section_lines:
            new_entry = entry.copy()
            new_entry[0] = str(int(entry[0]) + i * n_entries)

            for idx in range(2, len(new_entry)):
                new_entry[idx] = str(int(entry[idx]) + i * natoms)

            all_entries.append(new_entry)

    combined_file = f"combined_{file_name}"
    with open(combined_file, 'w') as f:
        f.writelines(header)
        for entry in all_entries:
            f.write(" ".join(entry) + "\n")

    print(f"{section_name} combined into {combined_file}")
    return combined_file

def generate_final_lammps_file():
    with open('crystal_2D.lammps-data', 'w') as f:
        with open('new_basic_info.dat', 'r') as basic:
            f.writelines(basic.readlines())

        f.write("\nAtoms\n\n")
        with open('combined_position_info.dat', 'r') as pos:
            f.writelines(pos.readlines()[2:])

        Files = ['bond_info.dat', 'angle_info.dat', 'dihedral_info.dat']
        Types = ['Bonds', 'Angles', 'Dihedrals']
        if polymer_type == 'PS':
            Files.append('improper_info.dat')
            Types.append('Impropers')

        for section, label in zip(Files, Types):
            f.write(f"\n{label}\n\n")
            with open(f'combined_{section}', 'r') as sec:
                f.writelines(sec.readlines()[2:])
    print("\nFinal LAMMPS data file generated: `crystal_2D.lammps-data`")

with open('position_info.dat', 'r') as f:
    lines = f.readlines()
    natoms = len(lines) - 2
    pos_lines = [line.strip().split() for line in lines[2:]]
    y_coords = np.array([float(line[5]) for line in pos_lines])
    z_coords = np.array([float(line[6]) for line in pos_lines])
    ymin = np.min(y_coords)
    zmin = np.min(z_coords)

with open('basic_info.dat', 'r') as f:
    lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 2 and parts[2] == 'xlo':
            Lx = float(parts[1]) - float(parts[0])


if polymer_type == 'PP' and tacticity == 1:
    a = 6.65
    b = 20.96
    nc = 4
    update_basic_info_iPP(nc, ymin, a, b)
    isotactic_PP(Lx, nc, a, b)

if polymer_type == 'PP' and tacticity == 0:
    a = 5.6
    b = 14.5
    nc = 8
    update_basic_info_sPP(nc, ymin, a)
    syndiotactic_PP(a,b)

if polymer_type == 'PS' and tacticity == 0:
    syndiotactic_PS()
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
        natoms = len(lines) - 2
    with open('combined_position_info.dat', 'r') as f:
        lines = f.readlines()
        natoms_total = len(lines) - 2
    nc = int(natoms_total / natoms)
    update_basic_info_sPS(nc)


if polymer_type == 'PS' and tacticity == 1:
    b = 21.9
    a = 21.9 * (3**0.5)
    nc = 12
    update_basic_info_iPS(nc, ymin, zmin, a, b)
    isotactic_PS(a,b)

replicate_section('bond_info.dat', 'Bonds', nc, natoms)
replicate_section('angle_info.dat', 'Angles', nc, natoms)
replicate_section('dihedral_info.dat', 'Dihedrals', nc, natoms)
if polymer_type == 'PS':
    replicate_section('improper_info.dat', 'Impropers', nc, natoms)

generate_final_lammps_file()
[os.remove(f) for f in glob.glob('*.dat')]
