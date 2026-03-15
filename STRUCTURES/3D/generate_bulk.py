import numpy as np
import os, glob

polymer_type = input('Enter PP or PS for PP/PS ')
tacticity = int(input('Enter 0 or 1 for s/i '))

if tacticity == 0 and polymer_type=='PP':
    z_shift = 14.5
    n_y = 8
    n_z = 2
    z_pos = np.loadtxt('position_info.dat',skiprows=2)[:,6]
    z_min = np.min(z_pos)
    z_max = z_min + (z_shift * n_z)

if tacticity == 1 and polymer_type=='PP':
    z_shift = 6.65
    n_y = 4
    n_z = 4
    z_pos = np.loadtxt('position_info.dat',skiprows=2)[:,6]
    z_min = np.min(z_pos)
    z_max = z_min + (z_shift * n_z)

if tacticity == 0 and polymer_type=='PS':
    circles = 6
    a = 26.26
    z_shift = a/2
    y_shift = (a*1.732)/3
    y_pos = np.loadtxt('position_info.dat',skiprows=2)[:,5]
    z_pos = np.loadtxt('position_info.dat',skiprows=2)[:,6]
    ymin = np.min(y_pos)
    zmin = np.min(z_pos)
    ymax = (circles/2)*y_shift + ymin
    zmax = zmin+a


if tacticity == 1 and polymer_type=='PS':
    n_y = 12
    z_shift = 0
    b = 21.9 
    n_z = 1
    z_pos = np.loadtxt('position_info.dat',skiprows=2)[:,6]
    z_min = np.min(z_pos)
    z_max = z_min+b

############################################### PP/iPS ###################################################

def update_basic_info(n_z):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * n_z} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'zlo':
            new_lines.append(f"{z_min:.6f} {z_max:.6f} zlo zhi\n")
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)

def translate_atoms(n_z, n_y, z_shift):
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]  # Preserve the header
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    for atom in atom_lines:
        all_atoms.append(atom)
    for i in range(1, n_z):
        for atom in atom_lines:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * natoms)
            new_atom[1] = str(int(atom[1]) + n_y*i)
            x, y, z = map(float, atom[4:7])
            z += i * z_shift
            new_atom[4:7] = f"{x:.15f}", f"{y:.15f}", f"{z:.15f}"
            all_atoms.append(new_atom)
    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in all_atoms:
            f.write(" ".join(atom) + "\n")

def replicate_section(file_name, section_name, nc_y, natoms):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    header = lines[:2]  # Preserve the header
    section_lines = [line.strip().split() for line in lines[2:]]
    n_entries = len(section_lines)
    all_entries = []
    for entry in section_lines:
        all_entries.append(entry)
    for i in range(1, n_z):
        for entry in section_lines:
            new_entry = entry.copy()
            new_entry[0] = str(int(entry[0]) + i * n_entries)
            for idx in range(2, len(new_entry)):
                new_entry[idx] = str(int(entry[idx]) + i * natoms)
            all_entries.append(new_entry)

    # Write the combined section data
    combined_file = f"combined_{file_name}"
    with open(combined_file, 'w') as f:
        f.writelines(header)
        for entry in all_entries:
            f.write(" ".join(entry) + "\n")

    print(f"{section_name} combined into {combined_file}")
    return combined_file

def generate_final_lammps_file():
    with open('bulk.lammps-data', 'w') as f:
        # Append the header
        with open('new_basic_info.dat', 'r') as basic:
            f.writelines(basic.readlines())

        # Append atom positions
        f.write("\nAtoms\n\n")
        with open('combined_position_info.dat', 'r') as pos:
            f.writelines(pos.readlines()[2:])  # Skip the original header

        Files = ['bond_info.dat', 'angle_info.dat', 'dihedral_info.dat']
        Types = ['Bonds', 'Angles', 'Dihedrals']

        for section, label in zip(Files, Types):
            f.write(f"\n{label}\n\n")
            with open(f'combined_{section}', 'r') as sec:
                f.writelines(sec.readlines()[2:])  # Skip the original header

    print("\nFinal LAMMPS data file generated: `bulk.lammps-data`")


####################################### sPS ##############################

def update_basic_info_sPS(ymin, ymax, circles):
    with open('basic_info.dat', 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) == 2 and parts[1] in {'atoms', 'bonds', 'angles', 'dihedrals', 'impropers'}:
            new_lines.append(f"{int(parts[0]) * circles} {parts[1]}\n")
        elif len(parts) > 2 and parts[2] == 'ylo':
            new_lines.append(f"{ymin:.6f} {ymax:.6f} ylo yhi\n")
        elif len(parts) > 2 and parts[2] == 'zlo':
            new_lines.append(f"{zmin:.6f} {zmax:.6f} zlo zhi\n")
        else:
            new_lines.append(line)
    with open('new_basic_info.dat', 'w') as f:
        f.writelines(new_lines)


def translate_atoms_sPS(circles, y_shift, z_shift):
    with open('position_info.dat', 'r') as f:
        lines = f.readlines()
    header = lines[:2]  # Preserve the header
    atom_lines = [line.strip().split() for line in lines[2:]]
    natoms = len(atom_lines)
    all_atoms = []
    for atom in atom_lines:
        all_atoms.append(atom)
    for i in range(1, circles):
        for atom in atom_lines:
            new_atom = atom.copy()
            new_atom[0] = str(int(atom[0]) + i * natoms)
            new_atom[1] = str(int(atom[1]) + i)
            x, y, z = map(float, atom[4:7])
            y += i*(y_shift/2)


            if i %2 != 0:
                z += z_shift
            new_atom[4:7] = f"{x:.15f}", f"{y:.15f}", f"{z:.15f}"
            all_atoms.append(new_atom)
    with open('combined_position_info.dat', 'w') as f:
        f.writelines(header)
        for atom in all_atoms:
            f.write(" ".join(atom) + "\n")

def replicate_section_sPS(file_name, section_name, circles, natoms):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    header = lines[:2]  # Preserve the header
    section_lines = [line.strip().split() for line in lines[2:]]
    n_entries = len(section_lines)
    all_entries = []
    for entry in section_lines:
        all_entries.append(entry)
    for i in range(1, circles):
        for entry in section_lines:
            new_entry = entry.copy()
            new_entry[0] = str(int(entry[0]) + i * n_entries)
            for idx in range(2, len(new_entry)):
                new_entry[idx] = str(int(entry[idx]) + i * natoms)
            all_entries.append(new_entry)

    # Write the combined section data
    combined_file = f"combined_{file_name}"
    with open(combined_file, 'w') as f:
        f.writelines(header)
        for entry in all_entries:
            f.write(" ".join(entry) + "\n")

    print(f"{section_name} combined into {combined_file}")
    return combined_file

def generate_final_lammps_file_PS():
    with open('bulk.lammps-data', 'w') as f:
        # Append the header
        with open('new_basic_info.dat', 'r') as basic:
            f.writelines(basic.readlines())

        # Append atom positions
        f.write("\nAtoms\n\n")
        with open('combined_position_info.dat', 'r') as pos:
            f.writelines(pos.readlines()[2:])  # Skip the original header

        Files = ['bond_info.dat', 'angle_info.dat', 'dihedral_info.dat', 'improper_info.dat']
        Types = ['Bonds', 'Angles', 'Dihedrals', 'Impropers']

        for section, label in zip(Files, Types):
            f.write(f"\n{label}\n\n")
            with open(f'combined_{section}', 'r') as sec:
                f.writelines(sec.readlines()[2:])  # Skip the original header

    print("\nFinal LAMMPS data file generated: `bulk.lammps-data`")


if polymer_type == 'PP':
    update_basic_info(n_z)

    with open('position_info.dat', 'r') as f:
        natoms = len(f.readlines()) - 2

    translate_atoms(n_z, n_y, z_shift)

    replicate_section('bond_info.dat', 'Bonds', n_z, natoms)
    replicate_section('angle_info.dat', 'Angles', n_z, natoms)
    replicate_section('dihedral_info.dat', 'Dihedrals', n_z, natoms)
    generate_final_lammps_file()


if polymer_type == 'PS' and tacticity == 0:
    update_basic_info_sPS(ymin, ymax, circles)

    with open('position_info.dat', 'r') as f:
        natoms = len(f.readlines()) - 2

    translate_atoms_sPS(circles, y_shift, z_shift)

    replicate_section_sPS('bond_info.dat', 'Bonds', circles, natoms)
    replicate_section_sPS('angle_info.dat', 'Angles', circles, natoms)
    replicate_section_sPS('dihedral_info.dat', 'Dihedrals', circles, natoms)
    replicate_section_sPS('improper_info.dat', 'Improperss', circles, natoms)
    generate_final_lammps_file_PS()


if polymer_type == 'PS' and tacticity == 1:
    update_basic_info(n_z)
    with open('position_info.dat', 'r') as f:
        natoms = len(f.readlines()) - 2
    translate_atoms(n_z, n_y, z_shift)
    replicate_section('bond_info.dat', 'Bonds', n_z, natoms)
    replicate_section('angle_info.dat', 'Angles', n_z, natoms)
    replicate_section('dihedral_info.dat', 'Dihedrals', n_z, natoms)
    replicate_section('improper_info.dat', 'Impropers', n_z, natoms)
    generate_final_lammps_file_PS()

[os.remove(f) for f in glob.glob('*.dat')]
