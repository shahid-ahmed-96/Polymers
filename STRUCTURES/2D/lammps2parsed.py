import sys

filename = sys.argv[1]

# Initialize variables to store line numbers and counts
natoms = nbonds = nangles = ndihedrals = nimpropers = 0
line_atom = line_vel = line_bond = line_angle = line_dih = line_imp = None

with open(filename) as f:
    lines = f.readlines()

# Parse the file to find section headers and counts
for i, line in enumerate(lines):  # Fixed typo: 'ennumerate' -> 'enumerate'
    if 'atoms' in line:
        natoms = int(line.strip().split()[0])  # Convert to integer
    if 'bonds' in line:
        nbonds = int(line.strip().split()[0])  # Convert to integer
    if 'angles' in line:
        nangles = int(line.strip().split()[0])  # Convert to integer
    if 'dihedrals' in line:
        ndihedrals = int(line.strip().split()[0])  # Convert to integer
    if 'impropers' in line:
        nimpropers = int(line.strip().split()[0])  # Convert to integer

    if line.startswith('Atoms'):
        line_atom = i
    if line.startswith('Velocities'):
        line_vel = i
    if line.startswith('Bonds'):
        line_bond = i
    if line.startswith('Angles'):
        line_angle = i
    if line.startswith('Dihedrals'):
        line_dih = i
    if line.startswith('Impropers'):
        line_imp = i

# Write basic info to 'basic_info.dat'
with open('basic_info.dat', 'w') as f1:
    for line in lines[:line_atom - 1]:  # Corrected slice range
        f1.write(line)

# Write atom info to 'position_info.dat'
with open('position_info.dat', 'w') as f1:
    for line in lines[line_atom :line_atom + natoms+2]:  # Skip header lines
        f1.write(line)

# Write velocity info to 'vel_info.dat' (if velocities exist)
if line_vel is not None:
    with open('vel_info.dat', 'w') as f1:
        for line in lines[line_vel :line_vel + natoms+2]:  # Skip header lines
            f1.write(line)

# Write bond info to 'bond_info.dat'
with open('bond_info.dat', 'w') as f1:
    for line in lines[line_bond :line_bond + nbonds+2]:  # Skip header lines
        f1.write(line)

# Write angle info to 'angle_info.dat'
with open('angle_info.dat', 'w') as f1:
    for line in lines[line_angle:line_angle + nangles+2]:  # Skip header lines
        f1.write(line)

# Write dihedral info to 'dihedral_info.dat'
with open('dihedral_info.dat', 'w') as f1:
    for line in lines[line_dih:line_dih + ndihedrals+2]:  # Skip header lines
        f1.write(line)

# Write improper info to 'improper_info.dat' (if impropers exist)
if line_imp is not None:
    with open('improper_info.dat', 'w') as f1:
        for line in lines[line_imp:line_imp + nimpropers+2]:  # Skip header lines
            f1.write(line)