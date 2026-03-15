import os
from radonpy.core import utils, poly
from radonpy.ff.gaff2_mod import GAFF2_mod
from radonpy.sim import qm, md
from radonpy.sim.preset import eq,tc
from rdkit import Chem
from rdkit.Chem import AllChem
from ase.io import write, read
import numpy as np

smiles = '*CC(*)C'
ter_smiles = '*C'
temp=300
press = 1.0
omp_psi4 = 1
mpi = 1
omp = 0
gpu = 0
mem = 2000
work_dir ='./'
ff = GAFF2_mod()

mol = utils.mol_from_smiles(smiles)        #This is rdkit mol object
print('MOL obtained')

qm.assign_charges(mol, charge='RESP', opt=False, work_dir=work_dir, omp=omp_psi4, memory=mem, log_name='monomer1')
print('Assigned charges')

ter = utils.mol_from_smiles(ter_smiles)
qm.assign_charges(ter, charge='RESP', opt=True, work_dir=work_dir, omp=omp_psi4, memory=mem, log_name='monomer1')

homopoly=poly.polymerize_mols(mol, 40, bond_length=1.5, dihedral=np.pi, random_rot=False, dih_type='monomer',confId=0, tacticity='atactic', atac_ratio=1.0, tac_array=None)
homopoly=poly.terminate_rw(homopoly, ter)
print('Polymerized')

result = ff.ff_assign(homopoly)
if not result:
    print('ERROR: Cannot assign force field')

ac = poly.amorphous_stc(homopoly,1,density=0.05)

print('Created amorphous cell')
eq_str = eq.INC_DEN(ac,  work_dir=work_dir)
ac = eq_str.exec(temp=temp, press=press, mpi=mpi, omp=omp, gpu=gpu) 
