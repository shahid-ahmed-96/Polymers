# Polymer Structure Generation Pipeline

This directory contains scripts to build polymer structures for molecular dynamics 
simulations using LAMMPS. The workflow progressively constructs structures from a 
**single polymer chain (1D)** to a **2D crystal** and finally a **3D bulk system**.
All generated structures are LAMMPS data files containing atomic coordinates and topology.

---

## Directory Structure
```
├── 1D/
│   ├── generate/       # Generate initial polymer chain
│   └── make_pbc/       # Convert chain to periodic boundary structure
├── 2D/                 # Build 2D crystal from periodic chain
└── 3D/                 # Build final 3D bulk structure
```

---

## Dependencies

This pipeline depends on **RadonPy** for the initial polymer chain generation (Step 1 only).

- **Repository:** https://github.com/RadonPy/RadonPy
- **Installation:** See RadonPy's README for conda/pip instructions
- **Required for:** `1D/generate/` only

### Citation

If you use this pipeline in your work, please cite RadonPy:

> Y. Hayashi, J. Shiomi, J. Morikawa, R. Yoshida, "RadonPy: Automated Physical 
> Property Calculation using All-atom Classical Molecular Dynamics Simulations for 
> Polymer Informatics," *npj Comput. Mater.* **8**, 222 (2022).  
> https://www.nature.com/articles/s41524-022-00906-4

---

## Workflow Overview

| Step | Directory | Input | Output |
|------|-----------|-------|--------|
| 1 | `1D/generate/` | SMILES string | `eq1.data` |
| 2 | `1D/make_pbc/` | `eq1.data` | `singlechain.lammps-data` |
| 3 | `2D/` | `singlechain.lammps-data` | `crystal_2D.lammps-data` |
| 4 | `3D/` | `crystal_2D.lammps-data` | `bulk.lammps-data` |

---

## Step 1 — Generate Polymer Chain
**Directory:** `1D/generate/`

Generates the initial polymer chain using **RadonPy**.

**Inputs (edit in `get_density.py`):**
- SMILES string (PP or PS)
- Degree of polymerization
- Tacticity

**Output:** `eq1.data` — LAMMPS data file with atomic coordinates, bonds, angles, dihedrals, and impropers.

---

## Step 2 — Convert to Periodic Chain
**Directory:** `1D/make_pbc/`

Aligns the chain along the **x-axis** and introduces vacuum in transverse directions.

1. Copy `eq1.data` from Step 1
2. Run:
```bash
python3 lammps2parsed.py eq1.data
python3 gen_pbc_chain.py
```

**Output:** `singlechain.lammps-data`

---

## Step 3 — Build 2D Crystal
**Directory:** `2D/`

Arranges multiple aligned chains into a **2D sheet**.

1. Copy `singlechain.lammps-data` from Step 2
2. Run:
```bash
python3 lammps2parsed.py singlechain.lammps-data
python3 generate_2D.py
```

**Output:** `crystal_2D.lammps-data`

---

## Step 4 — Build 3D Bulk Structure
**Directory:** `3D/`

Generates the final **3D bulk polymer structure**.

1. Copy `crystal_2D.lammps-data` from Step 3
2. Run:
```bash
python3 lammps2parsed.py crystal_2D.lammps-data
python3 generate_bulk.py
```

**Output:** `bulk.lammps-data`

---

## Notes

- `lammps2parsed.py` must be run before each generation step. Used to parse the .lammps-data files
- All outputs are LAMMPS-compatible structure files ready for simulation
