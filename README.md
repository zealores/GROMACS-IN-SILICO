# 🧬 GROMACS Protein Simulation Manual

This repository provides a step-by-step guide for preparing and running a **molecular dynamics (MD) simulation** of a protein system using **GROMACS**.

![GROMACS](https://www.gromacs.org/_static/gmx_logo_blue.png)


![Molecular Dynamics](https://img.shields.io/badge/Molecular-Dynamics-blue)

![GROMACS](https://img.shields.io/badge/GROMACS-2025-blue)

![Protein Simulation](https://img.shields.io/badge/Protein-Simulation-blue)

---

## ⚙️ Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```
pip install gromacs
```

## Protein Simulation 
#### All mdp files can be manually written or download on my github page.
### 1. Generate Topology for protein
#### Converts the PDB file to GROMACS format. Select force field and water model when prompted.
```
gmx pdb2gmx -f protein.pdb -o protein.gro 
# select 6 for force field and 1 for water model
```
- `-f` → input trajectory files
- `-o` → output trajectory files
- `-ignh` → ignore all hydrogens in pdb file

### 2. Define the Simulation Box
#### Creates a cubic box around the protein with 1.0nm between the protein and box edge.
```
gmx editconf -f protein.gro -o newbox.gro -bt cubic -d 1.0
```
- `-bt (box type)` set to cubic, meaning the system will be enclosed in a cubic box.
- `-d (distance)` between the protein and the edge of the box is 1.0 nanometer (nm).

### 3. Solvate the Protein
#### Fills the box with water (SPC216 model). Updates the topology file include water.
```
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```
### 4. Add Ions to Neutralize the System
#### Prepares the input (ions.tpr) and replaces water molecules with Na+ and Cl- ions. Keeps the system electrically neutral.
Note: when running 'genion' choose the group SOL when prompted.

```
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname Cl -neutral
# select 13 to replace SOL with ions
```
- `gmx grompp` → Stands for GROMACS preprocessor. It prepares a binary input file (.tpr) for running simulations or other GROMACS tools.

Options explained:

- `-f ions.mdp` → The MD parameter file. Contains settings for the simulation, like temperature, pressure, number of steps, constraints, and algorithms.

- `-c solv.gro` → The coordinate file of your system (usually solvated protein or molecule). `.gro` is a GROMACS structure format.

- `-p topol.top` → The topology file describing the molecular system, including bonds, angles, charges, and force field parameters.

- `-o ions.tpr` → The output file. This `.tpr` file contains all information needed for the next step (here, adding ions or running the simulation).

### 5. Energy Minimization
#### Minimizes potential energy to remove bad contacts. '-maxwarn 1' bypasses minor warning checks (use cautiously).
```
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -v -deffnm em
```
- `-deffnm em` → Default filename. GROMACS will use this prefix for all output files. For example:

  - `em.gro` → final structure

  - `em.edr` → energy file

  - `em.log` → log file

  - `em.trr` → trajectory file (if applicable)

### 6. Equilibration (NVT Ensemble)
#### Runs a short simulation at constant Number of particles (N), Volume (V), and Temperature (T). Stabilizes the temperature of the system.
```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
```
- `-c` → coordinate
- `-r` → reference

### 7. Run MD Simulation 
```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md
```
