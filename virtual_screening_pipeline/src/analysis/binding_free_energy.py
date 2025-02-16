import os
import subprocess
import numpy as np
import pandas as pd

def calculate_binding_free_energy(complex_pdb, ligand_pdb, receptor_pdb, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Prepare the input files for GROMACS
    subprocess.run(['gmx', 'grompp', '-f', 'binding_energy.mdp', '-c', complex_pdb, '-o', 'binding_energy.tpr'], check=True)

    # Run the binding free energy calculation
    subprocess.run(['gmx', 'energy', '-f', 'binding_energy.edr', '-o', 'binding_energy.xvg'], input='1\n', text=True, check=True)

    # Load the results
    data = pd.read_csv('binding_energy.xvg', delim_whitespace=True, comment='@', header=None)
    binding_energy = data.iloc[:, 1].mean()  # Assuming the second column contains the energy values

    # Save the results
    with open(os.path.join(output_dir, 'binding_free_energy.txt'), 'w') as f:
        f.write(f'Binding Free Energy: {binding_energy} kJ/mol\n')

    return binding_energy

def main():
    complex_pdb = 'path/to/complex.pdb'
    ligand_pdb = 'path/to/ligand.pdb'
    receptor_pdb = 'path/to/receptor.pdb'
    output_dir = 'output'

    binding_energy = calculate_binding_free_energy(complex_pdb, ligand_pdb, receptor_pdb, output_dir)
    print(f'Calculated Binding Free Energy: {binding_energy} kJ/mol')

if __name__ == "__main__":
    main()