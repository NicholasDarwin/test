import os
import subprocess
import numpy as np
import pandas as pd
import logging

def calculate_binding_free_energy(complex_pdb, ligand_pdb, receptor_pdb, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Prepare the input files for GROMACS (adjust as needed for your files)
    subprocess.run(['gmx', 'grompp', '-f', 'binding_energy.mdp', '-c', complex_pdb, '-o', 'binding_energy.tpr'], check=True)

    # Run the binding free energy calculation
    subprocess.run(['gmx', 'energy', '-f', 'binding_energy.edr', '-o', 'binding_energy.xvg'], input='1\n', text=True, check=True)

    # Load the results (assuming binding_energy.xvg is generated)
    data = pd.read_csv('binding_energy.xvg', delim_whitespace=True, comment='@', header=None)
    binding_energy = data.iloc[:, 1].mean()  # Assuming the second column contains the energy values

    # Save the results
    with open(os.path.join(output_dir, 'binding_free_energy.txt'), 'w') as f:
        print('Calculated Binding Free Energy: {} kJ/mol'.format(binding_energy))

    return binding_energy

def main():
    # Path to your files
    complex_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\1OXR_docked.pdbqt"
    ligand_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\1OXR_.pdbqt"
    receptor_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data"  # Adjust to the correct receptor file if needed
    output_dir = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\output"  # Define output directory

    # Calculate binding free energy
    binding_energy = calculate_binding_free_energy(complex_pdb, ligand_pdb, receptor_pdb, output_dir)
    print(f'Calculated Binding Free Energy: {binding_energy} kJ/mol')

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.debug("Calculating binding free energy...")

    # Define output directory
    DATA_DIR = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\output"  # Make sure this path is valid

    # Correct function call with output_dir argument
    complex_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\1OXR_docked.pdbqt"
    ligand_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data\1OXR_.pdbqt"
    receptor_pdb = r"C:\Users\Lenovo\Desktop\code\quantum_virus_simulation\data"

    binding_energy = calculate_binding_free_energy(complex_pdb, ligand_pdb, receptor_pdb, DATA_DIR)

    logging.debug(f"Binding free energy: {binding_energy}")
    main()
