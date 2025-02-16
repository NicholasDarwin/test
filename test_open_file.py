import subprocess
import os

def calculate_binding_free_energy(docked_pdbqt, data_dir):
    # Step 1: Convert the PDBQT file to PDB format using Open Babel
    pdb_file = os.path.join(data_dir, 'docked_complex.pdb')
    subprocess.run(['wsl', 'obabel', docked_pdbqt, '-O', pdb_file], check=True)

    # Step 2: Generate the topology files using GROMACS (protein and ligand preparation)
    # Assuming the receptor protein and ligand files are available
    protein_pdb = os.path.join(data_dir, 'protein.pdb')  # Modify this with the actual protein PDB file
    ligand_pdb = pdb_file  # The converted docked complex file will be treated as ligand
    subprocess.run(['wsl', 'gmx', 'pdb2gmx', '-f', protein_pdb, '-o', 'protein_processed.gro'], check=True)
    
    # Step 3: Prepare the complex (protein-ligand) for GROMACS
    # Here we're assuming the ligand is already in the proper format (from docked PDBQT)
    complex_pdb = os.path.join(data_dir, 'complex.pdb')
    subprocess.run(['wsl', 'gmx', 'editconf', '-f', ligand_pdb, '-o', complex_pdb], check=True)

    # Step 4: Energy minimization of the docked complex
    minim_mdp = os.path.join(data_dir, 'minim.mdp')  # Make sure the minim.mdp file is present
    subprocess.run(['wsl', 'gmx', 'grompp', '-f', minim_mdp, '-c', complex_pdb, '-o', 'minimized.tpr'], check=True)
    subprocess.run(['wsl', 'gmx', 'mdrun', '-v', '-deffnm', 'minimized'], check=True)

    # Step 5: Perform energy calculation
    binding_energy_mdp = os.path.join(data_dir, 'binding_energy.mdp')  # Ensure this file exists
    subprocess.run(['wsl', 'gmx', 'grompp', '-f', binding_energy_mdp, '-c', complex_pdb, '-o', 'binding_energy.tpr'], check=True)
    subprocess.run(['wsl', 'gmx', 'mdrun', '-v', '-deffnm', 'binding_energy'], check=True)
    
    # Step 6: Extract binding energy and other data from the output
    energy_file = 'binding_energy.edr'
    subprocess.run(['wsl', 'gmx', 'energy', '-f', energy_file, '-o', 'binding_energy.xvg'], check=True)
    
    # Parse the output .xvg file to extract data
    with open('binding_energy.xvg', 'r') as f:
        binding_data = f.readlines()
    
    # Step 7: Extract binding energy and other results
    binding_energy = None
    for line in binding_data:
        if 'G96' in line:  # Identify the line containing the binding energy data
            binding_energy = line.strip()
            break
    
    return binding_energy, binding_data

# Usage example:
docked_pdbqt = 'data1OXR_docked.pdbqt'  # Example input file path
data_dir = '/mnt/c/Users/Lenovo/Desktop/code/quantum_virus_simulation/data/'  # Your working directory
binding_energy, binding_data = calculate_binding_free_energy(docked_pdbqt, data_dir)

print(f"Binding Energy: {binding_energy}")
print(f"Other Data: {binding_data[:5]}")  # Print the first few lines of other data
