import os
import subprocess
import logging
import shutil

logging.basicConfig(level=logging.DEBUG)

def convert_path_to_wsl(path):
    """
    Convert a Windows path to a WSL path.
    
    :param path: Windows path.
    :return: WSL path.
    """
    if isinstance(path, str):
        return path.replace('C:\\', '/mnt/c/').replace('\\', '/')
    else:
        raise ValueError(f"Path must be a string, but got {type(path)}")

def run_gromacs_simulation(protein_file, ligand_file):
    """
    Run a series of GROMACS commands to perform molecular dynamics (MD) simulations.
    
    :param protein_file: Path to the protein file.
    :param ligand_file: Path to the ligand file.
    :return: Output of the GROMACS commands.
    """
    protein_file_wsl = convert_path_to_wsl(protein_file)
    ligand_file_wsl = convert_path_to_wsl(ligand_file)
    output_dir = "/mnt/c/Users/Lenovo/Desktop/code/quantum_virus_simulation/data/gromacs_sim"
    mdp_dir = "/mnt/c/Users/Lenovo/Desktop/code/quantum_virus_simulation/data/gromacs_mdp"

    # Clear the data directory before starting a new simulation
    data_dir = "C:\\Users\\Lenovo\\Desktop\\code\\quantum_virus_simulation\\data\\gromacs_sim"
    for filename in os.listdir(data_dir):
        file_path = os.path.join(data_dir, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            logging.error(f"Failed to delete {file_path}. Reason: {e}")

    commands = [
        f"gmx pdb2gmx -ignh -f {protein_file_wsl} -o {output_dir}/processed.gro -p {output_dir}/protein.top -i {output_dir}/posre.itp -water spce",  # Prepare protein structure in GROMACS format
        f"gmx editconf -f {output_dir}/processed.gro -o {output_dir}/protein_box.gro -c -d 1.0 -bt cubic",  # Define the simulation box and center the protein
        f"gmx solvate -cp {output_dir}/protein_box.gro -cs spc216.gro -o {output_dir}/protein_solvated.gro -p {output_dir}/protein.top",  # Add water molecules to the box
        f"gmx grompp -f {mdp_dir}/ions.mdp -c {output_dir}/protein_solvated.gro -p {output_dir}/protein.top -o {output_dir}/ions.tpr",  # Prepare the system for ion addition
        f"gmx genion -s {output_dir}/ions.tpr -o {output_dir}/protein_ions.gro -p {output_dir}/protein.top -pname NA -nname CL -neutral",  # Add ions to neutralize the system
        f"gmx make_ndx -f {output_dir}/protein_ions.gro -o {output_dir}/index.ndx",  # Generate index file
        f"gmx grompp -f {mdp_dir}/minim.mdp -c {output_dir}/protein_ions.gro -p {output_dir}/protein.top -n {output_dir}/index.ndx -o {output_dir}/em.tpr",  # Prepare the system for energy minimization
        f"gmx mdrun -v -deffnm {output_dir}/em",  # Perform energy minimization
        f"gmx grompp -f {mdp_dir}/nvt.mdp -c {output_dir}/em.gro -r {output_dir}/em.gro -p {output_dir}/protein.top -n {output_dir}/index.ndx -o {output_dir}/nvt.tpr",  # Prepare the system for NVT equilibration
        f"gmx mdrun -v -deffnm {output_dir}/nvt",  # Run NVT equilibration
        f"gmx grompp -f {mdp_dir}/md.mdp -c {output_dir}/nvt.gro -r {output_dir}/nvt.gro -p {output_dir}/protein.top -n {output_dir}/index.ndx -o {output_dir}/md.tpr -maxwarn 1",  # Prepare the system for MD simulation
        f"gmx mdrun -v -deffnm {output_dir}/md",  # Run the MD simulation
        f"gmx pdb2gmx -ignh -f {ligand_file_wsl} -o {output_dir}/ligand.gro -ff charmm27 -water spce",  # Convert the ligand's PDB file to GROMACS format
        f"gmx insert-molecules -f {output_dir}/protein_solvated.gro -ci {output_dir}/ligand.gro -o {output_dir}/protein_ligand.gro -p {output_dir}/protein.top -nmol 1"  # Insert the ligand into the system
    ]

    for command in commands:
        logging.debug(f"Running GROMACS command: {command}")
        try:
            if 'pdb2gmx' in command and 'ligand' not in command:  # Only for the pdb2gmx command for protein
                result = subprocess.run(f"wsl {command}", input="5\n", capture_output=True, text=True, shell=True, check=True)
            elif 'genion' in command:  # Only for the genion command
                result = subprocess.run(f"wsl {command}", input="13\n", capture_output=True, text=True, shell=True)
            elif 'make_ndx' in command:  # Only for the make_ndx command
                result = subprocess.run(f"wsl {command}", input="13\nq\n", capture_output=True, text=True, shell=True)
            elif 'pdb2gmx' in command and 'ligand' in command:  # Only for the pdb2gmx command for ligand
                result = subprocess.run(f"wsl {command}", input="5\n1\n", capture_output=True, text=True, shell=True)
            else:
                result = subprocess.run(f"wsl {command}", capture_output=True, text=True, shell=True)
        
            if result.returncode != 0:
                logging.error(f"Failed to run GROMACS command: {result.stderr}")
                return None
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running GROMACS command: {e.stderr}")
            return None
        
    return "GROMACS simulation completed successfully."

# Example usage:
if __name__ == "__main__":
    protein_file = "C:\\path\\to\\protein.pdb"
    ligand_file = "C:\\path\\to\\ligand.pdb"
    output_dir = "output"

    result = run_gromacs_simulation(protein_file, ligand_file)
    print(result)