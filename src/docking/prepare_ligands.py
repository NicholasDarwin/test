import os
import logging
import subprocess
import time

def prepare_ligand(ligand_pdb_path):
    """
    Prepare a ligand file for docking by converting it to PDBQT format.
    
    :param ligand_pdb_path: Path to the ligand PDB file.
    :return: Path to the prepared PDBQT ligand file.
    """
    logging.debug(f"Starting preparation of ligand: {ligand_pdb_path}")
    
    if not os.path.exists(ligand_pdb_path):
        raise FileNotFoundError(f"Ligand file not found: {ligand_pdb_path}")

    ligand_dir = os.path.dirname(ligand_pdb_path)
    ligand_filename = os.path.basename(ligand_pdb_path)
    ligand_pdbqt_filename = ligand_filename.replace(".pdb", ".pdbqt")
    ligand_pdbqt_path = os.path.join(ligand_dir, ligand_pdbqt_filename)

    prepare_ligand_script = (
        "C:\\Program Files (x86)\\MGLTools-1.5.7\\Lib\\site-packages\\AutoDockTools\\Utilities24\\prepare_ligand4.py"
    )
    
    cmd = (
        f'cd "{ligand_dir}" && '
        f'"C:\\Program Files (x86)\\MGLTools-1.5.7\\python.exe" "{prepare_ligand_script}" '
        f'-l "{ligand_filename}" -o "{ligand_pdbqt_filename}"'
    )
    
    logging.debug(f"Running command: {cmd}")
    
    try:
        # Open a new terminal window and run the command
        if os.name == 'nt':  # Windows
            subprocess.run(f'start cmd /c "{cmd}"', shell=True, check=True)
        else:  # Unix-based systems
            subprocess.run(f'gnome-terminal -- bash -c "{cmd}; exec bash"', shell=True, check=True)
        
        logging.info(f"Ligand prepared: {ligand_pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in ligand preparation: {e.stderr}")
        raise RuntimeError(f"Ligand preparation failed: {e.stderr}")
    
    # Add a short delay to ensure the file system is updated
    time.sleep(1)
    
    return ligand_pdbqt_path

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python prepare_ligands.py <ligand_pdb_path>")
        sys.exit(1)
    
    ligand_pdb_path = sys.argv[1]
    prepared_ligand_path = prepare_ligand(ligand_pdb_path)
    print(f"Prepared ligand saved to: {prepared_ligand_path}")