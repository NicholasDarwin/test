import os
import logging
import subprocess

def prepare_receptor(receptor_pdb_path):
    """
    Prepare a receptor file for docking by converting it to PDBQT format.
    
    :param receptor_pdb_path: Path to the receptor PDB file.
    :return: Path to the prepared PDBQT receptor file.
    """
    logging.debug(f"Starting preparation of receptor: {receptor_pdb_path}")
    
    if not os.path.exists(receptor_pdb_path):
        raise FileNotFoundError(f"Receptor file not found: {receptor_pdb_path}")

    receptor_dir = os.path.dirname(receptor_pdb_path)
    receptor_filename = os.path.basename(receptor_pdb_path)
    receptor_pdbqt_filename = receptor_filename.replace(".pdb", ".pdbqt")
    receptor_pdbqt_path = os.path.join(receptor_dir, receptor_pdbqt_filename)

    prepare_receptor_script = (
        "C:\\Program Files (x86)\\MGLTools-1.5.7\\Lib\\site-packages\\AutoDockTools\\Utilities24\\prepare_receptor4.py"
    )
    
    cmd = (
        f'cd "{receptor_dir}" && '
        f'"C:\\Program Files (x86)\\MGLTools-1.5.7\\python.exe" "{prepare_receptor_script}" '
        f'-r "{receptor_filename}" -o "{receptor_pdbqt_filename}"'
    )
    
    logging.debug(f"Running command: {cmd}")
    
    try:
        # Open a new terminal window and run the command
        if os.name == 'nt':  # Windows
            subprocess.run(f'start cmd /c "{cmd}"', shell=True, check=True)
        else:  # Unix-based systems
            subprocess.run(f'gnome-terminal -- bash -c "{cmd}; exec bash"', shell=True, check=True)
        
        logging.info(f"Receptor prepared: {receptor_pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in receptor preparation: {e.stderr}")
        raise RuntimeError(f"Receptor preparation failed: {e.stderr}")
    
    return receptor_pdbqt_path

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python prepare_receptors.py <receptor_pdb_path>")
        sys.exit(1)
    
    receptor_pdb_path = sys.argv[1]
    prepared_receptor_path = prepare_receptor(receptor_pdb_path)
    print(f"Prepared receptor saved to: {prepared_receptor_path}")