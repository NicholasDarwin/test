import os
import logging
import subprocess
import platform

def prepare_receptor(receptor_pdb_path):
    """
    Prepare a receptor file for docking by converting it to PDBQT format using OpenBabel.
    
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

    if platform.system() == "Windows":
        # Convert Windows path to WSL path manually
        receptor_pdb_path_wsl = receptor_pdb_path.replace("C:\\", "/mnt/c/").replace("\\", "/")
        receptor_pdbqt_path_wsl = receptor_pdbqt_path.replace("C:\\", "/mnt/c/").replace("\\", "/")
    else:
        receptor_pdb_path_wsl = subprocess.check_output(["wslpath", "-u", receptor_pdb_path]).decode().strip()
        receptor_pdbqt_path_wsl = subprocess.check_output(["wslpath", "-u", receptor_pdbqt_path]).decode().strip()

    # OpenBabel command to convert PDB to PDBQT in WSL
    cmd = (
        f"wsl obabel '{receptor_pdb_path_wsl}' -O '{receptor_pdbqt_path_wsl}'"
    )
    
    logging.debug(f"Running command: {cmd}")
    
    try:
        # Execute the OpenBabel command in WSL
        subprocess.run(cmd, shell=True, check=True)
        
        logging.info(f"Receptor prepared: {receptor_pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in receptor preparation: {e.stderr}")
        raise RuntimeError(f"Receptor preparation failed: {e.stderr}")
    
    return receptor_pdbqt_path

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python prepare_receptor.py <receptor_pdb_path>")
        sys.exit(1)
    
    receptor_pdb_path = sys.argv[1]
    prepared_receptor_path = prepare_receptor(receptor_pdb_path)
    print(f"Prepared receptor saved to: {prepared_receptor_path}")
