import os
import logging
import subprocess
import platform

def prepare_ligand(ligand_pdb_path):
    """
    Prepare a ligand file for docking by converting it to PDBQT format using OpenBabel.
    
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

    if platform.system() == "Windows":
        # Convert Windows path to WSL path manually
        ligand_pdb_path_wsl = ligand_pdb_path.replace("C:\\", "/mnt/c/").replace("\\", "/")
        ligand_pdbqt_path_wsl = ligand_pdbqt_path.replace("C:\\", "/mnt/c/").replace("\\", "/")
    else:
        ligand_pdb_path_wsl = subprocess.check_output(["wslpath", "-u", ligand_pdb_path]).decode().strip()
        ligand_pdbqt_path_wsl = subprocess.check_output(["wslpath", "-u", ligand_pdbqt_path]).decode().strip()

    # OpenBabel command to convert PDB to PDBQT in WSL
    cmd = (
        f"wsl obabel '{ligand_pdb_path_wsl}' -O '{ligand_pdbqt_path_wsl}'"
    )
    
    logging.debug(f"Running command: {cmd}")
    
    try:
        # Execute the OpenBabel command in WSL
        subprocess.run(cmd, shell=True, check=True)
        
        logging.info(f"Ligand prepared: {ligand_pdbqt_path}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in ligand preparation: {e.stderr}")
        raise RuntimeError(f"Ligand preparation failed: {e.stderr}")
    
    # Clean the PDBQT file by removing inappropriate tags
    with open(ligand_pdbqt_path, 'r') as file:
        lines = file.readlines()
    
    with open(ligand_pdbqt_path, 'w') as file:
        for line in lines:
            if not line.startswith("HEADER") and not line.startswith("TITLE"):
                file.write(line)
    
    return ligand_pdbqt_path

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python prepare_ligands.py <ligand_pdb_path>")
        sys.exit(1)
    
    ligand_pdb_path = sys.argv[1]
    prepared_ligand_path = prepare_ligand(ligand_pdb_path)
    print(f"Prepared ligand saved to: {prepared_ligand_path}")
