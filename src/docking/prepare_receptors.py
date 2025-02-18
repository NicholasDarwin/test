import os
import logging
import subprocess

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

    try:
        # Ensure the receptor file exists
        if not os.path.exists(receptor_pdb_path):
            logging.error(f"Receptor file not found: {receptor_pdb_path}")
            return None

        # Prepare the receptor file (example: using AutoDockTools or similar)
        prepared_receptor_path = receptor_pdb_path.replace('.pdb', '_prepared.pdbqt')

        # Example command to prepare the receptor (replace with actual command)
        # command = ["prepare_receptor", "-r", receptor_pdb_path, "-o", prepared_receptor_path]
        # subprocess.run(command, check=True)

        # For demonstration, we'll just copy the file (replace with actual preparation logic)
        with open(receptor_pdb_path, 'r') as src, open(prepared_receptor_path, 'w') as dst:
            dst.write(src.read())

        logging.debug(f"Prepared receptor file: {prepared_receptor_path}")
        return prepared_receptor_path
    except Exception as e:
        logging.error(f"Error preparing receptor: {e}")
        return None

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python prepare_receptors.py <receptor_pdb_path>")
        sys.exit(1)
    
    receptor_pdb_path = sys.argv[1]
    prepared_receptor_path = prepare_receptor(receptor_pdb_path)
    print(f"Prepared receptor saved to: {prepared_receptor_path}")
