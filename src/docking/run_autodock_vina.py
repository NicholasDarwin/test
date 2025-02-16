import subprocess
import os
import logging

logging.basicConfig(level=logging.DEBUG)

def convert_path_to_wsl(path):
    """
    Convert a Windows path or a list of paths to WSL paths.
    
    :param path: Windows path or list of paths.
    :return: WSL path or list of WSL paths.
    """
    if isinstance(path, str):
        return path.replace('C:\\', '/mnt/c/').replace('\\', '/')
    elif isinstance(path, list):
        return [convert_path_to_wsl(p) for p in path]
    else:
        raise ValueError(f"Path must be a string or a list of strings, but got {type(path)}")


def run_vina(receptor_file, ligand_files, output_dir):
    """
    Run AutoDock Vina for multiple ligands.

    :param receptor_file: Path to the receptor PDBQT file.
    :param ligand_files: List of ligand PDBQT file paths.
    :param output_dir: Directory to store output files.
    :return: List of docking results with output file paths.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    receptor_file_wsl = convert_path_to_wsl(receptor_file)
    ligand_files_wsl = convert_path_to_wsl(ligand_files)
    output_results = []

    for ligand_file, ligand_file_wsl in zip(ligand_files, ligand_files_wsl):
        ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
        output_file = os.path.join(output_dir, f"{ligand_name}_docked.pdbqt")
        output_file_wsl = convert_path_to_wsl(output_file)

        vina_command = [
            'wsl', 'vina',
            '--receptor', receptor_file_wsl,
            '--ligand', ligand_file_wsl,
            '--out', output_file_wsl,
            '--center_x', '0', '--center_y', '0', '--center_z', '0',
            '--size_x', '20', '--size_y', '20', '--size_z', '20'
        ]

        logging.debug(f"Running Vina command: {' '.join(vina_command)}")

        result = subprocess.run(vina_command, capture_output=True, text=True)

        if result.returncode != 0:
            logging.error(f"AutoDock Vina failed for {ligand_file}: {result.stderr}")
            continue

        logging.debug(f"AutoDock Vina output for {ligand_file}: {result.stdout}")
        output_results.append({"ligand": ligand_file, "output": output_file})

    return output_results


def run_docking(receptor_file, ligand_files):
    """
    Run docking simulations for multiple ligands.
    
    :param receptor_file: Path to the receptor PDBQT file.
    :param ligand_files: List of ligand PDBQT file paths.
    :return: List of docking results.
    """
    docking_results = []
    
    for ligand_file in ligand_files:
        if not isinstance(ligand_file, str):
            logging.error(f"Invalid ligand file type: {type(ligand_file)}. Expected string.")
            continue
        
        ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
        output_file = os.path.join(os.path.dirname(ligand_file), f"{ligand_name}_docked.pdbqt")
        
        try:
            output = run_vina(receptor_file, ligand_file, output_file)
            docking_results.append({"ligand": ligand_file, "output": output})
        except RuntimeError as e:
            logging.error(f"Docking failed for {ligand_file}: {e}")
    
    return docking_results


# Example usage:
if __name__ == "__main__":
    RECEPTOR = "C:\\path\\to\\receptor.pdbqt"
    LIGANDS = [
        "C:\\path\\to\\ligand1.pdbqt",
        "C:\\path\\to\\ligand2.pdbqt"
    ]

    docking_results = run_docking(RECEPTOR, LIGANDS)
    logging.info(f"Docking completed. Results: {docking_results}")
