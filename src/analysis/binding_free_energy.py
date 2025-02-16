import subprocess
import os

def convert_path_to_wsl(path):
    """
    Convert a Windows path to a WSL path.
    
    :param path: Windows path.
    :return: WSL path.
    """
    return path.replace('C:\\', '/mnt/c/').replace('\\', '/')

def calculate_binding_free_energy(docked_pdbqt, data_dir):
    """
    Calculate the binding free energy using AutoDock Vina.
    
    :param docked_pdbqt: Path to the docked PDBQT file.
    :param data_dir: Directory to store intermediate and output files.
    :return: Binding energy and binding data.
    """
    # Convert paths to WSL format
    docked_pdbqt_wsl = convert_path_to_wsl(docked_pdbqt)
    data_dir_wsl = convert_path_to_wsl(data_dir)
    pdb_file_wsl = os.path.join(data_dir_wsl, 'docked_complex.pdb').replace('\\', '/')
    receptor_pdbqt_wsl = os.path.join(data_dir_wsl, 'receptor.pdbqt').replace('\\', '/')
    output_pdbqt_wsl = os.path.join(data_dir_wsl, 'vina_output.pdbqt').replace('\\', '/')
    log_file_wsl = os.path.join(data_dir_wsl, 'vina_log.txt').replace('\\', '/')

    # Step 1: Convert the PDBQT file to PDB format using Open Babel
    subprocess.run(['wsl', 'obabel', docked_pdbqt_wsl, '-O', pdb_file_wsl], check=True)

    # Step 2: Prepare the receptor and ligand files for AutoDock Vina
    ligand_pdbqt_wsl = docked_pdbqt_wsl  # The docked complex file will be treated as ligand

    # Step 3: Run AutoDock Vina to calculate binding free energy
    vina_command = [
        'wsl', 'vina',
        '--receptor', receptor_pdbqt_wsl,
        '--ligand', ligand_pdbqt_wsl,
        '--out', output_pdbqt_wsl,
        '--log', log_file_wsl,
        '--center_x', '0', '--center_y', '0', '--center_z', '0',
        '--size_x', '20', '--size_y', '20', '--size_z', '20'
    ]
    subprocess.run(vina_command, check=True)

    # Step 4: Extract binding energy from the Vina log file
    binding_energy = None
    binding_data = []
    with open(log_file_wsl, 'r') as f:
        for line in f:
            binding_data.append(line.strip())
            if "REMARK VINA RESULT:" in line:
                binding_energy = line.split()[3]  # Extract the binding energy value

    return binding_energy, binding_data