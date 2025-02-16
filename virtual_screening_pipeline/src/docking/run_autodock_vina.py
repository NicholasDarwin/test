import subprocess
import pandas as pd

def run_vina(pdb_filename):
    # Define the command to run AutoDock Vina
    command = f"wsl vina --receptor data/{pdb_filename}_receptor.pdbqt --ligand data/{pdb_filename}_ligand.pdbqt --out data/{pdb_filename}_out.pdbqt --log data/{pdb_filename}_log.txt"
    
    # Run the command
    subprocess.run(command, shell=True, check=True)
    
    # Parse the results (assuming the log file contains the binding affinities)
    with open(f"data/{pdb_filename}_log.txt", "r") as log_file:
        lines = log_file.readlines()
    
    # Extract binding affinities (this is a simplified example, adjust as needed)
    affinities = []
    for line in lines:
        if "REMARK VINA RESULT" in line:
            parts = line.split()
            affinities.append(float(parts[3]))
    
    # Create a DataFrame with the results
    results_df = pd.DataFrame({"Ligand": [pdb_filename], "Affinity": affinities})
    
    return results_df

if __name__ == "__main__":
    receptor_path = "data/receptors/example_receptor.pdb"
    ligand_path = "data/ligands/example_ligand.pdb"
    output_path = "data/ligands/docking_results.pdb"

    try:
        output = run_vina("example")
        print("Docking completed successfully.")
        print(output)
    except Exception as e:
        print(e)