import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial import distance_matrix

def load_pdb(file_path):
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                atom = {
                    'name': parts[2],
                    'resname': parts[3],
                    'chain': parts[4],
                    'resid': int(parts[5]),
                    'x': float(parts[6]),
                    'y': float(parts[7]),
                    'z': float(parts[8])
                }
                atoms.append(atom)
    return atoms

def calculate_rmsd(trajectory, reference):
    assert len(trajectory) == len(reference), "Trajectory and reference must have the same number of atoms"
    diff = np.array([[atom['x'], atom['y'], atom['z']] for atom in trajectory]) - \
           np.array([[atom['x'], atom['y'], atom['z']] for atom in reference])
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

def analyze_rmsd(trajectory_files, reference_file):
    reference = load_pdb(reference_file)
    rmsd_values = []
    for traj_file in trajectory_files:
        trajectory = load_pdb(traj_file)
        rmsd = calculate_rmsd(trajectory, reference)
        rmsd_values.append(rmsd)
    return rmsd_values

def plot_rmsd(rmsd_values, output_file):
    plt.figure()
    plt.plot(rmsd_values, label='RMSD')
    plt.xlabel('Frame')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD over Time')
    plt.legend()
    plt.savefig(output_file)
    plt.close()

def analyze_binding_interactions(trajectory_files, ligand_selection):
    interaction_distances = []
    for traj_file in trajectory_files:
        trajectory = load_pdb(traj_file)
        ligand_atoms = [atom for atom in trajectory if atom['resname'] == ligand_selection]
        receptor_atoms = [atom for atom in trajectory if atom['resname'] != ligand_selection]
        ligand_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in ligand_atoms])
        receptor_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in receptor_atoms])
        distances = distance_matrix(ligand_coords, receptor_coords)
        interaction_distances.append(distances)
    return interaction_distances

def plot_binding_affinities(binding_affinities, output_file='binding_affinities.png'):
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Ligand', y='Binding Affinity (kcal/mol)', data=binding_affinities)
    plt.title('Binding Affinities of Ligands')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_heatmap(data, output_file='heatmap.png'):
    plt.figure(figsize=(10, 8))
    sns.heatmap(data, annot=True, fmt=".2f", cmap='coolwarm')
    plt.title('Heatmap of Interaction Energies')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def visualize_results(binding_free_energy_results):
    # Example visualization: Heatmap of binding free energy results
    plt.figure(figsize=(10, 8))
    sns.heatmap(binding_free_energy_results.pivot("Ligand", "Receptor", "FreeEnergy"), annot=True, cmap="coolwarm")
    plt.title("Binding Free Energy Heatmap")
    plt.savefig("data/binding_free_energy_heatmap.png")
    plt.close()

    # Example visualization: RMSD plot
    plt.figure(figsize=(10, 8))
    sns.lineplot(data=binding_free_energy_results, x="Time", y="RMSD", hue="Ligand")
    plt.title("RMSD Over Time")
    plt.savefig("data/rmsd_plot.png")
    plt.close()

    # Example visualization: Interaction map
    plt.figure(figsize=(10, 8))
    sns.heatmap(binding_free_energy_results.pivot("Ligand", "Interaction", "Frequency"), annot=True, cmap="viridis")
    plt.title("Interaction Map")
    plt.savefig("data/interaction_map.png")
    plt.close()

    # Return paths to the generated visualizations
    visualizations = {
        "heatmap": "data/binding_free_energy_heatmap.png",
        "rmsd_plot": "data/rmsd_plot.png",
        "interaction_map": "data/interaction_map.png"
    }

    return visualizations

def analyze_trajectory(trajectory_files, reference_file, ligand_selection, output_rmsd_file, output_interaction_file):
    rmsd_values = analyze_rmsd(trajectory_files, reference_file)
    plot_rmsd(rmsd_values, output_rmsd_file)

    interaction_distances = analyze_binding_interactions(trajectory_files, ligand_selection)
    np.save(output_interaction_file, interaction_distances)

if __name__ == "__main__":
    # Example usage
    analyze_trajectory(['path/to/trajectory1.pdb', 'path/to/trajectory2.pdb'], 'path/to/reference.pdb', 'LIG', 'rmsd_plot.png', 'binding_interactions.npy')