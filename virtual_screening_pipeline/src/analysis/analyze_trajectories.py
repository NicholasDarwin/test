import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def analyze_rmsd(trajectory_file, reference_structure):
    trajectory = md.load(trajectory_file)
    reference = md.load(reference_structure)

    # Align the trajectory to the reference structure
    aligned_trajectory = md.align(trajectory, reference)

    # Calculate RMSD
    rmsd = md.rmsd(aligned_trajectory, reference)

    return rmsd

def plot_rmsd(rmsd_values, output_file):
    plt.figure()
    plt.plot(rmsd_values, label='RMSD')
    plt.xlabel('Frame')
    plt.ylabel('RMSD (nm)')
    plt.title('RMSD over Time')
    plt.legend()
    plt.savefig(output_file)
    plt.close()

def analyze_binding_interactions(trajectory_file, ligand_selection):
    trajectory = md.load(trajectory_file)
    ligand_indices = trajectory.topology.select(ligand_selection)

    # Analyze distances or interactions
    distances = []
    for frame in trajectory:
        for index in ligand_indices:
            # Example: Calculate distance to a specific atom in the receptor
            receptor_atom = trajectory.topology.select('protein and name CA')[0]
            distance = md.compute_distances(frame, [[index, receptor_atom]])
            distances.append(distance)

    return np.array(distances)

def main(trajectory_file, reference_structure, ligand_selection, output_rmsd_file, output_interaction_file):
    rmsd_values = analyze_rmsd(trajectory_file, reference_structure)
    plot_rmsd(rmsd_values, output_rmsd_file)

    interaction_distances = analyze_binding_interactions(trajectory_file, ligand_selection)
    np.savetxt(output_interaction_file, interaction_distances)

if __name__ == "__main__":
    # Example usage
    main('path/to/trajectory.xtc', 'path/to/reference.pdb', 'resname LIG', 'rmsd_plot.png', 'binding_interactions.txt')