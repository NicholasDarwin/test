import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def plot_binding_affinities(binding_affinities, output_file='binding_affinities.png'):
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Ligand', y='Binding Affinity (kcal/mol)', data=binding_affinities)
    plt.title('Binding Affinities of Ligands')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def plot_rmsd(rmsd_data, output_file='rmsd_plot.png'):
    plt.figure(figsize=(10, 6))
    for ligand, rmsd_values in rmsd_data.items():
        plt.plot(rmsd_values, label=ligand)
    plt.title('RMSD of Protein-Ligand Complexes')
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.legend()
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