import os
from flask import Flask, request, jsonify
from docking.run_autodock_vina import run_docking
from dynamics.run_gromacs import run_md_simulation
from analysis.binding_free_energy import calculate_binding_free_energy
from analysis.visualize_results import visualize_results
from analysis.analyze_trajectories import analyze_trajectory

app = Flask(__name__)

@app.route('/run_pipeline', methods=['POST'])
def run_pipeline():
    data = request.json
    receptor_file = data.get('receptor_file')
    ligand_files = data.get('ligand_files')

    if not receptor_file or not ligand_files:
        return jsonify({"error": "Receptor file and ligand files are required."}), 400

    # Step 1: Run molecular docking
    docking_results = run_docking(receptor_file, ligand_files)

    # Step 2: Run molecular dynamics simulations for top ligands
    md_results = run_md_simulation(docking_results['top_ligands'])

    # Step 3: Calculate binding free energy
    binding_energy = calculate_binding_free_energy(md_results)

    # Step 4: Visualize results
    visualize_results(binding_energy)

    # Step 5: Analyze trajectories
    trajectory_analysis = analyze_trajectory(md_results)

    return jsonify({
        "docking_results": docking_results,
        "md_results": md_results,
        "binding_energy": binding_energy,
        "trajectory_analysis": trajectory_analysis
    })

if __name__ == "__main__":
    app.run(debug=True)