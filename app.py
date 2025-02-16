import os
import subprocess
import time
import pandas as pd
import logging
from flask import Flask, request, jsonify, render_template
from src.docking.prepare_receptors import prepare_receptor
from src.docking.prepare_ligands import prepare_ligand
from src.docking.run_autodock_vina import run_vina
from src.dynamics.run_gromacs import run_gromacs_simulation
from src.analysis.binding_free_energy import calculate_binding_free_energy
from src.analysis.visualize_results import visualize_results
from src.analysis.analyze_trajectories import analyze_trajectory

app = Flask(__name__)

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# Set up logging
logging.basicConfig(level=logging.DEBUG)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_pipeline', methods=['POST'])
def run_pipeline():
    try:
        receptor_file = request.form.get('receptor_file')
        ligand_files = request.form.get('ligand_files')

        if not receptor_file or not ligand_files:
            return jsonify({"error": "Receptor file and ligand files are required."}), 400

        # Append .pdb extension to the file names if not already present
        if not receptor_file.endswith('.pdb'):
            receptor_file += '.pdb'
        ligand_files = [file.strip() + '.pdb' if not file.strip().endswith('.pdb') else file.strip() for file in ligand_files.split(',')]

        receptor_file_path = os.path.join(DATA_DIR, receptor_file)
        ligand_files = [os.path.join(DATA_DIR, file) for file in ligand_files]

        logging.debug(f"Receptor file path: {receptor_file_path}")
        logging.debug(f"Ligand files: {ligand_files}")

        # Prepare receptor and ligands
        logging.debug("Preparing receptor...")
        try:
            prepared_receptor_file = prepare_receptor(receptor_file_path)
            logging.debug(f"Prepared receptor file: {prepared_receptor_file}")
        except Exception as e:
            logging.error(f"Error preparing receptor: {e}")
            return jsonify({"error": f"Error preparing receptor: {e}"}), 500

        logging.debug("Preparing ligands...")
        try:
            prepared_ligand_files = []
            for file in ligand_files:
                prepared_ligand_file = prepare_ligand(file)
                prepared_ligand_files.append(prepared_ligand_file)
                # Add a short delay to ensure the file system is updated
                time.sleep(10)
            logging.debug(f"Prepared ligand files: {prepared_ligand_files}")
        except Exception as e:
            logging.error(f"Error preparing ligands: {e}")
            return jsonify({"error": f"Error preparing ligands: {e}"}), 500

        # Run molecular docking
        logging.debug("Running molecular docking...")
        docking_results = run_vina(prepared_receptor_file, prepared_ligand_files, DATA_DIR)
        logging.debug(f"Docking results: {docking_results}")

        # Extract top ligands from docking results
        top_ligands = [result['pdb_output'] for result in docking_results]
        logging.debug(f"Top ligands: {top_ligands}")

        # Step 2: Run molecular dynamics simulations for top ligands
        logging.debug("Running molecular dynamics simulations...")
        md_results = run_gromacs_simulation(prepared_receptor_file, top_ligands, DATA_DIR)
        logging.debug(f"Molecular dynamics results: {md_results}")

        # Step 3: Calculate binding free energy
        logging.debug("Calculating binding free energy...")
        binding_energy_results = []
        for ligand_file in top_ligands:
            binding_energy, binding_data = calculate_binding_free_energy(ligand_file, DATA_DIR)
            binding_energy_results.append({
                "ligand_file": ligand_file,
                "binding_energy": binding_energy,
                "binding_data": binding_data
            })
        logging.debug(f"Binding free energy results: {binding_energy_results}")

        # Step 4: Visualize results
        logging.debug("Visualizing results...")
        visualizations = visualize_results(binding_energy_results)
        logging.debug(f"Visualizations: {visualizations}")

        # Step 5: Analyze trajectories
        logging.debug("Analyzing trajectories...")
        trajectory_analysis = analyze_trajectory(md_results)
        logging.debug(f"Trajectory analysis: {trajectory_analysis}")

        return jsonify({
            "docking_results": docking_results,
            "md_results": md_results,
            "binding_energy_results": binding_energy_results,
            "trajectory_analysis": trajectory_analysis,
            "visualizations": visualizations
        })
    except Exception as e:
        logging.error(f"Error in run_pipeline: {e}")
        return jsonify({"error": f"Error in run_pipeline: {e}"}), 500

if __name__ == "__main__":
    app.run(debug=True)
