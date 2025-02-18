import os
import time
import logging
import pandas as pd
import subprocess
import numpy as np
import scipy.integrate as spi
import sympy as sp
from flask import Flask, request, jsonify, render_template, send_file
from flask_socketio import SocketIO, emit
from rdkit import Chem
from rdkit.Chem import Descriptors
from src.docking.prepare_ligands import prepare_ligand
from src.docking.prepare_receptors import prepare_receptor
from src.dynamics.run_gromacs import run_gromacs_simulation
from src.docking.run_autodock_vina import run_vina

app = Flask(__name__)
socketio = SocketIO(app)

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
PYMOL_PATH = r"C:\Users\Lenovo\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\PyMOL (Anaconda3 (64-bit))\PyMOL + Console.lnk"  # Update this path to the actual PyMOL executable

# Set up logging
logging.basicConfig(level=logging.DEBUG)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/get_pdb/<filename>')
def get_pdb(filename):
    pdb_path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(pdb_path):
        return jsonify({"error": "File not found."}), 404
    return send_file(pdb_path, mimetype='chemical/x-pdb')

@app.route('/get_image/<filename>')
def get_image(filename):
    image_path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(image_path):
        return jsonify({"error": "File not found."}), 404
    return send_file(image_path, mimetype='image/png')

@app.route('/simulation_results')
def simulation_results():
    df = pd.read_csv(os.path.join(DATA_DIR, 'simulation_results.csv'))  # Load results
    return render_template('results.html', tables=[df.to_html(classes='data')])

def send_terminal_output(message):
    socketio.emit('terminal_output', message)

@app.route('/prepare_receptor', methods=['POST'])
def prepare_receptor_route():
    data = request.json
    receptor_file = data.get('receptor_file')
    if not receptor_file:
        return jsonify({"error": "Receptor file is required."}), 400

    receptor_file_path = os.path.join(DATA_DIR, receptor_file)
    logging.debug(f"Receptor file path: {receptor_file_path}")

    send_terminal_output("Preparing receptor...")
    prepared_receptor_file = prepare_receptor(receptor_file_path)
    if not prepared_receptor_file:
        return jsonify({"error": "Failed to prepare receptor file."}), 500
    send_terminal_output(f"Prepared receptor file: {prepared_receptor_file}")
    return jsonify({"message": "Receptor prepared successfully."})

@app.route('/prepare_ligands', methods=['POST'])
def prepare_ligands_route():
    data = request.json
    ligand_files = data.get('ligand_files')
    if not ligand_files:
        return jsonify({"error": "Ligand files are required."}), 400

    ligand_files = [file.strip() for file in ligand_files.split(',')]
    ligand_files = [os.path.join(DATA_DIR, file) for file in ligand_files]
    logging.debug(f"Ligand files: {ligand_files}")

    send_terminal_output("Preparing ligands...")
    prepared_ligand_files = []
    for file in ligand_files:
        prepared_ligand_file = prepare_ligand(file)
        if not prepared_ligand_file:
            return jsonify({"error": f"Failed to prepare ligand file: {file}"}), 500
        prepared_ligand_files.append(prepared_ligand_file)
        time.sleep(10)  # Ensure filesystem updates
    send_terminal_output(f"Prepared ligand files: {prepared_ligand_files}")
    return jsonify({"message": "Ligands prepared successfully."})

@app.route('/run_docking', methods=['POST'])
def run_docking_route():
    receptor_file = request.form.get('receptor_file')
    ligand_files = request.form.getlist('ligand_files')
    
    if not receptor_file:
        return jsonify({"error": "Receptor file is required."}), 400
    if not ligand_files:
        return jsonify({"error": "Ligand files are required."}), 400
    
    logging.debug(f"Receptor file path: {receptor_file}")
    logging.debug(f"Ligand files: {ligand_files}")
    
    prepared_receptor_file = prepare_receptor(receptor_file)
    logging.debug(f"Prepared receptor file: {prepared_receptor_file}")
    
    for file in ligand_files:
        prepared_ligand_file = prepare_ligand(file)
        logging.debug(f"Prepared ligand file: {prepared_ligand_file}")
        # ...existing code to run docking...

@app.route('/run_md', methods=['POST'])
def run_md_route():
    data = request.json
    receptor_file = data.get('receptor_file')
    ligand_files = data.get('ligand_files')
    if not receptor_file or not ligand_files:
        return jsonify({"error": "Receptor file and ligand files are required."}), 400

    receptor_file_path = os.path.join(DATA_DIR, receptor_file)
    ligand_files = [file.strip() for file in ligand_files.split(',')]
    ligand_files = [os.path.join(DATA_DIR, file) for file in ligand_files]
    logging.debug(f"Receptor file path: {receptor_file_path}")
    logging.debug(f"Ligand files: {ligand_files}")

    send_terminal_output("Running molecular dynamics simulations...")
    md_results = run_gromacs_simulation(receptor_file_path, ligand_files[0])
    if not md_results:
        return jsonify({"error": "Failed to run MD simulation."}), 500

    send_terminal_output(md_results)
    return jsonify({"message": "Molecular dynamics simulations completed successfully."})

@app.route('/open_pymol', methods=['POST'])
def open_pymol():
    try:
        receptor_file = request.json.get('receptor_file')
        ligand_files = request.json.get('ligand_files').split(',')
        receptor_file_path = os.path.join(DATA_DIR, receptor_file)
        ligand_file_paths = [os.path.join(DATA_DIR, file.strip()) for file in ligand_files]

        # Generate PyMOL script to load and visualize the molecules
        pymol_script_path = os.path.join(DATA_DIR, 'visualize.pml')
        with open(pymol_script_path, 'w') as script_file:
            script_file.write(f'load {receptor_file_path}, receptor\n')
            for i, ligand_file in enumerate(ligand_file_paths):
                script_file.write(f'load {ligand_file}, ligand{i+1}\n')
            script_file.write('show sticks, receptor\n')
            for i in range(len(ligand_file_paths)):
                script_file.write(f'show sticks, ligand{i+1}\n')
            script_file.write('zoom\n')

        # Open PyMOL with the generated script
        subprocess.run([PYMOL_PATH, '-r', pymol_script_path])
        return jsonify({"message": "PyMOL opened successfully."})
    except Exception as e:
        logging.error(f"Error opening PyMOL: {e}")
        return jsonify({"error": f"Failed to open PyMOL: {e}"}), 500

@app.route('/advanced_calculations')
def advanced_calculations():
    # Example advanced math calculation: Solving a differential equation
    t = sp.symbols('t')
    y = sp.Function('y')
    ode = sp.Eq(y(t).diff(t, t) - 2*y(t).diff(t) + y(t), sp.sin(t))
    solution = sp.dsolve(ode)

    # Example advanced chemistry calculation: Molecular weight of a molecule
    mol = Chem.MolFromSmiles('CCO')
    mol_weight = Descriptors.MolWt(mol)

    return render_template('advanced_calculations.html', solution=solution, mol_weight=mol_weight)

@socketio.on('connect')
def handle_connect():
    send_terminal_output("Connected to server.")

if __name__ == "__main__":
    socketio.run(app, debug=True)
