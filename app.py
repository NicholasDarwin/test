import os
from flask import Flask, render_template_string, request
import subprocess
import pandas as pd
import re
from virtual_screening_pipeline.src.docking.run_autodock_vina import run_vina
from virtual_screening_pipeline.src.dynamics.run_gromacs import run_gromacs
from virtual_screening_pipeline.src.analysis.binding_free_energy import calculate_binding_free_energy
from virtual_screening_pipeline.src.analysis.visualize_results import generate_visualizations

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def run_pipeline():
    # Ensure the 'data' folder exists
    if not os.path.exists("data"):
        os.makedirs("data")

    if request.method == 'POST':
        pdb_filename = request.form.get('pdb_filename')
        if pdb_filename:
            # Ensure that the filename doesn't have spaces or special characters
            pdb_filename = re.sub(r'[^a-zA-Z0-9_\-]', '', pdb_filename)

            # Run AutoDock Vina
            docking_results = run_vina(pdb_filename)

            # Select top-ranked ligands and run GROMACS MD simulations
            md_results = run_gromacs(docking_results)

            # Perform binding free energy calculations
            binding_free_energy_results = calculate_binding_free_energy(md_results)

            # Generate visualizations
            visualizations = generate_visualizations(binding_free_energy_results)

            # Convert the results to HTML for rendering
            docking_html = docking_results.to_html(classes='table table-striped')
            md_html = md_results.to_html(classes='table table-striped')
            binding_free_energy_html = binding_free_energy_results.to_html(classes='table table-striped')
            visualizations_html = visualizations.to_html(classes='table table-striped')

            # Return HTML table with results
            return render_template_string(""" 
                <h1>Virtual Screening Pipeline Results</h1>
                <h2>Docking Results:</h2>
                {{ docking_html|safe }}
                <h2>Molecular Dynamics Results:</h2>
                {{ md_html|safe }}
                <h2>Binding Free Energy Results:</h2>
                {{ binding_free_energy_html|safe }}
                <h2>Visualizations:</h2>
                {{ visualizations_html|safe }}
                <hr>
                <a href="/">Back</a>
            """, docking_html=docking_html, md_html=md_html, binding_free_energy_html=binding_free_energy_html, visualizations_html=visualizations_html)
        else:
            return "Please provide a PDB file name."
    
    return render_template_string(""" 
        <h1>Enter PDB Filename</h1>
        <form method="post">
            <label for="pdb_filename">PDB Filename (without extension):</label>
            <input type="text" id="pdb_filename" name="pdb_filename" required>
            <button type="submit">Submit</button>
        </form>
    """)

if __name__ == "__main__":
    app.run(debug=True)