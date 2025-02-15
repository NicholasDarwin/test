import os
from flask import Flask, render_template_string, request
import subprocess
import pandas as pd
import re

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def run_gmx_pdb2gmx():
    # Ensure the 'data' folder exists
    if not os.path.exists("data"):
        os.makedirs("data")

    if request.method == 'POST':
        pdb_filename = request.form.get('pdb_filename')
        if pdb_filename:
            # Ensure that the filename doesn't have spaces or special characters
            pdb_filename = re.sub(r'[^a-zA-Z0-9_\-]', '', pdb_filename)

            # Modify the command to run the gmx command in WSL
            command = f"wsl gmx pdb2gmx -f data/{pdb_filename}.pdb -o data/{pdb_filename}_processed.gro -water spc"

            try:
                result = subprocess.run(command, input="6\n", capture_output=True, text=True, shell=True)

                # Capture the command output and errors
                output = result.stdout
                errors = result.stderr

                # Extract relevant information using regular expressions
                data = {
                    'Force Field': [],
                    'Chain': [],
                    'Residues': [],
                    'Atoms': [],
                    'Mass (a.m.u.)': [],
                    'Charge (e)': []
                }

                # Force Field
                force_field_match = re.search(r"Using the (.*?) force field", output)
                if force_field_match:
                    force_field = force_field_match.group(1)
                    data['Force Field'].append(force_field)

                # Chain and residue information
                chain_matches = re.findall(r"chain  #res #atoms\n\s+(\d) '(.*?)'\s+(\d+)\s+(\d+)", output)
                for match in chain_matches:
                    data['Chain'].append(match[1])
                    data['Residues'].append(match[2])
                    data['Atoms'].append(match[3])

                # Mass and Charge (from the "Total mass" and "Total charge" lines)
                mass_match = re.search(r"Total mass (\d+\.\d+) a.m.u.", output)
                charge_match = re.search(r"Total charge (\d+\.\d+) e", output)

                if mass_match and charge_match:
                    total_mass = mass_match.group(1)
                    total_charge = charge_match.group(1)
                    data['Mass (a.m.u.)'].append(total_mass)
                    data['Charge (e)'].append(total_charge)

                # Create a DataFrame from the collected data
                df = pd.DataFrame(data)

                # Convert the dataframe to HTML for rendering
                output_html = df.to_html(classes='table table-striped')

                # Return HTML table with output and errors
                return render_template_string(""" 
                    <h1>GROMACS Command Output</h1>
                    <h2>Output:</h2>
                    {{ output_html|safe }}
                    <h2>Errors:</h2>
                    <pre>{{ errors }}</pre>
                    <hr>
                    <a href="/">Back</a>
                """, output_html=output_html, errors=errors)
            except Exception as e:
                return f"Error running command: {e}"
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
