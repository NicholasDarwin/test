import os
import subprocess

def prepare_system(receptor_pdb, ligand_pdb, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Solvate the protein-ligand complex
    solvate_command = f"gmx solvate -cp {receptor_pdb} -cs spc216.gro -o {output_dir}/solvated.gro -p {output_dir}/topol.top"
    subprocess.run(solvate_command, shell=True, check=True)

    # Add ions to the system
    ion_command = f"gmx grompp -f ions.mdp -c {output_dir}/solvated.gro -p {output_dir}/topol.top -o {output_dir}/ions.tpr"
    subprocess.run(ion_command, shell=True, check=True)

    # Generate the ion configuration
    subprocess.run("echo SOL | gmx genion -s {output_dir}/ions.tpr -o {output_dir}/ionized.gro -p {output_dir}/topol.top -neutral", shell=True, check=True)

    return os.path.join(output_dir, "ionized.gro")