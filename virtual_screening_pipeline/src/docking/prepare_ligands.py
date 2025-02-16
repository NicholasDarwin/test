import os
from rdkit import Chem

def prepare_ligand(ligand_file):
    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
    
    # Load the ligand using RDKit
    ligand_mol = Chem.MolFromPDBFile(ligand_file, removeHs=True)
    if ligand_mol is None:
        raise ValueError(f"Could not read ligand file: {ligand_file}")

    # Add hydrogens
    ligand_mol = Chem.AddHs(ligand_mol)

    # Generate 3D coordinates
    from rdkit.Chem import AllChem
    AllChem.EmbedMolecule(ligand_mol)
    AllChem.UFFOptimizeMolecule(ligand_mol)

    # Save the prepared ligand to a new PDB file
    prepared_ligand_file = os.path.join("data", "ligands", f"{ligand_name}_prepared.pdb")
    with open(prepared_ligand_file, 'w') as f:
        f.write(Chem.MolToPDBBlock(ligand_mol))

    return prepared_ligand_file

def prepare_ligands(ligand_directory):
    prepared_ligands = []
    for ligand_file in os.listdir(ligand_directory):
        if ligand_file.endswith('.pdb'):
            full_path = os.path.join(ligand_directory, ligand_file)
            prepared_ligand = prepare_ligand(full_path)
            prepared_ligands.append(prepared_ligand)
    return prepared_ligands