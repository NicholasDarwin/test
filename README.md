# Virtual Screening Pipeline

This project implements a virtual screening pipeline for molecular docking and dynamics simulations. It allows users to input a protein receptor and a set of ligand molecules, performing simulations to evaluate binding affinities and molecular dynamics.

## Project Structure

```
virtual_screening_pipeline
├── data
│   ├── ligands          # Directory for ligand PDB files
│   └── receptors        # Directory for protein receptor PDB files
├── src
│   ├── docking
│   │   ├── run_autodock_vina.py  # Main logic for AutoDock Vina docking
│   │   └── prepare_ligands.py     # Functions to prepare ligands for docking
│   ├── dynamics
│   │   ├── run_gromacs.py         # Logic to run GROMACS simulations
│   │   └── prepare_system.py       # Functions to prepare the simulation system
│   ├── analysis
│   │   ├── binding_free_energy.py  # Functions for binding free energy calculations
│   │   ├── visualize_results.py     # Functions for result visualizations
│   │   └── analyze_trajectories.py  # Functions to analyze MD trajectories
│   ├── utils
│   │   └── file_handling.py        # Utility functions for file handling
│   └── main.py                     # Entry point for the pipeline
├── requirements.txt                # Required Python packages and dependencies
└── README.md                       # Documentation for the project
```

## Installation

1. Clone the repository:
   ```
   git clone <repository-url>
   cd virtual_screening_pipeline
   ```

2. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

## Usage

1. Place your protein receptor PDB files in the `data/receptors` directory.
2. Place your ligand PDB files in the `data/ligands` directory.
3. Run the main pipeline:
   ```
   python src/main.py
   ```

## Workflow Overview

1. **Preparation**: Ligands are prepared for docking using `prepare_ligands.py`.
2. **Molecular Docking**: AutoDock Vina is executed via `run_autodock_vina.py` to perform docking simulations and retrieve binding affinities.
3. **Molecular Dynamics**: The top-ranked ligands undergo molecular dynamics simulations using GROMACS, managed by `run_gromacs.py`.
4. **Analysis**: Binding free energy calculations and trajectory analyses are performed using the analysis scripts.
5. **Visualization**: Results are visualized using tools like PyMOL or NGLView.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the LICENSE file for details.