# MetalloDock Streamlit App

MetalloDock is a Streamlit interface that automates the AutoDock4Zn workflow on Linux systems. It wraps the official AutoDock4Zn scripts and binaries to:

- Add tetrahedral zinc pseudo atoms to receptors (`zinc_pseudo.py`).
- Generate AD4 grid parameter files and affinity maps (`autogrid4`).
- Run AutoDock Vina with the AD4 scoring function (`vina --scoring ad4`).

## Repository Contents

- `metallodock_app.py` – main Streamlit application with Home, Demo, MetalloDock, and Documentation pages.
- `Files_for_GUI/` – bundled AutoGrid4/Vina executables, AD4Zn parameters, ligands, and scripts required by the workflow.
- `Carbonic_Anhydrase_I.pdbqt` – demo receptor used for validation.
- `PFOA_Test_Ligand.pdbqt` – demo ligand paired with the receptor above.
- `Demo_Gridbox.docx` – grid box dimensions and coordinate reference for the demo.
- `OFFICIAL - 2025 Docking Protocol.docx` – source protocol that inspired the documentation page.

## Setup (Linux)

1. **Install dependencies**
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

2. **Make bundled binaries executable**
   ```bash
   chmod +x Files_for_GUI/autogrid4
   chmod +x Files_for_GUI/autodock4
   chmod +x Files_for_GUI/vina
   ```

3. **Launch Streamlit**
   ```bash
   streamlit run metallodock_app.py
   ```

## Usage Highlights

- **Demo Page:** Reproduces the Carbonic Anhydrase I vs. PFOA docking with preset inputs and grid box (from `Demo_Gridbox.docx`).
- **MetalloDock Page:** Upload or reference your own receptor/ligand PDBQT files, adjust grid center & dimensions, and tune Vina parameters such as exhaustiveness, CPU count, and number of modes.
- **Outputs:** Each run is stored under `metallodock_runs/`, containing the TZ-adjusted receptor, `.gpf` file, AutoGrid log, Vina poses, and Vina log.

## Notes

- The packaged binaries are the Windows versions in the source project. Replace them with Linux builds of AutoGrid4, AutoDock4, and AutoDock Vina when deploying on Linux.
- Ensure Python can execute `Files_for_GUI/zinc_pseudo.py` (ships from the AutoDock-Vina examples).
- For batch docking, keep ligands in a dedicated folder and invoke the workflow repeatedly with different inputs and map prefixes.

