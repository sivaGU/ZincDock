# ZincDock GUI

ZincDock is a Streamlit-based interface that wraps AutoDock Vina and AutoDock4 (AD4) so metalloprotein docking runs are easy to configure, execute, and review.

Hosted app: https://metallodock-niv29aly8pujiceythjnx7.streamlit.app/

## Contents
- [Quick Start (Cloud)](#quick-start-cloud)
- [Run Locally](#run-locally)
- [Requirements](#requirements)
- [Project Structure](#project-structure)
- [Key Features](#key-features)
- [Documentation Overview](#documentation-overview)
- [Demo Tab Assets](#demo-tab-assets)
- [Tab-by-Tab Summary](#tab-by-tab-summary)
- [Usage Workflow](#usage-workflow)
- [Supported Formats](#supported-formats)
- [Configuration Files](#configuration-files)
- [Contact](#contact)

## Quick Start (Cloud)
1. Open the hosted Streamlit app.
2. Select the desired navigation tab (Demo, Standard AutoDock, or Metalloprotein Docking).
3. Follow the on-screen prompts to upload receptors, ligands, and configure grid parameters.
4. Use the *Documentation* tab for a detailed walkthrough of each step.

## Run Locally
1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd ZincDock-main
   ```
2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```
3. **Launch the app**
   ```bash
   streamlit run ZincDock.py
   ```
   *Windows:* double-click `run_docking.bat` to launch Streamlit.
4. Visit `http://localhost:8501` in your browser.

## Requirements
- Python 3.7+ (tested with 3.8–3.11)
- Streamlit, pandas (installed via `requirements.txt`)
- Windows-ready executables in `Files_for_GUI/`:
  - `vina.exe`
  - `autogrid4.exe`
  - `autodock4.exe`

*Linux/Mac:* compile or install equivalent binaries (without `.exe` extensions) and update paths in `ZincDock.py` or drop them inside `Files_for_GUI/`.

## Project Structure
```
ZincDock-main/
├── ZincDock.py          # Main Streamlit application
├── Files_for_GUI/          # Executables, parameters, sample ligands
│   ├── vina.exe
│   ├── autogrid4.exe
│   ├── autodock4.exe
│   ├── zinc_pseudo.py
│   ├── AD4_parameters.dat
│   ├── AD4Zn.dat
│   └── Ligands/
├── run_docking.bat         # Windows launcher
├── requirements.txt        # Python dependencies
└── README.md               # This document
```

## Key Features
- Unified Vina and AD4 workflows with zinc pseudoatom support
- Automatic executable discovery and parameter setup
- Grid box calculators, AD4 map generation, and progress tracking
- Downloads for per-ligand PDBQT/logs, CSV summaries, and ZIP archives
- Dedicated documentation tab mirroring the personalized vaccine pipeline

## Demo Tab Assets
The *Demo* tab is preconfigured for Carbonic Anhydrase I & II receptors. To use it:
1. **Download the demo assets** from the `Carbonic Anhydrase Receptor Files` and `18 PFAS Ligands` folders in the repository.
2. Place these folders alongside `ZincDock.py` (locally) or upload their contents to the Streamlit Cloud workspace under the same folder names.
3. In the app, choose either "Carbonic Anhydrase I" or "Carbonic Anhydrase II". Grid centers, sizes, spacing (0.375 Å), and docking parameters lock automatically.
4. Upload one of the bundled receptors (PDBQT) and select ligands from the PFAS set before running AD4 map building or docking.

Without these folders, the Demo tab cannot find receptors/ligands and will show missing-file warnings.

## Tab-by-Tab Summary
- **Home**: High-level overview and quick tips.
- **Documentation**: Full step-by-step instructions (see above) covering CSV ingestion, module selection, epitope features, VaxiJen workflow, aggregated file processing, and population coverage.
- **Demo**: AD4-only workflow with one-click presets for Carbonic Anhydrase I/II. Requires the demo receptor/ligand folders.
- **Standard AutoDock**: Vina-focused workflow with user-defined grid boxes.
- **Metalloprotein Docking**: AD4 workflow with full control over grid/map settings.

## Usage Workflow
1. **Upload receptor** (PDBQT) via file uploader or path.
2. **Prepare ligands** from a source folder or upload ready-made PDBQT files.
3. **Configure grid box** manually or via preset (Demo) / autodetect (Vina).
4. **Generate AD4 maps** when using the AD4 pipeline.
5. **Run docking**; monitor progress through the status table.
6. **Download results** (per-ligand logs/PDBQT, CSV summary, ZIP bundle).

## Supported Formats
- Receptors: PDBQT only
- Ligands: PDBQT, PDB, MOL2, SDF (converted to PDBQT during preparation)

## Configuration Files
Executables and parameters are auto-detected from `Files_for_GUI/`. Ensure permissions are set (`chmod +x` on Linux). Optional scripts:
- `zinc_pseudo.py` for inserting tetrahedral Zn pseudoatoms
- `AD4_parameters.dat` / `AD4Zn.dat` for AD4 scoring

## Contact
Questions, bug reports, or collaboration requests: **Dr. Sivanesan Dakshanamurthy** — sd233@georgetown.edu
---









