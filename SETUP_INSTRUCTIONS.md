# MetalBind Setup Instructions

## Overview
MetalBind is a Streamlit-based GUI for AutoDock Vina and AutoDock4 (AD4) metalloprotein docking. The Demo tab is configured for 8 Zinc Metal Proteins.

## Quick Start

### 1. Install Dependencies
```bash
pip install -r requirements.txt
```

### 2. Setup Demo Files (Optional)
For the Demo tab to work, you need to populate two folders:

- **`(DEMO) Zinc Metal Protein Receptors/`**: Place 8 receptor PDBQT files:
  - hACE.pdbqt
  - HDAC2.pdbqt
  - HDAC8_with_Hydroxamic_Acid.pdbqt
  - HDAC8_with_SAHA.pdbqt
  - HDAC10.pdbqt
  - Human_Neutral_Endopeptidase.pdbqt
  - Leukotriene.pdbqt
  - ADAMTS-5.pdbqt

- **`(DEMO) Ligands/`**: Place 8 endogenous ligand PDBQT files corresponding to each receptor.

### 3. Run the Application

**Windows:**
Double-click `run_docking.bat` or run:
```bash
streamlit run MetalBind.py
```

**Linux/Mac:**
```bash
streamlit run MetalBind.py
```

The application will open in your browser at `http://localhost:8501`

## Project Structure

```
FINAL METALBIND/
├── MetalBind.py                    # Main Streamlit application
├── README.md                       # Full documentation
├── requirements.txt                # Python dependencies
├── run_docking.bat                 # Windows launcher script
├── setup_executables.py            # Executable setup script
├── .gitignore                      # Git ignore file
├── Files_for_GUI/                  # Executables and parameters
│   ├── vina                        # AutoDock Vina executable
│   ├── autogrid4                   # AutoGrid4 executable
│   ├── autodock4                   # AutoDock4 executable
│   ├── AD4_parameters.dat          # AD4 parameters
│   ├── AD4Zn.dat                   # AD4 zinc parameters
│   ├── zinc_pseudo.py              # Zinc pseudoatom script
│   └── Ligands/                    # Sample ligands
├── (DEMO) Zinc Metal Protein Receptors/  # Demo receptors (populate with 8 PDBQT files)
└── (DEMO) Ligands/                 # Demo ligands (populate with 8 PDBQT files)
```

## Features

- **Demo Tab**: Pre-configured for 8 Zinc Metal Proteins with locked grid box parameters
- **Standard AutoDock Tab**: Vina-focused workflow with user-defined grid boxes
- **Metalloprotein Docking Tab**: AD4 workflow with full control over grid/map settings

## Notes

- The Demo tab requires the demo receptor and ligand folders to be populated
- Grid box coordinates are automatically set when selecting a zinc metal protein in the Demo tab
- All executables are included in `Files_for_GUI/` for Windows (Linux/Mac users need to compile or install equivalent binaries)

## Contact

Questions or issues: See README.md for contact information.
