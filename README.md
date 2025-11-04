# MetalloDock GUI

A user-friendly Streamlit-based GUI for molecular docking of ligands to metalloprotein receptors using AutoDock Vina and AutoDock4 (AD4).

Link: https://metallodock-9rfsp7szmdt7fevwbdb9vh.streamlit.app/

### Additional Option 1: Run Locally

1. **Clone this repository:**
   ```bash
   git clone <repository-url>
   cd "MetalloDock GUI GitHub"
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Launch the GUI:**
   ```bash
   streamlit run MetalloDock.py
   ```
   Or double-click `run_docking.bat` on Windows.

4. **Open your browser:**
   The GUI will automatically open at `http://localhost:8501`

## 📋 Requirements

- **Python 3.7+** (tested with Python 3.8, 3.9, 3.10, 3.11)
- **Streamlit** (installed via `requirements.txt`)
- **Pandas** (installed via `requirements.txt`)
- **Windows executables** (included in `Files_for_GUI/`):
  - `vina.exe` - AutoDock Vina
  - `autogrid4.exe` - AutoGrid4 for AD4 map generation
  - `autodock4.exe` - AutoDock4 for docking

## 📁 Project Structure

```
MetalloDock GUI GitHub/
├── MetalloDock.py          # Main GUI application
├── Files_for_GUI/          # Executables and parameter files
│   ├── vina.exe            # AutoDock Vina executable
│   ├── autogrid4.exe       # AutoGrid4 executable
│   ├── autodock4.exe       # AutoDock4 executable
│   ├── zinc_pseudo.py      # Zinc pseudoatom script
│   ├── AD4_parameters.dat  # AD4 base parameters
│   ├── AD4Zn.dat           # Zinc-specific parameters
│   └── Ligands/            # Sample ligand files (PDBQT format)
├── requirements.txt        # Python dependencies
├── run_docking.bat         # Windows launcher script
└── README.md               # This file
```

## 🎯 Features

- **Dual Docking Methods:**
  - **AutoDock Vina**: Fast docking for initial screening
  - **AutoDock4 (AD4)**: High-precision docking with zinc pseudoatom support for metalloproteins

- **Automatic Setup:**
  - Auto-detects executables from `Files_for_GUI/`
  - Auto-detects parameter files
  - Automatic receptor oxygen normalization (O→OA)

- **User-Friendly Interface:**
  - Drag-and-drop receptor upload
  - Batch ligand processing
  - Interactive grid box configuration
  - Real-time progress tracking
  - Results visualization and export

- **Advanced Options:**
  - Metal center auto-detection
  - Custom grid box settings
  - Configurable exhaustiveness and number of poses
  - AD4 map generation and reuse

## 🔧 Installation

### Windows

1. **Install Python:**
   - Download from [python.org](https://www.python.org/downloads/)
   - Make sure to check "Add Python to PATH" during installation

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the GUI:**
   - Double-click `run_docking.bat`, or
   - Run `streamlit run MetalloDock.py` in terminal

### Linux/Mac

**Note:** The included executables are Windows-only. For Linux/Mac, you'll need to:

1. Install AutoDock Vina and AutoGrid4 from source or package managers
2. Update the paths in `MetalloDock.py` or place executables in `Files_for_GUI/` with appropriate names (without `.exe`)

## 📖 Usage

1. **Upload Receptor:**
   - Drag and drop a PDBQT receptor file, or
   - Provide a local path to the receptor file

2. **Prepare Ligands:**
   - Click "Prepare Ligands" to select ligand files
   - Supports multiple formats (will be converted to PDBQT)

3. **Configure Grid Box:**
   - Set center coordinates (X, Y, Z)
   - Set box size (X, Y, Z)
   - Or use "Auto-detect metal center" for automatic positioning

4. **Select Docking Method:**
   - **Vina**: Fast docking
   - **AD4 Maps**: High-precision docking with pre-generated maps

5. **Run Docking:**
   - Click "Start Docking"
   - Monitor progress in real-time
   - View results when complete

6. **Export Results:**
   - Download individual PDBQT files
   - Export CSV summary
   - Download all results as ZIP

## 🔬 Supported Formats

- **Receptors:** PDBQT format (required)
- **Ligands:** PDB, PDBQT, MOL2, SDF (will be converted to PDBQT)

## ⚙️ Configuration

All executables and parameter files are automatically detected from the `Files_for_GUI/` folder. No manual configuration needed!

- Executables: `vina.exe`, `autogrid4.exe`, `autodock4.exe`
- Parameters: `AD4_parameters.dat`, `AD4Zn.dat`
- Scripts: `zinc_pseudo.py`

## 📧 Contact

For all questions and concerns contact Dr. Sivanesan Dakshanamurthy at sd233@georgetown.edu
---





