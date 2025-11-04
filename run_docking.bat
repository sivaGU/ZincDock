@echo off
echo MetalloDock - Quick Start
echo ========================

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo Error: Python not found. Please install Python from python.org
    echo.
    echo To fix this:
    echo 1. Install Python from https://python.org
    echo 2. Make sure to check "Add Python to PATH" during installation
    echo 3. Restart your command prompt after installation
    pause
    exit /b 1
)

REM Check if Streamlit is available
streamlit --version >nul 2>&1
if errorlevel 1 (
    echo Error: Streamlit not found. Installing...
    pip install streamlit pandas
    if errorlevel 1 (
        echo Error: Failed to install Streamlit. Please run manually:
        echo pip install streamlit pandas
        pause
        exit /b 1
    )
)

REM Check if required files exist (look in local directory)
if not exist "Files_for_GUI\vina.exe" (
    echo Error: vina.exe not found in Files_for_GUI directory
    echo Please ensure you are running this from the Final Code folder
    pause
    exit /b 1
)

if not exist "Files_for_GUI\7q0d_Zn_pp.pdbqt" (
    echo Error: Receptor file not found in Files_for_GUI directory
    echo Please ensure you are running this from the Final Code folder
    pause
    exit /b 1
)

if not exist "Files_for_GUI\Ligands" (
    echo Error: Ligands directory not found in Files_for_GUI directory
    echo Please ensure you are running this from the Final Code folder
    pause
    exit /b 1
)

if not exist "MetalloDock.py" (
    echo Error: MetalloDock.py not found in current directory
    echo Please ensure you are running this from the Final Code folder
    pause
    exit /b 1
)

echo Starting MetalloDock GUI...
echo.

REM Run the original GUI
echo Running: streamlit run MetalloDock.py
echo Opening browser at http://localhost:8501
echo.
streamlit run MetalloDock.py

pause
