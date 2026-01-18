#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import csv
import math
import shutil
import zipfile
import random
import subprocess
import platform
import stat
from pathlib import Path
from typing import List, Tuple, Optional, Set, Dict

REPO_ROOT = Path(__file__).resolve().parent

import streamlit as st
import pandas as pd
import argparse
import sys


# ---- ZincDock Theme Palette ----
LIGHT_POWDER_BLUE = "#6B9FC0"     # Darker Light Powder Blue
SOFT_SKY_BLUE = "#4A7FA8"         # Darker Soft Sky Blue
LIGHT_AZURE = "#2D6B94"           # Darker Light Azure
MEDIUM_STEEL_BLUE = "#1C5C8F"     # Darker Medium Steel Blue
DEEP_CERULEAN = "#0D4F7A"         # Darker Deep Cerulean
RICH_TEAL_BLUE = "#004566"        # Darker Rich Teal Blue
MIDNIGHT_AZURE = "#003A5F"        # Darker Midnight Azure
WHITE = "#FFFFFF"

# Demo preset defaults
DEMO_PRESETS = {
    "Carbonic Anhydrase I": {
        "center": (29.951, 0.420, -4.735),
        "size": (16.0, 18.0, 16.0),
    },
    "Carbonic Anhydrase II": {
        "center": (-6.421, 0.342, 17.256),
        "size": (21.0, 21.0, 21.0),
    },
}

DEMO_PARAM_DEFAULTS = {
    "base_exhaustiveness": 64,
    "base_num_modes": 10,
    "output_name": "PFAS_Docking_Results",
    "timeout_mode": "No timeout (recommended)",
    "timeout_s": 300,
    "max_retries": 2,
    "skip_existing": False,
    "exhaustiveness_retry": 1.50,
    "num_modes_retry": 1.25,
}

def render_home_page():
    st.header("Welcome to ZincDock")
    st.write(
        "ZincDock streamlines AutoDock4 map-based docking for metalloprotein projects. "
        "Use this interface to prepare receptors and ligands, build AD4 maps, and run either "
        "classical Vina or AD4 scoring workflows from a single place."
    )
    st.subheader("Highlights")
    st.markdown(
        "- Metalloprotein aware grid generation with optional zinc pseudo atoms\n"
        "- AD4 scoring workflow including intermolecular, internal, and torsional components\n"
        "- Auto-detected executables and reproducible working directory layout\n"
        "- Shared engine for both the Streamlit GUI and the CLI presets"
    )
    st.subheader("Workflow Overview")
    st.markdown(
        "1. Choose a working directory (ZincDock will create the required folders there).\n"
        "2. Load a receptor from file or path; normalization and optional pseudo atom insertion run automatically.\n"
        "3. Prepare ligands from a source folder or upload ready-to-dock PDBQT files.\n"
        "4. Build or update AD4 maps so that every ligand atom type is covered.\n"
        "5. Pick a docking mode (Vina box search or AD4 maps) and start the run.\n"
        "6. Review scores in the table, download the CSV summary, or fetch all pose PDBQTs as a ZIP."
    )
    st.subheader("Output Guide")
    st.markdown(
        "| Folder | Description |\n"
        "| --- | --- |\n"
        "| `prepared_ligands/ligands_no_hydrogens/` | Ligands copied or prepared for docking |\n"
        "| `ad4_maps/<prefix>/` | AutoGrid4 maps and supporting parameter files |\n"
        "| `<work_dir>/Vina_Docking_Results/` | Example Vina output folder (user-defined name in UI) |\n"
        "| `<work_dir>/AD4_Docking_Results/` | Example AD4 output folder (user-defined name in UI) |\n"
        "| `<work_dir>/<run>_results.csv` | Aggregated scores shown in the GUI |"
    )
    st.subheader("Automation Tip")
    st.write(
        "The CLI entry point (`python ZincDock.py --cli ...`) shares the same code paths as the GUI. "
        "Use it to script batch runs after maps are prepared."
    )

def render_documentation_page():
    st.header("ZincDock — Documentation")
    st.write(
        "Welcome to the documentation tab! Below is a detailed, step-by-step guide to the "
        "ZincDock workflow so you know what each section does and how to use it effectively."
    )

    st.subheader("① Prepare Your Workspace")
    st.markdown(
        "**1. Choose a working directory.**\\n"
        "ZincDock creates `prepared_ligands/`, `ad4_maps/`, and `outputs/` inside the folder you set. "
        "If you run on Streamlit Cloud, the directory defaults to `/mount/src/metallodock/`."
    )
    st.markdown(
        "**2. Review the navigation tabs.**\\n"
        "- *Demo*: AD4 workflow with carbonic anhydrase presets.\\n"
        "- *Standard AutoDock*: Vina box docking.\\n"
        "- *Metalloprotein Docking*: Manual AD4 configuration."
    )

    st.subheader("② Provide Receptor & Ligands")
    st.markdown(
        "**1. Upload a receptor (PDBQT).**\\n"
        "Use the uploader or a path. ZincDock normalizes oxygen labels (O → OA) and keeps coordinates intact."
    )
    st.markdown(
        "**2. Prepare ligands.**\\n"
        "- Supply a source folder and click *Prepare ligands* to convert to PDBQT.\\n"
        "- Or switch to *Upload now* to drop ready-made PDBQT files. Prepared ligands live under `prepared_ligands/ligands_no_hydrogens`."
    )

    st.subheader("③ Configure Grid & Backend")
    st.markdown(
        "**1. Select the docking backend.**\\n"
        "- *Vina (box)*: exhaustiveness-based sampling inside the defined box.\\n"
        "- *AD4 (maps)*: requires AutoGrid maps and supports component energy breakdown."
    )
    st.markdown(
        "**2. Set grid box parameters.**\\n"
        "Enter center (x, y, z), size (Å), and grid spacing. In the Demo tab these are locked to the preset you choose." 
    )
    st.markdown(
        "**3. Force extra atom types (optional).**\\n"
        "Add comma-separated atom symbols (e.g., `S,NA`) if your ligands contain uncommon types that need maps."
    )

    st.subheader("④ Generate AD4 Maps (when using AD4)")
    st.markdown(
        "ZincDock wraps AutoGrid4 to create or update map files. The workflow validates inputs before launching the executable:\\n"
        "1. Confirms `autogrid4` exists and has execute permissions.\\n"
        "2. Merges `AD4_parameters.dat` with optional `AD4Zn.dat`.\\n"
        "3. Runs `zinc_pseudo.py` (if present) to insert tetrahedral Zn pseudoatoms.\\n"
        "4. Normalizes receptor oxygen labels to OA.\\n"
        "5. Detects receptor & ligand atom types and unions them with forced types.\\n"
        "6. Builds the grid parameter file (GPF) and executes AutoGrid4."
    )
    st.warning(
        "Spacing must be greater than 0 Å. The Demo tab locks spacing at 0.375 Å to mimic published CA binding boxes."
    )

    st.subheader("⑤ Run Docking")
    st.markdown(
        "**1. Click *Run Docking*.**\\n"
        "ZincDock queues each ligand, calls the appropriate executable (Vina or AD4), captures stdout/stderr, and displays live status." 
    )
    st.markdown(
        "**2. Understand the results table.**\\n"
        "For every ligand you'll see binding affinity, pose counts, output/log paths, and status. AD4 runs also surface intermolecular, "
        "internal, torsional, and estimated free-energy components." 
    )

    st.subheader("⑥ Review & Export Outputs")
    st.markdown(
        "After docking completes you can:"
        "- Download individual ligand PDBQT and log files."
        "- Export a CSV summary of all ligands."
        "- Generate a ZIP archive with all outputs."
    )
    st.markdown(
        "Grid maps reside in `ad4_maps/<prefix>/` and are reused automatically if they already exist."
    )

    st.subheader("Demo Tab Notes")
    st.markdown(
        "The Demo tab is pre-populated for carbonic anhydrase receptors. Download the bundled folders (`Carbonic Anhydrase Receptor Files` and "
        "`18 PFAS Ligands`) from the repository so the tab can locate receptors and sample ligands. Switching between *Carbonic Anhydrase I* and "
        "*II* locks grid centers, box sizes, spacing (0.375 Å), and docking parameters accordingly."
    )

    st.info(
        "Tip: Use the *Tools → Test executables* button to confirm Vina, AutoGrid4, and AutoDock4 paths before starting long jobs."
    )

# ==============================
# Small helpers
# ==============================

def _save_uploaded_file(uploaded_file, dst_dir: Path) -> Path:
    dst_dir.mkdir(parents=True, exist_ok=True)
    out_path = dst_dir / uploaded_file.name
    with open(out_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return out_path

@st.cache_data(show_spinner=False)
def _cached_file_bytes(b: bytes) -> bytes:
    return b

def autodetect_metal_center(receptor_path: Path, metals=("ZN","MG","MN","FE","CU","CO","NI")) -> Optional[Tuple[float,float,float]]:
    try:
        with open(receptor_path, "r", errors="ignore") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                element = line[76:78].strip().upper() if len(line) >= 78 else ""
                if not element:
                    nm = line[12:16].strip().upper()
                    element = ''.join([c for c in nm if c.isalpha()])[:2]
                if element in metals:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    return (x, y, z)
    except Exception:
        pass
    return None

def parse_binding_affinity(pdbqt: Path) -> str:
    """Parse binding affinity from Vina output PDBQT file.
    Supports both Vina and AD4 scoring formats."""
    try:
        with open(pdbqt, "r", errors="ignore") as f:
            for ln in f:
                # Try Vina format first: "REMARK VINA RESULT:   -7.2   0.000   0.000"
                if ln.startswith("REMARK VINA RESULT:"):
                    parts = ln.split()
                    if len(parts) >= 4:
                        return parts[3]
                # Try AD4 format: "REMARK AD4 RESULT:   -7.2   0.000   0.000"
                elif ln.startswith("REMARK AD4 RESULT:"):
                    parts = ln.split()
                    if len(parts) >= 4:
                        return parts[3]
                # Try generic format: "REMARK RESULT:   -7.2   0.000   0.000"
                elif ln.startswith("REMARK RESULT:"):
                    parts = ln.split()
                    if len(parts) >= 3:
                        return parts[2]
                # Also check for estimated binding energy in AD4 verbose output
                elif "Estimated Free Energy of Binding" in ln:
                    parts = ln.split()
                    # Look for the number after "Binding"
                    for i, part in enumerate(parts):
                        if "Binding" in part and i + 1 < len(parts):
                            try:
                                return str(float(parts[i + 1]))
                            except (ValueError, IndexError):
                                pass
    except Exception:
        pass
    return "N/A"

def count_poses(pdbqt: Path) -> int:
    try:
        with open(pdbqt, "r", errors="ignore") as f:
            return f.read().count("MODEL")
    except Exception:
        return 0

def results_to_csv_bytes(rows: List[dict]) -> bytes:
    buf = io.StringIO()
    # Auto-detect fields from first row
    if rows and isinstance(rows[0], dict):
        fields = list(rows[0].keys())
    else:
        fields = ["Ligand","Binding_Affinity","Num_Poses","Output_File","Log_File","Status"]
    w = csv.DictWriter(buf, fieldnames=fields)
    w.writeheader(); w.writerows(rows)
    return buf.getvalue().encode("utf-8")

def zip_outputs(folder: Path) -> bytes:
    mem = io.BytesIO()
    with zipfile.ZipFile(mem, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in folder.rglob("*"):
            if p.is_file():
                zf.write(p, arcname=p.relative_to(folder))
    mem.seek(0)
    return mem.read()

# ==============================
# Atom-type utilities (diagnostics + normalization)
# ==============================

def read_types_from_pdbqt(pdbqt_path: Path) -> List[str]:
    """Collect unique AD4 atom types (last token) from ATOM/HETATM lines."""
    seen, types = set(), []
    try:
        with open(pdbqt_path, "r", errors="ignore") as f:
            for ln in f:
                if ln.startswith(("ATOM", "HETATM")):
                    toks = ln.split()
                    if toks:
                        t = toks[-1]
                        if t not in seen:
                            seen.add(t); types.append(t)
    except Exception:
        pass
    return types

def ligand_types_union(files: List[Path]) -> Set[str]:
    tset: Set[str] = set()
    for lf in files:
        tset.update(read_types_from_pdbqt(lf))
    return tset

def normalize_receptor_oxygen_to_OA(src: Path, dst: Path) -> None:
    """Copy receptor PDBQT, replacing bare 'O' AD4 atom type with 'OA'."""
    tmp = dst if dst != src else src.with_suffix(".tmp_norm.pdbqt")
    changed = False
    out_lines = []
    with open(src, "r", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                toks = ln.split()
                if toks and toks[-1] == "O":
                    toks[-1] = "OA"
                    ln = " ".join(toks) + "\n"
                    changed = True
            out_lines.append(ln)
    with open(tmp, "w", encoding="utf-8") as g:
        g.writelines(out_lines)
    if dst != src and tmp != dst:
        shutil.copy2(tmp, dst)
    elif dst == src and changed:
        src.write_text("".join(out_lines), encoding="utf-8")

# ==============================
# Ligand preparation
# ==============================

def prepare_ligands_from_folder(source_dir: Path, prepared_root: Path) -> List[Path]:
    """Copies *.pdbqt from source_dir → prepared_root/prepared_ligands/ligands_no_hydrogens as *_prepared.pdbqt."""
    if not source_dir.exists():
        raise FileNotFoundError(f"Ligand directory not found: {source_dir}")
    ligs = sorted(source_dir.glob("*.pdbqt"))
    if not ligs:
        raise FileNotFoundError(f"No .pdbqt ligands found in {source_dir}")

    out_dir = prepared_root / "prepared_ligands" / "ligands_no_hydrogens"
    out_dir.mkdir(parents=True, exist_ok=True)

    prepared = []
    for lig in ligs:
        dst = out_dir / f"{lig.stem}_prepared.pdbqt"
        shutil.copy2(lig, dst)
        prepared.append(dst)
    return prepared

# ==============================
# AD4Zn helpers (params merge, TZ, GPF, AutoGrid, map checks)
# ==============================

def merge_parameter_files(main_dat: Path, extra_dat: Optional[Path], out_dat: Path) -> None:
    out_dat.parent.mkdir(parents=True, exist_ok=True)
    with open(out_dat, "w", encoding="utf-8") as out:
        if main_dat and main_dat.exists():
            txt = main_dat.read_text(encoding="utf-8")
            out.write(txt)
            if not txt.endswith("\n"):
                out.write("\n")
        if extra_dat and extra_dat.exists():
            out.write("\n# --- Extra parameters appended ---\n")
            out.write(extra_dat.read_text(encoding="utf-8"))

def run_zinc_pseudo(python_exe: Path, script: Path, receptor_in: Path, receptor_tz_out: Path) -> subprocess.CompletedProcess:
    # Ensure receptor file exists and is accessible
    if not receptor_in.exists():
        raise FileNotFoundError(f"Receptor file not found: {receptor_in}")
    
    # Create output directory if it doesn't exist
    receptor_tz_out.parent.mkdir(parents=True, exist_ok=True)
    
    # Use the current Python interpreter instead of the broken one in Files_for_GUI
    import sys
    current_python = sys.executable
    
    # Run zinc_pseudo.py with absolute paths to avoid path issues
    return subprocess.run(
        [current_python, str(script), "-r", str(receptor_in.absolute()), "-o", str(receptor_tz_out.absolute())],
        capture_output=True, text=True, timeout=300
    )

def write_simple_gpf(
    gpf_path: Path,
    receptor_tz_filename: str,
    maps_prefix_basename: str,
    npts_xyz: Tuple[int,int,int],
    spacing: float,
    center_xyz: Tuple[float,float,float],
    receptor_types: List[str],
    ligand_types: List[str],
    parameter_file_rel: str,
) -> None:
    nx, ny, nz = npts_xyz
    cx, cy, cz = center_xyz
    rt = " ".join(receptor_types)
    lt = " ".join(ligand_types)
    lines = [
        f"npts {nx} {ny} {nz}",
        f"parameter_file {parameter_file_rel}",
        f"gridfld {maps_prefix_basename}.maps.fld",
        f"spacing {spacing:.3f}",
        f"receptor_types {rt}",
        f"ligand_types {lt}",
        f"receptor {receptor_tz_filename}",
        f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}",
        "smooth 0.5",
        f"elecmap {maps_prefix_basename}.e.map",
        f"dsolvmap {maps_prefix_basename}.d.map",
    ]
    for t in ligand_types:
        lines.append(f"map {maps_prefix_basename}.{t}.map")
    gpf_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def run_autogrid4(autogrid_exe: Path, work_dir: Path, gpf_path: Path, timeout_s: int = 1800) -> subprocess.CompletedProcess:
    """Run AutoGrid4 with proper error handling and permission checks."""
    if not autogrid_exe or not autogrid_exe.exists():
        raise FileNotFoundError(f"AutoGrid4 executable not found at: {autogrid_exe}")
    
    # Check if it's a Windows .exe on a non-Windows system FIRST (most important check)
    is_windows_os = platform.system() == "Windows"
    
    # Check filename for .exe extension (case-insensitive check)
    exe_name_lower = str(autogrid_exe.name).lower()
    has_exe_extension = exe_name_lower.endswith('.exe') or autogrid_exe.suffix.lower() == '.exe'
    
    if not is_windows_os and has_exe_extension:
        # Detect if running on Streamlit Cloud
        is_streamlit_cloud = os.environ.get("STREAMLIT_SERVER_URL", "").startswith("https://") or os.environ.get("STREAMLIT_SHARE", "") != ""
        cloud_context = "Streamlit Cloud runs on Linux servers" if is_streamlit_cloud else "This system runs on Linux"
        
        raise PermissionError(
            f"Windows executable (.exe) detected on Linux system.\n\n"
            f"**Why Linux?** {cloud_context}, so Windows .exe files cannot run here.\n\n"
            f"**Solution:**\n"
            f"1. Download Linux version of AutoGrid4 from AutoDock website\n"
            f"2. Place it in Files_for_GUI/ as 'autogrid4' (without .exe extension)\n"
            f"3. Commit and push to GitHub - Streamlit Cloud will automatically detect it\n"
            f"4. The Linux executable will work on Streamlit Cloud"
        )
    
    # Check if executable (on Unix-like systems) - only for Linux executables
    # This check should only run if we've confirmed it's NOT a Windows .exe
    if not is_windows_os and not has_exe_extension:
        if not os.access(autogrid_exe, os.X_OK):
            # Try to automatically fix execute permissions - more aggressive fix
            try:
                # First try the standard chmod
                current_stat = os.stat(autogrid_exe)
                new_mode = current_stat.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
                os.chmod(autogrid_exe, new_mode)
                
                # Verify it worked
                if os.access(autogrid_exe, os.X_OK):
                    # Successfully fixed, continue
                    pass
                else:
                    # Try setting permissions to 755 explicitly
                    os.chmod(autogrid_exe, 0o755)
                    if not os.access(autogrid_exe, os.X_OK):
                        raise PermissionError(
                            f"AutoGrid4 executable is not executable and could not be fixed: {autogrid_exe}\n"
                            f"The file exists but lacks execute permissions. This is a Streamlit Cloud limitation.\n"
                            f"Please ensure the files have execute permissions in Git before pushing to GitHub."
                        )
            except (OSError, PermissionError) as e:
                # Could not fix permissions automatically - this happens on some cloud platforms
                raise PermissionError(
                    f"AutoGrid4 executable is not executable: {autogrid_exe}\n"
                    f"Attempted to fix automatically but failed due to: {str(e)}\n\n"
                    f"**Solution:** Ensure files have execute permissions in Git before pushing:\n"
                    f"1. Run: git update-index --chmod=+x Files_for_GUI/autogrid4\n"
                    f"2. Commit and push to GitHub\n"
                    f"3. Streamlit Cloud will preserve the execute permissions"
                )
    
    try:
        return subprocess.run(
            [str(autogrid_exe), "-p", gpf_path.name, "-l", gpf_path.with_suffix(".glg").name],
            cwd=str(work_dir), capture_output=True, text=True, timeout=timeout_s
        )
    except PermissionError as e:
        raise PermissionError(
            f"Permission denied when trying to run AutoGrid4: {autogrid_exe}\n"
            f"On Linux/Mac, make sure the file has execute permissions: chmod +x {autogrid_exe}\n"
            f"Original error: {e}"
        )

def list_maps_present(maps_prefix: Path) -> Set[str]:
    """Return set of atom types that already have an affinity map file for this prefix."""
    present = set()
    folder = maps_prefix.parent
    base = maps_prefix.name
    for p in folder.glob(f"{base}.*.map"):
        # expecting base.<TYPE>.map
        t = p.suffixes[-2].lstrip(".") if len(p.suffixes) >= 2 else None
        if t:
            present.add(t)
    return present

# ==============================
# Docking (no-timeout option + retries/backoff + console prints)
# ==============================

def _vina_cmd(
    vina_exe: Path,
    mode: str,  # "ad4" or "vina"
    receptor_file: Path,
    ligand_file: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    exhaustiveness: int,
    num_modes: int,
    seed: Optional[int],
    maps_prefix: Optional[Path],
) -> List[str]:
    cx, cy, cz = center
    sx, sy, sz = size
    cmd = [str(vina_exe)]
    
    # Note: Standard Vina doesn't support --maps or --scoring ad4 options
    # If AD4 mode is requested but vina doesn't support it, fall back to regular vina with box
    # This ensures docking works even if AD4 maps aren't supported by the vina executable
    if mode == "ad4" and maps_prefix:
        # Try AD4 maps format first (for modified vina versions that support it)
        # But if this fails, the caller will detect it and can fall back
        cmd += ["--ligand", str(ligand_file), "--maps", str(maps_prefix), "--scoring", "ad4"]
    else:
        # Standard Vina format with receptor and search box
        cmd += ["--receptor", str(receptor_file),
                "--ligand", str(ligand_file),
                "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
                "--scoring", "vina"]
    
    cmd += ["--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    return cmd

def _parse_missing_map(stderr_txt: str) -> Optional[str]:
    # Vina message for AD4: 'Affinity map for atom type X is not present.'
    for line in (stderr_txt or "").splitlines():
        line = line.strip()
        if "Affinity map for atom type" in line and "not present" in line:
            # extract last token before 'is not present.'
            parts = line.replace("Affinity map for atom type", "").replace("is not present.", "").strip()
            # parts should be like: 'S'
            token = parts.strip(" :\"'").split()[-1] if parts else None
            return token
    return None

def _run_one(
    vina_exe: Path,
    mode: str,
    receptor_file: Path,
    ligand_file: Path,
    out_pdbqt: Path,
    log_file: Path,
    center: Tuple[float,float,float],
    size: Tuple[float,float,float],
    exhaustiveness: int,
    num_modes: int,
    seed: Optional[int],
    timeout_s: Optional[int],
    maps_prefix: Optional[Path],
) -> Tuple[bool, str, int, Optional[str], Optional[dict]]:
    """Return (ok, affinity, nposes, missing_atom_type, ad4_components). Writes log regardless."""
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = _vina_cmd(vina_exe, mode, receptor_file, ligand_file, center, size,
                    exhaustiveness, num_modes, seed, maps_prefix)
    cmd += ["--out", str(out_pdbqt)]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=None if (timeout_s is None or timeout_s == 0) else int(timeout_s)
        )
        with open(log_file, "w", encoding="utf-8") as lf:
            lf.write(f"Command: {' '.join(cmd)}\n")
            lf.write(f"Return code: {proc.returncode}\n")
            lf.write("\n---- STDOUT ----\n")
            if proc.stdout: lf.write(proc.stdout)
            lf.write("\n---- STDERR ----\n")
            if proc.stderr: lf.write(proc.stderr)

        if proc.returncode != 0:
            missing = _parse_missing_map(proc.stderr or "")
            error_msg = proc.stderr or proc.stdout or "Unknown error"
            
            # Check if the error is because vina doesn't support --maps option
            # If so, fall back to regular vina with box coordinates
            if "unknown option maps" in error_msg.lower() or "unknown option" in error_msg.lower():
                with open(log_file, "a", encoding="utf-8") as lf:
                    lf.write(f"\n---- WARNING: Vina doesn't support --maps option ----\n")
                    lf.write(f"Falling back to regular Vina with receptor and search box\n")
                
                # Retry with regular vina mode (receptor + box)
                if mode == "ad4" and maps_prefix:
                    # Fall back to vina mode with receptor file
                    cmd = _vina_cmd(vina_exe, "vina", receptor_file, ligand_file, center, size,
                                  exhaustiveness, num_modes, seed, None)
                    cmd += ["--out", str(out_pdbqt)]
                    
                    # Try again with regular vina mode
                    proc_fallback = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=None if (timeout_s is None or timeout_s == 0) else int(timeout_s)
                    )
                    
                    with open(log_file, "a", encoding="utf-8") as lf:
                        lf.write(f"\n---- FALLBACK ATTEMPT ----\n")
                        lf.write(f"Command: {' '.join(cmd)}\n")
                        lf.write(f"Return code: {proc_fallback.returncode}\n")
                        if proc_fallback.stdout: lf.write(f"STDOUT: {proc_fallback.stdout}\n")
                        if proc_fallback.stderr: lf.write(f"STDERR: {proc_fallback.stderr}\n")
                    
                    if proc_fallback.returncode == 0 and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                        # Success with fallback
                        aff = parse_binding_affinity(out_pdbqt)
                        nposes = count_poses(out_pdbqt)
                        if aff not in ("", "N/A") and nposes > 0:
                            with open(log_file, "a", encoding="utf-8") as lf:
                                lf.write(f"\n---- FALLBACK SUCCESS ----\n")
                                lf.write(f"Binding affinity: {aff}\n")
                                lf.write(f"Number of poses: {nposes}\n")
                            return (True, aff, nposes, None, None)
            
            # Log original error details
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR DETAILS ----\n")
                lf.write(f"Return code: {proc.returncode}\n")
                lf.write(f"Missing map: {missing}\n")
                lf.write(f"Error: {error_msg[:500]}\n")
            return (False, "", 0, missing, None)

        if not out_pdbqt.exists():
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR: Output file not created ----\n")
                lf.write(f"Expected: {out_pdbqt}\n")
            return (False, "", 0, None, None)

        if out_pdbqt.stat().st_size == 0:
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR: Output file is empty ----\n")
                lf.write(f"File: {out_pdbqt}\n")
            return (False, "", 0, None, None)

        aff = parse_binding_affinity(out_pdbqt)
        nposes = count_poses(out_pdbqt)
        ad4_components = None
        if mode == "ad4":
            try:
                ad4_components = parse_ad4_verbose_output(proc.stdout or "")
            except Exception:
                ad4_components = None
        
        # Log parsing results
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write(f"\n---- PARSING RESULTS ----\n")
            lf.write(f"Binding affinity: {aff}\n")
            lf.write(f"Number of poses: {nposes}\n")
            if ad4_components:
                for key, val in ad4_components.items():
                    lf.write(f"{key}: {val}\n")
        
        if aff in ("", "N/A") or nposes == 0:
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- WARNING: No valid scores found ----\n")
                lf.write(f"Affinity: '{aff}' (empty/N/A)\n")
                lf.write(f"Poses: {nposes}\n")
                # Try to read first few lines of output file for debugging
                try:
                    with open(out_pdbqt, "r", errors="ignore") as f:
                        first_lines = ''.join(f.readlines()[:20])
                        lf.write(f"\nFirst 20 lines of output:\n{first_lines}\n")
                except:
                    pass
            return (False, "", nposes, None, None)

        return (True, aff, nposes, None, ad4_components)

    except subprocess.TimeoutExpired as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- TIMEOUT ----\n")
            lf.write(str(e))
        return (False, "", 0, None, None)
    except Exception as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- EXCEPTION ----\n")
            lf.write(str(e))
        return (False, "", 0, None, None)

def run_vina_batch(
    vina_exe: Path,
    receptor_file: Path,
    ligand_files: List[Path],
    out_dir: Path,
    center: Tuple[float,float,float],
    size: Tuple[float,float,float],
    scoring: str,  # "vina" or "ad4"
    base_exhaustiveness: int,
    base_num_modes: int,
    timeout_mode: str,           # "no_timeout" or "soft_timeout"
    timeout_s: int,              # ignored if no_timeout
    max_retries: int,
    exhu_backoff: float,
    modes_backoff: float,
    progress_cb=None,
    maps_prefix: Optional[Path] = None,
    skip_if_output_exists: bool = False,
) -> List[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    mode = "ad4" if (scoring == "ad4" and maps_prefix) else "vina"
    rows = []

    for i, lig in enumerate(ligand_files, start=1):
        lig_name = lig.stem.replace("_prepared_no_h", "").replace("_prepared", "")
        suffix = "ad4_out" if mode == "ad4" else "vina_out"
        out_pdbqt = out_dir / f"{lig_name}_{suffix}.pdbqt"
        log_file = out_dir / f"{lig_name}.log"

        if skip_if_output_exists and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            aff = parse_binding_affinity(out_pdbqt); nposes = count_poses(out_pdbqt)
            rows.append({
                "Ligand": lig_name, "Binding_Affinity": aff, "Num_Poses": nposes,
                "Output_File": str(out_pdbqt), "Log_File": str(log_file), "Status": "Skipped (exists)"
            })
            if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Skipped (existing) | Score {aff}")
            continue

        tried = 0
        ok, aff, nposes = False, "", 0
        ex = base_exhaustiveness
        nm = base_num_modes
        per_try_timeout = None if timeout_mode == "no_timeout" else int(timeout_s)
        last_missing = None

        ad4_components = None
        while tried <= max_retries and not ok:
            seed = random.randint(1, 2**31-1)
            if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Running (try {tried+1}/{max_retries+1})")
            ok, aff, nposes, last_missing, ad4_components = _run_one(
                vina_exe, mode, receptor_file, lig, out_pdbqt, log_file,
                center, size, ex, nm, seed, per_try_timeout, maps_prefix
            )

            if ok:
                if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Success | Score {aff} ({nposes} poses)")
                break

            # report immediate cause if missing map
            if last_missing and progress_cb:
                progress_cb(i, len(ligand_files), lig_name, f"Failed - Missing map: {last_missing}")

            tried += 1
            ex = max(ex, int(math.ceil(ex * exhu_backoff)))
            nm = max(nm, int(math.ceil(nm * modes_backoff)))

        status = "Success" if ok else ("Failed - Timeout" if timeout_mode != "no_timeout" else f"Failed{(' - Missing '+last_missing) if last_missing else ''}")
        row = {
            "Ligand": lig_name,
            "Binding_Affinity": aff if ok else "",
            "Num_Poses": nposes if ok else 0,
            "Output_File": str(out_pdbqt if ok else ""),
            "Log_File": str(log_file),
            "Status": status,
        }
        if ok and ad4_components:
            row.update(ad4_components)
        rows.append(row)
        # final console line for this ligand
        if progress_cb:
            if ok:
                progress_cb(i, len(ligand_files), lig_name, f"Done | Score {aff}")
            else:
                if last_missing:
                    progress_cb(i, len(ligand_files), lig_name, f"FAILED | Missing map: {last_missing}")
                else:
                    progress_cb(i, len(ligand_files), lig_name, "FAILED")

    return rows

# ==============================
# GNINA Docking Functions
# ==============================

def find_smina_executable() -> Optional[Path]:
    """Find SMINA executable in PATH or conda environment."""
    # Check PATH
    exe = shutil.which("smina")
    if exe:
        return Path(exe)
    
    # Check conda environment
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "equibind_cpu")
    try:
        result = subprocess.run(
            ["conda", "run", "-n", conda_env, "which", "smina"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            return Path(result.stdout.strip())
    except:
        pass
    
    # Check common conda locations
    conda_base = os.environ.get("CONDA_PREFIX", "")
    if conda_base:
        is_windows_platform = platform.system() == "Windows"
        smina_path = Path(conda_base) / "Library" / "bin" / "smina.exe" if is_windows_platform else Path(conda_base) / "bin" / "smina"
        if smina_path.exists():
            return smina_path
    
    return None


def parse_gnina_affinities(stdout: str, num_modes: int = 10) -> List[float]:
    """Parse binding affinities from GNINA/SMINA stdout output."""
    import re
    affinities = []
    
    # Try parsing from stdout - SMINA outputs a table like:
    # mode |   affinity | dist from best mode
    #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
    # -----+------------+----------+----------
    # 1       -7.2       0.000      0.000
    for line in stdout.split('\n'):
        # Look for lines that start with a number (mode number) followed by affinity
        match = re.match(r'^\s*(\d+)\s+(-?\d+\.?\d*)', line)
        if match:
            try:
                affinity = float(match.group(2))
                affinities.append(affinity)
            except ValueError:
                pass
    
    # Remove duplicates and sort (best affinity first, most negative)
    affinities = sorted(set(affinities))[:num_modes]
    return affinities


def run_gnina_one(
    smina_exec: Path,
    receptor_file: Path,
    ligand_file: Path,
    out_pdbqt: Path,
    log_file: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    exhaustiveness: int,
    num_modes: int,
    seed: Optional[int],
    timeout_s: Optional[int],
    cnn_scoring: str = "none",  # "none", "rescore", "refinement", "all"
) -> Tuple[bool, str, int]:
    """Run GNINA (via SMINA) for a single ligand. Returns (ok, affinity, nposes)."""
    import re
    
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    
    # Build SMINA command
    cmd = [
        str(smina_exec),
        "-r", str(receptor_file),
        "-l", str(ligand_file),  # SMINA supports SDF, MOL2, PDBQT directly
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--num_modes", str(num_modes),
        "--exhaustiveness", str(exhaustiveness),
        "-o", str(out_pdbqt)
    ]
    
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    
    # CNN scoring options (if GNINA is available, otherwise SMINA will ignore)
    if cnn_scoring != "none":
        if cnn_scoring == "rescore":
            cmd.append("--cnn")
        elif cnn_scoring == "refinement":
            cmd.extend(["--cnn_scoring", "refinement"])
        elif cnn_scoring == "all":
            cmd.extend(["--cnn_scoring", "all"])
    
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace',
            timeout=None if (timeout_s is None or timeout_s == 0) else int(timeout_s)
        )
        
        with open(log_file, "w", encoding="utf-8") as lf:
            lf.write(f"Command: {' '.join(cmd)}\n")
            lf.write(f"Return code: {proc.returncode}\n")
            lf.write("\n---- STDOUT ----\n")
            if proc.stdout: lf.write(proc.stdout)
            lf.write("\n---- STDERR ----\n")
            if proc.stderr: lf.write(proc.stderr)
        
        if proc.returncode != 0:
            return (False, "", 0)
        
        if not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
            return (False, "", 0)
        
        # Parse binding affinity from output
        affinities = parse_gnina_affinities(proc.stdout or "", num_modes)
        
        # Also try parsing from PDBQT file
        if out_pdbqt.exists() and not affinities:
            with open(out_pdbqt, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                for line in content.split('\n'):
                    if 'REMARK' in line and 'Affinity' in line:
                        matches = re.findall(r'[-]?\d+\.?\d*', line)
                        for match in matches:
                            try:
                                val = float(match)
                                if val < 0:
                                    affinities.append(val)
                            except ValueError:
                                pass
        
        aff = affinities[0] if affinities else ""
        aff_str = f"{aff:.4f}" if aff != "" else ""
        nposes = count_poses(out_pdbqt) if out_pdbqt.exists() else len(affinities)
        
        return (True, aff_str, nposes)
        
    except subprocess.TimeoutExpired as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- TIMEOUT ----\n")
            lf.write(str(e))
        return (False, "", 0)
    except Exception as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- EXCEPTION ----\n")
            lf.write(str(e))
        return (False, "", 0)


def run_gnina_batch(
    smina_exec: Path,
    receptor_file: Path,
    ligand_files: List[Path],
    out_dir: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    base_exhaustiveness: int,
    base_num_modes: int,
    timeout_mode: str,  # "no_timeout" or "soft_timeout"
    timeout_s: int,
    max_retries: int,
    exhu_backoff: float,
    modes_backoff: float,
    progress_cb=None,
    skip_if_output_exists: bool = False,
    cnn_scoring: str = "none",
) -> List[dict]:
    """Run GNINA (via SMINA) docking for multiple ligands."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    
    for i, lig in enumerate(ligand_files, start=1):
        lig_name = lig.stem.replace("_prepared_no_h", "").replace("_prepared", "")
        out_pdbqt = out_dir / f"{lig_name}_gnina_out.pdbqt"
        log_file = out_dir / f"{lig_name}.log"
        
        if skip_if_output_exists and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            aff = parse_binding_affinity(out_pdbqt) if out_pdbqt.exists() else ""
            nposes = count_poses(out_pdbqt) if out_pdbqt.exists() else 0
            rows.append({
                "Ligand": lig_name, "Binding_Affinity": aff, "Num_Poses": nposes,
                "Output_File": str(out_pdbqt), "Log_File": str(log_file), "Status": "Skipped (exists)"
            })
            if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Skipped (existing) | Score {aff}")
            continue
        
        tried = 0
        ok, aff, nposes = False, "", 0
        ex = base_exhaustiveness
        nm = base_num_modes
        per_try_timeout = None if timeout_mode == "no_timeout" else int(timeout_s)
        
        while tried <= max_retries and not ok:
            seed = random.randint(1, 2**31-1)
            if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Running (try {tried+1}/{max_retries+1})")
            ok, aff, nposes = run_gnina_one(
                smina_exec, receptor_file, lig, out_pdbqt, log_file,
                center, size, ex, nm, seed, per_try_timeout, cnn_scoring
            )
            
            if ok:
                if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Success | Score {aff} ({nposes} poses)")
                break
            
            tried += 1
            ex = max(ex, int(math.ceil(ex * exhu_backoff)))
            nm = max(nm, int(math.ceil(nm * modes_backoff)))
        
        status = "Success" if ok else ("Failed - Timeout" if timeout_mode != "no_timeout" else "Failed")
        rows.append({
            "Ligand": lig_name,
            "Binding_Affinity": aff if ok else "",
            "Num_Poses": nposes if ok else 0,
            "Output_File": str(out_pdbqt if ok else ""),
            "Log_File": str(log_file),
            "Status": status,
        })
        
        if progress_cb:
            if ok:
                progress_cb(i, len(ligand_files), lig_name, f"Done | Score {aff}")
            else:
                progress_cb(i, len(ligand_files), lig_name, "FAILED")
    
    return rows

# ==============================
# AD4 Analysis Helpers
# ==============================

def extract_first_pose_simple(input_pdbqt: Path, output_pdbqt: Path) -> bool:
    """Extract first pose from multi-model PDBQT preserving PDBQT directives."""
    try:
        with open(input_pdbqt, "r") as f:
            lines = f.readlines()

        in_model = False
        pose_lines: List[str] = []

        for line in lines:
            line = line.replace("\x00", "")
            if line.startswith("MODEL"):
                if line.strip().startswith("MODEL 1"):
                    in_model = True
                    continue
                if in_model:
                    # next model started; stop after writing
                    break
                continue
            if line.startswith("ENDMDL"):
                break
            if in_model:
                if line.startswith("REMARK"):
                    continue
                pose_lines.append(line)

        if not pose_lines:
            return False

        with open(output_pdbqt, "w") as f:
            for ln in pose_lines:
                f.write(ln)
        return True
    except Exception:
        return False

def extract_all_poses(input_pdbqt: Path, output_dir: Path, max_poses: int = 10) -> List[Path]:
    """Extract all poses from multi-model PDBQT and save as separate files."""
    extracted_files: List[Path] = []
    try:
        with open(input_pdbqt, "r") as f:
            lines = f.readlines()

        current_pose = 0
        current_lines: List[str] = []
        in_model = False

        for line in lines:
            line = line.replace("\x00", "")
            if line.startswith("MODEL"):
                if in_model:
                    # Unexpected nested MODEL, flush existing lines first
                    if current_lines:
                        pose_file = output_dir / f"pose_{current_pose}.pdbqt"
                        with open(pose_file, "w") as pf:
                            pf.writelines(current_lines)
                        extracted_files.append(pose_file)
                        current_lines = []
                current_pose += 1
                if current_pose > max_poses:
                    break
                in_model = True
                current_lines = []
                continue

            if line.startswith("ENDMDL"):
                if in_model and current_lines:
                    pose_file = output_dir / f"pose_{current_pose}.pdbqt"
                    with open(pose_file, "w") as pf:
                        pf.writelines(current_lines)
                    extracted_files.append(pose_file)
                in_model = False
                current_lines = []
                continue

            if in_model:
                if line.startswith("REMARK"):
                    continue
                current_lines.append(line)

        return extracted_files
    except Exception:
        return []

def parse_ad4_verbose_output(stdout: str) -> dict:
    """Parse AD4 verbose output for energy components"""
    
    result = {
        'AD4_Affinity': None,
        'AD4_Intermolecular': None,
        'AD4_Internal': None,
        'AD4_Torsional': None
    }
    
    for line in stdout.split('\n'):
        if 'Estimated Free Energy of Binding' in line:
            result['AD4_Affinity'] = float(line.split(':')[1].split()[0])
        elif '(1) Final Intermolecular Energy' in line:
            result['AD4_Intermolecular'] = float(line.split(':')[1].split()[0])
        elif '(2) Final Total Internal Energy' in line:
            result['AD4_Internal'] = float(line.split(':')[1].split()[0])
        elif '(3) Torsional Free Energy' in line:
            result['AD4_Torsional'] = float(line.split(':')[1].split()[0])
    
    return result


# ==============================
# Endogenous Docking Presets (functions)
# ==============================

def _endogenous_presets() -> dict:
    root = Path("C:/Users/madas/Downloads/Protocols for Docking")
    return {
        "CA": {
            "receptor": root / "Receptor Files/Carbonic Anhydrase I/Carbonic_Anhydrase_Post_Processed.pdbqt",
            "lig_dir": root / "Endogenous Ligands/CA I",
            "center": (29.951, 0.420, -4.735),
            "size": (16.0, 18.0, 16.0),
        },
        "CA_CUSTOM": {
            "receptor": root / "Receptor Files/Carbonic Anhydrase I/Carbonic_Anhydrase_Post_Processed.pdbqt",
            "lig_dir": root / "ligands_no_hydrogens",
            "center": (29.951, 0.420, -4.735),
            "size": (16.0, 18.0, 16.0),
        },
        "CAI_PFAS": {
            "receptor": root / "Receptor Files/Carbonic Anhydrase I/Carbonic_Anhydrase_Post_Processed.pdbqt",
            "lig_dir": root / "18 PFAS",
            "center": (29.951, 0.420, -4.735),
            "size": (16.0, 18.0, 16.0),
        },
        "CAII_PFAS": {
            "receptor": root / "Receptor Files/Carbonic Anhydrase II/CA_2_pp.pdbqt",
            "lig_dir": root / "18 PFAS",
            "center": (-6.421, 0.342, 17.256),
            "size": (21.0, 21.0, 21.0),
        },
        "SOD1": {
            "receptor": root / "Receptor Files/SOD1 Receptor + Gridbox/5YTU_Cleaned.pdbqt",
            "lig_dir": root / "Endogenous Ligands/SOD1 Ligand",
            "center": (-77.967, 6.755, 1.13),
            "size": (18.0, 20.0, 20.0),
        },
        "HDAC1": {
            "receptor": root / "Receptor Files/HDAC Receptor + Gridobox/1HDAC.pdbqt",
            "lig_dir": root / "Endogenous Ligands/HDAC Ligands",
            "center": (205.989, 159.799, 161.458),
            "size": (20.0, 20.0, 20.0),
        },
        "HDAC2": {
            "receptor": root / "Receptor Files/HDAC Receptor + Gridobox/HDAC2.pdbqt",
            "lig_dir": root / "Endogenous Ligands/HDAC Ligands",
            "center": (19.980, 18.779, 0.606),
            "size": (18.0, 22.0, 20.0),
        },
        "HDAC3": {
            "receptor": root / "Receptor Files/HDAC Receptor + Gridobox/HDAC3.pdbqt",
            "lig_dir": root / "Endogenous Ligands/HDAC Ligands",
            "center": (2.976, 54.243, 24.699),
            "size": (22.0, 20.0, 20.0),
        },
        "HDAC4": {
            "receptor": root / "Receptor Files/HDAC Receptor + Gridobox/4HDAC.pdbqt",
            "lig_dir": root / "Endogenous Ligands/HDAC Ligands",
            "center": (20.225, 9.222, 1.723),
            "size": (22.0, 22.0, 22.0),
        },
        "HDAC6": {
            "receptor": root / "Receptor Files/HDAC Receptor + Gridobox/HDAC6.pdbqt",
            "lig_dir": root / "Endogenous Ligands/HDAC Ligands",
            "center": (-22.063, 20.237, 26.789),
            "size": (24.0, 20.0, 20.0),
        },
    }

def _collect_endogenous_ligands(lig_dir: Path, preset_key: str) -> List[Path]:
    if not lig_dir.exists():
        return []
    all_ligs = list(sorted(lig_dir.glob("*.pdbqt")))
    base_ligs = [p for p in all_ligs if not p.name.lower().endswith("_out.pdbqt")]
    if preset_key.startswith("HDAC"):
        x = preset_key[-1]
        specific = [p for p in base_ligs if f"HDAC{x}" in p.stem.upper()]
        return specific or base_ligs
    return base_ligs

def run_endogenous_preset_ad4(preset_key: str, headless: bool = False) -> List[dict]:
    presets = _endogenous_presets()
    if preset_key not in presets:
        if headless:
            print(f"ERROR: Unknown preset: {preset_key}")
        else:
            st.error("Unknown preset.")
        return []
    cfg = presets[preset_key]

    rec_path = Path(cfg["receptor"]).resolve()
    ligs = _collect_endogenous_ligands(Path(cfg["lig_dir"]).resolve(), preset_key)
    if not rec_path.exists():
        msg = f"Receptor not found: {rec_path}"
        if headless:
            print(f"ERROR: {msg}")
        else:
            st.error(msg)
        return []
    if not ligs:
        msg = f"No ligand PDBQTs found in: {cfg['lig_dir']}"
        if headless:
            print(f"ERROR: {msg}")
        else:
            st.error(msg)
        return []

    maps_dir = (work_dir / "ad4_maps" / preset_key)
    maps_dir.mkdir(parents=True, exist_ok=True)
    maps_prefix = maps_dir / f"{preset_key}_maps"

    if not (base_params and base_params.exists()):
        msg = "AD4_parameters.dat is missing. Set it in Configuration."
        if headless:
            print(f"ERROR: {msg}")
        else:
            st.error(msg)
        return []
    merged_params = maps_dir / "AD4_parameters_plus_ZnTZ.dat"
    merge_parameter_files(base_params, extra_params, merged_params)

    rec_tz = rec_path
    if normalize_OA:
        try:
            normalize_receptor_oxygen_to_OA(rec_tz, rec_tz)
        except Exception:
            pass

    # Filter out invalid AD4 receptor types (ions that AutoGrid4 doesn't support)
    invalid_types = {"K", "Na", "MG", "CA", "CL", "FE", "MN", "ZN", "CU", "CO", "NI"}
    
    # Create a filtered copy of receptor for map building (remove atoms with invalid types)
    rec_filtered = maps_dir / (rec_tz.stem + "_filtered.pdbqt")
    try:
        with open(rec_tz, "r", errors="ignore") as f_in:
            with open(rec_filtered, "w", encoding="utf-8") as f_out:
                for line in f_in:
                    if line.startswith(("ATOM", "HETATM")):
                        toks = line.split()
                        if toks and toks[-1] in invalid_types:
                            continue  # Skip atoms with invalid types
                    f_out.write(line)
        rec_tz = rec_filtered
        if not rec_tz.exists() or rec_tz.stat().st_size == 0:
            if headless:
                print(f"WARNING: Filtered receptor is empty, using original")
            rec_tz = rec_path
    except Exception as e:
        if headless:
            print(f"WARNING: Could not filter receptor: {e}, using original")
        rec_tz = rec_path

    rec_types_detected = read_types_from_pdbqt(rec_tz)
    # Also filter from detected types list
    rec_types_detected = [t for t in rec_types_detected if t not in invalid_types]
    if not rec_types_detected:
        msg = "No valid receptor types detected after filtering."
        if headless:
            print(f"ERROR: {msg}")
        else:
            st.error(msg)
        return []
    lig_types_detected = ligand_types_union(ligs)
    lig_types_full = sorted(set(lig_types_detected))

    spacing_val = float(spacing)
    sx, sy, sz = cfg["size"]
    nx = max(10, int(round(float(sx) / spacing_val)))
    ny = max(10, int(round(float(sy) / spacing_val)))
    nz = max(10, int(round(float(sz) / spacing_val)))
    gpf_out = maps_prefix.with_suffix(".gpf")
    write_simple_gpf(
        gpf_path=gpf_out,
        receptor_tz_filename=rec_tz.name,
        maps_prefix_basename=maps_prefix.name,
        npts_xyz=(nx, ny, nz),
        spacing=spacing_val,
        center_xyz=tuple(map(float, cfg["center"])),
        receptor_types=rec_types_detected,
        ligand_types=lig_types_full,
        parameter_file_rel=merged_params.name,
    )

    # Ensure receptor is in maps directory
    if rec_tz.parent != maps_dir:
        try:
            shutil.copy2(rec_tz, maps_dir / rec_tz.name)
            rec_tz = maps_dir / rec_tz.name
        except Exception:
            pass
    try:
        ex_lig = ligs[0]
        if ex_lig.parent != maps_dir:
            shutil.copy2(ex_lig, maps_dir / ex_lig.name)
    except Exception:
        pass

    if not (autogrid_exe and autogrid_exe.exists()):
        msg = "AutoGrid4 executable path is missing."
        if headless:
            print(f"ERROR: {msg}")
        else:
            st.error(msg)
        return []
    ag = run_autogrid4(autogrid_exe, maps_dir, gpf_out)
    if ag.returncode != 0:
        msg = "AutoGrid4 failed while building maps."
        if headless:
            print(f"ERROR: {msg}")
            print(f"AutoGrid4 stdout:\n{ag.stdout or ''}")
            print(f"AutoGrid4 stderr:\n{ag.stderr or ''}")
        else:
            st.error(msg)
            try:
                st.code((ag.stdout or '') + "\n" + (ag.stderr or ''))
            except Exception:
                pass
        return []

    out_dir = work_dir / "Endogenous_Results" / preset_key
    prog = None
    console = None
    def _cb(i, n, name, stat):
        if not headless:
            if prog is None:
                pass
            else:
                prog.progress(i / n, text=f"{i}/{n} {name} — {stat}")
                console.write(f"{i}/{n}  {name}: {stat}")

    rows = run_vina_batch(
        vina_exe=vina_exe,
        receptor_file=rec_tz,
        ligand_files=ligs,
        out_dir=out_dir,
        center=tuple(map(float, cfg["center"])),
        size=tuple(map(float, cfg["size"])),
        scoring="ad4",
        base_exhaustiveness=int(base_exhaustiveness),
        base_num_modes=int(base_num_modes),
        timeout_mode="no_timeout",
        timeout_s=int(timeout_s),
        max_retries=int(max_retries),
        exhu_backoff=float(exhu_backoff),
        modes_backoff=float(modes_backoff),
        progress_cb=_cb,
        maps_prefix=maps_prefix,
        skip_if_output_exists=bool(skip_exists),
    )
    return rows

# ==============================
# CLI mode (run presets without GUI)
# ==============================

def _run_cli():
    parser = argparse.ArgumentParser(description="ZincDock CLI (uses GUI code paths)")
    parser.add_argument("--cli", action="store_true", help="Run in CLI mode and skip Streamlit UI")
    parser.add_argument("--preset", type=str, default=None,
                        help="One of: CA, SOD1, HDAC1, HDAC2, HDAC3, HDAC4, HDAC6, or HDAC_ALL")
    parser.add_argument("--work-dir", type=str, default=None, help="Working directory for outputs and maps")
    parser.add_argument("--spacing", type=float, default=0.375, help="AD4 grid spacing (Å)")
    parser.add_argument("--exhaustiveness", type=int, default=64, help="Base exhaustiveness")
    parser.add_argument("--num-modes", type=int, default=10, help="Base num_modes")
    parser.add_argument("--timeout", type=int, default=300, help="Per-ligand soft timeout (s)")
    parser.add_argument("--retries", type=int, default=2, help="Max retries on failure")
    parser.add_argument("--skip-exists", action="store_true", help="Skip ligands with existing outputs")
    args, _ = parser.parse_known_args()

    if not args.cli:
        return False

    # Configure globals to mirror GUI defaults
    global work_dir, vina_exe, autogrid_exe, base_params, extra_params
    global spacing, base_exhaustiveness, base_num_modes, timeout_s, max_retries
    global exhu_backoff, modes_backoff, skip_exists, normalize_OA

    work_dir = Path(args.work_dir).expanduser().resolve() if args.work_dir else Path.cwd()
    work_dir.mkdir(parents=True, exist_ok=True)

    files_gui = (Path(__file__).resolve().parent / "Files_for_GUI")
    
    # Detect OS and use appropriate executable names
    is_windows_cli = platform.system() == "Windows"
    if is_windows_cli:
        vina_exe = (files_gui / "vina.exe").resolve()
        autogrid_exe = (files_gui / "autogrid4.exe").resolve()
    else:
        # Linux: try without .exe first
        vina_exe = (files_gui / "vina").resolve() if (files_gui / "vina").exists() else (files_gui / "vina.exe").resolve()
        autogrid_exe = (files_gui / "autogrid4").resolve() if (files_gui / "autogrid4").exists() else (files_gui / "autogrid4.exe").resolve()
    
    base_params = (files_gui / "AD4_parameters.dat").resolve()
    extra_params = (files_gui / "AD4Zn.dat").resolve() if (files_gui / "AD4Zn.dat").exists() else None

    spacing = float(args.spacing)
    base_exhaustiveness = int(args.exhaustiveness)
    base_num_modes = int(args.num_modes)
    timeout_s = int(args.timeout)
    max_retries = int(args.retries)
    exhu_backoff = 1.5
    modes_backoff = 1.25
    skip_exists = bool(args.skip_exists)
    normalize_OA = True

    if not vina_exe.exists():
        print(f"ERROR: vina.exe not found at {vina_exe}")
        sys.exit(1)
    if not autogrid_exe.exists():
        print(f"ERROR: autogrid4.exe not found at {autogrid_exe}")
        sys.exit(1)
    if not base_params.exists():
        print(f"ERROR: AD4_parameters.dat not found at {base_params}")
        sys.exit(1)

    targets: List[str]
    if args.preset is None:
        print("ERROR: --preset is required in --cli mode")
        sys.exit(2)
    if args.preset.upper() == "HDAC_ALL":
        targets = ["HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC6"]
    else:
        targets = [args.preset.upper()]

    for t in targets:
        print(f"\n=== Running preset: {t} ===")
        rows = run_endogenous_preset_ad4(t, headless=True)
        out_dir = work_dir / "Endogenous_Results" / t
        out_dir.mkdir(parents=True, exist_ok=True)
        if rows:
            csv_path = out_dir / f"{t}_results.csv"
            csv_path.write_bytes(results_to_csv_bytes(rows))
            print(f"Wrote: {csv_path}")
            succ = [r for r in rows if r.get("Status") == "Success"]
            print(f"Success: {len(succ)}/{len(rows)}")
        else:
            print("No results (see logs).")

    sys.exit(0)

if __name__ == "__main__":
    _run_cli()

# ==============================
# Setup: Ensure Linux executables have execute permissions
# ==============================
import stat
_setup_script = Path(__file__).parent / "setup_executables.py"
if _setup_script.exists():
    try:
        exec(open(_setup_script).read())
    except Exception:
        pass  # Silently fail if setup script has issues

# Also directly fix permissions on startup - more aggressive approach
_files_gui_setup = Path(__file__).parent / "Files_for_GUI"
if _files_gui_setup.exists():
    for exe_name in ["vina", "autogrid4", "autodock4"]:
        exe_path = _files_gui_setup / exe_name
        if exe_path.exists() and not exe_path.is_dir():
            try:
                # Try multiple methods to ensure execute permissions
                # Method 1: Add execute bits
                current_stat = exe_path.stat()
                new_mode = current_stat.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
                os.chmod(exe_path, new_mode)
                
                # Method 2: Set explicit 755 permissions
                try:
                    os.chmod(exe_path, 0o755)
                except:
                    pass
                
                # Method 3: Use subprocess to chmod if available
                try:
                    subprocess.run(["chmod", "+x", str(exe_path)], 
                                 capture_output=True, timeout=2, check=False)
                except:
                    pass
            except Exception:
                pass  # Silently fail if we can't set permissions

# ==============================
# Streamlit UI
# ==============================

st.set_page_config(
    page_title="ZincDock",
    layout="wide",
)

# Inject custom CSS to match blue theme
THEME_CSS = f"""
<style>
/* App background: light blue → dark blue gradient */
.stApp {{
    background: linear-gradient(135deg, {LIGHT_POWDER_BLUE} 0%, {SOFT_SKY_BLUE} 35%, {LIGHT_AZURE} 70%, {MIDNIGHT_AZURE} 100%);
}}

/* Main content container: white card on top of gradient */
.block-container {{
    background-color: rgba(255, 255, 255, 0.96);
    padding: 2rem 2rem 4rem 2rem;
    border-radius: 18px;
    box-shadow: 0 8px 24px rgba(0, 0, 0, 0.12);
}}

/* Force white text on all buttons in main content area */
.block-container .stButton > button,
.block-container button {{
    color: {WHITE} !important;
}}

.block-container .stButton > button *,
.block-container button * {{
    color: {WHITE} !important;
}}

/* Sidebar: solid color (no gradient) */
[data-testid="stSidebar"] {{
    background: {MIDNIGHT_AZURE};
    color: {WHITE};
}}

/* Sidebar text/icons stay light */
[data-testid="stSidebar"] * {{
    color: {WHITE} !important;
}}

/* Primary buttons (Run Docking, etc.) */
.stButton > button, .stDownloadButton > button {{
    background: {MEDIUM_STEEL_BLUE} !important;
    color: {WHITE} !important;
    border: none;
    border-radius: 999px;
    padding: 0.4rem 1.1rem;
    font-weight: 600;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
}}

.stButton > button:hover, .stDownloadButton > button:hover {{
    background: {DEEP_CERULEAN} !important;
    color: {WHITE} !important;
    transform: translateY(-1px);
    box-shadow: 0 6px 16px rgba(0, 0, 0, 0.18);
}}

/* Force white text on ALL Streamlit buttons - most aggressive selectors */
.stButton > button,
.stButton button,
button[data-baseweb="button"],
button.kind-primary,
button.kind-secondary,
.stButton > button[data-baseweb="button"],
.stButton > button.kind-primary,
.stButton > button.kind-secondary,
div[data-testid] button,
div[data-baseweb] button {{
    color: {WHITE} !important;
}}

/* Target button text/span elements inside buttons */
.stButton > button *,
.stButton button *,
button[data-baseweb="button"] *,
.stButton > button span,
.stButton > button p,
.stButton > button div,
button[data-baseweb="button"] span,
button[data-baseweb="button"] p,
button[data-baseweb="button"] div {{
    color: {WHITE} !important;
}}

/* Sidebar navigation buttons - distinct styling for better visibility */
[data-testid="stSidebar"] .stButton > button {{
    background: {LIGHT_AZURE} !important;
    color: {WHITE} !important;
    border: 2px solid {SOFT_SKY_BLUE} !important;
    border-radius: 8px !important;
    padding: 0.75rem 1rem !important;
    font-weight: 700 !important;
    font-size: 1rem !important;
    margin-bottom: 0.5rem !important;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.3), inset 0 1px 0 rgba(255, 255, 255, 0.1) !important;
    transition: all 0.2s ease !important;
}}

[data-testid="stSidebar"] .stButton > button * {{
    color: {WHITE} !important;
    font-weight: 700 !important;
    text-shadow: 0 1px 2px rgba(0, 0, 0, 0.3) !important;
}}

/* Selected/active sidebar button - darker and more prominent */
[data-testid="stSidebar"] .stButton > button.kind-primary,
[data-testid="stSidebar"] .stButton > button[aria-pressed="true"] {{
    background: {MEDIUM_STEEL_BLUE} !important;
    border: 2px solid {DEEP_CERULEAN} !important;
    box-shadow: 0 3px 12px rgba(0, 0, 0, 0.4), inset 0 1px 0 rgba(255, 255, 255, 0.15) !important;
}}

[data-testid="stSidebar"] .stButton > button:hover {{
    background: {SOFT_SKY_BLUE} !important;
    border-color: {LIGHT_AZURE} !important;
    color: {WHITE} !important;
    transform: translateX(4px) !important;
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.35), inset 0 1px 0 rgba(255, 255, 255, 0.15) !important;
}}

[data-testid="stSidebar"] .stButton > button:hover * {{
    color: {WHITE} !important;
    text-shadow: 0 1px 3px rgba(0, 0, 0, 0.4) !important;
}}

/* Final catch-all for any button with blue background - target by background color */
button[style*="background"],
button[style*="rgb(44, 124, 181)"],
button[style*="#2C7CB5"],
button[style*="#2c7cb5"] {{
    color: {WHITE} !important;
}}

button[style*="background"] *,
button[style*="rgb(44, 124, 181)"] *,
button[style*="#2C7CB5"] *,
button[style*="#2c7cb5"] * {{
    color: {WHITE} !important;
}}

/* File uploader components - revert text to black */
[data-testid="stFileUploader"],
[data-testid="stFileUploader"] *,
[data-testid="stFileUploader"] label,
[data-testid="stFileUploader"] span,
[data-testid="stFileUploader"] p,
[data-testid="stFileUploader"] div,
[data-testid="stFileUploader"] small,
.stFileUploader,
.stFileUploader *,
.stFileUploader label,
.stFileUploader span,
.stFileUploader p,
.stFileUploader div,
.stFileUploader small,
div[data-testid="stFileUploader"] label,
div[data-testid="stFileUploader"] span,
div[data-testid="stFileUploader"] p,
div[data-testid="stFileUploader"] small {{
    color: #000000 !important;
}}

/* File uploader button text should also be black - most specific selectors */
[data-testid="stFileUploader"] button,
[data-testid="stFileUploader"] button *,
[data-testid="stFileUploader"] button span,
[data-testid="stFileUploader"] button p,
[data-testid="stFileUploader"] button div,
[data-testid="stFileUploader"] .stButton button,
[data-testid="stFileUploader"] .stButton button *,
[data-testid="stFileUploader"] .stButton button span,
[data-testid="stFileUploader"] [data-baseweb="button"],
[data-testid="stFileUploader"] [data-baseweb="button"] *,
[data-testid="stFileUploader"] [data-baseweb="button"] span,
.stFileUploader button,
.stFileUploader button *,
.stFileUploader button span,
.stFileUploader .stButton button,
.stFileUploader .stButton button *,
.stFileUploader .stButton button span,
.stFileUploader [data-baseweb="button"],
.stFileUploader [data-baseweb="button"] *,
.stFileUploader [data-baseweb="button"] span {{
    color: #000000 !important;
}}

/* Target buttons inside file uploader containers with even higher specificity */
div[data-testid="stFileUploader"] > div button,
div[data-testid="stFileUploader"] > div button *,
div[data-testid="stFileUploader"] > div button span,
div[data-testid="stFileUploader"] div[data-baseweb="button"],
div[data-testid="stFileUploader"] div[data-baseweb="button"] *,
div[data-testid="stFileUploader"] div[data-baseweb="button"] span {{
    color: #000000 !important;
}}

/* Target the drag and drop area text specifically */
[data-testid="stFileUploader"] [data-baseweb="file-uploader"],
[data-testid="stFileUploader"] [data-baseweb="file-uploader"] *,
[data-testid="stFileUploader"] [data-baseweb="file-uploader"] span,
[data-testid="stFileUploader"] [data-baseweb="file-uploader"] p {{
    color: #000000 !important;
}}

/* Final override for ALL buttons within file uploader - must come after general button rules */
[data-testid="stFileUploader"] button[data-baseweb="button"],
[data-testid="stFileUploader"] button[data-baseweb="button"] *,
[data-testid="stFileUploader"] button[data-baseweb="button"] span,
[data-testid="stFileUploader"] button[data-baseweb="button"] p,
[data-testid="stFileUploader"] button[data-baseweb="button"] div,
.stFileUploader button[data-baseweb="button"],
.stFileUploader button[data-baseweb="button"] *,
.stFileUploader button[data-baseweb="button"] span {{
    color: #000000 !important;
    background-color: transparent !important;
}}

/* Radio + checkbox accent color */
input[type="radio"], input[type="checkbox"] {{
    accent-color: {MEDIUM_STEEL_BLUE};
}}

/* Progress bar to match palette */
[data-testid="stProgressBar"] > div > div {{
    background: linear-gradient(90deg, {LIGHT_POWDER_BLUE}, {SOFT_SKY_BLUE}, {LIGHT_AZURE});
}}

/* Tabs on a soft blue bar */
.stTabs [data-baseweb="tab-list"] {{
    background-color: rgba(153, 192, 222, 0.9);
    border-radius: 999px;
    padding: 0.25rem;
}}

.stTabs [data-baseweb="tab"] {{
    color: {MEDIUM_STEEL_BLUE};
    border-radius: 999px;
}}

.stTabs [data-baseweb="tab"][aria-selected="true"] {{
    background-color: {WHITE};
    color: {RICH_TEAL_BLUE};
}}

/* Headers in deep blue */
h1, h2, h3, h4, h5 {{
    color: {RICH_TEAL_BLUE};
}}

/* Metric labels, captions, small text slightly muted blue */
span, p, label {{
    color: {MIDNIGHT_AZURE};
}}

/* FINAL OVERRIDE: File uploader buttons MUST have black text - highest specificity */
[data-testid="stFileUploader"] .stButton > button,
[data-testid="stFileUploader"] .stButton > button *,
[data-testid="stFileUploader"] .stButton > button span,
[data-testid="stFileUploader"] .stButton > button p,
[data-testid="stFileUploader"] .stButton > button div,
[data-testid="stFileUploader"] button.stButton,
[data-testid="stFileUploader"] button.stButton *,
[data-testid="stFileUploader"] button.stButton span,
div[data-testid="stFileUploader"] .stButton > button,
div[data-testid="stFileUploader"] .stButton > button *,
div[data-testid="stFileUploader"] .stButton > button span,
div[data-testid="stFileUploader"] button,
div[data-testid="stFileUploader"] button *,
div[data-testid="stFileUploader"] button span {{
    color: #000000 !important;
}}
</style>
"""
st.markdown(THEME_CSS, unsafe_allow_html=True)

if "nav_open" not in st.session_state:
    st.session_state.nav_open = True
if "current_page" not in st.session_state:
    st.session_state.current_page = "ZincDock Demo"

nav_pages = [
    "Home",
    "Documentation",
    "ZincDock Demo",
    "Standard AutoDock",
    "Metalloprotein Docking",
    "GNINA ML Docking",
]

if st.session_state.current_page not in nav_pages:
    st.session_state.current_page = "ZincDock Demo"

with st.sidebar:
    toggle_label = "«" if st.session_state.nav_open else "»"
    if st.button(toggle_label):
        st.session_state.nav_open = not st.session_state.nav_open
        st.rerun()
    if st.session_state.nav_open:
        st.markdown("### Navigation")
        for nav_page in nav_pages:
            is_selected = nav_page == st.session_state.current_page
            button_label = f"**{nav_page}**" if is_selected else nav_page
            button_type = "primary" if is_selected else "secondary"
            if st.button(button_label, key=f"nav_{nav_page}", use_container_width=True, type=button_type):
                st.session_state.current_page = nav_page
                st.rerun()
    else:
        st.write("Navigation hidden")

page = st.session_state.current_page

if page == "Home":
    render_home_page()
    st.stop()

if page == "Documentation":
    render_documentation_page()
    st.stop()

page_mode = {
    "ZincDock Demo": "ad4",
    "Standard AutoDock": "vina",
    "Metalloprotein Docking": "ad4",
    "GNINA ML Docking": "gnina",
}.get(page, "generic")

state_prefix = "demo" if page == "ZincDock Demo" else page_mode

# Session state initialisation for docking workflow
if "docking_task" not in st.session_state:
    st.session_state.docking_task = None
if "stop_requested" not in st.session_state:
    st.session_state.stop_requested = False
if "docking_results" not in st.session_state:
    st.session_state.docking_results = None
if "docking_results_backend" not in st.session_state:
    st.session_state.docking_results_backend = None
if "docking_results_ad4" not in st.session_state:
    st.session_state.docking_results_ad4 = []
if "docking_results_out_dir" not in st.session_state:
    st.session_state.docking_results_out_dir = None
if "docking_status_message" not in st.session_state:
    st.session_state.docking_status_message = None

st.title(page)

# Working directory chooser
work_dir_input = st.text_input(
    "Working directory",
    value=str(Path.cwd()),
    help="All folders (prepared_ligands, ad4_maps, outputs) will be created here."
)
work_dir = Path(work_dir_input).expanduser().resolve()
work_dir.mkdir(parents=True, exist_ok=True)
st.caption(f"Using working directory: `{work_dir}`")

# Demo preset controls
_demo_default_center: Optional[Tuple[float, float, float]] = None
_demo_default_size: Optional[Tuple[float, float, float]] = None
_demo_default_spacing: Optional[float] = None
demo_selected_label: Optional[str] = None

if page == "ZincDock Demo":
    st.subheader("Choose Demo Receptor Preset")
    preset_cols = st.columns(len(DEMO_PRESETS))
    for idx, (label, settings) in enumerate(DEMO_PRESETS.items()):
        if preset_cols[idx].button(label, key=f"demo_preset_{label}"):
            st.session_state["demo_selected_preset"] = label

    demo_selected_label = st.session_state.get("demo_selected_preset", next(iter(DEMO_PRESETS)))
    demo_preset_info = DEMO_PRESETS[demo_selected_label]
    _demo_default_center = demo_preset_info["center"]
    _demo_default_size = demo_preset_info["size"]
    _demo_default_spacing = 0.375
else:
    demo_selected_label = None
    _demo_default_center = None
    _demo_default_size = None
    _demo_default_spacing = None

# Receptor and Ligand Setup
receptor_default_path = str((work_dir / "receptor.pdbqt").resolve())

st.subheader("Upload Receptor & Ligands")
upload_col1, upload_col2 = st.columns(2)

ligand_paths: List[Path] = []
receptor_path: Optional[Path] = None

with upload_col1:
    st.markdown("**Receptor**")
    # GNINA accepts PDB/PDBQT, others need PDBQT
    receptor_types = ["pdb", "pdbqt"] if page_mode == "gnina" else ["pdbqt"]
    receptor_label = "Upload receptor (PDB/PDBQT)" if page_mode == "gnina" else "Upload receptor (PDBQT)"
    receptor_upload = st.file_uploader(
        receptor_label,
        type=receptor_types,
        accept_multiple_files=False,
        key=f"{state_prefix}_receptor_upload"
    )
    receptor_store_dir = work_dir / f"{state_prefix}_receptor"
    receptor_store_dir.mkdir(parents=True, exist_ok=True)

    if receptor_upload is not None:
        saved_receptor = _save_uploaded_file(receptor_upload, receptor_store_dir)
        st.session_state[f"{state_prefix}_receptor_path"] = str(saved_receptor)
        receptor_path = saved_receptor

    stored_receptor_path = st.session_state.get(f"{state_prefix}_receptor_path")
    if receptor_path is None and stored_receptor_path:
        candidate = Path(stored_receptor_path)
        if candidate.exists():
            receptor_path = candidate

    if receptor_path:
        st.caption(f"Using receptor: {receptor_path.name}")
    else:
        st.info("Upload a receptor PDBQT to enable docking.")

with upload_col2:
    st.markdown("**Ligands**")
    # GNINA accepts MOL2/SDF/PDBQT, others need PDBQT
    ligand_types = ["mol2", "sdf", "pdbqt", "pdb"] if page_mode == "gnina" else ["pdbqt"]
    ligand_label = "Upload ligand files (MOL2/SDF/PDBQT/PDB)" if page_mode == "gnina" else "Upload ligand PDBQT files"
    ligand_uploads = st.file_uploader(
        ligand_label,
        type=ligand_types,
        accept_multiple_files=True,
        key=f"{state_prefix}_ligand_upload"
    )
    ligand_store_dir = work_dir / f"{state_prefix}_ligands"
    ligand_store_dir.mkdir(parents=True, exist_ok=True)

    stored_ligand_paths = st.session_state.get(f"{state_prefix}_ligand_paths", [])
    updated_paths: List[str] = list(stored_ligand_paths)

    if ligand_uploads:
        for uploaded in ligand_uploads:
            saved_lig = _save_uploaded_file(uploaded, ligand_store_dir)
            saved_str = str(saved_lig)
            if saved_str not in updated_paths:
                updated_paths.append(saved_str)
        st.session_state[f"{state_prefix}_ligand_paths"] = updated_paths

    for path_str in updated_paths:
        p = Path(path_str)
        if p.exists():
            ligand_paths.append(p)

    if ligand_paths:
        preview = ", ".join(p.name for p in ligand_paths[:5])
        if len(ligand_paths) > 5:
            preview += ", …"
        st.caption(f"{len(ligand_paths)} ligand(s) ready: {preview}")
    else:
        st.info("Upload one or more ligand PDBQT files to continue.")

# ---------------------------------------------

if page_mode == "gnina":
    # GNINA ML Docking page - uses SMINA backend
    allowed_backends = ["GNINA (ML)"]
    default_backend_label = "GNINA (ML)"
    grid_defaults = {
        "center": (0.0, 0.0, 0.0),
        "size": (20.0, 20.0, 20.0),
        "spacing": 0.375,
    }
    maps_prefix_default = str((work_dir / "ad4_maps" / "receptor_maps").resolve())
elif page_mode == "vina":
    allowed_backends = ["Vina (box)"]
    default_backend_label = "Vina (box)"
    grid_defaults = {
        "center": (0.0, 0.0, 0.0),
        "size": (0.0, 0.0, 0.0),
        "spacing": 0.375,
    }
    maps_prefix_default = str((work_dir / "ad4_maps" / "receptor_maps").resolve())
elif page == "ZincDock Demo":
    allowed_backends = ["AD4 (maps)"]
    default_backend_label = "AD4 (maps)"
    grid_defaults = {
        "center": _demo_default_center,
        "size": _demo_default_size,
        "spacing": _demo_default_spacing,
    }
    preset_slug = demo_selected_label.replace(" ", "_").lower()
    maps_prefix_default = str((work_dir / "ad4_maps" / preset_slug).resolve())
elif page_mode == "ad4":
    allowed_backends = ["AD4 (maps)"]
    default_backend_label = "AD4 (maps)"
    grid_defaults = {
        "center": (0.0, 0.0, 0.0),
        "size": (0.0, 0.0, 0.0),
        "spacing": 0.375,
    }
    maps_prefix_default = str((work_dir / "ad4_maps" / "receptor_maps").resolve())
else:
    allowed_backends = ["Vina (box)", "AD4 (maps)"]
    default_backend_label = "AD4 (maps)"
    grid_defaults = {
        "center": (-6.421, 0.342, 17.256),
        "size": (10.0, 10.0, 10.0),
        "spacing": 0.38,
    }
    maps_prefix_default = str((work_dir / "ad4_maps" / "receptor_maps").resolve())

saved_prefix = st.session_state.get(f"{state_prefix}_maps_prefix_saved")
if saved_prefix:
    st.session_state[f"{state_prefix}_maps_prefix"] = saved_prefix

build_maps_btn = False

with st.expander("Configuration", expanded=True):
    # For Standard AutoDock, only show Grid Box Settings
    if page == "Standard AutoDock":
        c2 = st.columns(1)[0]
        # Set backend and autodetect defaults for Standard AutoDock
        backend = allowed_backends[0] if len(allowed_backends) == 1 else default_backend_label
        autodetect = True
    else:
        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Executables & Scripts")
            if len(allowed_backends) == 1:
                backend = allowed_backends[0]
                st.markdown(f"**Docking backend:** `{backend}`")
            else:
                backend = st.radio(
                    "Docking backend",
                    allowed_backends,
                    index=allowed_backends.index(default_backend_label),
                    key=f"{state_prefix}_backend"
                )
            autodetect = False
            if backend == "Vina (box)":
                autodetect = st.checkbox(
                    "Auto-detect metal center (for Vina run)",
                    value=True,
                    key=f"{state_prefix}_autodetect"
                )
    with c2:
        st.subheader("Grid Box Settings")
        center_keys = {
            "x": f"{state_prefix}_center_x",
            "y": f"{state_prefix}_center_y",
            "z": f"{state_prefix}_center_z",
        }
        size_keys = {
            "x": f"{state_prefix}_size_x",
            "y": f"{state_prefix}_size_y",
            "z": f"{state_prefix}_size_z",
        }
        spacing_key = f"{state_prefix}_spacing"

        default_center = grid_defaults["center"]
        default_size = grid_defaults["size"]
        default_spacing = grid_defaults["spacing"]

        if page == "ZincDock Demo":
            for idx, axis in enumerate(["x", "y", "z"]):
                center_key = center_keys[axis]
                size_key = size_keys[axis]
                st.session_state[center_key] = default_center[idx]
                st.session_state[size_key] = default_size[idx]
            st.session_state[spacing_key] = default_spacing
            st.session_state[f"{state_prefix}_maps_prefix"] = maps_prefix_default
        else:
            for idx, axis in enumerate(["x", "y", "z"]):
                center_key = center_keys[axis]
                size_key = size_keys[axis]
                if center_key not in st.session_state:
                    st.session_state[center_key] = default_center[idx]
                if size_key not in st.session_state:
                    st.session_state[size_key] = default_size[idx]
            if spacing_key not in st.session_state:
                st.session_state[spacing_key] = default_spacing
            if f"{state_prefix}_maps_prefix" not in st.session_state:
                st.session_state[f"{state_prefix}_maps_prefix"] = maps_prefix_default

        grid_disabled = page == "ZincDock Demo"

        grid_c1, grid_c2, grid_c3 = st.columns(3)
        with grid_c1:
            center_x = st.number_input(
                "center_x",
                value=st.session_state[center_keys["x"]],
                format="%.3f",
                key=center_keys["x"],
                disabled=grid_disabled
            )
        with grid_c2:
            center_y = st.number_input(
                "center_y",
                value=st.session_state[center_keys["y"]],
                format="%.3f",
                key=center_keys["y"],
                disabled=grid_disabled
            )
        with grid_c3:
            center_z = st.number_input(
                "center_z",
                value=st.session_state[center_keys["z"]],
                format="%.3f",
                key=center_keys["z"],
                disabled=grid_disabled
            )

        sz1, sz2, sz3 = st.columns(3)
        with sz1:
            size_x = st.number_input(
                "size_x (Å)",
                value=st.session_state[size_keys["x"]],
                min_value=0.0,
                step=0.25,
                key=size_keys["x"],
                disabled=grid_disabled
            )
        with sz2:
            size_y = st.number_input(
                "size_y (Å)",
                value=st.session_state[size_keys["y"]],
                min_value=0.0,
                step=0.25,
                key=size_keys["y"],
                disabled=grid_disabled
            )
        with sz3:
            size_z = st.number_input(
                "size_z (Å)",
                value=st.session_state[size_keys["z"]],
                min_value=0.0,
                step=0.25,
                key=size_keys["z"],
                disabled=grid_disabled
            )

        spacing = st.number_input(
            "AD4 grid spacing (Å)",
            value=st.session_state[spacing_key],
            min_value=0.0,
            max_value=1.0,
            step=0.01,
            key=spacing_key,
            disabled=grid_disabled
        )

        if "AD4 (maps)" in allowed_backends:
            maps_prefix_input = st.text_input(
                "AD4 maps prefix (no extension)",
                value=st.session_state[f"{state_prefix}_maps_prefix"],
                help="Folder will be created if missing (receptor_maps.gpf, *.map, *.fld, etc.).",
                key=f"{state_prefix}_maps_prefix"
            )
            force_extra_types = st.text_input(
                "Force-include extra ligand atom types when building/patching maps (comma-separated)",
                value=st.session_state.get(f"{state_prefix}_force_types", "S,NA"),
                help="If you *know* you need maps like S or NA, list them here to guarantee creation.",
                key=f"{state_prefix}_force_types"
            )
            build_maps_btn = st.button(
                "Build/Update AD4 Maps",
                key=f"{state_prefix}_build_maps"
            )
        else:
            maps_prefix_input = st.session_state.get(f"{state_prefix}_maps_prefix", maps_prefix_default)
            force_extra_types = "S,NA"
            build_maps_btn = False

files_gui_dir = work_dir / "Files_for_GUI"
is_windows = platform.system() == "Windows"

# ==============================
# Endogenous Docking Presets (AD4 maps)
# ==============================

 

st.subheader("Docking Parameters")
p1, p2, p3, p4 = st.columns(4)
params_disabled = page == "ZincDock Demo"
with p1:
    if page_mode == "gnina":
        scoring = "GNINA (ML)"
        st.markdown(f"**Scoring function:** `{scoring}`")
    else:
        scoring = "ad4" if backend == "AD4 (maps)" else "vina"
        st.markdown(f"**Scoring function:** `{scoring}`")
with p2:
    base_exhaustiveness = st.number_input(
        "Base exhaustiveness",
        value=DEMO_PARAM_DEFAULTS["base_exhaustiveness"],
        min_value=1,
        step=1,
        key=f"{state_prefix}_base_exhaustiveness",
        disabled=params_disabled
    )
with p3:
    base_num_modes = st.number_input(
        "Base num_modes",
        value=DEMO_PARAM_DEFAULTS["base_num_modes"],
        min_value=1,
        step=1,
        key=f"{state_prefix}_base_num_modes",
        disabled=params_disabled
    )
with p4:
    out_dir_name = st.text_input(
        "Output folder name",
        value=DEMO_PARAM_DEFAULTS["output_name"],
        key=f"{state_prefix}_output_name",
        disabled=params_disabled
    )

timeout_options = ["No timeout (recommended)", "Soft timeout with retries"]
t_default_index = timeout_options.index(DEMO_PARAM_DEFAULTS["timeout_mode"])

t1, t2, t3, t4 = st.columns(4)
with t1:
    timeout_mode = st.selectbox(
        "Timeout mode",
        timeout_options,
        index=t_default_index,
        key=f"{state_prefix}_timeout_mode",
        disabled=params_disabled
    )
with t2:
    timeout_s = st.number_input(
        "Per-ligand timeout (s) if using soft timeout",
        value=DEMO_PARAM_DEFAULTS["timeout_s"],
        min_value=30,
        step=10,
        key=f"{state_prefix}_timeout_s",
        disabled=params_disabled
    )
with t3:
    max_retries = st.number_input(
        "Max retries on failure",
        value=DEMO_PARAM_DEFAULTS["max_retries"],
        min_value=0,
        step=1,
        key=f"{state_prefix}_max_retries",
        disabled=params_disabled
    )
with t4:
    skip_exists = st.checkbox(
        "Skip ligands with existing outputs",
        value=DEMO_PARAM_DEFAULTS["skip_existing"],
        key=f"{state_prefix}_skip_exists",
        disabled=params_disabled
    )

b1, b2 = st.columns(2)
with b1:
    exhu_backoff = st.number_input(
        "Exhaustiveness multiplier on retry",
        value=DEMO_PARAM_DEFAULTS["exhaustiveness_retry"],
        min_value=1.0,
        step=0.1,
        key=f"{state_prefix}_exhaustiveness_retry",
        disabled=params_disabled
    )
with b2:
    modes_backoff = st.number_input(
        "num_modes multiplier on retry",
        value=DEMO_PARAM_DEFAULTS["num_modes_retry"],
        min_value=1.0,
        step=0.05,
        key=f"{state_prefix}_num_modes_retry",
        disabled=params_disabled
    )

# Detect operating system
exe_ext = ".exe" if is_windows else ""

# Always enable receptor oxygen normalization (O→OA)
normalize_OA = True

# Auto-detect paths (fallback to defaults if Files_for_GUI doesn't exist)
# Try Windows .exe first, then Linux executable (no extension)
vina_exe = None
if is_windows:
    vina_exe = (files_gui_dir / "vina.exe").resolve() if (files_gui_dir / "vina.exe").exists() else None
else:
    # Linux: try without .exe extension
    vina_exe = (files_gui_dir / "vina").resolve() if (files_gui_dir / "vina").exists() else None
    if not vina_exe:
        # Fallback: try vina.exe in case it's there
        vina_exe = (files_gui_dir / "vina.exe").resolve() if (files_gui_dir / "vina.exe").exists() else None

autogrid_exe = None
if is_windows:
    autogrid_exe = (files_gui_dir / "autogrid4.exe").resolve() if (files_gui_dir / "autogrid4.exe").exists() else None
else:
    # Linux: try without .exe extension
    autogrid_exe = (files_gui_dir / "autogrid4").resolve() if (files_gui_dir / "autogrid4").exists() else None
    if not autogrid_exe:
        # Fallback: try autogrid4.exe in case it's there
        autogrid_exe = (files_gui_dir / "autogrid4.exe").resolve() if (files_gui_dir / "autogrid4.exe").exists() else None

python_exe = Path(sys.executable)
zinc_pseudo_py = (files_gui_dir / "zinc_pseudo.py").resolve() if (files_gui_dir / "zinc_pseudo.py").exists() else None
base_params = (files_gui_dir / "AD4_parameters.dat").resolve() if (files_gui_dir / "AD4_parameters.dat").exists() else None
extra_params = (files_gui_dir / "AD4Zn.dat").resolve() if (files_gui_dir / "AD4Zn.dat").exists() else None

if build_maps_btn:
    if backend != "AD4 (maps)":
        st.warning("Switch the docking backend to 'AD4 (maps)' before building affinity maps.")
    else:
        try:
            maps_prefix_str = (maps_prefix_input or "").strip()
            if not maps_prefix_str:
                raise ValueError("Enter a maps prefix path before building AutoGrid maps.")
            maps_prefix_path = Path(maps_prefix_str).expanduser()
            if not maps_prefix_path.is_absolute():
                maps_prefix_path = (work_dir / maps_prefix_path).resolve()
            else:
                maps_prefix_path = maps_prefix_path.resolve()

            if not maps_prefix_path.parts:
                maps_prefix_path = (work_dir / "ad4_maps" / "receptor_maps").resolve()

            force_types = set()
            if force_extra_types:
                force_types = {tok.strip().upper() for tok in force_extra_types.split(",") if tok.strip()}

            with st.spinner("Building AutoGrid4 maps…"):
                if receptor_path is None or not receptor_path.exists():
                    raise FileNotFoundError("Upload a receptor PDBQT before building maps.")
                if autogrid_exe is None or not autogrid_exe.exists():
                    raise FileNotFoundError("AutoGrid4 executable not found. Set it in the Executables section before building maps.")
                if base_params is None or not base_params.exists():
                    raise FileNotFoundError("AD4_parameters.dat is missing. Configure it in the Executables section before building maps.")

                spacing_val = float(spacing)
                if spacing_val <= 0.0:
                    raise ValueError("AD4 grid spacing must be greater than 0.0 Å before building maps.")
                size_vals = (float(size_x), float(size_y), float(size_z))
                if any(val <= 0.0 for val in size_vals):
                    raise ValueError("Grid box dimensions must all be greater than 0.0 Å before building maps.")

                maps_prefix_clean = maps_prefix_path.with_suffix("") if maps_prefix_path.suffix else maps_prefix_path
                maps_dir = maps_prefix_clean.parent
                maps_dir.mkdir(parents=True, exist_ok=True)

                merged_params = maps_dir / "AD4_parameters_plus_ZnTZ.dat"
                merge_parameter_files(base_params, extra_params, merged_params)

                receptor_copy = maps_dir / receptor_path.name
                shutil.copy2(receptor_path, receptor_copy)
                if any(ch in receptor_copy.name for ch in " ()"):
                    sanitized_name = receptor_copy.name.replace(" ", "_").replace("(", "").replace(")", "")
                    sanitized_path = receptor_copy.with_name(sanitized_name)
                    receptor_copy.rename(sanitized_path)
                    receptor_copy = sanitized_path
                if normalize_OA:
                    try:
                        normalize_receptor_oxygen_to_OA(receptor_copy, receptor_copy)
                    except Exception:
                        pass

                invalid_types = {"K", "NA", "MG", "CA", "CL", "FE", "MN", "ZN", "CU", "CO", "NI"}
                rec_filtered = maps_dir / f"{receptor_copy.stem}_filtered.pdbqt"
                try:
                    with open(receptor_copy, "r", errors="ignore") as fin, open(rec_filtered, "w", encoding="utf-8") as fout:
                        for line in fin:
                            if line.startswith(("ATOM", "HETATM")):
                                toks = line.split()
                                if toks and toks[-1] in invalid_types:
                                    continue
                            fout.write(line)
                    if rec_filtered.stat().st_size == 0:
                        rec_filtered.unlink(missing_ok=True)
                        rec_filtered = receptor_copy
                except Exception:
                    rec_filtered = receptor_copy
                rec_tz = rec_filtered

                rec_types = [t for t in read_types_from_pdbqt(rec_tz) if t not in invalid_types]
                if not rec_types:
                    rec_types = [t for t in read_types_from_pdbqt(receptor_copy) if t not in invalid_types]
                if not rec_types:
                    raise ValueError("No valid receptor atom types detected after filtering. Check the receptor PDBQT file.")

                ligand_types = ligand_types_union(ligand_paths) if ligand_paths else set()
                if force_types:
                    ligand_types.update(force_types)
                ligand_types = {t.strip().upper() for t in ligand_types if t}
                if not ligand_types:
                    ligand_types = {"C", "F", "O", "OA", "S", "NA"}
                ligand_types_sorted = sorted(ligand_types)

                npts = tuple(max(10, int(round(dim / spacing_val))) for dim in size_vals)
                gpf_out = maps_prefix_clean.with_suffix(".gpf")
                write_simple_gpf(
                    gpf_path=gpf_out,
                    receptor_tz_filename=rec_tz.name,
                    maps_prefix_basename=maps_prefix_clean.name,
                    npts_xyz=npts,
                    spacing=spacing_val,
                    center_xyz=(float(center_x), float(center_y), float(center_z)),
                    receptor_types=rec_types,
                    ligand_types=ligand_types_sorted,
                    parameter_file_rel=merged_params.name,
                )

                proc = run_autogrid4(autogrid_exe, maps_dir, gpf_out)
                if proc.returncode != 0:
                    debug_msg = [
                        "AutoGrid4 stderr:\n" + (proc.stderr or ""),
                        "AutoGrid4 stdout:\n" + (proc.stdout or ""),
                    ]
                    glg_file = gpf_out.with_suffix(".glg")
                    if glg_file.exists():
                        debug_msg.append("Generated .glg log:\n" + glg_file.read_text(encoding="utf-8", errors="ignore"))
                    st.error("AutoGrid4 failed while building maps. See console output below.")
                    with st.expander("AutoGrid4 error output", expanded=True):
                        for section in debug_msg:
                            st.code(section or "(empty)", language="text")
                    failure_log = maps_dir / "autogrid4_last_failure.log"
                    try:
                        failure_log.write_text("\n\n".join(debug_msg), encoding="utf-8")
                    except Exception:
                        pass
                    raise RuntimeError("AutoGrid4 failed while building maps.")

                map_details = {
                    "returncode": proc.returncode,
                    "stdout": proc.stdout or "",
                    "stderr": proc.stderr or "",
                    "maps_dir": maps_dir,
                    "maps_prefix": maps_prefix_clean,
                    "gpf": gpf_out,
                    "atom_types": ligand_types_sorted,
                    "map_files": sorted(str(p) for p in maps_dir.glob(f"{maps_prefix_clean.name}.*.map")),
                }

            atom_types = ", ".join(map_details.get("atom_types", [])) or "(none)"
            st.success(
                f"AutoGrid4 completed successfully. Maps stored under `{map_details['maps_dir']}` "
                f"for atom types: {atom_types}."
            )
            st.caption(f"GPF file: `{map_details['gpf'].name}`")
            st.session_state[f"{state_prefix}_maps_prefix_saved"] = str(map_details["maps_prefix"])

            maps_available = sorted(list_maps_present(map_details["maps_prefix"]))
            if not maps_available:
                st.warning("AutoGrid4 completed but no affinity map files were detected. Check the log output for details.")

            stdout = map_details.get("stdout")
            stderr = map_details.get("stderr")
            if stdout or stderr:
                with st.expander("AutoGrid4 console output", expanded=False):
                    if stdout:
                        st.code(stdout, language="text")
                    if stderr:
                        st.code(stderr, language="text")
            glg_file = map_details["gpf"].with_suffix(".glg")
            if glg_file.exists():
                with st.expander("AutoGrid4 .glg log", expanded=False):
                    st.code(glg_file.read_text(encoding="utf-8", errors="ignore"), language="text")
        except (FileNotFoundError, PermissionError, ValueError, RuntimeError) as err:
            st.error(str(err))
        except Exception as err:
            st.error(f"Unexpected error while building maps: {err}")

# Platform and executable status (silent detection - no warnings)

# Test executables button
st.subheader("Tools")
run_btn = st.button("Run Docking", type="primary")

rows: List[dict] = []
ad4_rows: List[dict] = []

if run_btn:
    if not vina_exe.exists():
        st.error("Vina executable not found.")
        st.stop()
    if receptor_path is None or not receptor_path.exists():
        st.error("Receptor file missing/invalid.")
        st.stop()
    if not ligand_paths:
        st.error("No ligand files found. Upload ligand PDBQT files first.")
        st.stop()

    cx, cy, cz = float(center_x), float(center_y), float(center_z)
    if autodetect and backend == "Vina (box)":
        auto_center = autodetect_metal_center(receptor_path)
        if auto_center:
            cx, cy, cz = auto_center
            st.info(f"Auto-detected metal center: {cx:.3f}, {cy:.3f}, {cz:.3f}")

    out_dir = (work_dir / out_dir_name).resolve()

    maps_prefix = None
    if backend == "AD4 (maps)":
        maps_prefix = Path(maps_prefix_input).expanduser().resolve()
        required_types = sorted(ligand_types_union(ligand_paths) or {"C", "F", "OA", "S", "NA"})
        have = list_maps_present(maps_prefix)
        base_req = [
            maps_prefix.parent / f"{maps_prefix.name}.maps.fld",
            maps_prefix.parent / f"{maps_prefix.name}.e.map",
            maps_prefix.parent / f"{maps_prefix.name}.d.map",
        ]
        missing_files = [p for p in base_req if not p.exists()]
        missing_types = [t for t in required_types if t not in have]
        if missing_files or missing_types:
            if missing_files:
                st.error("AD4 fld/e/d maps missing:\n" + "\n".join(str(p) for p in missing_files))
            if missing_types:
                st.error("Affinity maps missing for types:\n" + ", ".join(missing_types))
            st.stop()

    tm_mode_key = "no_timeout" if timeout_mode.startswith("No timeout") else "soft_timeout"

    prog = st.progress(0, text="Starting docking…")
    console = st.empty()

    def _cb(i, n, name, stat):
        prog.progress(i / n, text=f"{i}/{n} {name} — {stat}")
        console.write(f"{i}/{n}  {name}: {stat}")

    with st.spinner("Running docking…"):
        rows = run_vina_batch(
            vina_exe=vina_exe,
            receptor_file=receptor_path,
            ligand_files=ligand_paths,
            out_dir=out_dir,
            center=(cx, cy, cz),
            size=(float(size_x), float(size_y), float(size_z)),
            scoring="ad4" if backend == "AD4 (maps)" else "vina",
            base_exhaustiveness=int(base_exhaustiveness),
            base_num_modes=int(base_num_modes),
            timeout_mode=tm_mode_key,
            timeout_s=int(timeout_s),
            max_retries=int(max_retries),
            exhu_backoff=float(exhu_backoff),
            modes_backoff=float(modes_backoff),
            progress_cb=_cb,
            maps_prefix=maps_prefix,
            skip_if_output_exists=bool(skip_exists),
        )
        if backend == "AD4 (maps)":
            ad4_rows = rows

if rows:
    drop_cols = [c for c in rows[0].keys() if c.startswith("AD4_")]
    df = pd.DataFrame(rows)
    st.success("Docking complete.")
    st.dataframe(df[[c for c in df.columns if c not in drop_cols]], use_container_width=True)

    st.subheader("Result Files")
    for idx, row in enumerate(rows):
        ligand_name = row.get("Ligand", f"Ligand {idx+1}")
        out_path = row.get("Output_File")
        log_path = row.get("Log_File")
        if out_path and Path(out_path).exists():
            archive = io.BytesIO()
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(out_path, arcname=Path(out_path).name)
                if log_path and Path(log_path).exists():
                    zf.write(log_path, arcname=Path(log_path).name)
            archive.seek(0)
            st.download_button(
                label=f"Download {ligand_name} results (ZIP)",
                data=archive.getvalue(),
                file_name=f"{Path(out_path).stem}.zip",
                key=f"dl_zip_{idx}"
            )
        else:
            st.caption(f"No PDBQT/log available for {ligand_name}")

    if backend == "AD4 (maps)" and ad4_rows:
        st.subheader("AD4 Summary")
        ad4_success = [r for r in ad4_rows if r.get("Status") == "Success"]
        st.write(f"AD4 successes: {len(ad4_success)}/{len(ad4_rows)} ligands")
        comp_df = pd.DataFrame(ad4_success)
        display_cols = [
            "Ligand",
            "AD4_Affinity", "AD4_Intermolecular", "AD4_Internal", "AD4_Torsional",
            "Binding_Affinity", "Num_Poses"
        ]
        available_cols = [c for c in display_cols if c in comp_df.columns]
        if available_cols:
            st.write("Energy Components:")
            st.dataframe(comp_df[available_cols], use_container_width=True)
    else:
        st.subheader("Vina Summary")
        try:
            vina_affs = [float(r["Binding_Affinity"]) for r in rows if r.get("Binding_Affinity") not in ("", "N/A", None)]
            if vina_affs:
                st.write(f"Binding affinities range: {min(vina_affs):.1f} to {max(vina_affs):.1f} kcal/mol")
                st.write(f"Average binding affinity: {sum(vina_affs)/len(vina_affs):.1f} kcal/mol")
        except Exception:
            pass

    csv_bytes = _cached_file_bytes(results_to_csv_bytes(rows))
    st.download_button(
        "Download results CSV",
        data=csv_bytes,
        file_name="metallodock_results.csv",
        mime="text/csv",
    )
    if out_dir.exists():
        all_zip = _cached_file_bytes(zip_outputs(out_dir))
        st.download_button(
            "Download all output PDBQTs (ZIP)",
            data=all_zip,
            file_name=f"{out_dir.name}.zip",
            mime="application/zip",
        )
else:
    st.info("Run docking to see results.")

def _process_docking_task() -> None:
    """Legacy no-op to maintain compatibility with older session state."""
    st.session_state.docking_task = None
    st.session_state.stop_requested = False
    return

def list_maps_present(maps_prefix: Path) -> Set[str]:
    """Return set of atom types that already have an affinity map file for this prefix."""
    present = set()
    folder = maps_prefix.parent
    base = maps_prefix.name
    for p in folder.glob(f"{base}.*.map"):
        # expecting base.<TYPE>.map
        t = p.suffixes[-2].lstrip(".") if len(p.suffixes) >= 2 else None
        if t:
            present.add(t)
    return present

FORCE_DEFAULT_TYPES = {"C", "F", "O", "OA", "S", "NA"}


def build_ad4_maps(
    receptor: Path,
    ligands: List[Path],
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    spacing: float,
    maps_prefix: Path,
    autogrid_exe: Path,
    base_params: Optional[Path],
    extra_params: Optional[Path],
    normalize_oxygen: bool = True,
    force_types: Optional[Set[str]] = None,
) -> dict:
    if receptor is None or not receptor.exists():
        raise FileNotFoundError("Upload a receptor PDBQT before building maps.")
    if autogrid_exe is None or not autogrid_exe.exists():
        raise FileNotFoundError("AutoGrid4 executable not found. Set it in the Executables section before building maps.")
    if base_params is None or not base_params.exists():
        raise FileNotFoundError("AD4_parameters.dat is missing. Configure it in the Executables section before building maps.")

    spacing_val = float(spacing)
    if spacing_val <= 0.0:
        raise ValueError("AD4 grid spacing must be greater than 0.0 Å before building maps.")
    size_vals = tuple(float(v) for v in size)
    if any(val <= 0.0 for val in size_vals):
        raise ValueError("Grid box dimensions must all be greater than 0.0 Å before building maps.")

    prefix_clean = maps_prefix.with_suffix("") if maps_prefix.suffix else maps_prefix
    maps_dir = prefix_clean.parent
    maps_dir.mkdir(parents=True, exist_ok=True)

    merged_params = maps_dir / "AD4_parameters_plus_ZnTZ.dat"
    merge_parameter_files(base_params, extra_params, merged_params)

    receptor_copy = maps_dir / receptor.name
    shutil.copy2(receptor, receptor_copy)
    if normalize_oxygen:
        try:
            normalize_receptor_oxygen_to_OA(receptor_copy, receptor_copy)
        except Exception:
            pass

    invalid_types = {"K", "Na", "MG", "CA", "CL", "FE", "MN", "ZN", "CU", "CO", "NI"}
    rec_filtered = maps_dir / (receptor_copy.stem + "_filtered.pdbqt")
    try:
        with open(receptor_copy, "r", errors="ignore") as fin, open(rec_filtered, "w", encoding="utf-8") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    toks = line.split()
                    if toks and toks[-1] in invalid_types:
                        continue
                fout.write(line)
        if rec_filtered.stat().st_size == 0:
            rec_filtered.unlink(missing_ok=True)
            rec_filtered = receptor_copy
    except Exception:
        rec_filtered = receptor_copy
    rec_tz = rec_filtered

    rec_types = [t for t in read_types_from_pdbqt(rec_tz) if t not in invalid_types]
    if not rec_types:
        rec_types = [t for t in read_types_from_pdbqt(receptor_copy) if t not in invalid_types]
    if not rec_types:
        raise ValueError("No valid receptor atom types detected after filtering. Check the receptor PDBQT file.")

    ligand_types = ligand_types_union(ligands) if ligands else set()
    if force_types:
        ligand_types.update(force_types)
    ligand_types = {t.strip().upper() for t in ligand_types if t}
    if not ligand_types:
        ligand_types = set(FORCE_DEFAULT_TYPES)
    if not ligand_types:
        raise ValueError("Unable to detect ligand atom types. Upload ligands or list atom types in the force-include box.")
    ligand_types_sorted = sorted(ligand_types)

    npts = tuple(max(10, int(round(dim / spacing_val))) for dim in size_vals)
    gpf_out = prefix_clean.with_suffix(".gpf")
    write_simple_gpf(
        gpf_path=gpf_out,
        receptor_tz_filename=rec_tz.name,
        maps_prefix_basename=prefix_clean.name,
        npts_xyz=npts,
        spacing=spacing_val,
        center_xyz=tuple(float(v) for v in center),
        receptor_types=rec_types,
        ligand_types=ligand_types_sorted,
        parameter_file_rel=merged_params.name,
    )

    proc = run_autogrid4(autogrid_exe, maps_dir, gpf_out)
    if proc.returncode != 0:
        raise RuntimeError("AutoGrid4 failed while building maps.")

    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout or "",
        "stderr": proc.stderr or "",
        "maps_dir": maps_dir,
        "maps_prefix": prefix_clean,
        "gpf": gpf_out,
        "atom_types": ligand_types_sorted,
        "map_files": sorted(str(p) for p in maps_dir.glob(f"{prefix_clean.name}.*.map")),
    }


def build_ad4_maps_for_selection(*args, **kwargs):
    """Backward-compatible wrapper for legacy code paths."""
    return build_ad4_maps(*args, **kwargs)

