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

import streamlit as st
import pandas as pd
import argparse
import sys


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
    autodock4_exe: Optional[Path] = None,  # For AD4 docking
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
) -> Tuple[bool, str, int, Optional[str]]:
    """Return (ok, affinity, nposes, missing_atom_type). Writes log regardless."""
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
                                  exhaustiveness, num_modes, seed, None, autodock4_exe)
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
                            return (True, aff, nposes, None)
            
            # Log original error details
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR DETAILS ----\n")
                lf.write(f"Return code: {proc.returncode}\n")
                lf.write(f"Missing map: {missing}\n")
                lf.write(f"Error: {error_msg[:500]}\n")
            return (False, "", 0, missing)

        if not out_pdbqt.exists():
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR: Output file not created ----\n")
                lf.write(f"Expected: {out_pdbqt}\n")
            return (False, "", 0, None)

        if out_pdbqt.stat().st_size == 0:
            with open(log_file, "a", encoding="utf-8") as lf:
                lf.write(f"\n---- ERROR: Output file is empty ----\n")
                lf.write(f"File: {out_pdbqt}\n")
            return (False, "", 0, None)

        aff = parse_binding_affinity(out_pdbqt)
        nposes = count_poses(out_pdbqt)
        
        # Log parsing results
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write(f"\n---- PARSING RESULTS ----\n")
            lf.write(f"Binding affinity: {aff}\n")
            lf.write(f"Number of poses: {nposes}\n")
        
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
            return (False, "", nposes, None)

        return (True, aff, nposes, None)

    except subprocess.TimeoutExpired as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- TIMEOUT ----\n")
            lf.write(str(e))
        return (False, "", 0, None)
    except Exception as e:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("\n---- EXCEPTION ----\n")
            lf.write(str(e))
        return (False, "", 0, None)

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
    autodock4_exe: Optional[Path] = None,  # For AD4 docking
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

        while tried <= max_retries and not ok:
            seed = random.randint(1, 2**31-1)
            if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Running (try {tried+1}/{max_retries+1})")
            ok, aff, nposes, last_missing = _run_one(
                vina_exe, mode, receptor_file, lig, out_pdbqt, log_file,
                center, size, ex, nm, seed, per_try_timeout, maps_prefix
            )

            if ok:
                if progress_cb: progress_cb(i, len(ligand_files), lig_name, f"Success | Score {aff} ({nposes} poses)")
                break

            if last_missing and progress_cb:
                progress_cb(i, len(ligand_files), lig_name, f"Failed - Missing map: {last_missing}")

            tried += 1
            ex = max(ex, int(math.ceil(ex * exhu_backoff)))
            nm = max(nm, int(math.ceil(nm * modes_backoff)))

        status = "Success" if ok else ("Failed - Timeout" if timeout_mode != "no_timeout" else f"Failed{(' - Missing '+last_missing) if last_missing else ''}")
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
                if last_missing:
                    progress_cb(i, len(ligand_files), lig_name, f"FAILED | Missing map: {last_missing}")
                else:
                    progress_cb(i, len(ligand_files), lig_name, "FAILED")

    return rows

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
            "size": (20.0, 20.0, 20.0),
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
        progress_cb=_cb if not headless else None,
        maps_prefix=maps_prefix,
        skip_if_output_exists=bool(skip_exists),
    )
    return rows

# ==============================
# CLI mode (run presets without GUI)
# ==============================

def _run_cli():
    parser = argparse.ArgumentParser(description="MetalloDock CLI (uses GUI code paths)")
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

st.set_page_config(page_title="MetalloDock", layout="wide")
st.title("MetalloDock")

# Working directory chooser
work_dir_input = st.text_input(
    "Working directory",
    value=str(Path.cwd()),
    help="All folders (prepared_ligands, ad4_maps, outputs) will be created here."
)
work_dir = Path(work_dir_input).expanduser().resolve()
work_dir.mkdir(parents=True, exist_ok=True)
st.caption(f"Using working directory: `{work_dir}`")

# Receptor and Ligand Upload Section (Top)
st.subheader("📤 Upload Receptor & Ligands")
upload_col1, upload_col2 = st.columns(2)

with upload_col1:
    st.markdown("**Receptor**")
    receptor_input_mode = st.radio("Provide receptor via:", ["Upload file", "Local path"], index=0, horizontal=True)
    receptor_uploaded = None
    receptor_local_path = None
    if receptor_input_mode == "Upload file":
        receptor_uploaded = st.file_uploader("Upload receptor (PDBQT)", type=["pdbqt"], accept_multiple_files=False, key="receptor_upload")
    else:
        receptor_local_path = st.text_input("Receptor file path", value=str((work_dir / "receptor.pdbqt").resolve()), key="receptor_path")

with upload_col2:
    st.markdown("**Ligands**")
    lig_src = st.text_input("Ligand SOURCE folder (to prepare)", value=str((work_dir / "Files_for_GUI" / "Ligands").resolve()), key="lig_src")
    prep_btn = st.button("🧰 Prepare ligands from SOURCE → prepared_ligands", key="prep_btn")
    lig_mode = st.radio("Docking ligands come from:", ["prepared_ligands folder", "Upload now"], index=0, key="lig_mode")
    ligand_uploads = []
    if lig_mode == "Upload now":
        ligand_uploads = st.file_uploader("Upload ligand PDBQT files", type=["pdbqt"], accept_multiple_files=True, key="ligand_upload")

with st.expander("⚙️ Configuration", expanded=True):
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("Executables & Scripts")
        st.info("🗺️ **AD4 mode:** Uses AutoGrid4 maps and AutoDock4 (recommended for metalloproteins).")
        autodetect = st.checkbox("Auto-detect metal center", value=True)

        # Auto-detect executables and parameters from Files_for_GUI (no user input needed)
        files_gui_dir = work_dir / "Files_for_GUI"
        import sys

    with c2:
        st.subheader("Grid Box Settings")
        grid_c1, grid_c2, grid_c3 = st.columns(3)
        with grid_c1:
            center_x = st.number_input("center_x", value=24.654, format="%.3f")
        with grid_c2:
            center_y = st.number_input("center_y", value=-0.568, format="%.3f")
        with grid_c3:
            center_z = st.number_input("center_z", value=-1.090, format="%.3f")

        sz1, sz2, sz3 = st.columns(3)
        with sz1:
            size_x = st.number_input("size_x (Å)", value=28.0, min_value=1.0, step=0.25)
        with sz2:
            size_y = st.number_input("size_y (Å)", value=30.0, min_value=1.0, step=0.25)
        with sz3:
            size_z = st.number_input("size_z (Å)", value=28.0, min_value=1.0, step=0.25)

        spacing = st.number_input("AD4 grid spacing (Å)", value=0.375, min_value=0.2, max_value=1.0, step=0.025)

        st.markdown("**Maps**")
        maps_prefix_input = st.text_input(
            "AD4 maps prefix (no extension)",
            value=str((work_dir / "ad4_maps" / "receptor_maps").resolve()),
            help="Folder will be created if missing (receptor_maps.gpf, *.map, *.fld, etc.)."
        )
        force_extra_types = st.text_input(
            "Force-include extra ligand atom types when building/patching maps (comma-separated)",
            value="S,NA",
            help="If you *know* you need maps like S or NA, list them here to guarantee creation."
        )
        build_maps_btn = st.button("🗺️ Build/Update AD4 maps (auto-detect & include missing types)")

# ==============================
# Endogenous Docking Presets (AD4 maps)
# ==============================

 

st.subheader("Docking Parameters")
p1, p2, p3, p4 = st.columns(4)
with p1:
    scoring = "ad4"
    st.markdown("**Scoring function:** `ad4`")
with p2:
    base_exhaustiveness = st.number_input("Base exhaustiveness", value=64, min_value=1, step=1)
with p3:
    base_num_modes = st.number_input("Base num_modes", value=10, min_value=1, step=1)
with p4:
    out_dir_name = st.text_input("Output folder name", value="PFAS_Docking_Results")

t1, t2, t3, t4 = st.columns(4)
with t1:
    timeout_mode = st.selectbox("Timeout mode", ["No timeout (recommended)", "Soft timeout with retries"], index=0)
with t2:
    timeout_s = st.number_input("Per-ligand timeout (s) if using soft timeout", value=300, min_value=30, step=10)
with t3:
    max_retries = st.number_input("Max retries on failure", value=2, min_value=0, step=1)
with t4:
    skip_exists = st.checkbox("Skip ligands with existing outputs", value=False)

b1, b2 = st.columns(2)
with b1:
    exhu_backoff = st.number_input("Exhaustiveness multiplier on retry", value=1.5, min_value=1.0, step=0.1)
with b2:
    modes_backoff = st.number_input("num_modes multiplier on retry", value=1.25, min_value=1.0, step=0.05)

# Auto-detect executables and parameters from Files_for_GUI (no user input needed)
files_gui_dir = work_dir / "Files_for_GUI"
import sys

# Detect operating system
is_windows = platform.system() == "Windows"
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

# Platform and executable status (silent detection - no warnings)

# Test executables button
st.subheader("Tools")
test_btn = st.button("🔎 Test executables")

# Prepare ligand set for docking
ligand_paths: List[Path] = []
prepared_root = work_dir

if prep_btn:
    try:
        prepared = prepare_ligands_from_folder(Path(lig_src).expanduser().resolve(), prepared_root)
        st.success(f"Prepared {len(prepared)} ligands → {prepared[0].parent}")
    except Exception as e:
        st.error(str(e))

if lig_mode == "Upload now" and ligand_uploads:
    lig_dir = work_dir / "ligands_uploaded"
    lig_dir.mkdir(parents=True, exist_ok=True)
    for up in ligand_uploads:
        ligand_paths.append(_save_uploaded_file(up, lig_dir))
else:
    default_prepared_dir = prepared_root / "prepared_ligands" / "ligands_no_hydrogens"
    ligand_paths = sorted(default_prepared_dir.glob("*.pdbqt"))

# Resolve receptor path
receptor_path: Optional[Path] = None
if receptor_input_mode == "Upload file" and receptor_uploaded is not None:
    receptor_path = _save_uploaded_file(receptor_uploaded, work_dir / "receptor")
elif receptor_input_mode == "Local path" and receptor_local_path:
    receptor_path = Path(receptor_local_path).expanduser().resolve()

# Test executables
if test_btn:
    if vina_exe and vina_exe.exists():
        try:
            p = subprocess.run([str(vina_exe), "--help"], capture_output=True, text=True, timeout=10)
            st.success("Vina reachable.")
            st.code(p.stdout[:800] or p.stderr[:800])
        except Exception as e:
            st.error(f"Vina failed: {e}")
    else:
        st.error(f"Vina not found in Files_for_GUI: {files_gui_dir / 'vina.exe'}")

    for name, exe in [("AutoGrid4", autogrid_exe), ("Python", python_exe)]:
        if exe and exe.exists():
            st.info(f"{name} OK: {exe}")
        else:
            st.error(f"{name} not found in Files_for_GUI.")
    for name, pth in [("zinc_pseudo.py", zinc_pseudo_py), ("AD4_parameters.dat", base_params), ("AD4Zn.dat (extra)", extra_params)]:
        status = "✅" if (pth and pth.exists()) else "❌"
        loc = pth if pth else f"{files_gui_dir / name}"
        st.write(f"{name}: {loc} {status}")

# ==============================
# AD4 maps builder / updater (auto-detect all ligand types, force-include extras, patch missing)
# ==============================

if build_maps_btn:
    st.info("🔍 **Pre-validation checks...**")

    if receptor_path is None or not receptor_path.exists():
        st.error("❌ **Receptor file missing or invalid.**")
        st.error("**Solution:** Upload a valid PDBQT receptor file first.")
        st.stop()

    try:
        with open(receptor_path, 'r') as f:
            first_line = f.readline()
        if not first_line.startswith(('ATOM', 'HETATM', 'REMARK')):
            st.warning("⚠️ **Receptor file may not be in PDBQT format.**")
    except Exception as e:
        st.error(f"❌ **Cannot read receptor file:** {e}")
        st.stop()

    maps_prefix = Path(maps_prefix_input).expanduser().resolve()
    maps_dir = maps_prefix.parent
    maps_dir.mkdir(parents=True, exist_ok=True)

    if not autogrid_exe or not Path(autogrid_exe).exists():
        st.error("AutoGrid4 executable path is missing.")
    elif not ligand_paths:
        st.error("No ligands detected. Prepare or upload first.")
    elif not (base_params and base_params.exists()):
        st.error("AD4_parameters.dat is missing.")
    else:
        merged_params = maps_dir / "AD4_parameters_plus_ZnTZ.dat"
        merge_parameter_files(base_params, extra_params, merged_params)
        st.info(f"Using parameter file: {merged_params}")

        rec_tz = maps_dir / (Path(receptor_path).stem + "_tz.pdbqt")
        if zinc_pseudo_py and zinc_pseudo_py.exists() and python_exe.exists():
            with st.spinner("Adding tetrahedral Zn pseudoatom(s) (zinc_pseudo.py)…"):
                try:
                    zp = run_zinc_pseudo(python_exe, zinc_pseudo_py, receptor_path, rec_tz)
                    if zp.returncode != 0 or not rec_tz.exists():
                        st.error("❌ **zinc_pseudo.py failed.**")
                        st.code((zp.stdout or '') + "\n" + (zp.stderr or ''))
                        st.stop()
                    st.success(f"✅ Created {rec_tz.name}")
                except Exception as e:
                    st.error(f"❌ **Unexpected error:** {e}")
                    st.stop()
        else:
            rec_tz = Path(receptor_path)

        if normalize_OA:
            normalize_receptor_oxygen_to_OA(rec_tz, rec_tz)

        lig_types_detected = ligand_types_union(ligand_paths)
        forced = [t.strip() for t in force_extra_types.replace(",", " ").split() if t.strip()]
        lig_types_full = sorted(set(lig_types_detected).union(set(forced)))

        st.info(f"Ligand atom types detected: {' '.join(sorted(lig_types_detected)) if lig_types_detected else '(none)'}")
        if forced:
            st.warning(f"Force-including extra ligand types: {' '.join(forced)}")

        nx = max(10, int(round(float(size_x) / float(spacing))))
        ny = max(10, int(round(float(size_y) / float(spacing))))
        nz = max(10, int(round(float(size_z) / float(spacing))))
        gpf_out = maps_prefix.with_suffix(".gpf")

        write_simple_gpf(
            gpf_path=gpf_out,
            receptor_tz_filename=rec_tz.name,
            maps_prefix_basename=maps_prefix.name,
            npts_xyz=(nx, ny, nz),
            spacing=float(spacing),
            center_xyz=(float(center_x), float(center_y), float(center_z)),
            receptor_types=read_types_from_pdbqt(rec_tz),
            ligand_types=lig_types_full,
            parameter_file_rel=merged_params.name,
        )

        if rec_tz.parent != maps_dir:
            shutil.copy2(rec_tz, maps_dir / rec_tz.name)
        if ligand_paths:
            example_lig = ligand_paths[0]
            if example_lig.parent != maps_dir:
                try:
                    shutil.copy2(example_lig, maps_dir / example_lig.name)
                except Exception:
                    pass

        with st.spinner("Running AutoGrid4 to generate AD4 maps…"):
            try:
                ag = run_autogrid4(Path(autogrid_exe), maps_dir, gpf_out)
                if ag.returncode == 0:
                    st.success(f"Maps generated under {maps_dir}")
                else:
                    st.error("AutoGrid4 failed while building maps.")
                    st.code((ag.stdout or '') + "\n" + (ag.stderr or ''))
                    st.stop()
            except Exception as e:
                st.error(f"❌ **AutoGrid4 error:** {e}")
                st.stop()

        have = list_maps_present(maps_prefix)
        missing_after = [t for t in lig_types_full if t not in have]
        if missing_after:
            st.warning("Maps still missing for: " + ", ".join(missing_after))
        else:
            st.success("All requested ligand-type maps are present.")

# ==============================
# Run docking
# ==============================

run_btn = st.button("🚀 Run Docking", type="primary")

if run_btn:
    if not vina_exe or not Path(vina_exe).exists():
        st.error("Vina executable not found. Place the Linux binary in `Files_for_GUI/` and redeploy.")
        st.stop()
    if not autogrid_exe or not Path(autogrid_exe).exists():
        st.error("AutoGrid4 executable not found. Ensure the Linux binary is available in `Files_for_GUI/`.")
        st.stop()
    if receptor_path is None or not receptor_path.exists():
        st.error("Receptor file missing/invalid.")
        st.stop()
    if not ligand_paths:
        st.error("No ligand files found. Prepare or upload first.")
        st.stop()

    cx, cy, cz = float(center_x), float(center_y), float(center_z)
    if autodetect:
        auto = autodetect_metal_center(receptor_path)
        if auto:
            cx, cy, cz = auto
            st.info(f"Auto-detected metal center: {cx:.3f}, {cy:.3f}, {cz:.3f}")

    maps_prefix = Path(maps_prefix_input).expanduser().resolve()
    required_types = sorted(set(ligand_types_union(ligand_paths)).union({t.strip() for t in force_extra_types.replace(",", " ").split() if t.strip()}) or {"C", "F", "OA", "S", "NA"})
    have_maps = list_maps_present(maps_prefix)
    base_req = [
        maps_prefix.parent / f"{maps_prefix.name}.maps.fld",
        maps_prefix.parent / f"{maps_prefix.name}.e.map",
        maps_prefix.parent / f"{maps_prefix.name}.d.map",
    ]
    missing_files = [p for p in base_req if not p.exists()]
    missing_types = [t for t in required_types if t not in have_maps]
    if missing_files or missing_types:
        if missing_files:
            st.error("AD4 map headers missing:\n" + "\n".join(str(p) for p in missing_files))
        if missing_types:
            st.error("Affinity maps missing for atom types: " + ", ".join(missing_types))
        st.stop()

    out_dir = work_dir / out_dir_name
    prog = st.progress(0, text="Starting docking…")
    console = st.empty()

    def _cb(i, n, name, stat):
        prog.progress(i / n, text=f"{i}/{n} {name} — {stat}")
        console.write(f"{i}/{n}  {name}: {stat}")

    tm_mode_key = "no_timeout" if timeout_mode.startswith("No timeout") else "soft_timeout"

    with st.spinner("Running AD4 docking…"):
        rows = run_vina_batch(
            vina_exe=Path(vina_exe),
            receptor_file=receptor_path,
            ligand_files=ligand_paths,
            out_dir=out_dir,
            center=(cx, cy, cz),
            size=(float(size_x), float(size_y), float(size_z)),
            scoring="ad4",
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

    df = pd.DataFrame(rows)
    st.success("Docking complete.")
    st.dataframe(df, width="stretch")

    st.subheader("🧲 AD4 Summary")
    try:
        affinities = [float(r["Binding_Affinity"]) for r in rows if r.get("Binding_Affinity") not in ("", "N/A", None)]
        if affinities:
            st.write(f"Binding affinities range: {min(affinities):.1f} to {max(affinities):.1f} kcal/mol")
            st.write(f"Average binding affinity: {sum(affinities)/len(affinities):.1f} kcal/mol")
    except Exception:
        pass

    st.download_button(
        "⬇️ Download results CSV",
        data=_cached_file_bytes(results_to_csv_bytes(rows)),
        file_name="pfas_docking_results.csv",
        mime="text/csv",
    )
    if out_dir.exists():
        st.download_button(
            "⬇️ Download all output PDBQTs (ZIP)",
            data=_cached_file_bytes(zip_outputs(out_dir)),
            file_name=f"{out_dir.name}.zip",
            mime="application/zip",
        )

st.caption(
    """Tips:
• If you see “Affinity map for atom type X is not present”, click **Build/Update AD4 maps** with X in Force-include.
• The app now scans **all ligands** to decide which maps to make, and prints per-ligand **Score** or **missing map** in the console.
• Use **No timeout** for tough ligands; or enable soft timeouts with retries/backoff."""
)

