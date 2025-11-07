"""
MetalloDock Streamlit Application
---------------------------------

A Streamlit-based interface for running AutoDock4-style (AD4) docking workflows
backed by AutoDock Vina on Linux systems. The app provides:

- Home page with quick overview
- Documentation page derived from the 2025 docking protocol
- Demo page running a fixed ligand/receptor/gridbox example
- Full MetalloDock page allowing users to configure docking parameters

All docking runs reuse the AutoDock4Zn parameterisation and the bundled
AutoDock/AutoGrid/Vina binaries shipped with MetalloDock.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import textwrap
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import streamlit as st

# ------------------------------------------------------------------------------
# Configuration & Constants
# ------------------------------------------------------------------------------

APP_TITLE = "MetalloDock"
APP_ICON = "🧪"

BASE_DIR = Path(__file__).resolve().parent
FILES_DIR = BASE_DIR / "Files_for_GUI"
DEFAULT_PARAMETER_FILE = FILES_DIR / "AD4Zn.dat"
DEFAULT_BINS = {
    "autogrid": FILES_DIR / "autogrid4",
    "autodock": FILES_DIR / "autodock4",
    "vina": FILES_DIR / "vina",
}
ZINC_PSEUDO_SCRIPT = FILES_DIR / "zinc_pseudo.py"

DEMO_RECEPTOR = BASE_DIR / "Carbonic_Anhydrase_I.pdbqt"
DEMO_LIGAND = BASE_DIR / "PFOA_Test_Ligand.pdbqt"
DEMO_GRID_CENTER = (29.951, 0.420, -4.735)
DEMO_GRID_SIZE = (16, 18, 16)

DEFAULT_SPACING = 0.375
DEFAULT_EXHAUSTIVENESS = 16
DEFAULT_CPU = max(1, os.cpu_count() or 1)
DEFAULT_NUM_MODES = 9

OUTPUT_ROOT = BASE_DIR / "metallodock_runs"
OUTPUT_ROOT.mkdir(exist_ok=True)

# ------------------------------------------------------------------------------
# Streamlit Page Setup & Styling
# ------------------------------------------------------------------------------

st.set_page_config(
    page_title=APP_TITLE,
    page_icon=APP_ICON,
    layout="wide",
    initial_sidebar_state="expanded",
)

_CUSTOM_STYLE = """
<style>
.main-header {
    font-size: 2.6rem;
    font-weight: 700;
    color: #1f77b4;
    text-align: center;
    margin-bottom: 1.5rem;
}
.sub-header {
    font-size: 1.6rem;
    font-weight: 600;
    color: #1f77b4;
    margin-top: 1.2rem;
}
.section-card {
    background-color: #f8f9fa;
    border-radius: 12px;
    padding: 1.2rem 1.5rem;
    margin-bottom: 1.4rem;
    border: 1px solid #e1e5eb;
}
.stTabs [role="tablist"] {
    justify-content: center;
}
.stTabs [role="tab"] {
    font-weight: 600;
}
</style>
"""

st.markdown(_CUSTOM_STYLE, unsafe_allow_html=True)


# ------------------------------------------------------------------------------
# Dataclasses & Helpers
# ------------------------------------------------------------------------------

class DockingError(RuntimeError):
    """Raised when a docking step fails."""


@dataclass
class DockingInputs:
    receptor: Path
    ligand: Path
    parameter_file: Path
    grid_center: Tuple[float, float, float]
    grid_size: Tuple[int, int, int]
    spacing: float = DEFAULT_SPACING
    exhaustiveness: int = DEFAULT_EXHAUSTIVENESS
    cpu: int = DEFAULT_CPU
    num_modes: int = DEFAULT_NUM_MODES
    map_prefix: Optional[str] = None


@dataclass
class DockingOutputs:
    workdir: Path
    receptor_tz: Path
    gpf_file: Path
    grid_log: Optional[Path]
    vina_out: Optional[Path]
    vina_log: Optional[Path]
    metadata: Dict[str, str]


def _check_executable(path: Path) -> None:
    if not path.exists():
        raise DockingError(f"Executable not found: {path}")
    if not os.access(path, os.X_OK):
        raise DockingError(f"Executable is not marked as runnable: {path}")


def _ensure_binaries() -> Dict[str, Path]:
    bins = {}
    for key, path in DEFAULT_BINS.items():
        bins[key] = path
    for name, exe in bins.items():
        try:
            _check_executable(exe)
        except DockingError as exc:
            raise DockingError(
                f"{name} binary issue: {exc}. "
                "Ensure the Linux binaries are present and executable."
            ) from exc
    return bins


def _parse_pdbqt_atom_types(path: Path) -> List[str]:
    types: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                atom_type = line[77:79].strip().upper()
                if atom_type and atom_type not in types:
                    types.append(atom_type)
    return sorted(types)


def _run_subprocess(
    command: Sequence[str],
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    result = subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    return result


def _run_zinc_pseudo(receptor: Path, output: Path) -> subprocess.CompletedProcess:
    if not ZINC_PSEUDO_SCRIPT.exists():
        raise DockingError("zinc_pseudo.py script is missing.")
    command = [
        str(Path(shutil.which("python") or "python")),
        str(ZINC_PSEUDO_SCRIPT),
        "-r",
        str(receptor),
        "-o",
        str(output),
    ]
    return _run_subprocess(command, cwd=FILES_DIR)


def _build_gpf_file(
    inputs: DockingInputs,
    workdir: Path,
    map_prefix: str,
    receptor_tz: Path,
) -> Path:
    atom_types = _parse_pdbqt_atom_types(inputs.ligand)
    atom_types = [atype for atype in atom_types if atype] or ["C"]

    gpf_lines = [
        f"npts {inputs.grid_size[0]} {inputs.grid_size[1]} {inputs.grid_size[2]}",
        f"gridcenter {inputs.grid_center[0]:.3f} "
        f"{inputs.grid_center[1]:.3f} "
        f"{inputs.grid_center[2]:.3f}",
        f"spacing {inputs.spacing:.4f}",
        f"gridfld {map_prefix}.fld",
        f"receptor {receptor_tz.name}",
        f"ligand {inputs.ligand.name}",
    ]

    for atype in atom_types:
        gpf_lines.append(f"map {map_prefix}.{atype}.map")

    gpf_lines.extend(
        [
            f"elecmap {map_prefix}.e.map",
            f"dsolvmap {map_prefix}.d.map",
            "dielectric -0.1465",
            f"parameter_file {inputs.parameter_file}",
        ]
    )

    gpf_path = workdir / f"{map_prefix}.gpf"
    gpf_path.write_text("\n".join(gpf_lines) + "\n", encoding="utf-8")
    return gpf_path


def _parse_vina_log(log_path: Path) -> Dict[str, List[float]]:
    scores: List[float] = []
    with log_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.strip().startswith("-----+"):
                break
        for line in handle:
            if not line.strip():
                break
            parts = line.split()
            try:
                score = float(parts[1])
            except (IndexError, ValueError):
                continue
            scores.append(score)
    return {"binding_affinities": scores}


def run_docking(inputs: DockingInputs, demo_mode: bool = False) -> DockingOutputs:
    bins = _ensure_binaries()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_label = inputs.map_prefix or ("demo" if demo_mode else f"run_{timestamp}")
    workdir = OUTPUT_ROOT / run_label
    workdir.mkdir(parents=True, exist_ok=True)

    receptor_local = workdir / inputs.receptor.name
    ligand_local = workdir / inputs.ligand.name
    shutil.copy2(inputs.receptor, receptor_local)
    shutil.copy2(inputs.ligand, ligand_local)

    receptor_tz = workdir / f"{inputs.receptor.stem}_tz.pdbqt"
    zinc_proc = _run_zinc_pseudo(receptor_local, receptor_tz)
    if zinc_proc.returncode != 0:
        raise DockingError(
            "Failed to add TZ pseudo atoms:\n"
            f"{zinc_proc.stdout}\n{zinc_proc.stderr}"
        )

    map_prefix = inputs.map_prefix or f"{inputs.receptor.stem}_maps"
    gpf_path = _build_gpf_file(inputs, workdir, map_prefix, receptor_tz)

    grid_log = workdir / f"{map_prefix}.glg"
    autogrid_cmd = [
        str(bins["autogrid"]),
        "-p",
        str(gpf_path),
        "-l",
        str(grid_log),
    ]
    grid_proc = _run_subprocess(autogrid_cmd, cwd=workdir)
    if grid_proc.returncode != 0:
        raise DockingError(
            "AutoGrid failed. Review the log for details:\n"
            f"{grid_proc.stdout}\n{grid_proc.stderr}"
        )

    vina_out = workdir / f"{inputs.ligand.stem}_ad4_vina_out.pdbqt"
    vina_log = workdir / f"{inputs.ligand.stem}_ad4_vina.log"
    vina_cmd = [
        str(bins["vina"]),
        "--ligand",
        str(ligand_local),
        "--maps",
        map_prefix,
        "--scoring",
        "ad4",
        "--exhaustiveness",
        str(inputs.exhaustiveness),
        "--num_modes",
        str(inputs.num_modes),
        "--cpu",
        str(inputs.cpu),
        "--out",
        str(vina_out),
        "--log",
        str(vina_log),
    ]
    vina_proc = _run_subprocess(vina_cmd, cwd=workdir)
    if vina_proc.returncode != 0:
        raise DockingError(
            "AutoDock Vina (AD4 scoring) failed:\n"
            f"{vina_proc.stdout}\n{vina_proc.stderr}"
        )

    metadata = {
        "command_autogrid": " ".join(autogrid_cmd),
        "command_vina": " ".join(vina_cmd),
        "zinc_stdout": zinc_proc.stdout,
        "zinc_stderr": zinc_proc.stderr,
        "autogrid_stdout": grid_proc.stdout,
        "autogrid_stderr": grid_proc.stderr,
        "vina_stdout": vina_proc.stdout,
        "vina_stderr": vina_proc.stderr,
    }

    outputs = DockingOutputs(
        workdir=workdir,
        receptor_tz=receptor_tz,
        gpf_file=gpf_path,
        grid_log=grid_log if grid_log.exists() else None,
        vina_out=vina_out if vina_out.exists() else None,
        vina_log=vina_log if vina_log.exists() else None,
        metadata=metadata,
    )
    return outputs


# ------------------------------------------------------------------------------
# Documentation & Content Helpers
# ------------------------------------------------------------------------------

def _protocol_sections() -> List[Tuple[str, str]]:
    raw_steps = [
        (
            "Prerequisites",
            """
            - Prepare receptor and ligands as PDBQT files (MGLTools or OBabel).
            - Bundle all ligands within a single directory for batch docking.
            - Determine the grid box center using the endogenous/experimental ligand.
            - Keep `AD4Zn.dat` and the Vina/AutoDock binaries together.
            """,
        ),
        (
            "1. Add Tetrahedral Zn Pseudo Atoms",
            """
            Use `zinc_pseudo.py` to inject `TZ` pseudo-atoms around zinc sites and set the
            zinc formal charge to zero.
            """,
        ),
        (
            "2. Define Grid Maps",
            """
            Construct a `.gpf` file defining the grid dimensions (`npts`), Cartesian center,
            spacing, parameter file, and the atom maps to generate. Include additional map
            types such as `TZ`, `NA`, `SA`, or `F` whenever they appear in the ligand or
            receptor.
            """,
        ),
        (
            "3. Run AutoGrid4",
            """
            Execute `autogrid4 -p <protein>.gpf -l <protein>.glg` to build affinity,
            electrostatic, and desolvation maps used by AD4/Vina.
            """,
        ),
        (
            "4. Launch Docking",
            """
            Run AutoDock Vina with `--scoring ad4`, pointing to the generated maps via
            `--maps <prefix>` and providing ligand PDBQT, exhaustiveness, and CPU settings.
            """,
        ),
        (
            "5. Analyse Results",
            """
            Inspect the Vina log (`*.log`) for ranked binding affinities and
            convert/output poses from the generated `*_ad4_vina_out.pdbqt`.
            """,
        ),
    ]

    sections: List[Tuple[str, str]] = []
    for title, body in raw_steps:
        sections.append((title, textwrap.dedent(body)))
    return sections


def _grid_display(center: Tuple[float, float, float], size: Tuple[int, int, int]) -> str:
    return (
        f"Grid center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})\n"
        f"Dimensions (Å): {size[0]} × {size[1]} × {size[2]}"
    )


# ------------------------------------------------------------------------------
# Page Renderers
# ------------------------------------------------------------------------------

def render_home() -> None:
    st.markdown('<h1 class="main-header">MetalloDock</h1>', unsafe_allow_html=True)
    st.markdown(
        "Streamlined AutoDock4Zn workflow with AutoDock Vina scoring — "
        "optimized for metalloprotein-ligand docking in Linux environments."
    )

    col1, col2 = st.columns([2, 1])
    with col1:
        st.markdown(
            """
            ### Why MetalloDock?
            - One-click demo reproducing the official Carbonic Anhydrase I example.
            - Full flexibility to tweak grid boxes, exhaustiveness, CPU usage, and output paths.
            - Built-in Zn pseudo-atom handling using the official AutoDock4Zn toolkit.
            - Persists every run under `metallodock_runs/` for traceability and re-analysis.
            """
        )
    with col2:
        st.markdown(
            """
            ### Quick Facts
            - **Engine:** AutoDock Vina (AD4 scoring)
            - **Maps:** AutoGrid4 + AD4Zn parameters
            - **Outputs:** PDBQT poses & detailed logs
            - **Best for:** GPCRs and other zinc-rich receptors
            """
        )

    st.markdown(
        """
        #### Workflow Overview
        1. Upload or reference receptor & ligand PDBQT files.
        2. Define grid box using endogenous ligand coordinates.
        3. Auto-generate Zn pseudo-atoms and AD4 maps.
        4. Run docking via Vina and review ranked poses.
        """
    )


def render_documentation() -> None:
    st.markdown('<h1 class="main-header">Protocol Documentation</h1>', unsafe_allow_html=True)
    st.markdown(
        "Derived from the official 2025 MetalloDock workflow. Each section summarises "
        "the critical actions required for reproducible AD4 docking."
    )

    for title, body in _protocol_sections():
        with st.container():
            st.markdown(f'<div class="section-card"><div class="sub-header">{title}</div>', unsafe_allow_html=True)
            st.write(body)
            st.markdown("</div>", unsafe_allow_html=True)

    st.info(
        "Need a printable copy? Save this page as PDF from your browser's print dialog."
    )


def _render_run_summary(outputs: DockingOutputs) -> None:
    st.success("Docking completed successfully.")
    st.markdown(
        f"""
        **Work directory:** `{outputs.workdir}`

        **Generated files**
        - TZ receptor: `{outputs.receptor_tz.name}`
        - Grid params: `{outputs.gpf_file.name}`
        - Vina poses: `{outputs.vina_out.name if outputs.vina_out else 'n/a'}`
        - Vina log: `{outputs.vina_log.name if outputs.vina_log else 'n/a'}`
        """
    )

    if outputs.vina_log and outputs.vina_log.exists():
        scores = _parse_vina_log(outputs.vina_log)
        if scores["binding_affinities"]:
            st.markdown("**Top Binding Affinities (kcal/mol)**")
            st.table(
                {
                    "#": list(range(1, len(scores["binding_affinities"]) + 1)),
                    "Affinity": [f"{s:.2f}" for s in scores["binding_affinities"]],
                }
            )

        with outputs.vina_log.open("r", encoding="utf-8", errors="ignore") as fh:
            log_text = fh.read()
        st.markdown("**Raw Vina Log**")
        st.code(log_text, language="text")

    with st.expander("Execution Metadata"):
        st.json(outputs.metadata)


def render_demo() -> None:
    st.markdown('<h1 class="main-header">Demo Docking</h1>', unsafe_allow_html=True)
    st.markdown(
        "Reproduce the official demo: PFOA ligand docked into Carbonic Anhydrase I "
        "using the AD4 scoring function."
    )

    st.markdown(
        f"""
        **Ligand:** `{DEMO_LIGAND.name}`

        **Receptor:** `{DEMO_RECEPTOR.name}`

        { _grid_display(DEMO_GRID_CENTER, DEMO_GRID_SIZE) }
        """
    )

    run_placeholder = st.empty()
    if st.button("Run Demo Docking", type="primary"):
        with run_placeholder, st.spinner("Running demo docking..."):
            try:
                outputs = run_docking(
                    DockingInputs(
                        receptor=DEMO_RECEPTOR,
                        ligand=DEMO_LIGAND,
                        parameter_file=DEFAULT_PARAMETER_FILE,
                        grid_center=DEMO_GRID_CENTER,
                        grid_size=DEMO_GRID_SIZE,
                        map_prefix="demo_maps",
                    ),
                    demo_mode=True,
                )
            except DockingError as err:
                st.error(str(err))
            except FileNotFoundError as err:
                st.error(f"Input missing: {err}")
            else:
                _render_run_summary(outputs)
    else:
        st.info("Click **Run Demo Docking** to execute the full AD4 pipeline with the bundled example data.")


def _handle_file_upload(upload, destination_dir: Path) -> Optional[Path]:
    if upload is None:
        return None
    destination_dir.mkdir(parents=True, exist_ok=True)
    out_path = destination_dir / upload.name
    out_path.write_bytes(upload.getbuffer())
    return out_path


def render_metallodock() -> None:
    st.markdown('<h1 class="main-header">MetalloDock Workflow</h1>', unsafe_allow_html=True)
    st.markdown("Configure and launch customised docking runs using AD4 scoring.")

    with st.sidebar:
        st.markdown("### Docking Inputs")
        receptor_upload = st.file_uploader("Upload receptor (PDBQT)", type=["pdbqt"])
        ligand_upload = st.file_uploader("Upload ligand (PDBQT)", type=["pdbqt"])
        receptor_path_text = st.text_input("…or receptor path", value=str(DEMO_RECEPTOR))
        ligand_path_text = st.text_input("…or ligand path", value=str(DEMO_LIGAND))
        parameter_path = st.text_input("AD4 parameter file", value=str(DEFAULT_PARAMETER_FILE))

        st.markdown("### Grid Parameters")
        col_a, col_b, col_c = st.columns(3)
        with col_a:
            grid_x = st.number_input("Grid size X", min_value=4, max_value=128, value=DEMO_GRID_SIZE[0], step=1)
        with col_b:
            grid_y = st.number_input("Grid size Y", min_value=4, max_value=128, value=DEMO_GRID_SIZE[1], step=1)
        with col_c:
            grid_z = st.number_input("Grid size Z", min_value=4, max_value=128, value=DEMO_GRID_SIZE[2], step=1)

        col1, col2, col3 = st.columns(3)
        with col1:
            center_x = st.number_input("Center X", value=DEMO_GRID_CENTER[0], format="%.3f")
        with col2:
            center_y = st.number_input("Center Y", value=DEMO_GRID_CENTER[1], format="%.3f")
        with col3:
            center_z = st.number_input("Center Z", value=DEMO_GRID_CENTER[2], format="%.3f")

        spacing = st.number_input("Grid spacing (Å)", min_value=0.2, max_value=1.0, value=DEFAULT_SPACING, step=0.05, format="%.3f")

        st.markdown("### Vina Controls")
        exhaustiveness = st.number_input("Exhaustiveness", min_value=1, max_value=128, value=DEFAULT_EXHAUSTIVENESS, step=1)
        num_modes = st.number_input("Number of modes", min_value=1, max_value=50, value=DEFAULT_NUM_MODES, step=1)
        cpu = st.number_input("CPUs", min_value=1, max_value=DEFAULT_CPU, value=min(DEFAULT_CPU, 8), step=1)

        map_prefix = st.text_input("Map prefix (optional)", value="")

    work_container = st.container()
    with work_container:
        run_area = st.empty()
        if st.button("Run Docking", type="primary"):
            with st.spinner("Launching docking workflow..."):
                receptor_upload_path = _handle_file_upload(receptor_upload, OUTPUT_ROOT / "uploads")
                ligand_upload_path = _handle_file_upload(ligand_upload, OUTPUT_ROOT / "uploads")

                receptor_path = Path(receptor_upload_path or receptor_path_text).expanduser()
                ligand_path = Path(ligand_upload_path or ligand_path_text).expanduser()
                parameter_file_path = Path(parameter_path).expanduser()

                inputs = DockingInputs(
                    receptor=receptor_path,
                    ligand=ligand_path,
                    parameter_file=parameter_file_path,
                    grid_center=(center_x, center_y, center_z),
                    grid_size=(int(grid_x), int(grid_y), int(grid_z)),
                    spacing=float(spacing),
                    exhaustiveness=int(exhaustiveness),
                    cpu=int(cpu),
                    num_modes=int(num_modes),
                    map_prefix=map_prefix.strip() or None,
                )

                try:
                    outputs = run_docking(inputs)
                except DockingError as err:
                    st.error(str(err))
                except FileNotFoundError as err:
                    st.error(f"Input missing: {err}")
                else:
                    _render_run_summary(outputs)
        else:
            st.info("Upload/select your inputs and press **Run Docking** to begin.")


# ------------------------------------------------------------------------------
# Main Navigation
# ------------------------------------------------------------------------------

PAGES = {
    "Home": render_home,
    "Demo": render_demo,
    "MetalloDock": render_metallodock,
    "Documentation": render_documentation,
}


def main() -> None:
    st.sidebar.markdown('<div class="sidebar-header">Navigation</div>', unsafe_allow_html=True)
    selection = st.sidebar.radio("Navigate", list(PAGES.keys()), index=0)
    render_page = PAGES[selection]
    render_page()


if __name__ == "__main__":
    main()

