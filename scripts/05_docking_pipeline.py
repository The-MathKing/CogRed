"""
Batch molecular docking wrapper around AutoDock Vina / Webina (manuscript
Section 2.4). Handles ligand/receptor preparation via Open Babel, runs Vina
with explicit, logged parameters (exhaustiveness, seed, box), and parses
results into a ranked DataFrame.

Requires on PATH: `obabel` (Open Babel) and `vina` (AutoDock Vina) or
`webina` (the browser/WASM build has no CLI -- use native Vina for batch runs
and reserve Webina for interactive single-ligand demonstration figures).

Reproducibility note (critique doc Section 5): this script logs the exact
Vina version, seed, and exhaustiveness used for every run into the results
CSV so the manuscript's methods table can be generated directly from
`data/docking_results.csv` instead of re-typed by hand.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

VINA_BIN = "vina"
OBABEL_BIN = "obabel"


@dataclass
class DockingBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 20.0
    size_y: float = 20.0
    size_z: float = 20.0


def get_vina_version() -> str:
    out = subprocess.run([VINA_BIN, "--version"], capture_output=True, text=True, check=True)
    return out.stdout.strip()


def smiles_to_pdbqt(smiles: str, out_path: Path) -> Path:
    """3D-embed + protonate at pH 7.4 + convert to PDBQT via Open Babel."""
    smi_path = out_path.with_suffix(".smi")
    smi_path.write_text(smiles + "\n")
    subprocess.run(
        [OBABEL_BIN, str(smi_path), "-O", str(out_path),
         "--gen3d", "-p", "7.4", "--partialcharge", "gasteiger"],
        check=True, capture_output=True,
    )
    return out_path


def receptor_pdb_to_pdbqt(receptor_pdb: Path, out_path: Path) -> Path:
    """Strip waters/heteroatoms and convert receptor to PDBQT (rigid receptor)."""
    subprocess.run(
        [OBABEL_BIN, str(receptor_pdb), "-O", str(out_path),
         "-xr", "--partialcharge", "gasteiger"],
        check=True, capture_output=True,
    )
    return out_path


def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path, box: DockingBox,
             out_dir: Path, exhaustiveness: int = 16, num_modes: int = 9,
             seed: int = 42) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = out_dir / f"{ligand_pdbqt.stem}_docked.pdbqt"
    log_path = out_dir / f"{ligand_pdbqt.stem}_vina.log"

    cmd = [
        VINA_BIN,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box.center_x), "--center_y", str(box.center_y), "--center_z", str(box.center_z),
        "--size_x", str(box.size_x), "--size_y", str(box.size_y), "--size_z", str(box.size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    log_path.write_text(result.stdout)

    affinities = [float(m.group(1)) for m in re.finditer(
        r"^\s*\d+\s+(-?\d+\.\d+)", result.stdout, flags=re.MULTILINE)]
    best_affinity = min(affinities) if affinities else float("nan")

    return {
        "ligand": ligand_pdbqt.stem,
        "best_affinity_kcal_mol": best_affinity,
        "all_mode_affinities": affinities,
        "exhaustiveness": exhaustiveness,
        "seed": seed,
        "vina_version": get_vina_version(),
        "docked_pose_path": str(out_pdbqt),
    }


def dock_library(compound_csv: str, receptor_pdb: str, box: DockingBox,
                  target_name: str, work_dir: str = "data/docking",
                  n_replicates: int = 3) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Dock every SMILES in `compound_csv` against one receptor, `n_replicates`
    independent Vina runs per ligand (different seeds) for score stability.

    Returns (per-run results, per-compound summary)."""
    work = Path(work_dir) / target_name
    work.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = receptor_pdb_to_pdbqt(Path(receptor_pdb), work / "receptor.pdbqt")

    df = pd.read_csv(compound_csv)
    records = []
    for i, row in df.iterrows():
        ligand_pdbqt = smiles_to_pdbqt(row["smiles"], work / f"lig_{i:04d}.pdbqt")
        for rep in range(n_replicates):
            res = run_vina(receptor_pdbqt, ligand_pdbqt, box, work / "poses",
                            seed=1000 + rep)
            res.update({"compound_index": i, "smiles": row["smiles"],
                        "target": target_name, "replicate": rep})
            records.append(res)

    results = pd.DataFrame(records)
    summary = (results.groupby(["compound_index", "smiles", "target"])
               ["best_affinity_kcal_mol"].agg(["mean", "std", "min"])
               .reset_index()
               .rename(columns={"mean": "affinity_mean", "std": "affinity_std",
                                 "min": "affinity_best"}))
    return results, summary


if __name__ == "__main__":
    # Example call -- coordinates below are placeholders; derive the real box
    # from the co-crystallized ligand centroid of the chosen PDB structure
    # (e.g., via `obabel` centroid calc or a quick RDKit/Biopython script).
    box = DockingBox(center_x=0.0, center_y=0.0, center_z=0.0)
    results, summary = dock_library(
        compound_csv="data/developability_passed_compounds.csv",
        receptor_pdb="data/receptors/CHK1.pdb",
        box=box,
        target_name="CHK1",
    )
    results.to_csv("data/docking_results_raw.csv", index=False)
    summary.to_csv("data/docking_results_summary.csv", index=False)
    print(summary.sort_values("affinity_best").head(10))
