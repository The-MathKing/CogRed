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

Docking VALIDATION (added after external review round 2, manuscript §2.4):
`redock_native_ligand` extracts a receptor's own co-crystallized ligand and
docks it back into that receptor under the identical protocol, reporting the
RMSD to the native pose. This is a different question from replicate
agreement across seeds (does the search reproduce its own answer) -- redocking
answers whether the answer is *correct*. A receptor/box combination that
fails the pre-specified RMSD threshold should not be trusted for candidate
scoring; see manuscript §2.4 and §4.3 (Limitations).
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

VINA_BIN = "vina"
OBABEL_BIN = "obabel"
DEFAULT_EXHAUSTIVENESS = 32  # Vina's own default (8) is too low for a result
                              # reported as a reproducible benchmark score.


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
             out_dir: Path, exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
             num_modes: int = 9, seed: int = 42) -> dict:
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


def extract_cocrystal_ligand(receptor_pdb: Path, ligand_resname: str,
                              out_sdf: Path) -> Path:
    """Pull one HETATM residue (the co-crystallized ligand) out of a receptor
    PDB and write it as SDF. `ligand_resname` is the 3-letter PDB ligand code
    (e.g. the code for VE-822/berzosertib in PDB 9L40) -- look this up in the
    structure's PDB header before calling.

    Residue selection uses Biopython (bond orders are not reliable from raw
    PDB coordinates alone, so this is a starting geometry for redocking, not
    a substitute for the deposited chemical-component SDF/SMILES when RCSB
    provides one -- prefer that when available and fall back to this only
    when it is not).
    """
    from Bio.PDB import PDBIO, PDBParser, Select

    class _LigandSelect(Select):
        def accept_residue(self, residue):
            return residue.get_resname() == ligand_resname

    structure = PDBParser(QUIET=True).get_structure("receptor", str(receptor_pdb))
    tmp_pdb = out_sdf.with_suffix(".ligand.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(tmp_pdb), _LigandSelect())

    subprocess.run([OBABEL_BIN, str(tmp_pdb), "-O", str(out_sdf)],
                    check=True, capture_output=True)
    return out_sdf


def redock_native_ligand(receptor_pdb: Path, native_ligand_sdf: Path,
                          box: DockingBox, out_dir: Path,
                          rmsd_pass_threshold: float = 2.0,
                          exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
                          seed: int = 1) -> dict:
    """Redocking control (manuscript §2.4): dock a receptor's own
    co-crystallized ligand back into itself and report heavy-atom RMSD to the
    native pose against a pre-specified pass threshold.

    This is the accuracy check that replicate-seed agreement (`run_vina`
    called multiple times) cannot provide -- replicates only show the search
    is reproducible, not that the pose is correct.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt = receptor_pdb_to_pdbqt(receptor_pdb, out_dir / "receptor.pdbqt")

    native_mol = Chem.SDMolSupplier(str(native_ligand_sdf), removeHs=False)[0]
    if native_mol is None:
        raise ValueError(f"Could not parse native ligand from {native_ligand_sdf}")
    native_smiles = Chem.MolToSmiles(Chem.RemoveHs(native_mol))

    ligand_pdbqt = smiles_to_pdbqt(native_smiles, out_dir / "native_ligand.pdbqt")
    docked = run_vina(receptor_pdbqt, ligand_pdbqt, box, out_dir / "redock_poses",
                       exhaustiveness=exhaustiveness, seed=seed)

    docked_sdf = out_dir / "redocked_best_pose.sdf"
    subprocess.run(
        [OBABEL_BIN, docked["docked_pose_path"], "-O", str(docked_sdf), "-f", "1", "-l", "1"],
        check=True, capture_output=True,
    )
    docked_mol = Chem.SDMolSupplier(str(docked_sdf), removeHs=False)[0]

    rmsd = float("nan")
    if docked_mol is not None:
        try:
            rmsd = AllChem.GetBestRMS(Chem.RemoveHs(docked_mol), Chem.RemoveHs(native_mol))
        except (RuntimeError, ValueError):
            rmsd = float("nan")  # atom-ordering/connectivity mismatch; inspect manually

    return {
        "receptor": receptor_pdb.stem,
        "redocking_rmsd_angstrom": rmsd,
        "pass_threshold_angstrom": rmsd_pass_threshold,
        "passes_validation": bool(rmsd <= rmsd_pass_threshold) if rmsd == rmsd else False,
        "vina_version": docked["vina_version"],
        "exhaustiveness": exhaustiveness,
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
    # ATR default receptor per manuscript §2.1: PDB 9L40 (VE-822-bound kinase
    # domain, ~2.9 A) or 9L4B (RP-3500-bound) -- NOT 5YZ0 (4.7 A, unsuitable
    # for docking). CHK1 receptor PDB ID: [PENDING, see manuscript §2.1].
    box = DockingBox(center_x=0.0, center_y=0.0, center_z=0.0)

    # Validation MUST run before any candidate is scored (manuscript §2.4).
    ligand_sdf = extract_cocrystal_ligand(
        Path("data/receptors/CHK1.pdb"), ligand_resname="LIG",  # [PENDING: real 3-letter code]
        out_sdf=Path("data/receptors/CHK1_native_ligand.sdf"),
    )
    redock_result = redock_native_ligand(
        receptor_pdb=Path("data/receptors/CHK1.pdb"),
        native_ligand_sdf=ligand_sdf,
        box=box,
        out_dir=Path("data/docking/CHK1_validation"),
    )
    print("Redocking validation:", redock_result)
    if not redock_result["passes_validation"]:
        raise SystemExit(
            "Redocking RMSD exceeds the pre-specified threshold -- this "
            "receptor/box combination should not be used for candidate "
            "scoring without an explicit stated reason (manuscript §2.4)."
        )

    results, summary = dock_library(
        compound_csv="data/developability_passed_compounds.csv",
        receptor_pdb="data/receptors/CHK1.pdb",
        box=box,
        target_name="CHK1",
    )
    results.to_csv("data/docking_results_raw.csv", index=False)
    summary.to_csv("data/docking_results_summary.csv", index=False)
    print(summary.sort_values("affinity_best").head(10))
