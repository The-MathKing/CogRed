"""
MM-PBSA binding free energy estimation on top of the OpenMM MD trajectories
produced by `06_md_simulation_setup.py` (manuscript Section 2.8).

OpenMM has no native MM-PBSA implementation, and AmberTools' `MMPBSA.py`
remains the most widely cited, reviewer-recognized implementation, so this
script converts the OpenMM system/trajectory to Amber topology/coordinate
format via ParmEd, then drives `MMPBSA.py` as a subprocess and parses its
output. This requires AmberTools installed (`ambertools` conda package,
open-source/BSD, no separate MD engine required -- MMPBSA.py is used purely
as a post-processing tool here, the MD itself stays 100% OpenMM).

Entropy is neglected by default (interaction/MM-PBSA without normal-mode or
quasi-harmonic entropy) -- state this explicitly in Methods 2.8 as done here,
per the critique doc's reproducibility checklist. Add `--entropy nmode` to
the MMPBSA.py call if compute allows and cite the added assumption.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd
import parmed as pmd


def convert_openmm_to_amber(system, modeller, out_dir: Path,
                             complex_prefix: str = "complex") -> tuple[Path, Path]:
    """Use ParmEd to write Amber-format prmtop/inpcrd from an OpenMM System."""
    structure = pmd.openmm.load_topology(modeller.topology, system, xyz=modeller.positions)
    prmtop = out_dir / f"{complex_prefix}.prmtop"
    inpcrd = out_dir / f"{complex_prefix}.inpcrd"
    structure.save(str(prmtop), overwrite=True)
    structure.save(str(inpcrd), overwrite=True)
    return prmtop, inpcrd


def write_mmpbsa_input(out_path: Path, igb: int = 5, istrng: float = 0.15) -> Path:
    """Standard single-trajectory MM-GBSA input file for MMPBSA.py."""
    content = f"""Single-trajectory MM-GBSA input
&general
  startframe=1, endframe=9999, interval=1,
  verbose=2,
/
&gb
  igb={igb}, saltcon={istrng},
/
"""
    out_path.write_text(content)
    return out_path


def run_mmpbsa(prmtop_complex: Path, prmtop_receptor: Path, prmtop_ligand: Path,
               trajectory_dcd: Path, mmpbsa_in: Path, out_dir: Path) -> Path:
    out_dat = out_dir / "FINAL_RESULTS_MMPBSA.dat"
    cmd = [
        "MMPBSA.py",
        "-O",
        "-i", str(mmpbsa_in),
        "-cp", str(prmtop_complex),
        "-rp", str(prmtop_receptor),
        "-lp", str(prmtop_ligand),
        "-y", str(trajectory_dcd),
        "-o", str(out_dat),
    ]
    subprocess.run(cmd, cwd=out_dir, check=True, capture_output=True, text=True)
    return out_dat


def parse_mmpbsa_results(results_dat: Path) -> dict:
    """Extract DELTA TOTAL binding free energy (mean +/- std) from MMPBSA.py output."""
    text = results_dat.read_text()
    for line in text.splitlines():
        if line.strip().startswith("DELTA TOTAL"):
            parts = line.split()
            return {"delta_g_kcal_mol": float(parts[2]), "std_kcal_mol": float(parts[3])}
    raise ValueError(f"Could not find DELTA TOTAL line in {results_dat}")


def batch_mmpbsa(compound_dirs: list[Path]) -> pd.DataFrame:
    """Run MM-PBSA over every shortlisted compound's MD output directory.

    Each `compound_dir` is expected to already contain complex/receptor/ligand
    prmtop files (from `convert_openmm_to_amber`) and `production.dcd` (from
    the MD script).
    """
    rows = []
    for cdir in compound_dirs:
        mmpbsa_in = write_mmpbsa_input(cdir / "mmpbsa.in")
        results_dat = run_mmpbsa(
            prmtop_complex=cdir / "complex.prmtop",
            prmtop_receptor=cdir / "receptor.prmtop",
            prmtop_ligand=cdir / "ligand.prmtop",
            trajectory_dcd=cdir / "production.dcd",
            mmpbsa_in=mmpbsa_in,
            out_dir=cdir,
        )
        parsed = parse_mmpbsa_results(results_dat)
        parsed["compound"] = cdir.name
        rows.append(parsed)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    shortlist_dirs = sorted(Path("data/md").glob("*"))
    df = batch_mmpbsa(shortlist_dirs)
    df.to_csv("data/mmpbsa_results.csv", index=False)
    print(df.sort_values("delta_g_kcal_mol"))
