"""
MM-GBSA RELATIVE ENERGETIC ESTIMATION on top of the OpenMM MD trajectories
produced by `06_md_simulation_setup.py` (manuscript Methods 2.8).

Naming, fixed after external review round 2 (manuscript §2.8): the method
actually configured below (`&gb` block, `igb=5`) is MM-GBSA (Generalized
Born), not MM-PBSA (Poisson-Boltzmann). An earlier version of this docstring
used the two terms interchangeably, which is exactly the inconsistency a
computational-chemistry reviewer will catch fastest. `MMPBSA.py` remains the
correct name of the AmberTools *program* used below regardless of which
implicit-solvent model it is configured to run -- that is the program's own
name, not a description of the method. Every value produced here is used for
RELATIVE RANKING of ligand-associated energetics, never as a binding free
energy, an affinity, or a figure axis labelled ΔG. A value of, say,
-48 kcal/mol is not an experimentally meaningful affinity; it is a
model-dependent number whose usefulness is confined to ordering compounds
computed under identical settings.

Report alongside every value: frames used, frame correlation, dielectric
assumptions, entropy treatment, error estimation, and sensitivity to the
chosen trajectory window (`assess_window_sensitivity` below).

OpenMM has no native MM-GBSA implementation, and AmberTools' `MMPBSA.py`
remains the most widely cited, reviewer-recognized implementation of it, so
this script converts the OpenMM system/trajectory to Amber topology/coordinate
format via ParmEd, then drives `MMPBSA.py` as a subprocess and parses its
output. This requires AmberTools installed (`ambertools` conda package,
open-source/BSD, no separate MD engine required -- MMPBSA.py is used purely
as a post-processing tool here, the MD itself stays 100% OpenMM). If the
ParmEd conversion from an OpenMM/OpenFF system proves unreliable in practice
-- a known friction point for this specific toolchain combination --
`gmx_MMPBSA` is the fallback path.

Entropy is neglected by default (single-trajectory MM-GBSA without
normal-mode or quasi-harmonic entropy) -- state this explicitly in Methods
2.8 as done here, per the critique doc's reproducibility checklist. Add
`--entropy nmode` to the MMPBSA.py call if compute allows and cite the added
assumption.
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


def assess_window_sensitivity(cdir: Path, windows: tuple[tuple[int, int], ...] =
                               ((1, 3333), (3334, 6666), (6667, 9999))) -> pd.DataFrame:
    """Recompute the estimate over disjoint trajectory windows.

    If the estimate swings substantially between windows, the trajectory is not
    converged for this purpose and the ranking derived from it is unreliable --
    report the spread rather than the single whole-trajectory number.
    """
    rows = []
    for start, end in windows:
        mmpbsa_in = cdir / f"mmpbsa_{start}_{end}.in"
        mmpbsa_in.write_text(
            f"Window {start}-{end}\n&general\n  startframe={start}, endframe={end}, "
            f"interval=1, verbose=2,\n/\n&gb\n  igb=5, saltcon=0.15,\n/\n")
        results_dat = run_mmpbsa(
            prmtop_complex=cdir / "complex.prmtop",
            prmtop_receptor=cdir / "receptor.prmtop",
            prmtop_ligand=cdir / "ligand.prmtop",
            trajectory_dcd=cdir / "production.dcd",
            mmpbsa_in=mmpbsa_in, out_dir=cdir,
        )
        parsed = parse_mmpbsa_results(results_dat)
        parsed.update({"window_start": start, "window_end": end, "compound": cdir.name})
        rows.append(parsed)
    return pd.DataFrame(rows)


def parse_mmpbsa_results(results_dat: Path) -> dict:
    """Extract the DELTA TOTAL estimate (mean +/- std) from MMPBSA.py output.

    Returned as `mmpbsa_estimated_energy_kcal_mol` -- deliberately not named
    delta_g, so the column name itself resists being mislabelled downstream.
    """
    text = results_dat.read_text()
    for line in text.splitlines():
        if line.strip().startswith("DELTA TOTAL"):
            parts = line.split()
            return {"mmpbsa_estimated_energy_kcal_mol": float(parts[2]),
                    "std_kcal_mol": float(parts[3])}
    raise ValueError(f"Could not find DELTA TOTAL line in {results_dat}")


def batch_mmpbsa(compound_dirs: list[Path]) -> pd.DataFrame:
    """Run MM-GBSA over every shortlisted compound's MD output directory.

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
    print(df.sort_values("mmpbsa_estimated_energy_kcal_mol"))
