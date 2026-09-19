"""Symmetry-corrected RMSD (spyrmsd) for the redocking validation, addressing
the reviewer's objection that naive index-matched RMSD can overstate failure
for molecules with topological symmetry (e.g. berzosertib's two
para-disubstituted rings). Also runs the sanity controls the reviewer
requested: crystal-vs-itself (must be 0.00 A) and a deliberately perturbed
pose (must recover a known, non-zero displacement).
"""
import numpy as np
from spyrmsd import io, rmsd

NATIVE = {
    "9l40_ATR_VE822": "9l40_active_site_ligand_fixed.sdf",
    "9l4b_ATR_camonsertib": "9l4b_active_site_ligand_fixed.sdf",
    "2ym8_CHK1_YM8": "2ym8_active_site_ligand_fixed.sdf",
}
SEEDS = [1000, 1001, 1002]

print("=== Sanity control 1: crystal-vs-itself (expect 0.00 A) ===")
for label, native_path in NATIVE.items():
    mol = io.loadmol(native_path)
    r_naive = rmsd.rmsd(mol.coordinates, mol.coordinates,
                         mol.atomicnums, mol.atomicnums)
    r_symm = rmsd.symmrmsd(mol.coordinates, mol.coordinates,
                            mol.atomicnums, mol.atomicnums,
                            mol.adjacency_matrix, mol.adjacency_matrix)
    print(f"  {label}: naive={r_naive:.4f} A, symmetry-corrected={r_symm:.4f} A")

print("\n=== Sanity control 2: deliberately perturbed pose ===")
rng = np.random.default_rng(0)
for label, native_path in NATIVE.items():
    mol = io.loadmol(native_path)
    shift = np.array([1.0, 0.0, 0.0])  # exact 1.0 A rigid translation
    perturbed_coords = mol.coordinates + shift
    r_naive = rmsd.rmsd(mol.coordinates, perturbed_coords,
                         mol.atomicnums, mol.atomicnums)
    print(f"  {label}: known 1.0 A shift -> measured RMSD={r_naive:.4f} A "
          f"({'PASS' if abs(r_naive - 1.0) < 1e-6 else 'FAIL'})")

print("\n=== Symmetry-corrected redocking RMSD (3 seeded replicates) ===")
# rmsdwrapper takes Molecule objects and handles atom reordering via graph
# isomorphism internally (the independently-produced native and docked SDFs,
# from two separate Open Babel bond-perception passes, are not guaranteed to
# share atom order even though they share atom composition -- this is what
# broke the low-level rmsd.rmsd() call above and is exactly the kind of
# mismatch symmetry-aware matching is designed to handle).
results = {}
for label, native_path in NATIVE.items():
    native_mol = io.loadmol(native_path)
    symm_vals = []
    for seed in SEEDS:
        docked_path = f"{label}_seed{seed}_redocked.sdf"
        docked_mol = io.loadmol(docked_path)
        r_symm = rmsd.rmsdwrapper(native_mol, docked_mol, symmetry=True,
                                   minimize=False, strip=True)[0]
        symm_vals.append(r_symm)
        print(f"  {label} seed={seed}: symmetry-corrected={r_symm:.3f} A")
    results[label] = {
        "symm_mean": float(np.mean(symm_vals)), "symm_std": float(np.std(symm_vals)),
        "symm_values": [float(v) for v in symm_vals],
    }

print("\n=== Summary ===")
naive_from_prior_run = {  # from redocking_replicate_results.json (index-matched RMSD)
    "9l40_ATR_VE822": (2.9629, 0.0565),
    "9l4b_ATR_camonsertib": (0.6106, 0.0174),
    "2ym8_CHK1_YM8": (0.9004, 0.0066),
}
for label, r in results.items():
    naive_mean, naive_std = naive_from_prior_run[label]
    print(f"{label}: naive (index-matched) {naive_mean:.2f}+/-{naive_std:.2f} A -> "
          f"symmetry-corrected {r['symm_mean']:.2f}+/-{r['symm_std']:.2f} A")

import json
json.dump(results, open("symmrmsd_results.json", "w"), indent=2)
