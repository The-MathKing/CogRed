"""Real redocking validation run using the installed `vina` Python package
(AutoDock Vina 1.2.7, Forli lab) against the three real receptor structures
downloaded from the RCSB PDB AWS Open Data mirror (9L40, 9L4B, 2YM8).

This is the control described in manuscript section 2.4 / script
05_docking_pipeline.py's redock_native_ligand(): redock each receptor's own
co-crystallized ATP-competitive-site ligand and report RMSD to the native
pose against the pre-specified <=2.0 A pass threshold.
"""
import json
import numpy as np
from vina import Vina

BOX_SIZE = 22.0  # Angstrom, cubic box around the ligand centroid
EXHAUSTIVENESS = 32
SEED = 1000

centers = json.load(open("box_centers.json"))

RECEPTORS = {
    "9l40_ATR_VE822":       ("9l40_protein.pdbqt", "9l40_active_site_ligand_fixed.pdbqt", centers["9l40"]["center"]),
    "9l4b_ATR_camonsertib": ("9l4b_protein.pdbqt", "9l4b_active_site_ligand_fixed.pdbqt", centers["9l4b"]["center"]),
    "2ym8_CHK1_YM8":        ("2ym8_protein.pdbqt", "2ym8_active_site_ligand_fixed.pdbqt", centers["2ym8"]["center"]),
}


def parse_pdbqt_coords(path, pose_index=0):
    """Read heavy/all-atom coordinates from a PDBQT file, in file order,
    for exactly one pose (poses separated by ENDMDL in multi-pose output)."""
    coords = []
    names = []
    current_pose = 0
    with open(path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                current_pose += 1
                continue
            if current_pose != pose_index:
                continue
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
                names.append(line[12:16].strip())
    return np.array(coords), names


results = {}
for label, (receptor_pdbqt, ligand_pdbqt, center) in RECEPTORS.items():
    print(f"\n=== {label} ===")
    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=[BOX_SIZE] * 3)

    # Score the native (crystal) pose itself, before any search -- this is
    # the pose Vina's scoring function assigns to the deposited coordinates.
    native_energy = v.score()
    print(f"Native-pose Vina score: {native_energy}")

    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
    out_pdbqt = f"{label}_redocked.pdbqt"
    v.write_poses(out_pdbqt, n_poses=9, overwrite=True)
    docked_energies = v.energies(n_poses=9)

    native_coords, native_names = parse_pdbqt_coords(ligand_pdbqt, pose_index=0)
    best_coords, best_names = parse_pdbqt_coords(out_pdbqt, pose_index=0)

    same_order = (native_names == best_names)
    rmsd = float(np.nan)
    if same_order and len(native_coords) == len(best_coords):
        diff = native_coords - best_coords
        rmsd = float(np.sqrt((diff ** 2).sum(axis=1).mean()))
    else:
        print("  WARNING: atom order mismatch between native and docked pose "
              "files -- RMSD not computed by direct index correspondence.")

    results[label] = {
        "native_vina_score_kcal_mol": float(native_energy[0]),
        "best_pose_vina_score_kcal_mol": float(docked_energies[0][0]),
        "redocking_rmsd_angstrom": rmsd,
        "passes_2.0A_threshold": bool(rmsd <= 2.0) if rmsd == rmsd else False,
        "n_ligand_atoms": len(native_coords),
        "box_center": center,
        "box_size_angstrom": BOX_SIZE,
        "exhaustiveness": EXHAUSTIVENESS,
    }
    print(f"  Redocking RMSD: {rmsd:.3f} A -- "
          f"{'PASS' if results[label]['passes_2.0A_threshold'] else 'FAIL'} "
          f"(threshold 2.0 A)")

json.dump(results, open("redocking_validation_results.json", "w"), indent=2)
print("\n=== Summary ===")
for k, v in results.items():
    print(k, v)
