"""Replicate the redocking validation across 3 fixed, logged seeds per
receptor (manuscript protocol: >=3 independent seeded runs), to check
whether the 9L40 failure is a stable finding or a single-run fluke, and to
report search reproducibility (replicate agreement) separately from pose
accuracy (redocking RMSD) as manuscript section 2.4 requires.
"""
import json
import numpy as np
from vina import Vina

BOX_SIZE = 22.0
EXHAUSTIVENESS = 16
SEEDS = [1000, 1001, 1002]  # same seeds, different exhaustiveness

centers = json.load(open("box_centers.json"))
RECEPTORS = {
    "9l40_ATR_VE822":       ("9l40_protein.pdbqt", "9l40_active_site_ligand_fixed.pdbqt", centers["9l40"]["center"]),
    "9l4b_ATR_camonsertib": ("9l4b_protein.pdbqt", "9l4b_active_site_ligand_fixed.pdbqt", centers["9l4b"]["center"]),
    "2ym8_CHK1_YM8":        ("2ym8_protein.pdbqt", "2ym8_active_site_ligand_fixed.pdbqt", centers["2ym8"]["center"]),
}


def parse_pdbqt_coords(path, pose_index=0):
    coords, names = [], []
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


all_results = {}
for label, (receptor_pdbqt, ligand_pdbqt, center) in RECEPTORS.items():
    native_coords, native_names = parse_pdbqt_coords(ligand_pdbqt, pose_index=0)
    replicate_rmsds = []
    replicate_scores = []
    replicate_best_coords = []

    for seed in SEEDS:
        v = Vina(sf_name="vina", seed=seed)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=list(center), box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=9)
        out_pdbqt = f"{label}_seed{seed}_redocked.pdbqt"
        v.write_poses(out_pdbqt, n_poses=1, overwrite=True)
        energies = v.energies(n_poses=1)
        best_coords, best_names = parse_pdbqt_coords(out_pdbqt, pose_index=0)

        rmsd = float("nan")
        if best_names == native_names and len(best_coords) == len(native_coords):
            diff = native_coords - best_coords
            rmsd = float(np.sqrt((diff ** 2).sum(axis=1).mean()))

        replicate_rmsds.append(rmsd)
        replicate_scores.append(float(energies[0][0]))
        replicate_best_coords.append(best_coords)
        print(f"{label} seed={seed}: best score={energies[0][0]:.3f} kcal/mol, "
              f"RMSD to native={rmsd:.3f} A")

    # search reproducibility: pairwise RMSD between replicate best poses
    # (distinct question from accuracy -- does the search return the same
    # answer across seeds, regardless of whether that answer is correct)
    pairwise = []
    for i in range(len(replicate_best_coords)):
        for j in range(i + 1, len(replicate_best_coords)):
            d = replicate_best_coords[i] - replicate_best_coords[j]
            pairwise.append(float(np.sqrt((d ** 2).sum(axis=1).mean())))

    all_results[label] = {
        "seeds": SEEDS,
        "replicate_vina_scores_kcal_mol": replicate_scores,
        "replicate_redocking_rmsd_angstrom": replicate_rmsds,
        "mean_redocking_rmsd_angstrom": float(np.nanmean(replicate_rmsds)),
        "std_redocking_rmsd_angstrom": float(np.nanstd(replicate_rmsds)),
        "passes_2.0A_threshold_all_replicates": bool(np.all(np.array(replicate_rmsds) <= 2.0)),
        "search_reproducibility_pairwise_rmsd_angstrom": pairwise,
        "mean_search_reproducibility_rmsd": float(np.mean(pairwise)),
    }

json.dump(all_results, open("redocking_replicate_results_ex16.json", "w"), indent=2)
print("\n=== Final summary (mean +/- std RMSD across 3 seeded replicates) ===")
for k, v in all_results.items():
    status = "PASS" if v["passes_2.0A_threshold_all_replicates"] else "FAIL"
    print(f"{k}: {v['mean_redocking_rmsd_angstrom']:.3f} +/- "
          f"{v['std_redocking_rmsd_angstrom']:.3f} A -- {status} "
          f"(search reproducibility: {v['mean_search_reproducibility_rmsd']:.3f} A)")
