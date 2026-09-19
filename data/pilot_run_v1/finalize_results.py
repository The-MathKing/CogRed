"""Process the completed docking run into a ranked summary table and a
LaTeX table fragment ready to paste into the manuscript."""
import pandas as pd

raw = pd.read_csv("atr_docking_results_raw.csv")
summary = (raw.groupby(["compound_id", "source", "smiles"])["best_affinity_kcal_mol"]
           .agg(["mean", "std", "min"]).reset_index()
           .rename(columns={"mean": "affinity_mean", "std": "affinity_std", "min": "affinity_best"}))
summary = summary.sort_values("affinity_best").reset_index(drop=True)
summary.to_csv("atr_docking_results_summary_final.csv", index=False)

print(f"N compounds docked: {len(summary)}")
print(f"Affinity range: {summary['affinity_best'].min():.2f} to {summary['affinity_best'].max():.2f} kcal/mol")
seeds = summary[summary["source"] == "seed"]
derivs = summary[summary["source"] == "derivative"]
print(f"\nSeed compounds ranked:\n{seeds[['compound_id','affinity_mean','affinity_std','affinity_best']].to_string(index=False)}")
print(f"\nTop 10 derivatives:\n{derivs.head(10)[['compound_id','affinity_mean','affinity_std','affinity_best']].to_string(index=False)}")

n_better_than_best_seed = (summary["affinity_best"] < seeds["affinity_best"].min()).sum()
print(f"\nDerivatives with better (more negative) affinity than the best seed compound: {n_better_than_best_seed}")

# LaTeX table: top 10 overall + all 5 seeds for reference
top10 = summary.head(10).copy()
print("\n=== LaTeX table rows (top 10) ===")
for _, row in top10.iterrows():
    cid = row["compound_id"].replace("_", "\\_")
    src = "seed" if row["source"] == "seed" else "derivative"
    print(f"{cid} & {src} & {row['affinity_mean']:.2f} $\\pm$ {row['affinity_std']:.2f} & {row['affinity_best']:.2f} \\\\")
