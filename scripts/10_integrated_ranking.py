"""
Integrated multi-objective ranking with weight sensitivity analysis
(manuscript Methods 2.10 / Results 3.9). Added after external review.

The problem this addresses: a desirability function with hand-chosen weights
can produce whatever ranking the author wants. "Affinity 40%, selectivity 25%,
developability 20%, synthetic accessibility 15%" looks principled and is
arbitrary. A reviewer will ask why those numbers, and the honest answer is
usually that they were not derived from anything.

Two defenses, both implemented here:

  1. PRE-SPECIFICATION. Weights are declared in `PRESPECIFIED_WEIGHTS` with a
     written rationale, fixed before the ranking is computed and not adjusted
     afterwards.

  2. SENSITIVITY ANALYSIS. `weight_sensitivity_analysis` samples weight vectors
     from a Dirichlet distribution over the simplex and reports how often each
     compound lands in the top-k across thousands of alternative weightings.
     A candidate that is top-ranked under 5% of weightings is an artifact of
     the chosen weights; one that is top-ranked under 85% is a robust result.
     Report the stability frequency next to every candidate in the manuscript
     table -- it is the number that tells a reader whether the ranking means
     anything.

Note on correlated objectives: docking score, MM-GBSA estimate, and MD
stability are not independent evidence dimensions -- they share force-field
and pose assumptions. `objective_correlation_matrix` reports their empirical
correlation so the manuscript can state this rather than implying four
independent confirmations of one result.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Pre-specified before ranking. Rationale recorded here so it is auditable:
#   - predicted_engagement weighted highest because it is the only objective
#     directly tied to the target;
#   - selectivity weighted nearly as high because off-target liability is the
#     main practical failure mode for kinase-directed chemical matter;
#   - developability and synthetic accessibility weighted lower because they
#     are coarse heuristics, not measurements.
PRESPECIFIED_WEIGHTS = {
    "predicted_engagement": 0.35,
    "predicted_selectivity": 0.30,
    "conformational_stability": 0.15,
    "predicted_developability": 0.12,
    "synthetic_accessibility": 0.08,
}


def normalize_objectives(df: pd.DataFrame, objective_cols: list[str],
                          higher_is_better: dict[str, bool]) -> pd.DataFrame:
    """Min-max scale each objective to [0, 1], orienting all so higher = better."""
    out = pd.DataFrame(index=df.index)
    for col in objective_cols:
        series = df[col].astype(float)
        rng = series.max() - series.min()
        scaled = (series - series.min()) / rng if rng > 0 else pd.Series(0.5, index=series.index)
        out[col] = scaled if higher_is_better.get(col, True) else 1.0 - scaled
    return out


def desirability_score(normalized: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    """Weighted arithmetic mean of normalized objectives."""
    cols = [c for c in weights if c in normalized.columns]
    w = np.array([weights[c] for c in cols], dtype=float)
    w = w / w.sum()
    return pd.Series(normalized[cols].to_numpy() @ w, index=normalized.index)


def geometric_desirability(normalized: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    """Weighted geometric mean -- a compound scoring ~0 on any objective is
    penalized much harder than under the arithmetic mean. Reported alongside
    the arithmetic score because the choice between them is itself a modeling
    decision that can change the ranking."""
    cols = [c for c in weights if c in normalized.columns]
    w = np.array([weights[c] for c in cols], dtype=float)
    w = w / w.sum()
    vals = np.clip(normalized[cols].to_numpy(), 1e-9, None)
    return pd.Series(np.exp(np.log(vals) @ w), index=normalized.index)


def weight_sensitivity_analysis(normalized: pd.DataFrame,
                                 base_weights: dict[str, float] = PRESPECIFIED_WEIGHTS,
                                 n_samples: int = 5000, top_k: int = 5,
                                 concentration: float = 5.0,
                                 seed: int = 20260913) -> pd.DataFrame:
    """How often does each compound reach the top-k under perturbed weights?

    Weight vectors are drawn from Dirichlet(concentration * base_weights),
    which concentrates around the pre-specified weights while exploring the
    simplex. Lower `concentration` = wider exploration.
    """
    rng = np.random.default_rng(seed)
    cols = [c for c in base_weights if c in normalized.columns]
    alpha = np.array([base_weights[c] for c in cols], dtype=float)
    alpha = concentration * alpha / alpha.sum()

    counts = pd.Series(0, index=normalized.index, dtype=int)
    values = normalized[cols].to_numpy()
    for _ in range(n_samples):
        w = rng.dirichlet(alpha)
        scores = pd.Series(values @ w, index=normalized.index)
        counts.loc[scores.nlargest(top_k).index] += 1

    return pd.DataFrame({
        "top_k_frequency": counts / n_samples,
        "n_samples": n_samples,
        "top_k": top_k,
    }).sort_values("top_k_frequency", ascending=False)


def objective_correlation_matrix(normalized: pd.DataFrame) -> pd.DataFrame:
    """Spearman correlation between objectives.

    High correlation between docking score, MM-GBSA estimate, and MD stability
    means they are not independent lines of evidence. State this in Discussion
    rather than presenting them as mutual confirmation.
    """
    return normalized.corr(method="spearman")


def build_integrated_ranking(input_csv: str = "data/candidate_objectives.csv",
                              top_k: int = 5) -> pd.DataFrame:
    """Expects one row per candidate with the objective columns named in
    PRESPECIFIED_WEIGHTS (raw, un-normalized values)."""
    df = pd.read_csv(input_csv).set_index("smiles")

    higher_is_better = {
        "predicted_engagement": True,       # z-scored, higher = stronger
        "predicted_selectivity": True,      # selectivity margin in z units
        "conformational_stability": False,  # RMSD -- lower is better
        "predicted_developability": True,   # rules passed
        "synthetic_accessibility": False,   # SA score -- lower is easier
    }
    objective_cols = [c for c in PRESPECIFIED_WEIGHTS if c in df.columns]
    normalized = normalize_objectives(df, objective_cols, higher_is_better)

    ranking = pd.DataFrame({
        "desirability_arithmetic": desirability_score(normalized, PRESPECIFIED_WEIGHTS),
        "desirability_geometric": geometric_desirability(normalized, PRESPECIFIED_WEIGHTS),
    })
    stability = weight_sensitivity_analysis(normalized, top_k=top_k)
    ranking = ranking.join(stability["top_k_frequency"]).sort_values(
        "desirability_arithmetic", ascending=False)
    return ranking, objective_correlation_matrix(normalized)


if __name__ == "__main__":
    ranking, correlations = build_integrated_ranking()
    ranking.to_csv("data/integrated_ranking.csv")
    correlations.to_csv("data/objective_correlations.csv")

    print(ranking.head(10))
    print("\nObjective correlations (Spearman):")
    print(correlations)

    fragile = ranking[ranking["top_k_frequency"] < 0.5].index.tolist()
    if fragile:
        print(f"\n{len(fragile)} candidate(s) reach the top-{5} under fewer than "
              f"half of sampled weightings -- report their instability rather "
              f"than presenting them as robust hits.")
