"""
In silico developability and toxicity PREDICTION (manuscript Methods 2.6).

Renamed from "ADMET triage" after external review: nothing here is measured.
Every value this module produces is a prediction from a descriptor rule or a
QSAR web tool, and the manuscript must report it as such. In particular, a
predicted absence of toxicity is not evidence that a compound is non-toxic --
it is evidence that a model did not flag it.

Thresholds below are screening heuristics, not biological truth. The module
therefore emits a continuous score (`rules_passed_of_4`) alongside the boolean
gate so the manuscript can rank rather than merely filter.

SwissADME and ProTox-II have no official batch REST API, so a fully
reproducible open-science pipeline needs a documented, code-driven fallback:
this script computes the same descriptor-based rules SwissADME reports
(Lipinski, Ghose, Veber, Egan) directly with RDKit, plus a PAINS/Brenk
structural-alert screen.

For the two endpoints RDKit cannot approximate on its own (hERG liability,
hepatotoxicity), this script documents the manual step: submit
`data/developability_passed_compounds.csv` to SwissADME (http://www.swissadme.ch)
and ProTox-II (https://tox-new.charite.de/protox3/) and merge the downloaded
results with `merge_external_admet_results`. Record the access date and tool
version in the manuscript per the reproducibility checklist.
"""
from __future__ import annotations

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, Lipinski

RDLogger.DisableLog("rdApp.*")


def _rule_lipinski(mw, logp, hbd, hba) -> bool:
    return sum([mw > 500, logp > 5, hbd > 5, hba > 10]) <= 1


def _rule_veber(rotb, tpsa) -> bool:
    return rotb <= 10 and tpsa <= 140


def _rule_ghose(mw, logp, mr, n_atoms) -> bool:
    return (160 <= mw <= 480 and -0.4 <= logp <= 5.6
            and 40 <= mr <= 130 and 20 <= n_atoms <= 70)


def _rule_egan(logp, tpsa) -> bool:
    return tpsa <= 131.6 and logp <= 5.88


def compute_admet_descriptors(smiles: str, catalog: FilterCatalog.FilterCatalog) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "valid": False}

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    mr = Crippen.MolMR(mol)
    n_atoms = mol.GetNumHeavyAtoms()

    rules_passed = sum([
        _rule_lipinski(mw, logp, hbd, hba),
        _rule_veber(rotb, tpsa),
        _rule_ghose(mw, logp, mr, n_atoms),
        _rule_egan(logp, tpsa),
    ])

    return {
        "smiles": smiles,
        "valid": True,
        "mol_wt": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rotb,
        "tpsa": tpsa,
        "molar_refractivity": mr,
        "n_heavy_atoms": n_atoms,
        "rules_passed_of_4": rules_passed,      # continuous druglikeness score
        "structural_alert": catalog.HasMatch(mol),  # PAINS/Brenk flag
    }


def run_developability_prediction(input_csv: str = "data/derivative_library.csv",
                                   min_rules_passed: int = 3) -> pd.DataFrame:
    df_in = pd.read_csv(input_csv)
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)

    rows = [compute_admet_descriptors(smi, catalog) for smi in df_in["smiles"]]
    df = pd.DataFrame(rows)
    df["predicted_developability_pass"] = (
        (df["rules_passed_of_4"] >= min_rules_passed) & (~df["structural_alert"])
    )
    return df


def merge_external_admet_results(predictions_df: pd.DataFrame, swissadme_csv: str | None = None,
                                  protox_csv: str | None = None) -> pd.DataFrame:
    """Merge manually downloaded SwissADME / ProTox-II results by SMILES.

    Expected columns (rename to match your downloaded export headers):
      swissadme_csv: smiles, gi_absorption, bbb_permeant, cyp_inhibitor_flags, ...
      protox_csv:    smiles, predicted_ld50_mg_kg, tox_class, hepatotoxicity_prob, ...
    """
    merged = predictions_df
    if swissadme_csv:
        merged = merged.merge(pd.read_csv(swissadme_csv), on="smiles", how="left")
    if protox_csv:
        merged = merged.merge(pd.read_csv(protox_csv), on="smiles", how="left")
    return merged


if __name__ == "__main__":
    predictions = run_developability_prediction()
    predictions.to_csv("data/developability_predictions.csv", index=False)
    passed = predictions[predictions["predicted_developability_pass"]]
    passed[["smiles"]].to_csv("data/developability_passed_compounds.csv", index=False)
    print(f"{len(passed)}/{len(predictions)} compounds passed the descriptor-based "
          f"developability prediction. Submit data/developability_passed_compounds.csv to SwissADME "
          f"and ProTox-II for the two endpoints not covered here (hERG, "
          f"hepatotoxicity), then run merge_external_admet_results().")
