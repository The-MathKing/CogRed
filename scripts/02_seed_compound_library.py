"""
Build a reproducible seed compound library for ATR/CHK1 from PubChem, rather
than hand-typing SMILES (a common, hard-to-catch source of error in CADD
papers).

Strategy: resolve a curated list of known clinical/tool ATR & CHK1
modulators by name via the PubChem PUG-REST API, retrieve canonical SMILES,
CID, and molecular formula, then standardize with RDKit. This list is the
*named-compound* seed set; Section 2.2 of the manuscript outline also calls
for a broader ChEMBL bioactivity pull (see `fetch_chembl_actives` below) to
avoid seed-set bias toward only well-known clinical compounds.
"""
from __future__ import annotations

import time

import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

PUBCHEM_API = "https://pubchem.ncbi.gov/rest/pug"

# Known ATR/CHK1-pathway modulators (all inhibitors in current clinical/tool
# use -- see critique doc Section 2 on why this is the *starting* chemical
# space, not necessarily the target pharmacology).
SEED_COMPOUND_NAMES = {
    "ATR": ["berzosertib", "ceralasertib", "elimusertib", "gartisertib"],
    "CHK1": ["prexasertib", "SRA737", "rabusertib", "MK-8776"],
}


def fetch_pubchem_record(name: str) -> dict | None:
    """Resolve a compound name to CID + canonical SMILES via PubChem PUG-REST."""
    url = f"{PUBCHEM_API}/compound/name/{requests.utils.quote(name)}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula/JSON"
    resp = requests.get(url, timeout=30)
    if resp.status_code != 200:
        print(f"  [warn] PubChem lookup failed for '{name}' (HTTP {resp.status_code})")
        return None
    props = resp.json()["PropertyTable"]["Properties"][0]
    return {
        "name": name,
        "pubchem_cid": props["CID"],
        "smiles": props.get("IsomericSMILES", props["CanonicalSMILES"]),
        "formula": props["MolecularFormula"],
    }


def standardize_smiles(smiles: str) -> str | None:
    """RDKit standardization: strip salts, neutralize charges, canonicalize."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    lfc = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(lfc.choose(mol))
    return Chem.MolToSmiles(mol)


def build_named_seed_table(request_delay_s: float = 0.5) -> pd.DataFrame:
    rows = []
    for target, names in SEED_COMPOUND_NAMES.items():
        for name in names:
            rec = fetch_pubchem_record(name)
            time.sleep(request_delay_s)
            if rec is None:
                continue
            rec["target"] = target
            rec["standardized_smiles"] = standardize_smiles(rec["smiles"])
            rows.append(rec)
    return pd.DataFrame(rows)


def fetch_chembl_actives(target_chembl_id: str, ic50_nm_max: float = 1000.0,
                          limit: int = 1000) -> pd.DataFrame:
    """Pull bioactivity-annotated actives for a target from the ChEMBL REST API.

    target_chembl_id examples: ATR = CHEMBL2842, CHK1 = CHEMBL2996 (verify
    current IDs at https://www.ebi.ac.uk/chembl/ before running -- ChEMBL
    target IDs are stable but should be confirmed for the ChEMBL release used).
    """
    url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    params = {
        "target_chembl_id": target_chembl_id,
        "standard_type": "IC50",
        "standard_units": "nM",
        "standard_value__lte": ic50_nm_max,
        "limit": limit,
    }
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    activities = resp.json()["activities"]
    rows = [{
        "molecule_chembl_id": a.get("molecule_chembl_id"),
        "canonical_smiles": a.get("canonical_smiles"),
        "standard_value_nm": a.get("standard_value"),
        "assay_chembl_id": a.get("assay_chembl_id"),
    } for a in activities if a.get("canonical_smiles")]
    df = pd.DataFrame(rows).drop_duplicates("molecule_chembl_id")
    df["standardized_smiles"] = df["canonical_smiles"].apply(standardize_smiles)
    return df.dropna(subset=["standardized_smiles"])


if __name__ == "__main__":
    named_df = build_named_seed_table()
    named_df.to_csv("data/seed_named_compounds.csv", index=False)
    print(named_df[["name", "target", "pubchem_cid", "standardized_smiles"]])

    # Example broader pull -- confirm target_chembl_id before uncommenting.
    # atr_actives = fetch_chembl_actives("CHEMBL2842")
    # chk1_actives = fetch_chembl_actives("CHEMBL2996")
    # pd.concat([atr_actives.assign(target="ATR"),
    #            chk1_actives.assign(target="CHK1")]).to_csv(
    #     "data/chembl_actives.csv", index=False)
