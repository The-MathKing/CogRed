"""Real execution of the seed -> derivative -> developability pipeline
(scripts/02, 03, 04 logic) on the 5 verified seed compounds, run locally
with RDKit -- no network access needed for this stage.
"""
import sys
sys.path.insert(0, "/home/user/research-molecule-editing/scripts")

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Crippen, Descriptors, FilterCatalog, Lipinski

RDLogger.DisableLog("rdApp.*")

SEED_COMPOUNDS = {
    "berzosertib_VE822":     ("ATR",  "CNCc1ccc(-c2cc(-c3nc(-c4ccc(S(=O)(=O)C(C)C)cc4)cnc3N)on2)cc1"),
    "camonsertib_RP3500":    ("ATR",  "C[C@@H]1COCCN1c1cc([C@]2(O)C[C@H]3CC[C@@H](C2)O3)c2cnn(C3=NN=CC3)c2n1"),
    "CHK1_inhibitor_YM8":    ("CHK1", "C[C@H](CN(C)C)OC1=N/C(=N\\C2=NCc3c(Cl)cccc3C2)C=N[C@H]1C#N"),
    "prexasertib":           ("CHK1", "COc1cccc(OCCCN)c1-c1cc(Nc2cnc(C#N)cn2)n[nH]1"),
    "ceralasertib":          ("ATR",  "C[C@@H]1COCCN1c1cc(C2([S@](C)(=N)=O)CC2)nc(-c2cncc3[nH]ccc23)n1"),
}
# Provenance note: berzosertib/camonsertib/YM8 extracted directly from the
# co-crystallized active-site ligand of PDB 9L40/9L4B/2YM8 respectively
# (this session, via Biopython + Open Babel bond-order perception, verified
# against RDKit-parsed SMILES). prexasertib/ceralasertib SMILES obtained via
# WebSearch and cross-checked against independently reported molecular
# formula/exact mass (both matched exactly) -- PubChem/ChEMBL APIs were not
# reachable from this sandboxed session.


def standardize(smi):
    from rdkit.Chem.MolStandardize import rdMolStandardize
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    lfc = rdMolStandardize.LargestFragmentChooser()
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(lfc.choose(mol))
    return Chem.MolToSmiles(mol)


def build_filter_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


def passes_druglike_filters(smiles, catalog):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_ok = rotb <= 10 and tpsa <= 140
    if violations > 1 or not veber_ok:
        return False
    if catalog.HasMatch(mol):
        return False
    return True


def _rule_lipinski(mw, logp, hbd, hba):
    return sum([mw > 500, logp > 5, hbd > 5, hba > 10]) <= 1

def _rule_veber(rotb, tpsa):
    return rotb <= 10 and tpsa <= 140

def _rule_ghose(mw, logp, mr, n_atoms):
    return (160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 40 <= mr <= 130 and 20 <= n_atoms <= 70)

def _rule_egan(logp, tpsa):
    return tpsa <= 131.6 and logp <= 5.88


def compute_developability(smiles, catalog):
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
        _rule_lipinski(mw, logp, hbd, hba), _rule_veber(rotb, tpsa),
        _rule_ghose(mw, logp, mr, n_atoms), _rule_egan(logp, tpsa),
    ])
    return {
        "smiles": smiles, "valid": True, "mol_wt": mw, "logp": logp,
        "hbd": hbd, "hba": hba, "rotatable_bonds": rotb, "tpsa": tpsa,
        "molar_refractivity": mr, "n_heavy_atoms": n_atoms,
        "rules_passed_of_4": rules_passed,
        "structural_alert": catalog.HasMatch(mol),
    }


if __name__ == "__main__":
    # 1. Standardize seeds
    seed_rows = []
    for name, (target, smi) in SEED_COMPOUNDS.items():
        std = standardize(smi)
        seed_rows.append({"name": name, "target": target, "smiles": smi, "standardized_smiles": std})
    seeds_df = pd.DataFrame(seed_rows)
    seeds_df.to_csv("seed_compounds.csv", index=False)
    print(f"Seed library: {len(seeds_df)} compounds")
    print(seeds_df[["name", "target", "standardized_smiles"]].to_string(index=False))

    # 2. BRICS fragment pool + recombination
    seed_smiles = seeds_df["standardized_smiles"].dropna().tolist()
    fragments = set()
    for smi in seed_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fragments |= BRICS.BRICSDecompose(mol)
    print(f"\nBRICS fragment pool: {len(fragments)} unique fragments")

    # Sort for reproducibility -- Python set iteration order is not stable
    # across interpreter runs, which made BRICSBuild's enumeration (and thus
    # the first 5000 products under the cap) non-deterministic run-to-run.
    frag_mols = [Chem.MolFromSmiles(f) for f in sorted(fragments)]
    frag_mols = [m for m in frag_mols if m is not None]
    products = []
    for i, mol in enumerate(BRICS.BRICSBuild(frag_mols)):
        if i >= 5000:
            break
        try:
            Chem.SanitizeMol(mol)
            products.append(Chem.MolToSmiles(mol))
        except Exception:
            continue
    products = sorted(set(products))
    print(f"BRICS recombination candidates (pre-filter): {len(products)}")

    # 3. Pre-filter (Lipinski/Veber/PAINS/Brenk)
    catalog = build_filter_catalog()
    kept = [smi for smi in products if passes_druglike_filters(smi, catalog)]
    print(f"Candidates passing Lipinski/Veber/PAINS/Brenk pre-filter: {len(kept)}")

    derivative_rows = []
    for smi in kept:
        mol = Chem.MolFromSmiles(smi)
        derivative_rows.append({
            "smiles": smi, "mol_wt": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol), "tpsa": Descriptors.TPSA(mol),
        })
    derivatives_df = pd.DataFrame(derivative_rows).drop_duplicates("smiles").reset_index(drop=True)
    derivatives_df.to_csv("derivative_library.csv", index=False)

    # 4. Developability prediction on seeds + derivatives combined
    all_compounds = pd.concat([
        seeds_df[["standardized_smiles"]].rename(columns={"standardized_smiles": "smiles"}).assign(source="seed"),
        derivatives_df[["smiles"]].assign(source="derivative"),
    ], ignore_index=True)

    dev_rows = [compute_developability(smi, catalog) for smi in all_compounds["smiles"]]
    dev_df = pd.DataFrame(dev_rows)
    dev_df["source"] = all_compounds["source"].values
    dev_df["predicted_developability_pass"] = (
        (dev_df["rules_passed_of_4"] >= 3) & (~dev_df["structural_alert"])
    )
    dev_df.to_csv("developability_predictions.csv", index=False)

    passed = dev_df[dev_df["predicted_developability_pass"]]
    passed[["smiles", "source"]].to_csv("developability_passed_compounds.csv", index=False)

    print(f"\nDevelopability predictions: {len(passed)}/{len(dev_df)} compounds passed "
          f"(>=3/4 rules, no PAINS/Brenk alert)")
    print(dev_df.groupby("source")["predicted_developability_pass"].agg(["sum", "count"]))
