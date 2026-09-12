"""
Generate structural derivatives ("optimization" candidates) from the seed
compound library using RDKit: BRICS fragmentation/recombination and simple
R-group enumeration, followed by a cheap pre-filter (Lipinski/Veber + PAINS/
Brenk) so expensive docking is only run on plausible drug-like structures.

This corresponds to manuscript Section 2.3.
"""
from __future__ import annotations

import itertools

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import BRICS, Crippen, Descriptors, FilterCatalog, Lipinski

RDLogger.DisableLog("rdApp.*")


def build_brics_fragment_pool(seed_smiles: list[str]) -> set[str]:
    """Break each seed molecule into BRICS fragments (a fragment vocabulary)."""
    fragments: set[str] = set()
    for smi in seed_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fragments |= BRICS.BRICSDecompose(mol)
    return fragments


def enumerate_brics_derivatives(fragments: set[str], max_products: int = 5000) -> list[str]:
    """Recombine a fragment pool with BRICS.BRICSBuild to propose new molecules."""
    frag_mols = [Chem.MolFromSmiles(f) for f in fragments]
    frag_mols = [m for m in frag_mols if m is not None]
    products = []
    for i, mol in enumerate(BRICS.BRICSBuild(frag_mols)):
        if i >= max_products:
            break
        try:
            Chem.SanitizeMol(mol)
            products.append(Chem.MolToSmiles(mol))
        except Exception:
            continue
    return sorted(set(products))


def passes_druglike_filters(smiles: str, catalog: FilterCatalog.FilterCatalog) -> bool:
    """Lipinski Ro5 (relaxed: allow <=1 violation) + Veber + PAINS/Brenk screen."""
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


def build_filter_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


def generate_derivative_library(seed_csv: str = "data/seed_named_compounds.csv",
                                 max_products: int = 5000) -> pd.DataFrame:
    seeds = pd.read_csv(seed_csv)
    seed_smiles = seeds["standardized_smiles"].dropna().tolist()

    fragments = build_brics_fragment_pool(seed_smiles)
    candidates = enumerate_brics_derivatives(fragments, max_products=max_products)

    catalog = build_filter_catalog()
    kept = [smi for smi in candidates if passes_druglike_filters(smi, catalog)]

    df = pd.DataFrame({"smiles": kept})
    df["mol_wt"] = df["smiles"].apply(lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s)))
    df["logp"] = df["smiles"].apply(lambda s: Crippen.MolLogP(Chem.MolFromSmiles(s)))
    df["tpsa"] = df["smiles"].apply(lambda s: Descriptors.TPSA(Chem.MolFromSmiles(s)))
    return df.drop_duplicates("smiles").reset_index(drop=True)


if __name__ == "__main__":
    derivatives = generate_derivative_library()
    derivatives.to_csv("data/derivative_library.csv", index=False)
    print(f"Generated {len(derivatives)} drug-like derivative candidates "
          f"after BRICS recombination + PAINS/Brenk/Lipinski/Veber filtering.")
