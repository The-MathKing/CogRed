"""
Genomic mapping and feature extraction for common fragile sites (CFS).

Pulls gene coordinates for the classic CFS host genes (FHIT @ FRA3B, WWOX @
FRA16D, plus FRA7H/FRA1H as a comparison set) from the UCSC REST API, fetches
the underlying genomic sequence, and computes simple structural-instability
proxies (AT content, gene span, GC-skew) that ground the "fragile site"
narrative in actual genomics rather than treating it as a black box.

This is a supporting genomics analysis (manuscript Section 2.10 / 3.7) — it
does not feed directly into docking, but it is what lets the paper claim the
compounds are relevant to a well-defined genomic phenomenon.

Data source: UCSC Genome Browser REST API (https://api.genome.ucsc.edu),
hg38 assembly. No API key required.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import pandas as pd
import requests

UCSC_API = "https://api.genome.ucsc.edu"
ASSEMBLY = "hg38"

# Coordinates below are RefSeq gene spans (hg38), sourced from the UCSC
# Genome Browser refGene track. Verify against the current track before
# publication -- gene models are periodically revised.
FRAGILE_SITE_GENES = {
    # gene_symbol: (chrom, start, end, fragile_site_name)
    "FHIT":  ("chr3",  59_747_277,  61_251_452, "FRA3B"),
    "WWOX":  ("chr16", 78_073_305,  79_246_564, "FRA16D"),
    "CNTNAP2": ("chr7", 145_813_453, 148_120_022, "FRA7H"),
    "DMD":   ("chrX",  31_097_677,  33_339_609, "FRAXC (comparison, large gene)"),
}


@dataclass
class FragileSiteFeatures:
    gene: str
    fragile_site: str
    chrom: str
    start: int
    end: int
    span_bp: int
    at_content: float
    gc_skew: float


def fetch_sequence(chrom: str, start: int, end: int, assembly: str = ASSEMBLY) -> str:
    """Fetch genomic DNA sequence for a region via the UCSC REST API.

    UCSC coordinates in the DAS/REST API are 0-based half-open; refGene
    coordinates from the track table are already in that convention.
    """
    url = f"{UCSC_API}/getData/sequence"
    params = {"genome": assembly, "chrom": chrom, "start": start, "end": end}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    payload = resp.json()
    return payload["dna"].upper()


def compute_features(gene: str, chrom: str, start: int, end: int, fragile_site: str,
                      sequence: str) -> FragileSiteFeatures:
    span = end - start
    a, t, g, c = (sequence.count(b) for b in "ATGC")
    total = a + t + g + c
    at_content = (a + t) / total if total else float("nan")
    gc_skew = (g - c) / (g + c) if (g + c) else float("nan")
    return FragileSiteFeatures(gene, fragile_site, chrom, start, end, span,
                                at_content, gc_skew)


def build_feature_table(genes: dict[str, tuple[str, int, int, str]] = FRAGILE_SITE_GENES,
                         request_delay_s: float = 1.0) -> pd.DataFrame:
    rows = []
    for gene, (chrom, start, end, site) in genes.items():
        seq = fetch_sequence(chrom, start, end)
        feat = compute_features(gene, chrom, start, end, site, seq)
        rows.append(feat.__dict__)
        time.sleep(request_delay_s)  # be polite to the public API
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = build_feature_table()
    df.to_csv("data/fragile_site_features.csv", index=False)
    print(df.to_string(index=False))
