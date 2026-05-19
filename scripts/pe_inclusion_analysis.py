#!/usr/bin/env python3
"""
PE differential inclusion analysis driven by a long-read POISEN PE list.

Saba's spec: take Talha's POISEN poison-exon list as the authoritative event
catalog, map every short-read sample's STAR junction file to those PE coordinates,
and run differential inclusion across developmental timepoints.

This module is self-contained.  It reads a POISEN TSV (UPDATED_STRAND_POISEN_HITS
schema), parses each row's per-transcript exon chain to derive the three relevant
introns per PE (incl_donor, incl_acceptor, skip), looks up read counts in a per-sample
collection of STAR SJ.out.tab files, applies LC2-equivalent statistical hygiene,
computes per-sample inclusion PSI, then runs Spearman correlation of PSI vs hpf
across the developmental timeline.  Optional Dirichlet-multinomial confirmation is
gated behind ``--ds``.

Outputs:
  - pe_events.tsv               : canonical PE event catalog (one row per unique event)
  - pe_inclusion_psi.tsv.gz     : events x samples PSI matrix
  - pe_inclusion_counts.tsv.gz  : events x samples raw counts (incl_donor/incl_acceptor/skip)
  - pe_low_coverage.tsv         : audit trail of dropped (event, sample) pairs
  - pe_developmental_candidates.tsv : ranked candidate table (headline deliverable)
  - pe_developmental_candidates.json
  - summary.json
  - pe_psi_heatmap.png
  - pe_psi_vs_time_top.png
  - pe_spearman_distribution.png
  - pe_lc2_concordance.png      : when LC2 cross-reference is available
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Constants and small helpers
# ---------------------------------------------------------------------------

# Literature-standard hygiene thresholds for splicing differential inclusion.
# References:
#  - LeafCutter2 (Yang lab, 2025 bioRxiv): >=10 reads/junction for cross-tissue,
#    >=50 reads as "appreciable abundance" for developmental analysis,
#    intron-cluster usage in >=60% of samples for cross-condition tests.
#  - Mudge & Pritchard, Nat Genet 2024 ("Global impact of unproductive splicing"):
#    aggregate counts per condition, |dPSI| >= 0.10 effect-size bar.
#  - Leclair/Lareau, Mol Cell 2020 (canonical PE differentiation paper):
#    pre-filter PEs whose host-gene transcript is too lowly expressed to
#    accurately quantify (proxied here by total event reads across all samples).
DEFAULT_MIN_JUNCTION_READS = 3       # mirrors LC2 pool_junc_reads threshold
DEFAULT_MIN_LOCUS_READS = 10         # per-PE per-sample inclusion+skip floor
DEFAULT_MIN_OBSERVABILITY_FRAC = 0.60  # event must be observed in >=60% of dev samples
DEFAULT_MIN_TIMEPOINTS_OBSERVED = 3  # event must have data in >=3 distinct timepoints
DEFAULT_MIN_TOTAL_LOCUS_READS = 50   # LeafCutter2 "appreciable abundance" floor
DEFAULT_MIN_DPSI = 0.10              # biological effect-size bar
DEFAULT_MAX_INTRON_LEN = 500_000     # mirrors LC2 maxIntronLen default
DEFAULT_FDR = 0.05

# Developmental stage parsing — same regex idea as scripts/lc2_pipeline.py
_STAGE_PATTERN = re.compile(
    r"\d+(hpf|dpf|cell|somite)"
    r"|dome|shield|bud|sphere|oblong|gastrula|pharyngula|segmentation"
    r"|blastula|cleavage|epiboly|hatching|larval",
    re.IGNORECASE,
)

# Match parenthesized exon coord like "(NC_133177.1:10847256-10847464)"
_EXON_TOKEN = re.compile(r"\(?([^():\s]+):(\d+)-(\d+)\)?")


def _is_developmental_stage_label(label: str) -> bool:
    return bool(_STAGE_PATTERN.fullmatch(str(label).strip()))


def _timepoint_to_hours(label: str) -> float:
    """Convert a developmental-stage label to hours.  Kept identical in behaviour
    to scripts/lc2_pipeline.py:_timepoint_to_hours so cross-reference columns
    line up exactly with the existing LC2 candidate table."""
    s = str(label).strip().lower()
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*(hpf|dpf)", s)
    if m:
        val = float(m.group(1))
        return val if m.group(2) == "hpf" else val * 24.0
    return float("inf")


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> Tuple[Optional[float], Optional[float]]:
    """Spearman rho with two-tailed p-value, with proper tie handling.

    Uses average ranks for tied values (the same convention as
    ``scipy.stats.spearmanr``). When either input has zero variance after
    ranking (i.e. all values tied) the correlation is undefined and we return
    ``(None, None)`` instead of a spurious correlation. This is critical for
    PSI series that saturate at 0 or 1 across all timepoints — without proper
    tie handling, naive argsort-based ranks make a constant series perfectly
    correlated with the time axis.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 3:
        return None, None

    def _avg_rankdata(a: np.ndarray) -> np.ndarray:
        order = np.argsort(a, kind="mergesort")
        ranks = np.empty(n, dtype=float)
        i = 0
        while i < n:
            j = i + 1
            while j < n and a[order[j]] == a[order[i]]:
                j += 1
            avg = 0.5 * ((i + 1) + j)  # average of 1-based ranks i+1..j
            for k in range(i, j):
                ranks[order[k]] = avg
            i = j
        return ranks

    rx = _avg_rankdata(x)
    ry = _avg_rankdata(y)
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        return None, None

    mx, my = rx.mean(), ry.mean()
    dx, dy = rx - mx, ry - my
    denom = math.sqrt(float((dx ** 2).sum()) * float((dy ** 2).sum()))
    if denom == 0:
        return None, None
    rho = float((dx * dy).sum() / denom)
    rho = max(min(rho, 1.0), -1.0)
    if abs(rho) >= 1.0 - 1e-15:
        # Perfect correlation: t-stat blows up, p effectively 0. Cap at a
        # finite small value so downstream FDR maths stays sane.
        return rho, 1e-300
    t_stat = rho * math.sqrt((n - 2) / (1.0 - rho ** 2))
    df = n - 2
    x_beta = df / (df + t_stat ** 2)
    p = _betai(df / 2.0, 0.5, x_beta)
    return rho, p


def _betai(a: float, b: float, x: float) -> float:
    if x < 0 or x > 1:
        return 1.0
    if x == 0 or x == 1:
        return 1.0 - x
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1.0 - x) - lbeta) / a
    if x < (a + 1) / (a + b + 2):
        return front * _betacf(a, b, x)
    return 1.0 - math.exp(a * math.log(x) + b * math.log(1.0 - x) - lbeta) / b * _betacf(b, a, 1.0 - x)


def _betacf(a: float, b: float, x: float, max_iter: int = 200, eps: float = 3e-12) -> float:
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = max(1.0 - qab * x / qap, 1e-30)
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = max(1.0 + aa * d, 1e-30)
        c = max(1.0 + aa / c, 1e-30)
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = max(1.0 + aa * d, 1e-30)
        c = max(1.0 + aa / c, 1e-30)
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _t_two_tailed_pvalue(t_stat: float, df: float) -> float:
    """Two-tailed p-value for a t-statistic with ``df`` residual degrees of freedom.

    Uses the same regularized-incomplete-beta routine that backs ``_spearman_rho``.
    """
    if not math.isfinite(t_stat) or df <= 0:
        return 1.0
    x = df / (df + t_stat * t_stat)
    return _betai(df / 2.0, 0.5, x)


def _quasi_binomial_glm(
    x: np.ndarray, k: np.ndarray, n: np.ndarray, *, max_iter: int = 60, tol: float = 1e-8
) -> Optional[Dict]:
    """Fit ``logit(p_i) = b0 + b1 * x_i`` to (k_i, n_i) by IRLS, treating the
    response as binomial in the working step but inflating the standard errors
    by a Pearson-chi-squared overdispersion estimate (quasi-binomial).

    This is exactly the test used by ``glm(family = quasibinomial)`` in R and
    the count-based reference test for proportional splicing data (see
    LeafCutter2 supplementary methods, Mudge & Pritchard 2024 Methods).

    Returns ``None`` when the fit is not informative (degenerate response,
    insufficient degrees of freedom, no variation in the predictor, or singular
    Fisher information).
    """
    x = np.asarray(x, dtype=float)
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    if x.size != k.size or x.size != n.size:
        return None

    keep = np.isfinite(x) & np.isfinite(k) & np.isfinite(n) & (n > 0)
    if keep.sum() < 4:
        return None
    x = x[keep]
    k = k[keep]
    n = n[keep]

    if np.ptp(x) == 0:
        return None
    k_sum = float(k.sum())
    n_sum = float(n.sum())
    if k_sum <= 0 or k_sum >= n_sum:
        return None  # response is degenerate (always 0 or always 1) -> no slope info

    n_obs = x.size
    df_resid = n_obs - 2
    if df_resid <= 0:
        return None

    X = np.column_stack([np.ones(n_obs, dtype=float), x])

    # Initialise on the empirical proportions (Laplace-smoothed to keep eta finite).
    p = np.clip((k + 0.5) / (n + 1.0), 1e-6, 1.0 - 1e-6)
    eta = np.log(p / (1.0 - p))
    beta = np.array([eta.mean(), 0.0], dtype=float)

    for _ in range(max_iter):
        eta = X @ beta
        # Numerically stable sigmoid.
        mu = np.where(eta >= 0,
                      1.0 / (1.0 + np.exp(-eta)),
                      np.exp(eta) / (1.0 + np.exp(eta)))
        mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
        w = n * mu * (1.0 - mu)         # Fisher weight for binomial logit
        good = w > 1e-12
        if good.sum() < 3:
            return None
        # Working response: z = eta + (k - n*mu) / (n*mu*(1-mu))
        z = eta + (k - n * mu) / (n * mu * (1.0 - mu))
        sw = np.sqrt(w[good])
        Xw = X[good] * sw[:, None]
        zw = z[good] * sw
        try:
            beta_new, *_ = np.linalg.lstsq(Xw, zw, rcond=None)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(beta_new)):
            return None
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    else:
        # IRLS did not converge; treat as uninformative.
        return None

    # Final fitted values for SE / overdispersion.
    eta = X @ beta
    mu = np.where(eta >= 0,
                  1.0 / (1.0 + np.exp(-eta)),
                  np.exp(eta) / (1.0 + np.exp(eta)))
    mu = np.clip(mu, 1e-9, 1.0 - 1e-9)
    w = n * mu * (1.0 - mu)
    good = w > 1e-12
    if good.sum() < 3:
        return None
    XtWX = (X[good] * w[good][:, None]).T @ X[good]
    try:
        cov_unscaled = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return None

    # Pearson chi-squared / df_resid, floored at 1.0 (quasi-binomial convention).
    pearson = ((k[good] - n[good] * mu[good]) ** 2
               / (n[good] * mu[good] * (1.0 - mu[good])))
    phi = max(float(pearson.sum()) / float(df_resid), 1.0)

    se_b1 = math.sqrt(cov_unscaled[1, 1] * phi)
    if not math.isfinite(se_b1) or se_b1 <= 0:
        return None
    t_stat = float(beta[1]) / se_b1
    p_val = _t_two_tailed_pvalue(t_stat, df_resid)
    return {
        "beta0": float(beta[0]),
        "beta1": float(beta[1]),
        "se_beta1": float(se_b1),
        "t_stat": float(t_stat),
        "phi": float(phi),
        "n_obs": int(n_obs),
        "df_resid": int(df_resid),
        "p_value": float(p_val),
    }


def _bh_fdr(pvals: List[Optional[float]]) -> List[Optional[float]]:
    """Benjamini-Hochberg FDR.  Returns ``None`` for missing inputs in place."""
    idx = [i for i, p in enumerate(pvals) if p is not None and not math.isnan(p)]
    if not idx:
        return [None] * len(pvals)
    p_arr = np.array([pvals[i] for i in idx], dtype=float)
    m = len(p_arr)
    order = np.argsort(p_arr)
    ranks = np.empty(m, dtype=int)
    ranks[order] = np.arange(1, m + 1)
    adj = p_arr * m / ranks
    # Enforce monotonicity from largest to smallest p-value
    sorted_adj = adj[order]
    for j in range(m - 2, -1, -1):
        sorted_adj[j] = min(sorted_adj[j], sorted_adj[j + 1])
    monotone = np.empty(m, dtype=float)
    monotone[order] = np.minimum(sorted_adj, 1.0)
    out: List[Optional[float]] = [None] * len(pvals)
    for k, i in enumerate(idx):
        out[i] = float(monotone[k])
    return out


# ---------------------------------------------------------------------------
# Stage A — parse POISEN events, derive flanking introns
# ---------------------------------------------------------------------------


def _parse_genomic_coordinates(field: str) -> List[Tuple[str, int, int]]:
    """Parse the POISEN ``genomic_coordinates`` column.

    Format: semicolon-separated, parenthesised entries like
    ``(NC_133177.1:10847256-10847464);(NC_133177.1:10844621-10844750);...``.
    Entries are listed in transcript order (5' -> 3'), which means for + strand
    they are sorted ascending by genomic position and for - strand descending.
    """
    if not field:
        return []
    out: List[Tuple[str, int, int]] = []
    for tok in field.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        m = _EXON_TOKEN.search(tok)
        if not m:
            continue
        chrom = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        if start > end:
            start, end = end, start
        out.append((chrom, start, end))
    return out


def _parse_pe_coord(field: str, fallback_chrom: str = "") -> Optional[Tuple[str, int, int]]:
    """Parse a POISEN ``poison_exon_coordinate`` value (``NC_xxx.1:start-end``)."""
    if not field:
        return None
    m = _EXON_TOKEN.search(field)
    if not m:
        return None
    chrom = m.group(1) or fallback_chrom
    start = int(m.group(2))
    end = int(m.group(3))
    if start > end:
        start, end = end, start
    return chrom, start, end


def _derive_flanking_introns(
    exons_in_tx_order: List[Tuple[str, int, int]],
    pe_chrom: str,
    pe_start: int,
    pe_end: int,
) -> Optional[Dict[str, Tuple[str, int, int]]]:
    """Identify the PE within the transcript exon chain and derive the three
    relevant introns (always returned in genomic-low/high orientation, ready for
    STAR ``SJ.out.tab`` lookup).

    Returns ``None`` when:
      - the PE coordinate cannot be matched to any exon in the chain,
      - the PE is the first or last exon of the transcript (no skip junction),
      - the chrom mismatches the rest of the transcript.
    """
    if not exons_in_tx_order:
        return None

    # Locate PE in the chain.  Allow exact match first; then near-match (within 1nt)
    # to tolerate trivial 0/1-based off-by-one drift in the source TSV.
    pe_idx = -1
    for i, (chrom, s, e) in enumerate(exons_in_tx_order):
        if chrom != pe_chrom:
            continue
        if s == pe_start and e == pe_end:
            pe_idx = i
            break
    if pe_idx < 0:
        for i, (chrom, s, e) in enumerate(exons_in_tx_order):
            if chrom == pe_chrom and abs(s - pe_start) <= 1 and abs(e - pe_end) <= 1:
                pe_idx = i
                break
    if pe_idx < 0:
        return None

    if pe_idx == 0 or pe_idx == len(exons_in_tx_order) - 1:
        return None

    prev_exon = exons_in_tx_order[pe_idx - 1]   # upstream in transcript order
    next_exon = exons_in_tx_order[pe_idx + 1]   # downstream in transcript order
    if prev_exon[0] != pe_chrom or next_exon[0] != pe_chrom:
        return None

    # Identify the genomically-flanking exons (regardless of strand).
    if prev_exon[1] < pe_start:
        upstream_genomic = prev_exon
        downstream_genomic = next_exon
    else:
        upstream_genomic = next_exon
        downstream_genomic = prev_exon

    if not (upstream_genomic[2] <= pe_start < pe_end <= downstream_genomic[1]):
        return None

    # POISEN uses 0-based half-open BED-style coords (verified: pe_size == pe_end - pe_start).
    # STAR SJ.out.tab uses 1-based inclusive intron boundaries (first base, last base of intron).
    # Conversion: BED end == 1-based last base of the exon. So:
    #   upstream_intron 1-based: [u_e + 1, pe_start]
    #   downstream_intron 1-based: [pe_end + 1, d_s]
    #   skip_intron 1-based:     [u_e + 1, d_s]
    intron_5p_genomic = (pe_chrom, upstream_genomic[2] + 1, pe_start)
    intron_3p_genomic = (pe_chrom, pe_end + 1, downstream_genomic[1])
    intron_skip = (pe_chrom, upstream_genomic[2] + 1, downstream_genomic[1])

    return {
        "intron_5p_genomic": intron_5p_genomic,
        "intron_3p_genomic": intron_3p_genomic,
        "intron_skip": intron_skip,
        "upstream_exon": upstream_genomic,
        "downstream_exon": downstream_genomic,
        "pe_idx_in_tx": pe_idx,
        "n_exons_in_tx": len(exons_in_tx_order),
    }


def parse_pe_events(poisen_tsv: Path) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Stage A: read POISEN TSV, derive flanking introns, deduplicate to unique events.

    Returns the events DataFrame and a stats dict for summary.json.
    """
    if not poisen_tsv.exists():
        raise FileNotFoundError(f"POISEN PE list not found: {poisen_tsv}")

    df = pd.read_csv(poisen_tsv, sep="\t", dtype=str, keep_default_na=False)

    required = {
        "read_id", "strand", "genomic_coordinates",
        "pe_chrom", "pe_start", "pe_end",
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"POISEN TSV missing required columns: {sorted(missing)}")

    # Per-row PE event derivation
    derived_rows: List[Dict] = []
    n_total = len(df)
    n_no_chain = 0
    n_pe_not_found = 0
    n_terminal_exon = 0
    n_chrom_mismatch = 0

    for _, row in df.iterrows():
        chain = _parse_genomic_coordinates(row.get("genomic_coordinates", ""))
        if not chain:
            n_no_chain += 1
            continue

        pe_chrom = (row.get("pe_chrom") or "").strip()
        try:
            pe_start = int(row["pe_start"])
            pe_end = int(row["pe_end"])
        except (KeyError, ValueError, TypeError):
            n_pe_not_found += 1
            continue
        if pe_start > pe_end:
            pe_start, pe_end = pe_end, pe_start

        if not pe_chrom:
            parsed_pe = _parse_pe_coord(row.get("poison_exon_coordinate", ""))
            if parsed_pe is None:
                n_pe_not_found += 1
                continue
            pe_chrom, pe_start, pe_end = parsed_pe

        introns = _derive_flanking_introns(chain, pe_chrom, pe_start, pe_end)
        if introns is None:
            # Distinguish the failure mode for accounting
            if not any(c == pe_chrom for (c, _, _) in chain):
                n_chrom_mismatch += 1
                continue
            # Try to detect whether PE was first/last
            pe_pos = -1
            for i, (c, s, e) in enumerate(chain):
                if c == pe_chrom and s == pe_start and e == pe_end:
                    pe_pos = i
                    break
            if pe_pos == 0 or pe_pos == len(chain) - 1:
                n_terminal_exon += 1
            else:
                n_pe_not_found += 1
            continue

        gene_symbol = (
            (row.get("gene_symbol_annotated") or "").strip()
            or (row.get("gene_symbol") or "").strip()
        )
        gene_id = (
            (row.get("gene_id_annotated") or "").strip()
            or (row.get("gene_id") or "").strip()
        )

        derived_rows.append({
            "read_id": (row.get("read_id") or "").strip(),
            "gene_symbol": gene_symbol,
            "gene_id": gene_id,
            "chrom": pe_chrom,
            "strand": (row.get("strand") or "").strip(),
            "pe_start": pe_start,
            "pe_end": pe_end,
            "pe_size": pe_end - pe_start + 1,
            "poison_exon_type": (row.get("poison_exon_type") or "").strip(),
            "coord_match_ptc": (row.get("coord_match_ptc") or "").strip(),
            "frame_restored": (row.get("frame_restored") or "").strip(),
            "ptc_exon": (row.get("ptc_exon") or "").strip(),
            "intron_5p_genomic_start": introns["intron_5p_genomic"][1],
            "intron_5p_genomic_end": introns["intron_5p_genomic"][2],
            "intron_3p_genomic_start": introns["intron_3p_genomic"][1],
            "intron_3p_genomic_end": introns["intron_3p_genomic"][2],
            "intron_skip_start": introns["intron_skip"][1],
            "intron_skip_end": introns["intron_skip"][2],
            "upstream_exon_end": introns["upstream_exon"][2],
            "downstream_exon_start": introns["downstream_exon"][1],
            "pe_idx_in_tx": introns["pe_idx_in_tx"],
            "n_exons_in_tx": introns["n_exons_in_tx"],
        })

    if not derived_rows:
        raise SystemExit("[STAGE A] No usable PE events derived from POISEN TSV")

    derived = pd.DataFrame(derived_rows)

    def _mode_or_first(s: pd.Series) -> str:
        s = s.dropna().astype(str)
        s = s[s != ""]
        if s.empty:
            return ""
        c = Counter(s)
        return c.most_common(1)[0][0]

    grouped = derived.groupby(
        ["chrom", "strand", "pe_start", "pe_end"], sort=False, dropna=False
    )

    events: List[Dict] = []
    for (chrom, strand, pe_start, pe_end), g in grouped:
        gene_symbol = _mode_or_first(g["gene_symbol"])
        gene_id = _mode_or_first(g["gene_id"])
        pe_type = _mode_or_first(g["poison_exon_type"])
        coord_match_ptc = _mode_or_first(g["coord_match_ptc"])
        # Use the modal flanking-intron / skip coords (rare disagreements show up
        # for shared PEs whose donor/acceptor coordinates differ slightly across
        # supporting transcripts; a vote is the safest way to lock one set in)
        i5p_s = int(g["intron_5p_genomic_start"].mode().iloc[0])
        i5p_e = int(g["intron_5p_genomic_end"].mode().iloc[0])
        i3p_s = int(g["intron_3p_genomic_start"].mode().iloc[0])
        i3p_e = int(g["intron_3p_genomic_end"].mode().iloc[0])
        sk_s = int(g["intron_skip_start"].mode().iloc[0])
        sk_e = int(g["intron_skip_end"].mode().iloc[0])

        # Skip events whose flanking introns exceed maxIntronLen — these are
        # either gene-fusion-like artifacts or annotations that LC2's pool_junc_reads
        # would also discard.  We mark them but still emit so the audit trail is intact.
        max_intron_len = max(
            i5p_e - i5p_s + 1,
            i3p_e - i3p_s + 1,
            sk_e - sk_s + 1,
        )

        event_id = f"{chrom}:{pe_start}-{pe_end}:{strand or '?'}"
        events.append({
            "event_id": event_id,
            "chrom": chrom,
            "strand": strand,
            "pe_start": int(pe_start),
            "pe_end": int(pe_end),
            "pe_size": int(pe_end) - int(pe_start) + 1,
            "gene_symbol": gene_symbol,
            "gene_id": gene_id,
            "poison_exon_type": pe_type,
            "coord_match_ptc": coord_match_ptc,
            "intron_5p_genomic_start": i5p_s,
            "intron_5p_genomic_end": i5p_e,
            "intron_3p_genomic_start": i3p_s,
            "intron_3p_genomic_end": i3p_e,
            "intron_skip_start": sk_s,
            "intron_skip_end": sk_e,
            "upstream_exon_end": int(g["upstream_exon_end"].mode().iloc[0]),
            "downstream_exon_start": int(g["downstream_exon_start"].mode().iloc[0]),
            "max_intron_len": int(max_intron_len),
            "n_supporting_long_read_tx": int(len(g)),
            "supporting_read_ids": ",".join(sorted(set(g["read_id"]))[:64]),
        })

    events_df = pd.DataFrame(events).sort_values(
        ["chrom", "pe_start", "pe_end"]
    ).reset_index(drop=True)

    stats = {
        "n_input_rows": int(n_total),
        "n_unique_events": int(len(events_df)),
        "n_dropped_no_chain": int(n_no_chain),
        "n_dropped_pe_not_in_chain": int(n_pe_not_found),
        "n_dropped_terminal_exon": int(n_terminal_exon),
        "n_dropped_chrom_mismatch": int(n_chrom_mismatch),
    }
    return events_df, stats


# ---------------------------------------------------------------------------
# Stage B — per-sample PSI from STAR SJ.out.tab
# ---------------------------------------------------------------------------


def _strip_chr_prefix(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


# Optional global mapping {refseq_accession -> primary chrom name}, populated
# from --chrom-map by main(). Built once so _build_chrom_aliases stays cheap.
_CHROM_MAP: Dict[str, str] = {}


def _load_chrom_map(path: Path) -> Dict[str, str]:
    """Parse an NCBI assembly report (or simple 2-col TSV) into ``{refseq->name}``.

    NCBI ``assembly_report.txt`` layout (after the ``#`` comment header):
      ``Sequence-Name<TAB>...<TAB>RefSeq-Accn<TAB>...``  (RefSeq is col index 6)
    We accept either that layout or a simple ``alias<TAB>name`` TSV.
    """
    mapping: Dict[str, str] = {}
    if not path.exists():
        return mapping
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 7 and parts[6] and parts[6].lower() not in {"na", "none", ""}:
                seq_name = parts[0].strip()
                refseq = parts[6].strip()
                if refseq and seq_name:
                    mapping[refseq] = seq_name
                    mapping[refseq.split(".")[0]] = seq_name
            elif len(parts) >= 2:
                a, b = parts[0].strip(), parts[1].strip()
                if a and b:
                    mapping[a] = b
                    mapping[a.split(".")[0]] = b
    return mapping


def _build_chrom_aliases(chrom: str) -> List[str]:
    """Return all reasonable spellings of a chromosome.

    POISEN uses RefSeq accessions like ``NC_133177.1``; STAR runs are typically
    built against ``chrN``-style names. When ``--chrom-map`` (an NCBI assembly
    report) is provided, ``_CHROM_MAP`` lets us translate the accession to the
    canonical sequence name and we add all reasonable spellings as aliases.
    Without a map we still try the ``chr``/no-``chr`` permutations so a same-style
    run still matches transparently.
    """
    aliases = {chrom}
    no_chr = _strip_chr_prefix(chrom)
    aliases.add(no_chr)
    if not chrom.startswith("chr"):
        aliases.add(f"chr{chrom}")
    mapped = _CHROM_MAP.get(chrom) or _CHROM_MAP.get(chrom.split(".")[0])
    if mapped:
        aliases.add(mapped)
        aliases.add(_strip_chr_prefix(mapped))
        if not mapped.startswith("chr"):
            aliases.add(f"chr{mapped}")
    return list(aliases)


def _index_sj_file(sj_path: Path, chrom_set: set) -> Dict[Tuple[str, int, int], int]:
    """Read a STAR SJ.out.tab into a dict keyed by (chrom, intron_start, intron_end).

    Format (1-based inclusive intron boundaries):
      col1=chrom col2=intron_start col3=intron_end col4=strand col5=motif
      col6=annotated col7=unique_reads col8=multi_reads col9=max_overhang
    """
    out: Dict[Tuple[str, int, int], int] = {}
    if not sj_path.exists():
        return out
    opener = gzip.open if sj_path.suffix == ".gz" else open
    with opener(str(sj_path), "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            chrom = parts[0]
            if chrom_set and chrom not in chrom_set:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
                uniq = int(parts[6])
            except ValueError:
                continue
            if uniq <= 0:
                continue
            key = (chrom, start, end)
            # Multiple alignments can list the same junction (rare); sum
            out[key] = out.get(key, 0) + uniq
    return out


def _resolve_sj_files(
    sj_list: List[Path],
    samples_tsv: Optional[Path],
) -> List[Tuple[str, str, Path]]:
    """Return list of (sample_id, condition, sj_path) tuples.

    The sample_id is the basename of the SJ file with ``.SJ.out.tab`` stripped, which
    matches how scripts/zebrafish_quest_pipeline.py:write_samples_tsv writes the
    ``sample`` column.
    """
    cond_map: Dict[str, str] = {}
    if samples_tsv is not None and samples_tsv.exists():
        with open(samples_tsv, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                s = (row.get("sample") or "").strip()
                c = (row.get("condition") or "").strip()
                if s:
                    cond_map[s] = c
                    cond_map[Path(s).name] = c

    triples: List[Tuple[str, str, Path]] = []
    for p in sj_list:
        name = p.name
        for suf in [".SJ.out.tab.gz", ".SJ.out.tab"]:
            if name.endswith(suf):
                sample_id = name[: -len(suf)]
                break
        else:
            sample_id = p.stem
        cond = cond_map.get(sample_id) or cond_map.get(name) or p.parent.name
        triples.append((sample_id, cond, p))
    return triples


def compute_pe_psi(
    events_df: pd.DataFrame,
    sj_files: List[Tuple[str, str, Path]],
    *,
    min_junction_reads: int,
    min_locus_reads: int,
    max_intron_len: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Stage B.

    For every (event, sample) pair, look up the three relevant introns in the
    sample's STAR SJ.out.tab and compute inclusion PSI under LC2-equivalent hygiene.

    PSI formula:
        incl  = mean(reads_5p_intron, reads_3p_intron)   # tolerant of one-sided dropout
        total = incl + reads_skip_intron
        PSI   = incl / total when total >= min_locus_reads, else NaN

    Returns three DataFrames:
      - psi_df    : index=event_id, columns=samples, float PSI (NaN allowed)
      - counts_df : long-form (event_id, sample, incl_5p, incl_3p, skip)
      - low_cov   : audit rows for dropped (event, sample) pairs
    """
    keep_mask = events_df["max_intron_len"] <= max_intron_len
    n_filtered_long = int((~keep_mask).sum())
    events_use = events_df[keep_mask].reset_index(drop=True)

    chroms_needed = set(events_use["chrom"])
    chrom_alias_lookup: Dict[str, List[str]] = {
        c: _build_chrom_aliases(c) for c in chroms_needed
    }
    full_alias_set = set()
    for aliases in chrom_alias_lookup.values():
        full_alias_set.update(aliases)

    n_events = len(events_use)
    n_samples = len(sj_files)
    psi_arr = np.full((n_events, n_samples), np.nan, dtype=np.float64)
    incl_5p_arr = np.zeros((n_events, n_samples), dtype=np.int32)
    incl_3p_arr = np.zeros((n_events, n_samples), dtype=np.int32)
    skip_arr = np.zeros((n_events, n_samples), dtype=np.int32)

    low_cov_records: List[Dict] = []

    for j, (sample_id, cond, sj_path) in enumerate(sj_files):
        sj_idx = _index_sj_file(sj_path, full_alias_set)
        if not sj_idx:
            for i, evt in events_use.iterrows():
                low_cov_records.append({
                    "event_id": evt["event_id"],
                    "sample": sample_id,
                    "condition": cond,
                    "reason": "sj_file_empty_or_missing",
                })
            continue

        for i, evt in events_use.iterrows():
            chrom = evt["chrom"]
            aliases = chrom_alias_lookup[chrom]

            def _lookup(start: int, end: int) -> int:
                for a in aliases:
                    v = sj_idx.get((a, start, end))
                    if v is not None:
                        return v
                return 0

            r5 = _lookup(int(evt["intron_5p_genomic_start"]), int(evt["intron_5p_genomic_end"]))
            r3 = _lookup(int(evt["intron_3p_genomic_start"]), int(evt["intron_3p_genomic_end"]))
            rs = _lookup(int(evt["intron_skip_start"]), int(evt["intron_skip_end"]))

            r5 = r5 if r5 >= min_junction_reads else 0
            r3 = r3 if r3 >= min_junction_reads else 0
            rs = rs if rs >= min_junction_reads else 0

            incl_5p_arr[i, j] = r5
            incl_3p_arr[i, j] = r3
            skip_arr[i, j] = rs

            incl_pair = [v for v in (r5, r3) if v > 0]
            if not incl_pair and rs == 0:
                low_cov_records.append({
                    "event_id": evt["event_id"],
                    "sample": sample_id,
                    "condition": cond,
                    "reason": "no_supporting_reads",
                })
                continue
            incl_mean = float(np.mean(incl_pair)) if incl_pair else 0.0
            total = incl_mean + rs
            if total < min_locus_reads:
                low_cov_records.append({
                    "event_id": evt["event_id"],
                    "sample": sample_id,
                    "condition": cond,
                    "reason": f"locus_reads_below_{min_locus_reads}",
                })
                continue
            psi_arr[i, j] = incl_mean / total

    sample_cols = [s for s, _, _ in sj_files]

    psi_df = pd.DataFrame(psi_arr, columns=sample_cols)
    psi_df.insert(0, "event_id", events_use["event_id"].values)

    # Long-form counts table
    counts_records = []
    for i, evt_id in enumerate(events_use["event_id"].values):
        for j, sample_id in enumerate(sample_cols):
            counts_records.append({
                "event_id": evt_id,
                "sample": sample_id,
                "incl_5p_reads": int(incl_5p_arr[i, j]),
                "incl_3p_reads": int(incl_3p_arr[i, j]),
                "skip_reads": int(skip_arr[i, j]),
            })
    counts_df = pd.DataFrame.from_records(counts_records)

    low_cov_df = pd.DataFrame.from_records(low_cov_records)

    if n_filtered_long:
        # Note: events filtered out for max_intron_len don't appear in psi/counts
        for _, evt in events_df[~keep_mask].iterrows():
            low_cov_records.append({
                "event_id": evt["event_id"],
                "sample": "(all)",
                "condition": "(all)",
                "reason": f"intron_len_above_{max_intron_len}",
            })
    return psi_df, counts_df, low_cov_df


# ---------------------------------------------------------------------------
# Stage C — differential inclusion across development
# ---------------------------------------------------------------------------


def run_differential_inclusion(
    events_df: pd.DataFrame,
    psi_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    sj_files: List[Tuple[str, str, Path]],
    *,
    min_locus_reads: int = DEFAULT_MIN_LOCUS_READS,
    min_observability_frac: float = DEFAULT_MIN_OBSERVABILITY_FRAC,
    min_timepoints_observed: int = DEFAULT_MIN_TIMEPOINTS_OBSERVED,
    min_total_locus_reads: int = DEFAULT_MIN_TOTAL_LOCUS_READS,
    min_dpsi: float = DEFAULT_MIN_DPSI,
    fdr_threshold: float = DEFAULT_FDR,
    run_ds: bool = False,
    ds_workdir: Optional[Path] = None,
    leafcutter_repo: Optional[Path] = None,
) -> Tuple[pd.DataFrame, List[str], Dict]:
    """Stage C — differential PE inclusion across development.

    Primary test: per-event quasi-binomial GLM of (inclusion, skip) counts on
    developmental hpf (continuous predictor).  This is the count-based,
    overdispersion-aware reference test in the field (LeafCutter2 supplementary
    methods; Mudge & Pritchard 2024 Nat Genet methods; ``glm(family =
    quasibinomial)`` in R).  We report a Wald p-value, BH-FDR, and a literature-
    standard biological-significance bar of |dPSI| >= ``min_dpsi`` (default 0.10)
    measured between the earliest and latest *observed* timepoints for the event.

    Filters applied per literature standard before testing:
      * Per-(event,sample) coverage >= ``min_locus_reads`` (default 10).
      * Event observed in >= ``min_observability_frac`` of developmental samples
        (LeafCutter2 default 60%).
      * Event observed in >= ``min_timepoints_observed`` distinct timepoints.
      * Total locus reads across all samples >= ``min_total_locus_reads``
        (LeafCutter2 "appreciable abundance" floor, default 50).

    Secondary outputs (cross-check):
      * Spearman rho of per-sample PSI vs hpf (rank-based, distribution-free).
      * Dirichlet-multinomial test via ``leafcutter_ds.R`` first-vs-last
        timepoint contrast when ``run_ds=True`` and ``leafcutter_repo`` is set.

    Returns (results_df, ordered_timepoints, stage_stats).
    """
    sample_cond: Dict[str, str] = {s: c for s, c, _ in sj_files}
    sample_to_hours: Dict[str, float] = {s: _timepoint_to_hours(c) for s, c in sample_cond.items()}

    valid_dev_samples = [s for s, h in sample_to_hours.items() if math.isfinite(h)]
    n_dev_samples = len(valid_dev_samples)
    if not valid_dev_samples:
        print(
            "[STAGE C][WARN] No samples have parseable developmental timepoints; "
            "GLM test will be skipped.",
            flush=True,
        )

    cond_to_hours: Dict[str, float] = {}
    for c in {sample_cond[s] for s in valid_dev_samples}:
        cond_to_hours[c] = _timepoint_to_hours(c)
    ordered_timepoints = sorted(cond_to_hours.keys(), key=lambda t: cond_to_hours[t])
    timepoint_hours = np.array([cond_to_hours[t] for t in ordered_timepoints], dtype=float)
    cond_to_samples: Dict[str, List[str]] = defaultdict(list)
    for s, c in sample_cond.items():
        cond_to_samples[c].append(s)

    psi_indexed = psi_df.set_index("event_id")

    # Index per-(event,sample) counts for O(1) lookup
    counts_indexed = counts_df.set_index(["event_id", "sample"])
    counts_by_event = counts_df.groupby("event_id")[["incl_5p_reads", "incl_3p_reads", "skip_reads"]].sum()

    # Per-event accumulators
    glm_beta1_list: List[Optional[float]] = []
    glm_se_list: List[Optional[float]] = []
    glm_p_list: List[Optional[float]] = []
    glm_phi_list: List[Optional[float]] = []
    glm_n_obs_list: List[int] = []
    rho_list: List[Optional[float]] = []
    spearman_p_list: List[Optional[float]] = []
    n_samples_observed_list: List[int] = []
    n_timepoints_observed_list: List[int] = []
    per_tp_means: Dict[str, List[Optional[float]]] = {tp: [] for tp in ordered_timepoints}
    per_tp_count_psi: Dict[str, List[Optional[float]]] = {tp: [] for tp in ordered_timepoints}
    total_reads_list: List[int] = []
    dpsi_list: List[Optional[float]] = []
    first_observed_tp_list: List[Optional[str]] = []
    last_observed_tp_list: List[Optional[str]] = []
    pass_observability_list: List[bool] = []
    pass_total_reads_list: List[bool] = []
    pass_timepoints_list: List[bool] = []
    tested_list: List[bool] = []
    skip_reason_list: List[Optional[str]] = []

    min_observed_samples = max(1, int(math.ceil(min_observability_frac * n_dev_samples))) if n_dev_samples else 0

    for evt_id in events_df["event_id"].values:
        # Default no-data row
        if evt_id not in psi_indexed.index or n_dev_samples == 0:
            for tp in ordered_timepoints:
                per_tp_means[tp].append(None)
                per_tp_count_psi[tp].append(None)
            glm_beta1_list.append(None); glm_se_list.append(None); glm_p_list.append(None)
            glm_phi_list.append(None); glm_n_obs_list.append(0)
            rho_list.append(None); spearman_p_list.append(None)
            n_samples_observed_list.append(0); n_timepoints_observed_list.append(0)
            total_reads_list.append(0); dpsi_list.append(None)
            first_observed_tp_list.append(None); last_observed_tp_list.append(None)
            pass_observability_list.append(False); pass_total_reads_list.append(False)
            pass_timepoints_list.append(False); tested_list.append(False)
            skip_reason_list.append("no_psi_row")
            continue

        psi_row = psi_indexed.loc[evt_id]

        # Per-timepoint sample-mean PSI (for the candidate table & plots).
        tp_means: List[Optional[float]] = []
        for tp in ordered_timepoints:
            samples = cond_to_samples.get(tp, [])
            vals = np.array([psi_row.get(s, np.nan) for s in samples], dtype=float)
            vals = vals[~np.isnan(vals)]
            tp_means.append(float(np.mean(vals)) if len(vals) > 0 else None)
            per_tp_means[tp].append(tp_means[-1])

        # Build (hpf, k=incl, n=incl+skip) arrays for samples that pass
        # the per-(event,sample) coverage filter.
        per_sample_x: List[float] = []
        per_sample_k: List[int] = []
        per_sample_n: List[int] = []
        per_sample_psi: List[float] = []
        per_sample_tp: List[str] = []
        for s in valid_dev_samples:
            key = (evt_id, s)
            if key not in counts_indexed.index:
                continue
            r = counts_indexed.loc[key]
            r5 = int(r["incl_5p_reads"]); r3 = int(r["incl_3p_reads"])
            rs = int(r["skip_reads"])
            # Tolerate one-sided dropout: take the better-supported flank.
            incl = max(r5, r3)
            total = incl + rs
            if total < min_locus_reads:
                continue
            per_sample_x.append(sample_to_hours[s])
            per_sample_k.append(incl)
            per_sample_n.append(total)
            per_sample_psi.append(incl / total)
            per_sample_tp.append(sample_cond[s])

        n_obs = len(per_sample_x)
        n_tps_obs = len(set(per_sample_tp))
        n_samples_observed_list.append(n_obs)
        n_timepoints_observed_list.append(n_tps_obs)

        # Total locus reads (across all samples that passed coverage)
        total_event_reads = int(sum(per_sample_n))
        total_reads_list.append(total_event_reads)

        # Per-timepoint count-pooled PSI = sum(incl) / sum(total) within tp.
        # This is the literature-standard "condition PSI" used by LeafCutter2 and
        # Mudge 2024 -- properly read-weighted, unlike the mean-of-PSIs.
        tp_count_psi: Dict[str, Optional[float]] = {tp: None for tp in ordered_timepoints}
        tp_pool_k: Dict[str, int] = defaultdict(int)
        tp_pool_n: Dict[str, int] = defaultdict(int)
        for tp_label, k_, n_ in zip(per_sample_tp, per_sample_k, per_sample_n):
            tp_pool_k[tp_label] += k_
            tp_pool_n[tp_label] += n_
        for tp in ordered_timepoints:
            if tp_pool_n[tp] > 0:
                tp_count_psi[tp] = tp_pool_k[tp] / tp_pool_n[tp]
            per_tp_count_psi[tp].append(tp_count_psi[tp])

        # First / last *observed* timepoint and dPSI between them.
        observed_tps = [tp for tp in ordered_timepoints if tp_count_psi[tp] is not None]
        if observed_tps:
            first_tp = observed_tps[0]
            last_tp = observed_tps[-1]
            dpsi = tp_count_psi[last_tp] - tp_count_psi[first_tp]
        else:
            first_tp = last_tp = None
            dpsi = None
        first_observed_tp_list.append(first_tp)
        last_observed_tp_list.append(last_tp)
        dpsi_list.append(dpsi)

        # Filter pass flags
        passes_obs = (n_obs >= min_observed_samples)
        passes_tps = (n_tps_obs >= min_timepoints_observed)
        passes_reads = (total_event_reads >= min_total_locus_reads)
        pass_observability_list.append(passes_obs)
        pass_timepoints_list.append(passes_tps)
        pass_total_reads_list.append(passes_reads)

        if not (passes_obs and passes_tps and passes_reads):
            glm_beta1_list.append(None); glm_se_list.append(None); glm_p_list.append(None)
            glm_phi_list.append(None); glm_n_obs_list.append(n_obs)
            rho_list.append(None); spearman_p_list.append(None)
            tested_list.append(False)
            if not passes_reads:
                skip_reason_list.append(f"total_reads_below_{min_total_locus_reads}")
            elif not passes_obs:
                skip_reason_list.append(
                    f"observed_in_<{int(round(min_observability_frac*100))}pct_samples"
                )
            else:
                skip_reason_list.append(
                    f"observed_in_<{min_timepoints_observed}_timepoints"
                )
            continue

        # Primary: quasi-binomial GLM (count-based, hpf as continuous predictor)
        glm_fit = _quasi_binomial_glm(
            np.array(per_sample_x, dtype=float),
            np.array(per_sample_k, dtype=float),
            np.array(per_sample_n, dtype=float),
        )

        # Secondary: rank-based Spearman as a distribution-free cross-check
        rho, sp = _spearman_rho(
            np.array(per_sample_x, dtype=float),
            np.array(per_sample_psi, dtype=float),
        )

        if glm_fit is None:
            glm_beta1_list.append(None); glm_se_list.append(None); glm_p_list.append(None)
            glm_phi_list.append(None); glm_n_obs_list.append(n_obs)
            rho_list.append(rho); spearman_p_list.append(sp)
            tested_list.append(False)
            skip_reason_list.append("glm_uninformative")
            continue

        glm_beta1_list.append(glm_fit["beta1"])
        glm_se_list.append(glm_fit["se_beta1"])
        glm_p_list.append(glm_fit["p_value"])
        glm_phi_list.append(glm_fit["phi"])
        glm_n_obs_list.append(glm_fit["n_obs"])
        rho_list.append(rho); spearman_p_list.append(sp)
        tested_list.append(True)
        skip_reason_list.append(None)

    glm_fdrs = _bh_fdr(glm_p_list)
    spearman_fdrs = _bh_fdr(spearman_p_list)

    results = events_df.copy().set_index("event_id", drop=False)
    for tp in ordered_timepoints:
        results[f"psi_{tp}"] = per_tp_means[tp]
        results[f"psi_count_{tp}"] = per_tp_count_psi[tp]

    results["glm_beta1_per_hpf"] = glm_beta1_list
    results["glm_se_beta1"] = glm_se_list
    results["glm_p"] = glm_p_list
    results["glm_phi"] = glm_phi_list
    results["glm_n_samples_used"] = glm_n_obs_list
    results["glm_bh_fdr"] = glm_fdrs

    # Backward-compat alias: webapp / downstream consumers expect ``bh_fdr``.
    results["bh_fdr"] = glm_fdrs

    results["spearman_rho"] = rho_list
    results["spearman_p"] = spearman_p_list
    results["spearman_bh_fdr"] = spearman_fdrs

    results["n_samples_observed"] = n_samples_observed_list
    results["n_timepoints_observed"] = n_timepoints_observed_list
    results["total_event_reads"] = total_reads_list
    results["dpsi_first_to_last"] = dpsi_list
    results["first_observed_timepoint"] = first_observed_tp_list
    results["last_observed_timepoint"] = last_observed_tp_list
    results["passes_observability_filter"] = pass_observability_list
    results["passes_timepoint_filter"] = pass_timepoints_list
    results["passes_total_reads_filter"] = pass_total_reads_list
    results["tested"] = tested_list
    results["skip_reason"] = skip_reason_list

    # Combined biological + statistical significance bar.
    def _combined(row) -> bool:
        q = row["glm_bh_fdr"]
        d = row["dpsi_first_to_last"]
        if q is None or d is None:
            return False
        try:
            return (float(q) < fdr_threshold) and (abs(float(d)) >= min_dpsi)
        except (TypeError, ValueError):
            return False
    results["passes_significance"] = results.apply(_combined, axis=1)

    ds_stats: Dict = {}
    if run_ds:
        ds_stats = _run_leafcutter_ds(
            events_df, counts_df, sj_files, ordered_timepoints,
            cond_to_samples, ds_workdir, leafcutter_repo, results,
        )

    stage_stats = {
        # Filter funnel (literature-standard)
        "min_locus_reads": int(min_locus_reads),
        "min_observability_frac": float(min_observability_frac),
        "min_observed_samples_required": int(min_observed_samples),
        "min_timepoints_observed_required": int(min_timepoints_observed),
        "min_total_locus_reads": int(min_total_locus_reads),
        "min_dpsi": float(min_dpsi),
        "fdr_threshold": float(fdr_threshold),
        "n_events_passing_observability": int(sum(pass_observability_list)),
        "n_events_passing_timepoints": int(sum(pass_timepoints_list)),
        "n_events_passing_total_reads": int(sum(pass_total_reads_list)),
        "n_events_passing_all_filters": int(
            sum(o and t and r for o, t, r in zip(
                pass_observability_list, pass_timepoints_list, pass_total_reads_list))
        ),
        # Tests
        "n_events_tested_glm": int(sum(tested_list)),
        "n_significant_glm_q05": int(((results["glm_bh_fdr"].notna())
                                       & (results["glm_bh_fdr"] < fdr_threshold)).sum()),
        "n_significant_combined": int(results["passes_significance"].sum()),
        "n_significant_spearman_q05": int(((results["spearman_bh_fdr"].notna())
                                            & (results["spearman_bh_fdr"] < fdr_threshold)).sum()),
        # Backwards-compat keys
        "n_events_tested": int(sum(tested_list)),
        "n_significant_q05": int(results["passes_significance"].sum()),
        "n_significant_p05": int(((results["glm_p"].notna())
                                   & (results["glm_p"] < 0.05)).sum()),
        "ds": ds_stats,
        "ordered_timepoints": ordered_timepoints,
        "timepoint_hours": [float(h) for h in timepoint_hours],
    }
    return results.reset_index(drop=True), ordered_timepoints, stage_stats


def _run_leafcutter_ds(
    events_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    sj_files: List[Tuple[str, str, Path]],
    ordered_timepoints: List[str],
    cond_to_samples: Dict[str, List[str]],
    workdir: Optional[Path],
    leafcutter_repo: Optional[Path],
    results: pd.DataFrame,
) -> Dict:
    """Optional Dirichlet-multinomial confirmation via leafcutter_ds.R.

    Treats each PE event as a 2-junction pseudo-cluster {inclusion, skipping} and
    contrasts the first vs last developmental timepoint.  The results columns
    ``ds_log_effect`` and ``ds_p`` are appended to ``results`` in place.
    """
    if len(ordered_timepoints) < 2:
        return {"skipped": "fewer_than_two_timepoints"}
    if workdir is None:
        return {"skipped": "no_ds_workdir"}
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    if leafcutter_repo is None or not (Path(leafcutter_repo) / "scripts" / "leafcutter_ds.R").exists():
        return {"skipped": "leafcutter_repo_missing"}

    # Pseudo-cluster file: chrom:start:end:clu_<idx>:strand <count1> <count2> ...
    first_tp = ordered_timepoints[0]
    last_tp = ordered_timepoints[-1]
    samples_first = cond_to_samples.get(first_tp, [])
    samples_last = cond_to_samples.get(last_tp, [])
    contrast_samples = list(samples_first) + list(samples_last)
    if not samples_first or not samples_last:
        return {"skipped": "missing_samples_in_extreme_timepoints"}

    counts_indexed = counts_df.set_index(["event_id", "sample"])
    pseudo_path = workdir / "pe_pseudo_clusters.txt.gz"
    with gzip.open(str(pseudo_path), "wt") as fh:
        header = ["chrom"] + contrast_samples
        fh.write(" ".join(header) + "\n")
        for idx, evt in events_df.iterrows():
            evt_id = evt["event_id"]
            chrom = evt["chrom"]
            strand = evt["strand"] or "+"
            clu = f"clu_pe_{idx}"
            incl_id = f"{chrom}:{int(evt['intron_skip_start'])}:{int(evt['intron_skip_end'])}:{clu}_incl:{strand}"
            skip_id = f"{chrom}:{int(evt['intron_skip_start'])}:{int(evt['intron_skip_end'])}:{clu}_skip:{strand}"
            incl_row = [incl_id]
            skip_row = [skip_id]
            for s in contrast_samples:
                key = (evt_id, s)
                if key in counts_indexed.index:
                    r = counts_indexed.loc[key]
                    r5 = int(r["incl_5p_reads"])
                    r3 = int(r["incl_3p_reads"])
                    rs = int(r["skip_reads"])
                else:
                    r5 = r3 = rs = 0
                incl = max(r5, r3)
                incl_row.append(str(incl))
                skip_row.append(str(rs))
            fh.write(" ".join(incl_row) + "\n")
            fh.write(" ".join(skip_row) + "\n")

    groups_path = workdir / "pe_pseudo_groups.txt"
    with open(groups_path, "w") as fh:
        for s in samples_first:
            fh.write(f"{s}\t{first_tp}\n")
        for s in samples_last:
            fh.write(f"{s}\t{last_tp}\n")

    out_prefix = workdir / "pe_ds"
    cmd = [
        "Rscript",
        str(Path(leafcutter_repo) / "scripts" / "leafcutter_ds.R"),
        f"--num_threads=4",
        f"--output_prefix={out_prefix}",
        str(pseudo_path),
        str(groups_path),
    ]
    import subprocess
    print(f"[STAGE C][DS] Running leafcutter_ds: {' '.join(cmd)}", flush=True)
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[STAGE C][DS] Failed: {exc}", flush=True)
        return {"skipped": "leafcutter_ds_failed", "error": str(exc)}

    sig_path = Path(str(out_prefix) + "_cluster_significance.txt")
    if not sig_path.exists():
        return {"skipped": "no_ds_output"}
    sig_df = pd.read_csv(sig_path, sep="\t")
    cluster_to_p = dict(zip(sig_df["cluster"], sig_df["p"]))

    ds_p_col: List[Optional[float]] = []
    for idx in range(len(events_df)):
        clu_id = f"clu_pe_{idx}"
        # leafcutter_ds may key by chrom:cluster_id; try both
        match_keys = [k for k in cluster_to_p if k.endswith(clu_id)]
        if match_keys:
            ds_p_col.append(float(cluster_to_p[match_keys[0]]))
        else:
            ds_p_col.append(None)
    results["ds_p"] = ds_p_col
    results["ds_fdr"] = _bh_fdr(ds_p_col)
    return {
        "first_timepoint": first_tp,
        "last_timepoint": last_tp,
        "n_first": len(samples_first),
        "n_last": len(samples_last),
        "n_events_with_ds_p": int(sum(1 for v in ds_p_col if v is not None)),
    }


# ---------------------------------------------------------------------------
# Stage D — cross-reference annotation joins
# ---------------------------------------------------------------------------


def _load_longorf(longorf_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(longorf_tsv, sep="\t", dtype=str, keep_default_na=False)
    expected = {"transcript_id", "is_ptc", "last_ejc", "stop_genomic", "orf_length", "n_exons"}
    if not expected.issubset(df.columns):
        return pd.DataFrame()
    df = df.rename(columns={"transcript_id": "read_id"})

    def _last_ejc_distance(row: pd.Series) -> Optional[float]:
        """Return |stop_genomic - last_ejc| in nucleotides where parseable."""
        try:
            sg = float(row["stop_genomic"])
            lejc = float(row["last_ejc"])
            return abs(sg - lejc)
        except (ValueError, TypeError):
            return None

    df["longorf_supports_ptc"] = df["is_ptc"].astype(str).str.lower() == "true"
    df["longorf_last_ejc_distance"] = df.apply(_last_ejc_distance, axis=1)
    df["longorf_orf_length"] = pd.to_numeric(df["orf_length"], errors="coerce")
    df["longorf_n_exons"] = pd.to_numeric(df["n_exons"], errors="coerce")
    return df[[
        "read_id", "longorf_supports_ptc", "longorf_last_ejc_distance",
        "longorf_orf_length", "longorf_n_exons",
    ]]


def _load_lc2_classifications(class_path: Path) -> Dict[Tuple[str, int, int], str]:
    """Return ``{(chrom, intron_start, intron_end): classification}``.

    Tolerates the two flavours of LC2 classification format observed in this repo:
      a) gene<TAB>chrom:start-end<TAB>classification (lc2_pipeline.py convention)
      b) junction_id<TAB>classification<TAB>... where junction_id is
         ``chrom:start:end:cluster:strand`` (example_output convention)
    """
    out: Dict[Tuple[str, int, int], str] = {}
    if not class_path.exists():
        return out
    with open(class_path, "r") as fh:
        header_line = fh.readline()
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            # Try format (b): junction_id in col0 with embedded chrom:start:end
            j = parts[0]
            colon_count = j.count(":")
            cls = None
            chrom = start = end = None
            if colon_count >= 2:
                tok = j.split(":")
                try:
                    chrom = tok[0]
                    start = int(tok[1])
                    end = int(tok[2])
                    cls = parts[1].strip().upper() if len(parts) > 1 else None
                except (ValueError, IndexError):
                    chrom = start = end = None
            if chrom is None and len(parts) >= 3:
                # Fall back to format (a)
                coord_str = parts[1]
                m = re.match(r"([^:]+):(\d+)-(\d+)", coord_str)
                if m:
                    chrom = m.group(1)
                    start = int(m.group(2))
                    end = int(m.group(3))
                    cls = parts[2].strip().upper()
            if chrom is None or start is None or end is None or cls is None:
                continue
            out[(chrom, start, end)] = cls
    return out


def _load_lc2_cluster_ratios(
    cluster_ratios_path: Path,
    samples_tsv: Optional[Path],
) -> Tuple[Dict[Tuple[str, int, int], Dict[str, float]], List[str]]:
    """Parse LC2 cluster_ratios.gz and return ``{(chrom, start, end): {tp: mean PSI}}``.

    Same parsing strategy as scripts/lc2_pipeline.py:compute_poison_exon_analysis.
    """
    if not cluster_ratios_path.exists():
        return {}, []

    cond_map: Dict[str, str] = {}
    if samples_tsv is not None and samples_tsv.exists():
        with open(samples_tsv, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                s = (row.get("sample") or "").strip()
                c = (row.get("condition") or "").strip()
                if s:
                    cond_map[s] = c
                    cond_map[Path(s).stem] = c

    with gzip.open(str(cluster_ratios_path), "rt") as fh:
        header = fh.readline().strip().split()
        sample_cols = header[1:]
        sample_conds = [cond_map.get(c, cond_map.get(Path(c).stem, c)) for c in sample_cols]
        unique_tps = sorted(set(sample_conds), key=_timepoint_to_hours)

        out: Dict[Tuple[str, int, int], Dict[str, float]] = {}
        for line in fh:
            fields = line.strip().split()
            jid = fields[0]
            toks = jid.split(":")
            if len(toks) < 3:
                continue
            try:
                chrom = toks[0]
                start = int(toks[1])
                end = int(toks[2])
            except ValueError:
                continue

            psi_per_tp: Dict[str, List[float]] = defaultdict(list)
            for j in range(len(sample_cols)):
                token = fields[j + 1] if j + 1 < len(fields) else "0/0"
                if "/" in token:
                    num, den = token.split("/")
                    try:
                        num_f, den_f = float(num), float(den)
                    except ValueError:
                        continue
                    if den_f > 0:
                        psi_per_tp[sample_conds[j]].append(num_f / den_f)
            mean_per_tp = {tp: float(np.mean(vals)) for tp, vals in psi_per_tp.items() if vals}
            if mean_per_tp:
                out[(chrom, start, end)] = mean_per_tp
    return out, unique_tps


def annotate_results(
    results: pd.DataFrame,
    poisen_tsv: Path,
    longorf_tsv: Optional[Path],
    lc2_classifications: Optional[Path],
    lc2_cluster_ratios: Optional[Path],
    samples_tsv: Optional[Path],
) -> pd.DataFrame:
    """Stage D.

    Adds LongORF + LC2 cross-reference columns to the candidate table.  All
    joins are best-effort; missing sources just leave those columns blank.
    """
    # --- LongORF join via read_id (POISEN read_id == LongORF transcript_id) ---
    poisen_df = pd.read_csv(poisen_tsv, sep="\t", dtype=str, keep_default_na=False)
    if "read_id" not in poisen_df.columns:
        poisen_df["read_id"] = ""
    pe_to_reads = (
        poisen_df.groupby(["pe_chrom", "pe_start", "pe_end"], dropna=False)["read_id"]
        .agg(lambda s: set(s))
        .to_dict()
    )

    longorf_df = pd.DataFrame()
    if longorf_tsv is not None and longorf_tsv.exists():
        longorf_df = _load_longorf(longorf_tsv)

    longorf_supports = []
    longorf_ejc = []
    longorf_orflen = []
    longorf_nexons = []
    if not longorf_df.empty:
        longorf_idx = longorf_df.set_index("read_id")
        for _, row in results.iterrows():
            key = (row["chrom"], str(int(row["pe_start"])), str(int(row["pe_end"])))
            reads = pe_to_reads.get(key, set())
            matched = longorf_idx.index.intersection(reads)
            if len(matched) == 0:
                longorf_supports.append(None)
                longorf_ejc.append(None)
                longorf_orflen.append(None)
                longorf_nexons.append(None)
                continue
            sub = longorf_idx.loc[list(matched)]
            longorf_supports.append(bool(sub["longorf_supports_ptc"].any()))
            longorf_ejc.append(float(sub["longorf_last_ejc_distance"].dropna().mean()) if sub["longorf_last_ejc_distance"].notna().any() else None)
            longorf_orflen.append(float(sub["longorf_orf_length"].dropna().mean()) if sub["longorf_orf_length"].notna().any() else None)
            longorf_nexons.append(float(sub["longorf_n_exons"].dropna().mean()) if sub["longorf_n_exons"].notna().any() else None)
    else:
        longorf_supports = [None] * len(results)
        longorf_ejc = [None] * len(results)
        longorf_orflen = [None] * len(results)
        longorf_nexons = [None] * len(results)

    results["longorf_supports_ptc"] = longorf_supports
    results["longorf_last_ejc_distance"] = longorf_ejc
    results["longorf_orf_length"] = longorf_orflen
    results["longorf_n_exons"] = longorf_nexons

    # --- LC2 classifications: per-intron lookup ---
    lc2_class = (
        _load_lc2_classifications(lc2_classifications)
        if lc2_classifications is not None else {}
    )

    def _chrom_candidates(chrom: str) -> List[str]:
        """All reasonable spellings of ``chrom`` for LC2 lookups.

        Includes the assembly-report mapping (RefSeq -> chr*) so POISEN's
        ``NC_133180.1`` resolves to ``chr5`` etc.
        """
        out = [chrom, _strip_chr_prefix(chrom), f"chr{chrom}"]
        mapped = _CHROM_MAP.get(chrom) or _CHROM_MAP.get(chrom.split(".")[0])
        if mapped:
            out += [mapped, _strip_chr_prefix(mapped), f"chr{mapped}"]
        return out

    # LC2 stores intron coordinates as 0-based fully-closed (= STAR 1-based
    # first/last base of intron, each minus 1). Some forks instead use the
    # STAR-native 1-based form, so we probe both spellings.
    def _coord_candidates(start: int, end: int) -> List[Tuple[int, int]]:
        return [(start - 1, end - 1), (start, end), (start - 1, end)]

    def _lookup_class(chrom: str, start: int, end: int) -> Optional[str]:
        if not lc2_class:
            return None
        for c in _chrom_candidates(chrom):
            for s, e in _coord_candidates(start, end):
                v = lc2_class.get((c, s, e))
                if v is not None:
                    return v
        return None

    lc2_donor: List[Optional[str]] = []
    lc2_acceptor: List[Optional[str]] = []
    lc2_skip: List[Optional[str]] = []
    for _, row in results.iterrows():
        lc2_donor.append(_lookup_class(row["chrom"], int(row["intron_5p_genomic_start"]), int(row["intron_5p_genomic_end"])))
        lc2_acceptor.append(_lookup_class(row["chrom"], int(row["intron_3p_genomic_start"]), int(row["intron_3p_genomic_end"])))
        lc2_skip.append(_lookup_class(row["chrom"], int(row["intron_skip_start"]), int(row["intron_skip_end"])))
    results["lc2_class_5p_intron"] = lc2_donor
    results["lc2_class_3p_intron"] = lc2_acceptor
    results["lc2_class_skip_intron"] = lc2_skip
    results["lc2_seen_any_intron"] = [
        any(v is not None for v in (a, b, c))
        for a, b, c in zip(lc2_donor, lc2_acceptor, lc2_skip)
    ]
    results["lc2_called_up"] = [
        any(v == "UP" for v in (a, b, c))
        for a, b, c in zip(lc2_donor, lc2_acceptor, lc2_skip)
    ]

    # --- LC2 cluster ratios: per-timepoint PSI for the skip junction ---
    if lc2_cluster_ratios is not None and lc2_cluster_ratios.exists():
        ratios, lc2_tps = _load_lc2_cluster_ratios(lc2_cluster_ratios, samples_tsv)

        def _lookup_ratio(chrom: str, start: int, end: int) -> Dict[str, float]:
            if not ratios:
                return {}
            for c in _chrom_candidates(chrom):
                for s, e in _coord_candidates(start, end):
                    v = ratios.get((c, s, e))
                    if v is not None:
                        return v
            return {}

        for tp in lc2_tps:
            results[f"lc2_psi_{tp}"] = [
                _lookup_ratio(row["chrom"], int(row["intron_skip_start"]), int(row["intron_skip_end"])).get(tp)
                for _, row in results.iterrows()
            ]

    return results


# ---------------------------------------------------------------------------
# Stage E — outputs
# ---------------------------------------------------------------------------


def _abs_or_none(v) -> float:
    if v is None:
        return -1.0
    try:
        if math.isnan(v):
            return -1.0
        return abs(float(v))
    except (TypeError, ValueError):
        return -1.0


def write_outputs(
    *,
    outdir: Path,
    events_df: pd.DataFrame,
    psi_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    low_cov_df: pd.DataFrame,
    results: pd.DataFrame,
    ordered_timepoints: List[str],
    stage_a_stats: Dict,
    stage_c_stats: Dict,
    parameters: Dict,
) -> Dict[str, str]:
    """Stage E.

    Writes the canonical output tables, the headline candidate TSV, summary.json,
    and the four matplotlib figures.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    events_path = outdir / "pe_events.tsv"
    events_df.to_csv(events_path, sep="\t", index=False)

    psi_path = outdir / "pe_inclusion_psi.tsv.gz"
    psi_df.to_csv(psi_path, sep="\t", index=False, compression="gzip")

    counts_path = outdir / "pe_inclusion_counts.tsv.gz"
    counts_df.to_csv(counts_path, sep="\t", index=False, compression="gzip")

    low_cov_path = outdir / "pe_low_coverage.tsv"
    low_cov_df.to_csv(low_cov_path, sep="\t", index=False)

    candidates = results.copy()
    candidates["abs_dpsi"] = candidates["dpsi_first_to_last"].apply(_abs_or_none)
    # Primary sort key: combined-significance flag, then GLM FDR, then |dPSI|.
    candidates["_sig_rank"] = (~candidates["passes_significance"].fillna(False)).astype(int)
    candidates = candidates.sort_values(
        by=["_sig_rank", "glm_bh_fdr", "abs_dpsi"],
        ascending=[True, True, False],
        na_position="last",
    ).drop(columns=["_sig_rank", "abs_dpsi"])

    cand_path = outdir / "pe_developmental_candidates.tsv"
    candidates.to_csv(cand_path, sep="\t", index=False)

    cand_json_path = outdir / "pe_developmental_candidates.json"
    with open(cand_json_path, "w") as fh:
        json.dump({
            "timepoints": ordered_timepoints,
            "candidates": json.loads(candidates.to_json(orient="records")),
        }, fh, indent=2)

    heatmap_png = outdir / "pe_psi_heatmap.png"
    psi_vs_time_png = outdir / "pe_psi_vs_time_top.png"
    volcano_png = outdir / "pe_dpsi_vs_pvalue.png"
    spearman_dist_png = outdir / "pe_spearman_distribution.png"
    concordance_png = outdir / "pe_lc2_concordance.png"
    try:
        _plot_psi_heatmap(candidates, ordered_timepoints, heatmap_png)
        _plot_psi_vs_time_top(candidates, ordered_timepoints, psi_vs_time_png)
        _plot_dpsi_volcano(candidates, volcano_png,
                           min_dpsi=stage_c_stats.get("min_dpsi", DEFAULT_MIN_DPSI),
                           fdr_threshold=stage_c_stats.get("fdr_threshold", DEFAULT_FDR))
        _plot_spearman_distribution(candidates, spearman_dist_png)
        if "lc2_called_up" in candidates.columns and candidates["lc2_seen_any_intron"].any():
            _plot_lc2_concordance(candidates, concordance_png)
        else:
            concordance_png = None
    except Exception as exc:
        print(f"[STAGE E][WARN] Plot generation failed: {exc}", flush=True)
        volcano_png = volcano_png if volcano_png and volcano_png.exists() else None

    summary = {
        "stage": "stage3_pe_inclusion",
        "parameters": parameters,
        "stage_a": stage_a_stats,
        "stage_c": stage_c_stats,
        "n_candidates": int(len(candidates)),
        "n_events_in_psi": int(len(psi_df)),
        "ordered_timepoints": ordered_timepoints,
        "artifacts": {
            "pe_events": events_path.name,
            "pe_inclusion_psi": psi_path.name,
            "pe_inclusion_counts": counts_path.name,
            "pe_low_coverage": low_cov_path.name,
            "pe_developmental_candidates": cand_path.name,
            "pe_developmental_candidates_json": cand_json_path.name,
            "pe_psi_heatmap_png": heatmap_png.name if heatmap_png else None,
            "pe_psi_vs_time_top_png": psi_vs_time_png.name if psi_vs_time_png else None,
            "pe_dpsi_vs_pvalue_png": volcano_png.name if volcano_png else None,
            "pe_spearman_distribution_png": spearman_dist_png.name if spearman_dist_png else None,
            "pe_lc2_concordance_png": concordance_png.name if concordance_png else None,
        },
        "top_hits": json.loads(
            candidates.head(20)[[
                c for c in [
                    "event_id", "gene_symbol", "chrom", "pe_start", "pe_end",
                    "glm_beta1_per_hpf", "glm_p", "glm_bh_fdr",
                    "dpsi_first_to_last", "first_observed_timepoint",
                    "last_observed_timepoint", "passes_significance",
                    "spearman_rho", "spearman_p",
                    "longorf_supports_ptc", "lc2_called_up",
                ] if c in candidates.columns
            ]].to_json(orient="records")
        ),
    }
    summary_path = outdir / "summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    return {k: str(v) for k, v in {
        "summary": summary_path,
        "events": events_path,
        "psi": psi_path,
        "counts": counts_path,
        "low_cov": low_cov_path,
        "candidates": cand_path,
        "candidates_json": cand_json_path,
        "heatmap_png": heatmap_png,
        "psi_vs_time_png": psi_vs_time_png,
        "dpsi_volcano_png": volcano_png,
        "spearman_distribution_png": spearman_dist_png,
        "lc2_concordance_png": concordance_png,
    }.items() if v is not None}


def _plot_psi_heatmap(candidates: pd.DataFrame, timepoints: List[str], out_path: Path, top_n: int = 40) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    psi_cols = [f"psi_{tp}" for tp in timepoints]
    if not psi_cols:
        return
    sub = candidates.dropna(subset=psi_cols, how="all").head(top_n)
    if sub.empty:
        return
    mat = sub[psi_cols].astype(float).to_numpy()
    labels = [
        f"{r['gene_symbol'] or r['event_id']}" for _, r in sub.iterrows()
    ]

    fig, ax = plt.subplots(figsize=(max(5, 0.55 * len(timepoints) + 3), max(4, 0.18 * len(sub) + 1.5)), dpi=150)
    im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(timepoints)))
    ax.set_xticklabels(timepoints, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title(f"Top {len(sub)} POISEN PEs by FDR\nPSI across developmental timepoints",
                 fontsize=10, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Inclusion PSI", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_psi_vs_time_top(candidates: pd.DataFrame, timepoints: List[str], out_path: Path, top_n: int = 12) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    psi_cols = [f"psi_{tp}" for tp in timepoints]
    sub = candidates.dropna(subset=psi_cols, how="all").head(top_n)
    if sub.empty:
        return
    hours = np.array([_timepoint_to_hours(t) for t in timepoints], dtype=float)

    cols = 4
    rows = math.ceil(len(sub) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.0, rows * 2.4), dpi=150)
    axes = np.atleast_2d(axes).reshape(rows, cols)
    def _fmt(v, fmt):
        try:
            if v is None or (isinstance(v, float) and not math.isfinite(v)):
                return "n/a"
            return format(float(v), fmt)
        except (TypeError, ValueError):
            return "n/a"
    for k, (_, r) in enumerate(sub.iterrows()):
        ax = axes[k // cols][k % cols]
        ys = np.array([r.get(f"psi_{tp}", np.nan) for tp in timepoints], dtype=float)
        ys_count = np.array([r.get(f"psi_count_{tp}", np.nan) for tp in timepoints], dtype=float)
        # Plot count-pooled PSI (the literature-standard condition PSI) in colour
        # and the per-sample mean PSI as a faded reference.
        ax.plot(hours, ys, "o:", color="#aab7b8", lw=0.9, ms=3, label="mean of sample PSIs")
        ax.plot(hours, ys_count, "o-", color="#4e79a7", lw=1.4, ms=4.5, label="count-pooled PSI")
        ax.set_title(
            f"{r.get('gene_symbol') or r['event_id']}\n"
            f"\u0394PSI={_fmt(r.get('dpsi_first_to_last'), '+.2f')}, "
            f"q={_fmt(r.get('glm_bh_fdr'), '.2e')}",
            fontsize=8,
        )
        ax.set_ylim(0, 1)
        ax.set_xlabel("hpf", fontsize=7)
        ax.set_ylabel("PSI", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(alpha=0.25)
        if k == 0:
            ax.legend(fontsize=5, loc="best", framealpha=0.7)
    for k in range(len(sub), rows * cols):
        axes[k // cols][k % cols].axis("off")
    fig.suptitle("POISEN PE inclusion across development (top FDR)", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_spearman_distribution(candidates: pd.DataFrame, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rho = candidates["spearman_rho"].dropna().astype(float).to_numpy()
    if rho.size == 0:
        return
    fig, ax = plt.subplots(figsize=(6, 3.4), dpi=150)
    ax.hist(rho, bins=40, color="#59a14f", edgecolor="white", alpha=0.85)
    ax.axvline(0, color="#666666", linestyle="--", lw=0.8)
    n_pos = int((rho > 0).sum())
    n_neg = int((rho < 0).sum())
    ax.set_xlabel("Spearman \u03c1 (PSI vs hpf)", fontsize=9)
    ax.set_ylabel("# events", fontsize=9)
    ax.set_title(
        f"Distribution of PE Spearman \u03c1 across {len(rho)} events"
        f"\n({n_pos} positive, {n_neg} negative)",
        fontsize=10, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_dpsi_volcano(
    candidates: pd.DataFrame,
    out_path: Path,
    *,
    min_dpsi: float = DEFAULT_MIN_DPSI,
    fdr_threshold: float = DEFAULT_FDR,
) -> None:
    """Volcano: x = dPSI(last_observed - first_observed), y = -log10(GLM p)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sub = candidates.dropna(subset=["dpsi_first_to_last", "glm_p"]).copy()
    if sub.empty:
        return
    sub["neg_log10_p"] = -np.log10(np.clip(sub["glm_p"].astype(float), 1e-300, 1.0))
    dpsi_arr = sub["dpsi_first_to_last"].astype(float).to_numpy()
    p_arr = sub["neg_log10_p"].to_numpy()

    sig_mask = sub["passes_significance"].fillna(False).to_numpy().astype(bool)

    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=150)
    ax.scatter(dpsi_arr[~sig_mask], p_arr[~sig_mask],
               s=10, alpha=0.55, color="#bab0ac", label="not sig")
    if sig_mask.any():
        ax.scatter(dpsi_arr[sig_mask], p_arr[sig_mask],
                   s=18, alpha=0.95, color="#e15759",
                   label=f"FDR<{fdr_threshold} & |\u0394PSI|\u2265{min_dpsi}")
    ax.axhline(-math.log10(fdr_threshold), color="#666666", linestyle="--", lw=0.7)
    ax.axvline(min_dpsi, color="#666666", linestyle="--", lw=0.7)
    ax.axvline(-min_dpsi, color="#666666", linestyle="--", lw=0.7)
    ax.set_xlabel("\u0394PSI (last \u2212 first observed timepoint)", fontsize=9)
    ax.set_ylabel("\u2212log\u2081\u2080(GLM p)", fontsize=9)
    ax.set_title(
        f"Quasi-binomial GLM volcano (n={len(sub)} tested events; "
        f"{int(sig_mask.sum())} pass combined bar)",
        fontsize=10, fontweight="bold",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc="best", framealpha=0.8)

    # Label up to 8 most-significant points.
    label_sub = sub.assign(score=lambda d: d["neg_log10_p"] * d["dpsi_first_to_last"].abs())
    label_sub = label_sub.sort_values("score", ascending=False).head(8)
    for _, r in label_sub.iterrows():
        gene = r.get("gene_symbol") or r.get("event_id", "")
        ax.annotate(
            str(gene)[:18],
            xy=(float(r["dpsi_first_to_last"]), float(r["neg_log10_p"])),
            xytext=(3, 3), textcoords="offset points",
            fontsize=7, color="#333333",
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_lc2_concordance(candidates: pd.DataFrame, out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_total = len(candidates)
    n_seen = int(candidates["lc2_seen_any_intron"].sum()) if "lc2_seen_any_intron" in candidates.columns else 0
    n_up = int(candidates["lc2_called_up"].sum()) if "lc2_called_up" in candidates.columns else 0
    n_unseen = n_total - n_seen

    fig, ax = plt.subplots(figsize=(5, 3.4), dpi=150)
    bars = ax.bar(
        ["Not seen by LC2", "Seen, not UP", "Seen, called UP"],
        [n_unseen, max(0, n_seen - n_up), n_up],
        color=["#bab0ac", "#76b7b2", "#e15759"],
    )
    for b, v in zip(bars, [n_unseen, max(0, n_seen - n_up), n_up]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{v}",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("# POISEN PE events", fontsize=9)
    ax.set_title(f"Concordance of POISEN PEs with LC2 short-read calls (n={n_total})",
                 fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _resolve_sj_paths(sj_dir: Optional[Path], sj_list_file: Optional[Path]) -> List[Path]:
    paths: List[Path] = []
    if sj_list_file is not None and sj_list_file.exists():
        for line in sj_list_file.read_text().splitlines():
            line = line.strip()
            if line:
                paths.append(Path(line))
    if sj_dir is not None and sj_dir.exists():
        for ext in ("*.SJ.out.tab", "*.SJ.out.tab.gz"):
            paths.extend(sorted(sj_dir.rglob(ext)))
    seen = set()
    deduped: List[Path] = []
    for p in paths:
        key = str(p.resolve()) if p.exists() else str(p)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(p)
    return deduped


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="POISEN PE differential inclusion analysis (Saba's spec)",
    )
    ap.add_argument("--pe-list", required=True, type=Path,
                    help="POISEN PE TSV (UPDATED_STRAND_POISEN_HITS schema)")
    ap.add_argument("--sj-dir", type=Path, default=None,
                    help="Directory of STAR SJ.out.tab files (recursive)")
    ap.add_argument("--sj-list", type=Path, default=None,
                    help="Text file listing SJ.out.tab paths, one per line "
                         "(alternative to --sj-dir; both may be combined)")
    ap.add_argument("--samples-tsv", type=Path, default=None,
                    help="samples.tsv mapping sample -> condition (developmental "
                         "stage), as written by zebrafish_quest_pipeline.py")
    ap.add_argument("--longorf-tsv", type=Path, default=None,
                    help="LongORF_PTC+_fromClustered.tsv for cross-reference annotation")
    ap.add_argument("--lc2-classifications", type=Path, default=None,
                    help="LC2 *_lc2_junction_classifications.txt for cross-reference")
    ap.add_argument("--lc2-cluster-ratios", type=Path, default=None,
                    help="LC2 *_lc2.cluster_ratios.gz for per-timepoint PSI cross-reference")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="Output directory (created if missing)")
    ap.add_argument("--min-junction-reads", type=int, default=DEFAULT_MIN_JUNCTION_READS,
                    help="Per-junction read floor (mirrors LC2 pool_junc_reads)")
    ap.add_argument("--min-locus-reads", type=int, default=DEFAULT_MIN_LOCUS_READS,
                    help="Per-PE per-sample inclusion+skip read floor (LC2 default 10)")
    ap.add_argument("--min-observability-frac", type=float, default=DEFAULT_MIN_OBSERVABILITY_FRAC,
                    help="Event must be observed (above coverage) in at least this "
                         "fraction of developmental samples (LC2 default 0.60)")
    ap.add_argument("--min-timepoints-observed", type=int, default=DEFAULT_MIN_TIMEPOINTS_OBSERVED,
                    help="Event must have data in at least this many distinct "
                         "developmental timepoints to be tested")
    ap.add_argument("--min-total-locus-reads", type=int, default=DEFAULT_MIN_TOTAL_LOCUS_READS,
                    help="Aggregate inclusion+skip reads across all samples that "
                         "pass coverage; floor for testing (LC2 'appreciable "
                         "abundance' default 50)")
    ap.add_argument("--min-dpsi", type=float, default=DEFAULT_MIN_DPSI,
                    help="Biological-significance bar on |dPSI| between earliest "
                         "and latest *observed* timepoints (default 0.10)")
    ap.add_argument("--fdr-threshold", type=float, default=DEFAULT_FDR,
                    help="FDR threshold for combined-significance reporting")
    ap.add_argument("--max-intron-len", type=int, default=DEFAULT_MAX_INTRON_LEN,
                    help="Drop events whose flanking introns exceed this length")
    ap.add_argument("--ds", action="store_true",
                    help="Also run leafcutter Dirichlet-multinomial test on first vs "
                         "last developmental timepoint as a validation track. Auto-"
                         "enabled when --leafcutter-repo is provided.")
    ap.add_argument("--leafcutter-repo", type=Path, default=None,
                    help="Path to leafcutter repo (required for --ds)")
    ap.add_argument("--ds-workdir", type=Path, default=None,
                    help="Working directory for --ds intermediate files; defaults to <outdir>/_ds")
    ap.add_argument("--chrom-map", type=Path, default=None,
                    help="NCBI assembly_report.txt (or 2-col TSV) mapping POISEN's "
                         "RefSeq accessions (e.g. NC_133177.1) to STAR chrom names "
                         "(e.g. chr2). Required when POISEN and STAR use different "
                         "naming conventions.")
    args = ap.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.chrom_map is not None:
        global _CHROM_MAP
        _CHROM_MAP = _load_chrom_map(args.chrom_map)
        print(f"[INIT] Loaded {len(_CHROM_MAP)} chrom name mappings from {args.chrom_map}", flush=True)

    print(f"[STAGE A] Parsing POISEN PE list: {args.pe_list}", flush=True)
    events_df, stage_a_stats = parse_pe_events(args.pe_list)
    print(
        f"[STAGE A] {stage_a_stats['n_input_rows']} POISEN rows -> "
        f"{stage_a_stats['n_unique_events']} unique events "
        f"(dropped: {stage_a_stats['n_dropped_no_chain']} no-chain, "
        f"{stage_a_stats['n_dropped_pe_not_in_chain']} PE-not-in-chain, "
        f"{stage_a_stats['n_dropped_terminal_exon']} terminal-exon, "
        f"{stage_a_stats['n_dropped_chrom_mismatch']} chrom-mismatch)",
        flush=True,
    )

    sj_paths = _resolve_sj_paths(args.sj_dir, args.sj_list)
    if not sj_paths:
        print("[STAGE B][WARN] No SJ.out.tab files found; emitting events catalog only "
              "(skipping PSI computation and downstream stages).",
              flush=True)
        events_df.to_csv(args.outdir / "pe_events.tsv", sep="\t", index=False)
        with open(args.outdir / "summary.json", "w") as fh:
            json.dump({
                "stage": "stage3_pe_inclusion",
                "parameters": vars(args) | {"note": "events_only"},
                "stage_a": stage_a_stats,
                "n_sj_files": 0,
            }, fh, indent=2, default=str)
        return 0

    sj_files = _resolve_sj_files(sj_paths, args.samples_tsv)
    print(f"[STAGE B] Resolved {len(sj_files)} STAR SJ.out.tab files", flush=True)

    psi_df, counts_df, low_cov_df = compute_pe_psi(
        events_df, sj_files,
        min_junction_reads=args.min_junction_reads,
        min_locus_reads=args.min_locus_reads,
        max_intron_len=args.max_intron_len,
    )
    print(
        f"[STAGE B] PSI matrix: {psi_df.shape[0]} events x {psi_df.shape[1] - 1} samples; "
        f"{len(low_cov_df)} (event,sample) low-coverage records",
        flush=True,
    )

    print("[STAGE C] Running differential inclusion across development", flush=True)
    # Auto-enable DM validation track if a leafcutter repo is supplied.
    run_ds = bool(args.ds) or (args.leafcutter_repo is not None)
    results, ordered_timepoints, stage_c_stats = run_differential_inclusion(
        events_df, psi_df, counts_df, sj_files,
        min_locus_reads=args.min_locus_reads,
        min_observability_frac=args.min_observability_frac,
        min_timepoints_observed=args.min_timepoints_observed,
        min_total_locus_reads=args.min_total_locus_reads,
        min_dpsi=args.min_dpsi,
        fdr_threshold=args.fdr_threshold,
        run_ds=run_ds,
        ds_workdir=args.ds_workdir or (args.outdir / "_ds"),
        leafcutter_repo=args.leafcutter_repo,
    )
    print(
        f"[STAGE C] Filter funnel: "
        f"{stage_c_stats['n_events_passing_observability']} pass observability, "
        f"{stage_c_stats['n_events_passing_timepoints']} pass timepoints, "
        f"{stage_c_stats['n_events_passing_total_reads']} pass total-reads, "
        f"{stage_c_stats['n_events_passing_all_filters']} pass all",
        flush=True,
    )
    print(
        f"[STAGE C] GLM: {stage_c_stats['n_events_tested_glm']} tested; "
        f"{stage_c_stats['n_significant_glm_q05']} at FDR<{args.fdr_threshold}; "
        f"{stage_c_stats['n_significant_combined']} pass combined "
        f"(FDR<{args.fdr_threshold} AND |dPSI|>={args.min_dpsi}); "
        f"Spearman cross-check: {stage_c_stats['n_significant_spearman_q05']} at FDR<{args.fdr_threshold}",
        flush=True,
    )
    if stage_c_stats.get("ds"):
        print(f"[STAGE C][DS] {stage_c_stats['ds']}", flush=True)

    print("[STAGE D] Cross-reference annotation joins", flush=True)
    results = annotate_results(
        results, args.pe_list, args.longorf_tsv,
        args.lc2_classifications, args.lc2_cluster_ratios,
        args.samples_tsv,
    )

    print("[STAGE E] Writing outputs", flush=True)
    parameters = {
        "pe_list": str(args.pe_list),
        "samples_tsv": str(args.samples_tsv) if args.samples_tsv else None,
        "min_junction_reads": int(args.min_junction_reads),
        "min_locus_reads": int(args.min_locus_reads),
        "min_observability_frac": float(args.min_observability_frac),
        "min_timepoints_observed": int(args.min_timepoints_observed),
        "min_total_locus_reads": int(args.min_total_locus_reads),
        "min_dpsi": float(args.min_dpsi),
        "fdr_threshold": float(args.fdr_threshold),
        "max_intron_len": int(args.max_intron_len),
        "ds_enabled": bool(run_ds),
        "leafcutter_repo": str(args.leafcutter_repo) if args.leafcutter_repo else None,
    }
    artifacts = write_outputs(
        outdir=args.outdir,
        events_df=events_df,
        psi_df=psi_df,
        counts_df=counts_df,
        low_cov_df=low_cov_df,
        results=results,
        ordered_timepoints=ordered_timepoints,
        stage_a_stats=stage_a_stats,
        stage_c_stats=stage_c_stats,
        parameters=parameters,
    )
    print(f"[DONE] {len(results)} candidates -> {artifacts['candidates']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
