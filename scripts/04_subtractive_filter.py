#!/usr/bin/env python3
"""
Step 4: Subtractive Genomics & Target Prioritization
====================================================
Screens Core pangenome representative proteins through three subtractive filters:
1. DEG Essentiality Filter: Identifies indispensable genes (Database of Essential Genes).
2. Host (Human) Non-Homology Filter: Excludes host cross-reactive proteins (UniProt UP000005640).
3. Gut Microbiome Preservation Filter: Identifies commensal-sparing proteins (HRGM2 2026).

Academic Citations:
- MMseqs2: Steinegger M & Söding J. Nature Biotechnology (2017) 35:1026–1028; Nat Commun (2018) 9:1629.
- DEG: Luo H et al. Nucleic Acids Res (2014) 42:D574–D580; (2021) 49:D677–D684.
- Human Proteome: UniProt Consortium. Nucleic Acids Res (2025) 53:D501–D517.
- HRGM2: Ma J et al. Nature Microbiology (2026) 11:317–334 (PMID: 41345261).

Usage:
    python scripts/04_subtractive_filter.py [--config config.yaml]
"""

import os
import sys
import time
import gzip
import json
import shutil
import logging
import argparse
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, Set
import pandas as pd
import yaml

try:
    from scripts.utils_hardware import resolve_threads
except ImportError:
    from utils_hardware import resolve_threads

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("04_subtractive_filter")


def ensure_mmseqs_binary() -> str:
    """Finds MMseqs2 in PATH or project bin/ directory."""
    which_mmseqs = shutil.which("mmseqs")
    if which_mmseqs:
        return which_mmseqs

    project_root = Path(__file__).resolve().parent.parent
    local_bin = project_root / "bin" / "mmseqs"
    if local_bin.exists() and os.access(local_bin, os.X_OK):
        return str(local_bin)

    raise RuntimeError(
        "MMseqs2 binary not found. Please install mmseqs or ensure bin/mmseqs exists."
    )


def download_file_with_progress(url: str, output_path: Path, description: str):
    """Downloads a file from a URL with retry logic and progress reporting."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    logger.info("Downloading %s from %s...", description, url)

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Pangenome-Drug Pipeline)"})
    with urllib.request.urlopen(req, timeout=120) as resp, open(tmp_path, "wb") as out_f:
        shutil.copyfileobj(resp, out_f)

    tmp_path.rename(output_path)
    logger.info("Successfully downloaded %s (%s bytes).", description, f"{output_path.stat().st_size:,}")


def ensure_fasta_database(
    fasta_path: Path,
    download_url: Optional[str] = None,
    description: str = "",
    is_gzipped: bool = False,
) -> Path:
    """
    Ensures that a reference FASTA database is present and uncompressed.
    If missing, downloads and decompresses it automatically.
    """
    if fasta_path.exists() and fasta_path.stat().st_size > 0:
        logger.info("Found %s database at: %s (%s bytes)", description, fasta_path, f"{fasta_path.stat().st_size:,}")
        return fasta_path

    if not download_url:
        raise FileNotFoundError(
            f"Reference database {description} not found at {fasta_path} and no download URL was provided."
        )

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    if is_gzipped or download_url.endswith(".gz"):
        gz_path = fasta_path.with_suffix(fasta_path.suffix + ".gz")
        if not (gz_path.exists() and gz_path.stat().st_size > 0):
            download_file_with_progress(download_url, gz_path, f"{description} (compressed)")
        logger.info("Decompressing %s to %s...", gz_path.name, fasta_path.name)
        with gzip.open(gz_path, "rb") as f_in, open(fasta_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logger.info("Decompressed %s (%s bytes).", fasta_path.name, f"{fasta_path.stat().st_size:,}")
    else:
        download_file_with_progress(download_url, fasta_path, description)

    return fasta_path


def run_mmseqs_search(
    mmseqs_bin: str,
    query_fasta: Path,
    target_fasta: Path,
    output_tsv: Path,
    threads: int = 4,
    max_evalue: float = 1.0e-4,
    max_seqs: Optional[int] = None,
) -> Path:
    """
    Runs mmseqs easy-search query against target database with standard tabular output.
    """
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_tsv.parent / f"tmp_{output_tsv.stem}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    format_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov"
    )

    cmd = [
        mmseqs_bin,
        "easy-search",
        str(query_fasta),
        str(target_fasta),
        str(output_tsv),
        str(tmp_dir),
        "--threads",
        str(threads),
        "-e",
        str(max_evalue),
        "--format-output",
        format_fields,
        "-v",
        "1",
    ]

    if max_seqs is not None and int(max_seqs) > 0:
        cmd.extend(["--max-seqs", str(int(max_seqs))])

    logger.info("Executing MMseqs2 search: %s -> %s", query_fasta.name, target_fasta.name)
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    # Cleanup temp directory
    shutil.rmtree(tmp_dir, ignore_errors=True)

    if res.returncode != 0:
        logger.error("MMseqs2 search failed:\n%s\n%s", res.stdout, res.stderr)
        raise RuntimeError(f"MMseqs2 search failed with code {res.returncode}")

    logger.info("MMseqs2 search completed in %.2f seconds (saved %s).", elapsed, output_tsv.name)
    return output_tsv


def parse_significant_hits(
    search_tsv: Path,
    max_evalue: float,
    min_pident: float,
    min_bitscore: float,
    min_query_cov: float,
) -> Dict[str, Dict[str, Any]]:
    """
    Parses MMseqs2 tabular output and returns the best significant hit per query.
    Criteria:
      evalue <= max_evalue
      fident >= min_pident / 100.0
      bits >= min_bitscore
      qcov >= min_query_cov
    """
    best_hits = {}
    if not search_tsv.exists() or search_tsv.stat().st_size == 0:
        return best_hits

    cols = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen", "qcov"
    ]

    try:
        df = pd.read_csv(search_tsv, sep="\t", names=cols, low_memory=False)
    except pd.errors.EmptyDataError:
        return best_hits

    # Numeric conversion
    numeric_cols = ["fident", "alnlen", "evalue", "bits", "qlen", "tlen", "qcov"]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Filter significant hits
    sig_mask = (
        (df["evalue"] <= max_evalue) &
        (df["fident"] >= (min_pident / 100.0)) &
        (df["bits"] >= min_bitscore) &
        (df["qcov"] >= min_query_cov)
    )
    df_sig = df[sig_mask].copy()

    # Sort by bitscore descending, evalue ascending to pick best hit
    df_sig = df_sig.sort_values(by=["bits", "evalue"], ascending=[False, True])

    for _, row in df_sig.iterrows():
        qid = str(row["query"])
        if qid not in best_hits:
            best_hits[qid] = {
                "target": str(row["target"]),
                "evalue": float(row["evalue"]),
                "pident": float(row["fident"]) * 100.0,
                "bitscore": float(row["bits"]),
                "qcov": float(row["qcov"]) * 100.0,
            }

    return best_hits


def parse_gut_microbiome_hits(
    search_tsv: Path,
    max_evalue: float,
    min_pident: float,
    min_bitscore: float,
    min_query_cov: float,
    mode: str = "percentage",
    max_prevalence_pct: float = 25.0,
    max_hit_count: int = 50,
    total_gut_species: int = 4824,
) -> Dict[str, Dict[str, Any]]:
    """
    Parses MMseqs2 tabular output for gut microbiome preservation screening.
    Calculates the number of unique gut genomes matched and prevalence across gut species.
    
    A candidate protein is marked as 'is_homologous = True' (Tier 2 / broad-spectrum)
    only if its prevalence meets or exceeds the configured threshold.
    Otherwise, it is marked as 'is_homologous = False' (Gut-Safe / Tier 1 eligible).
    """
    gut_results = {}
    if not search_tsv.exists() or search_tsv.stat().st_size == 0:
        return gut_results

    cols = [
        "query", "target", "fident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen", "qcov"
    ]

    try:
        df = pd.read_csv(search_tsv, sep="\t", names=cols, low_memory=False)
    except pd.errors.EmptyDataError:
        return gut_results

    numeric_cols = ["fident", "alnlen", "evalue", "bits", "qlen", "tlen", "qcov"]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Filter significant hits based on stringency cutoffs
    sig_mask = (
        (df["evalue"] <= max_evalue) &
        (df["fident"] >= (min_pident / 100.0)) &
        (df["bits"] >= min_bitscore) &
        (df["qcov"] >= min_query_cov)
    )
    df_sig = df[sig_mask].copy()

    if df_sig.empty:
        return gut_results

    # Extract genome ID (e.g. 'GENOME124745' from 'GENOME124745_3603')
    df_sig["genome_id"] = df_sig["target"].astype(str).str.split("_").str[0]

    grouped = df_sig.groupby("query")

    for qid, group in grouped:
        unique_genomes = set(group["genome_id"].dropna().unique())
        unique_genomes_count = len(unique_genomes)

        denom = max(int(total_gut_species), 1)
        prevalence_pct = (unique_genomes_count / float(denom)) * 100.0

        # Determine homology status based on mode
        if mode == "percentage":
            is_homologous = (prevalence_pct >= float(max_prevalence_pct))
        elif mode == "count":
            is_homologous = (unique_genomes_count >= int(max_hit_count))
        elif mode == "strict":
            is_homologous = (unique_genomes_count >= 1)
        else:
            is_homologous = (prevalence_pct >= float(max_prevalence_pct))

        best_row = group.sort_values(by=["bits", "evalue"], ascending=[False, True]).iloc[0]

        gut_results[str(qid)] = {
            "target": str(best_row["target"]),
            "evalue": float(best_row["evalue"]),
            "pident": float(best_row["fident"]) * 100.0,
            "bitscore": float(best_row["bits"]),
            "qcov": float(best_row["qcov"]) * 100.0,
            "hit_count": unique_genomes_count,
            "prevalence_pct": round(prevalence_pct, 2),
            "is_homologous": is_homologous,
        }

    return gut_results


def resolve_total_gut_species(gut_cfg: dict, fasta_path: Path) -> Tuple[int, str]:
    """
    Dynamically determines the total number of species/genomes in the gut reference database.
    If 'total_gut_species' is an integer, returns it.
    If 'auto' (or omitted):
      1. Checks for metadata tables (e.g. HRGMv2_Cluster_metadata.tsv) in the database directory.
      2. Checks for a cached .species_count sidecar file.
      3. Performs a dynamic scan of FASTA headers (identifying unique taxon/genome IDs).
    """
    val = gut_cfg.get("total_gut_species", "auto")
    if val is not None and str(val).lower() != "auto":
        try:
            return int(val), f"user configuration ({val})"
        except ValueError:
            pass

    # 1. Check for metadata table in database directory
    db_dir = fasta_path.parent
    for meta_name in [
        "HRGMv2_Cluster_metadata.tsv",
        "cluster_metadata.tsv",
        "species_metadata.tsv",
        "metadata.tsv",
    ]:
        meta_file = db_dir / meta_name
        if meta_file.exists() and meta_file.stat().st_size > 0:
            try:
                df_meta = pd.read_csv(meta_file, sep="\t")
                for col in ["HRGMv2 Cluster", "Cluster", "Species", "species", "Representative Genome", "genome_id"]:
                    if col in df_meta.columns:
                        n_unique = int(df_meta[col].nunique())
                        return n_unique, f"metadata table {meta_file.name} ({n_unique:,} species/clusters)"
                return len(df_meta), f"metadata table {meta_file.name} ({len(df_meta):,} rows)"
            except Exception as e:
                logger.warning("Could not read %s for species count: %s", meta_file.name, e)

    # 2. Check cached count file
    cache_file = fasta_path.with_suffix(".species_count")
    if cache_file.exists():
        try:
            n_cached = int(cache_file.read_text().strip())
            return n_cached, f"cached index {cache_file.name} ({n_cached:,} species)"
        except Exception:
            pass

    # 3. Dynamic header scan of the database FASTA
    logger.info("Scanning FASTA headers in %s to dynamically detect total gut species...", fasta_path.name)
    species_set = set()
    with open(fasta_path, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                hdr = line[1:].strip()
                if "|" in hdr:
                    species_set.add(hdr.split("|")[0].strip())
                elif "_" in hdr and hdr.startswith("GENOME"):
                    species_set.add(hdr.split("_")[0].strip())
                else:
                    species_set.add(hdr.split()[0].strip())

    total_detected = max(len(species_set), 1)
    try:
        cache_file.write_text(str(total_detected))
    except Exception:
        pass

    return total_detected, f"dynamic FASTA header scan ({total_detected:,} unique species/genomes)"


def parse_fasta_headers(fasta_path: Path) -> Dict[str, Dict[str, str]]:
    """
    Parses FASTA file headers to extract sequence IDs, gene names, descriptions, and full sequences.
    """
    records = {}
    current_id = None
    current_header = ""
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    records[current_id] = {
                        "header": current_header,
                        "seq": "".join(current_seq),
                    }
                current_header = line[1:]
                current_id = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            records[current_id] = {
                "header": current_header,
                "seq": "".join(current_seq),
            }

    return records


def main():
    parser = argparse.ArgumentParser(description="Step 4: Subtractive Genomics & Target Prioritization")
    parser.add_argument("--config", default="config.yaml", help="Path to pipeline config YAML")
    args = parser.parse_args()

    # Load configuration
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    job_name = config.get("job_name", "default_job")
    results_dir = Path(config.get("results_dir", "results"))
    job_dir = results_dir / job_name

    step4_cfg = config.get("step4_subtractive_filter", {})
    if not step4_cfg.get("enabled", True):
        logger.info("Step 4 (Subtractive Filter) is disabled in config. Exiting.")
        return

    databases_dir = Path(step4_cfg.get("databases_dir", "databases"))
    databases_dir.mkdir(parents=True, exist_ok=True)
    raw_threads = step4_cfg.get("threads", "auto")
    threads = resolve_threads(raw_threads, config=config, default=32, task_name="subtractive_mmseqs_searches")

    # Thresholds
    max_evalue = float(step4_cfg.get("max_evalue", 1.0e-5))
    min_pident = float(step4_cfg.get("min_pident", 35.0))
    min_bitscore = float(step4_cfg.get("min_bitscore", 100.0))
    min_query_cov = float(step4_cfg.get("min_query_cov", 0.35))

    # Toggles
    enable_deg = bool(step4_cfg.get("enable_deg_essentiality_filter", True))
    enable_human = bool(step4_cfg.get("enable_human_non_homology_filter", True))
    enable_gut = bool(step4_cfg.get("enable_gut_microbiome_filter", True))

    logger.info("=" * 70)
    logger.info("Step 4: Subtractive Genomics Filter & Target Prioritization")
    logger.info("Job: %s", job_name)
    logger.info("Hardware Threads: %d allocated for MMseqs2 searches", threads)
    logger.info("Active Filters: DEG Essentiality=%s, Human Non-Homology=%s, HRGM2 Gut Microbiome=%s",
                enable_deg, enable_human, enable_gut)
    logger.info("Thresholds: E-value <= %g, Identity >= %.1f%%, Bitscore >= %.1f, Coverage >= %.1f%%",
                max_evalue, min_pident, min_bitscore, min_query_cov * 100)
    logger.info("=" * 70)

    # Input: Core Representatives FASTA from Step 3
    core_fasta = job_dir / "03_pangenome_analysis" / "fasta" / "core_representatives.fasta"
    if not core_fasta.exists():
        raise FileNotFoundError(
            f"Core representatives FASTA not found at {core_fasta}. Please run Step 3 first."
        )

    queries = parse_fasta_headers(core_fasta)
    total_core = len(queries)
    logger.info("Loaded %d Core representative proteins from %s.", total_core, core_fasta.name)

    mmseqs_bin = ensure_mmseqs_binary()
    logger.info("Using MMseqs2 binary: %s", mmseqs_bin)

    step4_dir = job_dir / "04_subtractive_filter"
    search_dir = step4_dir / "searches"
    tables_dir = step4_dir / "tables"
    fasta_dir = step4_dir / "fasta"
    summary_dir = step4_dir / "summary"

    for d in [search_dir, tables_dir, fasta_dir, summary_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------------------------
    # 1. Host (Human) Filter
    # --------------------------------------------------------------------------
    human_hits = {}
    if enable_human:
        human_cfg = step4_cfg.get("human_filter", {})
        human_path = Path(human_cfg.get("fasta_path", databases_dir / "human" / "human_UP000005640.fasta"))
        human_url = human_cfg.get(
            "url",
            f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3A{human_cfg.get('proteome_id', 'UP000005640')}%29"
        )
        ensure_fasta_database(human_path, human_url, "Human Proteome (UP000005640)", is_gzipped=True)
        human_search_tsv = search_dir / "human_search_results.tsv"
        run_mmseqs_search(mmseqs_bin, core_fasta, human_path, human_search_tsv, threads, max_evalue)
        human_hits = parse_significant_hits(human_search_tsv, max_evalue, min_pident, min_bitscore, min_query_cov)
        logger.info("Identified %d / %d core proteins with significant Human homology (Excluded).",
                    len(human_hits), total_core)

    # --------------------------------------------------------------------------
    # 2. DEG Essentiality Filter
    # --------------------------------------------------------------------------
    deg_hits = {}
    if enable_deg:
        deg_cfg = step4_cfg.get("deg_filter", {})
        deg_path = Path(deg_cfg.get("fasta_path", databases_dir / "deg" / "deg_bacteria.fasta"))
        deg_url = deg_cfg.get("deg_url", "http://tubic.org/deg/public/download/DEG10.aa.gz")
        ensure_fasta_database(deg_path, deg_url, "Database of Essential Genes (DEG)", is_gzipped=True)
        deg_search_tsv = search_dir / "deg_search_results.tsv"
        run_mmseqs_search(mmseqs_bin, core_fasta, deg_path, deg_search_tsv, threads, max_evalue)
        deg_hits = parse_significant_hits(deg_search_tsv, max_evalue, min_pident, min_bitscore, min_query_cov)
        logger.info("Identified %d / %d core proteins with significant DEG Essentiality matches.",
                    len(deg_hits), total_core)

    # --------------------------------------------------------------------------
    # 3. Commensal Gut Microbiome Filter (HRGM2)
    # --------------------------------------------------------------------------
    gut_hits = {}
    gut_homologous_count = 0
    gut_safe_count = total_core
    gut_mode = "percentage"
    max_prevalence_pct = 25.0
    max_hit_count = 50

    if enable_gut:
        gut_cfg = step4_cfg.get("gut_microbiome_filter", {})
        gut_path = Path(gut_cfg.get("fasta_path", databases_dir / "gut_microbiota" / "hrgm2_50_rep.fasta"))
        hrgm2_url = gut_cfg.get(
            "hrgm2_url",
            "https://www.decodebiome.org/HRGM/data/protein_catalog/0.HRGMv2_Proteins/6.HRGMv2_50_Proteins/HRGMv2_50_Proteins_rep_seq.faa.gz"
        )
        ensure_fasta_database(gut_path, hrgg2_url if (hrgg2_url := hrgm2_url) else None, "HRGM2 50% Representative Catalogue (Nature Micro 2026)", is_gzipped=True)
        gut_search_tsv = search_dir / "gut_search_results.tsv"
        max_seqs = gut_cfg.get("max_seqs", 5000)
        run_mmseqs_search(mmseqs_bin, core_fasta, gut_path, gut_search_tsv, threads, max_evalue, max_seqs=max_seqs)

        meta_path = databases_dir / "gut_microbiota" / "HRGMv2_Cluster_metadata.tsv"
        meta_url = gut_cfg.get(
            "metadata_url",
            "https://zenodo.org/api/records/19480672/files/HRGMv2_Cluster_metadata.tsv/content",
        )
        if not meta_path.exists() and "hrgm2" in str(gut_path).lower():
            try:
                download_file_with_progress(meta_url, meta_path, "HRGM2 Species Metadata (Zenodo)")
            except Exception as e:
                logger.warning("Could not download HRGM2 metadata: %s. Falling back to dynamic FASTA scan.", e)

        gut_mode = gut_cfg.get("mode", "percentage")
        max_prevalence_pct = float(gut_cfg.get("max_prevalence_pct", 25.0))
        max_hit_count = int(gut_cfg.get("max_hit_count", 50))
        total_gut_species, count_source = resolve_total_gut_species(gut_cfg, gut_path)
        logger.info("Total reference gut species resolved: %d (via %s)", total_gut_species, count_source)

        gut_hits = parse_gut_microbiome_hits(
            gut_search_tsv,
            max_evalue=max_evalue,
            min_pident=min_pident,
            min_bitscore=min_bitscore,
            min_query_cov=min_query_cov,
            mode=gut_mode,
            max_prevalence_pct=max_prevalence_pct,
            max_hit_count=max_hit_count,
            total_gut_species=total_gut_species,
        )

        gut_homologous_count = sum(1 for h in gut_hits.values() if h["is_homologous"])
        gut_safe_count = total_core - gut_homologous_count

        thresh_str = f"{max_prevalence_pct}% species prevalence" if gut_mode == "percentage" else f"{max_hit_count} genomes"
        logger.info(
            "Gut Microbiome Screening (Mode: %s | Threshold: >= %s): %d / %d core proteins Homologous (Shared across gut), %d Gut-Safe.",
            gut_mode,
            thresh_str,
            gut_homologous_count,
            total_core,
            gut_safe_count,
        )

    # --------------------------------------------------------------------------
    # 4. Integrate Annotations & Assign Target Tiers
    # --------------------------------------------------------------------------
    rows = []
    tier_counts = {
        "Tier 1: Pathogen-Selective (Microbiome-Safe)": 0,
        "Tier 2: Broad-Spectrum (Shared with Gut)": 0,
        "Tier 3: Non-Essential (Safe)": 0,
        "Excluded: Host Cross-Reactive (Cytotoxic)": 0,
    }

    tier1_records = []
    tier2_records = []
    essential_non_human_records = []
    non_human_core_records = []

    for qid, record in queries.items():
        header = record["header"]
        seq = record["seq"]

        parts = header.split(maxsplit=1)
        first_token = parts[0]
        desc = parts[1] if len(parts) > 1 else ""

        tokens = first_token.split("|")
        cluster_id = tokens[0] if len(tokens) >= 1 else qid
        rep_accession = tokens[2] if len(tokens) >= 3 else (tokens[1] if len(tokens) >= 2 else qid)

        gene_name = "N/A"
        for item in desc.split():
            if item.startswith("GN="):
                gene_name = item.split("=")[1]
                break

        is_human = qid in human_hits
        h_info = human_hits.get(qid, {})

        is_essential = qid in deg_hits
        d_info = deg_hits.get(qid, {})

        g_info = gut_hits.get(qid, {})
        is_gut_homologous = g_info.get("is_homologous", False)
        gut_hit_count = g_info.get("hit_count", 0)
        gut_prev_pct = g_info.get("prevalence_pct", 0.0)

        # Tier classification:
        # Tier 1 requires: Essential (DEG) + Non-Human + Gut-Safe (< threshold prevalence)
        # Tier 2: Essential (DEG) + Non-Human + Gut-Homologous (>= threshold prevalence)
        if enable_human and is_human:
            tier = "Excluded: Host Cross-Reactive (Cytotoxic)"
        elif enable_deg and is_essential:
            if enable_gut:
                if not is_gut_homologous:
                    tier = "Tier 1: Pathogen-Selective (Microbiome-Safe)"
                else:
                    tier = "Tier 2: Broad-Spectrum (Shared with Gut)"
            else:
                tier = "Tier 1: Pathogen-Selective (Microbiome-Safe)"
        else:
            tier = "Tier 3: Non-Essential (Safe)"

        tier_counts[tier] = tier_counts.get(tier, 0) + 1

        fasta_entry = f">{header} [Tier={tier}]\n{seq}\n"
        if tier == "Tier 1: Pathogen-Selective (Microbiome-Safe)":
            tier1_records.append(fasta_entry)
            essential_non_human_records.append(fasta_entry)
            non_human_core_records.append(fasta_entry)
        elif tier == "Tier 2: Broad-Spectrum (Shared with Gut)":
            tier2_records.append(fasta_entry)
            essential_non_human_records.append(fasta_entry)
            non_human_core_records.append(fasta_entry)
        elif tier == "Tier 3: Non-Essential (Safe)":
            non_human_core_records.append(fasta_entry)

        rows.append({
            "Cluster_ID": cluster_id,
            "Query_ID": qid,
            "Representative_Accession": rep_accession,
            "Gene_Name": gene_name,
            "Description": desc,
            "Sequence_Length": len(seq),
            "Target_Tier": tier,
            "DEG_Essential": is_essential,
            "DEG_Hit_ID": d_info.get("target", "None"),
            "DEG_Evalue": d_info.get("evalue", None),
            "DEG_Identity_Pct": d_info.get("pident", None),
            "DEG_Bitscore": d_info.get("bitscore", None),
            "Human_Homolog": is_human,
            "Human_Hit_ID": h_info.get("target", "None"),
            "Human_Evalue": h_info.get("evalue", None),
            "Human_Identity_Pct": h_info.get("pident", None),
            "Human_Bitscore": h_info.get("bitscore", None),
            "Gut_Microbiome_Homolog": is_gut_homologous,
            "Gut_Hit_Count": gut_hit_count,
            "Gut_Prevalence_Pct": gut_prev_pct,
            "Gut_Hit_ID": g_info.get("target", "None"),
            "Gut_Evalue": g_info.get("evalue", None),
            "Gut_Identity_Pct": g_info.get("pident", None),
            "Gut_Bitscore": g_info.get("bitscore", None),
        })

    df_master = pd.DataFrame(rows)

    # --------------------------------------------------------------------------
    # 5. Export Master and Tiered Tables
    # --------------------------------------------------------------------------
    master_tsv = tables_dir / "subtractive_master.tsv"
    master_csv = tables_dir / "subtractive_master.csv"
    df_master.to_csv(master_tsv, sep="\t", index=False)
    df_master.to_csv(master_csv, index=False)
    logger.info("Saved Subtractive Master Table (%d rows) to %s and %s.", len(df_master), master_tsv.name, master_csv.name)

    df_tier1 = df_master[df_master["Target_Tier"] == "Tier 1: Pathogen-Selective (Microbiome-Safe)"].copy()
    tier1_tsv = tables_dir / "tier1_golden_targets.tsv"
    tier1_csv = tables_dir / "tier1_golden_targets.csv"
    df_tier1.to_csv(tier1_tsv, sep="\t", index=False)
    df_tier1.to_csv(tier1_csv, index=False)
    logger.info("Saved Tier 1 Pathogen-Selective Targets (%d rows) to %s.", len(df_tier1), tier1_tsv.name)

    df_tier2 = df_master[df_master["Target_Tier"] == "Tier 2: Broad-Spectrum (Shared with Gut)"].copy()
    tier2_tsv = tables_dir / "tier2_broad_spectrum.tsv"
    tier2_csv = tables_dir / "tier2_broad_spectrum.csv"
    df_tier2.to_csv(tier2_tsv, sep="\t", index=False)
    df_tier2.to_csv(tier2_csv, index=False)

    # --------------------------------------------------------------------------
    # 6. Export Filtered FASTA Files
    # --------------------------------------------------------------------------
    tier1_fasta = fasta_dir / "tier1_golden_targets.fasta"
    with open(tier1_fasta, "w") as f:
        f.writelines(tier1_records)
    logger.info("Saved %d Tier 1 FASTA sequences to %s.", len(tier1_records), tier1_fasta.name)

    tier2_fasta = fasta_dir / "tier2_broad_spectrum.fasta"
    with open(tier2_fasta, "w") as f:
        f.writelines(tier2_records)

    essential_fasta = fasta_dir / "essential_non_human.fasta"
    with open(essential_fasta, "w") as f:
        f.writelines(essential_non_human_records)
    logger.info("Saved %d Essential Non-Human FASTA sequences to %s.", len(essential_non_human_records), essential_fasta.name)

    non_human_core_fasta = fasta_dir / "non_human_core.fasta"
    with open(non_human_core_fasta, "w") as f:
        f.writelines(non_human_core_records)

    # --------------------------------------------------------------------------
    # 7. Summary Metrics
    # --------------------------------------------------------------------------
    summary_data = {
        "job_name": job_name,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_core_proteins": total_core,
        "filters_enabled": {
            "deg_essentiality": enable_deg,
            "human_non_homology": enable_human,
            "gut_microbiome": enable_gut,
        },
        "thresholds": {
            "max_evalue": max_evalue,
            "min_pident": min_pident,
            "min_bitscore": min_bitscore,
            "min_query_cov": min_query_cov,
            "gut_filter_mode": gut_mode,
            "gut_prevalence_pct": max_prevalence_pct if gut_mode == "percentage" else None,
            "gut_hit_count_threshold": max_hit_count if gut_mode == "count" else None,
        },
        "filtering_metrics": {
            "human_homologs_count": len(human_hits),
            "human_non_homologous_count": total_core - len(human_hits),
            "deg_essential_count": len(deg_hits),
            "gut_microbiome_homologs_count": gut_homologous_count,
            "gut_microbiome_safe_count": gut_safe_count,
            "tier_counts": tier_counts,
        },
        "output_files": {
            "master_tsv": str(master_tsv),
            "tier1_tsv": str(tier1_tsv),
            "tier1_fasta": str(tier1_fasta),
            "essential_non_human_fasta": str(essential_fasta),
        },
    }

    summary_json = summary_dir / "subtractive_summary.json"
    with open(summary_json, "w") as f:
        json.dump(summary_data, f, indent=2)

    summary_tsv = summary_dir / "subtractive_summary.tsv"
    gut_label = (
        f"HRGM2 Gut Homologs (Prevalence >= {max_prevalence_pct}%)"
        if gut_mode == "percentage"
        else f"HRGM2 Gut Homologs (Hit Count >= {max_hit_count})"
    )
    summary_rows = [
        {"Metric": "Total Core Pangenome Proteins", "Count": total_core, "Percentage": 100.0},
        {"Metric": "Human Homologs (Excluded)", "Count": len(human_hits), "Percentage": (len(human_hits) / total_core * 100)},
        {"Metric": "Human Non-Homologous (Safe)", "Count": total_core - len(human_hits), "Percentage": ((total_core - len(human_hits)) / total_core * 100)},
        {"Metric": "DEG Essential Genes", "Count": len(deg_hits), "Percentage": (len(deg_hits) / total_core * 100)},
        {"Metric": gut_label, "Count": gut_homologous_count, "Percentage": (gut_homologous_count / total_core * 100)},
        {"Metric": "HRGM2 Gut Microbiome Safe (< threshold)", "Count": gut_safe_count, "Percentage": (gut_safe_count / total_core * 100)},
        {"Metric": "Tier 1: Pathogen-Selective (Microbiome-Safe)", "Count": tier_counts["Tier 1: Pathogen-Selective (Microbiome-Safe)"], "Percentage": (tier_counts["Tier 1: Pathogen-Selective (Microbiome-Safe)"] / total_core * 100)},
        {"Metric": "Tier 2: Broad-Spectrum (Shared with Gut)", "Count": tier_counts["Tier 2: Broad-Spectrum (Shared with Gut)"], "Percentage": (tier_counts["Tier 2: Broad-Spectrum (Shared with Gut)"] / total_core * 100)},
        {"Metric": "Tier 3: Non-Essential (Safe)", "Count": tier_counts["Tier 3: Non-Essential (Safe)"], "Percentage": (tier_counts["Tier 3: Non-Essential (Safe)"] / total_core * 100)},
    ]
    pd.DataFrame(summary_rows).to_csv(summary_tsv, sep="\t", index=False)
    logger.info("Saved Subtractive Filter Summary to %s.", summary_json.name)

    sentinel = step4_dir / ".subtractive_complete"
    sentinel.touch()

    logger.info("=" * 70)
    logger.info("Step 4 Subtractive Genomics successfully finished!")
    logger.info("Tier 1 Pathogen-Selective: %d (%.1f%% of Core)",
                tier_counts["Tier 1: Pathogen-Selective (Microbiome-Safe)"],
                tier_counts["Tier 1: Pathogen-Selective (Microbiome-Safe)"] / total_core * 100)
    logger.info("Tier 2 Broad-Spectrum: %d (%.1f%% of Core)",
                tier_counts["Tier 2: Broad-Spectrum (Shared with Gut)"],
                tier_counts["Tier 2: Broad-Spectrum (Shared with Gut)"] / total_core * 100)
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
