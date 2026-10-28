#!/usr/bin/env python3
"""
Step 6: Virulence Factor (VFDB) & Antimicrobial Resistance (CARD) Screening
==========================================================================
Screens prioritized drug targets against:
1. Virulence Factor Database (VFDB 2022/2026 Core Dataset):
   Identifies dual-action super-targets (bactericidal + anti-virulence) and
   antigenic/anti-virulence targets.
2. Comprehensive Antibiotic Resistance Database (CARD v3.3+ ARO Models):
   Identifies antibiotic resistance determinants, candidate resistance-breakers /
   adjuvants (e.g. efflux pumps, beta-lactamases), and de-risks targets from
   pre-existing clinical resistance alleles.

Academic Citations:
- VFDB: Liu B et al. Nucleic Acids Res (2022) 50:D704–D707 (PMID: 34850947).
- CARD: Alcock BP et al. Nucleic Acids Res (2023) 51:D690–D699 (PMID: 36263822).
- MMseqs2: Steinegger M & Söding J. Nature Biotechnology (2017) 35:1026–1028.
- Anti-Virulence: Dickey SW et al. Nat Rev Drug Discov (2017) 16:457–471.
- Adjuvants: Wright GD. Trends Microbiol (2016) 24:862–871.

Usage:
    python scripts/06_virulence_resistance_screening.py [--config config.yaml]
"""

import os
import re
import sys
import time
import json
import gzip
import tarfile
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
logger = logging.getLogger("06_virulence_resistance")


def ensure_mmseqs_binary() -> str:
    """Finds MMseqs2 binary in PATH or local project bin directory."""
    which_mmseqs = shutil.which("mmseqs")
    if which_mmseqs:
        return which_mmseqs

    project_root = Path(__file__).resolve().parent.parent
    local_bin = project_root / "bin" / "mmseqs"
    if local_bin.exists() and os.access(local_bin, os.X_OK):
        return str(local_bin)

    raise RuntimeError(
        "MMseqs2 binary not found. Please install MMseqs2 or ensure bin/mmseqs exists."
    )


def download_file_with_retry(
    url: str,
    output_path: Path,
    description: str,
    retries: int = 3,
    timeout: int = 180,
) -> Path:
    """Downloads a remote file with retries and temporary file swapping."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info("Found existing %s at: %s (%s bytes)", description, output_path, f"{output_path.stat().st_size:,}")
        return output_path

    tmp_path = output_path.with_suffix(output_path.suffix + ".download_tmp")
    logger.info("Downloading %s from %s...", description, url)

    headers = {"User-Agent": "Mozilla/5.0 (Pangenome-Drug Pipeline/1.0; Research)"}
    req = urllib.request.Request(url, headers=headers)

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp_path, "wb") as out_f:
                shutil.copyfileobj(resp, out_f)
            tmp_path.rename(output_path)
            logger.info("Successfully downloaded %s (%s bytes).", description, f"{output_path.stat().st_size:,}")
            return output_path
        except Exception as e:
            last_err = e
            logger.warning("Attempt %d/%d failed to download %s: %s", attempt, retries, description, e)
            if tmp_path.exists():
                tmp_path.unlink()
            time.sleep(attempt * 2)

    raise RuntimeError(f"Failed to download {description} after {retries} attempts: {last_err}")


def ensure_vfdb_database(vfdb_cfg: Dict[str, Any], databases_dir: Path) -> Path:
    """
    Ensures that the VFDB protein database is downloaded and uncompressed.
    Returns path to uncompressed FASTA.
    """
    default_url = "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
    url = vfdb_cfg.get("url", default_url)
    fasta_path_cfg = vfdb_cfg.get("fasta_path")
    
    if fasta_path_cfg:
        fasta_path = Path(fasta_path_cfg)
        if not fasta_path.is_absolute():
            fasta_path = Path.cwd() / fasta_path
    else:
        fasta_path = databases_dir / "vfdb" / "VFDB_setA_pro.fasta"

    if fasta_path.exists() and fasta_path.stat().st_size > 0:
        logger.info("Found VFDB database at: %s (%s bytes)", fasta_path, f"{fasta_path.stat().st_size:,}")
        return fasta_path

    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    gz_path = fasta_path.with_suffix(fasta_path.suffix + ".gz")

    download_file_with_retry(url, gz_path, "VFDB Core dataset")

    logger.info("Decompressing VFDB archive: %s -> %s", gz_path.name, fasta_path.name)
    with gzip.open(gz_path, "rb") as gz_in, open(fasta_path, "wb") as f_out:
        shutil.copyfileobj(gz_in, f_out)

    logger.info("VFDB database decompressed successfully (%s bytes).", f"{fasta_path.stat().st_size:,}")
    return fasta_path


def ensure_card_database(card_cfg: Dict[str, Any], databases_dir: Path) -> Tuple[Path, Path]:
    """
    Ensures that the CARD database archive is downloaded and extracted.
    Returns (card_fasta_path, aro_index_path).
    """
    default_url = "https://card.mcmaster.ca/latest/data"
    url = card_cfg.get("url", default_url)
    database_dir_cfg = card_cfg.get("database_dir")
    
    if database_dir_cfg:
        card_dir = Path(database_dir_cfg)
        if not card_dir.is_absolute():
            card_dir = Path.cwd() / card_dir
    else:
        card_dir = databases_dir / "card"

    card_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = card_dir / "protein_fasta_protein_homolog_model.fasta"
    aro_index_path = card_dir / "aro_index.tsv"

    if fasta_path.exists() and fasta_path.stat().st_size > 0 and aro_index_path.exists() and aro_index_path.stat().st_size > 0:
        logger.info("Found CARD database at: %s and aro_index at: %s", fasta_path, aro_index_path)
        return fasta_path, aro_index_path

    tar_path = card_dir / "card_data.tar.bz2"
    download_file_with_retry(url, tar_path, "CARD latest data archive")

    logger.info("Extracting CARD archive: %s -> %s", tar_path.name, card_dir)
    with tarfile.open(tar_path, "r:bz2") as tar:
        for member in tar.getmembers():
            base_name = Path(member.name).name
            if base_name in ("protein_fasta_protein_homolog_model.fasta", "aro_index.tsv", "aro_categories.tsv"):
                member_file = tar.extractfile(member)
                if member_file:
                    target_dest = card_dir / base_name
                    with open(target_dest, "wb") as out_f:
                        shutil.copyfileobj(member_file, out_f)
                    logger.info("Extracted CARD member: %s (%s bytes)", base_name, f"{target_dest.stat().st_size:,}")

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        raise RuntimeError(f"CARD FASTA {fasta_path} could not be extracted from {tar_path}")

    return fasta_path, aro_index_path


def parse_vfdb_metadata(vfdb_fasta: Path) -> Dict[str, Dict[str, str]]:
    """
    Parses headers from VFDB FASTA file into structured metadata.
    Example header:
    >VFG037176(gb|WP_001081735) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]
    """
    logger.info("Parsing VFDB headers from %s...", vfdb_fasta.name)
    metadata = {}
    
    with open(vfdb_fasta, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            first_token = header.split()[0]
            
            # Extract VFG id
            vfg_match = re.search(r"^(VFG\d+)", first_token)
            vfg_id = vfg_match.group(1) if vfg_match else first_token
            
            # Extract gene symbol (e.g. (plc1))
            gene_match = re.search(r"\(([a-zA-Z0-9_\-]+)\)\s+[a-zA-Z]", header)
            gene_name = gene_match.group(1) if gene_match else ""
            
            # Extract brackets
            brackets = re.findall(r"\[(.*?)\]", header)
            organism = brackets[-1] if brackets else ""
            
            vf_details = brackets[0] if len(brackets) >= 2 else (brackets[0] if brackets else "")
            
            vf_name = vf_details
            vf_class = "Unclassified Virulence Factor"
            if " - " in vf_details:
                parts = vf_details.split(" - ", 1)
                vf_name = parts[0].strip()
                vf_class = parts[1].strip()
            elif vf_details:
                vf_name = vf_details
                vf_class = vf_details

            # Extract protein description
            desc_match = re.search(r"\)\s+(.*?)\s+\[", header)
            description = desc_match.group(1).strip() if desc_match else header

            entry = {
                "vfg_id": vfg_id,
                "gene_name": gene_name,
                "vf_name": vf_name,
                "vf_class": vf_class,
                "organism": organism,
                "description": description,
                "raw_header": header
            }

            metadata[first_token] = entry
            metadata[vfg_id] = entry
            
            # Also index by accession if present inside (gb|...)
            gb_match = re.search(r"gb\|([a-zA-Z0-9_\.]+)", first_token)
            if gb_match:
                metadata[gb_match.group(1)] = entry

    logger.info("Indexed %d VFDB virulence factor identifiers.", len(metadata))
    return metadata


def parse_card_metadata(aro_index_path: Path, card_fasta: Path) -> Dict[str, Dict[str, str]]:
    """
    Parses CARD aro_index.tsv and protein FASTA headers into structured metadata
    mapped by ARO Accession, Model ID, Protein Accession, and FASTA sequence tokens.
    """
    logger.info("Parsing CARD ARO index from %s...", aro_index_path.name)
    metadata = {}
    if not aro_index_path.exists() or aro_index_path.stat().st_size == 0:
        logger.warning("CARD aro_index.tsv not found or empty.")
        return metadata

    df_aro = pd.read_csv(aro_index_path, sep="\t", dtype=str).fillna("")
    aro_by_id = {}
    
    for _, row in df_aro.iterrows():
        aro_acc = row.get("ARO Accession", "").strip()
        aro_name = row.get("ARO Name", "").strip()
        gene_family = row.get("AMR Gene Family", "").strip()
        drug_class = row.get("Drug Class", "").strip()
        mechanism = row.get("Resistance Mechanism", "").strip()
        prot_acc = row.get("Protein Accession", "").strip()
        model_name = row.get("Model Name", "").strip()

        entry = {
            "aro_accession": aro_acc,
            "aro_name": aro_name or model_name,
            "amr_gene_family": gene_family,
            "drug_class": drug_class,
            "resistance_mechanism": mechanism or "Antibiotic Resistance Determinant",
            "protein_accession": prot_acc,
        }

        if aro_acc:
            aro_by_id[aro_acc] = entry
            metadata[aro_acc] = entry
        if prot_acc:
            metadata[prot_acc] = entry
            if "." in prot_acc:
                metadata[prot_acc.split(".")[0]] = entry

    # Parse FASTA headers to map every accession/token to its ARO entry
    if card_fasta.exists():
        logger.info("Indexing CARD FASTA headers from %s...", card_fasta.name)
        with open(card_fasta, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                tokens = header.split()
                first_token = tokens[0]
                aro_match = re.search(r"ARO:(\d+)", header)
                if aro_match:
                    aro_acc = f"ARO:{aro_match.group(1)}"
                    if aro_acc in aro_by_id:
                        entry = aro_by_id[aro_acc]
                        metadata[first_token] = entry
                        for part in first_token.split("|"):
                            clean_part = part.strip()
                            if clean_part:
                                metadata[clean_part] = entry
                                if "." in clean_part:
                                    metadata[clean_part.split(".")[0]] = entry

    logger.info("Indexed %d CARD ARO resistance determinants.", len(metadata))
    return metadata


def run_mmseqs_search(
    mmseqs_bin: str,
    query_fasta: Path,
    target_fasta: Path,
    output_tsv: Path,
    threads: int,
    max_evalue: float,
    tmp_base: Path,
) -> Path:
    """Executes MMseqs2 easy-search with standard tabular alignment formatting."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = tmp_base / f"tmp_mmseqs_{output_tsv.stem}_{int(time.time()*1000)}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    format_fields = (
        "query,target,fident,alnlen,mismatch,gapopen,"
        "qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov"
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

    logger.info("Running MMseqs2: %s -> %s (threads: %d)", query_fasta.name, target_fasta.name, threads)
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    shutil.rmtree(tmp_dir, ignore_errors=True)

    if res.returncode != 0:
        logger.error("MMseqs2 easy-search failed:\n%s\n%s", res.stdout, res.stderr)
        raise RuntimeError(f"MMseqs2 easy-search failed with code {res.returncode}")

    logger.info("MMseqs2 search completed in %.2f s (output: %s).", elapsed, output_tsv.name)
    return output_tsv


def parse_best_hits(
    search_tsv: Path,
    max_evalue: float,
    min_pident: float,
    min_bitscore: float,
    min_query_cov: float,
) -> Dict[str, Dict[str, Any]]:
    """
    Parses MMseqs2 tabular search output and returns best significant hit per query.
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

    # Significant criteria
    sig_mask = (
        (df["evalue"] <= max_evalue) &
        (df["fident"] >= (min_pident / 100.0)) &
        (df["bits"] >= min_bitscore) &
        (df["qcov"] >= min_query_cov)
    )
    df_sig = df[sig_mask].copy()

    # Sort best hits: highest bitscore, lowest evalue
    df_sig = df_sig.sort_values(by=["bits", "evalue"], ascending=[False, True])

    for _, row in df_sig.iterrows():
        qid = str(row["query"])
        if qid not in best_hits:
            best_hits[qid] = {
                "target": str(row["target"]),
                "evalue": float(row["evalue"]),
                "fident": float(row["fident"]) * 100.0,
                "bits": float(row["bits"]),
                "qcov": float(row["qcov"]) * 100.0,
                "alnlen": int(row["alnlen"]),
            }

    return best_hits


def export_fasta_records(
    df: pd.DataFrame,
    source_fasta: Path,
    output_fasta: Path,
    id_col: str = "Query_ID"
):
    """Exports FASTA sequences for records present in the filtered DataFrame."""
    target_ids = set(df[id_col].dropna().astype(str).tolist())
    if "Representative_Accession" in df.columns:
        target_ids.update(df["Representative_Accession"].dropna().astype(str).tolist())
    if "Cluster_ID" in df.columns:
        target_ids.update(df["Cluster_ID"].dropna().astype(str).tolist())

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    count = 0

    with open(source_fasta, "r", encoding="utf-8") as f_in, open(output_fasta, "w", encoding="utf-8") as f_out:
        write_current = False
        for line in f_in:
            if line.startswith(">"):
                header = line[1:].strip()
                tokens = header.split()
                header_id = tokens[0]
                acc = header_id.split("|")[1] if "|" in header_id else header_id
                
                if (header_id in target_ids) or (acc in target_ids):
                    write_current = True
                    f_out.write(line)
                    count += 1
                else:
                    write_current = False
            elif write_current:
                f_out.write(line)

    logger.info("Exported %d sequences to %s.", count, output_fasta.name)


def main():
    parser = argparse.ArgumentParser(
        description="Step 6: Virulence Factor (VFDB) & Antimicrobial Resistance (CARD) Screening"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to pipeline configuration YAML file",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        logger.error("Configuration file not found: %s", config_path)
        sys.exit(1)

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    job_name = config.get("job_name", "klebsiella_pneumoniae")
    results_dir = Path(config.get("results_dir", "results")) / job_name
    databases_dir = Path(config.get("step4_subtractive_filter", {}).get("databases_dir", "databases"))
    if not databases_dir.is_absolute():
        databases_dir = Path.cwd() / databases_dir

    step6_cfg = config.get("step6_virulence_resistance", {})
    if not step6_cfg.get("enabled", True):
        logger.info("Step 6 (Virulence & Resistance Screening) is disabled in config. Skipping.")
        sys.exit(0)

    # Resolve hardware multithreading
    threads = resolve_threads(step6_cfg.get("threads", "auto"), config, task_name="Step 6")
    logger.info("================================================================================")
    logger.info("Step 6: Virulence Factor (VFDB) & Resistance (CARD) Screening")
    logger.info("Job Name: %s | Allocated Multithreading: %d cores", job_name, threads)
    logger.info("================================================================================")

    # Verify MMseqs2
    mmseqs_bin = ensure_mmseqs_binary()

    # Determine input dataset from Step 5
    input_dataset_name = step6_cfg.get("input_dataset", "prioritized_drug_candidates")
    step5_tables_dir = results_dir / "05_subcellular_localization" / "tables"
    step5_fasta_dir = results_dir / "05_subcellular_localization" / "fasta"

    input_tsv = step5_tables_dir / f"{input_dataset_name}.tsv"
    if not input_tsv.exists():
        input_csv = step5_tables_dir / f"{input_dataset_name}.csv"
        if input_csv.exists():
            input_tsv = input_csv
        else:
            raise FileNotFoundError(f"Step 5 input table not found: {input_tsv}")

    input_fasta = step5_fasta_dir / f"{input_dataset_name}.fasta"
    if not input_fasta.exists():
        fallback_fasta = step5_fasta_dir / "prioritized_drug_candidates.fasta"
        if fallback_fasta.exists():
            input_fasta = fallback_fasta
        else:
            raise FileNotFoundError(f"Step 5 input FASTA not found: {input_fasta}")

    logger.info("Input Candidates Table: %s", input_tsv)
    logger.info("Input Candidates FASTA: %s", input_fasta)

    # Load input table
    df_input = pd.read_csv(input_tsv, sep="\t" if str(input_tsv).endswith(".tsv") else ",", low_memory=False)
    logger.info("Loaded %d candidate targets from Step 5.", len(df_input))

    # Significance parameters
    max_evalue = float(step6_cfg.get("max_evalue", 1.0e-5))
    min_pident = float(step6_cfg.get("min_pident", 30.0))
    min_bitscore = float(step6_cfg.get("min_bitscore", 100.0))
    min_query_cov = float(step6_cfg.get("min_query_cov", 0.35))

    logger.info("Homology Cutoffs: E-value <= %s, Identity >= %.1f%%, Bitscore >= %.1f, Cov >= %.1f%%",
                max_evalue, min_pident, min_bitscore, min_query_cov * 100)

    # 1. Provision Databases
    vfdb_cfg = step6_cfg.get("vfdb", {})
    card_cfg = step6_cfg.get("card", {})

    vfdb_fasta = ensure_vfdb_database(vfdb_cfg, databases_dir)
    card_fasta, aro_index_tsv = ensure_card_database(card_cfg, databases_dir)

    # Parse metadata dictionaries
    vfdb_meta = parse_vfdb_metadata(vfdb_fasta)
    card_meta = parse_card_metadata(aro_index_tsv, card_fasta)

    # Output directories
    step6_out = results_dir / "06_virulence_resistance_screening"
    tables_dir = step6_out / "tables"
    fasta_dir = step6_out / "fasta"
    summary_dir = step6_out / "summary"
    raw_search_dir = step6_out / "raw_searches"

    for d in (tables_dir, fasta_dir, summary_dir, raw_search_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 2. Run MMseqs2 searches
    vfdb_search_tsv = raw_search_dir / "vfdb_search_results.tsv"
    card_search_tsv = raw_search_dir / "card_search_results.tsv"

    run_mmseqs_search(mmseqs_bin, input_fasta, vfdb_fasta, vfdb_search_tsv, threads, max_evalue, step6_out)
    run_mmseqs_search(mmseqs_bin, input_fasta, card_fasta, card_search_tsv, threads, max_evalue, step6_out)

    # 3. Parse Significant Hits
    vfdb_hits = parse_best_hits(vfdb_search_tsv, max_evalue, min_pident, min_bitscore, min_query_cov)
    card_hits = parse_best_hits(card_search_tsv, max_evalue, min_pident, min_bitscore, min_query_cov)

    logger.info("Found %d significant VFDB virulence factor hits.", len(vfdb_hits))
    logger.info("Found %d significant CARD antimicrobial resistance hits.", len(card_hits))

    # 4. Annotate and Stratify DataFrame
    annotated_rows = []
    
    count_dual_action = 0
    count_clean_essential = 0
    count_resistance_breaker = 0
    count_anti_virulence = 0
    count_other = 0

    vf_class_counts = {}
    card_mech_counts = {}
    card_drug_counts = {}

    for _, row in df_input.iterrows():
        r = dict(row)
        
        query_id = str(r.get("Query_ID", ""))
        rep_acc = str(r.get("Representative_Accession", ""))
        clust_id = str(r.get("Cluster_ID", ""))

        hit_vf = vfdb_hits.get(query_id) or vfdb_hits.get(rep_acc) or vfdb_hits.get(clust_id)
        hit_card = card_hits.get(query_id) or card_hits.get(rep_acc) or card_hits.get(clust_id)

        # Annotate VFDB
        if hit_vf:
            r["VFDB_Hit"] = True
            r["VFDB_Hit_ID"] = hit_vf["target"]
            r["VFDB_Evalue"] = hit_vf["evalue"]
            r["VFDB_Identity_Pct"] = round(hit_vf["fident"], 2)
            r["VFDB_Bitscore"] = hit_vf["bits"]
            r["VFDB_Query_Cov"] = round(hit_vf["qcov"], 2)

            target_str = hit_vf["target"]
            target_vfg = target_str.split()[0] if target_str else ""
            meta = vfdb_meta.get(target_vfg, {})
            if not meta and target_vfg.startswith("VFG"):
                meta = vfdb_meta.get(target_vfg[:9], {})

            vf_class = meta.get("vf_class", "Virulence Determinant")
            r["VFDB_VFG_ID"] = meta.get("vfg_id", target_vfg)
            r["VFDB_Gene_Name"] = meta.get("gene_name", "")
            r["VFDB_VF_Name"] = meta.get("vf_name", "")
            r["VFDB_Class"] = vf_class
            r["VFDB_Organism"] = meta.get("organism", "")

            vf_class_counts[vf_class] = vf_class_counts.get(vf_class, 0) + 1
        else:
            r["VFDB_Hit"] = False
            r["VFDB_Hit_ID"] = ""
            r["VFDB_Evalue"] = ""
            r["VFDB_Identity_Pct"] = ""
            r["VFDB_Bitscore"] = ""
            r["VFDB_Query_Cov"] = ""
            r["VFDB_VFG_ID"] = ""
            r["VFDB_Gene_Name"] = ""
            r["VFDB_VF_Name"] = ""
            r["VFDB_Class"] = ""
            r["VFDB_Organism"] = ""

        # Annotate CARD
        if hit_card:
            r["CARD_Hit"] = True
            r["CARD_Hit_ID"] = hit_card["target"]
            r["CARD_Evalue"] = hit_card["evalue"]
            r["CARD_Identity_Pct"] = round(hit_card["fident"], 2)
            r["CARD_Bitscore"] = hit_card["bits"]
            r["CARD_Query_Cov"] = round(hit_card["qcov"], 2)

            card_target = hit_card["target"]
            aro_match = re.search(r"ARO:(\d+)", card_target)
            aro_acc = f"ARO:{aro_match.group(1)}" if aro_match else ""
            
            prot_match = re.search(r"gb\|([a-zA-Z0-9_\.]+)", card_target)
            prot_acc = prot_match.group(1) if prot_match else ""

            meta = (
                card_meta.get(card_target)
                or card_meta.get(card_target.split(".")[0])
                or (card_meta.get(aro_acc) if aro_acc else {})
                or (card_meta.get(prot_acc) if prot_acc else {})
                or {}
            )

            aro_acc = meta.get("aro_accession", aro_acc)
            aro_name = meta.get("aro_name", card_target)
            gene_family = meta.get("amr_gene_family", "")
            drug_class = meta.get("drug_class", "")
            mechanism = meta.get("resistance_mechanism", "Antibiotic Resistance Determinant")

            r["CARD_ARO_Accession"] = aro_acc
            r["CARD_ARO_Name"] = aro_name
            r["CARD_AMR_Gene_Family"] = gene_family
            r["CARD_Drug_Class"] = drug_class
            r["CARD_Resistance_Mechanism"] = mechanism

            card_mech_counts[mechanism] = card_mech_counts.get(mechanism, 0) + 1
            if drug_class:
                for dc in [x.strip() for x in drug_class.split(";") if x.strip()]:
                    card_drug_counts[dc] = card_drug_counts.get(dc, 0) + 1
        else:
            r["CARD_Hit"] = False
            r["CARD_Hit_ID"] = ""
            r["CARD_Evalue"] = ""
            r["CARD_Identity_Pct"] = ""
            r["CARD_Bitscore"] = ""
            r["CARD_Query_Cov"] = ""
            r["CARD_ARO_Accession"] = ""
            r["CARD_ARO_Name"] = ""
            r["CARD_AMR_Gene_Family"] = ""
            r["CARD_Drug_Class"] = ""
            r["CARD_Resistance_Mechanism"] = ""

        # Pharmacological Stratification
        is_essential = (r.get("DEG_Essential") is True) or (str(r.get("DEG_Essential")).lower() == "true")
        has_vfdb = r["VFDB_Hit"] is True
        has_card = r["CARD_Hit"] is True

        if is_essential and has_vfdb:
            r["Therapeutic_Profile"] = "Dual-Action Super-Target (Lethal + Anti-Virulence)"
            count_dual_action += 1
        elif is_essential and not has_vfdb and not has_card:
            r["Therapeutic_Profile"] = "Clean Essential Drug Target (Low Resistance Risk)"
            count_clean_essential += 1
        elif has_card:
            r["Therapeutic_Profile"] = "Resistance-Breaker / Adjuvant Candidate"
            count_resistance_breaker += 1
        elif has_vfdb and not is_essential:
            r["Therapeutic_Profile"] = "Anti-Virulence Candidate (Disarming Pathogen)"
            count_anti_virulence += 1
        else:
            r["Therapeutic_Profile"] = "Secondary / Accessory Target"
            count_other += 1

        annotated_rows.append(r)

    df_master = pd.DataFrame(annotated_rows)

    # 5. Export Master and Subsets
    master_tsv = tables_dir / "virulence_resistance_master.tsv"
    master_csv = tables_dir / "virulence_resistance_master.csv"
    df_master.to_csv(master_tsv, sep="\t", index=False)
    df_master.to_csv(master_csv, sep=",", index=False)
    logger.info("Saved Master Step 6 table: %s (%d targets)", master_tsv.name, len(df_master))

    # Dual-Action Super-Targets
    df_dual = df_master[df_master["Therapeutic_Profile"] == "Dual-Action Super-Target (Lethal + Anti-Virulence)"]
    dual_tsv = tables_dir / "dual_action_super_targets.tsv"
    dual_csv = tables_dir / "dual_action_super_targets.csv"
    df_dual.to_csv(dual_tsv, sep="\t", index=False)
    df_dual.to_csv(dual_csv, sep=",", index=False)
    export_fasta_records(df_dual, input_fasta, fasta_dir / "dual_action_super_targets.fasta")

    # Clean Essential Targets
    df_clean = df_master[df_master["Therapeutic_Profile"] == "Clean Essential Drug Target (Low Resistance Risk)"]
    clean_tsv = tables_dir / "clean_essential_drug_targets.tsv"
    clean_csv = tables_dir / "clean_essential_drug_targets.csv"
    df_clean.to_csv(clean_tsv, sep="\t", index=False)
    df_clean.to_csv(clean_csv, sep=",", index=False)
    export_fasta_records(df_clean, input_fasta, fasta_dir / "clean_essential_drug_targets.fasta")

    # Resistance-Breakers
    df_res = df_master[df_master["Therapeutic_Profile"] == "Resistance-Breaker / Adjuvant Candidate"]
    res_tsv = tables_dir / "resistance_breaker_adjuvant_targets.tsv"
    res_csv = tables_dir / "resistance_breaker_adjuvant_targets.csv"
    df_res.to_csv(res_tsv, sep="\t", index=False)
    df_res.to_csv(res_csv, sep=",", index=False)
    if len(df_res) > 0:
        export_fasta_records(df_res, input_fasta, fasta_dir / "resistance_breaker_adjuvant_targets.fasta")

    # Anti-Virulence Candidates
    df_antivir = df_master[df_master["Therapeutic_Profile"] == "Anti-Virulence Candidate (Disarming Pathogen)"]
    antivir_tsv = tables_dir / "anti_virulence_targets.tsv"
    antivir_csv = tables_dir / "anti_virulence_targets.csv"
    df_antivir.to_csv(antivir_tsv, sep="\t", index=False)
    df_antivir.to_csv(antivir_csv, sep=",", index=False)
    if len(df_antivir) > 0:
        export_fasta_records(df_antivir, input_fasta, fasta_dir / "anti_virulence_targets.fasta")

    # 6. Generate Summary JSON & TXT
    summary = {
        "job_name": job_name,
        "input_dataset": input_dataset_name,
        "total_targets_evaluated": len(df_master),
        "vfdb_hits_count": int((df_master["VFDB_Hit"] == True).sum()),
        "card_hits_count": int((df_master["CARD_Hit"] == True).sum()),
        "therapeutic_profiles": {
            "dual_action_super_targets": count_dual_action,
            "clean_essential_drug_targets": count_clean_essential,
            "resistance_breaker_adjuvants": count_resistance_breaker,
            "anti_virulence_candidates": count_anti_virulence,
            "secondary_targets": count_other,
        },
        "top_virulence_classes": dict(sorted(vf_class_counts.items(), key=lambda x: x[1], reverse=True)[:10]),
        "top_resistance_mechanisms": dict(sorted(card_mech_counts.items(), key=lambda x: x[1], reverse=True)[:10]),
        "top_resistance_drug_classes": dict(sorted(card_drug_counts.items(), key=lambda x: x[1], reverse=True)[:10]),
        "parameters": {
            "max_evalue": max_evalue,
            "min_pident": min_pident,
            "min_bitscore": min_bitscore,
            "min_query_cov": min_query_cov,
            "threads": threads,
        },
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    summary_json = summary_dir / "step6_virulence_resistance_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    summary_txt = summary_dir / "step6_virulence_resistance_summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("================================================================================\n")
        f.write(f"Step 6: Virulence Factor (VFDB) & Resistance (CARD) Screening Summary\n")
        f.write(f"Job: {job_name} | Date: {summary['timestamp']}\n")
        f.write("================================================================================\n\n")
        f.write(f"Total Targets Evaluated: {summary['total_targets_evaluated']}\n")
        f.write(f"VFDB Virulence Hits:     {summary['vfdb_hits_count']}\n")
        f.write(f"CARD Resistance Hits:    {summary['card_hits_count']}\n\n")
        f.write("Therapeutic Profile Breakdown:\n")
        f.write(f"  - Dual-Action Super-Targets (Lethal + Virulence): {count_dual_action}\n")
        f.write(f"  - Clean Essential Drug Targets (Low Resistance):  {count_clean_essential}\n")
        f.write(f"  - Resistance-Breakers / Adjuvants:               {count_resistance_breaker}\n")
        f.write(f"  - Anti-Virulence Targets (Disarming):             {count_anti_virulence}\n")
        f.write(f"  - Secondary / Accessory Targets:                 {count_other}\n\n")
        f.write("Top Virulence Factor Categories:\n")
        for k, v in summary["top_virulence_classes"].items():
            f.write(f"  * {k}: {v}\n")
        f.write("\nTop Resistance Mechanisms:\n")
        for k, v in summary["top_resistance_mechanisms"].items():
            f.write(f"  * {k}: {v}\n")
        f.write("================================================================================\n")

    # Touch completion flag
    done_file = step6_out / ".virulence_resistance_complete"
    done_file.touch()

    logger.info("================================================================================")
    logger.info("Step 6 Completed Successfully!")
    logger.info("Dual-Action Super-Targets: %d", count_dual_action)
    logger.info("Clean Essential Targets:  %d", count_clean_essential)
    logger.info("Resistance-Breakers:      %d", count_resistance_breaker)
    logger.info("Anti-Virulence Targets:   %d", count_anti_virulence)
    logger.info("Summary saved to: %s", summary_json)
    logger.info("================================================================================")


if __name__ == "__main__":
    main()
