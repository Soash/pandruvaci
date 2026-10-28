#!/usr/bin/env python3
"""
Step 7: Target Druggability & Novelty Assessment (TTD 2026)
==========================================================
Cross-references prioritized antibacterial candidates against the Therapeutic
Target Database (TTD 2026 update) to evaluate clinical druggability, discover
drug-repurposing candidates, and identify completely novel, first-in-class
patentable antibacterial targets.

Academic Citations:
- TTD 2026: Zhang YT, Zhou Y, Xu HW, et al. Nucleic Acids Res (2026) 54:D1692–D1701.
- MMseqs2: Steinegger M & Söding J. Nature Biotechnology (2017) 35:1026–1028.

Usage:
    python scripts/07_target_druggability_ttd.py [--config config.yaml]
"""

import os
import re
import sys
import time
import json
import shutil
import logging
import argparse
import subprocess
import urllib.request
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, Set, List
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
logger = logging.getLogger("07_target_druggability_ttd")


def ensure_mmseqs_binary() -> str:
    """Finds MMseqs2 binary in PATH or local bin directory."""
    which_mmseqs = shutil.which("mmseqs")
    if which_mmseqs:
        return which_mmseqs

    project_root = Path(__file__).resolve().parent.parent
    local_bin = project_root / "bin" / "mmseqs"
    if local_bin.exists() and os.access(local_bin, os.X_OK):
        return str(local_bin)

    raise RuntimeError("MMseqs2 binary not found.")


def download_file_with_retry(
    url: str,
    output_path: Path,
    description: str,
    retries: int = 3,
    timeout: int = 180,
) -> Path:
    """Downloads a file with retries and temp-file atomic rename."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.stat().st_size > 0:
        logger.info("Found existing %s at: %s (%s bytes)", description, output_path.name, f"{output_path.stat().st_size:,}")
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


def provision_ttd_databases(ttd_cfg: Dict[str, Any], databases_dir: Path) -> Dict[str, Path]:
    """
    Downloads and prepares TTD 2026 database files.
    Returns dictionary of file paths.
    """
    ttd_dir = Path(ttd_cfg.get("database_dir", databases_dir / "ttd"))
    if not ttd_dir.is_absolute():
        ttd_dir = Path.cwd() / ttd_dir
    ttd_dir.mkdir(parents=True, exist_ok=True)

    urls = ttd_cfg.get("ttd_urls", {})
    files = {
        "target_download": (urls.get("target_download", "https://ttd.idrblab.cn/files/download/P1-01-TTD_target_download.txt"), "P1-01-TTD_target_download.txt"),
        "uniprot_all": (urls.get("uniprot_all", "https://ttd.idrblab.cn/files/download/P2-01-TTD_uniprot_all.txt"), "P2-01-TTD_uniprot_all.txt"),
        "sequence_all": (urls.get("sequence_all", "https://ttd.idrblab.cn/files/download/P2-06-TTD_sequence_all.txt"), "P2-06-TTD_sequence_all.txt"),
        "drug_download": (urls.get("drug_download", "https://ttd.idrblab.cn/files/download/P1-02-TTD_drug_download.txt"), "P1-02-TTD_drug_download.txt"),
    }

    paths = {}
    for key, (url, fname) in files.items():
        dest = ttd_dir / fname
        download_file_with_retry(url, dest, f"TTD {fname}")
        paths[key] = dest

    return paths


def parse_ttd_target_metadata(p1_01_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parses P1-01-TTD_target_download.txt into structured metadata dictionary.
    File format (tab-separated):
    TARGETID \t Field \t Value [\t Extra ...]
    """
    logger.info("Parsing TTD target metadata from %s...", p1_01_path.name)
    targets: Dict[str, Dict[str, Any]] = {}

    status_priority = {
        "approved": 10,
        "phase 3": 9,
        "phase 2/phase 3": 8,
        "phase 2": 7,
        "phase 1/phase 2": 6,
        "phase 1": 5,
        "investigational": 4,
        "preclinical": 3,
        "patented": 2,
        "research": 1,
        "terminated": 0,
        "withdrawn": 0,
    }

    with open(p1_01_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("-") or line.startswith("Download") or line.startswith("Title"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            tid = parts[0].strip()
            field = parts[1].strip().upper()
            val = parts[2].strip()

            if tid not in targets:
                targets[tid] = {
                    "target_id": tid,
                    "target_name": "",
                    "target_type": "Research target",
                    "uniprot_ids": set(),
                    "bioclass": "",
                    "ec_number": "",
                    "drugs": [],
                    "drug_names": set(),
                    "highest_status": "Novel / Unexplored",
                    "highest_score": -1,
                }

            t = targets[tid]

            if field == "TARGNAME":
                t["target_name"] = val
            elif field == "TARGTYPE":
                t["target_type"] = val
            elif field == "UNIPROID":
                t["uniprot_ids"].add(val)
                # If UniProt entry name like FOLL_ECOLI or accession
                if "_" in val:
                    t["uniprot_ids"].add(val.split("_")[0])
            elif field == "BIOCLASS":
                t["bioclass"] = val
            elif field == "ECNUMBER":
                t["ec_number"] = val
            elif field == "DRUGINFO":
                drug_id = val
                drug_name = parts[3].strip() if len(parts) > 3 else ""
                status = parts[4].strip() if len(parts) > 4 else "Investigational"
                
                t["drugs"].append({
                    "drug_id": drug_id,
                    "drug_name": drug_name,
                    "status": status,
                })
                if drug_name:
                    t["drug_names"].add(drug_name)

                score = status_priority.get(status.lower(), 1)
                if score > t["highest_score"]:
                    t["highest_score"] = score
                    t["highest_status"] = status

    # Clean highest status
    for tid, t in targets.items():
        if t["highest_status"] == "Novel / Unexplored":
            if "successful" in t["target_type"].lower():
                t["highest_status"] = "Approved"
            elif "clinical" in t["target_type"].lower():
                t["highest_status"] = "Clinical Trial"
            elif "patented" in t["target_type"].lower():
                t["highest_status"] = "Patented"
            else:
                t["highest_status"] = "Research"

    logger.info("Parsed %d structured TTD targets from P1-01.", len(targets))
    return targets


def parse_ttd_uniprot_crossmatch(p2_01_path: Path) -> Dict[str, str]:
    """
    Parses P2-01-TTD_uniprot_all.txt to build mapping from UniProt ID/Accession -> TTD TargetID.
    """
    logger.info("Parsing UniProt to TTD crossmatching from %s...", p2_01_path.name)
    uniprot_map: Dict[str, str] = {}

    with open(p2_01_path, "r", encoding="utf-8", errors="replace") as f:
        curr_tid = ""
        for line in f:
            line = line.strip()
            if not line or line.startswith("-") or line.startswith("Abbreviation"):
                continue

            parts = line.split("\t")
            if len(parts) >= 2:
                field = parts[0].strip().upper()
                val = parts[1].strip()
                if field == "TARGETID":
                    curr_tid = val
                elif field == "UNIPROID" and curr_tid:
                    uniprot_map[val] = curr_tid
                    if "_" in val:
                        uniprot_map[val.split("_")[0]] = curr_tid

    logger.info("Indexed %d UniProt-to-TTD accession cross-matches.", len(uniprot_map))
    return uniprot_map


def clean_ttd_fasta_database(raw_seq_path: Path, clean_fasta_path: Path) -> Path:
    """
    Cleans P2-06-TTD_sequence_all.txt by removing text headers and ensuring valid FASTA.
    """
    if clean_fasta_path.exists() and clean_fasta_path.stat().st_size > 0:
        return clean_fasta_path

    logger.info("Preparing clean FASTA database from %s -> %s...", raw_seq_path.name, clean_fasta_path.name)
    clean_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0

    with open(raw_seq_path, "r", encoding="utf-8", errors="replace") as f_in, open(clean_fasta_path, "w", encoding="utf-8") as f_out:
        in_fasta = False
        for line in f_in:
            if line.startswith(">"):
                in_fasta = True
                f_out.write(line)
                count += 1
            elif in_fasta:
                f_out.write(line)

    logger.info("Cleaned %d TTD FASTA sequences for MMseqs2 search.", count)
    return clean_fasta_path


def run_mmseqs_search(
    mmseqs_bin: str,
    query_fasta: Path,
    target_fasta: Path,
    output_tsv: Path,
    threads: int,
    max_evalue: float,
    tmp_base: Path,
) -> Path:
    """Executes MMseqs2 easy-search against TTD sequence database."""
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

    logger.info("Running MMseqs2 search: %s -> %s (threads: %d)", query_fasta.name, target_fasta.name, threads)
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    shutil.rmtree(tmp_dir, ignore_errors=True)

    if res.returncode != 0:
        logger.error("MMseqs2 search failed:\n%s\n%s", res.stdout, res.stderr)
        raise RuntimeError(f"MMseqs2 search failed with code {res.returncode}")

    logger.info("MMseqs2 search completed in %.2f s (output: %s).", elapsed, output_tsv.name)
    return output_tsv


def parse_best_hits(
    search_tsv: Path,
    max_evalue: float,
    min_pident: float,
    min_bitscore: float,
    min_query_cov: float,
) -> Dict[str, Dict[str, Any]]:
    """Parses MMseqs2 results and returns top significant hit per query."""
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

    numeric_cols = ["fident", "alnlen", "evalue", "bits", "qlen", "tlen", "qcov"]
    for c in numeric_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    sig_mask = (
        (df["evalue"] <= max_evalue) &
        (df["fident"] >= (min_pident / 100.0)) &
        (df["bits"] >= min_bitscore) &
        (df["qcov"] >= min_query_cov)
    )
    df_sig = df[sig_mask].copy()
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
    """Exports FASTA sequences for records in DataFrame."""
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
        description="Step 7: Target Druggability & Novelty Assessment (TTD 2026)"
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

    step7_cfg = config.get("step7_target_druggability_ttd", {})
    if not step7_cfg.get("enabled", True):
        logger.info("Step 7 (TTD Target Druggability) is disabled in config. Skipping.")
        sys.exit(0)

    threads = resolve_threads(step7_cfg.get("threads", "auto"), config, task_name="Step 7")
    logger.info("================================================================================")
    logger.info("Step 7: Target Druggability & Novelty Assessment (TTD 2026)")
    logger.info("Job Name: %s | Allocated Multithreading: %d cores", job_name, threads)
    logger.info("================================================================================")

    mmseqs_bin = ensure_mmseqs_binary()

    # Determine input dataset from Step 6
    step6_tables_dir = results_dir / "06_virulence_resistance_screening" / "tables"
    step6_master_tsv = step6_tables_dir / "virulence_resistance_master.tsv"
    if not step6_master_tsv.exists():
        step6_master_tsv = step6_tables_dir / "virulence_resistance_master.csv"
        if not step6_master_tsv.exists():
            raise FileNotFoundError(f"Step 6 master table not found: {step6_master_tsv}")

    # Query FASTA sequences
    step5_fasta = results_dir / "05_subcellular_localization" / "fasta" / "prioritized_drug_candidates.fasta"
    if not step5_fasta.exists():
        raise FileNotFoundError(f"Query FASTA not found at: {step5_fasta}")

    df_input = pd.read_csv(step6_master_tsv, sep="\t" if str(step6_master_tsv).endswith(".tsv") else ",", low_memory=False)
    logger.info("Loaded %d candidate targets from Step 6.", len(df_input))

    # Significance cutoffs
    max_evalue = float(step7_cfg.get("max_evalue", 1.0e-5))
    min_pident = float(step7_cfg.get("min_pident", 30.0))
    min_bitscore = float(step7_cfg.get("min_bitscore", 100.0))
    min_query_cov = float(step7_cfg.get("min_query_cov", 0.35))

    logger.info("Homology Cutoffs: E-value <= %s, Identity >= %.1f%%, Bitscore >= %.1f, Cov >= %.1f%%",
                max_evalue, min_pident, min_bitscore, min_query_cov * 100)

    # 1. Provision Databases
    ttd_paths = provision_ttd_databases(step7_cfg, databases_dir)
    ttd_dir = ttd_paths["target_download"].parent

    clean_ttd_fasta = clean_ttd_fasta_database(ttd_paths["sequence_all"], ttd_dir / "TTD_sequence_clean.fasta")
    ttd_targets_meta = parse_ttd_target_metadata(ttd_paths["target_download"])
    ttd_uniprot_map = parse_ttd_uniprot_crossmatch(ttd_paths["uniprot_all"])

    # Output directories
    step7_out = results_dir / "07_target_druggability_ttd"
    tables_dir = step7_out / "tables"
    fasta_dir = step7_out / "fasta"
    summary_dir = step7_out / "summary"
    raw_search_dir = step7_out / "raw_searches"

    for d in (tables_dir, fasta_dir, summary_dir, raw_search_dir):
        d.mkdir(parents=True, exist_ok=True)

    # 2. Run MMseqs2 Search against TTD Clean Sequence Database
    ttd_search_tsv = raw_search_dir / "ttd_search_results.tsv"
    run_mmseqs_search(mmseqs_bin, step5_fasta, clean_ttd_fasta, ttd_search_tsv, threads, max_evalue, step7_out)

    # 3. Parse Significant MMseqs2 Hits
    ttd_hits = parse_best_hits(ttd_search_tsv, max_evalue, min_pident, min_bitscore, min_query_cov)
    logger.info("Found %d significant MMseqs2 hits against TTD targets.", len(ttd_hits))

    # 4. Annotate & Stratify Targets
    annotated_rows = []

    count_approved = 0
    count_clinical = 0
    count_patented = 0
    count_research = 0
    count_novel = 0
    count_novel_dual = 0
    count_novel_clean = 0
    count_novel_adjuvant = 0

    clinical_status_counts = {}

    for _, row in df_input.iterrows():
        r = dict(row)

        qid = str(r.get("Query_ID", ""))
        rep_acc = str(r.get("Representative_Accession", ""))
        clust_id = str(r.get("Cluster_ID", ""))
        gene_name = str(r.get("Gene_Name", ""))
        therap_prof = str(r.get("Therapeutic_Profile", ""))

        hit_seq = ttd_hits.get(qid) or ttd_hits.get(rep_acc) or ttd_hits.get(clust_id)

        # Check exact UniProt crossmatch
        exact_tid = ttd_uniprot_map.get(rep_acc) or ttd_uniprot_map.get(qid)
        if not exact_tid and gene_name:
            exact_tid = ttd_uniprot_map.get(gene_name)

        matched_tid = ""
        matched_meta = {}

        if hit_seq:
            raw_target = hit_seq["target"]
            # Target ID is first token (e.g. T47101)
            tid_cand = raw_target.split()[0]
            matched_tid = tid_cand
            matched_meta = ttd_targets_meta.get(matched_tid, {})
        elif exact_tid:
            matched_tid = exact_tid
            matched_meta = ttd_targets_meta.get(matched_tid, {})

        if matched_tid and matched_meta:
            r["TTD_Hit"] = True
            r["TTD_Target_ID"] = matched_tid
            r["TTD_Target_Name"] = matched_meta.get("target_name", "")
            r["TTD_Target_Type"] = matched_meta.get("target_type", "")
            r["TTD_Clinical_Status"] = matched_meta.get("highest_status", "Investigational")
            r["TTD_Bioclass"] = matched_meta.get("bioclass", "")
            r["TTD_EC_Number"] = matched_meta.get("ec_number", "")
            
            # Format top drugs
            drugs_list = matched_meta.get("drugs", [])
            drug_strs = [f"{d['drug_name']} ({d['status']})" for d in drugs_list[:5] if d.get("drug_name")]
            r["TTD_Associated_Drugs"] = "; ".join(drug_strs) if drug_strs else "None reported"

            if hit_seq:
                r["TTD_Evalue"] = hit_seq["evalue"]
                r["TTD_Identity_Pct"] = round(hit_seq["fident"], 2)
                r["TTD_Bitscore"] = hit_seq["bits"]
                r["TTD_Query_Cov"] = round(hit_seq["qcov"], 2)
                r["TTD_Match_Mode"] = "Sequence Homology (MMseqs2)"
            else:
                r["TTD_Evalue"] = 0.0
                r["TTD_Identity_Pct"] = 100.0
                r["TTD_Bitscore"] = 999.0
                r["TTD_Query_Cov"] = 100.0
                r["TTD_Match_Mode"] = "Exact UniProt Identifier"

            status_clean = r["TTD_Clinical_Status"]
            clinical_status_counts[status_clean] = clinical_status_counts.get(status_clean, 0) + 1

            if "approved" in status_clean.lower():
                r["Clinical_Druggability_Category"] = "Approved Drug Target (Repurposing Opportunity)"
                count_approved += 1
            elif "phase" in status_clean.lower() or "clinical" in status_clean.lower():
                r["Clinical_Druggability_Category"] = "Clinical Trial Target (Phase I-III Validated)"
                count_clinical += 1
            elif "patented" in status_clean.lower():
                r["Clinical_Druggability_Category"] = "Patented / Preclinical Target"
                count_patented += 1
            else:
                r["Clinical_Druggability_Category"] = "Research Target (Literature Precedent)"
                count_research += 1
        else:
            # Completely Novel Target
            r["TTD_Hit"] = False
            r["TTD_Target_ID"] = ""
            r["TTD_Target_Name"] = ""
            r["TTD_Target_Type"] = "Novel / Unexplored"
            r["TTD_Clinical_Status"] = "Novel / Untapped"
            r["TTD_Bioclass"] = ""
            r["TTD_EC_Number"] = ""
            r["TTD_Associated_Drugs"] = ""
            r["TTD_Evalue"] = ""
            r["TTD_Identity_Pct"] = ""
            r["TTD_Bitscore"] = ""
            r["TTD_Query_Cov"] = ""
            r["TTD_Match_Mode"] = "None (Novel)"

            count_novel += 1
            if "dual-action" in therap_prof.lower():
                r["Clinical_Druggability_Category"] = "Novel Dual-Action Super-Target (Pristine First-in-Class)"
                count_novel_dual += 1
            elif "clean" in therap_prof.lower():
                r["Clinical_Druggability_Category"] = "Novel Clean Essential Target (First-in-Class Metabolic Enzyme)"
                count_novel_clean += 1
            elif "resistance-breaker" in therap_prof.lower():
                r["Clinical_Druggability_Category"] = "Novel Resistance-Breaker / Adjuvant Candidate"
                count_novel_adjuvant += 1
            else:
                r["Clinical_Druggability_Category"] = "Novel Antibacterial Target"

        annotated_rows.append(r)

    df_master = pd.DataFrame(annotated_rows)

    # 5. Export Master Table and Distinct Clinical Subsets
    master_tsv = tables_dir / "ttd_druggability_master.tsv"
    master_csv = tables_dir / "ttd_druggability_master.csv"
    df_master.to_csv(master_tsv, sep="\t", index=False)
    df_master.to_csv(master_csv, sep=",", index=False)
    logger.info("Saved Master Step 7 table: %s (%d targets)", master_tsv.name, len(df_master))

    # Approved Targets (Repurposing)
    df_app = df_master[df_master["Clinical_Druggability_Category"] == "Approved Drug Target (Repurposing Opportunity)"]
    app_tsv = tables_dir / "approved_drug_targets.tsv"
    app_csv = tables_dir / "approved_drug_targets.csv"
    df_app.to_csv(app_tsv, sep="\t", index=False)
    df_app.to_csv(app_csv, sep=",", index=False)
    if len(df_app) > 0:
        export_fasta_records(df_app, step5_fasta, fasta_dir / "approved_drug_targets.fasta")

    # Clinical Trial Targets
    df_clin = df_master[df_master["Clinical_Druggability_Category"] == "Clinical Trial Target (Phase I-III Validated)"]
    clin_tsv = tables_dir / "clinical_trial_targets.tsv"
    clin_csv = tables_dir / "clinical_trial_targets.csv"
    df_clin.to_csv(clin_tsv, sep="\t", index=False)
    df_clin.to_csv(clin_csv, sep=",", index=False)
    if len(df_clin) > 0:
        export_fasta_records(df_clin, step5_fasta, fasta_dir / "clinical_trial_targets.fasta")

    # Novel First-in-Class Targets
    df_novel = df_master[df_master["TTD_Hit"] == False]
    novel_tsv = tables_dir / "novel_first_in_class_targets.tsv"
    novel_csv = tables_dir / "novel_first_in_class_targets.csv"
    df_novel.to_csv(novel_tsv, sep="\t", index=False)
    df_novel.to_csv(novel_csv, sep=",", index=False)
    export_fasta_records(df_novel, step5_fasta, fasta_dir / "novel_first_in_class_targets.fasta")

    # Novel Dual-Action Targets
    df_novel_dual = df_master[df_master["Clinical_Druggability_Category"] == "Novel Dual-Action Super-Target (Pristine First-in-Class)"]
    novel_dual_tsv = tables_dir / "novel_dual_action_targets.tsv"
    novel_dual_csv = tables_dir / "novel_dual_action_targets.csv"
    df_novel_dual.to_csv(novel_dual_tsv, sep="\t", index=False)
    df_novel_dual.to_csv(novel_dual_csv, sep=",", index=False)
    if len(df_novel_dual) > 0:
        export_fasta_records(df_novel_dual, step5_fasta, fasta_dir / "novel_dual_action_targets.fasta")

    # Novel Clean Essential Targets
    df_novel_clean = df_master[df_master["Clinical_Druggability_Category"] == "Novel Clean Essential Target (First-in-Class Metabolic Enzyme)"]
    novel_clean_tsv = tables_dir / "novel_clean_essential_targets.tsv"
    novel_clean_csv = tables_dir / "novel_clean_essential_targets.csv"
    df_novel_clean.to_csv(novel_clean_tsv, sep="\t", index=False)
    df_novel_clean.to_csv(novel_clean_csv, sep=",", index=False)
    if len(df_novel_clean) > 0:
        export_fasta_records(df_novel_clean, step5_fasta, fasta_dir / "novel_clean_essential_targets.fasta")

    # 6. Generate Summary JSON & TXT
    summary = {
        "job_name": job_name,
        "total_targets_evaluated": len(df_master),
        "ttd_hits_count": int((df_master["TTD_Hit"] == True).sum()),
        "ttd_novel_targets_count": count_novel,
        "clinical_druggability_breakdown": {
            "approved_drug_targets": count_approved,
            "clinical_trial_targets": count_clinical,
            "patented_preclinical_targets": count_patented,
            "research_literature_targets": count_research,
            "novel_first_in_class_targets": count_novel,
        },
        "novel_targets_stratification": {
            "novel_dual_action_super_targets": count_novel_dual,
            "novel_clean_essential_targets": count_novel_clean,
            "novel_resistance_breaker_targets": count_novel_adjuvant,
        },
        "top_clinical_statuses": dict(sorted(clinical_status_counts.items(), key=lambda x: x[1], reverse=True)),
        "parameters": {
            "max_evalue": max_evalue,
            "min_pident": min_pident,
            "min_bitscore": min_bitscore,
            "min_query_cov": min_query_cov,
            "threads": threads,
        },
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    summary_json = summary_dir / "step7_ttd_druggability_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    summary_txt = summary_dir / "step7_ttd_druggability_summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("================================================================================\n")
        f.write(f"Step 7: Target Druggability & Novelty Assessment (TTD 2026) Summary\n")
        f.write(f"Job: {job_name} | Date: {summary['timestamp']}\n")
        f.write("================================================================================\n\n")
        f.write(f"Total Targets Evaluated:          {summary['total_targets_evaluated']}\n")
        f.write(f"Targets with Known TTD Matches:   {summary['ttd_hits_count']}\n")
        f.write(f"Completely Novel / Untapped:      {summary['ttd_novel_targets_count']}\n\n")
        f.write("Clinical Development Status Breakdown:\n")
        f.write(f"  - Approved / Successful Targets (Repurposing): {count_approved}\n")
        f.write(f"  - Clinical Trial Targets (Phase I-III):       {count_clinical}\n")
        f.write(f"  - Patented / Preclinical Targets:             {count_patented}\n")
        f.write(f"  - Research / Literature Precedent:            {count_research}\n")
        f.write(f"  - Completely Novel / First-in-Class:          {count_novel}\n\n")
        f.write("Novel Targets Stratification (First-in-Class Potential):\n")
        f.write(f"  * Novel Dual-Action Super-Targets:            {count_novel_dual}\n")
        f.write(f"  * Novel Clean Essential Metabolic Targets:    {count_novel_clean}\n")
        f.write(f"  * Novel Resistance-Breakers / Adjuvants:      {count_novel_adjuvant}\n")
        f.write("================================================================================\n")

    done_file = step7_out / ".ttd_druggability_complete"
    done_file.touch()

    logger.info("================================================================================")
    logger.info("Step 7 Completed Successfully!")
    logger.info("Approved Drug Targets (Repurposing): %d", count_approved)
    logger.info("Clinical Trial Targets (Phase I-III): %d", count_clinical)
    logger.info("Completely Novel Targets:            %d", count_novel)
    logger.info("  * Novel Dual-Action Super-Targets: %d", count_novel_dual)
    logger.info("  * Novel Clean Essential Targets:   %d", count_novel_clean)
    logger.info("Summary saved to: %s", summary_json)
    logger.info("================================================================================")


if __name__ == "__main__":
    main()
