#!/usr/bin/env python3
"""
Step 1: Data Acquisition & Local Filtering Script
=================================================
Fetches proteomes metadata from UniProt REST API using a broad query:
    ((taxonomy_id:<id>) NOT (proteome_type:EXCLUDED))

Then applies flexible local filtering on the resulting metadata table:
  - Preserves Reference proteomes (even when CPD is Unknown).
  - Filters Non-Reference proteomes by CPD status (e.g. Standard / cpd:1).
  - Streams FASTA proteome sequences for the selected proteomes.

Usage:
    python scripts/01_data_acquisition.py [--config config.yaml]
"""

import os
import sys
import time
import json
import logging
import argparse
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
import requests
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
logger = logging.getLogger("01_data_acquisition")

# CPD status code mappings
CPD_CODE_TO_STATUS = {
    1: "Standard",
    2: "Close to standard (high value)",
    3: "Close to standard (low value)",
    4: "Outlier (high value)",
    5: "Outlier (low value)",
    6: "Unknown",
}


def load_config(config_path: str) -> dict:
    """Load and validate the YAML configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at {config_path}")
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


def construct_uniprot_query(params: dict) -> str:
    """
    Constructs a clean, broad query string for UniProt.
    Example: ((taxonomy_id:573) AND (proteome_type:REFERENCE))
    """
    clauses = []

    # Taxonomy ID
    tax_id = params.get("taxonomy_id")
    if tax_id:
        clauses.append(f"(taxonomy_id:{tax_id})")

    # Reference Proteome Only filter at query level
    if params.get("reference_only", False):
        clauses.append("(proteome_type:REFERENCE)")
    elif params.get("exclude_excluded", True):
        clauses.append("NOT (proteome_type:EXCLUDED)")

    query_str = " AND ".join(clauses)
    return query_str


def get_uniprot_metadata_url(config: dict) -> str:
    """Returns the full UniProt REST API streaming URL for metadata in JSON format."""
    mode = config.get("query_mode", "builder")

    if mode == "custom_url":
        custom_url = config.get("custom_url", "").strip()
        if custom_url:
            logger.info("Using custom UniProt URL from config.")
            return custom_url

    # Builder mode: retrieve proteomes for the taxon via JSON stream
    params = config.get("query_params", {}).copy()
    if config.get("taxonomy_id"):
        params["taxonomy_id"] = config["taxonomy_id"]
    filter_opts = config.get("filter_options", {})
    # If reference_only is enabled in filter_options or top level, propagate to query params
    if filter_opts.get("reference_only", False) or config.get("reference_only", False):
        params["reference_only"] = True

    query = construct_uniprot_query(params)
    logger.info("Constructed UniProt Query: %s", query)

    base_url = "https://rest.uniprot.org/proteomes/stream"
    url = f"{base_url}?format=json&query={urllib.parse.quote(query)}"
    return url


def fetch_and_parse_proteomes(url: str, timeout: int = 90) -> pd.DataFrame:
    """
    Fetches proteomes from UniProt stream API and parses into a structured DataFrame.
    """
    logger.info("Fetching proteome records from UniProt...")
    logger.debug("Request URL: %s", url)

    response = requests.get(url, timeout=timeout)
    if response.status_code != 200:
        raise RuntimeError(
            f"Failed to fetch metadata from UniProt (Status Code: {response.status_code})\n"
            f"Response: {response.text[:500]}"
        )

    # If response is TSV (e.g. from a custom TSV URL)
    content_type = response.headers.get("Content-Type", "")
    if "tsv" in content_type or "tab-separated" in content_type or response.text.startswith("Proteome Id\t") or response.text.startswith("upid\t"):
        import io
        logger.info("Parsing TSV response from UniProt...")
        return pd.read_csv(io.StringIO(response.text), sep="\t")

    # Parse JSON stream
    logger.info("Parsing JSON stream response from UniProt...")
    try:
        data = response.json()
    except Exception as e:
        raise RuntimeError(f"Error parsing UniProt response as JSON: {e}")

    results = data.get("results", [])
    logger.info("Received %d total proteome entries from UniProt.", len(results))

    rows = []
    for entry in results:
        upid = entry.get("id")
        tax_info = entry.get("taxonomy", {})
        organism = tax_info.get("scientificName", "")
        taxon_id = tax_info.get("taxonId", "")
        strain = entry.get("strain", "")
        proteome_type = entry.get("proteomeType", "")

        # Completeness & CPD report
        completeness = entry.get("proteomeCompletenessReport", {})
        cpd_report = completeness.get("cpdReport", {})
        cpd_status = cpd_report.get("status", "Unknown")

        busco_report = completeness.get("buscoReport", {})
        busco_score = busco_report.get("score", None)

        # Protein count and UniProtKB sequence statistics
        stats = entry.get("proteomeStatistics", {})
        reviewed = stats.get("reviewedProteinCount", 0) or 0
        unreviewed = stats.get("unreviewedProteinCount", 0) or 0
        uniprot_seq_count = reviewed + unreviewed
        genome_protein_count = entry.get("proteinCount", 0) or 0

        # Effective protein count: use UniProtKB count if statistics provided, else genome count
        effective_protein_count = uniprot_seq_count if stats else genome_protein_count
        genome_assembly = entry.get("genomeAssembly", {}).get("assemblyId", "")

        rows.append({
            "Proteome Id": upid,
            "Organism": organism,
            "Taxon Id": taxon_id,
            "Strain": strain,
            "Proteome Type": proteome_type,
            "CPD": cpd_status,
            "Protein Count": effective_protein_count,
            "UniProt Sequences": uniprot_seq_count,
            "Genome Assembly ID": genome_assembly,
            "BUSCO Score": busco_score,
        })

    df = pd.DataFrame(rows)
    return df


def apply_local_filtering(df: pd.DataFrame, filter_opts: dict) -> pd.DataFrame:
    """
    Applies user-defined selection rules:
    - If reference_only is True, selects ONLY Reference proteomes (excluding all non-reference).
    - Otherwise:
      - Retains Reference proteomes (even if CPD is Unknown) if include_reference is True.
      - Filters Non-Reference proteomes by specified CPD levels (if include_non_reference is True).
    - Requires verified protein sequence entries in UniProtKB (skip 0-sequence unindexed assemblies).
    - Applies optional min/max protein count filters.
    """
    if df.empty:
        return df

    ref_only = filter_opts.get("reference_only", False)
    include_ref = filter_opts.get("include_reference", True)
    include_non_ref = filter_opts.get("include_non_reference", not ref_only)
    if ref_only:
        include_non_ref = False
        logger.info("Local filter rule: REFERENCE ONLY mode active (All non-reference proteomes excluded).")

    non_ref_cpd_levels = filter_opts.get("non_reference_cpd_levels", [1])
    min_count = filter_opts.get("min_protein_count")
    max_count = filter_opts.get("max_protein_count")

    # Map requested numeric CPD levels to text statuses
    target_non_ref_statuses = set()
    no_cpd_filter = (0 in non_ref_cpd_levels or "0" in non_ref_cpd_levels)

    if not no_cpd_filter:
        for lvl in non_ref_cpd_levels:
            if lvl in CPD_CODE_TO_STATUS:
                target_non_ref_statuses.add(CPD_CODE_TO_STATUS[lvl])
            elif str(lvl).lower() == "unknown":
                target_non_ref_statuses.add("Unknown")

    selected_mask = pd.Series(False, index=df.index)

    for idx, row in df.iterrows():
        p_type = str(row.get("Proteome Type", "")).lower()
        cpd_val = str(row.get("CPD", "Unknown"))
        is_ref = ("reference" in p_type) and ("non" not in p_type)

        # 1. Reference Proteome check
        if is_ref:
            if include_ref:
                selected_mask.loc[idx] = True
            continue

        # 2. Non-Reference Proteome check
        if not is_ref:
            if include_non_ref and not ref_only:
                if no_cpd_filter or (cpd_val in target_non_ref_statuses):
                    selected_mask.loc[idx] = True

    # 3. Require verified protein sequence entries in UniProtKB (skip 0-sequence unindexed assemblies)
    if "UniProt Sequences" in df.columns:
        has_uniprot_seqs = df["UniProt Sequences"] > 0
    elif "Protein Count" in df.columns:
        has_uniprot_seqs = df["Protein Count"] > 0
    else:
        has_uniprot_seqs = pd.Series(True, index=df.index)

    selected_mask = selected_mask & has_uniprot_seqs

    # 4. Protein count boundaries (if configured)
    if min_count is not None and "Protein Count" in df.columns:
        selected_mask = selected_mask & (df["Protein Count"] >= min_count)
    if max_count is not None and "Protein Count" in df.columns:
        selected_mask = selected_mask & (df["Protein Count"] <= max_count)

    df["Selected"] = selected_mask
    return df


def count_fasta_sequences(fasta_path: Path) -> int:
    """Counts the number of protein records (headers starting with '>') in a FASTA file."""
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return 0
    count = 0
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def download_single_fasta(
    upid: str,
    output_dir: Path,
    expected_count: int = None,
    skip_existing: bool = True,
    timeout: int = 90,
    retry_attempts: int = 3,
) -> dict:
    """Downloads a single proteome FASTA file directly from UniProtKB."""
    target_file = output_dir / f"{upid}.fasta"
    result = {
        "upid": upid,
        "file": str(target_file),
        "source": "uniprot",
        "status": "pending",
        "seq_count": 0,
        "expected_count": expected_count,
        "error": None,
    }

    if skip_existing and target_file.exists() and target_file.stat().st_size > 0:
        seq_count = count_fasta_sequences(target_file)
        if seq_count > 0:
            result["status"] = "skipped"
            result["seq_count"] = seq_count
            return result
        else:
            try:
                target_file.unlink()
            except OSError:
                pass

    uniprot_url = f"https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=proteome:{upid}"

    for attempt in range(1, retry_attempts + 1):
        try:
            resp = requests.get(uniprot_url, timeout=timeout)
            if resp.status_code == 200:
                count = resp.text.count(">")
                if count > 0:
                    with open(target_file, "w", encoding="utf-8") as f:
                        f.write(resp.text)
                    result["status"] = "success"
                    result["seq_count"] = count
                    return result
                else:
                    # UniProt returned empty sequence content for this proteome
                    result["error"] = "Empty sequence response from UniProtKB"
                    result["status"] = "failed"
                    return result
            else:
                result["error"] = f"HTTP {resp.status_code}"
        except Exception as e:
            result["error"] = str(e)

        if attempt < retry_attempts:
            time.sleep(2 * attempt)

    result["status"] = "failed"
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Step 1: Download proteomes metadata from UniProt, filter locally, and download FASTA files."
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML configuration file (default: config.yaml)",
    )
    args = parser.parse_args()

    # 1. Load Configuration
    logger.info("Loading configuration from %s", args.config)
    config = load_config(args.config)

    # Resolve output directory structure: results/{job_name}/01_data_acquisition/
    job_name = config.get("job_name", "default_job")
    results_base = Path(config.get("results_dir", "results")) / job_name
    step1_dir = results_base / "01_data_acquisition"

    output_cfg = config.get("output", {})
    metadata_tsv = Path(output_cfg.get("metadata_tsv", step1_dir / "proteomes_metadata.tsv"))
    selected_tsv = Path(step1_dir / "selected_proteomes.tsv")
    fasta_dir = Path(output_cfg.get("fasta_dir", step1_dir / "fasta"))
    summary_file = Path(output_cfg.get("summary_report", step1_dir / "acquisition_summary.json"))
    download_opts = config.get("download_options", {})
    filter_opts = config.get("filter_options", {})

    fasta_dir.mkdir(parents=True, exist_ok=True)

    # 2. Fetch All Proteomes Metadata from UniProt
    url = get_uniprot_metadata_url(config)
    df = fetch_and_parse_proteomes(url, timeout=download_opts.get("timeout_seconds", 90))

    if df.empty:
        raise ValueError(
            f"No proteomes found matching the query criteria!\n"
            f"Query URL used: {url}\n"
            f"Please check your taxonomy_id and query_params in {args.config}."
        )

    # 3. Apply Local Filtering
    logger.info("Applying local filtering on retrieved proteome metadata...")
    df = apply_local_filtering(df, filter_opts)

    # Save full metadata TSV (with 'Selected' flag)
    metadata_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(metadata_tsv, sep="\t", index=False)
    logger.info("Saved complete metadata table (%d proteomes) to %s", len(df), metadata_tsv)

    # Filter down to selected proteomes
    selected_df = df[df["Selected"] == True].copy()

    # Cap number of proteomes if max_proteomes is configured (e.g. For testing)
    max_proteomes = filter_opts.get("max_proteomes")
    if max_proteomes is not None and int(max_proteomes) > 0:
        max_p = int(max_proteomes)
        is_ref = selected_df["Proteome Type"].str.lower().apply(lambda x: 1 if "reference" in x and "non" not in x else 0)
        selected_df = selected_df.assign(_is_ref=is_ref).sort_values(by=["_is_ref", "Protein Count"], ascending=[False, False]).head(max_p).drop(columns=["_is_ref"])
        logger.info("Capped selection to %d proteomes (max_proteomes = %d for testing)", len(selected_df), max_p)

    selected_df.to_csv(selected_tsv, sep="\t", index=False)
    logger.info("Saved selected proteomes table (%d proteomes) to %s", len(selected_df), selected_tsv)

    ref_count = len(selected_df[selected_df["Proteome Type"].str.lower().str.contains("reference") & ~selected_df["Proteome Type"].str.lower().str.contains("non")])
    non_ref_count = len(selected_df[selected_df["Proteome Type"].str.lower().str.contains("non")])
    logger.info(
        "Selection breakdown: %d Reference proteome(s), %d Non-Reference proteome(s).",
        ref_count,
        non_ref_count,
    )

    if selected_df.empty:
        if filter_opts.get("reference_only", False) or config.get("query_params", {}).get("reference_only", False):
            raise ValueError(
                "Local filtering resulted in 0 selected proteomes! No Reference proteomes met the criteria for this taxon. "
                "Consider setting reference_only: false to include high-quality non-reference proteomes."
            )
        raise ValueError("Local filtering resulted in 0 selected proteomes! Check filter_options in config.yaml.")

    # 4. Download FASTA files for Selected Proteomes
    results = []
    if download_opts.get("download_fasta", True):
        raw_workers = download_opts.get("max_workers", "auto")
        max_workers = resolve_threads(raw_workers, config=config, default=8, task_name="fasta_downloads")
        # Cap network download threads at 16 to avoid remote server throttling/socket exhaustion
        max_workers = min(max_workers, 16)
        timeout = download_opts.get("timeout_seconds", 90)
        retries = download_opts.get("retry_attempts", 3)
        skip_existing = download_opts.get("skip_existing", True)

        upids_to_download = selected_df["Proteome Id"].dropna().astype(str).tolist()
        logger.info(
            "Starting FASTA downloads for %d selected proteomes into '%s' using %d worker(s)...",
            len(upids_to_download),
            fasta_dir,
            max_workers,
        )

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_upid = {}
            for idx, row in selected_df.iterrows():
                u = str(row["Proteome Id"])
                expected = int(row["Protein Count"]) if "Protein Count" in selected_df.columns and pd.notna(row["Protein Count"]) else None
                f = executor.submit(
                    download_single_fasta,
                    upid=u,
                    output_dir=fasta_dir,
                    expected_count=expected,
                    skip_existing=skip_existing,
                    timeout=timeout,
                    retry_attempts=retries,
                )
                future_to_upid[f] = u

            completed_count = 0
            for future in as_completed(future_to_upid):
                res = future.result()
                results.append(res)
                completed_count += 1
                if completed_count % 10 == 0 or completed_count == len(upids_to_download):
                    logger.info(
                        "Progress: %d/%d proteomes processed (Latest: %s [%s, %d seqs])",
                        completed_count,
                        len(upids_to_download),
                        res["upid"],
                        res["status"],
                        res["seq_count"],
                    )

    # 5. Rigorous FASTA Integrity Check: ensure no 0-byte or 0-sequence files
    logger.info("Verifying FASTA file integrity for all %d selected proteomes...", len(selected_df))
    empty_or_missing = []
    verified_proteomes = []

    for idx, row in selected_df.iterrows():
        u = str(row["Proteome Id"])
        fasta_file = fasta_dir / f"{u}.fasta"
        if not fasta_file.exists() or fasta_file.stat().st_size == 0:
            empty_or_missing.append(u)
        else:
            n_seqs = count_fasta_sequences(fasta_file)
            if n_seqs == 0:
                empty_or_missing.append(u)
                try:
                    fasta_file.unlink()  # remove empty/corrupted file
                except OSError:
                    pass
            else:
                verified_proteomes.append((u, n_seqs))

    auto_prune = filter_opts.get("auto_prune_empty", True)

    if empty_or_missing:
        if auto_prune:
            logger.warning("==================================================")
            logger.warning("AUTO-PRUNE TRIGGERED: Found %d proteome(s) with empty or missing FASTA:", len(empty_or_missing))
            for u in empty_or_missing:
                logger.warning("  --> Pruned Proteome: %s (No sequence data retrieved from UniProtKB)", u)
            logger.warning("These proteome(s) will be logged and excluded from downstream analysis.")
            logger.warning("==================================================")

            # Exclude pruned proteomes from selected table
            selected_df = selected_df[~selected_df["Proteome Id"].isin(empty_or_missing)].copy()
            df.loc[df["Proteome Id"].isin(empty_or_missing), "Selected"] = False

            # Update TSV files on disk
            selected_df.to_csv(selected_tsv, sep="\t", index=False)
            df.to_csv(metadata_tsv, sep="\t", index=False)

            if selected_df.empty:
                raise RuntimeError("All selected proteomes were pruned! No valid FASTA files available to proceed.")

            logger.info("Auto-pruning applied. Proceeding with %d verified proteome(s).", len(selected_df))
        else:
            logger.error(
                "FASTA Validation Failed! The following %d proteome(s) have empty or missing FASTA files: %s",
                len(empty_or_missing),
                empty_or_missing,
            )
            raise RuntimeError(
                f"FASTA integrity check failed: {len(empty_or_missing)} proteome(s) are empty or missing: {empty_or_missing}"
            )

    logger.info(
        "FASTA Integrity Check Passed: All %d selected FASTA files contain verified sequences (Total: %d seqs).",
        len(verified_proteomes),
        sum(s[1] for s in verified_proteomes),
    )

    # 6. Generate Summary Report
    total_seqs = sum(s[1] for s in verified_proteomes)
    successful = len(verified_proteomes)
    failed = [r for r in results if r.get("status") == "failed"]

    summary = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config_file": args.config,
        "query_url": url,
        "total_proteomes_retrieved": len(df),
        "total_proteomes_selected": len(selected_df),
        "selected_reference_proteomes": ref_count,
        "selected_non_reference_proteomes": non_ref_count,
        "total_fasta_downloaded": successful,
        "total_fasta_failed": len(failed),
        "total_fasta_pruned": len(empty_or_missing),
        "pruned_proteomes": empty_or_missing,
        "total_protein_sequences": total_seqs,
        "fasta_directory": str(fasta_dir),
        "metadata_file": str(metadata_tsv),
        "selected_file": str(selected_tsv),
        "failed_upids": [f["upid"] for f in failed],
    }

    summary_file.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    if download_opts.get("download_fasta", True):
        (fasta_dir / ".download_complete").touch()

    logger.info("==================================================")
    logger.info("Step 1 Acquisition Completed Successfully!")
    logger.info("Total Retrieved: %d | Total Selected: %d", len(df), len(selected_df))
    logger.info("Reference Proteomes Kept: %d", ref_count)
    logger.info("Non-Reference Proteomes Kept: %d", non_ref_count)
    if empty_or_missing:
        logger.info("Auto-Pruned Empty Proteomes: %d (%s)", len(empty_or_missing), empty_or_missing)
    logger.info("Successfully Downloaded & Verified: %d", successful)
    logger.info("Total Protein Sequences: %d", total_seqs)
    if failed:
        logger.warning("Failed downloads: %d (%s)", len(failed), [f['upid'] for f in failed])
    logger.info("Summary report saved to %s", summary_file)
    logger.info("==================================================")


if __name__ == "__main__":
    main()
