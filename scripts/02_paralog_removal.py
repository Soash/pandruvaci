#!/usr/bin/env python3
"""
Step 2: Redundancy Reduction & Paralog Removal (MMseqs2)
======================================================
Removes paralogous / redundant sequences from each strain's proteome
using MMseqs2 easy-cluster / easy-linclust.

If MMseqs2 is not found in PATH, the script automatically installs it
into bin/mmseqs using the official precompiled AVX2/SSE4 static binary.

Usage:
    python scripts/02_paralog_removal.py [--config config.yaml]
"""

import os
import sys
import time
import json
import shutil
import logging
import argparse
import tarfile
import tempfile
import urllib.request
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import yaml

try:
    from scripts.utils_hardware import optimize_parallel_workers_and_threads, resolve_threads
except ImportError:
    from utils_hardware import optimize_parallel_workers_and_threads, resolve_threads

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("02_paralog_removal")


def ensure_mmseqs_binary() -> str:
    """
    Checks if mmseqs is available in PATH or project bin/.
    If not found, automatically downloads and installs the official static binary.
    """
    # 1. Check existing PATH
    which_mmseqs = shutil.which("mmseqs")
    if which_mmseqs:
        logger.debug("Found MMseqs2 in PATH: %s", which_mmseqs)
        return which_mmseqs

    # 2. Check project bin/mmseqs
    project_root = Path(__file__).resolve().parent.parent
    local_bin = project_root / "bin" / "mmseqs"
    if local_bin.exists() and os.access(local_bin, os.X_OK):
        logger.debug("Found MMseqs2 in local bin: %s", local_bin)
        return str(local_bin)

    # 3. Auto-download static binary
    logger.info("MMseqs2 not found. Automatically downloading official static binary...")
    bin_dir = project_root / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    tar_path = bin_dir / "mmseqs-linux-avx2.tar.gz"

    download_url = "https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz"
    try:
        urllib.request.urlretrieve(download_url, str(tar_path))
        with tarfile.open(str(tar_path), "r:gz") as tar:
            for member in tar.getmembers():
                if member.name.endswith("bin/mmseqs"):
                    member.name = "mmseqs"
                    tar.extract(member, path=str(bin_dir))
        if tar_path.exists():
            tar_path.unlink()
        local_bin.chmod(0o755)
        logger.info("Successfully installed MMseqs2 to %s", local_bin)
        return str(local_bin)
    except Exception as e:
        logger.error("Failed to auto-install MMseqs2: %s", e)
        raise RuntimeError(
            f"MMseqs2 binary is required. Please install via 'conda install -c conda-forge -c bioconda mmseqs2' or check network."
        )


def load_config(config_path: str) -> dict:
    """Load configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def count_fasta_sequences(fasta_path: Path) -> int:
    """Count records starting with '>' in a FASTA file."""
    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        return 0
    count = 0
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def cluster_proteome_mmseqs(
    mmseqs_bin: str,
    input_fasta: Path,
    output_np_fasta: Path,
    output_cluster_tsv: Path,
    min_seq_id: float = 0.6,
    coverage: float = 0.8,
    cov_mode: int = 0,
    cluster_mode: str = "easy-cluster",
    threads: int = 4,
) -> dict:
    """
    Executes MMseqs2 clustering on a single proteome to remove paralogs.
    """
    upid = input_fasta.stem
    orig_count = count_fasta_sequences(input_fasta)

    result = {
        "upid": upid,
        "input_fasta": str(input_fasta),
        "output_np_fasta": str(output_np_fasta),
        "output_cluster_tsv": str(output_cluster_tsv),
        "original_count": orig_count,
        "non_paralog_count": 0,
        "paralogs_removed": 0,
        "reduction_pct": 0.0,
        "status": "pending",
        "error": None,
    }

    if orig_count == 0:
        result["status"] = "failed"
        result["error"] = "Input FASTA is empty"
        return result

    # Check if already processed
    if output_np_fasta.exists() and output_cluster_tsv.exists() and output_np_fasta.stat().st_size > 0:
        np_count = count_fasta_sequences(output_np_fasta)
        if np_count > 0:
            result["status"] = "skipped"
            result["non_paralog_count"] = np_count
            result["paralogs_removed"] = orig_count - np_count
            result["reduction_pct"] = round((orig_count - np_count) / orig_count * 100, 2) if orig_count else 0.0
            return result

    with tempfile.TemporaryDirectory(prefix=f"mmseqs_{upid}_") as tmpdir:
        tmp_path = Path(tmpdir)
        out_prefix = tmp_path / "clust"

        cmd = [
            mmseqs_bin,
            cluster_mode,
            str(input_fasta),
            str(out_prefix),
            str(tmp_path / "tmp"),
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
            "--threads", str(threads),
            "-v", "1",
        ]

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            result["status"] = "failed"
            result["error"] = proc.stderr.strip() or proc.stdout.strip()
            return result

        rep_fasta = tmp_path / "clust_rep_seq.fasta"
        clu_tsv = tmp_path / "clust_cluster.tsv"

        if not rep_fasta.exists() or rep_fasta.stat().st_size == 0:
            result["status"] = "failed"
            result["error"] = "MMseqs2 did not generate representative FASTA"
            return result

        # Move outputs to target locations
        output_np_fasta.parent.mkdir(parents=True, exist_ok=True)
        output_cluster_tsv.parent.mkdir(parents=True, exist_ok=True)

        shutil.copyfile(rep_fasta, output_np_fasta)
        if clu_tsv.exists():
            shutil.copyfile(clu_tsv, output_cluster_tsv)

        np_count = count_fasta_sequences(output_np_fasta)
        result["status"] = "success"
        result["non_paralog_count"] = np_count
        result["paralogs_removed"] = orig_count - np_count
        result["reduction_pct"] = round((orig_count - np_count) / orig_count * 100, 2) if orig_count else 0.0

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Step 2: Remove paralogs / redundant sequences per strain using MMseqs2."
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML configuration file (default: config.yaml)",
    )
    args = parser.parse_args()

    # 1. Ensure MMseqs2 binary
    mmseqs_bin = ensure_mmseqs_binary()

    # 2. Load Configuration
    config = load_config(args.config)
    job_name = config.get("job_name", "default_job")
    results_base = Path(config.get("results_dir", "results")) / job_name

    # Step 1 input paths
    step1_dir = results_base / "01_data_acquisition"
    step1_fasta_dir = step1_dir / "fasta"
    selected_tsv = step1_dir / "selected_proteomes.tsv"

    # Step 2 output paths
    step2_dir = results_base / "02_paralog_removal"
    step2_fasta_dir = step2_dir / "fasta"
    step2_cluster_dir = step2_dir / "clusters"
    step2_summary_dir = step2_dir / "summary"

    step2_fasta_dir.mkdir(parents=True, exist_ok=True)
    step2_cluster_dir.mkdir(parents=True, exist_ok=True)
    step2_summary_dir.mkdir(parents=True, exist_ok=True)

    stats_tsv = step2_summary_dir / "paralog_removal_stats.tsv"
    summary_json = step2_summary_dir / "paralog_removal_summary.json"
    done_flag = step2_fasta_dir / ".paralog_removal_complete"

    # Step 2 parameters
    p_opts = config.get("step2_paralog_removal", {})
    min_seq_id = float(p_opts.get("min_seq_id", 0.6))
    coverage = float(p_opts.get("coverage", 0.8))
    cov_mode = int(p_opts.get("cov_mode", 0))
    cluster_mode = str(p_opts.get("cluster_mode", "easy-cluster"))

    # 3. Identify input proteomes
    if selected_tsv.exists():
        df_sel = pd.read_csv(selected_tsv, sep="\t")
        upids = df_sel["Proteome Id"].dropna().astype(str).tolist()
    else:
        upids = [f.stem for f in step1_fasta_dir.glob("*.fasta")]

    if not upids:
        raise RuntimeError(f"No proteome FASTA files found in {step1_fasta_dir}")

    raw_threads = p_opts.get("threads", "auto")
    max_workers, threads_per_job = optimize_parallel_workers_and_threads(
        total_items=len(upids),
        step_threads=raw_threads,
        config=config,
        preferred_threads_per_worker=8,
    )

    logger.info("==================================================")
    logger.info("Step 2: MMseqs2 Paralog Removal Starting")
    logger.info("Proteomes to process: %d", len(upids))
    logger.info("Hardware Allocation: %d concurrent worker(s) x %d threads each (Total: %d cores)",
                max_workers, threads_per_job, max_workers * threads_per_job)
    logger.info("Parameters: min_seq_id=%.2f, coverage=%.2f, cov_mode=%d, mode=%s",
                min_seq_id, coverage, cov_mode, cluster_mode)
    logger.info("==================================================")

    # 4. Run clustering
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_upid = {}
        for u in upids:
            in_fasta = step1_fasta_dir / f"{u}.fasta"
            out_np = step2_fasta_dir / f"{u}_np.fasta"
            out_clu = step2_cluster_dir / f"{u}_cluster.tsv"

            if not in_fasta.exists():
                logger.warning("FASTA file missing for %s: %s", u, in_fasta)
                continue

            f = executor.submit(
                cluster_proteome_mmseqs,
                mmseqs_bin=mmseqs_bin,
                input_fasta=in_fasta,
                output_np_fasta=out_np,
                output_cluster_tsv=out_clu,
                min_seq_id=min_seq_id,
                coverage=coverage,
                cov_mode=cov_mode,
                cluster_mode=cluster_mode,
                threads=threads_per_job,
            )
            future_to_upid[f] = u

        for future in as_completed(future_to_upid):
            res = future.result()
            results.append(res)
            logger.info(
                "Proteome %s: %d -> %d seqs (-%d paralogs, -%.1f%%) [%s]",
                res["upid"],
                res["original_count"],
                res["non_paralog_count"],
                res["paralogs_removed"],
                res["reduction_pct"],
                res["status"],
            )

    # 5. Integrity verification
    empty_outputs = [r for r in results if r["status"] == "failed" or r["non_paralog_count"] == 0]
    if empty_outputs:
        failed_upids = [r["upid"] for r in empty_outputs]
        logger.error("Failed to generate non-paralogous FASTAs for: %s", failed_upids)
        raise RuntimeError(f"Step 2 failed for proteomes: {failed_upids}")

    # 6. Save statistics and summary
    stats_df = pd.DataFrame(results)[
        ["upid", "original_count", "non_paralog_count", "paralogs_removed", "reduction_pct", "status"]
    ]
    stats_df.rename(columns={
        "upid": "Proteome Id",
        "original_count": "Original Proteins",
        "non_paralog_count": "Non-Paralogous Proteins",
        "paralogs_removed": "Paralogs Removed",
        "reduction_pct": "Reduction (%)",
    }, inplace=True)
    stats_df.to_csv(stats_tsv, sep="\t", index=False)

    total_orig = stats_df["Original Proteins"].sum()
    total_np = stats_df["Non-Paralogous Proteins"].sum()
    total_removed = stats_df["Paralogs Removed"].sum()
    overall_red = round((total_removed / total_orig) * 100, 2) if total_orig else 0.0

    summary = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config_file": args.config,
        "mmseqs_binary": mmseqs_bin,
        "parameters": {
            "min_seq_id": min_seq_id,
            "coverage": coverage,
            "cov_mode": cov_mode,
            "cluster_mode": cluster_mode,
            "threads": threads_per_job,
        },
        "total_proteomes_processed": len(results),
        "total_original_sequences": int(total_orig),
        "total_non_paralogous_sequences": int(total_np),
        "total_paralogs_removed": int(total_removed),
        "overall_reduction_pct": overall_red,
        "output_fasta_dir": str(step2_fasta_dir),
        "output_cluster_dir": str(step2_cluster_dir),
        "stats_tsv": str(stats_tsv),
    }

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    done_flag.touch()

    logger.info("==================================================")
    logger.info("Step 2: Paralog Removal Completed Successfully!")
    logger.info("Proteomes Processed: %d", len(results))
    logger.info("Total Original Sequences: %d", total_orig)
    logger.info("Total Non-Paralogous Sequences: %d", total_np)
    logger.info("Total Paralogs Filtered Out: %d (%.2f%% reduction)", total_removed, overall_red)
    logger.info("Summary saved to %s", summary_json)
    logger.info("Stats TSV saved to %s", stats_tsv)
    logger.info("==================================================")


if __name__ == "__main__":
    main()
