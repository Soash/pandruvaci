#!/usr/bin/env python3
"""
Step 3: Pangenome Analysis & Partitioning (MMseqs2)
===================================================
Clusters non-paralogous proteome sequences across all strains using MMseqs2,
computes the strain presence/absence matrix, partitions gene clusters into
Core, Shell, and Cloud genomes, and extracts Core Representative Proteins
for downstream drug target discovery.

Outputs:
  - matrix/pangenome_matrix.tsv
  - clusters/pangenome_clusters.tsv
  - fasta/core_representatives.fasta
  - fasta/core_all_members.fasta
  - fasta/shell_representatives.fasta
  - fasta/cloud_representatives.fasta
  - fasta/all_pangenome_representatives.fasta
  - summary/pangenome_summary.tsv
  - summary/pangenome_summary.json
  - summary/pangenome_distribution.txt
  - .pangenome_complete

Usage:
  python scripts/03_pangenome_analysis.py [--config config.yaml]
"""

import os
import sys
import re
import time
import json
import shutil
import logging
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict, OrderedDict
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
logger = logging.getLogger("03_pangenome_analysis")


def ensure_mmseqs_binary() -> str:
    """Finds or confirms mmseqs binary from PATH or project bin/."""
    which_mmseqs = shutil.which("mmseqs")
    if which_mmseqs:
        return which_mmseqs

    project_root = Path(__file__).resolve().parent.parent
    local_bin = project_root / "bin" / "mmseqs"
    if local_bin.exists() and os.access(local_bin, os.X_OK):
        return str(local_bin)

    raise RuntimeError(
        "MMseqs2 binary not found in PATH or project bin/. Please ensure bin/mmseqs exists."
    )


def load_config(config_path: str) -> dict:
    """Load configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def parse_fasta_header(header_line: str):
    """
    Parses a UniProt FASTA header line.
    Example: >tr|A0A7U3Y4P7|A0A7U3Y4P7_ACIB5 RNA polymerase sigma-54 factor OS=Acinetobacter baumannii OX=480119 GN=AB57_0729 PE=3 SV=1
    Returns:
      seq_id: 'tr|A0A7U3Y4P7|A0A7U3Y4P7_ACIB5'
      accession: 'A0A7U3Y4P7'
      gene_name: 'AB57_0729' (or '')
      protein_name: 'RNA polymerase sigma-54 factor'
    """
    clean_line = header_line.lstrip(">").strip()
    parts = clean_line.split(maxsplit=1)
    seq_id = parts[0]
    rest = parts[1] if len(parts) > 1 else ""

    # Accession extraction (between pipes if present)
    pipes = seq_id.split("|")
    accession = pipes[1] if len(pipes) >= 2 else seq_id

    # Gene name extraction
    gn_match = re.search(r"\bGN=([^\s]+)", rest)
    gene_name = gn_match.group(1) if gn_match else ""

    # Protein name extraction (up to OS= or first metadata tag)
    tag_match = re.search(r"\s+[A-Z]{2}=", rest)
    if tag_match:
        protein_name = rest[:tag_match.start()].strip()
    else:
        protein_name = rest.strip()

    return seq_id, accession, gene_name, protein_name


def pool_non_paralogous_fastas(
    np_fasta_files: list, pooled_fasta_path: Path
) -> tuple:
    """
    Combines all non-paralogous FASTA files into one pooled file,
    tagging each header with {upid}___{accession} so that MMseqs2
    unambiguously preserves the strain association and unique ID.

    Returns:
      seq_metadata: dict mapping pooled_id -> dict(upid, raw_id, accession, gene_name, protein_name, original_header)
      seq_records: dict mapping pooled_id -> sequence_str
      strain_seq_counts: dict mapping upid -> count
    """
    seq_metadata = {}
    seq_records = {}
    strain_seq_counts = defaultdict(int)

    logger.info("Pooling %d non-paralogous FASTA files...", len(np_fasta_files))
    pooled_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    with open(pooled_fasta_path, "w", encoding="utf-8") as out_f:
        for fasta_file in np_fasta_files:
            upid = fasta_file.name.replace("_np.fasta", "").replace(".fasta", "")
            current_pooled_id = None
            current_header = ""
            current_seq_chunks = []

            with open(fasta_file, "r", encoding="utf-8", errors="replace") as in_f:
                for line in in_f:
                    line_s = line.strip()
                    if not line_s:
                        continue
                    if line_s.startswith(">"):
                        # Flush previous sequence
                        if current_pooled_id is not None:
                            full_seq = "".join(current_seq_chunks)
                            seq_records[current_pooled_id] = full_seq
                            out_f.write(f">{current_pooled_id}\n")
                            out_f.write(f"{full_seq}\n")
                            strain_seq_counts[upid] += 1

                        current_header = line_s.lstrip(">")
                        seq_id, acc, gn, pn = parse_fasta_header(line_s)
                        current_pooled_id = f"{upid}___{acc}"
                        seq_metadata[current_pooled_id] = {
                            "upid": upid,
                            "raw_id": seq_id,
                            "accession": acc,
                            "gene_name": gn,
                            "protein_name": pn,
                            "original_header": current_header,
                        }
                        current_seq_chunks = []
                    else:
                        current_seq_chunks.append(line_s)

                # Flush trailing sequence
                if current_pooled_id is not None:
                    full_seq = "".join(current_seq_chunks)
                    seq_records[current_pooled_id] = full_seq
                    out_f.write(f">{current_pooled_id}\n")
                    out_f.write(f"{full_seq}\n")
                    strain_seq_counts[upid] += 1

    logger.info(
        "Successfully pooled %d sequences across %d strains into %s",
        len(seq_records),
        len(strain_seq_counts),
        pooled_fasta_path,
    )
    return seq_metadata, seq_records, strain_seq_counts


def run_mmseqs_clustering(
    mmseqs_bin: str,
    pooled_fasta: Path,
    out_prefix: Path,
    tmp_dir: Path,
    min_seq_id: float = 0.6,
    coverage: float = 0.8,
    cov_mode: int = 0,
    cluster_mode: str = "easy-cluster",
    threads: int = 4,
) -> Path:
    """Runs MMseqs2 clustering on the pooled FASTA."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cluster_tsv = Path(f"{out_prefix}_cluster.tsv")

    cmd = [
        mmseqs_bin,
        cluster_mode,
        str(pooled_fasta),
        str(out_prefix),
        str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
        "--threads", str(threads),
        "-v", "1",
    ]

    logger.info("Executing MMseqs2 clustering:")
    logger.info("Command: %s", " ".join(cmd))
    start_time = time.time()

    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        logger.error("MMseqs2 error: %s", res.stderr)
        raise RuntimeError(f"MMseqs2 clustering failed: {res.stderr}")

    duration = time.time() - start_time
    logger.info("MMseqs2 clustering completed in %.2f seconds.", duration)

    if not cluster_tsv.exists():
        raise FileNotFoundError(f"Expected MMseqs2 cluster TSV missing: {cluster_tsv}")

    return cluster_tsv


def compute_pangenome_matrix_and_partitions(
    cluster_tsv_path: Path,
    seq_metadata: dict,
    seq_records: dict,
    strain_list: list,
    core_threshold_pct: float = 100.0,
    shell_threshold_pct: float = 15.0,
):
    """
    Parses MMseqs2 clustering results, computes the presence/absence matrix,
    and categorizes each cluster into Core, Shell, or Cloud.
    """
    total_strains = len(strain_list)

    # 1. Parse clusters: rep_id -> list of member_ids
    clusters = defaultdict(list)
    with open(cluster_tsv_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rep_id, member_id = parts[0], parts[1]
                clusters[rep_id].append(member_id)

    logger.info(
        "Parsed %d orthologous clusters from MMseqs2 output.", len(clusters)
    )

    # 2. Build cluster records and matrix rows
    matrix_rows = []
    cluster_tidy_rows = []
    category_counts = {"Core": 0, "Shell": 0, "Cloud": 0}
    category_proteins = {"Core": 0, "Shell": 0, "Cloud": 0}

    # Partitioned sequence lists
    rep_seqs_by_cat = defaultdict(list)  # cat -> list of (header, seq)
    all_seqs_by_cat = defaultdict(list)  # cat -> list of (header, seq)
    all_rep_seqs = []                    # list of (header, seq)

    # Sort clusters deterministically (by size descending, then rep_id)
    sorted_reps = sorted(
        clusters.keys(),
        key=lambda r: (len(clusters[r]), r),
        reverse=True,
    )

    for cluster_num, rep_id in enumerate(sorted_reps, start=1):
        cluster_id = f"CL_{cluster_num:06d}"
        members = clusters[rep_id]
        rep_meta = seq_metadata.get(rep_id, {})
        rep_acc = rep_meta.get("accession", rep_id)
        rep_upid = rep_meta.get("upid", "Unknown")
        gene_name = rep_meta.get("gene_name", "")
        protein_name = rep_meta.get("protein_name", "")

        # Strain distribution
        strain_members = defaultdict(list)
        for m_id in members:
            m_upid = seq_metadata.get(m_id, {}).get("upid", m_id.split("___")[0] if "___" in m_id else "Unknown")
            strain_members[m_upid].append(m_id)

        present_strains = set(strain_members.keys())
        strain_count = len(present_strains)
        total_seqs = len(members)
        presence_pct = (strain_count / total_strains) * 100.0

        # Classification
        if presence_pct >= core_threshold_pct:
            category = "Core"
        elif presence_pct >= shell_threshold_pct:
            category = "Shell"
        else:
            category = "Cloud"

        category_counts[category] += 1
        category_proteins[category] += total_seqs

        # Build matrix row
        matrix_row = OrderedDict([
            ("cluster_id", cluster_id),
            ("representative_id", rep_id),
            ("accession", rep_acc),
            ("gene_name", gene_name),
            ("protein_name", protein_name),
            ("representative_strain", rep_upid),
            ("category", category),
            ("strain_count", strain_count),
            ("presence_pct", round(presence_pct, 2)),
            ("total_sequences", total_seqs),
        ])

        # Strain copy counts
        for s in strain_list:
            matrix_row[s] = len(strain_members.get(s, []))

        matrix_rows.append(matrix_row)

        # Sequences for FASTA export: standard header starting with accession for downstream tool compatibility
        rep_seq = seq_records.get(rep_id, "")
        rep_header = f"{rep_acc} cluster={cluster_id} category={category} strain={rep_upid} strains_present={strain_count}/{total_strains} presence_pct={presence_pct:.1f}% GN={gene_name} {protein_name}".strip()
        all_rep_seqs.append((rep_header, rep_seq))
        rep_seqs_by_cat[category].append((rep_header, rep_seq))

        # Tidy cluster mappings
        for m_id in members:
            m_meta = seq_metadata.get(m_id, {})
            m_upid = m_meta.get("upid", m_id.split("___")[0] if "___" in m_id else "Unknown")
            m_acc = m_meta.get("accession", m_id)
            m_header = m_meta.get("original_header", m_id)
            is_rep = (m_id == rep_id)
            cluster_tidy_rows.append({
                "cluster_id": cluster_id,
                "representative_id": rep_id,
                "representative_accession": rep_acc,
                "member_id": m_id,
                "member_accession": m_acc,
                "strain_id": m_upid,
                "is_representative": is_rep,
                "category": category,
                "gene_name": m_meta.get("gene_name", ""),
                "protein_name": m_meta.get("protein_name", ""),
            })

            m_seq = seq_records.get(m_id, "")
            all_seqs_by_cat[category].append((m_header, m_seq))

    return (
        matrix_rows,
        cluster_tidy_rows,
        category_counts,
        category_proteins,
        rep_seqs_by_cat,
        all_seqs_by_cat,
        all_rep_seqs,
    )


def write_fasta_file(filepath: Path, seq_tuples: list):
    """Writes list of (header, seq) tuples to a FASTA file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        for header, seq in seq_tuples:
            f.write(f">{header}\n")
            # Write sequence wrapped at 80 chars
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def generate_ascii_distribution(
    total_clusters: int, category_counts: dict, total_proteins: int, category_proteins: dict
) -> str:
    """Generates an ASCII bar chart summarizing pangenome categories."""
    lines = []
    lines.append("=" * 70)
    lines.append("                  PANGENOME PARTITION DISTRIBUTION")
    lines.append("=" * 70)
    lines.append(f"{'Category':<10} | {'Clusters':<10} | {'% of Clusters':<14} | {'Proteins':<10} | {'% of Proteins':<14}")
    lines.append("-" * 70)

    for cat in ["Core", "Shell", "Cloud"]:
        c_cnt = category_counts.get(cat, 0)
        c_pct = (c_cnt / total_clusters * 100.0) if total_clusters > 0 else 0.0
        p_cnt = category_proteins.get(cat, 0)
        p_pct = (p_cnt / total_proteins * 100.0) if total_proteins > 0 else 0.0
        bar = "█" * int(round(c_pct / 2.5))
        lines.append(f"{cat:<10} | {c_cnt:<10} | {c_pct:>6.2f}%        | {p_cnt:<10} | {p_pct:>6.2f}%")
        lines.append(f"  [{bar:<40}]")

    lines.append("-" * 70)
    lines.append(f"{'TOTAL':<10} | {total_clusters:<10} | 100.00%        | {total_proteins:<10} | 100.00%")
    lines.append("=" * 70)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Step 3: Pangenome Analysis & Partitioning using MMseqs2"
    )
    parser.add_argument(
        "--config", default="config.yaml", help="Path to configuration YAML"
    )
    args = parser.parse_args()

    # Load configuration
    cfg = load_config(args.config)
    job_name = cfg.get("job_name", "default_job")
    results_dir = Path(cfg.get("results_dir", "results"))
    job_dir = results_dir / job_name

    step2_dir = job_dir / "02_paralog_removal"
    step2_fasta_dir = step2_dir / "fasta"

    step3_dir = job_dir / "03_pangenome_analysis"
    intermediate_dir = step3_dir / "intermediate"
    matrix_dir = step3_dir / "matrix"
    clusters_dir = step3_dir / "clusters"
    fasta_dir = step3_dir / "fasta"
    summary_dir = step3_dir / "summary"

    for d in [intermediate_dir, matrix_dir, clusters_dir, fasta_dir, summary_dir]:
        d.mkdir(parents=True, exist_ok=True)

    step3_cfg = cfg.get("step3_pangenome", {})
    if not step3_cfg.get("enabled", True):
        logger.info("Step 3 is disabled in config.yaml. Exiting.")
        sys.exit(0)

    min_seq_id = float(step3_cfg.get("min_seq_id", 0.6))
    coverage = float(step3_cfg.get("coverage", 0.8))
    cov_mode = int(step3_cfg.get("cov_mode", 0))
    cluster_mode = str(step3_cfg.get("cluster_mode", "easy-cluster"))
    raw_threads = step3_cfg.get("threads", "auto")
    threads = resolve_threads(raw_threads, config=cfg, default=32, task_name="pangenome_clustering")
    core_threshold_pct = float(step3_cfg.get("core_threshold_pct", 100.0))
    shell_threshold_pct = float(step3_cfg.get("shell_threshold_pct", 15.0))

    logger.info("=" * 70)
    logger.info("Step 3: Pangenome Analysis & Partitioning (MMseqs2)")
    logger.info("=" * 70)
    logger.info("Job Name:             %s", job_name)
    logger.info("Hardware Threads:     %d threads allocated", threads)
    logger.info("Identity Cutoff:      %.2f", min_seq_id)
    logger.info("Coverage Cutoff:      %.2f (Mode %d)", coverage, cov_mode)
    logger.info("Clustering Mode:      %s", cluster_mode)
    logger.info("Core Threshold:       >= %.1f%% of strains", core_threshold_pct)
    logger.info("Shell Threshold:      >= %.1f%% and < %.1f%%", shell_threshold_pct, core_threshold_pct)
    logger.info("Cloud Threshold:      < %.1f%%", shell_threshold_pct)
    logger.info("=" * 70)

    mmseqs_bin = ensure_mmseqs_binary()
    logger.info("Using MMseqs2 executable: %s", mmseqs_bin)

    # 1. Discover input non-paralogous FASTAs
    np_fastas = sorted(list(step2_fasta_dir.glob("*_np.fasta")))
    if not np_fastas:
        raise FileNotFoundError(
            f"No non-paralogous FASTA files found in {step2_fasta_dir}. "
            "Please ensure Step 2 (paralog_removal) has completed."
        )

    # Extract strain IDs
    strain_list = [f.name.replace("_np.fasta", "") for f in np_fastas]
    logger.info("Found %d non-paralogous proteome files: %s", len(strain_list), strain_list)

    # 2. Pool non-paralogous sequences with strain-tagged IDs
    pooled_fasta = intermediate_dir / "pooled_non_paralogous.fasta"
    seq_metadata, seq_records, strain_seq_counts = pool_non_paralogous_fastas(
        np_fastas, pooled_fasta
    )

    # 3. Execute MMseqs2 cross-strain clustering
    mmseqs_out_prefix = intermediate_dir / "mmseqs_pan"
    mmseqs_tmp = intermediate_dir / "mmseqs_tmp"

    cluster_tsv = run_mmseqs_clustering(
        mmseqs_bin=mmseqs_bin,
        pooled_fasta=pooled_fasta,
        out_prefix=mmseqs_out_prefix,
        tmp_dir=mmseqs_tmp,
        min_seq_id=min_seq_id,
        coverage=coverage,
        cov_mode=cov_mode,
        cluster_mode=cluster_mode,
        threads=threads,
    )

    # 4. Compute Pangenome Matrix & Partitions
    (
        matrix_rows,
        cluster_tidy_rows,
        category_counts,
        category_proteins,
        rep_seqs_by_cat,
        all_seqs_by_cat,
        all_rep_seqs,
    ) = compute_pangenome_matrix_and_partitions(
        cluster_tsv_path=cluster_tsv,
        seq_metadata=seq_metadata,
        seq_records=seq_records,
        strain_list=strain_list,
        core_threshold_pct=core_threshold_pct,
        shell_threshold_pct=shell_threshold_pct,
    )

    total_clusters = len(matrix_rows)
    total_proteins = len(seq_records)

    # 5. Export Pangenome Matrix
    matrix_df = pd.DataFrame(matrix_rows)
    matrix_tsv_path = matrix_dir / "pangenome_matrix.tsv"
    matrix_df.to_csv(matrix_tsv_path, sep="\t", index=False)
    logger.info("Saved pangenome matrix TSV: %s (%d rows)", matrix_tsv_path, len(matrix_df))

    # 6. Export Tidy Cluster Protein Table
    cluster_df = pd.DataFrame(cluster_tidy_rows)
    cluster_tsv_path = clusters_dir / "pangenome_clusters.tsv"
    cluster_df.to_csv(cluster_tsv_path, sep="\t", index=False)
    logger.info("Saved tidy cluster TSV: %s (%d rows)", cluster_tsv_path, len(cluster_df))

    # 7. Export Partitioned FASTAs
    core_rep_path = fasta_dir / "core_representatives.fasta"
    core_all_path = fasta_dir / "core_all_members.fasta"
    shell_rep_path = fasta_dir / "shell_representatives.fasta"
    cloud_rep_path = fasta_dir / "cloud_representatives.fasta"
    all_rep_path = fasta_dir / "all_pangenome_representatives.fasta"

    write_fasta_file(core_rep_path, rep_seqs_by_cat.get("Core", []))
    write_fasta_file(core_all_path, all_seqs_by_cat.get("Core", []))
    write_fasta_file(shell_rep_path, rep_seqs_by_cat.get("Shell", []))
    write_fasta_file(cloud_rep_path, rep_seqs_by_cat.get("Cloud", []))
    write_fasta_file(all_rep_path, all_rep_seqs)

    logger.info("Saved core representatives:  %s (%d seqs)", core_rep_path, len(rep_seqs_by_cat.get("Core", [])))
    logger.info("Saved core all members:      %s (%d seqs)", core_all_path, len(all_seqs_by_cat.get("Core", [])))
    logger.info("Saved shell representatives: %s (%d seqs)", shell_rep_path, len(rep_seqs_by_cat.get("Shell", [])))
    logger.info("Saved cloud representatives: %s (%d seqs)", cloud_rep_path, len(rep_seqs_by_cat.get("Cloud", [])))

    # 8. Export Summary Files
    summary_stats = [
        {"Category": "Core", "Clusters": category_counts["Core"], "Cluster_Pct": round(category_counts["Core"] / total_clusters * 100, 2), "Proteins": category_proteins["Core"], "Protein_Pct": round(category_proteins["Core"] / total_proteins * 100, 2)},
        {"Category": "Shell", "Clusters": category_counts["Shell"], "Cluster_Pct": round(category_counts["Shell"] / total_clusters * 100, 2), "Proteins": category_proteins["Shell"], "Protein_Pct": round(category_proteins["Shell"] / total_proteins * 100, 2)},
        {"Category": "Cloud", "Clusters": category_counts["Cloud"], "Cluster_Pct": round(category_counts["Cloud"] / total_clusters * 100, 2), "Proteins": category_proteins["Cloud"], "Protein_Pct": round(category_proteins["Cloud"] / total_proteins * 100, 2)},
        {"Category": "Total", "Clusters": total_clusters, "Cluster_Pct": 100.0, "Proteins": total_proteins, "Protein_Pct": 100.0},
    ]
    summary_df = pd.DataFrame(summary_stats)
    summary_tsv_path = summary_dir / "pangenome_summary.tsv"
    summary_df.to_csv(summary_tsv_path, sep="\t", index=False)

    summary_json_path = summary_dir / "pangenome_summary.json"
    summary_dict = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config_file": args.config,
        "mmseqs_binary": mmseqs_bin,
        "parameters": {
            "min_seq_id": min_seq_id,
            "coverage": coverage,
            "cov_mode": cov_mode,
            "cluster_mode": cluster_mode,
            "core_threshold_pct": core_threshold_pct,
            "shell_threshold_pct": shell_threshold_pct,
            "threads": threads,
        },
        "total_strains": len(strain_list),
        "total_input_proteins": total_proteins,
        "total_orthologous_clusters": total_clusters,
        "category_counts": category_counts,
        "category_proteins": category_proteins,
        "strains": {s: strain_seq_counts[s] for s in strain_list},
        "outputs": {
            "matrix_tsv": str(matrix_tsv_path),
            "clusters_tsv": str(cluster_tsv_path),
            "core_representatives_fasta": str(core_rep_path),
            "core_all_fasta": str(core_all_path),
            "shell_representatives_fasta": str(shell_rep_path),
            "cloud_representatives_fasta": str(cloud_rep_path),
            "summary_tsv": str(summary_tsv_path),
        },
    }
    with open(summary_json_path, "w", encoding="utf-8") as jf:
        json.dump(summary_dict, jf, indent=2)

    # 9. Generate ASCII distribution
    ascii_dist = generate_ascii_distribution(
        total_clusters, category_counts, total_proteins, category_proteins
    )
    dist_path = summary_dir / "pangenome_distribution.txt"
    with open(dist_path, "w", encoding="utf-8") as df:
        df.write(ascii_dist + "\n")

    print("\n" + ascii_dist + "\n")

    # 10. Write sentinel file
    sentinel_path = step3_dir / ".pangenome_complete"
    with open(sentinel_path, "w", encoding="utf-8") as f:
        f.write(f"Completed at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    logger.info("Step 3 completed successfully! Sentinel created: %s", sentinel_path)


if __name__ == "__main__":
    main()
