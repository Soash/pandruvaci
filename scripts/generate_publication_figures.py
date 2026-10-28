#!/usr/bin/env python3
"""
Publication-Ready Visualizations & CSV Exports for Steps 1, 2, and 3
===================================================================
Generates high-resolution publication-quality figures (300 DPI PNG and vector PDF)
and standard CSV files for:
  - Step 1: Data Acquisition (CPD & Proteome distributions)
  - Step 2: Paralog Removal (Redundancy reduction per strain)
  - Step 3: Pangenome Partitioning (Core/Shell/Cloud distributions, strain frequency, matrix heatmap)
  - Master Multi-panel Figure: Complete overview of Steps 1-3

Usage:
  python scripts/generate_publication_figures.py [--config config.yaml]
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("generate_figures")

# Publication style parameters (Nature / Cell style)
PUB_PALETTE = {
    "Core": "#1b7837",      # Forest Green
    "Shell": "#e66101",     # Warm Amber / Orange
    "Cloud": "#5e3c99",     # Slate Purple
    "Reference": "#0571b0",  # Royal Blue
    "Non-Reference": "#92c5de",  # Sky Blue
    "Original": "#4575b4",   # Steel Blue
    "Non-Paralog": "#74add1", # Soft Blue
    "Paralogs": "#d73027",    # Crimson Red
}


def apply_publication_theme():
    """Applies clean publication-grade typography and layout settings."""
    sns.set_theme(style="ticks", font_scale=1.05)
    plt.rcParams.update({
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "axes.edgecolor": "#333333",
        "axes.linewidth": 1.1,
        "axes.labelcolor": "#111111",
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.labelsize": 11,
        "axes.labelweight": "medium",
        "xtick.color": "#222222",
        "ytick.color": "#222222",
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "grid.color": "#e0e0e0",
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "#cccccc",
        "legend.fontsize": 10,
        "figure.titlesize": 14,
        "figure.titleweight": "bold",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.12,
    })


def export_figure(fig, out_prefix: Path):
    """Saves figure in both 300 DPI PNG and vector PDF formats."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_prefix.with_suffix(".png")
    pdf_path = out_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    logger.info("Exported figures: %s and %s", png_path.name, pdf_path.name)


# ==============================================================================
# Step 1: Data Acquisition Figures & CSVs
# ==============================================================================
def process_step1_visuals(step1_dir: Path):
    """Generates publication charts and CSV exports for Step 1."""
    meta_tsv = step1_dir / "proteomes_metadata.tsv"
    sel_tsv = step1_dir / "selected_proteomes.tsv"
    fig_dir = step1_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not meta_tsv.exists():
        logger.warning("Step 1 metadata TSV not found at %s. Skipping Step 1 figures.", meta_tsv)
        return

    df_meta = pd.read_csv(meta_tsv, sep="\t")
    # Export CSVs
    meta_csv = step1_dir / "proteomes_metadata.csv"
    df_meta.to_csv(meta_csv, index=False)
    logger.info("Saved Step 1 metadata CSV: %s", meta_csv)

    if sel_tsv.exists():
        df_sel = pd.read_csv(sel_tsv, sep="\t")
        sel_csv = step1_dir / "selected_proteomes.csv"
        df_sel.to_csv(sel_csv, index=False)
        logger.info("Saved Step 1 selected proteomes CSV: %s", sel_csv)

    # Figure 1: Proteome Type & CPD Assessment Distribution
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.5))

    # Panel A: Proteome Type Distribution
    type_counts = df_meta["Proteome Type"].value_counts()
    colors_type = [PUB_PALETTE["Non-Reference"], PUB_PALETTE["Reference"]]
    bars = axes[0].bar(
        type_counts.index,
        type_counts.values,
        color=colors_type[:len(type_counts)],
        edgecolor="#333333",
        linewidth=1.0,
        width=0.55,
    )
    axes[0].set_title("A. Proteome Classification", loc="left")
    axes[0].set_ylabel("Number of Proteomes")
    axes[0].set_yscale("log")
    axes[0].set_ylim(1, 45000)
    axes[0].grid(axis="y", alpha=0.5)
    for bar in bars:
        h = bar.get_height()
        axes[0].text(
            bar.get_x() + bar.get_width() / 2,
            h * 1.25,
            f"{int(h):,}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    # Panel B: CPD Assessment Breakdown (Clean Donut with Side Legend)
    cpd_counts = df_meta["CPD"].value_counts()
    total_proteomes = len(df_meta)
    cpd_palette = {
        "Unknown": "#7f7f7f",
        "Standard": "#2ca02c",
        "Close to standard (low value)": "#66c2a4",
        "Close to standard (high value)": "#b2e2e2",
        "Outlier (high value)": "#ff7f0e",
        "Outlier (low value)": "#d62728",
    }
    colors_cpd = [cpd_palette.get(k, "#cccccc") for k in cpd_counts.index]

    wedges, _ = axes[1].pie(
        cpd_counts.values,
        startangle=140,
        colors=colors_cpd,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5),
    )
    # Donut center text
    axes[1].text(
        0, 0,
        f"Total\n{total_proteomes:,}",
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="bold",
        color="#333333",
    )
    # External legend to prevent any text overlap on small slices
    legend_labels = [
        f"{k}: {v:,}  ({v / total_proteomes * 100:.2f}%)"
        for k, v in cpd_counts.items()
    ]
    axes[1].legend(
        wedges,
        legend_labels,
        title="CPD Completeness Category",
        title_fontsize=10.5,
        loc="center left",
        bbox_to_anchor=(0.92, 0.5),
        frameon=True,
        fontsize=9.5,
    )
    axes[1].set_title("B. Complete Proteome Detector (CPD) Status", loc="left")

    plt.tight_layout()
    export_figure(fig, fig_dir / "step1_proteome_composition")

    # Figure 2: Protein Count Distribution
    if "Protein Count" in df_meta.columns and not df_meta["Protein Count"].dropna().empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        valid_counts = df_meta["Protein Count"].dropna()
        sns.histplot(
            valid_counts,
            kde=True,
            ax=ax,
            color="#2b8cbe",
            edgecolor="#1c5a7d",
            bins=40,
            alpha=0.65,
        )
        # Reference mean line
        mean_val = valid_counts.mean()
        median_val = valid_counts.median()
        ax.axvline(mean_val, color="#d95f02", linestyle="--", linewidth=1.5, label=f"Mean: {mean_val:.0f} aa")
        ax.axvline(median_val, color="#7570b3", linestyle=":", linewidth=1.5, label=f"Median: {median_val:.0f} aa")
        ax.set_title("Distribution of Proteome Protein Counts", loc="left")
        ax.set_xlabel("Protein Count per Proteome")
        ax.set_ylabel("Frequency")
        ax.legend(loc="upper right")
        ax.grid(axis="y", alpha=0.5)
        plt.tight_layout()
        export_figure(fig, fig_dir / "step1_protein_count_distribution")


# ==============================================================================
# Step 2: Paralog Removal Figures & CSVs
# ==============================================================================
def process_step2_visuals(step2_dir: Path):
    """Generates publication charts and CSV exports for Step 2."""
    stats_tsv = step2_dir / "summary" / "paralog_removal_stats.tsv"
    fig_dir = step2_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not stats_tsv.exists():
        logger.warning("Step 2 stats TSV not found at %s. Skipping Step 2 figures.", stats_tsv)
        return

    df_stats = pd.read_csv(stats_tsv, sep="\t")
    # Export CSV
    stats_csv = step2_dir / "summary" / "paralog_removal_stats.csv"
    df_stats.to_csv(stats_csv, index=False)
    logger.info("Saved Step 2 stats CSV: %s", stats_csv)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(df_stats))
    width = 0.38

    b1 = ax.bar(
        x - width / 2,
        df_stats["Original Proteins"],
        width,
        label="Original Proteins",
        color=PUB_PALETTE["Original"],
        edgecolor="#333333",
        linewidth=0.9,
    )
    b2 = ax.bar(
        x + width / 2,
        df_stats["Non-Paralogous Proteins"],
        width,
        label="Retained Non-Paralogous",
        color=PUB_PALETTE["Non-Paralog"],
        edgecolor="#333333",
        linewidth=0.9,
    )

    # Reduction annotations on top of bars
    for idx, row in df_stats.iterrows():
        red_pct = row["Reduction (%)"]
        ax.text(
            x[idx],
            row["Original Proteins"] + 40,
            f"-{red_pct:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
            color=PUB_PALETTE["Paralogs"],
        )

    ax.set_ylabel("Number of Proteins")
    ax.set_xlabel("Strain / Proteome Identifier")
    ax.set_title("Intra-Strain Redundancy Reduction via MMseqs2", loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels(df_stats["Proteome Id"], rotation=25, ha="right", fontweight="medium")
    ax.set_ylim(0, df_stats["Original Proteins"].max() * 1.15)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True)
    ax.grid(axis="y", alpha=0.5)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step2_paralog_reduction_barplot")


# ==============================================================================
# Step 3: Pangenome Partitioning Figures & CSVs
# ==============================================================================
def process_step3_visuals(step3_dir: Path):
    """Generates publication charts and CSV exports for Step 3."""
    matrix_tsv = step3_dir / "matrix" / "pangenome_matrix.tsv"
    clusters_tsv = step3_dir / "clusters" / "pangenome_clusters.tsv"
    summary_tsv = step3_dir / "summary" / "pangenome_summary.tsv"
    fig_dir = step3_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not matrix_tsv.exists():
        logger.warning("Step 3 matrix TSV not found at %s. Skipping Step 3 figures.", matrix_tsv)
        return

    df_matrix = pd.read_csv(matrix_tsv, sep="\t")
    # Export CSVs
    matrix_csv = step3_dir / "matrix" / "pangenome_matrix.csv"
    df_matrix.to_csv(matrix_csv, index=False)
    logger.info("Saved Step 3 matrix CSV: %s", matrix_csv)

    if clusters_tsv.exists():
        df_clusters = pd.read_csv(clusters_tsv, sep="\t")
        clusters_csv = step3_dir / "clusters" / "pangenome_clusters.csv"
        df_clusters.to_csv(clusters_csv, index=False)
        logger.info("Saved Step 3 clusters CSV: %s", clusters_csv)

    if summary_tsv.exists():
        df_summary = pd.read_csv(summary_tsv, sep="\t")
        summary_csv = step3_dir / "summary" / "pangenome_summary.csv"
        df_summary.to_csv(summary_csv, index=False)
        logger.info("Saved Step 3 summary CSV: %s", summary_csv)

    # --------------------------------------------------------------------------
    # Figure 3A: Pangenome Partition Nested Donut Chart (Clusters vs Proteins)
    # --------------------------------------------------------------------------
    cat_order = ["Core", "Shell", "Cloud"]
    cat_clusters = df_matrix["category"].value_counts().reindex(cat_order).fillna(0)
    cat_proteins = df_matrix.groupby("category")["total_sequences"].sum().reindex(cat_order).fillna(0)

    # Filter out categories with 0 for clean pie rendering
    valid_cats = [c for c in cat_order if cat_clusters[c] > 0]
    colors = [PUB_PALETTE[c] for c in valid_cats]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # Panel 1: Cluster distribution
    axes[0].pie(
        cat_clusters[valid_cats],
        labels=[f"{c}\n({int(cat_clusters[c]):,} clusters)" for c in valid_cats],
        autopct="%1.1f%%",
        pctdistance=0.72,
        startangle=140,
        colors=colors,
        wedgeprops=dict(width=0.42, edgecolor="#ffffff", linewidth=1.5),
    )
    axes[0].set_title("A. Orthologous Gene Clusters", loc="left")

    # Panel 2: Total proteins distribution
    axes[1].pie(
        cat_proteins[valid_cats],
        labels=[f"{c}\n({int(cat_proteins[c]):,} proteins)" for c in valid_cats],
        autopct="%1.1f%%",
        pctdistance=0.72,
        startangle=140,
        colors=colors,
        wedgeprops=dict(width=0.42, edgecolor="#ffffff", linewidth=1.5),
    )
    axes[1].set_title("B. Total Protein Sequences", loc="left")

    plt.tight_layout()
    export_figure(fig, fig_dir / "step3_pangenome_partition_donut")

    # --------------------------------------------------------------------------
    # Figure 3B: Strain Frequency Distribution (U-shaped pangenome curve)
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    strain_counts = df_matrix["strain_count"].value_counts().sort_index()

    # Determine colors for each bar based on category
    bar_colors = []
    for sc in strain_counts.index:
        sample_cat = df_matrix[df_matrix["strain_count"] == sc]["category"].iloc[0]
        bar_colors.append(PUB_PALETTE.get(sample_cat, "#6c757d"))

    bars = ax.bar(
        strain_counts.index,
        strain_counts.values,
        color=bar_colors,
        edgecolor="#333333",
        linewidth=1.0,
        width=0.65,
    )

    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + max(strain_counts.values) * 0.015,
            f"{int(h):,}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=9.5,
        )

    ax.set_ylim(0, max(strain_counts.values) * 1.20)
    ax.set_xlabel("Number of Strains Sharing Gene Cluster")
    ax.set_ylabel("Number of Orthologous Clusters")
    ax.set_title("Gene Cluster Frequency Distribution Across Strains", loc="left")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.grid(axis="y", alpha=0.5)

    # Custom legend outside the plot border
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PUB_PALETTE[c], edgecolor="#333333", label=c)
        for c in valid_cats
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        title="Pangenome Category",
        frameon=True,
        alignment="left",
    )

    plt.tight_layout()
    export_figure(fig, fig_dir / "step3_strain_frequency_distribution")

    # --------------------------------------------------------------------------
    # Figure 3C: Strain Proteome Composition (Stacked Barplot)
    # --------------------------------------------------------------------------
    strain_cols = [c for c in df_matrix.columns if c.startswith("UP")]
    if strain_cols:
        strain_cat_data = defaultdict(lambda: defaultdict(int))
        for _, row in df_matrix.iterrows():
            cat = row["category"]
            for s in strain_cols:
                if row[s] > 0:
                    strain_cat_data[s][cat] += row[s]

        df_strain_comp = pd.DataFrame(strain_cat_data).T.reindex(columns=valid_cats).fillna(0)

        fig, ax = plt.subplots(figsize=(11.5, 5.5))
        bottom = np.zeros(len(df_strain_comp))

        for cat in valid_cats:
            values = df_strain_comp[cat].values
            ax.bar(
                df_strain_comp.index,
                values,
                bottom=bottom,
                label=cat,
                color=PUB_PALETTE[cat],
                edgecolor="#ffffff",
                linewidth=0.8,
                width=0.55,
            )
            bottom += values

        ax.set_ylabel("Number of Proteins")
        ax.set_xlabel("Strain / Proteome Identifier")
        ax.set_title("Strain Proteome Architecture by Pangenome Category", loc="left")
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            title="Category",
            frameon=True,
            alignment="left",
        )
        ax.grid(axis="y", alpha=0.5)
        plt.xticks(rotation=20, ha="right", fontweight="medium")
        plt.tight_layout()
        export_figure(fig, fig_dir / "step3_strain_proteome_composition")

    # --------------------------------------------------------------------------
    # Figure 3D: Pangenome Matrix Heatmap (Presence/Absence)
    # --------------------------------------------------------------------------
    if strain_cols and len(df_matrix) > 0:
        # Sample or subset clusters for clear visualization
        plot_df = df_matrix.copy()
        # Sort by category (Core -> Shell -> Cloud) then by strain_count descending
        plot_df["cat_rank"] = plot_df["category"].map({"Core": 0, "Shell": 1, "Cloud": 2})
        plot_df = plot_df.sort_values(by=["cat_rank", "strain_count"], ascending=[True, False])

        # If too many clusters, downsample evenly for the visual heatmap
        max_heatmap_clusters = 500
        if len(plot_df) > max_heatmap_clusters:
            sample_step = len(plot_df) // max_heatmap_clusters
            sampled_df = plot_df.iloc[::sample_step]
        else:
            sampled_df = plot_df

        heatmap_data = (sampled_df[strain_cols] > 0).astype(int).T

        fig, ax = plt.subplots(figsize=(12, 5.2))
        sns.heatmap(
            heatmap_data,
            cmap=["#f0f0f0", "#1b7837"],
            cbar_kws={"label": "Gene Presence (1) / Absence (0)", "ticks": [0, 1]},
            linewidths=0.0,
            ax=ax,
        )
        ax.tick_params(axis="y", rotation=0, pad=8)
        ax.set_title(
            f"Pangenome Binary Presence/Absence Matrix ({len(sampled_df):,} Sampled Clusters)",
            loc="left",
        )
        ax.set_xlabel("Gene Clusters (Ordered by Core -> Shell -> Cloud)")
        ax.set_ylabel("Strain Proteome", labelpad=16, fontweight="bold")
        plt.tight_layout()
        export_figure(fig, fig_dir / "step3_pangenome_heatmap")


# ==============================================================================
# Step 4 Visualizations: Subtractive Genomics Filter
# ==============================================================================
def process_step4_visuals(step4_dir: Path):
    """Generates publication-quality charts for Step 4 Subtractive Genomics."""
    master_tsv = step4_dir / "tables" / "subtractive_master.tsv"
    summary_json = step4_dir / "summary" / "subtractive_summary.json"
    fig_dir = step4_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not (master_tsv.exists() and summary_json.exists()):
        logger.info("Step 4 outputs not found; skipping Step 4 figures.")
        return

    df_master = pd.read_csv(master_tsv, sep="\t")
    with open(summary_json, "r") as f:
        summary_data = json.load(f)

    metrics = summary_data["filtering_metrics"]
    tier_counts = metrics["tier_counts"]
    total_core = summary_data["total_core_proteins"]

    # --------------------------------------------------------------------------
    # Figure 4A: Subtractive Genomics Attrition Waterfall / Funnel
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 5.2))
    
    stages = [
        "Core Pangenome",
        "Human Non-Homologous\n(Safe)",
        "DEG Essential\n(Indispensable)",
        "Tier 1: Pathogen-Selective\n(Microbiome-Safe)",
    ]
    t1_count = tier_counts.get("Tier 1: Pathogen-Selective (Microbiome-Safe)", tier_counts.get("Tier 1: Golden Target", 0))
    counts = [
        total_core,
        metrics["human_non_homologous_count"],
        metrics["deg_essential_count"],
        t1_count,
    ]
    colors = ["#2b5c8f", "#4192d9", "#f39c12", "#27ae60"]

    bars = ax.bar(stages, counts, color=colors, edgecolor="#333333", linewidth=1.0, width=0.55)
    for bar, count in zip(bars, counts):
        h = bar.get_height()
        pct = (count / total_core) * 100
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + total_core * 0.015,
            f"{count:,}\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=9.5,
        )

    ax.set_ylim(0, total_core * 1.22)
    ax.set_ylabel("Number of Candidate Proteins")
    ax.set_title("Subtractive Genomics Funnel: Core Target Attrition", loc="left")
    ax.grid(axis="y", alpha=0.5)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step4_subtractive_genomics_funnel")

    # --------------------------------------------------------------------------
    # Figure 4B: Target Prioritization Tier Distribution (Donut Chart)
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5.2))
    
    tier_palette = {
        "Tier 1: Pathogen-Selective (Microbiome-Safe)": "#27ae60",
        "Tier 2: Broad-Spectrum (Shared with Gut)": "#2980b9",
        "Tier 3: Non-Essential (Safe)": "#f39c12",
        "Excluded: Host Cross-Reactive (Cytotoxic)": "#c0392b",
        "Tier 1: Golden Target": "#27ae60",
        "Tier 2: Broad-Spectrum Essential": "#2980b9",
        "Tier 3: Non-Essential Safe": "#f39c12",
        "Excluded (Host Homolog)": "#c0392b",
    }
    
    labels = []
    values = []
    wedge_colors = []
    for k, v in tier_counts.items():
        if v > 0:
            labels.append(f"{k}: {v:,} ({v / total_core * 100:.1f}%)")
            values.append(v)
            wedge_colors.append(tier_palette.get(k, "#7f7f7f"))

    wedges, _ = ax.pie(
        values,
        startangle=140,
        colors=wedge_colors,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5),
    )
    ax.text(
        0, 0,
        f"Total Core\n{total_core:,}",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="#333333",
    )
    ax.legend(
        wedges,
        labels,
        title="Target Classification",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
        alignment="left",
    )
    ax.set_title("Pangenome Target Prioritization by Safety & Essentiality", loc="left")
    plt.tight_layout()
    export_figure(fig, fig_dir / "step4_target_tier_distribution")


# ==============================================================================
# Step 5 Visualizations: Subcellular Localization & TMH Druggability
# ==============================================================================
def process_step5_visuals(step5_dir: Path):
    """Generates publication-quality charts for Step 5: Subcellular Localization & TM Topology."""
    master_tsv = step5_dir / "tables" / "subcellular_master.tsv"
    fig_dir = step5_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not master_tsv.exists():
        logger.info("Step 5 master TSV not found; skipping Step 5 visuals.")
        return

    df = pd.read_csv(master_tsv, sep="\t")
    total_targets = len(df)

    # --------------------------------------------------------------------------
    # Figure 5A: Subcellular Compartment Distribution by Target Tier
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    # Order compartments logically from inside out
    compartment_order = ["Cytoplasmic", "Cytoplasmic Membrane", "Periplasmic", "Outer Membrane", "Extracellular"]
    tier_order = [
        "Tier 1: Pathogen-Selective (Microbiome-Safe)",
        "Tier 2: Broad-Spectrum (Shared with Gut)",
        "Tier 3: Non-Essential (Safe)",
        "Excluded: Host Cross-Reactive (Cytotoxic)"
    ]
    tier_palette = {
        "Tier 1: Pathogen-Selective (Microbiome-Safe)": "#27ae60",
        "Tier 2: Broad-Spectrum (Shared with Gut)": "#2980b9",
        "Tier 3: Non-Essential (Safe)": "#f39c12",
        "Excluded: Host Cross-Reactive (Cytotoxic)": "#c0392b"
    }

    # Cross-tabulate
    ct = pd.crosstab(df["Predicted_Compartment"], df["Target_Tier"])
    ct = ct.reindex(index=[c for c in compartment_order if c in ct.index], fill_value=0)
    valid_tiers = [t for t in tier_order if t in ct.columns]
    ct = ct[valid_tiers]

    if len(valid_tiers) == 1:
        # Focused Single Tier Mode (e.g. Tier 1 only)
        tier_name = valid_tiers[0]
        comp_palette = {
            "Cytoplasmic": "#27ae60",
            "Cytoplasmic Membrane": "#2980b9",
            "Periplasmic": "#f39c12",
            "Outer Membrane": "#8e44ad",
            "Extracellular": "#16a085"
        }
        bar_colors = [comp_palette.get(c, "#7f7f7f") for c in ct.index]
        bars = ax.bar(
            range(len(ct)),
            ct[tier_name],
            color=bar_colors,
            edgecolor="#333333",
            linewidth=0.9,
            width=0.55
        )
        totals = ct[tier_name]
        for idx, total in enumerate(totals):
            if total > 0:
                pct = (total / total_targets) * 100
                ax.text(
                    idx,
                    total + total_targets * 0.015,
                    f"{total:,}\n({pct:.1f}%)",
                    ha="center",
                    va="bottom",
                    fontweight="bold",
                    fontsize=9.5
                )
        ax.set_xticks(range(len(ct)))
        ax.set_xticklabels(ct.index, rotation=15, ha="right", fontweight="medium")
        ax.set_title(f"Subcellular Localization: {tier_name}", loc="left")
    else:
        # Multi-tier stacked barplot
        ct.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=[tier_palette.get(t, "#7f7f7f") for t in valid_tiers],
            edgecolor="#333333",
            linewidth=0.9,
            width=0.55
        )
        totals = ct.sum(axis=1)
        for idx, (comp, total) in enumerate(totals.items()):
            if total > 0:
                pct = (total / total_targets) * 100
                ax.text(
                    idx,
                    total + total_targets * 0.015,
                    f"{total:,}\n({pct:.1f}%)",
                    ha="center",
                    va="bottom",
                    fontweight="bold",
                    fontsize=9.5
                )
        ax.set_title("Subcellular Localization Distribution by Target Safety Tier", loc="left")
        ax.legend(title="Target Prioritization Tier", loc="upper right", frameon=True, fontsize=8.5)

    ax.set_ylim(0, totals.max() * 1.22)
    ax.set_xlabel("Bacterial Subcellular Compartment", fontweight="medium")
    ax.set_ylabel("Number of Target Candidates", fontweight="medium")
    plt.xticks(rotation=15, ha="right", fontweight="medium")
    ax.grid(axis="y", alpha=0.5)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step5_subcellular_compartment_distribution")

    # --------------------------------------------------------------------------
    # Figure 5B: Transmembrane Topology & Druggability Feasibility (Dual Panel)
    # --------------------------------------------------------------------------
    fig, (ax_donut, ax_hist) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1.2]})

    # Panel 1: Feasible vs Deprioritized Donut
    feas_counts = df["Druggability_Feasibility"].value_counts()
    feas_colors = {
        "Feasible (TMH <= 1)": "#27ae60",
        "Deprioritized (TMH >= 2)": "#e74c3c"
    }
    wedge_colors = [feas_colors.get(k, "#95a5a6") for k in feas_counts.index]
    wedges, _ = ax_donut.pie(
        feas_counts.values,
        startangle=140,
        colors=wedge_colors,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5)
    )
    center_label = f"Tier 1\n{total_targets:,}" if len(valid_tiers) == 1 else f"Total\n{total_targets:,}"
    ax_donut.text(
        0, 0,
        center_label,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="#333333"
    )
    legend_labels = [f"{k}: {v:,} ({v / total_targets * 100:.1f}%)" for k, v in feas_counts.items()]
    ax_donut.legend(
        wedges,
        legend_labels,
        title="Assay Feasibility",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        frameon=True,
        fontsize=9
    )
    ax_donut.set_title("A. Druggability Screening (≤ 1 TMH)", loc="left")

    # Panel 2: TMH Count Frequency Distribution
    tmh_counts = df["TMH_Count"].value_counts().sort_index()
    # Group >= 8 TMHs
    bins_data = {}
    for tmh, cnt in tmh_counts.items():
        key = str(tmh) if tmh < 8 else "≥8"
        bins_data[key] = bins_data.get(key, 0) + cnt

    bin_keys = [str(i) for i in range(8)] + ["≥8"]
    bin_values = [bins_data.get(k, 0) for k in bin_keys]
    bar_colors = ["#27ae60" if k in ["0", "1"] else "#e74c3c" for k in bin_keys]

    bars = ax_hist.bar(bin_keys, bin_values, color=bar_colors, edgecolor="#333333", linewidth=0.9, width=0.6)
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax_hist.text(
                bar.get_x() + bar.get_width() / 2,
                h + max(bin_values) * 0.015,
                f"{h:,}",
                ha="center",
                va="bottom",
                fontweight="bold",
                fontsize=8.5
            )

    ax_hist.set_ylim(0, max(bin_values) * 1.18)
    ax_hist.set_xlabel("Predicted Transmembrane Helices (TMH Count)", fontweight="medium")
    ax_hist.set_ylabel("Number of Candidate Proteins", fontweight="medium")
    ax_hist.set_title("B. Transmembrane Helix Distribution Across Pangenome", loc="left")
    ax_hist.grid(axis="y", alpha=0.5)

    # Custom legend for panel B
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#27ae60", edgecolor="#333", label="Soluble / Single-Pass (TMH ≤ 1)"),
        Patch(facecolor="#e74c3c", edgecolor="#333", label="Multi-Pass Integral Membrane (TMH ≥ 2)")
    ]
    ax_hist.legend(handles=legend_elements, loc="upper right", frameon=True, fontsize=8.5)

    plt.tight_layout()
    export_figure(fig, fig_dir / "step5_tmh_druggability_distribution")


# ==============================================================================
# Step 6: Virulence Factor (VFDB) & Antimicrobial Resistance (CARD) Visuals
# ==============================================================================
def process_step6_visuals(step6_dir: Path):
    """Generates publication-quality charts for Step 6 (Virulence & Resistance)."""
    master_tsv = step6_dir / "tables" / "virulence_resistance_master.tsv"
    fig_dir = step6_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not master_tsv.exists():
        logger.warning("Step 6 master table not found at: %s. Skipping Step 6 figures.", master_tsv)
        return

    logger.info("Generating Step 6 publication figures from %s...", master_tsv.name)
    df = pd.read_csv(master_tsv, sep="\t", low_memory=False)
    total_targets = len(df)

    # --------------------------------------------------------------------------
    # Figure 6A: Therapeutic Profile Stratification Donut & Breakdown
    # --------------------------------------------------------------------------
    profile_counts = df["Therapeutic_Profile"].value_counts()
    profile_palette = {
        "Dual-Action Super-Target (Lethal + Anti-Virulence)": "#27ae60",   # Emerald Green
        "Clean Essential Drug Target (Low Resistance Risk)": "#2980b9",   # Royal Blue
        "Resistance-Breaker / Adjuvant Candidate": "#d35400",             # Rust Amber
        "Anti-Virulence Candidate (Disarming Pathogen)": "#8e44ad",       # Amethyst Purple
        "Secondary / Accessory Target": "#7f8c8d"                         # Gray
    }
    
    fig, ax = plt.subplots(figsize=(9, 6.2))
    colors_prof = [profile_palette.get(k, "#95a5a6") for k in profile_counts.index]
    wedges, _ = ax.pie(
        profile_counts.values,
        startangle=140,
        colors=colors_prof,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5)
    )
    ax.text(0, 0, f"Prioritized\nTargets\n{total_targets:,}", ha="center", va="center", fontsize=11, fontweight="bold", color="#333333")
    legend_labels = [f"{k}: {v:,} ({v / total_targets * 100:.1f}%)" for k, v in profile_counts.items()]
    ax.legend(
        wedges,
        legend_labels,
        title="Target Therapeutic Profile",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        frameon=True,
        fontsize=9
    )
    ax.set_title("Pharmacological Stratification of Prioritized Targets (VFDB & CARD)", loc="left", pad=12)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step6_target_therapeutic_profiles")

    # --------------------------------------------------------------------------
    # Figure 6B: Virulence Factor Functional Classes (if any hits)
    # --------------------------------------------------------------------------
    vf_df = df[df["VFDB_Hit"] == True]
    if len(vf_df) > 0 and "VFDB_Class" in vf_df.columns:
        vf_classes = vf_df["VFDB_Class"].value_counts().head(10)
        fig, ax = plt.subplots(figsize=(10, 5.5))
        bars = ax.barh(vf_classes.index[::-1], vf_classes.values[::-1], color="#27ae60", edgecolor="#333333", linewidth=0.9, height=0.65)
        for bar in bars:
            w = bar.get_width()
            ax.text(w + max(vf_classes.values) * 0.015, bar.get_y() + bar.get_height() / 2, f"{w:,}", ha="left", va="center", fontweight="bold", fontsize=9)
        ax.set_xlim(0, max(vf_classes.values) * 1.15)
        ax.set_xlabel("Number of Targets", fontweight="medium")
        ax.set_ylabel("Virulence Factor Classification", fontweight="medium")
        ax.set_title(f"Virulence Factor Diversity Across Prioritized Targets (n={len(vf_df)})", loc="left")
        ax.grid(axis="x", alpha=0.5)
        plt.tight_layout()
        export_figure(fig, fig_dir / "step6_virulence_classes_distribution")

    # --------------------------------------------------------------------------
    # Figure 6C: CARD Resistance Mechanisms (if any hits)
    # --------------------------------------------------------------------------
    card_df = df[df["CARD_Hit"] == True]
    if len(card_df) > 0 and "CARD_Resistance_Mechanism" in card_df.columns:
        mech_counts = card_df["CARD_Resistance_Mechanism"].value_counts().head(8)
        fig, ax = plt.subplots(figsize=(10, 5.2))
        bars = ax.barh(mech_counts.index[::-1], mech_counts.values[::-1], color="#d35400", edgecolor="#333333", linewidth=0.9, height=0.6)
        for bar in bars:
            w = bar.get_width()
            ax.text(w + max(mech_counts.values) * 0.02, bar.get_y() + bar.get_height() / 2, f"{w:,}", ha="left", va="center", fontweight="bold", fontsize=9)
        ax.set_xlim(0, max(mech_counts.values) * 1.18)
        ax.set_xlabel("Number of Targets", fontweight="medium")
        ax.set_ylabel("Molecular Resistance Mechanism", fontweight="medium")
        ax.set_title(f"Antimicrobial Resistance Mechanisms Identified via CARD (n={len(card_df)})", loc="left")
        ax.grid(axis="x", alpha=0.5)
        plt.tight_layout()
        export_figure(fig, fig_dir / "step6_card_resistance_mechanisms")


# ==============================================================================
# Step 7: Target Druggability & Novelty Assessment (TTD 2026) Visuals
# ==============================================================================
def process_step7_visuals(step7_dir: Path):
    """Generates publication-quality charts for Step 7 (TTD Druggability & Novelty)."""
    master_tsv = step7_dir / "tables" / "ttd_druggability_master.tsv"
    fig_dir = step7_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not master_tsv.exists():
        logger.warning("Step 7 master table not found at: %s. Skipping Step 7 figures.", master_tsv)
        return

    logger.info("Generating Step 7 publication figures from %s...", master_tsv.name)
    df = pd.read_csv(master_tsv, sep="\t", low_memory=False)
    total_targets = len(df)

    # --------------------------------------------------------------------------
    # Figure 7A: Clinical Development Status & Target Novelty Donut Chart
    # --------------------------------------------------------------------------
    status_counts = df["Clinical_Druggability_Category"].value_counts()
    status_palette = {
        "Approved Drug Target (Repurposing Opportunity)": "#27ae60",               # Emerald Green
        "Clinical Trial Target (Phase I-III Validated)": "#2980b9",                # Royal Blue
        "Patented / Preclinical Target": "#f39c12",                                # Amber
        "Research Target (Literature Precedent)": "#8e44ad",                       # Purple
        "Novel Dual-Action Super-Target (Pristine First-in-Class)": "#16a085",     # Dark Teal
        "Novel Clean Essential Target (First-in-Class Metabolic Enzyme)": "#2c3e50", # Midnight Navy
        "Novel Resistance-Breaker / Adjuvant Candidate": "#d35400",                # Rust Orange
        "Novel Antibacterial Target": "#7f8c8d"                                    # Gray
    }

    fig, ax = plt.subplots(figsize=(10, 6.5))
    colors_stat = [status_palette.get(k, "#95a5a6") for k in status_counts.index]
    wedges, _ = ax.pie(
        status_counts.values,
        startangle=140,
        colors=colors_stat,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5)
    )
    novel_count = int((df["TTD_Hit"] == False).sum())
    novel_pct = (novel_count / total_targets) * 100
    ax.text(
        0, 0,
        f"Total: {total_targets:,}\nNovel: {novel_pct:.1f}%",
        ha="center", va="center",
        fontsize=11, fontweight="bold", color="#333333"
    )
    legend_labels = [f"{k}: {v:,} ({v / total_targets * 100:.1f}%)" for k, v in status_counts.items()]
    ax.legend(
        wedges,
        legend_labels,
        title="Clinical Status & Novelty (TTD 2026)",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        frameon=True,
        fontsize=8.5
    )
    ax.set_title("Clinical Druggability & Novelty Assessment of Prioritized Targets (TTD 2026)", loc="left", pad=12)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step7_clinical_druggability_novelty")

    # --------------------------------------------------------------------------
    # Figure 7B: Stratification of Completely Novel / First-in-Class Targets
    # --------------------------------------------------------------------------
    novel_df = df[df["TTD_Hit"] == False]
    if len(novel_df) > 0:
        novel_counts = novel_df["Clinical_Druggability_Category"].value_counts()
        fig, ax = plt.subplots(figsize=(9, 4.8))
        colors_bar = [status_palette.get(k, "#16a085") for k in novel_counts.index]
        bars = ax.barh(novel_counts.index[::-1], novel_counts.values[::-1], color=colors_bar[::-1], edgecolor="#333333", linewidth=0.9, height=0.6)
        for bar in bars:
            w = bar.get_width()
            ax.text(w + max(novel_counts.values) * 0.015, bar.get_y() + bar.get_height() / 2, f"{w:,} ({w/len(novel_df)*100:.1f}%)", ha="left", va="center", fontweight="bold", fontsize=9)
        ax.set_xlim(0, max(novel_counts.values) * 1.25)
        ax.set_xlabel("Number of Novel Targets", fontweight="medium")
        ax.set_ylabel("First-in-Class Candidate Subcategory", fontweight="medium")
        ax.set_title(f"First-in-Class Novel Antibacterial Target Landscape (n={len(novel_df)})", loc="left")
        ax.grid(axis="x", alpha=0.5)
        plt.tight_layout()
        export_figure(fig, fig_dir / "step7_novel_targets_stratification")


# ==============================================================================
# Step 8 Visualizations: Metabolic Pathway & Choke-Point Analysis
# ==============================================================================
def process_step8_visuals(step8_dir: Path):
    """Generates publication-quality charts for Step 8 (Metabolic Choke-Points & Pathways)."""
    master_tsv = step8_dir / "tables" / "metabolic_master.tsv"
    fig_dir = step8_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not master_tsv.exists():
        logger.warning("Step 8 master table not found at: %s. Skipping Step 8 figures.", master_tsv)
        return

    logger.info("Generating Step 8 publication figures from %s...", master_tsv.name)
    df = pd.read_csv(master_tsv, sep="\t", low_memory=False)
    total_targets = len(df)

    # --------------------------------------------------------------------------
    # Figure 8A: Choke-Point Vulnerability Status Donut Chart
    # --------------------------------------------------------------------------
    choke_counts = df["Choke_Point_Status"].value_counts()
    choke_palette = {
        "Double Choke-Point (Bactericidal + Starvation)": "#c0392b",
        "Single Choke-Point (Unique Substrate Consumer)": "#d35400",
        "Single Choke-Point (Unique Product Synthesizer)": "#e67e22",
        "Non-Choke-Point (Metabolic Bypass Available)": "#2980b9",
        "Non-Enzymatic / Structural Target": "#7f8c8d"
    }

    fig, ax = plt.subplots(figsize=(10, 6.5))
    colors_choke = [choke_palette.get(k, "#95a5a6") for k in choke_counts.index]
    wedges, _ = ax.pie(
        choke_counts.values,
        startangle=140,
        colors=colors_choke,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5)
    )
    total_choke = int(df["Is_Choke_Point"].sum())
    double_choke = int(df["Is_Double_Choke_Point"].sum())
    ax.text(
        0, 0,
        f"Bottlenecks\n{total_choke} / {total_targets}\n({total_choke/total_targets*100:.1f}%)\nDouble: {double_choke}",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color="#333333"
    )
    labels = [f"{k}: {v:,} ({v / total_targets * 100:.1f}%)" for k, v in choke_counts.items()]
    ax.legend(
        wedges,
        labels,
        title="Metabolic Choke-Point Classification",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
        alignment="left"
    )
    ax.set_title("Metabolic Network Choke-Point Vulnerability Profile", loc="left")
    plt.tight_layout()
    export_figure(fig, fig_dir / "step8_chokepoint_vulnerability_distribution")

    # --------------------------------------------------------------------------
    # Figure 8B: Definitive Target Priority & Pathogen-Exclusive Pathways
    # --------------------------------------------------------------------------
    prio_counts = df["Definitive_Target_Priority"].value_counts()
    prio_palette = {
        "Priority 1: Gold-Standard Dual-Action Double Choke-Point": "#c0392b",
        "Priority 2: Pathogen-Exclusive Double Choke-Point": "#e74c3c",
        "Priority 3: Core Metabolic Double Choke-Point": "#e67e22",
        "Priority 4: Pathogen-Exclusive Single Choke-Point": "#f39c12",
        "Priority 5: Core Metabolic Single Choke-Point": "#3498db",
        "Clinically Approved Repurposing Benchmark": "#27ae60",
        "Deprioritized / Non-Choke-Point Enzyme": "#95a5a6",
    }

    fig, ax = plt.subplots(figsize=(10, 5.2))
    colors_p = [prio_palette.get(k, "#7f8c8d") for k in prio_counts.index]
    bars = ax.barh(prio_counts.index[::-1], prio_counts.values[::-1], color=colors_p[::-1], edgecolor="#333333", linewidth=0.9, height=0.6)
    for bar in bars:
        w = bar.get_width()
        ax.text(w + max(prio_counts.values) * 0.015, bar.get_y() + bar.get_height() / 2, f"{w:,} ({w/total_targets*100:.1f}%)", ha="left", va="center", fontweight="bold", fontsize=9)
    ax.set_xlim(0, max(prio_counts.values) * 1.25)
    ax.set_xlabel("Number of Candidate Targets", fontweight="medium")
    ax.set_ylabel("Definitive Metabolic Priority Tier", fontweight="medium")
    ax.set_title(f"Definitive Lead Target Prioritization (Total n={total_targets})", loc="left")
    ax.grid(axis="x", alpha=0.5)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step8_definitive_lead_prioritization")


# ==============================================================================
# Step 9 Visualizations: PPI Interactome & Network Centrality Analysis
# ==============================================================================
def process_step9_visuals(step9_dir: Path):
    """Generates publication-quality charts for Step 9 (Network Hub-Bottleneck & Centrality)."""
    master_tsv = step9_dir / "tables" / "network_centrality_master.tsv"
    summary_json = step9_dir / "summary" / "step9_network_summary.json"
    fig_dir = step9_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not master_tsv.exists():
        logger.warning("Step 9 master table not found at: %s. Skipping Step 9 figures.", master_tsv)
        return

    logger.info("Generating Step 9 publication figures from %s...", master_tsv.name)
    df = pd.read_csv(master_tsv, sep="\t", low_memory=False)
    total_targets = len(df)

    # Load summary metrics for exact cutoffs
    deg_cutoff = 42.65
    bet_cutoff = 0.004638
    total_nodes = 4126
    total_edges = 33685
    if summary_json.exists():
        try:
            with open(summary_json, "r") as f:
                sdata = json.load(f)
                inter_stats = sdata.get("proteome_interactome_statistics", {})
                deg_cutoff = inter_stats.get("degree_cutoff_hub", deg_cutoff)
                bet_cutoff = inter_stats.get("betweenness_cutoff_bottleneck", bet_cutoff)
                total_nodes = inter_stats.get("total_proteins_in_graph", total_nodes)
                total_edges = inter_stats.get("total_high_confidence_edges", total_edges)
        except Exception as e:
            logger.warning("Error reading Step 9 summary JSON: %s", e)

    # --------------------------------------------------------------------------
    # Figure 9A: Degree vs. Betweenness Centrality (Hub-Bottleneck Quadrants)
    # --------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10.5, 6.8))

    df_mapped = df[df["STRING_Mapped"] == True].copy()

    topo_palette = {
        "Hub-Bottleneck (HB)": "#c0392b",
        "Non-Hub Bottleneck (NH-B)": "#e67e22",
        "Hub Non-Bottleneck (H-NB)": "#2980b9",
        "Non-Hub Non-Bottleneck (NH-NB)": "#7f8c8d"
    }

    for tclass, color in topo_palette.items():
        sub = df_mapped[df_mapped["Topological_Class"] == tclass]
        if len(sub) > 0:
            ax.scatter(
                sub["STRING_Degree"],
                sub["Betweenness_Centrality"],
                c=color,
                label=f"{tclass} (n={len(sub)})",
                alpha=0.65,
                s=45,
                edgecolors="none"
            )

    ax.axvline(x=deg_cutoff, color="#34495e", linestyle="--", linewidth=1.2, alpha=0.85, label=f"Hub Cutoff (k ≥ {deg_cutoff:.1f})")
    ax.axhline(y=bet_cutoff, color="#34495e", linestyle=":", linewidth=1.2, alpha=0.85, label=f"Bottleneck Cutoff (BC ≥ {bet_cutoff:.5f})")

    highlight_genes = {
        "entB": ("#c0392b", "entB (Priority 1A Ultimate Lead)"),
        "acpP": ("#8e44ad", "acpP (Fatty Acid Hub)"),
        "rpoA": ("#2c3e50", "rpoA (RNA Pol Hub)"),
        "glmM": ("#d35400", "glmM (Cell Wall Hub)"),
        "ugd": ("#16a085", "ugd (Capsule/LPS Lead)"),
        "lpxD": ("#27ae60", "lpxD (Lipid A Lead)"),
        "murF": ("#2980b9", "murF (Peptidoglycan Lead)")
    }

    max_bet = max(df_mapped["Betweenness_Centrality"].max() if len(df_mapped) > 0 else 0.01, bet_cutoff * 1.5)

    for gene, (color, annot_text) in highlight_genes.items():
        sub = df[df["Gene_Name"] == gene]
        if len(sub) > 0:
            r = sub.iloc[0]
            if r["STRING_Mapped"]:
                x_val = r["STRING_Degree"]
                y_val = r["Betweenness_Centrality"]
                ax.scatter(x_val, y_val, c=color, s=120, edgecolors="#ffffff", linewidth=1.5, zorder=5)
                y_offset = (max_bet * 0.035) if y_val < (max_bet * 0.5) else (-max_bet * 0.035)
                ax.annotate(
                    annot_text,
                    xy=(x_val, y_val),
                    xytext=(x_val + 3, y_val + y_offset),
                    fontsize=8.5,
                    fontweight="bold",
                    color=color,
                    arrowprops=dict(arrowstyle="->", color=color, lw=0.8),
                    zorder=6
                )

    ax.set_xlabel("Degree Centrality (k, Interacting Partners in High-Confidence PPI)", fontweight="medium")
    ax.set_ylabel("Betweenness Centrality (Shortest Path Fraction)", fontweight="medium")
    ax.set_title("Interactome Network Centrality & Hub-Bottleneck Quadrant Stratification", loc="left")
    ax.legend(title="Topological Classifications", loc="upper left", frameon=True, fontsize=8.5)
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    export_figure(fig, fig_dir / "step9_degree_vs_betweenness_centrality")

    # --------------------------------------------------------------------------
    # Figure 9B: Network Topological Class Distribution Donut Chart
    # --------------------------------------------------------------------------
    topo_counts = df["Topological_Class"].value_counts()
    fig, ax = plt.subplots(figsize=(9.5, 6.0))

    class_colors = {
        "Hub-Bottleneck (HB)": "#c0392b",
        "Non-Hub Bottleneck (NH-B)": "#e67e22",
        "Hub Non-Bottleneck (H-NB)": "#2980b9",
        "Non-Hub Non-Bottleneck (NH-NB)": "#95a5a6",
        "Unmapped / Peripheral": "#bdc3c7"
    }

    colors = [class_colors.get(k, "#7f8c8d") for k in topo_counts.index]
    wedges, _ = ax.pie(
        topo_counts.values,
        startangle=140,
        colors=colors,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5)
    )

    hb_count = int((df["Topological_Class"] == "Hub-Bottleneck (HB)").sum())

    ax.text(
        0, 0,
        f"Interactome\n{total_nodes:,} Nodes\n{total_edges:,} Edges\n\nHub-Bottlenecks:\n{hb_count} ({hb_count/total_targets*100:.1f}%)",
        ha="center",
        va="center",
        fontsize=10.5,
        fontweight="bold",
        color="#333333"
    )

    labels = [f"{k}: {v:,} ({v / total_targets * 100:.1f}%)" for k, v in topo_counts.items()]
    ax.legend(
        wedges,
        labels,
        title="Topological Classification",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
        alignment="left"
    )
    ax.set_title("Target Stratification by Network Centrality & Lethality Profile", loc="left")
    plt.tight_layout()
    export_figure(fig, fig_dir / "step9_topological_class_distribution")


# ==============================================================================
# Master Overview Figure: Complete Steps 1-3 Summary Panel
# ==============================================================================
def process_master_figure(job_dir: Path):
    """Creates a unified 4-panel publication figure summarizing Steps 1, 2, and 3."""
    step1_meta = job_dir / "01_data_acquisition" / "proteomes_metadata.tsv"
    step2_stats = job_dir / "02_paralog_removal" / "summary" / "paralog_removal_stats.tsv"
    step3_matrix = job_dir / "03_pangenome_analysis" / "matrix" / "pangenome_matrix.tsv"
    master_fig_dir = job_dir / "figures"
    master_fig_dir.mkdir(parents=True, exist_ok=True)

    if not (step1_meta.exists() and step2_stats.exists() and step3_matrix.exists()):
        logger.info("Not all steps complete yet; skipping master multi-panel figure.")
        return

    df_meta = pd.read_csv(step1_meta, sep="\t")
    df_stats = pd.read_csv(step2_stats, sep="\t")
    df_matrix = pd.read_csv(step3_matrix, sep="\t")

    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32)

    # --------------------------------------------------------------------------
    # Panel A: Step 1 CPD & Taxon Acquisition
    # --------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    cpd_counts = df_meta["CPD"].value_counts()
    total_meta = len(df_meta)
    cpd_palette = {
        "Unknown": "#7f7f7f",
        "Standard": "#2ca02c",
        "Close to standard (low value)": "#66c2a4",
        "Close to standard (high value)": "#b2e2e2",
        "Outlier (high value)": "#ff7f0e",
        "Outlier (low value)": "#d62728",
    }
    colors_a = [cpd_palette.get(k, "#cccccc") for k in cpd_counts.index]
    wedges, _ = ax_a.pie(
        cpd_counts.values,
        startangle=140,
        colors=colors_a,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5),
    )
    ax_a.text(0, 0, f"Total\n{total_meta:,}", ha="center", va="center", fontsize=10.5, fontweight="bold", color="#333333")
    legend_a = [f"{k}: {v:,} ({v / total_meta * 100:.1f}%)" for k, v in cpd_counts.items()]
    ax_a.legend(
        wedges,
        legend_a,
        title="CPD Completeness",
        title_fontsize=9.5,
        loc="upper center",
        bbox_to_anchor=(0.0, -1.10),
        bbox_transform=ax_a.transData,
        ncol=2,
        frameon=True,
        fontsize=8.5,
    )
    ax_a.set_title("A. Step 1: UniProt Proteome Completeness (CPD)", loc="left")

    # --------------------------------------------------------------------------
    # Panel B: Step 2 Intra-strain Paralog Removal
    # --------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    x = np.arange(len(df_stats))
    width = 0.36
    ax_b.bar(
        x - width / 2,
        df_stats["Original Proteins"],
        width,
        label="Original Proteins",
        color=PUB_PALETTE["Original"],
        edgecolor="#333333",
        linewidth=0.9,
    )
    ax_b.bar(
        x + width / 2,
        df_stats["Non-Paralogous Proteins"],
        width,
        label="Non-Paralogous",
        color=PUB_PALETTE["Non-Paralog"],
        edgecolor="#333333",
        linewidth=0.9,
    )
    for idx, row in df_stats.iterrows():
        ax_b.text(
            x[idx],
            row["Original Proteins"] + 30,
            f"-{row['Reduction (%)']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=PUB_PALETTE["Paralogs"],
        )
    ax_b.set_ylabel("Proteins")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(df_stats["Proteome Id"], rotation=20, ha="right", fontweight="medium")
    ax_b.set_ylim(0, 4750)
    ax_b.legend(loc="upper right", frameon=True, fontsize=9)
    ax_b.grid(axis="y", alpha=0.5)
    ax_b.set_title("B. Step 2: Paralog Filtering (MMseqs2)", loc="left")

    # --------------------------------------------------------------------------
    # Panel C: Step 3 Pangenome Partition
    # --------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    cat_order = ["Core", "Shell", "Cloud"]
    cat_clusters = df_matrix["category"].value_counts().reindex(cat_order).fillna(0)
    valid_cats = [c for c in cat_order if cat_clusters[c] > 0]
    p_colors = [PUB_PALETTE[c] for c in valid_cats]

    wedges, texts, autotexts = ax_c.pie(
        cat_clusters[valid_cats],
        labels=[f"{c} ({int(cat_clusters[c]):,} clusters)" for c in valid_cats],
        autopct="%1.1f%%",
        pctdistance=0.72,
        startangle=140,
        colors=p_colors,
        wedgeprops=dict(width=0.45, edgecolor="#ffffff", linewidth=1.5),
    )
    for at in autotexts:
        at.set_fontsize(9)
        at.set_weight("bold")
    ax_c.set_title("C. Step 3: Pangenome Partitioning", loc="left")

    # --------------------------------------------------------------------------
    # Panel D: Step 3 Strain Sharing Distribution
    # --------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    strain_counts = df_matrix["strain_count"].value_counts().sort_index()
    d_colors = [
        PUB_PALETTE.get(df_matrix[df_matrix["strain_count"] == sc]["category"].iloc[0], "#6c757d")
        for sc in strain_counts.index
    ]
    bars = ax_d.bar(
        strain_counts.index,
        strain_counts.values,
        color=d_colors,
        edgecolor="#333333",
        linewidth=1.0,
        width=0.6,
    )
    for bar in bars:
        h = bar.get_height()
        ax_d.text(
            bar.get_x() + bar.get_width() / 2,
            h + max(strain_counts.values) * 0.015,
            f"{int(h):,}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=9,
        )
    ax_d.set_ylim(0, max(strain_counts.values) * 1.18)
    ax_d.set_xlabel("Number of Strains")
    ax_d.set_ylabel("Gene Clusters")
    ax_d.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax_d.grid(axis="y", alpha=0.5)
    ax_d.set_title("D. Step 3: Gene Frequency Across Strains", loc="left")

    fig.suptitle(
        f"Pangenome-Guided Drug Target Discovery Pipeline: Steps 1–3 ({df_meta['Organism'].iloc[0]})",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )

    export_figure(fig, master_fig_dir / "pangenome_pipeline_overview")


def main():
    parser = argparse.ArgumentParser(
        description="Generate publication-ready figures and CSV files for Steps 1-3"
    )
    parser.add_argument(
        "--config", default="config.yaml", help="Path to configuration YAML"
    )
    args = parser.parse_args()

    apply_publication_theme()

    # Load configuration
    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    job_name = cfg.get("job_name", "default_job")
    results_dir = Path(cfg.get("results_dir", "results"))
    job_dir = results_dir / job_name

    logger.info("Generating publication-ready figures & CSVs for %s...", job_name)

    # Process each step
    step1_dir = job_dir / "01_data_acquisition"
    step2_dir = job_dir / "02_paralog_removal"
    step3_dir = job_dir / "03_pangenome_analysis"
    step4_dir = job_dir / "04_subtractive_filter"
    step5_dir = job_dir / "05_subcellular_localization"
    step6_dir = job_dir / "06_virulence_resistance_screening"
    step7_dir = job_dir / "07_target_druggability_ttd"
    step8_dir = job_dir / "08_metabolic_chokepoint"
    step9_dir = job_dir / "09_network_biology"

    process_step1_visuals(step1_dir)
    process_step2_visuals(step2_dir)
    process_step3_visuals(step3_dir)
    process_step4_visuals(step4_dir)
    process_step5_visuals(step5_dir)
    process_step6_visuals(step6_dir)
    process_step7_visuals(step7_dir)
    process_step8_visuals(step8_dir)
    process_step9_visuals(step9_dir)
    process_master_figure(job_dir)

    logger.info("All publication figures and CSVs generated successfully!")


if __name__ == "__main__":
    import yaml
    main()
