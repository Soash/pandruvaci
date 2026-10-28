# ==============================================================================
# Snakefile: Pangenome-Guided Drug Target Identification Pipeline
# ==============================================================================

import os
import sys
import yaml

configfile: "config.yaml"

from scripts.utils_hardware import detect_available_cores, resolve_threads
SYSTEM_CORES = detect_available_cores()

PYTHON = sys.executable

# ------------------------------------------------------------------------------
# Target Pathogen Validation & Parameter Resolution
# ------------------------------------------------------------------------------
JOB = config.get("job_name")
TAX_ID = config.get("taxonomy_id") or config.get("query_params", {}).get("taxonomy_id")
GRAM = config.get("gram_stain") or config.get("step5_subcellular_localization", {}).get("gram_stain")

# If required parameters are missing, display clear instruction guide
if (
    not JOB or str(JOB).lower() in ("null", "none", "") or
    not TAX_ID or str(TAX_ID).lower() in ("null", "none", "") or
    not GRAM or str(GRAM).lower() in ("null", "none", "")
):
    sys.stderr.write(
        "\n"
        "================================================================================\n"
        "  [PANGENOME-DRUG PIPELINE] ERROR: Missing Required Pathogen Parameters!\n"
        "================================================================================\n"
        "The pipeline is organism-agnostic. Please specify your target pathogen via\n"
        "Snakemake's --config command-line flags:\n\n"
        "  snakemake --cores all --config \\\n"
        "      job_name=\"<organism_name>\" \\\n"
        "      taxonomy_id=<ncbi_taxonomy_id> \\\n"
        "      gram_stain=\"negative\" # or \"positive\"\n\n"
        "Examples:\n"
        "  - Acinetobacter baumannii:\n"
        "      snakemake --cores all --config job_name=\"acinetobacter_baumannii\" taxonomy_id=470 gram_stain=\"negative\"\n\n"
        "  - Pseudomonas aeruginosa:\n"
        "      snakemake --cores all --config job_name=\"pseudomonas_aeruginosa\" taxonomy_id=287 gram_stain=\"negative\"\n\n"
        "  - Staphylococcus aureus (Gram-Positive):\n"
        "      snakemake --cores all --config job_name=\"staphylococcus_aureus\" taxonomy_id=1280 gram_stain=\"positive\"\n\n"
        "  - Klebsiella pneumoniae:\n"
        "      snakemake --cores all --config job_name=\"klebsiella_pneumoniae\" taxonomy_id=573 gram_stain=\"negative\"\n"
        "================================================================================\n\n"
    )
    sys.exit(1)

config["job_name"] = str(JOB).strip()
config["taxonomy_id"] = int(TAX_ID)
config["gram_stain"] = str(GRAM).strip().lower()

if config["gram_stain"] not in ("negative", "positive"):
    sys.stderr.write(
        f"\n[PANGENOME-DRUG PIPELINE] ERROR: Invalid gram_stain '{config['gram_stain']}'. "
        "Must be either 'negative' or 'positive'.\n\n"
    )
    sys.exit(1)

if "step5_subcellular_localization" in config and isinstance(config["step5_subcellular_localization"], dict):
    config["step5_subcellular_localization"]["gram_stain"] = config["gram_stain"]

RESULTS_DIR = config.get("results_dir", "results")
JOB_DIR = f"{RESULTS_DIR}/{config['job_name']}"

# Propagate active runtime config (including all CLI --config overrides) to child scripts
os.makedirs(JOB_DIR, exist_ok=True)
ACTIVE_CONFIG = f"{JOB_DIR}/active_config.yaml"
with open(ACTIVE_CONFIG, "w") as _cfg_out:
    yaml.dump(dict(config), _cfg_out, default_flow_style=False)

# Step 1 Paths
STEP1_DIR = f"{JOB_DIR}/01_data_acquisition"
STEP1_METADATA = f"{STEP1_DIR}/proteomes_metadata.tsv"
STEP1_SELECTED = f"{STEP1_DIR}/selected_proteomes.tsv"
STEP1_FASTA_DIR = f"{STEP1_DIR}/fasta"
STEP1_SUMMARY = f"{STEP1_DIR}/acquisition_summary.json"
STEP1_DONE = f"{STEP1_FASTA_DIR}/.download_complete"

# Step 2 Paths
STEP2_DIR = f"{JOB_DIR}/02_paralog_removal"
STEP2_FASTA_DIR = f"{STEP2_DIR}/fasta"
STEP2_STATS = f"{STEP2_DIR}/summary/paralog_removal_stats.tsv"
STEP2_SUMMARY = f"{STEP2_DIR}/summary/paralog_removal_summary.json"
STEP2_DONE = f"{STEP2_FASTA_DIR}/.paralog_removal_complete"

# Step 3 Paths
STEP3_DIR = f"{JOB_DIR}/03_pangenome_analysis"
STEP3_MATRIX = f"{STEP3_DIR}/matrix/pangenome_matrix.tsv"
STEP3_CLUSTERS = f"{STEP3_DIR}/clusters/pangenome_clusters.tsv"
STEP3_CORE_REP = f"{STEP3_DIR}/fasta/core_representatives.fasta"
STEP3_SUMMARY = f"{STEP3_DIR}/summary/pangenome_summary.json"
STEP3_DONE = f"{STEP3_DIR}/.pangenome_complete"

# Step 4 Paths
STEP4_DIR = f"{JOB_DIR}/04_subtractive_filter"
STEP4_MASTER = f"{STEP4_DIR}/tables/subtractive_master.tsv"
STEP4_TIER1_TSV = f"{STEP4_DIR}/tables/tier1_golden_targets.tsv"
STEP4_TIER1_FASTA = f"{STEP4_DIR}/fasta/tier1_golden_targets.fasta"
STEP4_SUMMARY = f"{STEP4_DIR}/summary/subtractive_summary.json"
STEP4_DONE = f"{STEP4_DIR}/.subtractive_complete"

# Step 5 Paths
STEP5_DIR = f"{JOB_DIR}/05_subcellular_localization"
STEP5_MASTER = f"{STEP5_DIR}/tables/subcellular_master.tsv"
STEP5_CYTO_TSV = f"{STEP5_DIR}/tables/tier1_druggable_cytoplasmic.tsv"
STEP5_CYTO_FASTA = f"{STEP5_DIR}/fasta/tier1_druggable_cytoplasmic.fasta"
STEP5_SURFACE_TSV = f"{STEP5_DIR}/tables/tier1_surface_vaccine_targets.tsv"
STEP5_SURFACE_FASTA = f"{STEP5_DIR}/fasta/tier1_surface_vaccine_targets.fasta"
STEP5_PRIO_TSV = f"{STEP5_DIR}/tables/prioritized_drug_candidates.tsv"
STEP5_PRIO_FASTA = f"{STEP5_DIR}/fasta/prioritized_drug_candidates.fasta"
STEP5_SUMMARY = f"{STEP5_DIR}/summary/subcellular_summary.json"
STEP5_DONE = f"{STEP5_DIR}/.subcellular_complete"

# Step 6 Paths
STEP6_DIR = f"{JOB_DIR}/06_virulence_resistance_screening"
STEP6_MASTER = f"{STEP6_DIR}/tables/virulence_resistance_master.tsv"
STEP6_DUAL_TSV = f"{STEP6_DIR}/tables/dual_action_super_targets.tsv"
STEP6_DUAL_FASTA = f"{STEP6_DIR}/fasta/dual_action_super_targets.fasta"
STEP6_CLEAN_TSV = f"{STEP6_DIR}/tables/clean_essential_drug_targets.tsv"
STEP6_CLEAN_FASTA = f"{STEP6_DIR}/fasta/clean_essential_drug_targets.fasta"
STEP6_SUMMARY = f"{STEP6_DIR}/summary/step6_virulence_resistance_summary.json"
STEP6_DONE = f"{STEP6_DIR}/.virulence_resistance_complete"

# Step 7 Paths
STEP7_DIR = f"{JOB_DIR}/07_target_druggability_ttd"
STEP7_MASTER = f"{STEP7_DIR}/tables/ttd_druggability_master.tsv"
STEP7_NOVEL_TSV = f"{STEP7_DIR}/tables/novel_first_in_class_targets.tsv"
STEP7_NOVEL_FASTA = f"{STEP7_DIR}/fasta/novel_first_in_class_targets.fasta"
STEP7_APPROVED_TSV = f"{STEP7_DIR}/tables/approved_drug_targets.tsv"
STEP7_SUMMARY = f"{STEP7_DIR}/summary/step7_ttd_druggability_summary.json"
STEP7_DONE = f"{STEP7_DIR}/.ttd_druggability_complete"

# Step 8 Paths
STEP8_DIR = f"{JOB_DIR}/08_metabolic_chokepoint"
STEP8_MASTER = f"{STEP8_DIR}/tables/metabolic_master.tsv"
STEP8_DOUBLE_TSV = f"{STEP8_DIR}/tables/double_chokepoint_targets.tsv"
STEP8_GOLD_TSV = f"{STEP8_DIR}/tables/definitive_gold_standard_targets.tsv"
STEP8_GOLD_FASTA = f"{STEP8_DIR}/fasta/definitive_gold_standard_targets.fasta"
STEP8_SUMMARY = f"{STEP8_DIR}/summary/step8_metabolic_summary.json"
STEP8_DONE = f"{STEP8_DIR}/.metabolic_complete"

# Step 9 Paths
STEP9_DIR = f"{JOB_DIR}/09_network_biology"
STEP9_MASTER = f"{STEP9_DIR}/tables/network_centrality_master.tsv"
STEP9_HB_TSV = f"{STEP9_DIR}/tables/hub_bottleneck_targets.tsv"
STEP9_LEADS_TSV = f"{STEP9_DIR}/tables/definitive_interactome_leads.tsv"
STEP9_EDGES_TSV = f"{STEP9_DIR}/tables/target_ppi_edges.tsv"
STEP9_HB_FASTA = f"{STEP9_DIR}/fasta/hub_bottleneck_targets.fasta"
STEP9_LEADS_FASTA = f"{STEP9_DIR}/fasta/definitive_interactome_leads.fasta"
STEP9_SUMMARY = f"{STEP9_DIR}/summary/step9_network_summary.json"
STEP9_DONE = f"{STEP9_DIR}/.network_complete"

# Publication Reports
JOB_METHODOLOGY = f"{JOB_DIR}/METHODOLOGY.md"
JOB_RESULTS = f"{JOB_DIR}/RESULTS.md"
MASTER_FIG = f"{JOB_DIR}/figures/pangenome_pipeline_overview.png"

LOG_DIR = f"{JOB_DIR}/logs"

# ------------------------------------------------------------------------------
# Target Rule: All
# ------------------------------------------------------------------------------
rule all:
    input:
        STEP1_METADATA,
        STEP1_SELECTED,
        STEP1_SUMMARY,
        STEP1_DONE,
        STEP2_STATS,
        STEP2_SUMMARY,
        STEP2_DONE,
        STEP3_MATRIX,
        STEP3_CLUSTERS,
        STEP3_CORE_REP,
        STEP3_SUMMARY,
        STEP3_DONE,
        STEP4_MASTER,
        STEP4_TIER1_TSV,
        STEP4_TIER1_FASTA,
        STEP4_SUMMARY,
        STEP4_DONE,
        STEP5_MASTER,
        STEP5_CYTO_TSV,
        STEP5_CYTO_FASTA,
        STEP5_SURFACE_TSV,
        STEP5_SURFACE_FASTA,
        STEP5_PRIO_TSV,
        STEP5_PRIO_FASTA,
        STEP5_SUMMARY,
        STEP5_DONE,
        STEP6_MASTER,
        STEP6_DUAL_TSV,
        STEP6_DUAL_FASTA,
        STEP6_CLEAN_TSV,
        STEP6_CLEAN_FASTA,
        STEP6_SUMMARY,
        STEP6_DONE,
        STEP7_MASTER,
        STEP7_NOVEL_TSV,
        STEP7_NOVEL_FASTA,
        STEP7_APPROVED_TSV,
        STEP7_SUMMARY,
        STEP7_DONE,
        STEP8_MASTER,
        STEP8_DOUBLE_TSV,
        STEP8_GOLD_TSV,
        STEP8_GOLD_FASTA,
        STEP8_SUMMARY,
        STEP8_DONE,
        STEP9_MASTER,
        STEP9_HB_TSV,
        STEP9_LEADS_TSV,
        STEP9_EDGES_TSV,
        STEP9_HB_FASTA,
        STEP9_LEADS_FASTA,
        STEP9_SUMMARY,
        STEP9_DONE,
        JOB_METHODOLOGY,
        JOB_RESULTS

# ------------------------------------------------------------------------------
# Step 1: Data Acquisition
# ------------------------------------------------------------------------------
rule data_acquisition:
    """
    Step 1: Queries UniProt REST API using configuration parameters,
    downloads proteome metadata TSV and individual proteome FASTA files
    into results/{job_name}/01_data_acquisition/
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/01_data_acquisition.py"
    output:
        metadata=STEP1_METADATA,
        selected=STEP1_SELECTED,
        summary=STEP1_SUMMARY,
        done=STEP1_DONE
    threads:
        resolve_threads(config.get("download_options", {}).get("max_workers"), config)
    params:
        job_dir=JOB_DIR,
        fasta_dir=STEP1_FASTA_DIR
    log:
        f"{LOG_DIR}/01_data_acquisition.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 2: Redundancy Reduction & Paralog Removal (MMseqs2)
# ------------------------------------------------------------------------------
rule paralog_removal:
    """
    Step 2: Removes paralogous / redundant sequences per strain using MMseqs2,
    outputting non-paralogous FASTA files and cluster mapping tables.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/02_paralog_removal.py",
        step1_done=STEP1_DONE,
        selected=STEP1_SELECTED
    output:
        stats=STEP2_STATS,
        summary=STEP2_SUMMARY,
        done=STEP2_DONE
    threads:
        resolve_threads(config.get("step2_paralog_removal", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/02_paralog_removal.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 3: Pangenome Partitioning (Core / Shell / Cloud) via MMseqs2
# ------------------------------------------------------------------------------
rule pangenome_analysis:
    """
    Step 3: Clusters non-paralogous sequences across all strains, calculates
    the pangenome matrix, classifies into Core/Shell/Cloud, and extracts
    core representative sequences for downstream drug targeting.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/03_pangenome_analysis.py",
        step2_done=STEP2_DONE
    output:
        matrix=STEP3_MATRIX,
        clusters=STEP3_CLUSTERS,
        core_rep=STEP3_CORE_REP,
        summary=STEP3_SUMMARY,
        done=STEP3_DONE
    threads:
        resolve_threads(config.get("step3_pangenome", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/03_pangenome_analysis.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 4: Subtractive Genomics & Target Prioritization (DEG, Human, HRGM2)
# ------------------------------------------------------------------------------
rule subtractive_filter:
    """
    Step 4: Screens Core representative proteins through DEG Essentiality,
    Human Non-Homology, and HRGM2 Gut Microbiome preservation filters.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/04_subtractive_filter.py",
        core_rep=STEP3_CORE_REP,
        step3_done=STEP3_DONE
    output:
        master_tsv=STEP4_MASTER,
        tier1_tsv=STEP4_TIER1_TSV,
        tier1_fasta=STEP4_TIER1_FASTA,
        summary=STEP4_SUMMARY,
        done=STEP4_DONE
    threads:
        resolve_threads(config.get("step4_subtractive_filter", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/04_subtractive_filter.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 5: Subcellular Localization & Transmembrane Topology
# ------------------------------------------------------------------------------
rule subcellular_localization:
    """
    Step 5: Predicts subcellular compartments and transmembrane topology
    using a dual-engine consensus (UniProtKB curated + Kyte-Doolittle hydropathy
    and von Heijne signal peptide / positive-inside model).
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/05_subcellular_localization.py",
        step4_done=STEP4_DONE,
        tier1_tsv=STEP4_TIER1_TSV,
        tier1_fasta=STEP4_TIER1_FASTA,
        master_tsv=STEP4_MASTER,
        core_fasta=STEP3_CORE_REP
    output:
        master_tsv=STEP5_MASTER,
        cyto_tsv=STEP5_CYTO_TSV,
        cyto_fasta=STEP5_CYTO_FASTA,
        surface_tsv=STEP5_SURFACE_TSV,
        surface_fasta=STEP5_SURFACE_FASTA,
        prio_tsv=STEP5_PRIO_TSV,
        prio_fasta=STEP5_PRIO_FASTA,
        summary=STEP5_SUMMARY,
        done=STEP5_DONE
    threads:
        resolve_threads(config.get("step5_subcellular_localization", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/05_subcellular_localization.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 6: Virulence Factor (VFDB) & Antimicrobial Resistance (CARD) Screening
# ------------------------------------------------------------------------------
rule virulence_resistance_screening:
    """
    Step 6: Screens prioritized targets against VFDB (Virulence Factors) and
    CARD (Comprehensive Antibiotic Resistance Database) to identify dual-action
    super-targets, antibiotic adjuvants, and clean essential targets.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/06_virulence_resistance_screening.py",
        step5_done=STEP5_DONE,
        step5_prio_tsv=STEP5_PRIO_TSV,
        step5_prio_fasta=STEP5_PRIO_FASTA
    output:
        master=STEP6_MASTER,
        dual_tsv=STEP6_DUAL_TSV,
        dual_fasta=STEP6_DUAL_FASTA,
        clean_tsv=STEP6_CLEAN_TSV,
        clean_fasta=STEP6_CLEAN_FASTA,
        summary=STEP6_SUMMARY,
        done=STEP6_DONE
    threads:
        resolve_threads(config.get("step6_virulence_resistance", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/06_virulence_resistance_screening.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 7: Target Druggability & Novelty Assessment (TTD 2026)
# ------------------------------------------------------------------------------
rule target_druggability_ttd:
    """
    Step 7: Evaluates prioritized targets against Therapeutic Target Database
    (TTD 2026 update) to discover drug-repurposing candidates and identify
    completely novel, first-in-class patentable antibacterial targets.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/07_target_druggability_ttd.py",
        step6_done=STEP6_DONE,
        step6_master=STEP6_MASTER
    output:
        master=STEP7_MASTER,
        novel_tsv=STEP7_NOVEL_TSV,
        novel_fasta=STEP7_NOVEL_FASTA,
        approved_tsv=STEP7_APPROVED_TSV,
        summary=STEP7_SUMMARY,
        done=STEP7_DONE
    threads:
        resolve_threads(config.get("step7_target_druggability_ttd", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/07_target_druggability_ttd.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 8: Metabolic Pathway & Choke-Point Analysis (KEGG / BioCyc)
# ------------------------------------------------------------------------------
rule metabolic_chokepoint_analysis:
    """
    Step 8: Maps prioritized targets into the complete stoichiometric metabolic
    network of the pathogen to identify unbypassable metabolic choke-points
    and pathogen-exclusive biochemical pathways absent in humans.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/08_metabolic_chokepoint_analysis.py",
        step7_done=STEP7_DONE,
        step7_master=STEP7_MASTER
    output:
        master=STEP8_MASTER,
        double_tsv=STEP8_DOUBLE_TSV,
        gold_tsv=STEP8_GOLD_TSV,
        gold_fasta=STEP8_GOLD_FASTA,
        summary=STEP8_SUMMARY,
        done=STEP8_DONE
    threads:
        resolve_threads(config.get("step8_metabolic_chokepoint", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/08_metabolic_chokepoint_analysis.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Step 9: Protein-Protein Interactome & Network Hub Analysis (STRING v12)
# ------------------------------------------------------------------------------
rule network_biology_string:
    """
    Step 9: Models the pathogen's high-confidence proteome interactome via
    STRING v12, computes NetworkX Degree, Betweenness, and Closeness centralities,
    and identifies master Hub-Bottleneck regulators per the Centrality-Lethality Rule.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/09_network_biology_string.py",
        step8_done=STEP8_DONE,
        step8_master=STEP8_MASTER
    output:
        master=STEP9_MASTER,
        hb_tsv=STEP9_HB_TSV,
        leads_tsv=STEP9_LEADS_TSV,
        edges_tsv=STEP9_EDGES_TSV,
        hb_fasta=STEP9_HB_FASTA,
        leads_fasta=STEP9_LEADS_FASTA,
        summary=STEP9_SUMMARY,
        done=STEP9_DONE
    threads:
        resolve_threads(config.get("step9_network_biology_string", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/09_network_biology_string.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Publication-Ready Documentation (METHODOLOGY.md & RESULTS.md)
# ------------------------------------------------------------------------------
rule generate_reports:
    """
    Generates publication-ready METHODOLOGY.md and RESULTS.md reports
    directly inside the job's results folder, customized with empirical multi-omics data.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/generate_job_reports.py",
        step1_summary=STEP1_SUMMARY,
        step2_summary=STEP2_SUMMARY,
        step3_summary=STEP3_SUMMARY,
        step4_summary=STEP4_SUMMARY,
        step5_summary=STEP5_SUMMARY,
        step6_summary=STEP6_SUMMARY,
        step7_summary=STEP7_SUMMARY,
        step8_summary=STEP8_SUMMARY,
        step9_summary=STEP9_SUMMARY,
        step9_done=STEP9_DONE
    output:
        methodology=JOB_METHODOLOGY,
        results=JOB_RESULTS
    threads: 1
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/generate_reports.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Publication-Ready Visualizations & CSV Exports (Steps 1 through 9)
# ------------------------------------------------------------------------------
rule publication_figures:
    """
    Generates publication-quality charts (PNG 300 DPI, vector PDF)
    and standard CSV files across Steps 1 through 9.
    """
    input:
        config=ACTIVE_CONFIG,
        script="scripts/generate_publication_figures.py",
        step1_meta=STEP1_METADATA,
        step2_stats=STEP2_STATS,
        step3_matrix=STEP3_MATRIX,
        step4_done=STEP4_DONE,
        step5_done=STEP5_DONE,
        step6_done=STEP6_DONE,
        step7_done=STEP7_DONE,
        step8_done=STEP8_DONE,
        step9_done=STEP9_DONE
    output:
        overview_fig=MASTER_FIG
    threads:
        resolve_threads(config.get("hardware", {}).get("threads"), config)
    params:
        job_dir=JOB_DIR
    log:
        f"{LOG_DIR}/publication_figures.log"
    shell:
        """
        mkdir -p $(dirname {log})
        {PYTHON} {input.script} --config {input.config} 2>&1 | tee {log}
        """

# ------------------------------------------------------------------------------
# Downstream Roadmap (Steps 2-11 will attach here as we build them)
# Step 2: Redundancy Reduction & Paralog Removal (CD-HIT)
# Step 3: Pangenome Analysis & Partitioning (Core / Shell / Cloud)
# Step 4: Metabolic Pathway & Choke-point Prioritization
# Step 5: Pathogenicity Screening (VFDB & CARD)
# Step 6: Target Pool Integration
# Step 7: Essentiality Screening (DEG)
# Step 8: Host & Gut Flora Non-Homology Filtering
# Step 9: Subcellular Localization & Druggability
# Step 10: Target Novelty Assessment (TTD)
# Step 11: Network Systems Biology (STRING PPI & Hubs)
# ------------------------------------------------------------------------------
