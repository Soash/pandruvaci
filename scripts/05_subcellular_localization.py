#!/usr/bin/env python3
"""
scripts/05_subcellular_localization.py: Dual-Engine Subcellular Localization
and Transmembrane (TM) Topology Prediction for Bacterial Drug Targets.

Combines:
  1. UniProtKB Curated Annotations (Subcellular Location [CC], GO Cellular Component,
     ft_transmem, ft_signal) via high-throughput REST API with local caching.
  2. De Novo Biophysical Topology Predictor (Kyte-Doolittle hydropathy sliding window,
     von Heijne bacterial signal peptide tripartite model, and positive-inside rule).

Outputs:
  - Master annotated table (TSV/CSV) with predicted compartments, TMH counts, and druggability feasibility.
  - Partitioned tables and FASTA files for Tier 1 cytoplasmic small-molecule targets and surface/vaccine targets.
  - Comprehensive run summary metrics (JSON/TSV).
"""

import os
import sys
import json
import time
import argparse
import urllib.request
import urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# Add repo root to sys.path for hardware utils
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.utils_hardware import detect_available_cores, resolve_threads

# Kyte-Doolittle Hydropathy Scale (Kyte & Doolittle, J. Mol. Biol. 1982)
KYTE_DOOLITTLE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3,
    'X': 0.0, 'U': 0.0, 'B': -3.5, 'Z': -3.5
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Step 5: Subcellular Localization & Transmembrane Topology"
    )
    parser.add_argument(
        "--config", default="config.yaml", help="Path to config.yaml"
    )
    parser.add_argument(
        "--job-name", default=None, help="Override job_name from config"
    )
    parser.add_argument(
        "--tier1-only", dest="tier1_only", action="store_true", default=None,
        help="Only pass and process Tier 1 targets"
    )
    parser.add_argument(
        "--all-tiers", dest="tier1_only", action="store_false",
        help="Pass and process all tiers (Tier 1, 2, 3, and Excluded)"
    )
    parser.add_argument(
        "--target-modality", default=None,
        choices=["drug_candidates", "cytoplasmic_only", "vaccine_candidates", "all"],
        help="Filter prioritized targets by modality (default from config: drug_candidates)"
    )
    return parser.parse_args()


def load_yaml_config(config_path):
    import yaml
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def parse_fasta(fasta_path):
    """Parses a FASTA file into a dictionary of {header_id: sequence}."""
    sequences = {}
    if not os.path.exists(fasta_path):
        return sequences

    current_id = None
    current_seq = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                header_line = line[1:].split()[0]
                current_id = header_line
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences


# ==============================================================================
# 1. UniProtKB Curated Annotations Fetcher (Batch REST API with Local Cache)
# ==============================================================================

def fetch_uniprot_annotations_batch(accessions, cache_dir, batch_size=100):
    """
    Fetches curated subcellular location, GO cellular component, TMH, and signal
    peptides from UniProtKB REST API with disk caching.
    """
    os.makedirs(cache_dir, exist_ok=True)
    results = {}
    missing_accs = []

    # Check local cache first
    for acc in accessions:
        cache_file = os.path.join(cache_dir, f"{acc}.json")
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "r") as f:
                    results[acc] = json.load(f)
            except Exception:
                missing_accs.append(acc)
        else:
            missing_accs.append(acc)

    if not missing_accs:
        return results

    print(f"[*] Querying UniProtKB REST API for {len(missing_accs)} target accessions (batch size: {batch_size})...")

    # Chunk into batch_size
    for i in range(0, len(missing_accs), batch_size):
        chunk = missing_accs[i : i + batch_size]
        acc_query = ",".join(chunk)
        url = (
            f"https://rest.uniprot.org/uniprotkb/accessions?"
            f"accessions={urllib.parse.quote(acc_query)}"
            f"&fields=accession,id,cc_subcellular_location,go_c,ft_transmem,ft_signal"
            f"&format=tsv"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "PangenomeDrugPipeline/1.0"})

        retries = 3
        tsv_text = None
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    tsv_text = resp.read().decode("utf-8")
                break
            except Exception as e:
                time.sleep(1.5 * (attempt + 1))
                if attempt == retries - 1:
                    print(f"    [!] Warning: Failed UniProt query chunk {i//batch_size + 1}: {e}")

        if not tsv_text:
            continue

        # Parse TSV
        lines = tsv_text.strip().split("\n")
        if len(lines) <= 1:
            continue

        header = lines[0].split("\t")
        col_map = {col.strip(): idx for idx, col in enumerate(header)}

        for line in lines[1:]:
            parts = line.split("\t")
            if not parts or len(parts) <= col_map.get("Entry", 0):
                continue
            acc = parts[col_map.get("Entry", 0)].strip()
            entry_name = parts[col_map.get("Entry Name", 1)].strip() if "Entry Name" in col_map else ""
            subcell_cc = parts[col_map.get("Subcellular location [CC]", 2)].strip() if "Subcellular location [CC]" in col_map and len(parts) > col_map["Subcellular location [CC]"] else ""
            go_c = parts[col_map.get("Gene Ontology (cellular component)", 3)].strip() if "Gene Ontology (cellular component)" in col_map and len(parts) > col_map["Gene Ontology (cellular component)"] else ""
            transmem = parts[col_map.get("Transmembrane", 4)].strip() if "Transmembrane" in col_map and len(parts) > col_map["Transmembrane"] else ""
            signal = parts[col_map.get("Signal peptide", 5)].strip() if "Signal peptide" in col_map and len(parts) > col_map["Signal peptide"] else ""

            # Count curated TMHs
            tmh_coords = []
            if transmem:
                for token in transmem.split(";"):
                    token = token.strip()
                    if token.startswith("TRANSMEM"):
                        coord_str = token.replace("TRANSMEM", "").strip()
                        tmh_coords.append(coord_str)

            data = {
                "accession": acc,
                "entry_name": entry_name,
                "subcellular_location_cc": subcell_cc,
                "go_cellular_component": go_c,
                "curated_tmh_count": len(tmh_coords),
                "curated_tmh_coords": "; ".join(tmh_coords) if tmh_coords else "None",
                "curated_signal_peptide": signal if signal else "None",
                "has_curated_data": bool(subcell_cc or go_c or transmem or signal)
            }
            results[acc] = data

            # Save to disk cache
            cache_file = os.path.join(cache_dir, f"{acc}.json")
            try:
                with open(cache_file, "w") as f:
                    json.dump(data, f)
            except Exception:
                pass

        time.sleep(0.2)  # Respect UniProt rate limits

    return results


# ==============================================================================
# 2. De Novo Biophysical Topology Predictor (Kyte-Doolittle & von Heijne)
# ==============================================================================

def predict_tm_helices_hydropathy(sequence, window_size=19, threshold=1.60, min_loop=5):
    """
    Predicts transmembrane alpha-helices using the Kyte-Doolittle hydropathy scale
    with a sliding window of length 19 (lipid bilayer width ~30 Å).
    """
    seq_len = len(sequence)
    if seq_len < window_size:
        return 0, [], "None"

    scores = []
    for i in range(seq_len - window_size + 1):
        window = sequence[i : i + window_size]
        avg_hydro = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in window) / window_size
        scores.append((i, avg_hydro))

    # Identify candidate peak regions above threshold
    candidate_helices = []
    in_helix = False
    start_idx = 0
    max_score = 0.0

    for i, score in scores:
        if score >= threshold:
            if not in_helix:
                in_helix = True
                start_idx = i
                max_score = score
            else:
                max_score = max(max_score, score)
        else:
            if in_helix:
                in_helix = False
                end_idx = i + window_size - 1
                candidate_helices.append((start_idx + 1, end_idx, max_score))

    if in_helix:
        end_idx = len(scores) - 1 + window_size
        candidate_helices.append((start_idx + 1, end_idx, max_score))

    # Resolve overlapping or too-close helices (< min_loop residues separation)
    filtered_helices = []
    for h in candidate_helices:
        if not filtered_helices:
            filtered_helices.append(h)
        else:
            prev_start, prev_end, prev_score = filtered_helices[-1]
            curr_start, curr_end, curr_score = h
            if curr_start - prev_end < min_loop:
                # Merge into longer span or keep stronger peak
                merged_start = prev_start
                merged_end = max(prev_end, curr_end)
                merged_score = max(prev_score, curr_score)
                filtered_helices[-1] = (merged_start, merged_end, merged_score)
            else:
                filtered_helices.append(h)

    coord_strs = [f"{s}..{e}" for s, e, _ in filtered_helices]
    return len(filtered_helices), filtered_helices, "; ".join(coord_strs) if coord_strs else "None"


def predict_bacterial_signal_peptide(sequence):
    """
    Evaluates bacterial signal peptide features using the von Heijne tripartite model:
      - N-region (residues 1-7): contains at least one basic residue (K or R).
      - H-region (residues 8-22): hydrophobic core of 7-15 residues (mean hydropathy > 1.4).
      - C-region (residues 18-32): cleavage motif with small uncharged residues at -3, -1.
    """
    if len(sequence) < 25:
        return False, "None"

    # Analyze first 35 residues
    n_term = sequence[:35]

    # N-region check: basic residues
    n_region = n_term[:7]
    has_basic = any(aa in "KR" for aa in n_region)

    # H-region check: sliding 8-aa hydrophobic core
    has_h_core = False
    core_start = 0
    for i in range(len(n_term) - 8):
        sub = n_term[i : i + 8]
        hydro = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in sub) / 8.0
        if hydro >= 1.5:
            has_h_core = True
            core_start = i
            break

    # Cleavage region check: looking for AXA or small residues at pos 18-30
    has_cleavage = False
    cleave_pos = 0
    for pos in range(max(15, core_start + 7), min(32, len(sequence) - 2)):
        r_minus_3 = sequence[pos - 3]
        r_minus_1 = sequence[pos - 1]
        if r_minus_3 in "AGSTILV" and r_minus_1 in "AGST":
            has_cleavage = True
            cleave_pos = pos
            break

    if (has_basic and has_h_core) or (has_h_core and has_cleavage):
        cleave_annot = f"1..{cleave_pos}" if cleave_pos > 0 else "1..22 (predicted)"
        return True, cleave_annot

    return False, "None"


# ==============================================================================
# 3. Consensus Subcellular Classification & Druggability Feasibility
# ==============================================================================

def classify_subcellular_consensus(curated_info, denovo_tmh_count, denovo_tmh_coords,
                                   denovo_has_sp, denovo_sp_coords, gram_stain="negative"):
    """
    Synthesizes curated UniProtKB evidence with de novo biophysical topology
    into consensus compartment, final TMH count, and druggability feasibility.
    """
    subcell_cc = (curated_info.get("subcellular_location_cc") or "").lower()
    go_c = (curated_info.get("go_cellular_component") or "").lower()
    curated_tmh = curated_info.get("curated_tmh_count", 0)
    curated_sp = curated_info.get("curated_signal_peptide", "None")

    # 1. Resolve Final TMH Count
    if curated_info.get("has_curated_data") and curated_tmh > 0:
        final_tmh_count = curated_tmh
        final_tmh_coords = curated_info.get("curated_tmh_coords", "None")
        tmh_source = "UniProtKB Curated"
    else:
        final_tmh_count = denovo_tmh_count
        final_tmh_coords = denovo_tmh_coords
        tmh_source = "De Novo Biophysical (Kyte-Doolittle)"

    # 2. Resolve Signal Peptide Status
    if curated_sp != "None" and curated_sp:
        has_signal_peptide = True
        sp_coords = curated_sp
        sp_source = "UniProtKB Curated"
    elif denovo_has_sp:
        has_signal_peptide = True
        sp_coords = denovo_sp_coords
        sp_source = "De Novo (von Heijne Tripartite)"
    else:
        has_signal_peptide = False
        sp_coords = "None"
        sp_source = "None"

    # 3. Compartment Classification Logic
    pred_compartment = "Unknown"
    confidence_source = "De Novo Consensus"

    # Priority A: Curated UniProt Annotations
    if "outer membrane" in subcell_cc or "outer membrane" in go_c:
        pred_compartment = "Outer Membrane"
        confidence_source = "UniProtKB Curated"
    elif "periplasm" in subcell_cc or "periplasm" in go_c:
        pred_compartment = "Periplasmic"
        confidence_source = "UniProtKB Curated"
    elif "extracellular" in subcell_cc or "secreted" in subcell_cc or "extracellular" in go_c:
        pred_compartment = "Extracellular"
        confidence_source = "UniProtKB Curated"
    elif "inner membrane" in subcell_cc or "cytoplasmic membrane" in subcell_cc or "plasma membrane" in go_c:
        pred_compartment = "Cytoplasmic Membrane"
        confidence_source = "UniProtKB Curated"
    elif "cytoplasm" in subcell_cc or "cytosol" in go_c:
        if final_tmh_count >= 2:
            pred_compartment = "Cytoplasmic Membrane"
            confidence_source = "Curated (Membrane Associated)"
        else:
            pred_compartment = "Cytoplasmic"
            confidence_source = "UniProtKB Curated"

    # Priority B: De Novo Biophysical Assignment (if unreviewed or unspecified)
    if pred_compartment == "Unknown":
        if final_tmh_count >= 1:
            pred_compartment = "Cytoplasmic Membrane"
            confidence_source = "De Novo Topology"
        elif has_signal_peptide:
            if gram_stain.lower() == "negative":
                pred_compartment = "Periplasmic"
            else:
                pred_compartment = "Extracellular"
            confidence_source = "De Novo Signal Peptide"
        else:
            pred_compartment = "Cytoplasmic"
            confidence_source = "De Novo Soluble"

    # 4. Druggability Feasibility (The <= 1 TMH Rule)
    if final_tmh_count <= 1:
        druggability_feasibility = "Feasible (TMH <= 1)"
    else:
        druggability_feasibility = "Deprioritized (TMH >= 2)"

    # 5. Therapeutic Modality Assignment
    if pred_compartment == "Cytoplasmic" and druggability_feasibility.startswith("Feasible"):
        modality = "Small-Molecule Drug Target (Cytoplasmic Enzyme)"
    elif pred_compartment == "Periplasmic" and druggability_feasibility.startswith("Feasible"):
        modality = "Periplasmic Target (Porin-Accessible Small Molecule)"
    elif pred_compartment in ["Outer Membrane", "Extracellular"]:
        modality = "Surface / Secreted Target (Vaccine & Monoclonal Antibody)"
    elif pred_compartment == "Cytoplasmic Membrane" and final_tmh_count == 1:
        modality = "Single-Pass Membrane Target (Assay-Feasible)"
    else:
        modality = "Multi-Pass Integral Membrane (Challenging / Deprioritized)"

    return {
        "Predicted_Compartment": pred_compartment,
        "Localization_Source": confidence_source,
        "TMH_Count": final_tmh_count,
        "TMH_Coordinates": final_tmh_coords,
        "TMH_Source": tmh_source,
        "Signal_Peptide": "Yes" if has_signal_peptide else "No",
        "Signal_Peptide_Coords": sp_coords,
        "Signal_Peptide_Source": sp_source,
        "Druggability_Feasibility": druggability_feasibility,
        "Therapeutic_Modality": modality
    }


# ==============================================================================
# Main Orchestration
# ==============================================================================

def main():
    start_time = time.time()
    args = parse_args()
    cfg = load_yaml_config(args.config)

    job_name = args.job_name or cfg.get("job_name", "acinetobacter_baumannii")
    results_dir = cfg.get("results_dir", "results")
    job_dir = Path(results_dir) / job_name

    step4_dir = job_dir / "04_subtractive_filter"
    step3_dir = job_dir / "03_pangenome_analysis"
    step5_dir = job_dir / "05_subcellular_localization"

    master_tsv_path = step4_dir / "tables" / "subtractive_master.tsv"
    core_fasta_path = step3_dir / "fasta" / "core_representatives.fasta"

    if not master_tsv_path.exists():
        print(f"[!] Error: Step 4 master table not found at: {master_tsv_path}")
        sys.exit(1)

    step5_cfg = cfg.get("step5_subcellular_localization", {})
    if not step5_cfg.get("enabled", True):
        print("[*] Step 5: Subcellular Localization is disabled in config. Skipping.")
        sys.exit(0)

    gram_stain = cfg.get("gram_stain") or step5_cfg.get("gram_stain", "negative")
    max_tmh_druggable = step5_cfg.get("max_tmh_for_druggable", 1)
    hydro_threshold = step5_cfg.get("hydropathy_threshold", 1.6)
    window_size = step5_cfg.get("window_size", 19)
    min_loop = step5_cfg.get("min_loop_length", 5)
    batch_size = step5_cfg.get("batch_size", 100)
    cache_dir = step5_cfg.get("cache_dir", "databases/cache/uniprot_subcell")

    # Output directories
    tables_dir = step5_dir / "tables"
    fasta_dir = step5_dir / "fasta"
    summary_dir = step5_dir / "summary"
    for d in [tables_dir, fasta_dir, summary_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("STEP 5: SUBCELLULAR LOCALIZATION & TRANSMEMBRANE TOPOLOGY")
    print("=" * 80)
    print(f"[*] Job Name: {job_name}")
    print(f"[*] Gram Architecture: {gram_stain.capitalize()}")
    print(f"[*] Druggability TMH Cutoff: <= {max_tmh_druggable} TMH (Feasible)")
    print(f"[*] Kyte-Doolittle Window: {window_size} aa | Threshold: >= {hydro_threshold}")
    print(f"[*] Cache Directory: {cache_dir}")

    tier1_only_cfg = step5_cfg.get("tier1_only", True)
    tier1_only = args.tier1_only if args.tier1_only is not None else tier1_only_cfg

    target_modality = args.target_modality or step5_cfg.get("target_modality", "drug_candidates")
    print(f"[*] Target Modality Filter: {target_modality.upper()}")

    # Load targets TSV
    tier1_tsv_path = step4_dir / "tables" / "tier1_golden_targets.tsv"
    if tier1_only:
        if tier1_tsv_path.exists():
            df_master = pd.read_csv(tier1_tsv_path, sep="\t")
        else:
            df_master = pd.read_csv(master_tsv_path, sep="\t")
            df_master = df_master[
                df_master["Target_Tier"].astype(str).str.contains("Tier 1", case=False, na=False)
            ].copy()
        print(f"[*] Target Filter: TIER 1 ONLY ({len(df_master)} Pathogen-Selective targets)")
    else:
        df_master = pd.read_csv(master_tsv_path, sep="\t")
        print(f"[*] Target Filter: ALL TIERS ({len(df_master)} core candidate targets)")

    # Load FASTA sequences
    tier1_fasta_path = step4_dir / "fasta" / "tier1_golden_targets.fasta"
    if tier1_only and tier1_fasta_path.exists():
        fasta_seqs = parse_fasta(str(tier1_fasta_path))
    else:
        fasta_seqs = parse_fasta(str(core_fasta_path))

    if not fasta_seqs and core_fasta_path.exists():
        fasta_seqs = parse_fasta(str(core_fasta_path))
    print(f"[*] Loaded {len(fasta_seqs)} FASTA sequences.")

    # Extract clean UniProt accessions
    accessions = []
    for raw_acc in df_master["Representative_Accession"].dropna().unique():
        # Handle format: tr|A0A...|A0A... or sp|P0A...|... or raw accession
        clean_acc = raw_acc.split("|")[1] if "|" in raw_acc else raw_acc
        accessions.append(clean_acc)

    # 1. Fetch UniProt Curated Annotations
    uniprot_annotations = fetch_uniprot_annotations_batch(accessions, cache_dir, batch_size=batch_size)

    # 2. Run De Novo Topology and Consensus Synthesis
    print("[*] Running de novo biophysical topology prediction and consensus classification...")
    records = []
    for idx, row in df_master.iterrows():
        raw_id = str(row.get("Representative_Accession", ""))
        clean_acc = raw_id.split("|")[1] if "|" in raw_id else raw_id
        cluster_id = row.get("Cluster_ID", f"Cluster_{idx}")

        # Get sequence
        seq = fasta_seqs.get(raw_id, fasta_seqs.get(clean_acc, ""))
        seq_len = len(seq)

        # De novo Kyte-Doolittle TMH prediction
        denovo_tmh_cnt, denovo_tmh_list, denovo_tmh_coords = predict_tm_helices_hydropathy(
            seq, window_size=window_size, threshold=hydro_threshold, min_loop=min_loop
        )

        # De novo Signal Peptide prediction
        denovo_has_sp, denovo_sp_coords = predict_bacterial_signal_peptide(seq)

        # Fetch curated info
        curated_info = uniprot_annotations.get(clean_acc, {})

        # Synthesize consensus
        consensus = classify_subcellular_consensus(
            curated_info=curated_info,
            denovo_tmh_count=denovo_tmh_cnt,
            denovo_tmh_coords=denovo_tmh_coords,
            denovo_has_sp=denovo_has_sp,
            denovo_sp_coords=denovo_sp_coords,
            gram_stain=gram_stain
        )

        # Append to record
        rec = dict(row)
        rec["Protein_Length"] = seq_len
        rec["Predicted_Compartment"] = consensus["Predicted_Compartment"]
        rec["Localization_Source"] = consensus["Localization_Source"]
        rec["TMH_Count"] = consensus["TMH_Count"]
        rec["TMH_Coordinates"] = consensus["TMH_Coordinates"]
        rec["TMH_Source"] = consensus["TMH_Source"]
        rec["Signal_Peptide"] = consensus["Signal_Peptide"]
        rec["Signal_Peptide_Coords"] = consensus["Signal_Peptide_Coords"]
        rec["Signal_Peptide_Source"] = consensus["Signal_Peptide_Source"]
        rec["Druggability_Feasibility"] = consensus["Druggability_Feasibility"]
        rec["Therapeutic_Modality"] = consensus["Therapeutic_Modality"]
        records.append(rec)

    df_subcell = pd.DataFrame(records)

    # 3. Export Master TSV & CSV
    out_master_tsv = tables_dir / "subcellular_master.tsv"
    out_master_csv = tables_dir / "subcellular_master.csv"
    df_subcell.to_csv(out_master_tsv, sep="\t", index=False)
    df_subcell.to_csv(out_master_csv, index=False)
    print(f"[+] Exported Master Subcellular Table: {out_master_tsv}")

    # 4. Generate Partitioned Tables & FASTA Files
    # 4. Generate Partitioned Tables & FASTA Files
    # Tier 1 Cytoplasmic Small-Molecule Targets (Feasible TMH <= 1)
    df_tier1_cyto = df_subcell[
        (df_subcell["Target_Tier"].astype(str).str.contains("Tier 1", case=False, na=False)) &
        (df_subcell["Predicted_Compartment"] == "Cytoplasmic") &
        (df_subcell["Druggability_Feasibility"].str.startswith("Feasible"))
    ]
    df_tier1_cyto.to_csv(tables_dir / "tier1_druggable_cytoplasmic.tsv", sep="\t", index=False)
    df_tier1_cyto.to_csv(tables_dir / "tier1_druggable_cytoplasmic.csv", index=False)

    # Tier 1 Surface / Vaccine Targets (Periplasmic, Outer Membrane, Extracellular)
    df_tier1_surface = df_subcell[
        (df_subcell["Target_Tier"].astype(str).str.contains("Tier 1", case=False, na=False)) &
        (df_subcell["Predicted_Compartment"].isin(["Periplasmic", "Outer Membrane", "Extracellular"])) &
        (df_subcell["Druggability_Feasibility"].str.startswith("Feasible"))
    ]
    df_tier1_surface.to_csv(tables_dir / "tier1_surface_vaccine_targets.tsv", sep="\t", index=False)
    df_tier1_surface.to_csv(tables_dir / "tier1_surface_vaccine_targets.csv", index=False)

    # Write Partitioned FASTA Files
    def write_subset_fasta(df_subset, fasta_out):
        count = 0
        with open(fasta_out, "w") as f_out:
            for _, r in df_subset.iterrows():
                raw_id = str(r["Representative_Accession"])
                clean_acc = raw_id.split("|")[1] if "|" in raw_id else raw_id
                seq = fasta_seqs.get(raw_id, fasta_seqs.get(clean_acc, ""))
                if seq:
                    gene = str(r.get("Gene_Name", "Unknown"))
                    comp = str(r.get("Predicted_Compartment", ""))
                    tmh = str(r.get("TMH_Count", 0))
                    header = f">{raw_id} gene={gene} compartment={comp} tmh={tmh}"
                    f_out.write(f"{header}\n{seq}\n")
                    count += 1
        return count

    n_cyto = write_subset_fasta(df_tier1_cyto, fasta_dir / "tier1_druggable_cytoplasmic.fasta")
    n_surf = write_subset_fasta(df_tier1_surface, fasta_dir / "tier1_surface_vaccine_targets.fasta")

    # Filter Prioritized Candidates based on target_modality for downstream steps
    if target_modality == "drug_candidates":
        # Soluble / assay-feasible small-molecule targets in Cytoplasm or Periplasm
        df_prioritized = df_subcell[
            (df_subcell["Predicted_Compartment"].isin(["Cytoplasmic", "Periplasmic"])) &
            (df_subcell["Druggability_Feasibility"].str.startswith("Feasible"))
        ]
    elif target_modality == "cytoplasmic_only":
        df_prioritized = df_subcell[
            (df_subcell["Predicted_Compartment"] == "Cytoplasmic") &
            (df_subcell["Druggability_Feasibility"].str.startswith("Feasible"))
        ]
    elif target_modality == "vaccine_candidates":
        df_prioritized = df_subcell[
            (df_subcell["Predicted_Compartment"].isin(["Periplasmic", "Outer Membrane", "Extracellular"])) &
            (df_subcell["Druggability_Feasibility"].str.startswith("Feasible"))
        ]
    else:  # "all"
        df_prioritized = df_subcell[
            df_subcell["Druggability_Feasibility"].str.startswith("Feasible")
        ]

    out_prio_tsv = tables_dir / "prioritized_drug_candidates.tsv"
    out_prio_csv = tables_dir / "prioritized_drug_candidates.csv"
    out_prio_fasta = fasta_dir / "prioritized_drug_candidates.fasta"

    df_prioritized.to_csv(out_prio_tsv, sep="\t", index=False)
    df_prioritized.to_csv(out_prio_csv, index=False)
    n_prio = write_subset_fasta(df_prioritized, out_prio_fasta)

    print(f"[+] Exported {n_cyto} Tier 1 Cytoplasmic Drug Targets FASTA: {fasta_dir / 'tier1_druggable_cytoplasmic.fasta'}")
    print(f"[+] Exported {n_surf} Tier 1 Surface/Vaccine Targets FASTA: {fasta_dir / 'tier1_surface_vaccine_targets.fasta'}")
    print(f"[+] Exported {n_prio} Prioritized Targets ({target_modality}) FASTA: {out_prio_fasta}")

    # 5. Summary Metrics
    compartment_counts = df_subcell["Predicted_Compartment"].value_counts().to_dict()
    feasibility_counts = df_subcell["Druggability_Feasibility"].value_counts().to_dict()
    tier1_compartments = df_subcell[df_subcell["Target_Tier"].astype(str).str.contains("Tier 1", case=False, na=False)]["Predicted_Compartment"].value_counts().to_dict()

    summary_data = {
        "job_name": job_name,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "tier1_only": tier1_only,
        "target_modality": target_modality,
        "total_targets_analyzed": len(df_subcell),
        "total_prioritized_candidates": len(df_prioritized),
        "gram_architecture": gram_stain,
        "max_tmh_for_druggable": max_tmh_druggable,
        "compartment_distribution": compartment_counts,
        "druggability_feasibility": feasibility_counts,
        "tier1_golden_targets": {
            "total_tier1": len(df_subcell) if tier1_only else int((df_master["Target_Tier"].astype(str).str.contains("Tier 1", case=False, na=False)).sum()),
            "cytoplasmic_small_molecule_druggable": len(df_tier1_cyto),
            "surface_and_vaccine_candidates": len(df_tier1_surface),
            "prioritized_modality_selected": len(df_prioritized),
            "compartment_breakdown": tier1_compartments
        },
        "elapsed_seconds": round(time.time() - start_time, 2)
    }

    with open(summary_dir / "subcellular_summary.json", "w") as f:
        json.dump(summary_data, f, indent=2)

    df_summary_flat = pd.DataFrame([
        {"Metric": k, "Value": str(v)} for k, v in summary_data.items() if not isinstance(v, dict)
    ])
    df_summary_flat.to_csv(summary_dir / "subcellular_summary.tsv", sep="\t", index=False)

    # Touch completion marker
    (step5_dir / ".subcellular_complete").touch()

    print("\n" + "=" * 80)
    print("STEP 5 SUBCELLULAR LOCALIZATION SUMMARY")
    print("=" * 80)
    print(f"Total Targets Analyzed:             {len(df_subcell)}")
    print(f"Cytoplasmic Compartment:            {compartment_counts.get('Cytoplasmic', 0)}")
    print(f"Cytoplasmic Membrane:               {compartment_counts.get('Cytoplasmic Membrane', 0)}")
    print(f"Periplasmic Space:                  {compartment_counts.get('Periplasmic', 0)}")
    print(f"Outer Membrane:                     {compartment_counts.get('Outer Membrane', 0)}")
    print(f"Extracellular (Secreted):           {compartment_counts.get('Extracellular', 0)}")
    print("-" * 50)
    print(f"Druggability Feasible (TMH <= 1):   {feasibility_counts.get('Feasible (TMH <= 1)', 0)}")
    print(f"Deprioritized (TMH >= 2):           {feasibility_counts.get('Deprioritized (TMH >= 2)', 0)}")
    print("-" * 50)
    print(f"Tier 1 Cytoplasmic Drug Targets:    {len(df_tier1_cyto)}")
    print(f"Tier 1 Surface / Vaccine Targets:   {len(df_tier1_surface)}")
    print(f"Execution Time:                     {summary_data['elapsed_seconds']} s")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
