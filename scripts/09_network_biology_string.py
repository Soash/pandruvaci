#!/usr/bin/env python3
"""
Step 9: Protein-Protein Interactome & Network Hub-Bottleneck Analysis (STRING v12)
=================================================================================
Models the whole-proteome high-confidence interactome of the bacterial pathogen,
computes graph topological centrality metrics (Degree, Betweenness, Closeness,
and Clustering Coefficient), and classifies candidate targets into four quadrants:
  - Hub-Bottleneck (HB): Master topological regulators (Tier 1 Lethal Vulnerabilities)
  - Non-Hub Bottleneck (NH-B): Inter-modular signaling connectors
  - Hub Non-Bottleneck (H-NB): Intra-modular complex stabilizers
  - Non-Hub Non-Bottleneck (NH-NB): Peripheral effectors

Academic Citations:
- Centrality-Lethality Rule: Barabási AL, Oltvai ZN. Nat Rev Genet (2004) 5:101–113.
- Network Bottlenecks: Yu H, et al. Science (2007) 318:1458–1461.
- STRING Database v12: Szklarczyk D, et al. Nucleic Acids Res (2023) 51:D638–D646.

Usage:
    python scripts/09_network_biology_string.py [--config config.yaml]
"""

import os
import re
import sys
import gzip
import time
import json
import logging
import argparse
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, Set, List
import pandas as pd
import numpy as np
import networkx as nx
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
logger = logging.getLogger("09_network_biology")


def download_file_cached(url: str, dest_path: Path, description: str, retries: int = 3, timeout: int = 60) -> Path:
    """Downloads a file with caching, progress logging, and retries."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if dest_path.exists() and dest_path.stat().st_size > 0:
        logger.info("Using cached %s (%s, %.2f MB)", description, dest_path.name, dest_path.stat().st_size / (1024 * 1024))
        return dest_path

    logger.info("Downloading %s from %s...", description, url)
    temp_path = dest_path.with_suffix(dest_path.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": "PangenomeDrugPipeline/1.0 (Academic Research)"})
    last_err = None

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(temp_path, "wb") as out_f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out_f.write(chunk)
            temp_path.rename(dest_path)
            logger.info("Successfully downloaded %s (%.2f MB)", description, dest_path.stat().st_size / (1024 * 1024))
            return dest_path
        except Exception as e:
            last_err = e
            logger.warning("Attempt %d/%d failed downloading %s: %s", attempt, retries, description, e)
            if temp_path.exists():
                temp_path.unlink()
            time.sleep(attempt * 2)

    raise RuntimeError(f"Failed to download {description} from {url}: {last_err}")


def resolve_string_species_taxon_id(cfg: Dict[str, Any], cache_dir: Path) -> Tuple[int, str]:
    """
    Dynamically resolves the NCBI Taxon ID used by STRING v12 for the pathogen.
    Works for any bacterial species without hardcoding.
    """
    tax_id = cfg.get("taxonomy_id") or cfg.get("query_params", {}).get("taxonomy_id", "")
    job_name = cfg.get("job_name", "klebsiella_pneumoniae")
    species_name = job_name.replace("_", " ").strip()
    cache_file = cache_dir / f"string_taxon_resolution_{job_name}.json"

    if cache_file.exists():
        try:
            with open(cache_file, "r") as f:
                cached = json.load(f)
                if cached.get("resolved_taxon_id"):
                    logger.info("Loaded cached STRING taxon resolution: %s (%s)", cached["resolved_taxon_id"], cached.get("taxon_name"))
                    return int(cached["resolved_taxon_id"]), cached.get("taxon_name", "")
        except Exception:
            pass

    logger.info("Auto-resolving STRING species taxon ID for '%s' (tax_id=%s)...", species_name, tax_id)

    probe_genes = ["gyrA", "rpoB", "dnaA", "recA", "atpD"]
    query_str = "%0d".join(probe_genes)
    url = f"https://string-db.org/api/json/get_string_ids?identifiers={query_str}&species={tax_id}&limit=1"
    req = urllib.request.Request(url, headers={"User-Agent": "PangenomeDrugPipeline/1.0"})

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8", errors="replace"))
        if data and isinstance(data, list) and len(data) > 0:
            first_hit = data[0]
            resolved_tax_id = int(first_hit.get("ncbiTaxonId", tax_id))
            taxon_name = first_hit.get("taxonName", species_name)
            logger.info("STRING resolved taxon ID: %d (%s)", resolved_tax_id, taxon_name)

            cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_file, "w") as f:
                json.dump({
                    "original_tax_id": tax_id,
                    "species_name": species_name,
                    "resolved_taxon_id": resolved_tax_id,
                    "taxon_name": taxon_name,
                    "resolved_at": time.strftime("%Y-%m-%d %H:%M:%S")
                }, f, indent=2)

            return resolved_tax_id, taxon_name
    except Exception as e:
        logger.warning("STRING API resolution query encountered an issue: %s. Using taxonomy_id %s directly.", e, tax_id)

    return int(tax_id), species_name


def load_string_protein_info(info_gz_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parses STRING protein.info.v12.0.txt.gz to map STRING IDs to gene names and annotations.
    """
    logger.info("Loading STRING protein information from %s...", info_gz_path.name)
    protein_info = {}
    with gzip.open(info_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
        col_map = {col.lstrip("#"): idx for idx, col in enumerate(header)}

        id_idx = col_map.get("string_protein_id", 0)
        name_idx = col_map.get("preferred_name", 1)
        size_idx = col_map.get("protein_size", 2)
        desc_idx = col_map.get("annotation", 3)

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(id_idx, name_idx):
                continue
            str_id = parts[id_idx]
            pref_name = parts[name_idx]
            prot_size = int(parts[size_idx]) if len(parts) > size_idx and parts[size_idx].isdigit() else 0
            annotation = parts[desc_idx] if len(parts) > desc_idx else ""

            protein_info[str_id] = {
                "string_id": str_id,
                "preferred_name": pref_name,
                "protein_size": prot_size,
                "annotation": annotation
            }

    logger.info("Indexed %d proteins in STRING proteome info.", len(protein_info))
    return protein_info


def load_string_aliases(aliases_gz_path: Path) -> Dict[str, str]:
    """
    Parses STRING protein.aliases.v12.0.txt.gz to construct an extensive lookup
    mapping UniProt accessions, gene synonyms, locus tags, and RefSeq IDs to STRING IDs.
    """
    logger.info("Loading STRING protein aliases from %s...", aliases_gz_path.name)
    alias_map = {}
    total_aliases = 0
    with gzip.open(aliases_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                str_id = parts[0]
                alias = parts[1].strip()
                if alias:
                    total_aliases += 1
                    alias_lower = alias.lower()
                    if alias not in alias_map:
                        alias_map[alias] = str_id
                    if alias_lower not in alias_map:
                        alias_map[alias_lower] = str_id

    logger.info("Indexed %d aliases across STRING proteome.", total_aliases)
    return alias_map


def build_interactome_graph(links_gz_path: Path, min_score: int = 700) -> Tuple[nx.Graph, int, int]:
    """
    Loads STRING protein.links.v12.0.txt.gz, filters edges by min_score,
    and constructs an undirected NetworkX Graph.
    """
    logger.info("Building organism interactome graph from %s (score >= %d)...", links_gz_path.name, min_score)
    G = nx.Graph()
    raw_edges = 0
    passed_edges = 0

    with gzip.open(links_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            raw_edges += 1
            score = int(parts[2])
            if score >= min_score:
                passed_edges += 1
                u, v = parts[0], parts[1]
                if not G.has_edge(u, v):
                    G.add_edge(u, v, weight=score / 1000.0, score=score)

    logger.info("Filtered %d high-confidence interactions into graph: %d nodes, %d unique undirected edges.",
                passed_edges, G.number_of_nodes(), G.number_of_edges())
    return G, raw_edges, passed_edges


def compute_network_centralities(G: nx.Graph) -> Dict[str, Dict[str, float]]:
    """
    Computes whole-interactome centrality metrics:
    - Degree Centrality (normalized and raw degree)
    - Betweenness Centrality (Brandes' shortest paths algorithm)
    - Closeness Centrality
    - Clustering Coefficient
    """
    total_nodes = G.number_of_nodes()
    logger.info("Computing Degree Centrality across %d nodes...", total_nodes)
    deg_dict = dict(G.degree())
    deg_norm_dict = nx.degree_centrality(G)

    logger.info("Computing exact Betweenness Centrality (Brandes algorithm) across %d nodes...", total_nodes)
    t0 = time.time()
    betweenness_dict = nx.betweenness_centrality(G, normalized=True)
    logger.info("Betweenness Centrality completed in %.2f seconds.", time.time() - t0)

    logger.info("Computing Closeness Centrality across %d nodes...", total_nodes)
    t0 = time.time()
    closeness_dict = nx.closeness_centrality(G)
    logger.info("Closeness Centrality completed in %.2f seconds.", time.time() - t0)

    logger.info("Computing Clustering Coefficients across %d nodes...", total_nodes)
    clustering_dict = nx.clustering(G)

    node_metrics = {}
    for node in G.nodes():
        node_metrics[node] = {
            "degree": deg_dict.get(node, 0),
            "degree_centrality": deg_norm_dict.get(node, 0.0),
            "betweenness_centrality": betweenness_dict.get(node, 0.0),
            "closeness_centrality": closeness_dict.get(node, 0.0),
            "clustering_coefficient": clustering_dict.get(node, 0.0)
        }

    return node_metrics


def build_gene_to_string_id_map(protein_info: Dict[str, Dict[str, Any]], alias_map: Dict[str, str]) -> Dict[str, str]:
    """
    Builds comprehensive mapping from gene name / locus tag / synonym to STRING protein ID.
    Includes bacterial LPS / cell wall standard synonym conventions:
      - rfa* <-> waa*
      - kdsC <-> yrbI
    """
    combined_map = dict(alias_map)

    for str_id, info in protein_info.items():
        pref = info["preferred_name"].strip()
        if pref:
            combined_map[pref.lower()] = str_id
            combined_map[pref] = str_id

        if "." in str_id:
            locus = str_id.split(".", 1)[1]
            combined_map[locus.lower()] = str_id
            combined_map[locus] = str_id

    # Bacterial nomenclature standard synonyms
    if "yrbi" in combined_map:
        combined_map["kdsc"] = combined_map["yrbi"]
        combined_map["kdsC"] = combined_map["yrbi"]

    # rfa -> waa synonyms
    rfa_waa_pairs = [("rfaC", "waaC"), ("rfaF", "waaF"), ("rfaQ", "waaQ"),
                     ("rfaP", "waaP"), ("rfaY", "waaY"), ("rfaL", "waaL"),
                     ("rfaI", "waaI"), ("rfaJ", "waaJ")]
    for rfa, waa in rfa_waa_pairs:
        if waa.lower() in combined_map:
            combined_map[rfa] = combined_map[waa.lower()]
            combined_map[rfa.lower()] = combined_map[waa.lower()]
        elif rfa.lower() in combined_map:
            combined_map[waa] = combined_map[rfa.lower()]
            combined_map[waa.lower()] = combined_map[rfa.lower()]

    return combined_map


def classify_topological_class(
    degree: float,
    betweenness: float,
    deg_cutoff: float,
    bet_cutoff: float,
    is_mapped: bool
) -> Tuple[bool, bool, str]:
    """
    Classifies a node into one of four topological quadrants per Yu et al. Science (2007):
      - Hub-Bottleneck (HB): High Degree AND High Betweenness
      - Non-Hub Bottleneck (NH-B): Low Degree AND High Betweenness
      - Hub Non-Bottleneck (H-NB): High Degree AND Low Betweenness
      - Non-Hub Non-Bottleneck (NH-NB): Low Degree AND Low Betweenness
    """
    if not is_mapped:
        return False, False, "Unmapped / Peripheral"

    is_hub = degree >= deg_cutoff
    is_bottleneck = betweenness >= bet_cutoff

    if is_hub and is_bottleneck:
        return True, True, "Hub-Bottleneck (HB)"
    elif (not is_hub) and is_bottleneck:
        return False, True, "Non-Hub Bottleneck (NH-B)"
    elif is_hub and (not is_bottleneck):
        return True, False, "Hub Non-Bottleneck (H-NB)"
    else:
        return False, False, "Non-Hub Non-Bottleneck (NH-NB)"


def main():
    parser = argparse.ArgumentParser(description="Step 9: Protein-Protein Interactome & Network Hub-Bottleneck Analysis")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    job_name = config.get("job_name", "klebsiella_pneumoniae")
    results_dir = Path(config.get("results_dir", "results")) / job_name
    step9_cfg = config.get("step9_network_biology_string", {})

    if not step9_cfg.get("enabled", True):
        logger.info("Step 9 is disabled in config.yaml. Skipping.")
        return

    confidence_score = int(step9_cfg.get("confidence_score", 700))
    hub_std_mult = float(step9_cfg.get("hub_std_multiplier", 1.0))
    bet_std_mult = float(step9_cfg.get("bottleneck_std_multiplier", 1.0))
    cache_dir = Path(step9_cfg.get("cache_dir", "databases/cache/string_network"))
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Output directories
    out_dir = results_dir / "09_network_biology"
    tables_dir = out_dir / "tables"
    fasta_dir = out_dir / "fasta"
    summary_dir = out_dir / "summary"

    for d in (tables_dir, fasta_dir, summary_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Input file from Step 8
    step8_master_path = results_dir / "08_metabolic_chokepoint" / "tables" / "metabolic_master.tsv"
    if not step8_master_path.exists():
        raise FileNotFoundError(f"Step 8 master table not found: {step8_master_path}. Please execute Step 8 first.")

    logger.info("================================================================================")
    logger.info("Step 9: Protein-Protein Interactome & Network Hub-Bottleneck Analysis (STRING v12)")
    logger.info("Organism Job: %s", job_name)
    logger.info("Confidence Score Cutoff: >= %d", confidence_score)
    logger.info("================================================================================")

    # 1. Resolve STRING taxon ID dynamically
    resolved_tax_id, taxon_name = resolve_string_species_taxon_id(config, cache_dir)

    # 2. Download and cache STRING whole-organism datasets (info, aliases, links)
    info_url = f"https://stringdb-downloads.org/download/protein.info.v12.0/{resolved_tax_id}.protein.info.v12.0.txt.gz"
    aliases_url = f"https://stringdb-downloads.org/download/protein.aliases.v12.0/{resolved_tax_id}.protein.aliases.v12.0.txt.gz"
    links_url = f"https://stringdb-downloads.org/download/protein.links.v12.0/{resolved_tax_id}.protein.links.v12.0.txt.gz"

    info_path = cache_dir / f"{resolved_tax_id}.protein.info.v12.0.txt.gz"
    aliases_path = cache_dir / f"{resolved_tax_id}.protein.aliases.v12.0.txt.gz"
    links_path = cache_dir / f"{resolved_tax_id}.protein.links.v12.0.txt.gz"

    download_file_cached(info_url, info_path, f"STRING protein info for {taxon_name}")
    download_file_cached(aliases_url, aliases_path, f"STRING protein aliases for {taxon_name}")
    download_file_cached(links_url, links_path, f"STRING protein links for {taxon_name}")

    # 3. Load protein info, aliases, and construct interactome graph
    protein_info = load_string_protein_info(info_path)
    alias_map = load_string_aliases(aliases_path)
    gene_to_string_id = build_gene_to_string_id_map(protein_info, alias_map)

    G, raw_edges_count, passed_edges_count = build_interactome_graph(links_path, min_score=confidence_score)

    # 4. Compute network centrality metrics across all nodes in the proteome
    node_metrics = compute_network_centralities(G)

    # 5. Calculate proteome-wide baseline statistical thresholds
    all_degrees = [m["degree"] for m in node_metrics.values()]
    all_betweenness = [m["betweenness_centrality"] for m in node_metrics.values()]
    all_closeness = [m["closeness_centrality"] for m in node_metrics.values()]

    deg_mean = float(np.mean(all_degrees))
    deg_std = float(np.std(all_degrees))
    deg_p90 = float(np.percentile(all_degrees, 90))

    bet_mean = float(np.mean(all_betweenness))
    bet_std = float(np.std(all_betweenness))
    bet_p90 = float(np.percentile(all_betweenness, 90))

    deg_cutoff = deg_mean + (hub_std_mult * deg_std)
    bet_cutoff = bet_mean + (bet_std_mult * bet_std)

    logger.info("--------------------------------------------------------------------------------")
    logger.info("Proteome Interactome Baseline Statistics (%d proteins, %d interactions):", len(G.nodes()), G.number_of_edges())
    logger.info("  Degree: Mean = %.2f, Std = %.2f, Cutoff (Mean + %.1f*Std) = %.2f (90th Pct = %.2f)",
                deg_mean, deg_std, hub_std_mult, deg_cutoff, deg_p90)
    logger.info("  Betweenness: Mean = %.6f, Std = %.6f, Cutoff (Mean + %.1f*Std) = %.6f (90th Pct = %.6f)",
                bet_mean, bet_std, bet_std_mult, bet_cutoff, bet_p90)
    logger.info("--------------------------------------------------------------------------------")

    # 6. Load candidate targets from Step 8
    df_targets = pd.read_csv(step8_master_path, sep="\t")
    logger.info("Loaded %d prioritized candidate targets from Step 8.", len(df_targets))

    # 7. Map targets to STRING and assign centrality & topological classes
    mapped_count = 0
    results_rows = []

    target_string_ids = set()

    for idx, row in df_targets.iterrows():
        row_dict = dict(row)
        gene_name = str(row.get("Gene_Name", "")).strip() if pd.notna(row.get("Gene_Name")) else ""
        accession = str(row.get("Representative_Accession", "")).strip()
        query_id = str(row.get("Query_ID", "")).strip()
        description = str(row.get("Description", "")).strip()

        # Resolution cascade:
        str_id = None
        for candidate_key in [gene_name, gene_name.lower(), accession, accession.lower(), query_id, query_id.lower()]:
            if candidate_key and candidate_key in gene_to_string_id:
                str_id = gene_to_string_id[candidate_key]
                break

        # If not resolved by direct match, attempt regex extraction from description
        if not str_id and description:
            gn_match = re.search(r'\bGN=([A-Za-z0-9_-]+)', description)
            if gn_match:
                candidate_gn = gn_match.group(1).strip()
                if candidate_gn in gene_to_string_id:
                    str_id = gene_to_string_id[candidate_gn]
                elif candidate_gn.lower() in gene_to_string_id:
                    str_id = gene_to_string_id[candidate_gn.lower()]

        if str_id and str_id in node_metrics:
            mapped_count += 1
            metrics = node_metrics[str_id]
            deg = metrics["degree"]
            deg_cent = metrics["degree_centrality"]
            bet_cent = metrics["betweenness_centrality"]
            clo_cent = metrics["closeness_centrality"]
            clust_coeff = metrics["clustering_coefficient"]

            is_hub, is_bottleneck, topo_class = classify_topological_class(
                deg, bet_cent, deg_cutoff, bet_cutoff, is_mapped=True
            )
            str_annot = protein_info.get(str_id, {}).get("annotation", "")
            target_string_ids.add(str_id)
        else:
            deg = 0
            deg_cent = 0.0
            bet_cent = 0.0
            clo_cent = 0.0
            clust_coeff = 0.0
            is_hub, is_bottleneck, topo_class = classify_topological_class(
                deg, bet_cent, deg_cutoff, bet_cutoff, is_mapped=False
            )
            str_annot = ""

        # Step 9 annotations
        row_dict["STRING_Species_Taxon_ID"] = resolved_tax_id
        row_dict["STRING_Protein_ID"] = str_id if str_id else "Unmapped"
        row_dict["STRING_Mapped"] = bool(str_id is not None)
        row_dict["STRING_Degree"] = deg
        row_dict["Degree_Centrality"] = round(deg_cent, 6)
        row_dict["Betweenness_Centrality"] = round(bet_cent, 8)
        row_dict["Closeness_Centrality"] = round(clo_cent, 6)
        row_dict["Clustering_Coefficient"] = round(clust_coeff, 4)
        row_dict["Is_Network_Hub"] = is_hub
        row_dict["Is_Network_Bottleneck"] = is_bottleneck
        row_dict["Topological_Class"] = topo_class
        row_dict["STRING_Annotation"] = str_annot

        # Multi-Omics Lead Integration: Double Choke-Point + Hub-Bottleneck
        is_double_cp = str(row_dict.get("Is_Double_Choke_Point", "False")).lower() in ["true", "1"]
        is_choke = str(row_dict.get("Is_Choke_Point", "False")).lower() in ["true", "1"]
        is_dual_action = "dual" in str(row_dict.get("Therapeutic_Profile", "")).lower()

        if is_double_cp and is_hub and is_bottleneck:
            if is_dual_action:
                interactome_lead_rank = "Priority 1A: Dual-Action Choke-Point Hub-Bottleneck (Ultimate Vulnerability)"
            else:
                interactome_lead_rank = "Priority 1B: Double Choke-Point Hub-Bottleneck Master Target"
        elif is_double_cp and (is_hub or is_bottleneck):
            interactome_lead_rank = "Priority 2: Double Choke-Point Network Coordinator"
        elif is_double_cp:
            interactome_lead_rank = "Priority 3: Double Choke-Point Modular Enzyme"
        elif is_choke and is_hub and is_bottleneck:
            interactome_lead_rank = "Priority 4: Single Choke-Point Hub-Bottleneck Target"
        elif is_hub and is_bottleneck:
            interactome_lead_rank = "Priority 5: Non-Metabolic Hub-Bottleneck Master Regulator"
        else:
            interactome_lead_rank = "Standard Candidate"

        row_dict["Interactome_Lead_Priority"] = interactome_lead_rank
        results_rows.append(row_dict)

    df_master = pd.DataFrame(results_rows)

    logger.info("Successfully mapped %d of %d targets (%.1f%%) to STRING v12 interactome.",
                mapped_count, len(df_targets), (mapped_count / len(df_targets)) * 100)

    # 8. Extract High-Priority Subsets
    df_hb = df_master[df_master["Topological_Class"] == "Hub-Bottleneck (HB)"]
    df_definitive_leads = df_master[df_master["Interactome_Lead_Priority"].str.startswith("Priority 1")]

    logger.info("Identified %d Hub-Bottleneck (HB) Master Regulators.", len(df_hb))
    logger.info("Identified %d Definitive Interactome Leads (Double Choke-Point + Hub-Bottleneck).", len(df_definitive_leads))

    # 9. Extract Target PPI Subgraph Edges (Target-Target and Target-Neighbor)
    target_edges = []
    for u, v, data in G.edges(data=True):
        if u in target_string_ids or v in target_string_ids:
            score = data.get("score", 0)
            u_name = protein_info.get(u, {}).get("preferred_name", u)
            v_name = protein_info.get(v, {}).get("preferred_name", v)
            u_is_target = u in target_string_ids
            v_is_target = v in target_string_ids
            edge_type = "Target-Target" if (u_is_target and v_is_target) else "Target-InteractionMode"
            target_edges.append({
                "Source_STRING_ID": u,
                "Target_STRING_ID": v,
                "Source_Gene": u_name,
                "Target_Gene": v_name,
                "Interaction_Score": score,
                "Edge_Type": edge_type,
                "Source_Is_Candidate": u_is_target,
                "Target_Is_Candidate": v_is_target
            })

    df_edges = pd.DataFrame(target_edges)
    logger.info("Extracted %d high-confidence interactome edges involving candidate targets.", len(df_edges))

    # 10. Save Output Tables (TSV & CSV)
    master_tsv = tables_dir / "network_centrality_master.tsv"
    master_csv = tables_dir / "network_centrality_master.csv"
    df_master.to_csv(master_tsv, sep="\t", index=False)
    df_master.to_csv(master_csv, index=False)

    hb_tsv = tables_dir / "hub_bottleneck_targets.tsv"
    hb_csv = tables_dir / "hub_bottleneck_targets.csv"
    df_hb.to_csv(hb_tsv, sep="\t", index=False)
    df_hb.to_csv(hb_csv, index=False)

    leads_tsv = tables_dir / "definitive_interactome_leads.tsv"
    leads_csv = tables_dir / "definitive_interactome_leads.csv"
    df_definitive_leads.to_csv(leads_tsv, sep="\t", index=False)
    df_definitive_leads.to_csv(leads_csv, index=False)

    edges_tsv = tables_dir / "target_ppi_edges.tsv"
    edges_csv = tables_dir / "target_ppi_edges.csv"
    df_edges.to_csv(edges_tsv, sep="\t", index=False)
    df_edges.to_csv(edges_csv, index=False)

    # 11. Export FASTA for Hub-Bottleneck Targets
    step4_tier1_fasta = results_dir / "04_subtractive_filter" / "fasta" / "tier1_golden_targets.fasta"
    fasta_seqs = {}
    if step4_tier1_fasta.exists():
        with open(step4_tier1_fasta, "r") as f:
            cur_id = None
            cur_seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id and cur_seq:
                        fasta_seqs[cur_id] = "".join(cur_seq)
                    cur_id = line[1:].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line)
            if cur_id and cur_seq:
                fasta_seqs[cur_id] = "".join(cur_seq)

    hb_fasta_path = fasta_dir / "hub_bottleneck_targets.fasta"
    exported_seqs = 0
    with open(hb_fasta_path, "w") as out_fa:
        for _, r in df_hb.iterrows():
            q_id = str(r["Representative_Accession"])
            if q_id in fasta_seqs:
                g_name = r.get("Gene_Name", "unnamed")
                prio = r.get("Interactome_Lead_Priority", "")
                header = f">{q_id} gene={g_name} topo=Hub-Bottleneck degree={r['STRING_Degree']} betweenness={r['Betweenness_Centrality']} priority={prio}"
                out_fa.write(f"{header}\n{fasta_seqs[q_id]}\n")
                exported_seqs += 1

    logger.info("Exported %d Hub-Bottleneck target sequences to %s.", exported_seqs, hb_fasta_path.name)

    # 11b. Export FASTA for Definitive Interactome Leads (Final Candidates)
    leads_fasta_path = fasta_dir / "definitive_interactome_leads.fasta"
    exported_leads = 0
    with open(leads_fasta_path, "w") as out_leads_fa:
        for _, r in df_definitive_leads.iterrows():
            q_id = str(r["Representative_Accession"])
            if q_id in fasta_seqs:
                g_name = r.get("Gene_Name", "unnamed")
                prio = r.get("Interactome_Lead_Priority", "")
                cp_status = r.get("Choke_Point_Status", "")
                deg = r.get("STRING_Degree", 0)
                bet = r.get("Betweenness_Centrality", 0.0)
                desc = r.get("Description", "")
                header = f">{q_id} GN={g_name} ChokePoint='{cp_status}' Degree={deg} Betweenness={bet:.6f} Priority='{prio}' Desc='{desc}'"
                out_leads_fa.write(f"{header}\n{fasta_seqs[q_id]}\n")
                exported_leads += 1

    logger.info("Exported %d Definitive Interactome Lead (Final Candidate) sequences to %s.", exported_leads, leads_fasta_path.name)

    # 12. Write Step 9 Summary Metrics
    topo_counts = df_master["Topological_Class"].value_counts().to_dict()
    lead_counts = df_master["Interactome_Lead_Priority"].value_counts().to_dict()

    summary_data = {
        "step": "09_network_biology_string",
        "job_name": job_name,
        "resolved_species_taxon_id": resolved_tax_id,
        "taxon_name": taxon_name,
        "string_confidence_score_cutoff": confidence_score,
        "proteome_interactome_statistics": {
            "total_proteins_in_graph": G.number_of_nodes(),
            "total_high_confidence_edges": G.number_of_edges(),
            "degree_mean": round(deg_mean, 2),
            "degree_std": round(deg_std, 2),
            "degree_cutoff_hub": round(deg_cutoff, 2),
            "degree_90th_percentile": round(deg_p90, 2),
            "betweenness_mean": round(bet_mean, 8),
            "betweenness_std": round(bet_std, 8),
            "betweenness_cutoff_bottleneck": round(bet_cutoff, 8),
            "betweenness_90th_percentile": round(bet_p90, 8)
        },
        "target_pool_analysis": {
            "total_targets_evaluated": len(df_master),
            "targets_mapped_to_string": mapped_count,
            "target_mapping_rate_pct": round((mapped_count / len(df_master)) * 100, 2),
            "topological_class_counts": topo_counts,
            "interactome_lead_priority_counts": lead_counts,
            "total_hub_bottleneck_targets": len(df_hb),
            "total_definitive_interactome_leads": len(df_definitive_leads)
        },
        "top_leads_interactome_status": [
            {
                "Gene_Name": r.get("Gene_Name"),
                "Accession": r.get("Representative_Accession"),
                "Topological_Class": r.get("Topological_Class"),
                "Degree": int(r.get("STRING_Degree", 0)),
                "Betweenness_Centrality": float(r.get("Betweenness_Centrality", 0.0)),
                "Choke_Point_Status": r.get("Choke_Point_Status"),
                "Priority": r.get("Interactome_Lead_Priority")
            }
            for _, r in df_definitive_leads.head(10).iterrows()
        ]
    }

    summary_json_path = summary_dir / "step9_network_summary.json"
    with open(summary_json_path, "w") as f:
        json.dump(summary_data, f, indent=2)

    # 13. Create completion sentinel
    done_file = out_dir / ".network_complete"
    with open(done_file, "w") as f:
        f.write(f"Completed Step 9 at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    logger.info("================================================================================")
    logger.info("Step 9 Completed Successfully!")
    logger.info("  Master TSV: %s", master_tsv)
    logger.info("  Hub-Bottleneck TSV: %s (%d targets)", hb_tsv, len(df_hb))
    logger.info("  Definitive Interactome Leads: %s (%d targets)", leads_tsv, len(df_definitive_leads))
    logger.info("  Candidate PPI Edges: %s (%d edges)", edges_tsv, len(df_edges))
    logger.info("  Summary Metrics: %s", summary_json_path)
    logger.info("================================================================================")


if __name__ == "__main__":
    main()
