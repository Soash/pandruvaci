#!/usr/bin/env python3
"""
Step 8: Metabolic Pathway & Choke-Point Analysis (KEGG / BioCyc)
===============================================================
Maps prioritized antibacterial candidates into the complete biochemical network
of the pathogen, identifying unbypassable metabolic bottlenecks (choke-points)
and pathogen-exclusive biochemical pathways absent in humans.

Academic Citations:
- Choke-Point Theory: Yeh I, et al. Genome Res (2004) 14:917–924.
- KEGG Pathways: Kanehisa M, et al. Nucleic Acids Res (2023) 51:D587–D592.
- Rhea Reactions: Bansal P, et al. Nucleic Acids Res (2022) 50:D693–D700.

Usage:
    python scripts/08_metabolic_chokepoint_analysis.py [--config config.yaml]
"""

import os
import re
import sys
import time
import json
import logging
import argparse
import urllib.request
import urllib.parse
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
logger = logging.getLogger("08_metabolic_chokepoint")


DEFAULT_CURRENCY_METABOLITES = {
    "h2o", "water", "atp", "adp", "amp", "nadh", "nad+", "nadph", "nadp+",
    "phosphate", "diphosphate", "orthophosphate", "pyrophosphate", "pi", "ppi",
    "h+", "proton", "co2", "carbon dioxide", "coa", "coenzyme a", "oxygen", "o2",
    "ammonia", "nh3", "nh4+", "na+", "k+", "cl-", "h2o2", "hydrogen peroxide",
    "glutathione", "oxidized glutathione", "reduced glutathione",
    "thioredoxin", "oxidized thioredoxin", "reduced thioredoxin"
}


def download_url_cached(url: str, cache_path: Path, description: str, retries: int = 3, timeout: int = 30) -> str:
    """Fetches URL with local caching and retries."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    if cache_path.exists() and cache_path.stat().st_size > 0:
        logger.debug("Loading cached %s from %s", description, cache_path.name)
        with open(cache_path, "r", encoding="utf-8", errors="replace") as f:
            return f.read()

    logger.info("Downloading %s from %s...", description, url)
    req = urllib.request.Request(url, headers={"User-Agent": "PangenomeDrugPipeline/1.0 (Research)"})
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            with open(cache_path, "w", encoding="utf-8") as out_f:
                out_f.write(text)
            return text
        except Exception as e:
            last_err = e
            logger.warning("Attempt %d/%d failed for %s: %s", attempt, retries, description, e)
            time.sleep(attempt * 2)

    raise RuntimeError(f"Failed to fetch {description} from {url}: {last_err}")


def resolve_kegg_organism_code(cfg: Dict[str, Any], cache_dir: Path) -> str:
    """
    Dynamically resolves KEGG organism code from config or KEGG REST API.
    Handles any bacterial pathogen dynamically without hardcoding.
    """
    step8_cfg = cfg.get("step8_metabolic_chokepoint", {})
    code = step8_cfg.get("kegg_organism_code", "auto")
    if code and code != "auto":
        logger.info("Using configured KEGG organism code: '%s'", code)
        return code

    tax_id = cfg.get("taxonomy_id") or cfg.get("query_params", {}).get("taxonomy_id", "")
    job_name = cfg.get("job_name", "klebsiella_pneumoniae")
    # Search by full species name (e.g. "klebsiella pneumoniae" or "pseudomonas aeruginosa")
    species_query = job_name.replace("_", " ").strip()
    cache_path = cache_dir / f"kegg_genome_{job_name}.txt"
    logger.info("Auto-resolving KEGG organism code for species='%s' (tax_id=%s)...", species_query, tax_id)

    url = f"http://rest.kegg.jp/find/genome/{urllib.parse.quote(species_query)}"
    try:
        raw = download_url_cached(url, cache_path, f"KEGG genome search ({species_query})")
        lines = [l.strip() for l in raw.strip().split("\n") if l.strip()]
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 2:
                org_part = parts[1].split(";")[0].strip()
                if org_part and len(org_part) <= 5 and org_part.isalnum():
                    logger.info("Dynamically matched KEGG organism code: '%s' (%s)", org_part, parts[1])
                    return org_part
    except Exception as e:
        logger.warning("Could not auto-resolve via species query (%s), attempting tax_id fallback...", e)

    # Fallback to searching by organism name
    clean_name = job_name.replace("_", " ").split()[0]
    fallback_cache = cache_dir / f"kegg_genome_{clean_name}.txt"
    try:
        url_fb = f"http://rest.kegg.jp/find/genome/{urllib.parse.quote(clean_name)}"
        raw_fb = download_url_cached(url_fb, fallback_cache, f"KEGG genome fallback ({clean_name})")
        lines = [l.strip() for l in raw_fb.strip().split("\n") if l.strip()]
        if lines:
            parts = lines[0].split("\t")
            if len(parts) >= 2:
                org_part = parts[1].split(";")[0].strip()
                logger.info("Fallback matched KEGG organism code: '%s'", org_part)
                return org_part
    except Exception as e:
        logger.error("Failed KEGG fallback resolution: %s", e)

    default_code = "kpn" if "klebsiella" in job_name.lower() else "eco"
    logger.warning("Falling back to default KEGG organism code: '%s'", default_code)
    return default_code


def fetch_pathway_catalog(kegg_org: str, human_code: str, cache_dir: Path) -> Tuple[Dict[str, str], Set[str], Set[str]]:
    """
    Downloads organism pathways and human pathways, and identifies pathogen-exclusive pathways.
    Returns: (organism_pathways {id: name}, pathogen_exclusive_ids, all_human_path_ids)
    """
    org_cache = cache_dir / f"pathways_{kegg_org}.txt"
    hsa_cache = cache_dir / f"pathways_{human_code}.txt"

    raw_org = download_url_cached(f"http://rest.kegg.jp/list/pathway/{kegg_org}", org_cache, f"{kegg_org} pathways")
    raw_hsa = download_url_cached(f"http://rest.kegg.jp/list/pathway/{human_code}", hsa_cache, f"{human_code} pathways")

    org_pathways: Dict[str, str] = {}
    for line in raw_org.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            num_match = re.search(r"\d{5}", parts[0])
            if num_match:
                num_id = num_match.group()
                pname = parts[1].split(" - ")[0].strip()
                org_pathways[num_id] = pname

    human_nums: Set[str] = set()
    for line in raw_hsa.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 1:
            num_match = re.search(r"\d{5}", parts[0])
            if num_match:
                human_nums.add(num_match.group())

    exclusive_ids = {pid for pid in org_pathways if pid not in human_nums}
    logger.info("Pathogen Pathways: %d | Human Pathways: %d | Pathogen-Exclusive Pathways: %d",
                len(org_pathways), len(human_nums), len(exclusive_ids))

    return org_pathways, exclusive_ids, human_nums


def fetch_kegg_reactions_and_enzymes(kegg_org: str, cache_dir: Path) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Downloads and parses KEGG reaction definitions, enzyme-to-reaction mappings, and gene-to-pathway mappings.
    """
    rxn_cache = cache_dir / "kegg_reaction_list.txt"
    ec_rxn_cache = cache_dir / "kegg_ec_to_reaction.txt"
    gene_path_cache = cache_dir / f"kegg_{kegg_org}_gene_pathway.txt"
    gene_ec_cache = cache_dir / f"kegg_{kegg_org}_gene_enzyme.txt"

    raw_rxn = download_url_cached("http://rest.kegg.jp/list/reaction", rxn_cache, "KEGG reaction list")
    raw_ec_rxn = download_url_cached("http://rest.kegg.jp/link/reaction/enzyme", ec_rxn_cache, "EC-reaction mapping")
    raw_gene_path = download_url_cached(f"http://rest.kegg.jp/link/pathway/{kegg_org}", gene_path_cache, f"{kegg_org} gene pathways")
    raw_gene_ec = download_url_cached(f"http://rest.kegg.jp/link/enzyme/{kegg_org}", gene_ec_cache, f"{kegg_org} gene enzymes")

    reactions: Dict[str, Dict[str, Any]] = {}
    for line in raw_rxn.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            rid = parts[0].strip()
            rest = parts[1].strip()
            eq = rest.split(";")[-1].strip()
            if "<=>" in eq:
                left, right = eq.split("<=>", 1)
            elif "=>" in eq:
                left, right = eq.split("=>", 1)
            else:
                continue

            def parse_metabolites(text_side: str) -> List[str]:
                mets = []
                for item in text_side.split(" + "):
                    m = re.sub(r"^\d+\s+", "", item.strip())
                    m = re.sub(r"^\([a-z0-9+-]+\)\s*", "", m)
                    m = m.strip()
                    if m:
                        mets.append(m)
                return mets

            reactions[rid] = {
                "definition": rest,
                "substrates": parse_metabolites(left),
                "products": parse_metabolites(right),
            }

    ec_to_rxn: Dict[str, List[str]] = {}
    for line in raw_ec_rxn.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            ec_id = parts[0].replace("ec:", "").strip()
            rxn_id = parts[1].replace("rn:", "").strip()
            ec_to_rxn.setdefault(ec_id, []).append(rxn_id)

    gene_to_paths: Dict[str, List[str]] = {}
    for line in raw_gene_path.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            gid = parts[0].replace(f"{kegg_org}:", "").strip().lower()
            num_match = re.search(r"\d{5}", parts[1])
            if num_match:
                gene_to_paths.setdefault(gid, []).append(num_match.group())

    for line in raw_gene_ec.strip().split("\n"):
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            gid = parts[0].replace(f"{kegg_org}:", "").strip().lower()
            ec_id = parts[1].replace("ec:", "").strip()
            for rxn_id in ec_to_rxn.get(ec_id, []):
                gene_to_paths.setdefault(f"ec_{ec_id}", []).append(rxn_id)

    logger.info("Loaded %d KEGG reactions, %d EC-reaction links.", len(reactions), len(ec_to_rxn))
    return reactions, ec_to_rxn, gene_to_paths


def query_uniprot_annotations(accessions: List[str], cache_dir: Path, batch_size: int = 100) -> Dict[str, Dict[str, Any]]:
    """
    Fetches curated catalytic activity, EC numbers, and pathways from UniProtKB REST API.
    """
    up_cache_dir = cache_dir / "uniprot_metabolic"
    up_cache_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, Dict[str, Any]] = {}
    missing: List[str] = []

    for acc in accessions:
        cp = up_cache_dir / f"{acc}.json"
        if cp.exists():
            try:
                with open(cp, "r", encoding="utf-8") as f:
                    results[acc] = json.load(f)
            except Exception:
                missing.append(acc)
        else:
            missing.append(acc)

    if not missing:
        logger.info("All %d UniProt metabolic records loaded from cache.", len(accessions))
        return results

    logger.info("Querying UniProtKB REST API for %d accessions (batch size: %d)...", len(missing), batch_size)
    for i in range(0, len(missing), batch_size):
        chunk = missing[i : i + batch_size]
        acc_str = ",".join(chunk)
        url = (
            f"https://rest.uniprot.org/uniprotkb/accessions?"
            f"accessions={urllib.parse.quote(acc_str)}"
            f"&fields=accession,id,gene_names,protein_name,cc_catalytic_activity,cc_pathway,xref_kegg"
            f"&format=tsv"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "PangenomeDrugPipeline/1.0"})
        tsv_text = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    tsv_text = resp.read().decode("utf-8")
                break
            except Exception as e:
                time.sleep(1.5 * (attempt + 1))
                if attempt == 2:
                    logger.warning("Failed UniProt chunk %d: %s", (i // batch_size) + 1, e)

        if not tsv_text:
            continue

        lines = tsv_text.strip().split("\n")
        if len(lines) <= 1:
            continue

        header = lines[0].split("\t")
        col_map = {col.strip(): idx for idx, col in enumerate(header)}

        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= col_map.get("Entry", 0):
                continue

            acc = parts[col_map["Entry"]].strip()
            prot_name = parts[col_map.get("Protein names", 1)] if "Protein names" in col_map and len(parts) > col_map["Protein names"] else ""
            cat_act = parts[col_map.get("Catalytic activity", 2)] if "Catalytic activity" in col_map and len(parts) > col_map["Catalytic activity"] else ""
            pathway = parts[col_map.get("Pathway", 3)] if "Pathway" in col_map and len(parts) > col_map["Pathway"] else ""
            kegg_xref = parts[col_map.get("KEGG", 4)] if "KEGG" in col_map and len(parts) > col_map["KEGG"] else ""

            ec_matches = set(re.findall(r"EC[= ](\d+\.\d+\.\d+\.[\d-]+)", cat_act + " " + prot_name))
            rhea_matches = set(re.findall(r"Rhea:RHEA:(\d+)", cat_act))

            rxn_strings = []
            if "Reaction=" in cat_act:
                for match in re.finditer(r"Reaction=([^;]+);", cat_act):
                    rxn_strings.append(match.group(1).strip())

            record = {
                "accession": acc,
                "protein_name": prot_name,
                "ec_numbers": list(ec_matches),
                "rhea_ids": list(rhea_matches),
                "reaction_strings": rxn_strings,
                "uniprot_pathway": pathway,
                "kegg_xref": kegg_xref,
            }

            results[acc] = record
            with open(up_cache_dir / f"{acc}.json", "w", encoding="utf-8") as f:
                json.dump(record, f)

    return results


def is_currency_metabolite(met: str, currency_set: Set[str]) -> bool:
    """Checks if a metabolite matches currency molecule list."""
    clean = met.strip().lower()
    if clean in currency_set:
        return True
    for c in currency_set:
        if len(c) > 2 and (clean == c or clean.startswith(c + " ") or clean.endswith(" " + c)):
            return True
    return False


def build_metabolic_bipartite_network(
    candidate_records: Dict[str, Dict[str, Any]],
    kegg_reactions: Dict[str, Dict[str, Any]],
    ec_to_rxn: Dict[str, List[str]],
    currency_set: Set[str],
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, List[Dict[str, Any]]]]:
    """
    Constructs the metabolic bipartite network across all active reactions in the organism.
    """
    consuming_reactions: Dict[str, Set[str]] = {}
    producing_reactions: Dict[str, Set[str]] = {}
    target_reaction_map: Dict[str, List[Dict[str, Any]]] = {}

    for acc, rec in candidate_records.items():
        active_rxns: List[Dict[str, Any]] = []

        # 1. Parse UniProt explicit reaction strings
        for rxn_str in rec.get("reaction_strings", []):
            if "=" in rxn_str:
                parts = rxn_str.split("=")
                left, right = parts[0], parts[1]
                subs = [s.strip() for s in left.split(" + ") if s.strip()]
                prods = [p.strip() for p in right.split(" + ") if p.strip()]
                active_rxns.append({
                    "source": "UniProt",
                    "id": f"RHEA_{acc}",
                    "substrates": subs,
                    "products": prods,
                })

        # 2. Add mapped KEGG reactions via EC numbers
        for ec in rec.get("ec_numbers", []):
            for rxn_id in ec_to_rxn.get(ec, []):
                if rxn_id in kegg_reactions:
                    k_rxn = kegg_reactions[rxn_id]
                    active_rxns.append({
                        "source": "KEGG",
                        "id": rxn_id,
                        "ec": ec,
                        "substrates": k_rxn["substrates"],
                        "products": k_rxn["products"],
                    })

        target_reaction_map[acc] = active_rxns

        # Add to global bipartite reaction graph
        for rxn in active_rxns:
            rxn_key = f"{acc}_{rxn['id']}"
            for s in rxn["substrates"]:
                if not is_currency_metabolite(s, currency_set):
                    s_clean = s.strip().lower()
                    consuming_reactions.setdefault(s_clean, set()).add(rxn_key)

            for p in rxn["products"]:
                if not is_currency_metabolite(p, currency_set):
                    p_clean = p.strip().lower()
                    producing_reactions.setdefault(p_clean, set()).add(rxn_key)

    # Also incorporate global KEGG organism background reactions into network
    for rxn_id, k_rxn in kegg_reactions.items():
        rxn_key = f"KEGG_{rxn_id}"
        for s in k_rxn["substrates"]:
            if not is_currency_metabolite(s, currency_set):
                s_clean = s.strip().lower()
                consuming_reactions.setdefault(s_clean, set()).add(rxn_key)
        for p in k_rxn["products"]:
            if not is_currency_metabolite(p, currency_set):
                p_clean = p.strip().lower()
                producing_reactions.setdefault(p_clean, set()).add(rxn_key)

    return consuming_reactions, producing_reactions, target_reaction_map


def evaluate_chokepoints(
    target_reaction_map: Dict[str, List[Dict[str, Any]]],
    consuming_reactions: Dict[str, Set[str]],
    producing_reactions: Dict[str, Set[str]],
    currency_set: Set[str],
) -> Dict[str, Dict[str, Any]]:
    """
    Executes the Yeh et al. Choke-Point Algorithm for each target.
    """
    results: Dict[str, Dict[str, Any]] = {}

    for acc, rxns in target_reaction_map.items():
        if not rxns:
            results[acc] = {
                "chokepoint_status": "Non-Enzymatic / Structural Target",
                "is_chokepoint": False,
                "is_double_chokepoint": False,
                "unique_consumed_substrates": [],
                "unique_produced_products": [],
                "reaction_count": 0,
            }
            continue

        unique_consumed: Set[str] = set()
        unique_produced: Set[str] = set()

        for rxn in rxns:
            for s in rxn["substrates"]:
                if not is_currency_metabolite(s, currency_set):
                    s_clean = s.strip().lower()
                    consumers = consuming_reactions.get(s_clean, set())
                    target_rxn_keys = {k for k in consumers if k.startswith(acc + "_")}
                    if len(consumers) == len(target_rxn_keys) and len(consumers) > 0:
                        unique_consumed.add(s)

            for p in rxn["products"]:
                if not is_currency_metabolite(p, currency_set):
                    p_clean = p.strip().lower()
                    producers = producing_reactions.get(p_clean, set())
                    target_rxn_keys = {k for k in producers if k.startswith(acc + "_")}
                    if len(producers) == len(target_rxn_keys) and len(producers) > 0:
                        unique_produced.add(p)

        is_consumer_choke = len(unique_consumed) > 0
        is_producer_choke = len(unique_produced) > 0

        if is_consumer_choke and is_producer_choke:
            status = "Double Choke-Point (Bactericidal + Starvation)"
            is_choke = True
            is_double = True
        elif is_consumer_choke:
            status = "Single Choke-Point (Unique Substrate Consumer)"
            is_choke = True
            is_double = False
        elif is_producer_choke:
            status = "Single Choke-Point (Unique Product Synthesizer)"
            is_choke = True
            is_double = False
        else:
            status = "Non-Choke-Point (Metabolic Bypass Available)"
            is_choke = False
            is_double = False

        results[acc] = {
            "chokepoint_status": status,
            "is_chokepoint": is_choke,
            "is_double_chokepoint": is_double,
            "unique_consumed_substrates": sorted(list(unique_consumed)),
            "unique_produced_products": sorted(list(unique_produced)),
            "reaction_count": len(rxns),
        }

    return results


def run_step8(config_path: str = "config.yaml"):
    """Main execution function for Step 8: Metabolic Pathway & Choke-Point Analysis."""
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    job_name = cfg.get("job_name", "klebsiella_pneumoniae")
    results_dir = Path(cfg.get("results_dir", "results"))
    threads = resolve_threads(cfg.get("hardware", {}))
    step8_cfg = cfg.get("step8_metabolic_chokepoint", {})

    logger.info("=" * 80)
    logger.info("Step 8: Metabolic Pathway & Choke-Point Analysis (KEGG / BioCyc)")
    logger.info("Job Name: %s | Allocated Multithreading: %d cores", job_name, threads)
    logger.info("=" * 80)

    step7_master = results_dir / job_name / "07_target_druggability_ttd" / "tables" / "ttd_druggability_master.tsv"
    if not step7_master.exists():
        raise FileNotFoundError(f"Step 7 master table not found: {step7_master}")

    df_master = pd.read_csv(step7_master, sep="\t", low_memory=False)
    logger.info("Loaded %d prioritized targets from Step 7 master table.", len(df_master))

    step8_dir = results_dir / job_name / "08_metabolic_chokepoint"
    tables_dir = step8_dir / "tables"
    fasta_dir = step8_dir / "fasta"
    summary_dir = step8_dir / "summary"
    cache_dir = Path(step8_cfg.get("cache_dir", "databases/cache/kegg_metabolic"))
    if not cache_dir.is_absolute():
        cache_dir = Path.cwd() / cache_dir

    for d in [tables_dir, fasta_dir, summary_dir, cache_dir]:
        d.mkdir(parents=True, exist_ok=True)

    user_currency = step8_cfg.get("currency_metabolites", [])
    currency_set = set(DEFAULT_CURRENCY_METABOLITES)
    if user_currency:
        currency_set.update(c.strip().lower() for c in user_currency)
    filter_currency = step8_cfg.get("filter_currency_metabolites", True)

    kegg_org = resolve_kegg_organism_code(cfg, cache_dir)
    human_code = step8_cfg.get("human_kegg_code", "hsa")

    org_pathways, exclusive_pathways, _ = fetch_pathway_catalog(kegg_org, human_code, cache_dir)
    kegg_rxns, ec_to_rxn, gene_to_paths = fetch_kegg_reactions_and_enzymes(kegg_org, cache_dir)

    accessions = df_master["Representative_Accession"].dropna().unique().tolist()
    uniprot_records = query_uniprot_annotations(accessions, cache_dir)

    for _, row in df_master.iterrows():
        acc = str(row["Representative_Accession"]).strip()
        ttd_ec = str(row.get("TTD_EC_Number", "")).strip()
        if ttd_ec and ttd_ec.lower() != "nan":
            ec_clean = ttd_ec.replace("EC ", "").strip()
            if acc in uniprot_records:
                if ec_clean not in uniprot_records[acc]["ec_numbers"]:
                    uniprot_records[acc]["ec_numbers"].append(ec_clean)

    consuming_rxns, producing_rxns, target_rxns = build_metabolic_bipartite_network(
        uniprot_records,
        kegg_rxns,
        ec_to_rxn,
        currency_set if filter_currency else set(),
    )

    chokepoint_results = evaluate_chokepoints(
        target_rxns,
        consuming_rxns,
        producing_rxns,
        currency_set if filter_currency else set(),
    )

    master_rows = []
    for _, row in df_master.iterrows():
        acc = str(row["Representative_Accession"]).strip()
        gname = str(row.get("Gene_Name", "")).strip().lower()
        up_rec = uniprot_records.get(acc, {})
        cp_rec = chokepoint_results.get(acc, {
            "chokepoint_status": "Non-Enzymatic / Structural Target",
            "is_chokepoint": False,
            "is_double_chokepoint": False,
            "unique_consumed_substrates": [],
            "unique_produced_products": [],
            "reaction_count": 0,
        })

        associated_paths: Set[str] = set()
        if up_rec.get("uniprot_pathway"):
            associated_paths.add(up_rec["uniprot_pathway"])
        kegg_pids = gene_to_paths.get(gname, [])
        for pid in kegg_pids:
            if pid in org_pathways:
                associated_paths.add(f"KEGG:{pid} {org_pathways[pid]}")

        is_pathogen_exclusive = any(
            any(ex in p for ex in exclusive_pathways)
            for p in associated_paths
        )
        exclusive_keywords = ["peptidoglycan", "lipid a", "lipopolysaccharide", "shikimate", "d-alanine", "o-antigen", "teichoic"]
        for kw in exclusive_keywords:
            if any(kw in p.lower() for p in associated_paths) or any(kw in str(row.get("Description", "")).lower() for p in [1]):
                is_pathogen_exclusive = True
                break

        ecs = up_rec.get("ec_numbers", [])
        ec_str = "; ".join(ecs) if ecs else ""

        rxn_count = cp_rec["reaction_count"]
        choke_status = cp_rec["chokepoint_status"]
        consumed_str = "; ".join(cp_rec["unique_consumed_substrates"])
        produced_str = "; ".join(cp_rec["unique_produced_products"])
        pathway_str = "; ".join(sorted(list(associated_paths))) if associated_paths else "Uncharacterized Pathway"

        is_double = cp_rec["is_double_chokepoint"]
        is_choke = cp_rec["is_chokepoint"]
        is_dual_action = "Dual-Action" in str(row.get("Therapeutic_Profile", ""))
        is_novel = "Novel" in str(row.get("Clinical_Druggability_Category", ""))
        is_approved = "Approved" in str(row.get("Clinical_Druggability_Category", ""))

        if is_double and is_pathogen_exclusive and is_dual_action:
            definitive_tier = "Priority 1: Gold-Standard Dual-Action Double Choke-Point"
        elif is_double and is_pathogen_exclusive:
            definitive_tier = "Priority 2: Pathogen-Exclusive Double Choke-Point"
        elif is_double:
            definitive_tier = "Priority 3: Core Metabolic Double Choke-Point"
        elif is_choke and is_pathogen_exclusive:
            definitive_tier = "Priority 4: Pathogen-Exclusive Single Choke-Point"
        elif is_choke:
            definitive_tier = "Priority 5: Core Metabolic Single Choke-Point"
        elif is_approved:
            definitive_tier = "Clinically Approved Repurposing Benchmark"
        else:
            definitive_tier = "Deprioritized / Non-Choke-Point Enzyme"

        row_dict = row.to_dict()
        row_dict.update({
            "KEGG_Organism_Code": kegg_org,
            "Enzyme_Commission_EC": ec_str,
            "Active_Reactions_Count": rxn_count,
            "Choke_Point_Status": choke_status,
            "Is_Choke_Point": is_choke,
            "Is_Double_Choke_Point": is_double,
            "Unique_Consumed_Substrates": consumed_str,
            "Unique_Produced_Products": produced_str,
            "Associated_Pathways": pathway_str,
            "Pathogen_Exclusive_Pathway": is_pathogen_exclusive,
            "Definitive_Target_Priority": definitive_tier,
        })
        master_rows.append(row_dict)

    df_out = pd.DataFrame(master_rows)

    master_tsv = tables_dir / "metabolic_master.tsv"
    master_csv = tables_dir / "metabolic_master.csv"
    df_out.to_csv(master_tsv, sep="\t", index=False)
    df_out.to_csv(master_csv, index=False)
    logger.info("Saved Step 8 Master Table: %s (%d targets)", master_tsv.name, len(df_out))

    df_double = df_out[df_out["Is_Double_Choke_Point"] == True]
    df_double.to_csv(tables_dir / "double_chokepoint_targets.tsv", sep="\t", index=False)
    df_double.to_csv(tables_dir / "double_chokepoint_targets.csv", index=False)

    df_pathogen_excl = df_out[df_out["Pathogen_Exclusive_Pathway"] == True]
    df_pathogen_excl.to_csv(tables_dir / "pathogen_exclusive_pathway_targets.tsv", sep="\t", index=False)
    df_pathogen_excl.to_csv(tables_dir / "pathogen_exclusive_pathway_targets.csv", index=False)

    df_gold = df_out[df_out["Definitive_Target_Priority"].str.startswith("Priority 1") | df_out["Definitive_Target_Priority"].str.startswith("Priority 2")]
    df_gold.to_csv(tables_dir / "definitive_gold_standard_targets.tsv", sep="\t", index=False)
    df_gold.to_csv(tables_dir / "definitive_gold_standard_targets.csv", index=False)

    query_fasta = results_dir / job_name / "05_subcellular_localization" / "fasta" / "prioritized_drug_candidates.fasta"
    fasta_seqs = {}
    if query_fasta.exists():
        with open(query_fasta, "r") as f:
            cur_hdr = ""
            cur_seq = []
            for line in f:
                if line.startswith(">"):
                    if cur_hdr:
                        fasta_seqs[cur_hdr] = "".join(cur_seq)
                    cur_hdr = line.strip()[1:].split()[0]
                    cur_seq = []
                else:
                    cur_seq.append(line.strip())
            if cur_hdr:
                fasta_seqs[cur_hdr] = "".join(cur_seq)

    def export_cohort_fasta(df_cohort: pd.DataFrame, out_path: Path):
        count = 0
        with open(out_path, "w", encoding="utf-8") as f:
            for _, r in df_cohort.iterrows():
                qid = str(r["Representative_Accession"]).strip()
                if qid in fasta_seqs:
                    f.write(f">{qid} GN={r.get('Gene_Name', '')} ChokePoint={r.get('Choke_Point_Status', '')} Priority={r.get('Definitive_Target_Priority', '')}\n")
                    f.write(f"{fasta_seqs[qid]}\n")
                    count += 1
        logger.info("Exported %d sequences to %s.", count, out_path.name)

    export_cohort_fasta(df_double, fasta_dir / "double_chokepoint_targets.fasta")
    export_cohort_fasta(df_pathogen_excl, fasta_dir / "pathogen_exclusive_pathway_targets.fasta")
    export_cohort_fasta(df_gold, fasta_dir / "definitive_gold_standard_targets.fasta")

    choke_breakdown = df_out["Choke_Point_Status"].value_counts().to_dict()
    priority_breakdown = df_out["Definitive_Target_Priority"].value_counts().to_dict()

    summary_data = {
        "job_name": job_name,
        "kegg_organism_code": kegg_org,
        "human_kegg_code": human_code,
        "total_targets_evaluated": len(df_out),
        "double_chokepoints_count": int(df_out["Is_Double_Choke_Point"].sum()),
        "total_chokepoints_count": int(df_out["Is_Choke_Point"].sum()),
        "pathogen_exclusive_targets_count": int(df_out["Pathogen_Exclusive_Pathway"].sum()),
        "gold_standard_candidates_count": len(df_gold),
        "chokepoint_breakdown": choke_breakdown,
        "definitive_priority_breakdown": priority_breakdown,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    summary_json = summary_dir / "step8_metabolic_summary.json"
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2)

    flag = step8_dir / ".metabolic_complete"
    flag.touch()

    logger.info("=" * 80)
    logger.info("Step 8 Completed Successfully!")
    logger.info("Total Evaluated Targets:          %d", len(df_out))
    logger.info("Double Choke-Point Targets:       %d", summary_data["double_chokepoints_count"])
    logger.info("Total Choke-Point Bottlenecks:    %d", summary_data["total_chokepoints_count"])
    logger.info("Pathogen-Exclusive Targets:       %d", summary_data["pathogen_exclusive_targets_count"])
    logger.info("Definitive Gold-Standard Leads:   %d", summary_data["gold_standard_candidates_count"])
    logger.info("Summary JSON: %s", summary_json)
    logger.info("=" * 80)


def main():
    parser = argparse.ArgumentParser(description="Step 8: Metabolic Pathway & Choke-Point Analysis")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    args = parser.parse_args()
    run_step8(args.config)


if __name__ == "__main__":
    main()
