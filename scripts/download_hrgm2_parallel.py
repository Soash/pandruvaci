#!/usr/bin/env python3
"""
Parallel Multi-part Downloader for HRGM2 Gut Microbiome Catalogue
================================================================
Downloads HRGMv2_50_Proteins_rep_seq.faa.gz (~652 MB) in parallel
byte-range chunks, reassembles them, and decompresses into FASTA format.
"""

import os
import sys
import time
import gzip
import shutil
import urllib.request
from pathlib import Path
from typing import Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

URL = "https://www.decodebiome.org/HRGM/data/protein_catalog/0.HRGMv2_Proteins/6.HRGMv2_50_Proteins/HRGMv2_50_Proteins_rep_seq.faa.gz"


def download_chunk(url: str, start: int, end: int, part_file: Path, retries: int = 5) -> Tuple[int, int]:
    expected_size = end - start + 1
    if part_file.exists() and part_file.stat().st_size == expected_size:
        return part_file.stat().st_size

    tmp_file = part_file.with_suffix(".tmp")
    headers = {
        "User-Agent": "Mozilla/5.0 (HRGM2 Downloader; Bioinformatics)",
        "Range": f"bytes={start}-{end}",
    }

    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=35) as resp, open(tmp_file, "wb") as f:
                while True:
                    buf = resp.read(65536)
                    if not buf:
                        break
                    f.write(buf)
            
            if tmp_file.stat().st_size == expected_size:
                tmp_file.rename(part_file)
                return part_file.stat().st_size
        except Exception as e:
            if attempt == retries - 1:
                raise e
            time.sleep(2 * (attempt + 1))
    
    raise RuntimeError(f"Failed to download range {start}-{end}")


def main():
    target_dir = Path("databases/gut_microbiota")
    target_dir.mkdir(parents=True, exist_ok=True)
    out_gz = target_dir / "HRGMv2_50_Proteins_rep_seq.faa.gz"
    out_fa = target_dir / "hrgm2_50_rep.fasta"

    if out_fa.exists() and out_fa.stat().st_size > 500_000_000:
        print(f"HRGM2 FASTA already exists and is ready: {out_fa} ({out_fa.stat().st_size:,} bytes)")
        return

    # Check Content-Length
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        total_len = int(resp.headers.get("Content-Length"))

    print(f"Total HRGM2 archive size: {total_len:,} bytes (~{total_len/1e6:.1f} MB)")

    parts_dir = target_dir / "parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    num_parts = 32
    workers = 16
    chunk_size = (total_len + num_parts - 1) // num_parts

    tasks = []
    for i in range(num_parts):
        start = i * chunk_size
        end = min(total_len - 1, (i + 1) * chunk_size - 1)
        part_file = parts_dir / f"part_{i:02d}.bin"
        tasks.append((i, start, end, part_file))

    print(f"Downloading {num_parts} parts using {workers} concurrent workers...")
    t0 = time.time()
    completed = 0

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futs = {
            executor.submit(download_chunk, URL, start, end, pfile): i
            for i, start, end, pfile in tasks
        }
        for f in as_completed(futs):
            idx = futs[f]
            try:
                size = f.result()
                completed += 1
                elapsed = time.time() - t0
                pct = (completed / num_parts) * 100
                print(f"[{completed:02d}/{num_parts}] ({pct:5.1f}%) Part {idx:02d} complete ({size:,} bytes) in {elapsed:.1f}s")
            except Exception as e:
                print(f"Error on part {idx}: {e}")
                raise e

    print(f"All {num_parts} parts downloaded! Merging into {out_gz}...")
    with open(out_gz, "wb") as outfile:
        for i in range(num_parts):
            pfile = parts_dir / f"part_{i:02d}.bin"
            with open(pfile, "rb") as infile:
                shutil.copyfileobj(infile, outfile)

    print(f"Assembled {out_gz} ({out_gz.stat().st_size:,} bytes). Decompressing to {out_fa}...")
    with gzip.open(out_gz, "rb") as f_in, open(out_fa, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    print(f"HRGM2 FASTA ready: {out_fa} ({out_fa.stat().st_size:,} bytes)")
    # Cleanup parts
    shutil.rmtree(parts_dir, ignore_errors=True)
    print("Cleaned up temporary part files. Done!")


if __name__ == "__main__":
    main()
