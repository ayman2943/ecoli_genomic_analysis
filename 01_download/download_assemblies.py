#!/usr/bin/env python3
"""
Download E. coli genome assemblies listed in an EnteroBase metadata Excel file.

Parallel downloader with rate limiting. Uses the NCBI datasets tool (FTP first,
datasets fallback). Requires the "datasets" CLI from NCBI (conda install -c conda-forge ncbi-datasets-cli).

Usage:
    python 01_download/download_assemblies.py metadata/ST69_filtered.xlsx ST69

Credentials:
    The NCBI API key is read from the NCBI_API_KEY environment variable (or a
    .env file in the repo root). It is NEVER hardcoded. Without a key the
    script runs at the slower 3 req/sec rate.

Output:
    ST69/Escherichia_coli_<strain>.fna.gz   (one assembly FASTA per isolate)
    ST69/download_log_*.txt
    ST69/failed_accessions.txt
"""

import os
import re
import shutil
import subprocess
import sys
import pandas as pd
from multiprocessing import Pool
from datetime import datetime
import time
from functools import lru_cache
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# ================= SETTINGS =================
def parse_args(argv):
    if len(argv) < 2:
        sys.exit("Usage: download_assemblies.py <input_metadata.xlsx> <outdir> [species]")
    return argv[0], argv[1], argv[2] if len(argv) > 2 else "Escherichia_coli"

INPUT_FILE, OUTDIR, SPECIES = parse_args(sys.argv[1:])
MAX_WORKERS = int(os.environ.get("DOWNLOAD_WORKERS", "16"))
MAX_RETRIES = 3

# ================= NCBI API SETTINGS =================
# Key is provided via env var only - never commit a key to the repo.
def load_api_key():
    key = os.environ.get("NCBI_API_KEY", "").strip()
    if not key:
        env_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".env")
        if os.path.isfile(env_file):
            with open(env_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("NCBI_API_KEY="):
                        key = line.split("=", 1)[1].strip().strip("'\"")
                        break
    return key

NCBI_API_KEY = load_api_key()

if NCBI_API_KEY:
    os.environ["NCBI_API_KEY"] = NCBI_API_KEY
    REQUEST_DELAY = 0.11   # ~9 req/sec - safe margin under 10/sec limit
else:
    REQUEST_DELAY = 0.35   # 3 req/sec without key

# ================= OUTPUT FILES =================
os.makedirs(OUTDIR, exist_ok=True)
LOG_FILE = f"{OUTDIR}/download_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
FAILED_FILE = f"{OUTDIR}/failed_accessions.txt"

# Rate limiting with token bucket
class RateLimiter:
    def __init__(self, rate=10):
        self.rate = rate
        self.tokens = rate
        self.last_update = time.time()

    def acquire(self):
        now = time.time()
        elapsed = now - self.last_update
        self.tokens = min(self.rate, self.tokens + elapsed * self.rate)
        self.last_update = now

        if self.tokens >= 1:
            self.tokens -= 1
            return
        else:
            sleep_time = (1 - self.tokens) / self.rate
            time.sleep(sleep_time)
            self.tokens = 0

rate_limiter = RateLimiter(rate=9 if NCBI_API_KEY else 3)

# ────────────────────────────────────────────────
def log_message(message, also_print=False):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as f:
        f.write(f"{timestamp} - {message}\n")
    if also_print:
        print(f"{timestamp} - {message}")

# ────────────────────────────────────────────────
def find_all_sra_accessions(row):
    ncbi_sra, ddbj_sra = [], []
    for value in row:
        if pd.notna(value):
            ncbi_sra += re.findall(r"\b([SE]RR\d{6,})\b", str(value))
            ddbj_sra += re.findall(r"\b(DRR\d{6,})\b", str(value))
    return list(dict.fromkeys(ncbi_sra)), list(dict.fromkeys(ddbj_sra))

def find_assembly_accession(row):
    for value in row:
        if pd.notna(value):
            m = re.search(r"\b(GC[FA]_\d+(?:\.\d+)?)\b", str(value))
            if m:
                return m.group(1)
    return None

def find_biosample_accession(row):
    for value in row:
        if pd.notna(value):
            m = re.search(r"\b(SAM[NDE][A-Z]?\d+)\b", str(value))
            if m:
                return m.group(1)
    return None

def get_strain_name(row, name_col_idx):
    if name_col_idx is not None:
        v = row.iloc[name_col_idx]
        if pd.notna(v):
            return str(v).strip()
    return ""

# ────────────────────────────────────────────────
def run_ncbi_command_batch(cmd, retries=MAX_RETRIES, timeout=120):
    """Run NCBI command with rate limiting"""
    for attempt in range(retries):
        try:
            rate_limiter.acquire()
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            if result.returncode == 0:
                out = result.stdout.strip()
                if out:
                    return out
            else:
                if attempt == retries - 1:
                    log_message(f"Command failed: {result.stderr.strip()[:200]}")
        except subprocess.TimeoutExpired:
            log_message(f"Timeout after {timeout}s")
        except Exception as e:
            log_message(f"Exception: {str(e)}")

        if attempt < retries - 1:
            time.sleep(2 ** attempt)  # exponential backoff
    return None

# ────────────────────────────────────────────────
@lru_cache(maxsize=10000)
def get_biosample_from_sra_cached(sra):
    cmd = (
        f'esearch -db sra -query "{sra}" | '
        f'elink -target biosample | '
        f'efetch -format docsum | '
        f'xtract -pattern DocumentSummary -element Accession'
    )
    out = run_ncbi_command_batch(cmd, timeout=30)
    if out and out.startswith("SAM"):
        return out.splitlines()[0].strip()
    return None

@lru_cache(maxsize=10000)
def get_assembly_from_biosample_cached(biosample):
    if not biosample:
        return None

    cmd = (
        f'esearch -db biosample -query "{biosample}" | '
        f'elink -target assembly | '
        f'efetch -format docsum | '
        f'xtract -pattern DocumentSummary -if AssemblyAccession -element AssemblyAccession'
    )
    out = run_ncbi_command_batch(cmd, timeout=30)
    if out:
        for line in out.splitlines():
            a = line.strip()
            if a.startswith("GCF_") or a.startswith("GCA_"):
                return a
    return None

@lru_cache(maxsize=5000)
def verify_assembly_exists_cached(assembly):
    if not assembly:
        return False
    cmd = f'esearch -db assembly -query "{assembly}"'
    out = run_ncbi_command_batch(cmd, retries=2, timeout=20)
    return bool(out and "Count: 0" not in out)

# ────────────────────────────────────────────────
def download_assembly_ftp(assembly, out_prefix):
    """Download via FTP - much faster than datasets command"""
    final_fna = f"{out_prefix}.fna.gz"

    if os.path.exists(final_fna) and os.path.getsize(final_fna) > 10000:
        return True, "already exists"

    cmd = (
        f'esearch -db assembly -query "{assembly}" | '
        f'efetch -format docsum | '
        f'xtract -pattern DocumentSummary -element FtpPath_GenBank'
    )

    ftp_path = run_ncbi_command_batch(cmd, timeout=30)
    if not ftp_path or "ftp://" not in ftp_path:
        return download_assembly_datasets(assembly, out_prefix)

    https_path = ftp_path.replace("ftp://", "https://")
    base_name = os.path.basename(https_path)
    genome_url = f"{https_path}/{base_name}_genomic.fna.gz"

    try:
        response = requests.get(genome_url, stream=True, timeout=300)
        if response.status_code == 200:
            with open(final_fna, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            log_message(f"Downloaded via FTP: {assembly}")
            return True, "downloaded (FTP)"
        else:
            return download_assembly_datasets(assembly, out_prefix)
    except Exception as e:
        log_message(f"FTP download failed: {e}")
        return download_assembly_datasets(assembly, out_prefix)

def download_assembly_datasets(assembly, out_prefix):
    """Fallback: use datasets command"""
    final_fna = f"{out_prefix}.fna.gz"
    zipf = f"{out_prefix}.zip"
    tmpd = f"{out_prefix}_tmp"

    for attempt in range(MAX_RETRIES):
        try:
            cmd = [
                "datasets", "download", "genome", "accession", assembly,
                "--include", "genome", "--filename", zipf, "--no-progressbar"
            ]
            if NCBI_API_KEY:
                cmd += ["--api-key", NCBI_API_KEY]

            subprocess.run(cmd, check=True, capture_output=True, timeout=600)

            os.makedirs(tmpd, exist_ok=True)
            subprocess.run(["unzip", "-q", "-o", zipf, "-d", tmpd], check=True)

            for root, _, files in os.walk(tmpd):
                for f in files:
                    if f.endswith((".fna", ".fna.gz", ".fa", ".fasta")):
                        src = os.path.join(root, f)
                        dst = final_fna if f.endswith(".gz") else f"{out_prefix}.fna"
                        shutil.move(src, dst)
                        return True, "downloaded (datasets)"

        except Exception as e:
            if attempt == MAX_RETRIES - 1:
                log_message(f"Download failed: {str(e)}")
        finally:
            for f in [zipf]:
                if os.path.exists(f):
                    os.unlink(f)
            shutil.rmtree(tmpd, ignore_errors=True)

    return False, "download failed"

# ────────────────────────────────────────────────
def process_sample(args):
    i, total, acc_type, acc, prefix, strain = args
    label = strain or (acc[0] if isinstance(acc, list) else acc)
    full_prefix = os.path.join(OUTDIR, f"{SPECIES}_{prefix.replace(' ', '_')}")

    print(f"[{i}/{total}] Processing {label} ({acc_type})")

    if acc_type == "assembly":
        ok, detail = download_assembly_ftp(acc, full_prefix)
        status = "success" if ok else "failed"
        return (status, f"[{i}/{total}] {label}", acc, detail)

    elif acc_type == "biosample":
        asm = get_assembly_from_biosample_cached(acc)
        if not asm:
            return ("no_assembly", f"[{i}/{total}] {label}", acc, "no assembly")
        ok, detail = download_assembly_ftp(asm, full_prefix)
        status = "success" if ok else "failed"
        return (status, f"[{i}/{total}] {label}", acc, f"{detail} ({asm})")

    elif acc_type == "sra":
        for sra in acc:
            if sra.startswith("DRR"):
                continue
            bio = get_biosample_from_sra_cached(sra)
            if not bio:
                continue
            asm = get_assembly_from_biosample_cached(bio)
            if not asm:
                continue
            ok, detail = download_assembly_ftp(asm, full_prefix)
            if ok:
                return ("success", f"[{i}/{total}] {label}", sra, f"downloaded ({asm})")
        return ("no_assembly", f"[{i}/{total}] {label}", acc[0] if acc else "?", "no assembly")

    return ("failed", f"[{i}/{total}] {label}", acc, "unknown")

# ────────────────────────────────────────────────
def main():
    print("=" * 70)
    print(f"E. coli Genome Downloader (OPTIMIZED) — {SPECIES}")
    print("=" * 70)
    print(f"Input: {INPUT_FILE}")
    print(f"Outdir: {OUTDIR}")
    print(f"API key: {'YES' if NCBI_API_KEY else 'NO'}")
    print(f"Workers: {MAX_WORKERS} | Rate: {9 if NCBI_API_KEY else 3} req/sec")

    try:
        df = pd.read_excel(INPUT_FILE)
    except Exception as e:
        print(f"Error reading Excel: {e}")
        sys.exit(1)

    name_col_idx = None
    for i, c in enumerate(df.columns):
        if str(c).strip().lower() in ["name", "strain", "isolate"]:
            name_col_idx = i
            break

    targets = []
    for idx, row in df.iterrows():
        strain = get_strain_name(row, name_col_idx)
        asm = find_assembly_accession(row)
        if asm:
            prefix = f"{strain or asm}".replace(" ", "_").replace("/", "_")
            targets.append(("assembly", asm, prefix, strain))
            continue

        bio = find_biosample_accession(row)
        sra_list, _ = find_all_sra_accessions(row)

        if bio:
            prefix = f"{strain or bio}".replace(" ", "_").replace("/", "_")
            targets.append(("biosample", bio, prefix, strain))
        elif sra_list:
            prefix = f"{strain or sra_list[0]}".replace(" ", "_").replace("/", "_")
            targets.append(("sra", sra_list, prefix, strain))

    print(f"Found {len(targets)} targets")
    if not targets:
        print("No valid targets found.")
        sys.exit(0)

    if os.environ.get("DOWNLOAD_ASSUME_YES", "0") != "1":
        if input("Continue? (y/n): ").strip().lower() != "y":
            print("Aborted.")
            sys.exit(0)

    start_time = time.time()

    with Pool(MAX_WORKERS) as pool:
        results = pool.imap_unordered(
            process_sample,
            [(i+1, len(targets), *t) for i, t in enumerate(targets)]
        )

        success_count = 0
        for status, msg, acc, detail in results:
            symbol = "OK" if status == "success" else "XX"
            print(f"{symbol} {msg} — {detail}")
            if status == "success":
                success_count += 1

    elapsed = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"Finished in {elapsed/60:.1f} minutes")
    print(f"Success: {success_count}/{len(targets)}")
    print(f"Speed: {len(targets)/elapsed*60:.1f} genomes/hour")
    print(f"Log: {LOG_FILE}")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
