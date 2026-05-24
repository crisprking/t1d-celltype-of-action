"""02_ld_clumping: LD-clump T1D associations with PLINK 1.9 → 145 independent loci.

Extracted from notebook cell 16. Part of the t1d-celltype-of-action pipeline.
Run modules in numeric order — each consumes outputs from earlier steps.
"""

# Switched from Dropbox-hosted Phase 3 (returned HTML error from Kaggle IPs)
# to the MAGMA pre-built EUR panel at ctg.cncr.nl.
# That panel:
#   - already EUR-only (no PLINK2 conversion step)
#   - already PLINK 1.9 bed/bim/fam
#   - already MAF>0.01 filtered
#   - ~600 MB zipped
#   - hosted by Center for Neurogenomics, Vrije Universiteit Amsterdam
#   - genome build GRCh37; we clump by rsID so build doesn't matter
#
# References:
#   https://ctg.cncr.nl/software/magma          (panel home)
#   biostars.org/p/329901                       (community docs)

import os, subprocess, tempfile, zipfile
from pathlib import Path
import pandas as pd
import httpx
from tqdm.auto import tqdm

PROJECT = Path("/kaggle/working/t1d_mech")
RAW = PROJECT / "data" / "raw"
INTERIM = PROJECT / "data" / "interim"
PROCESSED = PROJECT / "data" / "processed"
TOOLS = PROJECT / "tools"
REF = PROJECT / "reference"
for d in (INTERIM, PROCESSED, TOOLS, REF):
    d.mkdir(parents=True, exist_ok=True)

# --- PLINK 1.9 install (skip if from previous run) -------------------------
PLINK = TOOLS / "plink"
if not PLINK.exists():
    log.info("Installing PLINK 1.9...")
    plink_zip = TOOLS / "plink_linux_x86_64.zip"
    url = "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
    with httpx.stream("GET", url, follow_redirects=True, timeout=120) as r:
        r.raise_for_status()
        with plink_zip.open("wb") as f:
            for chunk in r.iter_bytes(1 << 16):
                f.write(chunk)
    with zipfile.ZipFile(plink_zip) as zf:
        zf.extractall(TOOLS)
    PLINK.chmod(0o755)
ver = subprocess.run([str(PLINK), "--version"], capture_output=True, text=True)
print(ver.stdout.strip() or ver.stderr.strip())

# --- Probe both candidate URLs for the EUR panel before committing ---------
MAGMA_CANDIDATES = [
    "https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download",  # SURFsara mirror (newer)
    "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",         # original CNCR
]
EUR_ZIP = REF / "g1000_eur.zip"
EUR_PREFIX = REF / "g1000_eur"
EUR_BED = EUR_PREFIX.with_suffix(".bed")

def probe(url, timeout=15.0):
    try:
        with httpx.Client(timeout=timeout, follow_redirects=True) as c:
            r = c.head(url)
            if r.status_code == 405:
                r = c.get(url, headers={"Range": "bytes=0-0"})
            return r.status_code, int(r.headers.get("content-length") or 0) or None
    except httpx.HTTPError as e:
        return -1, None

if not EUR_BED.exists():
    log.info("Probing EUR panel mirrors...")
    chosen = None
    for url in MAGMA_CANDIDATES:
        status, size = probe(url)
        print(f"  [{status:>4}]  {size}  {url}")
        if status in (200, 206) and chosen is None:
            chosen = url
    if chosen is None:
        raise RuntimeError(
            "Neither MAGMA EUR panel mirror responded. Try a browser to confirm "
            "https://ctg.cncr.nl/software/magma is up; if it is, paste the current "
            "download link manually."
        )
    log.info("Using: %s", chosen)

    # Stream download
    fd, tmp_name = tempfile.mkstemp(prefix="g1000_eur.", suffix=".part", dir=REF)
    os.close(fd); tmp = Path(tmp_name)
    try:
        with httpx.stream("GET", chosen, follow_redirects=True, timeout=900) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length") or 0)
            written = 0
            with tmp.open("wb") as f:
                bar = tqdm(total=total or None, unit="B", unit_scale=True, desc="g1000_eur.zip")
                for chunk in r.iter_bytes(1 << 16):
                    f.write(chunk); bar.update(len(chunk)); written += len(chunk)
                bar.close()
        # Sanity: must be ≥100 MB, must be a real zip
        if tmp.stat().st_size < 100_000_000:
            head = tmp.read_bytes()[:200]
            raise RuntimeError(
                f"Too-small download ({tmp.stat().st_size:,} bytes). Probably HTML error. "
                f"Head: {head!r}"
            )
        with open(tmp, "rb") as f:
            magic = f.read(4)
        if magic[:2] != b"PK":
            raise RuntimeError(f"Not a zip file. Magic bytes: {magic!r}")
        tmp.replace(EUR_ZIP)
        log_artifact(EUR_ZIP, source=chosen, notes="MAGMA pre-built 1000G EUR PLINK panel")
    except Exception:
        if tmp.exists(): tmp.unlink(missing_ok=True)
        raise

    # Extract
    log.info("Extracting %s ...", EUR_ZIP.name)
    with zipfile.ZipFile(EUR_ZIP) as zf:
        names = zf.namelist()
        print(f"  Members: {names}")
        zf.extractall(REF)

    # Find the bed/bim/fam triplet and rename to a canonical prefix
    beds = list(REF.rglob("*.bed"))
    if not beds:
        raise RuntimeError("No .bed file in extracted archive.")
    actual = beds[0].with_suffix("")  # strip .bed
    log.info("Extracted prefix: %s", actual)

    # If extracted prefix isn't already "g1000_eur", symlink/rename .bed/.bim/.fam
    if actual != EUR_PREFIX:
        for ext in (".bed", ".bim", ".fam"):
            src = actual.with_suffix(ext)
            dst = EUR_PREFIX.with_suffix(ext)
            if src.exists() and not dst.exists():
                src.rename(dst)
        log.info("Renamed triplet to %s.*", EUR_PREFIX)

# Verify
for ext in (".bed", ".bim", ".fam"):
    p = EUR_PREFIX.with_suffix(ext)
    print(f"  {p.name}: {p.stat().st_size/1e6:.1f} MB" if p.exists() else f"  MISSING: {p.name}")

# How many SNPs and samples?
n_snps = sum(1 for _ in EUR_PREFIX.with_suffix(".bim").open())
n_samples = sum(1 for _ in EUR_PREFIX.with_suffix(".fam").open())
print(f"\nPanel: {n_samples} EUR samples, {n_snps:,} SNPs")

# --- Prepare clumping input -----------------------------------------------
t1d_gw = pd.read_csv(RAW / "gwas_catalog_t1d_gwsig.tsv", sep="\t", low_memory=False, dtype=str)
t1d_gw["P_NUM"] = pd.to_numeric(t1d_gw["P-VALUE"], errors="coerce")

clump_input = INTERIM / "t1d_clump_input.tsv"
best = (
    t1d_gw.dropna(subset=["SNPS", "P_NUM"])
          .sort_values("P_NUM")
          .drop_duplicates("SNPS", keep="first")
          [["SNPS", "P_NUM"]]
          .rename(columns={"SNPS": "SNP", "P_NUM": "P"})
)
best = best[~best["SNP"].str.contains(",", na=False)]
best = best[best["SNP"].str.startswith("rs", na=False)]
log.info("SNPs going into clumping: %d", len(best))
best.to_csv(clump_input, sep="\t", index=False)

# --- Run PLINK --clump ----------------------------------------------------
clump_out = INTERIM / "t1d_clumped"
cmd = [
    str(PLINK),
    "--bfile", str(EUR_PREFIX),
    "--clump", str(clump_input),
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", "5e-8",
    "--clump-p2", "1e-5",
    "--clump-r2", "0.1",
    "--clump-kb", "1000",
    "--out", str(clump_out),
]
print("Running:", " ".join(cmd))
result = subprocess.run(cmd, capture_output=True, text=True)
print(result.stdout[-2000:])
if result.returncode != 0:
    print("STDERR:", result.stderr[-2000:])
    raise RuntimeError("PLINK clump failed")

# --- Parse and save ------------------------------------------------------
clumped_file = Path(str(clump_out) + ".clumped")
if not clumped_file.exists():
    raise RuntimeError(f"Expected clump output not found: {clumped_file}")

clumped = pd.read_csv(clumped_file, sep=r"\s+", engine="python")
log.info("Independent loci after clumping: %d", len(clumped))

loci_path = PROCESSED / "t1d_independent_loci.tsv"
clumped.to_csv(loci_path, sep="\t", index=False)
log_artifact(loci_path,
             source="PLINK --clump (r2<0.1, 1Mb) on T1D GW-sig SNPs vs MAGMA 1000G EUR",
             notes=f"{len(clumped)} independent T1D loci")

clumped["n_absorbed"] = clumped["SP2"].fillna("").apply(
    lambda s: 0 if s in ("", "NONE") else len(s.split(","))
)

print(f"\n=== T1D locus inventory (post-clumping) ===")
print(f"  Independent loci:         {len(clumped):,}")
print(f"  Median p of lead SNP:     {clumped['P'].median():.2e}")
print(f"  Chromosomes with signal:  {clumped['CHR'].nunique()}")

print(f"\nTop 15 most-loaded loci (expect HLA/INS/PTPN22 to dominate):")
top = clumped.nlargest(15, "n_absorbed")[["CHR", "SNP", "BP", "P", "n_absorbed"]]
print(top.to_string(index=False))

print(f"\nFirst 20 independent loci (ordered as PLINK reports them):")
print(clumped[["CHR", "SNP", "BP", "P", "n_absorbed"]].head(20).to_string(index=False))

print(f"\nSaved: {loci_path}")
