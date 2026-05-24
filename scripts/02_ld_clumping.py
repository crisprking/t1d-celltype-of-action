"""Step 02 — LD-clump T1D GW-sig SNPs against the 1000G EUR panel.

PLINK 1.9's ``--clump`` collapses tagging SNPs into independent loci. We
use the MAGMA pre-built EUR PLINK panel because it ships MAF-filtered
and PLINK-ready, avoiding a PLINK2 → PLINK1 conversion step on Kaggle.

Output: ``data/processed/t1d_independent_loci.tsv`` (≈145 loci).
"""

from __future__ import annotations

import os
import subprocess
import tempfile
import zipfile
from pathlib import Path

import httpx
import pandas as pd
from tqdm.auto import tqdm

from t1d_coa.config import INTERIM, PROCESSED, RAW, REFERENCE, TOOLS
from t1d_coa.provenance import log_artifact


PLINK_URL = (
    "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip"
)
EUR_MIRRORS = (
    "https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download",
    "https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip",
)
CLUMP_PARAMS = {
    "--clump-p1": "5e-8",
    "--clump-p2": "1e-5",
    "--clump-r2": "0.1",
    "--clump-kb": "1000",
}


def install_plink() -> Path:
    plink = TOOLS / "plink"
    if plink.exists():
        return plink
    TOOLS.mkdir(parents=True, exist_ok=True)
    zpath = TOOLS / "plink.zip"
    with httpx.stream("GET", PLINK_URL, follow_redirects=True, timeout=120) as r:
        r.raise_for_status()
        with zpath.open("wb") as f:
            for chunk in r.iter_bytes(1 << 16):
                f.write(chunk)
    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(TOOLS)
    plink.chmod(0o755)
    return plink


def download_eur_panel(prefix: Path) -> None:
    if prefix.with_suffix(".bed").exists():
        return
    REFERENCE.mkdir(parents=True, exist_ok=True)

    chosen = None
    for url in EUR_MIRRORS:
        try:
            with httpx.Client(timeout=15, follow_redirects=True) as c:
                r = c.head(url)
                if r.status_code in (200, 206):
                    chosen = url
                    break
        except httpx.HTTPError:
            continue
    if chosen is None:
        raise RuntimeError("No MAGMA EUR mirror responded; check upstream.")

    zpath = REFERENCE / "g1000_eur.zip"
    fd, tmp_name = tempfile.mkstemp(prefix="g1000_eur.", suffix=".part",
                                    dir=REFERENCE)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        with httpx.stream("GET", chosen, follow_redirects=True, timeout=900) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length") or 0)
            bar = tqdm(total=total or None, unit="B", unit_scale=True,
                       desc="g1000_eur.zip")
            with tmp.open("wb") as f:
                for chunk in r.iter_bytes(1 << 16):
                    f.write(chunk)
                    bar.update(len(chunk))
            bar.close()
        if tmp.stat().st_size < 100_000_000:
            raise RuntimeError(f"Suspiciously small download: {tmp.stat().st_size:,} bytes")
        if tmp.open("rb").read(2) != b"PK":
            raise RuntimeError("Downloaded file is not a ZIP archive")
        tmp.replace(zpath)
        log_artifact(zpath, source=chosen,
                     notes="MAGMA pre-built 1000G EUR PLINK panel")
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)

    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(REFERENCE)

    # Locate the extracted .bed and standardize the prefix.
    extracted = next(REFERENCE.rglob("*.bed")).with_suffix("")
    if extracted != prefix:
        for ext in (".bed", ".bim", ".fam"):
            src = extracted.with_suffix(ext)
            dst = prefix.with_suffix(ext)
            if src.exists() and not dst.exists():
                src.rename(dst)


def prepare_clump_input() -> Path:
    gw = pd.read_csv(RAW / "gwas_catalog_t1d_gwsig.tsv", sep="\t",
                     low_memory=False, dtype=str)
    gw["P_NUM"] = pd.to_numeric(gw["P-VALUE"], errors="coerce")
    best = (
        gw.dropna(subset=["SNPS", "P_NUM"])
        .sort_values("P_NUM")
        .drop_duplicates("SNPS", keep="first")[["SNPS", "P_NUM"]]
        .rename(columns={"SNPS": "SNP", "P_NUM": "P"})
    )
    # Drop multi-SNP rows and any non-rsID identifiers.
    best = best[~best["SNP"].str.contains(",", na=False)]
    best = best[best["SNP"].str.startswith("rs", na=False)]

    INTERIM.mkdir(parents=True, exist_ok=True)
    path = INTERIM / "t1d_clump_input.tsv"
    best.to_csv(path, sep="\t", index=False)
    return path


def main() -> None:
    plink = install_plink()
    eur_prefix = REFERENCE / "g1000_eur"
    download_eur_panel(eur_prefix)
    clump_in = prepare_clump_input()

    clump_out = INTERIM / "t1d_clumped"
    cmd = [
        str(plink),
        "--bfile", str(eur_prefix),
        "--clump", str(clump_in),
        "--clump-snp-field", "SNP",
        "--clump-field", "P",
        *[v for pair in CLUMP_PARAMS.items() for v in pair],
        "--out", str(clump_out),
    ]
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout[-1500:])
    if result.returncode != 0:
        print("STDERR:", result.stderr[-1500:])
        raise RuntimeError("PLINK clump failed")

    clumped = pd.read_csv(f"{clump_out}.clumped", sep=r"\s+", engine="python")
    loci_path = PROCESSED / "t1d_independent_loci.tsv"
    clumped.to_csv(loci_path, sep="\t", index=False)
    log_artifact(
        loci_path,
        source="PLINK --clump (r2<0.1, 1Mb) on T1D GW-sig vs MAGMA 1000G EUR",
        notes=f"{len(clumped)} independent T1D loci",
    )
    print(f"\nIndependent loci: {len(clumped)}  →  {loci_path}")


if __name__ == "__main__":
    main()
