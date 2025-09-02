#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Descubrimiento de FASTQ (paired/single-end), generación de manifiesto y FastQC crudo.
Respeta la prioridad de rutas: MAG_PROJECT_DIR > MAGENTA_DIR > cwd.
"""

from __future__ import annotations

# --- Importaciones ---
from pathlib import Path
import os, re, sys, shutil, subprocess, csv, json, gzip
from datetime import datetime
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import glob
from typing import Dict, List, Optional

# --- Configuración y rutas (MAGENTA_DIR congruente) ---

# Raíz del proyecto: MAG_PROJECT_DIR tiene prioridad, si no MAGENTA_DIR, si no cwd.
PROJECT_DIR = Path(
    os.environ.get("MAG_PROJECT_DIR", os.environ.get("MAGENTA_DIR", Path.cwd()))
).resolve()

# Estructura base tipo MAGs dentro del proyecto
MAGS_DIR    = PROJECT_DIR / "mags"
DATA_DIR    = MAGS_DIR / "data"            # (no se usa como default de RAW)
RESULTS_DIR = MAGS_DIR / "results"
SCRIPTS_DIR = MAGS_DIR / "scripts"

# Subcarpetas de resultados
QC_DIR      = RESULTS_DIR / "01.qc"
TRIM_DIR    = RESULTS_DIR / "02.trimmed"
ASM_DIR     = RESULTS_DIR / "03.assembly"
ASM_LOG_DIR = ASM_DIR / "logs"

# Metadatos / reportes
PIPE_META          = RESULTS_DIR / "pipeline_meta"
PIPE_META.mkdir(parents=True, exist_ok=True)
PIPE_MANIFEST_CSV  = PIPE_META / "pipeline_manifest.csv"
TRIM_REPORT_CSV    = PIPE_META / "trim_report.csv"
ASM_RESUMEN_CSV    = ASM_DIR   / "resumen_ensamblaje.csv"
VALID_SAMPLES_CSV  = PIPE_META / "muestras_validas.csv"

# Donde 02_descargar_y_convertir_mangrove.py deja los FASTQ
DEFAULT_RAW_SRC = PROJECT_DIR / "rawdata" / "convertidos"
RAW_SRC = Path(os.environ.get("MAG_RAW_DIR", DEFAULT_RAW_SRC))

# Crear estructura necesaria
for d in [MAGS_DIR, DATA_DIR, RESULTS_DIR, SCRIPTS_DIR, QC_DIR, TRIM_DIR, ASM_DIR, ASM_LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Parámetros
N_THREADS_FASTQC   = int(os.environ.get("MAG_QC_THREADS",        "12"))
N_THREADS_FASTP    = int(os.environ.get("MAG_FASTP_THREADS",     "12"))
MAX_WORKERS_FASTP  = int(os.environ.get("MAG_FASTP_MAX_WORKERS", "4"))
MAX_READS_VALIDATE = int(os.environ.get("MAG_VALIDATE_READS",    "10000"))
MEGAHIT_THREADS    = int(os.environ.get("MAG_MEGAHIT_THREADS",   "40"))
SAMPLE_REGEX       = os.environ.get(
    "MAG_SAMPLE_REGEX",
    r"^([A-Za-z0-9_\-]+)_([12])\.fastq(?:\.gz)?$"
)

def have_bin(b: str) -> bool:
    return shutil.which(b) is not None

def run(cmd, **kwargs) -> int:
    printable = " ".join(map(str, cmd)) if isinstance(cmd, (list,tuple)) else str(cmd)
    print("[CMD]", printable)
    return subprocess.call(cmd, shell=isinstance(cmd, str), **kwargs)

print("[PORTABLE] PROJECT_DIR =", PROJECT_DIR)
print("[PORTABLE] RAW_SRC     =", RAW_SRC)
print("[PORTABLE] RESULTS_DIR  =", RESULTS_DIR)

# --- Descubrimiento de muestras + Manifiesto + FastQC (crudo) ---

# 1) Detectar pares (y single-end) a partir de SAMPLE_REGEX
def list_pairs(raw_dir: Path, sample_regex: str) -> Dict[str, Dict[str, Optional[Path]]]:
    rx = re.compile(sample_regex)
    buckets: Dict[str, Dict[str, Optional[Path]]] = {}
    # Paired-end con sufijos _1 / _2
    for p in sorted(raw_dir.rglob("*.fastq*")):
        m = rx.match(p.name)
        if not m:
            continue
        sample, read = m.group(1), m.group(2)
        d = buckets.setdefault(sample, {"R1": None, "R2": None})
        if read == "1" and d["R1"] is None:
            d["R1"] = p
        elif read == "2" and d["R2"] is None:
            d["R2"] = p
    # Single-end con nombres <sample>.fastq(.gz)
    for p in sorted(raw_dir.rglob("*.fastq*")):
        if rx.match(p.name):
            continue
        base = p.name.replace(".gz", "")
        if base.endswith(".fastq"):
            sample = base[:-6]  # quita ".fastq"
            buckets.setdefault(sample, {"R1": None, "R2": None})
            if buckets[sample]["R1"] is None:
                buckets[sample]["R1"] = p
    return buckets

def write_manifest(pairs: Dict[str, Dict[str, Optional[Path]]]) -> pd.DataFrame:
    rows = []
    for s, rr in sorted(pairs.items()):
        rows.append({"sample": s, "R1": str(rr.get("R1") or ""), "R2": str(rr.get("R2") or "")})
    df = pd.DataFrame(rows)
    df.to_csv(PIPE_MANIFEST_CSV, index=False)
    print(f"[META] Manifiesto muestras -> {PIPE_MANIFEST_CSV}")
    return df

# 2) FastQC crudo (paralelo por lotes para no saturar la CLI)
def run_fastqc(inputs: List[Path], outdir: Path, threads: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    if not inputs:
        print("[FASTQC] Sin archivos de entrada.")
        return
    batch: List[str] = []
    for p in inputs:
        batch.append(str(p))
        if len(batch) >= 64:
            run(["fastqc", "-t", str(threads), "-o", str(outdir)] + batch)
            batch = []
    if batch:
        run(["fastqc", "-t", str(threads), "-o", str(outdir)] + batch)

def main() -> int:
    # Validaciones iniciales
    if not RAW_SRC.exists():
        print(f"[ERROR] No existe RAW_SRC: {RAW_SRC}", file=sys.stderr)
        return 2

    pairs = list_pairs(RAW_SRC, SAMPLE_REGEX)
    if not pairs:
        print(f"[ERROR] No se detectaron FASTQ en {RAW_SRC}", file=sys.stderr)
        return 3

    manifest_df = write_manifest(pairs)

    raw_fastqs: List[Path] = []
    for rr in pairs.values():
        if rr.get("R1"): raw_fastqs.append(Path(rr["R1"]))  # type: ignore[arg-type]
        if rr.get("R2"): raw_fastqs.append(Path(rr["R2"]))  # type: ignore[arg-type]

    if have_bin("fastqc"):
        outdir = QC_DIR / "raw_fastqc"
        run_fastqc(raw_fastqs, outdir, N_THREADS_FASTQC)
        if have_bin("multiqc"):
            run(["multiqc", str(outdir), "-o", str(outdir)])
    else:
        print("[AVISO] 'fastqc' no está en PATH. Se omite QC inicial.")

    print("[OK] Proceso completado.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

