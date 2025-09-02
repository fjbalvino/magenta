#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
06_assembly_serial.py 
"""

from __future__ import annotations
import os, sys, csv, shutil, subprocess, argparse
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional

# ==========
# Rutas base
# ==========
PROJECT_DIR = Path(os.environ.get("MAG_PROJECT_DIR", os.environ.get("MAGENTA_DIR", Path.cwd()))).resolve()
MAGS_DIR    = PROJECT_DIR / "mags"
RESULTS_DIR = MAGS_DIR / "results"

TRIM_DIR    = RESULTS_DIR / "02.trimmed"
ASM_DIR     = RESULTS_DIR / "03.assembly"
ASM_LOG_DIR = ASM_DIR / "logs"
RESUMEN_CSV = ASM_DIR / "resumen_ensamblaje.csv"

ASM_LOG_DIR.mkdir(parents=True, exist_ok=True)

# ==========
# Utilidades
# ==========
def have_bin(b: str) -> bool:
    from shutil import which
    return which(b) is not None

def run(cmd, log_file: Path | None = None) -> int:
    printable = " ".join(map(str, cmd))
    print("[CMD]", printable, flush=True)
    if log_file:
        with open(log_file, "a") as log:
            log.write(f"[CMD] {printable}\n")
            return subprocess.call(cmd, stdout=log, stderr=log)
    else:
        return subprocess.call(cmd)

def size_sum(f1: Path, f2: Path) -> int:
    try:
        return f1.stat().st_size + f2.stat().st_size
    except Exception:
        return 0

# =========================
# Descubrimiento de pares
# =========================
R1_GLOBS_DEFAULT = [
    "*_R1*.fastq.gz", "*_1*.fastq.gz", "*R1*.fastq.gz", "*1*.fastq.gz",
    "*_R1*.fq.gz",    "*_1*.fq.gz",    "*R1*.fq.gz",    "*1*.fq.gz",
]
TOKEN_PAIRS = [("_R1", "_R2"), ("_1", "_2"), ("R1", "R2"), ("1", "2")]

def _try_match_r2_by_tokens(r1: Path) -> Optional[Path]:
    """Intenta encontrar R2 reemplazando tokens comunes en el nombre de R1, en el mismo directorio."""
    name = r1.name
    for t1, t2 in TOKEN_PAIRS:
        if t1 in name:
            cand = r1.with_name(name.replace(t1, t2, 1))
            if cand.exists():
                return cand
    return None

def _best_effort_r2_in_dir(r1: Path) -> Optional[Path]:
    """
    Si el reemplazo directo falla, intenta localizar un R2 plausible en el mismo directorio:
    - Algún archivo que contenga un token R2 y termine en .f*q.gz
    - Preferir el que comparte mayor prefijo con R1.
    """
    dirp = r1.parent
    candidates: List[Path] = []
    for patt in ["*_R2*.fastq.gz", "*_2*.fastq.gz", "*R2*.fastq.gz", "*2*.fastq.gz",
                 "*_R2*.fq.gz",   "*_2*.fq.gz",   "*R2*.fq.gz",   "*2*.fq.gz"]:
        candidates.extend(dirp.glob(patt))
    if not candidates:
        return None

    def common_prefix_len(a: str, b: str) -> int:
        n = min(len(a), len(b)); i = 0
        while i < n and a[i] == b[i]:
            i += 1
        return i

    best = max(candidates, key=lambda p: common_prefix_len(p.name, r1.name))
    return best

def _collect_sample_dirs(root: Path) -> List[Path]:
    """Devuelve subdirectorios de primer nivel si existen; si no, devuelve [root] para fallback plano."""
    subdirs = [d for d in sorted(root.iterdir()) if d.is_dir()]
    return subdirs if subdirs else [root]

def find_pairs(trim_root: Path, debug: bool=False) -> List[Tuple[str, Path, Path]]:
    """
    Busca pares R1/R2 en estructura:
      trim_root/<sample>/**/{R1,R2}*.f*q.gz
    - Explora recursivamente debajo de cada subcarpeta de muestra.
    - Acepta múltiples patrones de nombres.
    Devuelve: lista de (sample, R1, R2) con sample = nombre del subdirectorio de primer nivel
              o, en fallback plano, el prefijo derivado del archivo R1.
    """
    pairs: List[Tuple[str, Path, Path]] = []
    sample_dirs = _collect_sample_dirs(trim_root)

    for sample_dir in sample_dirs:
        if sample_dir == trim_root and sample_dir.is_dir() and not any(sample_dir.iterdir()):
            continue

        sample_name = sample_dir.name if sample_dir != trim_root else None

        r1_candidates: List[Path] = []
        for patt in R1_GLOBS_DEFAULT:
            r1_candidates.extend(sample_dir.rglob(patt))
        r1_candidates = sorted(set(r1_candidates))

        if debug:
            print(f"[DEBUG] Dir muestra: {sample_dir}")
            for r1 in r1_candidates:
                print(f"[DEBUG]   R1 cand: {r1}")

        used: set[Path] = set()
        for r1 in r1_candidates:
            if r1 in used:
                continue
            r2 = _try_match_r2_by_tokens(r1)
            if r2 is None or not r2.exists():
                r2 = _best_effort_r2_in_dir(r1)
            if r2 and r2.exists() and r2 not in used:
                used.add(r1); used.add(r2)
                if sample_dir == trim_root:
                    sname = r1.name
                    for suf in ["_R1.fastq.gz","_1.fastq.gz","R1.fastq.gz","1.fastq.gz",
                                "_R1.fq.gz","_1.fq.gz","R1.fq.gz","1.fq.gz"]:
                        if sname.endswith(suf):
                            sname = sname[:-len(suf)]
                            break
                else:
                    sname = sample_name
                pairs.append((sname, r1, r2))

        if debug and not pairs:
            print(f"[DEBUG]   No se encontraron pares en {sample_dir}")

    return pairs

# =========================
# Normalización de salidas
# =========================
def normalize_and_prefix_outputs_megahit(sample: str, out_dir: Path, log: Path):
    """
    MEGAHIT produce final.contigs.fa.
    Crear contigs.fasta y renombrar a <sample>_contigs.fasta para uniformidad.
    """
    src = out_dir / "final.contigs.fa"
    mid = out_dir / "contigs.fasta"
    if src.exists():
        try:
            shutil.copyfile(src, mid)
        except Exception as e:
            with open(log, "a") as lg:
                lg.write(f"[WARN] No se pudo copiar a contigs.fasta: {e}\n")
    for fname in ["contigs.fasta", "scaffolds.fasta"]:
        source = out_dir / fname
        if source.exists():
            renamed = out_dir / f"{sample}_{fname}"
            if not renamed.exists():
                source.rename(renamed)
            print(f"[RENAME] {renamed.name}")

def prefix_outputs_metaspades(sample: str, out_dir: Path):
    """
    MetaSPAdes produce contigs.fasta y scaffolds.fasta.
    Renombrar con prefijo de muestra.
    """
    for fname in ["contigs.fasta", "scaffolds.fasta"]:
        source = out_dir / fname
        if source.exists():
            renamed = out_dir / f"{sample}_{fname}"
            if not renamed.exists():
                source.rename(renamed)
            print(f"[RENAME] {renamed.name}")

# ==========
# Ensamblaje
# ==========
def assemble_one(sample: str, R1: Path, R2: Path, assembler: str, threads: int, mem_gb: int,
                 overwrite: bool=False) -> tuple[str, str, str]:
    sample_outdir = ASM_DIR / f"{sample}_assembled"
    log_file = ASM_LOG_DIR / f"{sample}.log"

    # Manejo del directorio de salida
    if sample_outdir.exists():
        if overwrite:
            try:
                shutil.rmtree(sample_outdir)
            except Exception as e:
                return (sample, "ERROR", f"No se pudo borrar salida previa: {e}")
        else:
            return (sample, "SKIPPED", "Ya existe (use --overwrite para rehacer)")

    try:
        if assembler == "megahit":
            if not have_bin("megahit"):
                return (sample, "ERROR", "megahit no encontrado en PATH")
            cmd = [
                "megahit",
                "-1", str(R1),
                "-2", str(R2),
                "-o", str(sample_outdir),   # dejar que megahit cree la carpeta
                "-t", str(threads),
                "--presets", "meta-sensitive",
                # Correcciones añadidas (flags extra de ejemplo):
                "--min-contig-len", "1000",
                "--k-list", "21,29,39,59,79,99,119",
            ]
            code = run(cmd, log_file)
            if code != 0:
                try:
                    print(f"[LOG TAIL] {log_file}")
                    os.system(f"tail -n 200 {log_file}")
                except Exception:
                    pass
                return (sample, "ERROR", f"megahit exit code {code}")
            normalize_and_prefix_outputs_megahit(sample, sample_outdir, log_file)

        elif assembler == "metaspades":
            exe = "metaspades.py" if have_bin("metaspades.py") else ("spades.py" if have_bin("spades.py") else None)
            if exe is None:
                return (sample, "ERROR", "metaspades.py no encontrado en PATH")
            cmd = [exe, "-1", str(R1), "-2", str(R2), "-o", str(sample_outdir), "-t", str(threads)]
            if mem_gb and int(mem_gb) > 0:
                cmd += ["-m", str(int(mem_gb))]
            # Correcciones añadidas (flag extra de ejemplo):
            cmd += ["--only-assembler"]  # si ya tienes corrección hecha
            # cmd += ["--meta"]  # (metaspades.py ya implica modo meta)
            code = run(cmd, log_file)
            if code != 0:
                try:
                    print(f"[LOG TAIL] {log_file}")
                    os.system(f"tail -n 200 {log_file}")
                except Exception:
                    pass
                return (sample, "ERROR", f"metaspades exit code {code}")
            prefix_outputs_metaspades(sample, sample_outdir)

        else:
            return (sample, "ERROR", f"Assembler no soportado: {assembler}")

        return (sample, "OK", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    except Exception as e:
        try:
            print(f"[LOG TAIL] {log_file}")
            os.system(f"tail -n 200 {log_file}")
        except Exception:
            pass
        return (sample, "ERROR", str(e))

def obtener_muestras_pendientes(pairs: List[Tuple[str, Path, Path]]) -> List[str]:
    """
    Filtra muestras sin carpeta de ensamblaje y las ordena por tamaño total ascendente.
    """
    info = []
    for sample, r1, r2 in pairs:
        outfolder = ASM_DIR / f"{sample}_assembled"
        if not outfolder.exists():
            sz = size_sum(r1, r2)
            info.append((sample, sz))
    return [m for m, _ in sorted(info, key=lambda x: x[1])]

# ==========
# Main
# ==========
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Ensamblaje metagenómico en serie (MEGAHIT/MetaSPAdes) con 02.trimmed/<sample>/**/R1,R2.")
    parser.add_argument("--assembler", choices=["megahit","metaspades"],
                        default=os.environ.get("MAG_ASSEMBLER", "megahit"),
                        help="Algoritmo de ensamblaje (default: megahit).")
    parser.add_argument("--threads", type=int,
                        default=int(os.environ.get("MAG_MEGAHIT_THREADS", os.environ.get("MAG_SPADES_THREADS", "40"))),
                        help="Hilos por muestra (default: env o 40).")
    parser.add_argument("--mem-gb", type=int,
                        default=int(os.environ.get("MAG_SPADES_MEM_GB", "0")),
                        help="Memoria (GB) para MetaSPAdes; 0=auto.")
    parser.add_argument("--debug", action="store_true",
                        help="Imprime información de depuración del descubrimiento de pares.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Si existe el directorio de salida de la muestra, lo elimina y vuelve a ensamblar.")
    args = parser.parse_args(argv)

    if not TRIM_DIR.exists():
        print(f"[ERROR] No existe carpeta de trimmed: {TRIM_DIR}", file=sys.stderr)
        return 2

    pairs = find_pairs(TRIM_DIR, debug=args.debug)
    if not pairs:
        print(f"[ERROR] No se detectaron pares R1/R2 dentro de subdirectorios en {TRIM_DIR}", file=sys.stderr)
        if not args.debug:
            print("Sugerencia: ejecute nuevamente con --debug para imprimir candidatos encontrados.", file=sys.stderr)
        return 3

    pendientes = obtener_muestras_pendientes(pairs)
    print(f" Muestras detectadas: {len(pairs)}")
    print(f" Muestras a ensamblar (ordenadas por tamaño): {len(pendientes)}")
    print(f" Iniciando ensamblaje en serie con {args.threads} hilos por muestra usando '{args.assembler}'...\n")

    idx = {s: (r1, r2) for s, r1, r2 in pairs}

    resultados = []
    for sample in pendientes:
        R1, R2 = idx[sample]
        print(f" Ensamblando: {sample}")
        res = assemble_one(sample, R1, R2, args.assembler, args.threads, args.mem_gb, overwrite=args.overwrite)
        resultados.append(res)
        print(f" Resultado: {res[1]} - {res[2]}\n")

    ASM_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESUMEN_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Muestra", "Estado", "Detalle"])
        writer.writerows(resultados)

    print(f"\n Logs por muestra en: {ASM_LOG_DIR}")
    print(f" Contigs generados en: {ASM_DIR}")
    print(f" Resumen CSV: {RESUMEN_CSV}")
    print("\n Ensamblajes completados.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

