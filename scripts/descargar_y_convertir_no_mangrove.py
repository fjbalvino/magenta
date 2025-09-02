#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
descargar_y_convertir.py
Lee selected_metagenomes_full.tsv (salida de selección) y:
- Si hay .fastq.gz (ENA) → descarga directa (aria2c o wget).
- Si no, intenta NCBI/SRA: fasterq-dump remoto por Run; o convierte .sra/.lite.1 si están enlazados.
- Soporta CA corporativa, exclusiones por EXCLUDE_RUNS, y normaliza nombres *_1/_2.
"""

import os
import sys
import re
import subprocess
import datetime
from pathlib import Path
from urllib.parse import unquote
from shutil import which
import pandas as pd

# ---------------------------------------------------------------------
# Rutas fijas del proyecto
# ---------------------------------------------------------------------
PROJECT_DIR = Path("/nfs/testing/.jbalvino/MAGENTA/MAGENTA_DIR").resolve()

# ⚠️ Cambiado: ahora leemos selected_metagenomes_full.tsv
META_PATH   = PROJECT_DIR / "metadata" / "selected_metagenomes_full.tsv"

ENA_DIR     = PROJECT_DIR / "rawdata" / "ena"
CONV_DIR    = PROJECT_DIR / "rawdata" / "convertidos"
TMP_FQD     = PROJECT_DIR / "tmp_fqd"
LOGS_DIR    = PROJECT_DIR / "logs"
AUX_DIR     = PROJECT_DIR / "aux"
SRA_IMG     = "docker://ncbi/sra-tools:3.1.1"  # fallback en contenedor (si existe)

for d in (ENA_DIR, CONV_DIR, TMP_FQD, LOGS_DIR, AUX_DIR):
    d.mkdir(parents=True, exist_ok=True)

if not META_PATH.exists():
    print(f"[ERROR] No se encuentra {META_PATH}")
    sys.exit(1)

RUNLOG = (LOGS_DIR / f"descarga_all_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log").open("w")

def logprint(*args):
    msg = " ".join(str(a) for a in args)
    print(msg)
    try:
        RUNLOG.write(msg + "\n"); RUNLOG.flush()
    except Exception:
        pass

def run_cmd(cmd_list, check=True, env=None):
    """Ejecuta comando mostrando stdout/stderr en el log. Devuelve True/False."""
    logprint("[CMD]", " ".join(cmd_list))
    try:
        res = subprocess.run(
            cmd_list, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env
        )
        out = (res.stdout.decode(errors="ignore") + res.stderr.decode(errors="ignore")).strip()
        if out:
            logprint("[OUT]", out)
        return True
    except subprocess.CalledProcessError as e:
        out = ((e.stdout.decode(errors="ignore") if e.stdout else "") +
               (e.stderr.decode(errors="ignore") if e.stderr else "")).strip()
        if out:
            logprint("[ERR]", out)
        return False

# ---------------------------------------------------------------------
# EXCLUSIONES (opcional)
# ---------------------------------------------------------------------
EXCLUDE_RUNS = set(r.strip() for r in os.environ.get("EXCLUDE_RUNS", "").split(",") if r.strip())

# ---------------------------------------------------------------------
# CA corporativa / TLS
# ---------------------------------------------------------------------
def apply_corporate_ca():
    corp_ca = os.environ.get("CORP_CA", "").strip()
    sys_cas = [
        "/etc/pki/tls/certs/ca-bundle.crt",
        "/etc/ssl/certs/ca-certificates.crt",
        "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem",
    ]
    vdbc = which("vdb-config")
    def vdb_set(k, v):
        if vdbc:
            run_cmd([vdbc, "-s", f"{k}={v}"], check=False)

    if vdbc:
        run_cmd([vdbc, "--restore-defaults"], check=False)

    if corp_ca and Path(corp_ca).exists():
        os.environ["SSL_CERT_FILE"] = corp_ca
        os.environ["REQUESTS_CA_BUNDLE"] = corp_ca
        os.environ["CURL_CA_BUNDLE"] = corp_ca
        vdb_set("/tls/verify-peer", "true")
        vdb_set("/tls/use-system-ca-cert", "false")
        vdb_set("/tls/ca-file", corp_ca)
        logprint(f"[TLS] CA corporativa aplicada: {corp_ca}")
        return corp_ca

    for p in sys_cas:
        if Path(p).exists():
            os.environ.setdefault("SSL_CERT_FILE", p)
            os.environ.setdefault("REQUESTS_CA_BUNDLE", p)
            os.environ.setdefault("CURL_CA_BUNDLE", p)
            vdb_set("/tls/verify-peer", "true")
            vdb_set("/tls/use-system-ca-cert", "yes")
            vdb_set("/tls/ca-file", p)
            logprint(f"[TLS] CA del sistema aplicada: {p}")
            return p

    vdb_set("/tls/verify-peer", "true")
    vdb_set("/tls/use-system-ca-cert", "yes")
    logprint("[TLS] Sin CA explícita; configuración por defecto.")
    return None

CA_IN_USE = apply_corporate_ca()

def wget_base_args():
    base = [which("wget") or "wget", "-c", "--tries=10", "--read-timeout=30"]
    if os.environ.get("ALLOW_INSECURE_WGET") == "1":
        base.insert(1, "--no-check-certificate")
    elif CA_IN_USE and Path(CA_IN_USE).exists():
        base.extend(["--ca-certificate", CA_IN_USE])
    return base

# ---------------------------------------------------------------------
# Utilidades generales
# ---------------------------------------------------------------------
def sra_tools_version():
    fqd = which("fasterq-dump")
    if not fqd:
        return (0, 0, 0)
    ok = subprocess.run([fqd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    m = re.search(r"fasterq-dump\.(\d+)\.(\d+)\.(\d+)", ok.stdout.decode())
    if m:
        return tuple(int(x) for x in m.groups())
    return (0, 0, 0)

SRA_VER = sra_tools_version()
logprint(f"[SRA-TOOLS] fasterq-dump versión detectada: {'.'.join(map(str,SRA_VER))}")

def parse_semicolon_urls(s: str):
    if not isinstance(s, str):
        return []
    return [u.strip() for u in s.split(";") if u and u.strip()]

def list_fastq_for_run(base_dir: Path, run_id: str):
    outs = []
    for root, _, files in os.walk(base_dir):
        for f in files:
            if f.endswith((".fastq", ".fastq.gz")) and (f.startswith(run_id) or f.startswith(f"{run_id}.sra")):
                outs.append(Path(root) / f)
    return sorted(outs)

def have_apptainer():
    return bool(which("apptainer") or which("singularity"))

def digits_after_prefix(run_id: str) -> str:
    # p.ej. SRR29002871 -> "29002871"
    return re.sub(r"^[A-Za-z]+", "", run_id)

def ena_dir2(run_id: str) -> str | None:
    # Regla ENA: 7 dígitos -> "00x"; 8 -> "0xx"; 9+ -> "xxx"; <=6 -> sin subcarpeta
    d = digits_after_prefix(run_id)
    L = len(d)
    if L <= 6: return None
    if L == 7:  return f"00{d[-1]}"
    if L == 8:  return f"0{d[-2:]}"
    return d[-3:]  # 9+

def ena_fastq_candidates(run_id: str):
    pre6 = run_id[:6]
    d2 = ena_dir2(run_id)
    parts = []
    if d2:
        parts.append((pre6, d2, run_id))
    parts.append((pre6, run_id))  # fallback sin d2
    urls = []
    for scheme in ("https", "http", "ftp"):
        for p in parts:
            base = f"{scheme}://ftp.sra.ebi.ac.uk/vol1/fastq/" + "/".join(p)
            urls += [f"{base}/{run_id}_1.fastq.gz", f"{base}/{run_id}_2.fastq.gz"]
    return urls

def ena_sra_candidates(run_id: str):
    pre6 = run_id[:6]
    d2 = ena_dir2(run_id)
    tuples = []
    if d2:
        tuples.append((pre6, d2, run_id))
    tuples.append((pre6, run_id))
    urls = []
    for scheme in ("https", "http", "ftp"):
        for p in tuples:
            base = f"{scheme}://ftp.sra.ebi.ac.uk/vol1/srr/" + "/".join(p)
            urls.append(f"{base}/{run_id}.sra")
    return urls

def wget_spider(url: str) -> bool:
    cmd = wget_base_args() + ["--spider", url]
    ok = run_cmd(cmd, check=False)
    return ok

def download_from_url(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        logprint(f"[DL] Ya existe: {dest.name}")
        return True
    cmd = wget_base_args() + [url, "-O", str(dest)]
    ok = run_cmd(cmd, check=False)
    return ok and dest.exists() and dest.stat().st_size > 0

def normalize_fastq_names(run_id: str):
    """Ajusta nombres post-dump: SRR.fastq -> SRR_1.fastq; gestiona single-end y opcional _2 vacío."""
    outs = list_fastq_for_run(CONV_DIR, run_id)
    if not outs:
        return

    # Caso común con input "archivo.sra": queda "SRRxxxxxx.sra.fastq"
    for p in list(outs):
        if p.name.endswith(".sra.fastq"):
            newp = p.with_name(p.name.replace(".sra.fastq", ".fastq"))
            p.rename(newp)
            logprint(f"[FIX] Renombrado {p.name} -> {newp.name}")
    outs = list_fastq_for_run(CONV_DIR, run_id)

    # Si hay exactamente 1 archivo y no es *_1/_2, renombrar a _1
    only = [p for p in outs if p.name == f"{run_id}.fastq" or p.name == f"{run_id}.fastq.gz"]
    if len(outs) == 1 and only:
        p = outs[0]
        suf = (p.suffixes[-2] + p.suffixes[-1]) if "".join(p.suffixes).endswith(".fastq.gz") else p.suffix
        target = p.with_name(f"{run_id}_1{suf}")
        p.rename(target)
        logprint(f"[FIX] Dataset single-end: {p.name} -> {target.name}")
        if os.environ.get("TOUCH_EMPTY_R2") == "1":
            r2 = target.with_name(f"{run_id}_2{suf}")
            if r2.suffix == ".gz":
                run_cmd(["bash", "-lc", f": | gzip -c > {r2}"], check=False)
            else:
                r2.touch()
            logprint(f"[FIX] Creado _2 vacío: {r2.name}")

# ---------------------------------------------------------------------
# Carga metadatos y detección de columnas
# ---------------------------------------------------------------------
logprint(f"[INFO] Cargando {META_PATH}")
df = pd.read_csv(META_PATH, sep="\t", dtype=str, low_memory=False).fillna("")
logprint(f"[INFO] Filas: {len(df):,}")

# Soportar nombres variados de columnas (mayúsculas/minúsculas/espacios)
cols_lower = {c.lower(): c for c in df.columns}
COL_RUN = (
    cols_lower.get("run") or
    cols_lower.get("run_accession") or
    cols_lower.get("run id") or
    cols_lower.get("run_id")
)

# Enlaces potenciales (coinciden con salidas de tu fetcher y selecciones)
LINK_KEYS = ("ftp_link", "download_path", "fastq_ftp", "fastq_http", "sra_link", "sra_ftp", "sra_http")
COL_LINKS = [cols_lower[k] for k in LINK_KEYS if k in cols_lower]

if not COL_RUN:
    logprint("[ERROR] No se encontró columna Run/run_accession en el TSV.")
    sys.exit(1)

def row_get_links(row):
    links = []
    for col in (COL_LINKS or []):
        links += parse_semicolon_urls(str(row.get(col, "")))
    if not links:
        # Plan C: si alguien dejó URLs separadas por ';' en otra columna
        big = " ".join(map(str, row.values))
        if ";" in big:
            links = parse_semicolon_urls(big)
    return links

def is_ena_row(row):
    return any(".fastq.gz" in u.lower() for u in row_get_links(row))

def is_ncbi_row(row):
    L = [u.lower() for u in row_get_links(row)]
    return any(".lite.1" in u or u.endswith(".sra") for u in L)

# ---------------------------------------------------------------------
# ENA directa (si hay)
# ---------------------------------------------------------------------
def download_ena_row(row):
    run = str(row.get(COL_RUN, "")).strip()
    urls = [u for u in row_get_links(row) if ".fastq.gz" in u.lower()]
    if not urls:
        # generar candidatos estándar ENA si no venían en el TSV
        urls = ena_fastq_candidates(run)

    got = 0
    seen = set()
    aria2 = which("aria2c")
    if aria2 and urls:
        lst = AUX_DIR / f"ena_{run}.txt"
        with open(lst, "w") as f:
            for u in urls:
                if u not in seen:
                    f.write(u + "\n"); seen.add(u)
        cmd = [aria2, "--continue=true",
               "--max-concurrent-downloads=4", "--max-connection-per-server=8",
               "--retry-wait=5", "--max-tries=10", "--file-allocation=none",
               f"--input-file={str(lst)}", f"--dir={str(ENA_DIR)}"]
        run_cmd(cmd, check=False)
        for p in ENA_DIR.glob(f"{run}*.fastq.gz"):
            logprint(f"[ENA] descargado: {p.name}"); got += 1
    else:
        for u in urls:
            if u in seen: continue
            seen.add(u)
            fname = ENA_DIR / os.path.basename(unquote(u))
            if wget_spider(u) and download_from_url(u, fname):
                logprint(f"[ENA] descargado: {fname.name}")
                got += 1
    return got

# ---------------------------------------------------------------------
# Cadena NCBI (igual que el original corregido)
# ---------------------------------------------------------------------
def fasterq_dump_remote(run_id: str) -> bool:
    fqd = which("fasterq-dump")
    if not fqd:
        return False
    threads = os.environ.get("THREADS", "8")
    cmd = [fqd, run_id, "--split-files", "-O", str(CONV_DIR),
           "--threads", threads, "--temp", str(TMP_FQD)]
    ok = run_cmd(cmd, check=False, env=os.environ.copy())
    normalize_fastq_names(run_id)
    outs = list_fastq_for_run(CONV_DIR, run_id)
    return bool(outs)

def fasterq_dump_local(sra_path: Path, run_id: str) -> bool:
    fqd = which("fasterq-dump")
    if not fqd:
        return False
    threads = os.environ.get("THREADS", "8")
    cmd = [fqd, str(sra_path), "--outdir", str(CONV_DIR),
           "--split-files", "--threads", threads, "--temp", str(TMP_FQD)]
    ok = run_cmd(cmd, check=False)
    normalize_fastq_names(run_id)
    outs = list_fastq_for_run(CONV_DIR, run_id)
    if outs and os.environ.get("KEEP_SRA", "0") != "1":
        try: sra_path.unlink(); logprint(f"[NCBI] Eliminado {sra_path.name}")
        except Exception: pass
    return bool(outs)

def fastq_dump_local(sra_path: Path, run_id: str) -> bool:
    fd = which("fastq-dump")
    if not fd:
        return False
    cmd = [fd, str(sra_path), "--outdir", str(CONV_DIR),
           "--split-files", "--gzip", "--skip-technical", "--readids", "--dumpbase", "--clip"]
    ok = run_cmd(cmd, check=False)
    normalize_fastq_names(run_id)
    outs = list_fastq_for_run(CONV_DIR, run_id)
    if outs and os.environ.get("KEEP_SRA", "0") != "1":
        try: sra_path.unlink(); logprint(f"[NCBI] Eliminado {sra_path.name}")
        except Exception: pass
    return bool(outs)

def fasterq_dump_container(sra_path: Path, run_id: str) -> bool:
    runner = which("apptainer") or which("singularity")
    if not runner:
        return False
    bind_arg = f"{str(PROJECT_DIR)}:/work"
    in_path  = f"/work/{sra_path.relative_to(PROJECT_DIR)}"
    out_path = f"/work/{CONV_DIR.relative_to(PROJECT_DIR)}"
    tmp_path = f"/work/{TMP_FQD.relative_to(PROJECT_DIR)}"
    threads = os.environ.get("THREADS", "8")
    cmd = [runner, "exec", "--cleanenv", "--bind", bind_arg,
           SRA_IMG, "fasterq-dump", in_path,
           "--outdir", out_path, "--split-files", "--threads", threads, "--temp", tmp_path]
    ok = run_cmd(cmd, check=False)
    normalize_fastq_names(run_id)
    outs = list_fastq_for_run(CONV_DIR, run_id)
    if outs and os.environ.get("KEEP_SRA", "0") != "1":
        try: sra_path.unlink(); logprint(f"[NCBI] Eliminado {sra_path.name}")
        except Exception: pass
    return bool(outs)

def download_sra_from_ena(run_id: str) -> Path | None:
    for u in ena_sra_candidates(run_id):
        if wget_spider(u):
            dest = CONV_DIR / f"{run_id}.sra"
            if download_from_url(u, dest):
                logprint(f"[ENA-SRA] descargado: {dest.name}")
                return dest
    logprint("[ENA-SRA] No se encontró .sra en ENA.")
    return None

def download_fastq_from_ena(run_id: str) -> int:
    got = 0
    for u in ena_fastq_candidates(run_id):
        dest = CONV_DIR / os.path.basename(unquote(u))
        if wget_spider(u) and download_from_url(u, dest):
            logprint(f"[ENA-FASTQ] descargado: {dest.name}")
            got += 1
    return got

def download_ncbi_row_and_convert(row):
    run = str(row.get(COL_RUN, "")).strip()
    links = row_get_links(row)

    logprint(f"[NCBI] Intento remoto por accesión: {run}")
    if fasterq_dump_remote(run):
        for p in list_fastq_for_run(CONV_DIR, run):
            logprint(f"[NCBI] generado (remoto): {p.name}")
        return

    # Preferir SRA completo si el TSV lo trae
    sra_link  = next((u for u in links if u.lower().endswith(".sra")), None)
    lite_link = next((u for u in links if u.lower().endswith(".lite.1")), None)
    sra_path = None

    if sra_link:
        base = os.path.basename(sra_link)
        sra_path = CONV_DIR / base
        if not sra_path.exists():
            logprint(f"[NCBI] Descargando archivo del TSV: {base}")
            if not download_from_url(sra_link, sra_path):
                logprint("[NCBI] Descarga .sra desde NCBI falló.")
                sra_path = None
        else:
            logprint(f"[NCBI] Ya existe archivo: {sra_path.name}")

    if sra_path and sra_path.exists() and sra_path.stat().st_size > 0:
        if fasterq_dump_local(sra_path, run) or fastq_dump_local(sra_path, run) or fasterq_dump_container(sra_path, run):
            for p in list_fastq_for_run(CONV_DIR, run):
                logprint(f"[NCBI] generado (.sra): {p.name}")
            return
        logprint("[NCBI] Conversión falló con .sra del TSV].")

    # ENA: primero FASTQ (si existen), luego SRA
    logprint("[NCBI] Intento obtener FASTQ directos desde ENA.")
    if download_fastq_from_ena(run) > 0:
        for p in list_fastq_for_run(CONV_DIR, run):
            logprint(f"[NCBI] generado (ENA FASTQ): {p.name}")
        return

    logprint("[NCBI] Intento obtener .sra completo desde ENA.")
    ena_sra = download_sra_from_ena(run)
    if ena_sra and ena_sra.exists() and ena_sra.stat().st_size > 0:
        if fasterq_dump_local(ena_sra, run) or fastq_dump_local(ena_sra, run) or fasterq_dump_container(ena_sra, run):
            for p in list_fastq_for_run(CONV_DIR, run):
                logprint(f"[NCBI] generado (ENA .sra): {p.name}")
            return

    # Último recurso: .lite.1 + contenedor 3.x
    if lite_link:
        base = os.path.basename(lite_link)
        lite_path = CONV_DIR / base
        if not lite_path.exists():
            logprint(f"[NCBI] Descargando archivo .lite.1 del TSV: {base}")
            if not download_from_url(lite_link, lite_path):
                lite_path = None
        if lite_path and lite_path.exists() and lite_path.stat().st_size > 0:
            if fasterq_dump_container(lite_path, run):
                for p in list_fastq_for_run(CONV_DIR, run):
                    logprint(f"[NCBI] generado (.lite.1 → contenedor 3.x): {p.name}")
                return

    logprint("[NCBI] No fue posible generar FASTQ para la corrida.")

# ---------------------------------------------------------------------
# EJECUCIÓN PARA TODAS LAS FILAS
# ---------------------------------------------------------------------
procesadas = set()
ok_ena = ok_ncbi = 0

for i, row in df.iterrows():
    run = str(row.get(COL_RUN, "")).strip()
    if not run:
        continue
    if run in EXCLUDE_RUNS:
        logprint(f"[SKIP] EXCLUDE_RUNS: {run}")
        continue
    if run in procesadas:
        continue
    procesadas.add(run)

    try:
        if is_ena_row(row):
            logprint(f"\n[PROCESS][ENA] idx={i} Run={run}")
            got = download_ena_row(row)
            ok_ena += 1 if got > 0 else 0
            # Si el TSV trae .fastq.gz ya no intentamos NCBI para evitar duplicados
            continue

        if is_ncbi_row(row):
            logprint(f"\n[PROCESS][NCBI] idx={i} Run={run}")
            download_ncbi_row_and_convert(row)
            ok_ncbi += 1
        else:
            logprint(f"[WARN] Fila sin enlaces útiles: idx={i} Run={run}")

    except Exception as e:
        logprint(f"[ERROR] Fallo procesando idx={i} Run={run}: {e}")

# ---------------------------------------------------------------------
# Verificaciones simples
# ---------------------------------------------------------------------
try:
    ena_files = [f for f in os.listdir(ENA_DIR) if f.endswith(".fastq.gz")]
    print(f"[CHK][ENA] archivos (directos): {len(ena_files)}")
except Exception as e:
    print("[CHK][ENA] error:", e)

try:
    conv_files = [f for f in os.listdir(CONV_DIR) if f.endswith(".fastq") or f.endswith(".fastq.gz")]
    print(f"[CHK][NCBI] archivos: {len(conv_files)}")
except Exception as e:
    print("[CHK][NCBI] error:", e)

print(f"Procesadas ENA: {ok_ena}, NCBI: {ok_ncbi}, Total RUNs únicos: {len(procesadas)}")
print("Proceso ENA + NCBI para TODAS las muestras finalizado.")

