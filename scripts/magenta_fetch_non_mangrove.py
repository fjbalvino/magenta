#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MAGENTA: Fetch non-mangrove coastal metagenomes (WGS + Illumina)
+ Enriquecimiento y filtro geográfico:
  - ENA read_run: lat/lon/location -> latitude/longitude/geographic_location
  - ENA read_sample (plan B): lat/lon/location por sample_accession
  - Parseo de geographic_location (N/S/E/W y 'lat, lon')
  - Queries SRA sin comodines explosivos
  - Filtro final: coordenadas válidas
"""

import os, sys, shlex, subprocess, csv, re, io
from pathlib import Path
from urllib.parse import urlencode, quote_plus

import pandas as pd
import numpy as np
import requests
from tqdm import tqdm
import certifi

# ---------------------------------------------------------------------
# CONFIG Y DIRECTORIOS
# ---------------------------------------------------------------------
os.environ.setdefault("MAGENTA_DIR", "/nfs/testing/.jbalvino/MAGENTA/MAGENTA_DIR")

PROJECT_DIR = Path(os.environ.get("MAGENTA_DIR", Path.cwd())).resolve()
RAW_DIR  = PROJECT_DIR / "rawdata" / "fastq"
META_DIR = PROJECT_DIR / "metadata"
OUT_DIR  = PROJECT_DIR / "outputs"
LOGS_DIR = PROJECT_DIR / "logs"
AUX_DIR  = PROJECT_DIR / "aux"
for d in [RAW_DIR, META_DIR, OUT_DIR, LOGS_DIR, AUX_DIR]:
    d.mkdir(parents=True, exist_ok=True)

META_BASE = "metadatos_enriquecidos"
META_FILE = OUT_DIR / f"{META_BASE}.csv"

def save_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)

def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str)

def ensure_float(df: pd.DataFrame, cols):
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def km_to_deg(km: float) -> float:
    return km / 111.32

print("PROJECT_DIR:", PROJECT_DIR)

# ---------------------------------------------------------------------
# ENTORNO SEGURO PARA curl/EDirect Y requests
# ---------------------------------------------------------------------
CA_BUNDLE = certifi.where()
SAFE_ENV = dict(os.environ)
SAFE_ENV["CURL_CA_BUNDLE"] = CA_BUNDLE
SAFE_ENV["SSL_CERT_FILE"]  = CA_BUNDLE
SAFE_ENV.setdefault("TERM", "dumb")
SAFE_ENV.setdefault("EUTILS_TOOL", "magenta_edirect")
SAFE_ENV.setdefault("EUTILS_EMAIL", "jbalvino@masternew")
os.environ["REQUESTS_CA_BUNDLE"] = CA_BUNDLE

# ---------------------------------------------------------------------
# QUERIES (SRA/EDirect) — WGS + ILLUMINA, sin wildcard-explosion
# ---------------------------------------------------------------------
SRA_QUERIES = {
    "salt_marsh": (
        '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
        '('
        ' ( "salt marsh"[All Fields] OR marsh[All Fields] OR marshes[All Fields] '
        '   OR "Spartina alterniflora"[All Fields] OR cordgrass[All Fields] ) '
        ' AND ( sediment[All Fields] OR sediments[All Fields] '
        '       OR rhizosphere[All Fields] OR rhizospheric[All Fields] '
        '       OR root[All Fields] OR roots[All Fields] OR "root-associated"[All Fields] ) '
        ' AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) '
        ') NOT ( mangrove[All Fields] OR mangroves[All Fields] )'
    ),

    "seagrass": (
        '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
        '('
        ' ( seagrass[All Fields] OR seagrasses[All Fields] '
        '   OR "Zostera"[All Fields] OR "Zostera marina"[All Fields] '
        '   OR "Posidonia"[All Fields] OR "Thalassia"[All Fields] ) '
        ' AND ( sediment[All Fields] OR sediments[All Fields] '
        '       OR rhizosphere[All Fields] OR rhizospheric[All Fields] '
        '       OR root[All Fields] OR roots[All Fields] OR "root-associated"[All Fields] ) '
        ' AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) '
        ') NOT ( mangrove[All Fields] OR mangroves[All Fields] )'
    ),

    "estuarine": (
        '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
        '('
        ' ( estuary[All Fields] OR estuarine[All Fields] OR "tidal flat"[All Fields] '
        '   OR "tidal flats"[All Fields] OR mudflat[All Fields] OR mudflats[All Fields] '
        '   OR intertidal[All Fields] ) '
        ' AND ( sediment[All Fields] OR sediments[All Fields] ) '
        ' AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) '
        ') NOT ( mangrove[All Fields] OR mangroves[All Fields] )'
    ),

    "microbial_mats": (
        '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
        '('
        ' ( "microbial mat"[All Fields] OR "microbial mats"[All Fields] '
        '   OR hypersaline[All Fields] OR "solar saltern"[All Fields] OR salterns[All Fields] ) '
        ' AND ( intertidal[All Fields] OR coastal[All Fields] OR lagoon[All Fields] OR lagoons[All Fields] OR sabkha[All Fields] ) '
        ' AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) '
        ') NOT ( mangrove[All Fields] OR mangroves[All Fields] )'
    ),

    "halophyte_rhizo": (
        '("WGS"[Strategy] AND "ILLUMINA"[Platform]) AND '
        '('
        ' ( halophyte[All Fields] OR halophytes[All Fields] OR Salicornia[All Fields] OR Suaeda[All Fields] ) '
        ' AND ( rhizosphere[All Fields] OR rhizospheric[All Fields] '
        '       OR root[All Fields] OR roots[All Fields] OR "root-associated"[All Fields] ) '
        ' AND ( metagenome[All Fields] OR metagenomic[All Fields] OR metagenomics[All Fields] ) '
        ') NOT ( mangrove[All Fields] OR mangroves[All Fields] )'
    ),
}

# ---------------------------------------------------------------------
# ENA QUERIES — WGS + ILLUMINA + NOT mangrove (read_run)
#   Campos de coordenadas en ENA read_run: lat, lon, location
#   Plataforma: instrument_platform
# ---------------------------------------------------------------------
def ena_terms_for(ecosystem: str) -> str:
    if ecosystem == "salt_marsh":
        return '( "salt marsh" OR marshes OR "Spartina alterniflora" OR cordgrass ) AND (sediment OR sediments OR rhizosphere OR rhizospheric OR root OR roots OR "root-associated")'
    if ecosystem == "seagrass":
        return '( seagrass OR seagrasses OR "Zostera" OR "Zostera marina" OR "Posidonia" OR "Thalassia" ) AND (sediment OR sediments OR rhizosphere OR rhizospheric OR root OR roots OR "root-associated")'
    if ecosystem == "estuarine":
        return '( estuary OR estuarine OR "tidal flat" OR "tidal flats" OR mudflat OR mudflats OR intertidal ) AND (sediment OR sediments)'
    if ecosystem == "microbial_mats":
        return '( "microbial mat" OR "microbial mats" OR hypersaline OR "solar saltern" OR salterns ) AND (intertidal OR coastal OR lagoon OR lagoons OR sabkha)'
    if ecosystem == "halophyte_rhizo":
        return '( halophyte OR halophytes OR Salicornia OR Suaeda ) AND (rhizosphere OR rhizospheric OR root OR roots OR "root-associated")'
    return ''

def ena_query_for(ecosystem: str) -> str:
    terms = ena_terms_for(ecosystem)
    return f'(library_strategy="WGS" AND instrument_platform="ILLUMINA" AND ({terms})) NOT mangrove'

ENA_FIELDS = [
    "run_accession","sample_accession","study_accession",
    "library_strategy","library_layout","instrument_platform",
    "collection_date","country",
    "location","lat","lon",
    "fastq_ftp","fastq_http","fastq_md5"
]
ENA_BASE = "https://www.ebi.ac.uk/ena/portal/api/search"

# ---------------------------------------------------------------------
# SUBPROCESOS Y RETRIES
# ---------------------------------------------------------------------
def run_cmd(cmd: str, retries: int = 3) -> int:
    print(">>", cmd)
    last = 1
    for i in range(1, retries+1):
        last = subprocess.call(cmd, shell=True, env=SAFE_ENV)
        if last == 0:
            return 0
        print(f"[run_cmd][retry {i}/{retries}] exit={last}")
    return last

def have_edirect() -> bool:
    for tool in ("esearch", "efetch"):
        code = subprocess.call(f"command -v {tool} >/dev/null 2>&1", shell=True, env=SAFE_ENV)
        if code != 0:
            return False
    return True

# ---------------------------------------------------------------------
# NORMALIZACIÓN
# ---------------------------------------------------------------------
def postfilter_wgs_illumina_csv(csv_in: Path, csv_out: Path):
    """Filtra RunInfo por LibraryStrategy=WGS y Platform=ILLUMINA (case-insensitive)."""
    if not csv_in.exists() or csv_in.stat().st_size == 0:
        return
    with open(csv_in, newline='', encoding='utf-8', errors="ignore") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        if not rows:
            csv_out.write_text("") ; return
        fields = reader.fieldnames
    def ok(row):
        ls = (row.get("LibraryStrategy") or "").upper()
        pf = (row.get("Platform") or "").upper()
        return ("WGS" in ls) and ("ILLUMINA" in pf)
    kept = [r for r in rows if ok(r)]
    with open(csv_out, "w", newline='', encoding="utf-8") as g:
        w = csv.DictWriter(g, fieldnames=fields)
        w.writeheader()
        for r in kept:
            w.writerow(r)

def normalize_ena(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns={
        "run_accession":      "Run",
        "sample_accession":   "Sample",
        "study_accession":    "Study",
        "library_layout":     "LibraryLayout",
        "library_strategy":   "LibraryStrategy",
        "instrument_platform":"Platform",
        "lat":                "latitude",
        "lon":                "longitude",
        "location":           "geographic_location",
    })
    return df

def normalize_sra(df: pd.DataFrame) -> pd.DataFrame:
    ren = {}
    if "lat" in df.columns: ren["lat"] = "latitude"
    if "lon" in df.columns: ren["lon"] = "longitude"
    return df.rename(columns=ren)

# ---------------------------------------------------------------------
# GEO: ENRIQUECIMIENTO Y PARSEO
# ---------------------------------------------------------------------
def _parse_ns_ew(val: str):
    """Interpreta N/S/E/W: 23.5S -> -23.5, 100.1W -> -100.1. Devuelve (lat, lon) si encuentra dos números."""
    if not isinstance(val, str):
        return None, None
    tokens = re.findall(r'([+-]?\d+(?:\.\d+)?)([NnSsEeWw]?)', val.strip())
    nums = []
    for num, hemi in tokens:
        x = float(num)
        if hemi.upper() == 'S':
            x = -abs(x)
        elif hemi.upper() == 'N':
            x = abs(x)
        elif hemi.upper() == 'W':
            x = -abs(x)
        elif hemi.upper() == 'E':
            x = abs(x)
        nums.append(x)
    if len(nums) >= 2:
        return nums[0], nums[1]
    return None, None

def _extract_lat_lon_from_text(s: str):
    if not isinstance(s, str):
        return np.nan, np.nan
    lat, lon = _parse_ns_ew(s)
    if lat is not None and lon is not None:
        return lat, lon
    # patrón genérico "lat, lon"
    nums = re.findall(r'([+-]?\d+(?:\.\d+)?)', s)
    if len(nums) >= 2:
        return float(nums[0]), float(nums[1])
    return np.nan, np.nan

def _valid_lat_lon(lat, lon):
    try:
        if pd.isna(lat) or pd.isna(lon):
            return False
        lat = float(lat); lon = float(lon)
        return np.isfinite(lat) and np.isfinite(lon) and -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0
    except Exception:
        return False

def enrich_coords_from_ena_by_runs(run_ids, chunk_size=200):
    """Devuelve DataFrame con columnas: Run, latitude, longitude, geographic_location (desde ENA read_run)."""
    records = []
    run_ids = [r for r in run_ids if isinstance(r, str) and r.strip()]
    for i in range(0, len(run_ids), chunk_size):
        chunk = run_ids[i:i+chunk_size]
        ors = " OR ".join([f'run_accession="{r}"' for r in chunk])
        params = {
            "result": "read_run",
            "query":  ors,
            "fields": ",".join(["run_accession", "lat", "lon", "location"]),
            "format": "tsv",
            "limit":  0,
        }
        url = f"{ENA_BASE}?{urlencode(params)}"
        try:
            r = requests.get(url, timeout=180, verify=CA_BUNDLE)
            r.raise_for_status()
            if r.content:
                dfc = pd.read_csv(io.StringIO(r.content.decode("utf-8")), sep="\t", dtype=str)
                if not dfc.empty:
                    dfc = dfc.rename(columns={
                        "run_accession": "Run",
                        "lat": "latitude",
                        "lon": "longitude",
                        "location": "geographic_location",
                    })
                    records.append(dfc)
        except Exception as e:
            print(f"[ENA][enrich run][WARN] chunk {i//chunk_size}: {e}")

    if not records:
        return pd.DataFrame(columns=["Run","latitude","longitude","geographic_location"])
    return pd.concat(records, ignore_index=True, sort=False)

def enrich_coords_from_ena_by_samples(sample_ids, chunk_size=200):
    """Plan B: devuelve DataFrame con columnas: Sample, latitude, longitude, geographic_location (desde ENA read_sample)."""
    records = []
    sample_ids = [s for s in sample_ids if isinstance(s, str) and s.strip()]
    for i in range(0, len(sample_ids), chunk_size):
        chunk = sample_ids[i:i+chunk_size]
        ors = " OR ".join([f'sample_accession="{s}"' for s in chunk])
        params = {
            "result": "sample",
            "query":  ors,
            "fields": ",".join(["sample_accession","lat","lon","location","country"]),
            "format": "tsv",
            "limit":  0,
        }
        url = f"{ENA_BASE}?{urlencode(params)}"
        try:
            r = requests.get(url, timeout=180, verify=CA_BUNDLE)
            r.raise_for_status()
            if r.content:
                dfc = pd.read_csv(io.StringIO(r.content.decode("utf-8")), sep="\t", dtype=str)
                if not dfc.empty:
                    dfc = dfc.rename(columns={
                        "sample_accession": "Sample",
                        "lat": "latitude",
                        "lon": "longitude",
                        "location": "geographic_location",
                    })
                    records.append(dfc)
        except Exception as e:
            print(f"[ENA][enrich sample][WARN] chunk {i//chunk_size}: {e}")

    if not records:
        return pd.DataFrame(columns=["Sample","latitude","longitude","geographic_location"])
    return pd.concat(records, ignore_index=True, sort=False)

def fill_and_filter_geo(all_df: pd.DataFrame) -> pd.DataFrame:
    """
    Rellena lat/lon faltantes:
      1) ENA read_run por Run
      2) ENA read_sample por Sample
      3) Parseo de geographic_location
    Luego filtra filas con coordenadas válidas.
    """
    # Asegura columnas
    for col in ["latitude","longitude","geographic_location","Sample","Run"]:
        if col not in all_df.columns:
            all_df[col] = np.nan

    # --- (1) Completar con ENA por Run ---
    need_geo = all_df["Run"].notna() & (all_df["latitude"].isna() | all_df["longitude"].isna())
    runs_to_enrich = all_df.loc[need_geo, "Run"].dropna().astype(str).unique().tolist()
    if runs_to_enrich:
        print(f"[GEO] Enriqueciendo coordenadas desde ENA read_run por Run (n={len(runs_to_enrich)})…")
        enr = enrich_coords_from_ena_by_runs(runs_to_enrich, chunk_size=200)
        if not enr.empty:
            all_df = all_df.merge(enr[["Run","latitude","longitude","geographic_location"]],
                                  on="Run", how="left", suffixes=("", "_ena"))
            for col in ["latitude","longitude","geographic_location"]:
                all_df[col] = all_df[col].fillna(all_df[f"{col}_ena"])
                if f"{col}_ena" in all_df.columns:
                    all_df.drop(columns=[f"{col}_ena"], inplace=True)

    # --- (2) Completar con ENA por Sample (plan B) ---
    still_missing = (all_df["latitude"].isna() | all_df["longitude"].isna()) & all_df["Sample"].notna()
    samples_to_enrich = all_df.loc[still_missing, "Sample"].dropna().astype(str).unique().tolist()
    if samples_to_enrich:
        print(f"[GEO] Enriqueciendo coordenadas desde ENA sample por Sample (n={len(samples_to_enrich)})…")
        ens = enrich_coords_from_ena_by_samples(samples_to_enrich, chunk_size=200)
        if not ens.empty:
            all_df = all_df.merge(ens[["Sample","latitude","longitude","geographic_location"]],
                                  on="Sample", how="left", suffixes=("", "_sena"))
            for col in ["latitude","longitude","geographic_location"]:
                all_df[col] = all_df[col].fillna(all_df[f"{col}_sena"])
                if f"{col}_sena" in all_df.columns:
                    all_df.drop(columns=[f"{col}_sena"], inplace=True)

    # --- (3) Parseo de geographic_location ---
    still_missing2 = all_df["latitude"].isna() | all_df["longitude"].isna()
    if still_missing2.any():
        print(f"[GEO] Parseando geographic_location para {int(still_missing2.sum())} filas…")
        latlon = all_df.loc[still_missing2, "geographic_location"].apply(_extract_lat_lon_from_text)
        lat_parsed = latlon.map(lambda t: t[0])
        lon_parsed = latlon.map(lambda t: t[1])
        all_df.loc[still_missing2, "latitude"]  = all_df.loc[still_missing2, "latitude"].fillna(lat_parsed)
        all_df.loc[still_missing2, "longitude"] = all_df.loc[still_missing2, "longitude"].fillna(lon_parsed)

    # Convertir a numérico y validar
    all_df["latitude"]  = pd.to_numeric(all_df["latitude"], errors="coerce")
    all_df["longitude"] = pd.to_numeric(all_df["longitude"], errors="coerce")
    valid_mask = all_df.apply(lambda r: _valid_lat_lon(r["latitude"], r["longitude"]), axis=1)
    before_n = len(all_df)
    all_df = all_df[valid_mask].copy()
    after_n  = len(all_df)
    print(f"[GEO] Filtrado por coordenadas válidas: {before_n} -> {after_n}")

    return all_df

# ---------------------------------------------------------------------
# FETCHERS
# ---------------------------------------------------------------------
def fetch_ena_for_ecosystem(ecosystem: str) -> Path:
    params = {
        "result": "read_run",
        "query":  ena_query_for(ecosystem),
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit":  0
    }
    url = f"{ENA_BASE}?{urlencode(params)}"
    print(f"[ENA] {ecosystem} -> {url}")
    out_tsv = META_DIR / f"ena_{ecosystem}.tsv"
    r = requests.get(url, timeout=180, verify=CA_BUNDLE)
    r.raise_for_status()
    with open(out_tsv, "wb") as f:
        f.write(r.content)
    try:
        with open(out_tsv, 'r', encoding='utf-8', errors='ignore') as fr:
            n = sum(1 for _ in fr) - 1
        print(f"[ENA] filas -> {max(n,0)} :: {out_tsv}")
    except Exception:
        pass
    return out_tsv

def sra_runinfo_fallback(query: str, out_csv: Path) -> bool:
    """Fallback: endpoint CSV de RunInfo (sin WebEnv)."""
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?term=" + quote_plus(query)
    print("[SRA][fallback] GET", url)
    try:
        r = requests.get(url, timeout=180, verify=CA_BUNDLE)
        r.raise_for_status()
        txt = r.text
        if not txt or "Run,ReleaseDate,LoadDate" not in txt.splitlines()[0]:
            print("[SRA][fallback][WARN] Respuesta no luce como RunInfo CSV.")
            return False
        out_csv.write_text(txt, encoding="utf-8")
        return True
    except Exception as e:
        print("[SRA][fallback][ERROR]", e)
        return False

def fetch_sra_for_ecosystem(ecosystem: str, query: str) -> Path:
    raw_csv = META_DIR / f"sra_{ecosystem}_raw_runinfo.csv"
    out_csv = META_DIR / f"sra_{ecosystem}_wgs_illumina.csv"

    cmd = 'esearch -db sra -query {q} | efetch -format runinfo > {out}'.format(
        q=shlex.quote(query),
        out=shlex.quote(str(raw_csv))
    )
    print(f"[SRA] {ecosystem} (WGS+ILLUMINA) — EDirect")
    ret = run_cmd(cmd)
    need_fallback = (ret != 0) or (not raw_csv.exists()) or (raw_csv.stat().st_size == 0)

    if need_fallback:
        print(f"[SRA][WARN] EDirect falló o sin filas para {ecosystem}. Intentando fallback RunInfo CSV…")
        if sra_runinfo_fallback(query, raw_csv):
            print(f"[SRA][fallback] OK -> {raw_csv}")
        else:
            print(f"[SRA][fallback] Sin datos para {ecosystem}.")

    postfilter_wgs_illumina_csv(raw_csv, out_csv)
    kept = 0
    if out_csv.exists() and out_csv.stat().st_size > 0:
        try:
            with open(out_csv, "r", encoding="utf-8", errors="ignore") as fr:
                kept = sum(1 for _ in fr) - 1
        except Exception:
            pass
    print(f"[SRA] {ecosystem} -> {max(kept,0)} filas (post-filter) :: {out_csv}")
    return out_csv

# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------
def main():
    print("\n==> Iniciando descarga por ecosistema (ENA + SRA)")
    have_ed = have_edirect()
    if not have_ed:
        print("[WARN] EDirect no está en PATH o no es ejecutable. Se usará solo ENA y el fallback de RunInfo.")

    ecosystems = list(SRA_QUERIES.keys())

    # ENA
    ena_paths = []
    for eco in tqdm(ecosystems, desc="ENA", unit="eco"):
        try:
            ena_paths.append(fetch_ena_for_ecosystem(eco))
        except Exception as e:
            print(f"[ENA][ERROR] {eco}: {e}")

    # SRA
    sra_paths = []
    for eco in tqdm(ecosystems, desc="SRA", unit="eco"):
        try:
            if have_ed:
                sra_paths.append(fetch_sra_for_ecosystem(eco, SRA_QUERIES[eco]))
            else:
                raw_csv = META_DIR / f"sra_{eco}_raw_runinfo.csv"
                if sra_runinfo_fallback(SRA_QUERIES[eco], raw_csv):
                    out_csv = META_DIR / f"sra_{eco}_wgs_illumina.csv"
                    postfilter_wgs_illumina_csv(raw_csv, out_csv)
                    sra_paths.append(out_csv)
                    print(f"[SRA][no-edirect] {eco} -> {out_csv}")
                else:
                    print(f"[SRA][no-edirect] Sin datos para {eco}.")
        except Exception as e:
            print(f"[SRA][ERROR] {eco}: {e}")

    # -----------------------------------------------------------------
    # NORMALIZAR Y UNIFICAR
    # -----------------------------------------------------------------
    print("\n==> Unificando metadatos ENA + SRA")
    frames = []

    # ENA
    for p in ena_paths:
        if p.exists() and p.stat().st_size > 0:
            try:
                df = pd.read_csv(p, sep="\t", dtype=str)
                df = normalize_ena(df)
                df["source"] = "ENA"
                df["ecosystem"] = p.stem.replace("ena_", "")
                frames.append(df)
            except Exception as e:
                print(f"[ENA][READ][WARN] {p.name}: {e}")

    # SRA
    for p in sra_paths:
        if p.exists() and p.stat().st_size > 0:
            try:
                df = pd.read_csv(p, dtype=str)
                df = normalize_sra(df)
                df["source"] = "SRA"
                df["ecosystem"] = (p.stem
                                   .replace("sra_", "")
                                   .replace("_wgs_illumina","")
                                   .replace("_raw_runinfo",""))
                frames.append(df)
            except Exception as e:
                print(f"[SRA][READ][WARN] {p.name}: {e}")

    if not frames:
        print("[WARN] No se obtuvieron tablas para unificar.")
        return

    all_df = pd.concat(frames, ignore_index=True, sort=False)

    # Columnas mínimas
    for col in ["Run","Sample","Study","LibraryLayout","LibraryStrategy","Platform",
                "latitude","longitude","geographic_location","ecosystem","source"]:
        if col not in all_df.columns:
            all_df[col] = np.nan

    # Filtro final WGS + ILLUMINA
    mask = (
        all_df["LibraryStrategy"].fillna("").str.upper().str.contains("WGS") &
        all_df["Platform"].fillna("").str.upper().str.contains("ILLUMINA")
    )
    all_df = all_df[mask].copy()

    # Dedupe por Run
    all_df["Run"] = all_df["Run"].astype(str)
    all_df = all_df.drop_duplicates(subset=["Run"], keep="first")

    # -----------------------------------------------------------------
    # ENRIQUECER COORDENADAS Y FILTRAR GEO VÁLIDO
    # -----------------------------------------------------------------
    print("\n==> Enriqueciendo/validando coordenadas (read_run → read_sample → parseo)")
    all_df = fill_and_filter_geo(all_df)
    if all_df.empty:
        print("[WARN] Tras filtrar por coordenadas válidas no quedaron corridas.")
        # Aun así, escribe archivos vacíos para trazabilidad
        for eco in ["salt_marsh","seagrass","estuarine","microbial_mats","halophyte_rhizo"]:
            out_eco = OUT_DIR / f"{eco}_metadata.tsv"
            save_tsv(pd.DataFrame(), out_eco)
        unified_tsv = META_DIR / "metadatos_unificados.tsv"
        save_tsv(pd.DataFrame(), unified_tsv)
        return

    # -----------------------------------------------------------------
    # EXPORTAR RESULTADOS
    # -----------------------------------------------------------------
    # Salvar por ecosistema
    for eco in sorted(all_df["ecosystem"].dropna().unique()):
        df_eco = all_df[all_df["ecosystem"] == eco].copy()
        out_eco = OUT_DIR / f"{eco}_metadata.tsv"
        save_tsv(df_eco, out_eco)
        print(f"[OK] {eco}: {len(df_eco)} filas -> {out_eco}")

    # Unificado total
    unified_tsv = META_DIR / "metadatos_unificados.tsv"
    save_tsv(all_df, unified_tsv)
    print(f"\n[OK] Metadatos unificados (geo-válidos) -> {unified_tsv}  ({len(all_df)} corridas)")

    # Resumen
    summary = (
        all_df.groupby(["ecosystem","source"])["Run"]
        .nunique()
        .reset_index()
        .rename(columns={"Run":"n_runs"})
        .sort_values(["ecosystem","source"])
    )
    log_path = LOGS_DIR / "summary_non_mangrove_wgs_illumina.tsv"
    save_tsv(summary, log_path)
    print(f"[OK] Resumen -> {log_path}")

# ---------------------------------------------------------------------
# ENTRYPOINT
# ---------------------------------------------------------------------
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrumpido por usuario.")
        sys.exit(130)

