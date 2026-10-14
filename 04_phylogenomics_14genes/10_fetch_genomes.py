#!/usr/bin/env python3
"""Descarga (NCBI datasets) proteoma, CDS y genoma de los 16 genomas del análisis de 14 genes
+ prepara H. cephalosporioides CH-1.1 (GCA_050941815.1; proteínas/CDS por miniprot).
Genomas en carpeta local persistente (no Drive, no /tmp)."""
import os, subprocess, zipfile, glob, shutil
from pathlib import Path
GDIR = str(Path.home() / "hongo_comparativa" / "genomes14")
os.makedirs(GDIR, exist_ok=True)
DS = "/Users/pablo/.local/share/mamba/envs/filogenia/bin/datasets"

ASSEMBLIES = {
 "Neurospora_crassa":            ("GCF_000182925.2", "Sordariales"),   # OR74A / Neucr2 (outgroup)
 "Chaetomium_globosum":          ("GCF_000143365.1", "Sordariales"),
 "Botrytis_cinerea":             ("GCF_000143535.2", "Helotiales"),
 "Sclerotinia_sclerotiorum":     ("GCF_000146945.2", "Helotiales"),
 "Glarea_lozoyensis":            ("GCF_000409485.1", "Helotiales"),
 "Amorphotheca_resinae":         ("GCF_003019875.1", "Helotiales"),
 "Mollisia_scopiformis":         ("GCF_001500285.1", "Helotiales"),
 "Cadophora_sp":                 ("GCA_030347535.1", "Helotiales"),
 "Chlorociboria_aeruginascens":  ("GCA_002276475.2", "Helotiales"),
 "Hymenoscyphus_fraxineus":      ("GCA_911649665.1", "Helotiales"),
 "Rhynchosporium_commune":       ("GCA_046529785.1", "Helotiales"),
 "Oidiodendron_maius":           ("GCA_000827325.1", "Myxotrichaceae"),
 "Blumeria_graminis":            ("GCA_905067625.1", "Erysiphales"),
 "Erysiphe_necator":             ("GCA_024703715.1", "Erysiphales"),
 "Erysiphe_pulchra":             ("GCA_002918395.1", "Erysiphales"),
 "Pseudogymnoascus_australis":   ("GCA_048593815.1", "Thelebolales"),
}

for tax, (acc, order) in ASSEMBLIES.items():
    d = os.path.join(GDIR, tax); os.makedirs(d, exist_ok=True)
    if all(os.path.exists(os.path.join(d, f)) for f in ("protein.faa", "cds.fna", "genomic.fna")):
        print(f"{tax:30s} cache"); continue
    z = os.path.join(d, "tmp.zip")
    for attempt in range(4):
        subprocess.run([DS, "download", "genome", "accession", acc, "--include", "protein,cds,genome", "--filename", z],
                       capture_output=True)
        if os.path.exists(z) and zipfile.is_zipfile(z): break
    else:
        print(f"{tax:30s} FALLO descarga {acc}"); continue
    with zipfile.ZipFile(z) as zf: zf.extractall(d)
    base = glob.glob(os.path.join(d, "ncbi_dataset/data/GC*"))[0]
    for src, dst in (("protein.faa", "protein.faa"), ("cds_from_genomic.fna", "cds.fna")):
        if os.path.exists(os.path.join(base, src)): shutil.move(os.path.join(base, src), os.path.join(d, dst))
    shutil.move(glob.glob(os.path.join(base, "*_genomic.fna"))[0], os.path.join(d, "genomic.fna"))
    shutil.rmtree(os.path.join(d, "ncbi_dataset")); os.remove(z)
    print(f"{tax:30s} {acc}  " + " ".join(f"{f}:{os.path.exists(os.path.join(d, f))}" for f in ("protein.faa", "cds.fna", "genomic.fna")))
