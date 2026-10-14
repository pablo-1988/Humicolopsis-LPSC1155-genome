#!/usr/bin/env python3
"""
Colector NCBI para la filogenia v2 (revision del colega).
Para cada (especie x marcador) busca en GenBank prefiriendo material tipo,
descarga la mejor secuencia y arma un FASTA por marcador + matriz de cobertura.

15 marcadores de Johnston et al. 2019:
  rDNA: SSU, ITS(5.8S), LSU  | Sanger multigen: mtSSU, TUB2, TUB1, EF1a, MCM7, RPB1, RPB2
  solo-genoma: RPA1, RPA2, RPC2, SF3B1, TFB4
"""
import os, re, time, sys
from collections import defaultdict
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.org"   # NCBI requiere un correo de contacto
BASE = os.path.dirname(os.path.abspath(__file__))
TAXTAB = os.path.join(BASE, "00_design", "taxon_sampling.tsv")
OUT = os.path.join(BASE, "02_ncbi_raw")
CACHE = os.path.join(OUT, "_cache")
os.makedirs(CACHE, exist_ok=True)

# marcador -> (terminos Title, len_min, len_max)
MARKERS = {
    "SSU":  ('("small subunit ribosomal"[Title] OR "18S ribosomal"[Title] OR "18S rRNA"[Title]) NOT mitochond*[Title]', 500, 2500),
    "LSU":  ('("large subunit ribosomal"[Title] OR "28S ribosomal"[Title] OR "28S rRNA"[Title] OR "LSU"[Title]) NOT mitochond*[Title]', 400, 2500),
    "ITS":  ('("internal transcribed spacer"[Title] OR "5.8S ribosomal"[Title] OR "ITS1"[Title])', 300, 1200),
    "mtSSU":('("small subunit ribosomal"[Title] AND mitochond*[Title]) OR "mtSSU"[Title]', 300, 1500),
    "TUB2": ('("beta-tubulin"[Title] OR "beta tubulin"[Title] OR "tub2"[Title] OR "TUB2"[Title] OR "benA"[Title])', 250, 1800),
    "TUB1": ('("alpha-tubulin"[Title] OR "alpha tubulin"[Title] OR "tub1"[Title] OR "tubA"[Title])', 250, 1800),
    "EF1a": ('("elongation factor 1-alpha"[Title] OR "elongation factor 1 alpha"[Title] OR "TEF1"[Title] OR "tef1"[Title] OR "EF-1alpha"[Title] OR "EF1-alpha"[Title])', 200, 1800),
    "MCM7": ('("MCM7"[Title] OR "minichromosome maintenance"[Title] OR "DNA replication licensing"[Title])', 250, 1800),
    "RPB1": ('("RPB1"[Title] OR "rpb1"[Title] OR "polymerase II largest"[Title] OR "polymerase II subunit 1"[Title])', 300, 2200),
    "RPB2": ('("RPB2"[Title] OR "rpb2"[Title] OR "polymerase II second"[Title] OR "polymerase II subunit 2"[Title])', 300, 2200),
    "RPA1": ('("RPA1"[Title] OR "replication protein A 70"[Title] OR "replication protein A1"[Title])', 300, 2500),
    "RPA2": ('("RPA2"[Title] OR "replication protein A 32"[Title] OR "replication protein A2"[Title])', 250, 2000),
    "RPC2": ('("RPC2"[Title] OR "polymerase III 127"[Title] OR "polymerase III second"[Title])', 300, 2500),
    "SF3B1":('("SF3B1"[Title] OR "splicing factor 3B"[Title])', 300, 3000),
    "TFB4": ('("TFB4"[Title] OR "transcription and DNA repair factor IIH subunit"[Title] OR "GTF2H3"[Title])', 250, 2500),
}

def read_taxa():
    taxa = []
    for line in open(TAXTAB):
        line = line.rstrip("\n")
        if not line or line.startswith("#") or line.startswith("order\t"):
            continue
        p = line.split("\t")
        if len(p) < 5:
            continue
        taxa.append({"order": p[0], "family": p[1], "species": p[2], "note": p[3], "role": p[4]})
    return taxa

def esearch(term, retmax=30):
    for attempt in range(3):
        try:
            h = Entrez.esearch(db="nucleotide", term=term, retmax=retmax)
            r = Entrez.read(h); h.close()
            time.sleep(0.35)
            return r["IdList"]
        except Exception as e:
            time.sleep(1.5)
    return []

def efetch(ids):
    if not ids:
        return []
    for attempt in range(3):
        try:
            h = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            recs = list(SeqIO.parse(h, "genbank")); h.close()
            time.sleep(0.4)
            return recs
        except Exception as e:
            time.sleep(1.5)
    return []

def is_type(rec):
    txt = (rec.description + " " + str(rec.annotations)).lower()
    for feat in rec.features:
        for k, v in feat.qualifiers.items():
            if "type_material" in k or any("type material" in str(x).lower() for x in v):
                return True
    return "type material" in txt or "type strain" in txt

def pick_best(recs, len_min, len_max):
    """Prefiere material tipo, luego mayor longitud dentro de rango."""
    cands = [r for r in recs if len_min <= len(r.seq) <= len_max]
    if not cands:
        cands = recs
    if not cands:
        return None
    cands.sort(key=lambda r: (is_type(r), len(r.seq)), reverse=True)
    return cands[0]

def collect():
    taxa = read_taxa()
    by_marker = defaultdict(list)          # marker -> [(species, rec)]
    coverage = defaultdict(dict)           # species -> marker -> accession
    for t in taxa:
        sp = t["species"]
        if t["role"] == "focal" and "T1" in sp:
            continue  # T1 se agrega desde archivos locales / genoma
        sp_clean = re.sub(r"\s+", "_", sp)
        org = sp.replace("-", " ")
        # el nombre de organismo para NCBI (quitar sufijos)
        org_q = re.sub(r"\bcf\.\b|\baff\.\b", "", org).strip()
        print(f"\n### {sp}  [{t['order']}]")
        for mk, (terms, lmin, lmax) in MARKERS.items():
            cache_f = os.path.join(CACHE, f"{sp_clean}__{mk}.fasta")
            if os.path.exists(cache_f):
                r = list(SeqIO.parse(cache_f, "fasta"))
                if r:
                    by_marker[mk].append((sp_clean, r[0]))
                    coverage[sp_clean][mk] = r[0].id.split("|")[-1]
                    print(f"   {mk:5s} cache {r[0].id.split('|')[-1]}")
                    continue
                else:
                    print(f"   {mk:5s} cache (vacio)")
                    continue
            term = f'"{org_q}"[Organism] AND ({terms})'
            ids = esearch(term, retmax=25)
            if not ids:
                open(cache_f, "w").close()  # marca vacio
                print(f"   {mk:5s} --- (0 hits)")
                continue
            recs = efetch(ids[:12])
            best = pick_best(recs, lmin, lmax)
            if best is None:
                open(cache_f, "w").close()
                print(f"   {mk:5s} --- (sin secuencia usable)")
                continue
            acc = best.id
            typ = "TYPE" if is_type(best) else ""
            best.id = f"{sp_clean}|{acc}"
            best.description = ""
            SeqIO.write([best], cache_f, "fasta")
            by_marker[mk].append((sp_clean, best))
            coverage[sp_clean][mk] = acc
            print(f"   {mk:5s} {acc} {len(best.seq)}bp {typ}")
    return by_marker, coverage, taxa

def write_outputs(by_marker, coverage, taxa):
    # FASTA por marcador
    for mk, items in by_marker.items():
        with open(os.path.join(OUT, f"{mk}.fasta"), "w") as fh:
            for sp, rec in items:
                rec.id = sp
                rec.description = ""
                SeqIO.write([rec], fh, "fasta")
    # matriz de cobertura
    mks = list(MARKERS.keys())
    with open(os.path.join(OUT, "coverage_matrix.tsv"), "w") as fh:
        fh.write("species\torder\t" + "\t".join(mks) + "\tn_markers\n")
        for t in taxa:
            if "T1" in t["species"]:
                continue
            sp = re.sub(r"\s+", "_", t["species"])
            cov = coverage.get(sp, {})
            row = [("1" if mk in cov else "0") for mk in mks]
            n = sum(int(x) for x in row)
            fh.write(f"{sp}\t{t['order']}\t" + "\t".join(row) + f"\t{n}\n")
    print("\n=== FASTA y coverage_matrix.tsv escritos en 02_ncbi_raw/ ===")

if __name__ == "__main__":
    bm, cov, taxa = collect()
    write_outputs(bm, cov, taxa)
