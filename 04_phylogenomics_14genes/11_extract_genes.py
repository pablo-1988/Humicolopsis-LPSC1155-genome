#!/usr/bin/env python3
"""Extrae 14 genes de cada genoma: 10 proteicos (blastp ref SwissProt -> proteoma ->
CDS nt) + 4 rDNA (blastn ref T1 -> genomic). Incluye T1 (proteoma/transcritos/rDNA local).
Versión F2: + H. cephalosporioides CH-1.1 (genomes14/, miniprot) y CBS 645.76 ex-tipo (taxón parcial:
EF1a, MCM7, RPB1, RPB2 CDS + SSU, LSU, ITS, mtSSU de GenBank)."""
import os, re, glob, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = os.path.dirname(os.path.abspath(__file__))
from pathlib import Path
os.environ["PATH"] = "/Users/pablo/.local/share/mamba/envs/filogenia/bin:" + os.environ["PATH"]
GDIR = str(Path.home() / "hongo_comparativa" / "genomes14")
REF = os.path.join(BASE, "refs")
OUT = os.path.join(BASE, "genes", "raw")
os.makedirs(OUT, exist_ok=True)
SCRATCH = str(Path.home() / "hongo_comparativa" / "tmp14"); os.makedirs(SCRATCH, exist_ok=True)
PROJ = Path(BASE).parents[2]
PRED = PROJ / "humicolopsis_annot/anotacion_final/predict_results"

PROT_GENES = ["TUB2","EF1a","MCM7","RPB1","RPB2","RPA1","RPA2","RPC2","SF3B1","TFB4"]
RDNA_GENES = ["SSU","LSU","mtSSU","ITS"]
seqs = {g: {} for g in PROT_GENES + RDNA_GENES}  # gene -> {taxon: SeqRecord(nt)}

def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True)

def cds_index(cds_path):
    """protein_id (o record.id) -> CDS nt."""
    idx = {}
    for r in SeqIO.parse(cds_path, "fasta"):
        m = re.search(r"\[protein_id=([^\]]+)\]", r.description)
        key = m.group(1) if m else r.id
        idx[key] = str(r.seq)
        idx.setdefault(r.id, str(r.seq))
    return idx

def extract_proteome(taxon, proteome, cds_path):
    db = os.path.join(SCRATCH, f"db_{taxon}")
    sh(f'makeblastdb -in "{proteome}" -dbtype prot -out "{db}"')
    idx = cds_index(cds_path)
    for g in PROT_GENES:
        r = sh(f'blastp -query "{REF}/{g}.faa" -db "{db}" -outfmt "6 sseqid pident length bitscore" -max_target_seqs 5 -evalue 1e-20')
        lines = [l for l in r.stdout.strip().split("\n") if l]
        if not lines: continue
        best = max(lines, key=lambda l: float(l.split("\t")[3]))
        sid = best.split("\t")[0]
        cds = idx.get(sid) or idx.get(sid.split("|")[-1])
        if not cds:
            # a veces sseqid trae prefijo lcl| o ref|
            for k in idx:
                if sid in k or k in sid: cds = idx[k]; break
        if cds and len(cds) >= 150:
            seqs[g][taxon] = SeqRecord(Seq(cds), id=taxon, description="")

def extract_rdna(taxon, genomic):
    db = os.path.join(SCRATCH, f"gdb_{taxon}")
    sh(f'makeblastdb -in "{genomic}" -dbtype nucl -out "{db}"')
    genome = None
    for g in RDNA_GENES:
        q = os.path.join(REF, f"T1__{g}.fasta")
        if not os.path.exists(q): continue
        r = sh(f'blastn -query "{q}" -db "{db}" -outfmt "6 sseqid sstart send pident length bitscore" -max_target_seqs 5 -evalue 1e-30')
        lines = [l for l in r.stdout.strip().split("\n") if l]
        if not lines: continue
        best = max(lines, key=lambda l: float(l.split("\t")[5]))
        sid,s,e = best.split("\t")[0], int(best.split("\t")[1]), int(best.split("\t")[2])
        if genome is None: genome = SeqIO.to_dict(SeqIO.parse(genomic, "fasta"))
        sub = genome[sid].seq[min(s,e)-1:max(s,e)]
        if s > e: sub = sub.reverse_complement()
        if len(sub) >= 200:
            seqs[g][taxon] = SeqRecord(Seq(str(sub)), id=taxon, description="")

# ── genomas NCBI ──
for d in sorted(glob.glob(os.path.join(GDIR, "*"))):
    if not os.path.isdir(d): continue
    taxon = os.path.basename(d)
    prot, cds, gen = [os.path.join(d,f) for f in ("protein.faa","cds.fna","genomic.fna")]
    print(f"### {taxon}")
    if os.path.exists(prot) and os.path.exists(cds):
        extract_proteome(taxon, prot, cds)
    if os.path.exists(gen):
        extract_rdna(taxon, gen)
    got = [g for g in seqs if taxon in seqs[g]]
    print(f"   {len(got)}/14: {' '.join(got)}")

# ── T1 (local) ──
print("### Humicolopsis_cephalosporioides_T1")
t1 = "Humicolopsis_cephalosporioides_T1"
import shutil  # makeblastdb no acepta rutas con espacios (Drive): copiar a SCRATCH
shutil.copy(PRED / "Humicolopsis_cephalosporioides_T1.proteins.fa", os.path.join(SCRATCH, "T1_proteins.fa"))
shutil.copy(PRED / "Humicolopsis_cephalosporioides_T1.cds-transcripts.fa", os.path.join(SCRATCH, "T1_cds.fa"))
extract_proteome(t1, os.path.join(SCRATCH, "T1_proteins.fa"), os.path.join(SCRATCH, "T1_cds.fa"))
for g in RDNA_GENES:
    q = os.path.join(REF, f"T1__{g}.fasta")
    if os.path.exists(q):
        r = list(SeqIO.parse(q,"fasta"))[0]
        seqs[g][t1] = SeqRecord(r.seq, id=t1, description="")
got=[g for g in seqs if t1 in seqs[g]]
print(f"   {len(got)}/14: {' '.join(got)}")

# ── CBS 645.76 (ex-tipo, sin genoma): marcadores GenBank ──
cbs = "Humicolopsis_cephalosporioides_CBS645.76"
for r in SeqIO.parse(PROJ / "Paper/analisis/F0_CBS645_marcadores/CBS645_cds.fna", "fasta"):
    seqs[r.id][cbs] = SeqRecord(r.seq, id=cbs, description="")
CACHE = PROJ / "filogenia_v2_revision/02_ncbi_raw/_cache"
for g in RDNA_GENES:
    f = CACHE / f"Humicolopsis_cephalosporioides__{g}.fasta"
    if f.exists() and f.stat().st_size > 0:
        r = next(SeqIO.parse(f, "fasta"))
        seqs[g][cbs] = SeqRecord(Seq(str(r.seq).replace("-", "").upper()), id=cbs, description="")
got = [g for g in seqs if cbs in seqs[g]]
print(f"### {cbs}\n   {len(got)}/14: {' '.join(got)}")

# ── escribir ──
for g in PROT_GENES + RDNA_GENES:
    recs = list(seqs[g].values())
    SeqIO.write(recs, os.path.join(OUT, f"{g}.fna"), "fasta")
print("\n=== cobertura por gen ===")
for g in PROT_GENES + RDNA_GENES:
    print(f"  {g:6s} {len(seqs[g])} taxa")
