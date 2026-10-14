#!/usr/bin/env python3
"""Extrae de T1: genes solo-genoma (TUB1,RPA1,RPA2,RPC2,SF3B1,TFB4) por tblastn con
proteinas de referencia de Neurospora, e ITS desde la region intergenica SSU-LSU."""
import os, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.join(BASE, "_scratch")   # genoma de LPSC 1155 (T1_genome.fasta) y base BLAST
GENOME = os.path.join(SCRATCH, "T1_genome.fasta")
DB = os.path.join(SCRATCH, "T1db")
REF = os.path.join(BASE, "02_ncbi_raw", "_refprot")
OUT = os.path.join(BASE, "02_ncbi_raw", "_T1")

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)

genome = SeqIO.to_dict(SeqIO.parse(GENOME, "fasta"))

# ── genes proteicos por tblastn ──
for g in ["TUB1", "RPA1", "RPA2", "RPC2", "SF3B1", "TFB4"]:
    q = os.path.join(REF, f"{g}.faa")
    if not os.path.exists(q):
        print(f"{g}: sin ref"); continue
    r = sh(f'tblastn -query "{q}" -db "{DB}" -outfmt "6 sseqid sstart send pident length evalue bitscore" -max_target_seqs 5 -evalue 1e-15')
    lines = [l for l in r.stdout.strip().split("\n") if l]
    if not lines:
        print(f"{g}: sin hit"); continue
    # agrupar HSPs por scaffold (mejor scaffold por suma de bitscore), luego rango min-max
    best_sid = max(set(l.split("\t")[0] for l in lines),
                   key=lambda s: sum(float(l.split("\t")[6]) for l in lines if l.split("\t")[0]==s))
    hsps = [l.split("\t") for l in lines if l.split("\t")[0]==best_sid]
    coords = []
    for h in hsps:
        s, e = int(h[1]), int(h[2]); coords += [s, e]
    lo, hi = min(coords), max(coords)
    strand_minus = sum(1 for h in hsps if int(h[1])>int(h[2])) > len(hsps)/2
    if hi - lo > 8000:  # evitar rango disparatado
        # usar solo el HSP de mayor bitscore
        h = max(hsps, key=lambda x: float(x[6])); s,e=int(h[1]),int(h[2]); lo,hi=min(s,e),max(s,e); strand_minus = s>e
    sub = genome[best_sid].seq[lo-1:hi]
    if strand_minus:
        sub = sub.reverse_complement()
    rec = SeqRecord(Seq(str(sub)), id="Humicolopsis_cephalosporioides_T1", description="")
    SeqIO.write([rec], os.path.join(OUT, f"T1__{g}.fasta"), "fasta")
    print(f"{g}: {best_sid}:{lo}-{hi} {'(-)' if strand_minus else '(+)'} {len(sub)}bp")

# ── ITS: region intergenica SSU-LSU en scaffold 12 ──
sc = "JBPDIK010000012.1_RagTag"
if sc in genome:
    its = genome[sc].seq[1211696-1:1212888].reverse_complement()  # entre LSU(fin) y SSU(inicio)
    rec = SeqRecord(Seq(str(its)), id="Humicolopsis_cephalosporioides_T1", description="")
    SeqIO.write([rec], os.path.join(OUT, "T1__ITS.fasta"), "fasta")
    print(f"ITS: {sc}:1211696-1212888 (-) {len(its)}bp")
