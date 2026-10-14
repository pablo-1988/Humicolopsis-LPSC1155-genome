#!/usr/bin/env python3
"""
Extrae los marcadores propios de la cepa T1 (LPSC 1155) desde su genoma ensamblado,
para que T1 sea un terminal INDEPENDIENTE (no una copia de los depositos de CBS 645.76).

Estrategia:
 - Markers estandar (LSU, SSU, mtSSU, TUB2, EF1a, MCM7, RPB1, RPB2): blastn de la
   secuencia de referencia CBS 645.76 (deposito GenBank, en 01_query_T1) vs genoma T1,
   se extrae la region homologa del genoma de T1.
 - Genes solo-genoma (TUB1, RPA1, RPA2, RPC2, SF3B1, TFB4) e ITS: se extraen en el
   paso 02b con proteinas de referencia (tblastn).
"""
import os, subprocess, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.join(BASE, "_scratch")   # genoma de LPSC 1155 (T1_genome.fasta) y base BLAST
GENOME = os.path.join(SCRATCH, "T1_genome.fasta")
QDIR = os.path.join(BASE, "01_query_T1")
OUT = os.path.join(BASE, "02_ncbi_raw", "_T1")
os.makedirs(OUT, exist_ok=True)
DB = os.path.join(SCRATCH, "T1db")

# archivo query local (deposito CBS645) -> marcador destino
QMAP = {
    "SSU_CBS645.fasta":  "SSU",
    "ITS_CBS645.fasta":  "LSU",   # este archivo es 28S/LSU
    "mtSSU_CBS645.fasta":"mtSSU",
    "Btub_CBS645.fasta": "TUB2",
    "TEF1_CBS645.fasta": "EF1a",
    "MCM7_CBS645.fasta": "MCM7",
    "RPB1_CBS645.fasta": "RPB1",
    "RPB2_CBS645.fasta": "RPB2",
}

def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)

def main():
    if not os.path.exists(DB + ".nin"):
        print("makeblastdb genoma T1...")
        sh(f'makeblastdb -in "{GENOME}" -dbtype nucl -out "{DB}"')
    genome = SeqIO.to_dict(SeqIO.parse(GENOME, "fasta"))
    for qf, mk in QMAP.items():
        qpath = os.path.join(QDIR, qf)
        if not os.path.exists(qpath):
            print(f"{mk}: falta query {qf}"); continue
        # blastn query vs genoma T1, mejor hit
        r = sh(f'blastn -query "{qpath}" -db "{DB}" -outfmt "6 sseqid sstart send pident length evalue bitscore" -max_target_seqs 5 -evalue 1e-20')
        lines = [l for l in r.stdout.strip().split("\n") if l]
        if not lines:
            print(f"{mk}: sin hit en genoma T1"); continue
        # tomar el hit de mayor bitscore
        best = max(lines, key=lambda l: float(l.split("\t")[6]))
        sid, s, e, pid, length, ev, bs = best.split("\t")
        s, e = int(s), int(e)
        seq = genome[sid].seq
        if s <= e:
            sub = seq[s-1:e]
        else:
            sub = seq[e-1:s].reverse_complement()
        rec = SeqRecord(Seq(str(sub)), id=f"Humicolopsis_cephalosporioides_T1", description="")
        SeqIO.write([rec], os.path.join(OUT, f"T1__{mk}.fasta"), "fasta")
        print(f"{mk}: {sid}:{s}-{e}  {length}bp pident={pid}  ({len(sub)}bp extraido)")

if __name__ == "__main__":
    main()
