#!/usr/bin/env python3
"""Alinea (codon-consciente para genes proteicos, nt para rDNA), concatena e infiere
el arbol filogenomico de 14 genes. IQ-TREE2 MFP+MERGE + UFBoot 1000."""
import os, glob, subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "genes", "raw")
ALN = os.path.join(BASE, "align"); TREE = os.path.join(BASE, "tree")
CONC = os.path.join(BASE, "align", "concat")
for d in (ALN, TREE, CONC): os.makedirs(d, exist_ok=True)

PROT = ["TUB2","EF1a","MCM7","RPB1","RPB2","RPA1","RPA2","RPC2","SF3B1","TFB4"]
RDNA = ["SSU","LSU","mtSSU","ITS"]
MIN_TAXA = 4; MIN_MARK_PER_TAXON = 3

def env(t): return f"mamba run -n filogenia {t}"
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True)

def mafft(inp, out, aa=False):
    sh(f'{env("mafft")} --auto --anysymbol "{inp}" > "{out}" 2>/dev/null')

def codon_align(gene, recs):
    """translate -> MAFFT AA -> backtranslate a codones."""
    nt = {}; aas = []
    for r in recs:
        s = str(r.seq).upper().replace("-","")
        s = s[: len(s)//3*3]
        if len(s) < 150: continue
        aa = str(Seq(s).translate(table=1)).replace("*","X")
        nt[r.id] = s
        aas.append(SeqRecord(Seq(aa), id=r.id, description=""))
    if len(aas) < MIN_TAXA: return None
    fa = os.path.join(ALN, f"{gene}_aa.fasta"); al = os.path.join(ALN, f"{gene}_aa_aln.fasta")
    SeqIO.write(aas, fa, "fasta"); mafft(fa, al, aa=True)
    out = {}
    for r in SeqIO.parse(al, "fasta"):
        s = nt[r.id]; i = 0; codons = []
        for a in str(r.seq):
            if a == "-": codons.append("---")
            else: codons.append(s[i:i+3] if i+3 <= len(s) else "NNN"); i += 3
        out[r.id] = "".join(codons)
    return out

def nt_align(gene, recs):
    recs = [r for r in recs if len(str(r.seq).replace("-","")) >= 200]
    if len(recs) < MIN_TAXA: return None
    fa = os.path.join(ALN, f"{gene}_nt.fasta"); al = os.path.join(ALN, f"{gene}_nt_aln.fasta")
    SeqIO.write(recs, fa, "fasta"); mafft(fa, al)
    return {r.id: str(r.seq) for r in SeqIO.parse(al, "fasta")}

def trimal(gene, aln_dict):
    fa = os.path.join(ALN, f"{gene}_pre.fasta"); tr = os.path.join(ALN, f"{gene}_trim.fasta")
    with open(fa,"w") as fh:
        for k,v in aln_dict.items(): fh.write(f">{k}\n{v}\n")
    sh(f'{env("trimal")} -in "{fa}" -out "{tr}" -automated1')
    recs = list(SeqIO.parse(tr, "fasta"))
    if not recs or len(recs[0].seq) < 60: return None
    return {r.id: str(r.seq) for r in recs}

def main():
    kept = {}
    for g in PROT + RDNA:
        f = os.path.join(RAW, f"{g}.fna")
        if not os.path.exists(f): continue
        recs = list(SeqIO.parse(f, "fasta"))
        aln = codon_align(g, recs) if g in PROT else nt_align(g, recs)
        if not aln: print(f"[skip] {g}"); continue
        tr = trimal(g, aln)
        if tr: kept[g] = tr; print(f"[ok] {g}: {len(tr)} taxa, {len(next(iter(tr.values())))} bp")
    # concat
    taxa = set().union(*[set(kept[g]) for g in kept])
    counts = {t: sum(t in kept[g] for g in kept) for t in taxa}
    taxa = sorted(t for t in taxa if counts[t] >= MIN_MARK_PER_TAXON or t.startswith("Humicolopsis"))
    order = [g for g in PROT+RDNA if g in kept]
    L = {g: len(next(iter(kept[g].values()))) for g in order}
    fa = os.path.join(CONC, "concat.fasta")
    with open(fa,"w") as fh:
        for t in taxa:
            fh.write(">"+t+"\n"+"".join(kept[g].get(t,"-"*L[g]) for g in order)+"\n")
    pos=1
    with open(os.path.join(CONC,"partition.txt"),"w") as fh:
        for g in order: fh.write(f"DNA, {g} = {pos}-{pos+L[g]-1}\n"); pos+=L[g]
    with open(os.path.join(CONC,"taxa_matrix.tsv"),"w") as fh:
        fh.write("taxon\t"+"\t".join(order)+"\tn\n")
        for t in taxa:
            row=["1" if t in kept[g] else "0" for g in order]
            fh.write(t+"\t"+"\t".join(row)+f"\t{sum(map(int,row))}\n")
    print(f"\nMatriz: {len(taxa)} taxa x {pos-1} bp, {len(order)} genes")
    pre=os.path.join(TREE,"Humicolopsis_phylogenomics")
    print("IQ-TREE...")
    r=sh(f'{env("iqtree2")} -s "{fa}" -p "{os.path.join(CONC,"partition.txt")}" -m MFP+MERGE -B 1000 -T AUTO --prefix "{pre}" -redo')
    print(r.stdout[-800:] if r.stdout else r.stderr[-800:])

if __name__ == "__main__":
    main()
