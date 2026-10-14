"""Prepara H. cephalosporioides CH-1.1 (GCA_050941815.1, sin anotación) para el pipeline de 14 genes:
proteínas y CDS a partir del GFF de miniprot (mejor alineamiento por locus, sin solapamiento)."""
import os, shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
W = Path.home() / "hongo_comparativa"
D = W / "genomes14" / "Humicolopsis_cephalosporioides_CH-1.1"; D.mkdir(parents=True, exist_ok=True)
G = W / "fna" / "HumcCH11.fna"
genome = SeqIO.to_dict(SeqIO.parse(G, "fasta"))
mrna, cds, sta = {}, {}, None
for l in open(W / "miniprot" / "HumcCH11.gff"):
    if l.startswith("##STA"): sta = l.split("\t")[1].strip(); continue
    f = l.rstrip().split("\t")
    if len(f) < 9: continue
    if f[2] == "mRNA":
        i = f[8].split("ID=")[1].split(";")[0]
        mrna[i] = (f[0], int(f[3]), int(f[4]), f[6], float(f[5]), sta.replace("*", ""))
    elif f[2] == "CDS":
        cds.setdefault(f[8].split("Parent=")[1].split(";")[0], []).append((int(f[3]), int(f[4])))
order = sorted(mrna, key=lambda i: -mrna[i][4]); occ, kept = {}, []
for i in order:
    sq, s, e, st, sc, aa = mrna[i]
    if any(not (e < a or s > b) for a, b in occ.get((sq, st), [])): continue
    occ.setdefault((sq, st), []).append((s, e)); kept.append(i)
with open(D / "protein.faa", "w") as fp, open(D / "cds.fna", "w") as fc:
    for i in kept:
        sq, s, e, st, sc, aa = mrna[i]
        segs = sorted(cds[i])
        nt = "".join(str(genome[sq].seq[a - 1:b]) for a, b in segs)
        if st == "-": nt = str(Seq(nt).reverse_complement())
        fp.write(f">CH11_{i}\n{aa}\n"); fc.write(f">CH11_{i}_cds [protein_id=CH11_{i}]\n{nt}\n")
shutil.copy(G, D / "genomic.fna")
print(len(kept), "loci ->", D)
