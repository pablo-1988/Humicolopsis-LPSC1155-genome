#!/usr/bin/env python3
"""
Alinea (MAFFT) -> recorta (TrimAl) -> concatena -> IQ-TREE2, para la filogenia v2.
Arbol principal: 9 marcadores estandar y de amplio muestreo.
Fuente de secuencias: 02_ncbi_raw/_cache/  (NCBI) + 02_ncbi_raw/_T1/ (genoma T1).
Outgroup: Neurospora crassa (como Quijada et al. 2022).
"""
import os, glob, subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "02_ncbi_raw", "_cache")
T1DIR = os.path.join(BASE, "02_ncbi_raw", "_T1")
ALN = os.path.join(BASE, "03_alignments")
TRIM = os.path.join(BASE, "04_trimmed")
CONC = os.path.join(BASE, "05_concatenated")
TREE = os.path.join(BASE, "06_tree")
for d in (ALN, TRIM, CONC, TREE): os.makedirs(d, exist_ok=True)

def env(tool): return f"mamba run -n filogenia {tool}"
def sh(cmd): return subprocess.run(cmd, shell=True, capture_output=True, text=True)

# marcador -> longitud maxima permitida (filtra mRNA/genomicos gigantes)
MARKERS = {"SSU":2500,"LSU":2500,"ITS":1200,"mtSSU":1600,
           "TUB2":1800,"EF1a":1800,"MCM7":1800,"RPB1":2400,"RPB2":2400}
MIN_TAXA_PER_MARKER = 6   # descarta particiones casi vacias
MIN_MARKERS_PER_TAXON = 2 # un taxon necesita >=2 marcadores

def gather():
    """marker -> {taxon: SeqRecord}"""
    data = defaultdict(dict)
    for mk, lmax in MARKERS.items():
        # NCBI cache
        for f in glob.glob(os.path.join(CACHE, f"*__{mk}.fasta")):
            if os.path.getsize(f) == 0: continue
            recs = list(SeqIO.parse(f, "fasta"))
            if not recs: continue
            r = recs[0]
            sp = r.id.split("|")[0]
            seq = str(r.seq).replace("-", "").upper()
            if len(seq) > lmax or len(seq) < 150: continue
            data[mk][sp] = SeqRecord(Seq(seq), id=sp, description="")
        # T1 (genoma)
        t1f = os.path.join(T1DIR, f"T1__{mk}.fasta")
        if os.path.exists(t1f):
            r = list(SeqIO.parse(t1f, "fasta"))[0]
            seq = str(r.seq).replace("-", "").upper()
            if len(seq) <= lmax:
                data[mk]["Humicolopsis_cephalosporioides_T1"] = SeqRecord(Seq(seq), id="Humicolopsis_cephalosporioides_T1", description="")
    return data

def align_trim(data):
    kept = {}
    for mk in MARKERS:
        recs = list(data[mk].values())
        if len(recs) < MIN_TAXA_PER_MARKER:
            print(f"[skip] {mk}: solo {len(recs)} taxa (<{MIN_TAXA_PER_MARKER})"); continue
        raw = os.path.join(ALN, f"{mk}_raw.fasta")
        SeqIO.write(recs, raw, "fasta")
        aln = os.path.join(ALN, f"{mk}_aln.fasta")
        sh(f'{env("mafft")} --auto --adjustdirection --anysymbol "{raw}" > "{aln}" 2>/dev/null')
        # limpiar prefijo _R_ que agrega adjustdirection
        seqs = list(SeqIO.parse(aln, "fasta"))
        for s in seqs: s.id = s.id.replace("_R_",""); s.description=""
        SeqIO.write(seqs, aln, "fasta")
        tr = os.path.join(TRIM, f"{mk}_trim.fasta")
        sh(f'{env("trimal")} -in "{aln}" -out "{tr}" -automated1')
        trecs = list(SeqIO.parse(tr, "fasta"))
        if trecs and len(trecs[0].seq) >= 80:
            kept[mk] = {s.id: s for s in trecs}
            print(f"[ok] {mk}: {len(trecs)} taxa, {len(trecs[0].seq)} bp (post-trim)")
        else:
            print(f"[skip] {mk}: alineamiento vacio tras trim")
    return kept

def concat(kept):
    taxa = set()
    for mk in kept: taxa |= set(kept[mk].keys())
    # filtrar taxa con pocos marcadores (excepto focal)
    counts = {t: sum(1 for mk in kept if t in kept[mk]) for t in taxa}
    taxa = sorted(t for t in taxa if counts[t] >= MIN_MARKERS_PER_TAXON
                  or t.startswith("Humicolopsis"))
    lengths = {mk: len(next(iter(kept[mk].values())).seq) for mk in kept}
    order = [mk for mk in MARKERS if mk in kept]
    seqs = {}
    for t in taxa:
        s = ""
        for mk in order:
            s += str(kept[mk][t].seq) if t in kept[mk] else "-"*lengths[mk]
        seqs[t] = s
    fa = os.path.join(CONC, "concatenated.fasta")
    with open(fa, "w") as fh:
        for t in taxa: fh.write(f">{t}\n{seqs[t]}\n")
    # particiones
    pos = 1; parts = []
    with open(os.path.join(CONC, "partition.txt"), "w") as fh:
        for mk in order:
            L = lengths[mk]; fh.write(f"DNA, {mk} = {pos}-{pos+L-1}\n"); parts.append((mk,pos,pos+L-1)); pos += L
    total = pos-1
    # reporte cobertura
    with open(os.path.join(CONC, "taxa_matrix.tsv"), "w") as fh:
        fh.write("taxon\t"+"\t".join(order)+"\tn\n")
        for t in taxa:
            row=["1" if t in kept[mk] else "0" for mk in order]
            fh.write(t+"\t"+"\t".join(row)+f"\t{sum(int(x) for x in row)}\n")
    print(f"\nMatriz: {len(taxa)} taxa x {total} bp, {len(order)} particiones -> {fa}")
    return fa, total, len(taxa)

def run_iqtree(fa):
    part = os.path.join(CONC, "partition.txt")
    pre = os.path.join(TREE, "Humicolopsis_v2")
    # NB: no se usa -o (crashea con MFP+MERGE en 2.4.0); se enraiza en Neurospora
    # en el paso de visualizacion (04_visualize.py).
    cmd = (f'{env("iqtree2")} -s "{fa}" -p "{part}" -m MFP+MERGE -B 1000 '
           f'-T AUTO --prefix "{pre}" -redo')
    print("IQ-TREE2:", cmd)
    r = sh(cmd)
    print(r.stdout[-1500:] if r.stdout else r.stderr[-1500:])

if __name__ == "__main__":
    data = gather()
    print("=== recuento crudo por marcador ===")
    for mk in MARKERS: print(f"  {mk}: {len(data[mk])} taxa")
    kept = align_trim(data)
    fa, total, ntax = concat(kept)
    run_iqtree(fa)
