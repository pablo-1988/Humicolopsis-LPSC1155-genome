"""A1 paso 2: matriz de ortólogos BUSCO de copia única -> árboles.

1. Proteínas de loci 'Single' de compleasm (leotiomycetes_odb10) por taxón.
2. Filtro de ocupación (>= OCC de los taxa) y presencia obligatoria de LPSC 1155.
3. MAFFT L-INS-i + trimAl -automated1 por locus (paralelo).
4. (b) Árboles génicos IQ-TREE3 (MFP) -> ASTRAL-IV (coalescente)
   (a) Concatenado IQ-TREE3 (modelo por BIC en 300 loci al azar, UFBoot 1000, -mem 20G)
   (c) Factores de concordancia gCF/sCF sobre el árbol concatenado.
Excluye referencias funcionales (S. cerevisiae, A. niger). Grupo externo: N. crassa.
Incluye taxones parciales sin genoma (partial_taxa.tsv): CBS 645.76 ex-tipo con MCM7, RPB1 y EF1a.
La ocupación se calcula solo sobre los taxones con genoma.
"""
import csv, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

W = Path.home() / "hongo_comparativa"
A = "/Users/pablo/.local/share/mamba/envs/comparativa_arm/bin"
OUT = W / "A1_phylogenomics"
OCC = float(sys.argv[1]) if len(sys.argv) > 1 else 0.9
EXCL = {"Sacccere", "Aspenige"}
FOCAL, OUTG = "HumcLPSC", "Neurcras"
CPU = 11

taxa = [r["code"] for r in csv.DictReader(open(W / "taxa_map.tsv"), delimiter="\t")
        if r["code"] not in EXCL and (W / "compleasm" / r["code"]).exists()]
# Taxones parciales sin genoma (p. ej. ex-tipo CBS 645.76: marcadores Sanger asignados a su BUSCO)
PARTIAL = {}
if (W / "partial_taxa.tsv").exists():
    for r in csv.DictReader(open(W / "partial_taxa.tsv"), delimiter="\t"):
        PARTIAL.setdefault(r["taxon"], {})[r["busco"]] = r["seq"]
for d in ("loci_raw", "loci_aln", "loci_trim", "genetrees"):
    (OUT / d).mkdir(parents=True, exist_ok=True)

# 1. recolectar secuencias
loci = {}
for t in taxa:
    lin = W / "compleasm" / t / "leotiomycetes_odb10"
    single = {r["Gene"] for r in csv.DictReader(open(lin / "full_table.tsv"), delimiter="\t") if r["Status"] == "Single"}
    name, seq = None, []
    for line in list(open(lin / "gene_marker.fasta")) + [">END"]:
        if line.startswith(">"):
            if name in single:
                loci.setdefault(name, {})[t] = "".join(seq).replace("*", "")
            name, seq = line[1:].strip().rsplit("_", 1)[0], []
        else:
            seq.append(line.strip())

# 2. ocupación
keep = sorted(b for b, d in loci.items() if len(d) >= OCC * len(taxa) and FOCAL in d and OUTG in d)
for pt, d in PARTIAL.items():
    for b, sq in d.items():
        if b in loci: loci[b][pt] = sq
genome_taxa = list(taxa)
taxa = taxa + list(PARTIAL)
print(f"{len(loci)} loci totales; {len(keep)} con ocupación >= {OCC:.0%} (+ LPSC 1155 y grupo externo)")
for b in keep:
    (OUT / "loci_raw" / f"{b}.faa").write_text("".join(f">{t}\n{s}\n" for t, s in loci[b].items()))

# 3. alinear + recortar
def aln(b):
    a, t = OUT / "loci_aln" / f"{b}.faa", OUT / "loci_trim" / f"{b}.faa"
    if not t.exists():
        with a.open("w") as fh:
            subprocess.run([f"{A}/mafft", "--localpair", "--maxiterate", "1000", "--thread", "1",
                            OUT / "loci_raw" / f"{b}.faa"], stdout=fh, stderr=subprocess.DEVNULL, check=True)
        subprocess.run([f"{A}/trimal", "-in", a, "-out", t, "-automated1"], check=True)
with ThreadPoolExecutor(CPU) as ex:
    list(ex.map(aln, keep))

# matriz concatenada + particiones
def read(p):
    d, n = {}, None
    for l in open(p):
        if l.startswith(">"): n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}
cat, parts, pos = {t: [] for t in taxa}, [], 1
for b in keep:
    s = read(OUT / "loci_trim" / f"{b}.faa")
    L = len(next(iter(s.values())))
    if L < 50: continue
    for t in taxa: cat[t].append(s.get(t, "-" * L))
    parts.append(f"AA, {b} = {pos}-{pos + L - 1}"); pos += L
(OUT / "concat.faa").write_text("".join(f">{t}\n{''.join(v)}\n" for t, v in cat.items()))
(OUT / "partitions.txt").write_text("\n".join(parts) + "\n")
print(f"Matriz: {len(taxa)} taxa x {pos - 1} aa, {len(parts)} loci")

# 4b. árboles génicos + ASTRAL
def gtree(b):
    p = OUT / "genetrees" / b
    if not Path(f"{p}.treefile").exists():
        subprocess.run([f"{A}/iqtree3", "-s", OUT / "loci_trim" / f"{b}.faa", "-m", "MFP", "-T", "1",
                        "--prefix", p, "--quiet"], check=False)
with ThreadPoolExecutor(CPU) as ex:
    list(ex.map(gtree, [p.split(" = ")[0].split(", ")[1] for p in parts]))
with (OUT / "genetrees.nwk").open("w") as fh:
    for f in sorted((OUT / "genetrees").glob("*.treefile")): fh.write(f.read_text())
subprocess.run([f"{A}/astral4", "-t", str(CPU), "-i", OUT / "genetrees.nwk", "-o", OUT / "astral.nwk"],
               check=True, stderr=subprocess.DEVNULL)
# 4a. concatenado. La matriz completa (~1.5 M aa) no permite ModelFinder en 26 GB de RAM:
# se elige el modelo con 300 loci al azar y se aplica a la matriz completa con -mem.
import random
random.seed(1)
sub = random.sample([p.split(" = ")[0].split(", ")[1] for p in parts], 300)
subcat = {t: [] for t in taxa}
for b in sub:
    s = read(OUT / "loci_trim" / f"{b}.faa"); L = len(next(iter(s.values())))
    for t in taxa: subcat[t].append(s.get(t, "-" * L))
# taxones sin datos en la submuestra (p. ej. CBS 645.76 parcial) se omiten: IQ-TREE no acepta secuencias solo-gaps
(OUT / "concat_sub300.faa").write_text("".join(f">{t}\n{''.join(v)}\n" for t, v in subcat.items()
                                               if set("".join(v)) - {"-"}))
subprocess.run([f"{A}/iqtree3", "-s", OUT / "concat_sub300.faa", "-m", "MF", "-mset", "LG,WAG,JTT,Q.pfam,Q.yeast",
                "-T", str(CPU), "--prefix", OUT / "modelsel_sub300", "-redo", "--quiet"], check=True)
best = [l.split(":", 1)[1].strip() for l in open(OUT / "modelsel_sub300.iqtree") if l.startswith("Best-fit model according to BIC")][0]
print("Modelo (BIC, 300 loci):", best)
if not (OUT / "concat.treefile").exists():
    subprocess.run([f"{A}/iqtree3", "-s", OUT / "concat.faa", "-m", best, "-B", "1000", "-T", str(CPU), "-mem", "20G",
                    "--prefix", OUT / "concat", "-o", OUTG, "-redo", "--quiet"], check=True)
# 4c. concordancia
subprocess.run([f"{A}/iqtree3", "-t", OUT / "concat.treefile", "--gcf", OUT / "genetrees.nwk",
                "-s", OUT / "concat.faa", "--scf", "100", "-T", str(CPU), "--prefix", OUT / "concord", "-redo", "--quiet"],
               check=True)
print("listo:", OUT)
