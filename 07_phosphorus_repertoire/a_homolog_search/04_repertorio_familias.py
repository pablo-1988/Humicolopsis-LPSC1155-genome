"""A3 paso 4: repertorio de genes de P y ácidos orgánicos por especie.

Familia = unión de los ortogrupos (OrthoFinder) que contienen el gen de referencia caracterizado.
 - Ancla de cada referencia: su propio gen en el proteoma de origen si está en el set
   (S. cerevisiae, N. crassa, A. niger, S. sclerotiorum, B. cinerea; identidad >= 90%);
   si no, el mejor hit global.
 - Se agrega el ortogrupo del mejor hit por especie (cobertura >= 60%) para no perder ortólogos
   que OrthoFinder separa del ortogrupo de levadura (p. ej. PHO81, PHO87/90/91).
 - Además: conteo por dominio Pfam (hmmsearch --cut_ga) por especie.
Salidas (resultados/): familias_copias.tsv, familias_genes.tsv, pfam_copias.tsv
"""
import csv, re
from collections import defaultdict
from pathlib import Path
import os
SUF = os.environ.get("PROT_SET", "")   # "" = proteomas originales de A3; "_A7prot" = proteomas reanotados de A7 (carpetas nuevas)

W = Path.home() / "hongo_comparativa"
HERE = Path(__file__).parent
OF = sorted((W / f"A3{SUF}/prot/OrthoFinder").glob("Results_*"))[-1]
OUT = HERE / f"resultados{SUF}"; OUT.mkdir(exist_ok=True)
SRC = {"YEAST": "Sacccere", "NEUCR": "Neurcras", "ASPNG": "Aspenige", "ASPNA": "Aspenige",
       "SCLS1": "Sclescle", "BOTFB": "Botrcine", "BOTFU": "Botrcine"}

# ortogrupos
og_of, og_genes = {}, {}
rows = list(csv.reader(open(OF / "Orthogroups/Orthogroups.tsv"), delimiter="\t"))
species = rows[0][1:]
for r in rows[1:]:
    og_genes[r[0]] = defaultdict(list)
    for sp, cell in zip(species, r[1:]):
        for g in filter(None, (x.strip() for x in cell.split(","))):
            og_of[g] = r[0]; og_genes[r[0]][sp].append(g)
for r in csv.reader(open(OF / "Orthogroups/Orthogroups_UnassignedGenes.tsv"), delimiter="\t"):
    if r[0] == "Orthogroup": continue
    for sp, cell in zip(species, r[1:]):
        if cell.strip():
            og_of[cell.strip()] = r[0]; og_genes[r[0]] = defaultdict(list); og_genes[r[0]][sp].append(cell.strip())

# referencias
meta = {r["gen"]: r for r in csv.DictReader(open(HERE / "refs/referencias.tsv"), delimiter="\t")}
hits = defaultdict(list)
for l in open(W / f"A3{SUF}/refs_vs_all.tsv"):
    q, s, pid, ln, ql, sl, ev, bs = l.rstrip().split("\t")
    hits[q].append((s, float(pid), float(bs), int(ln) / int(ql)))
fam_ogs, anchor = defaultdict(set), {}
for q, hs in hits.items():
    _, gene, acc, _ = q.split("|")
    m = meta[gene]
    org = re.search(r"_([A-Z0-9]+) ", m["header"]).group(1) if m["header"] else ""
    src = SRC.get(org)
    cand = [h for h in hs if src and h[0].startswith(src + "|") and h[1] >= 90 and h[3] >= 0.8]
    best = max(cand or hs, key=lambda h: h[2])
    anchor[gene] = (best[0], best[1], "propio" if cand else "mejor_hit")
    fam_ogs[(m["categoria"], m["familia"])].add(og_of.get(best[0], "NA"))
    # + ortogrupo del mejor hit en cada especie (cobertura >= 60% de la referencia):
    # evita falsas ausencias cuando el ortólogo de Leotiomycetes cae en otro ortogrupo que el de levadura.
    per_sp = {}
    for h in hs:
        sp = h[0].split("|")[0]
        if h[3] >= 0.6 and (sp not in per_sp or h[2] > per_sp[sp][2]): per_sp[sp] = h
    for h in per_sp.values():
        fam_ogs[(m["categoria"], m["familia"])].add(og_of.get(h[0], "NA"))

with open(OUT / "familias_copias.tsv", "w", newline="") as fc, open(OUT / "familias_genes.tsv", "w", newline="") as fg:
    wc, wg = csv.writer(fc, delimiter="\t"), csv.writer(fg, delimiter="\t")
    wc.writerow(["categoria", "familia", "ortogrupos", "referencias"] + species)
    wg.writerow(["categoria", "familia", "especie", "gen", "ortogrupo"])
    for (cat, fam), ogs in sorted(fam_ogs.items()):
        refs = [g for g in meta if (meta[g]["categoria"], meta[g]["familia"]) == (cat, fam) and g in anchor]
        counts = {sp: 0 for sp in species}
        for og in ogs:
            if og == "NA": continue
            for sp, gs in og_genes[og].items():
                counts[sp] += len(gs)
                for g in gs: wg.writerow([cat, fam, sp, g.split("|", 1)[1], og])
        wc.writerow([cat, fam, ",".join(sorted(ogs)), ",".join(refs)] + [counts[sp] for sp in species])

# dominios Pfam por especie
pf = defaultdict(lambda: defaultdict(set))
for l in open(W / f"A3{SUF}/pfam_familias.domtbl"):
    if l.startswith("#"): continue
    f = l.split()
    pf[f[3]][f[0].split("|")[0]].add(f[0])
with open(OUT / "pfam_copias.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["pfam"] + species)
    for d in sorted(pf): w.writerow([d] + [len(pf[d][sp]) for sp in species])
print("OrthoFinder:", OF.name, "| familias:", len(fam_ogs))
for g, (s, pid, how) in sorted(anchor.items()): print(f"  {g:10s} -> {s:32s} {pid:5.1f}% {how}")
