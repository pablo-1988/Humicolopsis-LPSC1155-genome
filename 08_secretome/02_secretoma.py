"""A5 paso 2 — Secretoma predicho y comparación entre genomas.

Proteína secretada (definición clásica, sin herramientas con licencia):
  (1) péptido señal según DeepSig (modelo eucarionte), y
  (2) sin hélices ni barriles transmembrana según TMbed fuera del péptido señal, y
  (3) sin señal C-terminal de retención en el retículo (HDEL/KDEL).
No se predice anclaje GPI (no hay un predictor abierto y confiable; NetGPI requiere licencia): el secretoma puede incluir
proteínas de pared ancladas por GPI.
Salidas (resultados/):
  secretoma_por_genoma.tsv       proteoma, péptido señal, secretadas (%), pequeñas (< 300 aa), CAZymes secretadas
  secretadas_P_por_genoma.tsv    proteínas secretadas en cada familia de P (ortogrupos de A7 de los genes de LPSC)
  secretadas_LPSC.tsv            lista de proteínas secretadas de LPSC 1155 con anotación eggNOG
  expandidas_secretadas.txt      ¿las familias expandidas en Humicolopsis están enriquecidas en proteínas secretadas?
  control_concordancia.txt       concordancia DeepSig vs TMbed para el péptido señal
Uso: python 02_secretoma.py   (env `comparativa`: scipy)"""
import csv, re, glob, collections
from pathlib import Path
from scipy.stats import fisher_exact

HERE = Path(__file__).parent; OUT = HERE / "resultados"; OUT.mkdir(exist_ok=True)
W = Path.home() / "hongo_comparativa"; A = W / "A5"
PROJ = HERE.parents[1]
REFS = {"Sacccere", "Aspenige"}

# proteomas
seq = {}
for f in glob.glob(str(W / "A7/prot/*.faa")):
    n = None
    for l in open(f):
        if l.startswith(">"): n = l[1:].split()[0]; seq[n] = []
        else: seq[n].append(l.strip())
seq = {k: "".join(v) for k, v in seq.items()}
sp_of = lambda g: g.split("|")[0]
species = sorted({sp_of(g) for g in seq})

# DeepSig
sp_end = {}
for f in glob.glob(str(A / "deepsig/*.gff3")):
    for l in open(f):
        r = l.rstrip("\n").split("\t")
        if len(r) > 4 and r[2] == "Signal peptide": sp_end[r[0]] = int(r[4])

# TMbed (formato 0: id / secuencia / etiquetas; B,b barril; H,h hélice; S péptido señal; . resto)
tm = {}
for f in glob.glob(str(A / "tmbed/chunk[0-9][0-9].pred")):
    L = open(f).read().split("\n")
    for i in range(0, len(L) - 2, 3):
        if L[i].startswith(">"): tm[L[i][1:].split()[0]] = L[i + 2]
def n_tm(lab, sp_fin):  # segmentos transmembrana (corridas de H/h/B/b) que empiezan después del péptido señal de DeepSig
    # (si TMbed no reconoce el péptido señal, su núcleo hidrofóbico puede aparecer como hélice: se ignora)
    return sum(1 for m in re.finditer(r"[HhBb]+", lab) if m.start() >= sp_fin)

secr, motivos = set(), collections.Counter()
for g, e in sp_end.items():
    lab = tm.get(g)
    if lab is None: motivos["sin TMbed"] += 1; continue
    if n_tm(lab, e) > 0: motivos["con TM"] += 1; continue
    if re.search(r"[HK]DEL$", seq[g].rstrip("*")): motivos["retención RE"] += 1; continue
    secr.add(g)
with open(OUT / "control_concordancia.txt", "w") as fh:
    both = sum(1 for g in sp_end if g in tm and tm[g].startswith("S"))
    fh.write(f"proteínas con péptido señal (DeepSig): {len(sp_end)}\n")
    fh.write(f"con predicción TMbed: {sum(1 for g in sp_end if g in tm)}\n")
    fh.write(f"TMbed también predice péptido señal: {both} ({both / max(1, len(sp_end)):.1%})\n")
    fh.write(f"descartadas: {dict(motivos)}\nsecretadas: {len(secr)}\n")

# anotación eggNOG
ann = {}
for l in open(W / "A7/eggnog/todos.emapper.annotations"):
    if l.startswith("#"): continue
    r = l.rstrip("\n").split("\t"); ann[r[0]] = {"name": r[8], "desc": r[7], "cazy": r[18], "pfam": r[20]}

# por genoma
rows = []
for s in species:
    P = [g for g in seq if sp_of(g) == s]
    S = [g for g in P if g in sp_end]; X = [g for g in P if g in secr]
    small = [g for g in X if len(seq[g].rstrip("*")) < 300]
    caz = [g for g in X if ann.get(g, {}).get("cazy", "-") not in ("-", "")]
    rows.append([s, len(P), len(S), len(X), round(100 * len(X) / len(P), 2), len(small), len(caz)])
with open(OUT / "secretoma_por_genoma.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genoma", "proteinas", "peptido_senal", "secretadas", "pct_secretadas", "secretadas_<300aa", "CAZymes_secretadas"]); w.writerows(rows)

# familias de P: ortogrupos de A7 que contienen los genes de LPSC de cada familia (Table S5b)
og_of, og_mem = {}, collections.defaultdict(list)
for r in list(csv.reader(open(W / "A7/prot/OrthoFinder/Results_A7/Orthogroups/Orthogroups.tsv"), delimiter="\t"))[1:]:
    for c in r[1:]:
        for g in filter(None, (x.strip() for x in c.split(","))): og_of[g] = r[0]; og_mem[r[0]].append(g)
fam_genes = collections.defaultdict(set)
for r in csv.DictReader(open(PROJ / "04_suplementario/Table_S5b_P_repertoire.tsv"), delimiter="\t"):
    fam_genes[(r["Category"], r["Gene family"])].add("HumcLPSC|" + r["LPSC 1155 gene ID"])
prow = []
for (cat, fam), gs in fam_genes.items():
    ogs = {og_of[g] for g in gs if g in og_of}
    mem = [g for o in ogs for g in og_mem[o]]
    c = collections.Counter(sp_of(g) for g in mem if g in secr)
    if sum(c.values()) == 0: continue
    prow.append([cat, fam, ",".join(sorted(ogs))] + [c.get(s, 0) for s in species])
with open(OUT / "secretadas_P_por_genoma.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["categoria", "familia", "ortogrupos"] + species); w.writerows(prow)

# lista de LPSC
with open(OUT / "secretadas_LPSC.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["gen", "largo_aa", "fin_peptido_senal", "ortogrupo", "nombre", "descripcion", "CAZy", "Pfam"])
    for g in sorted(x for x in secr if sp_of(x) == "HumcLPSC"):
        a = ann.get(g, {})
        w.writerow([g.split("|")[1], len(seq[g].rstrip("*")), sp_end[g], og_of.get(g, ""), a.get("name", ""), a.get("desc", ""), a.get("cazy", ""), a.get("pfam", "")])

# familias expandidas en el tronco de Humicolopsis vs fondo (proteínas de LPSC + CH-1.1)
A7R = PROJ / "analisis/A7_expansion_enriquecimiento/resultados"
exp = set(open(A7R / "enriquecimiento/ogs_Humicolopsis_tronco_expandidas.txt").read().split())
hum = [g for g in seq if sp_of(g) in ("HumcLPSC", "HumcCH11")]
inE = [g for g in hum if og_of.get(g) in exp]; outE = [g for g in hum if og_of.get(g) not in exp]
a, b = sum(g in secr for g in inE), len(inE) - sum(g in secr for g in inE)
c, d = sum(g in secr for g in outE), len(outE) - sum(g in secr for g in outE)
odds, p = fisher_exact([[a, b], [c, d]])
fams_sec = sorted({og_of[g] for g in inE if g in secr})
with open(OUT / "expandidas_secretadas.txt", "w") as fh:
    fh.write(f"proteínas de Humicolopsis en familias expandidas: {len(inE)}; secretadas {a} ({a / max(1, len(inE)):.1%})\n")
    fh.write(f"resto del proteoma: {len(outE)}; secretadas {c} ({c / max(1, len(outE)):.1%})\n")
    fh.write(f"Fisher: OR = {odds:.2f}, p = {p:.3g}\n")
    fh.write(f"familias expandidas con al menos una proteína secretada: {len(fams_sec)} de {len(exp)}\n")
    for o in fams_sec:
        ex = next(g for g in og_mem[o] if g in secr and sp_of(g) in ("HumcLPSC", "HumcCH11"))
        fh.write(f"  {o}\t{ann.get(ex, {}).get('name', '')}\t{ann.get(ex, {}).get('desc', '')[:80]}\n")
print(open(OUT / "control_concordancia.txt").read())
for r in rows: print("\t".join(map(str, r)))
print(open(OUT / "expandidas_secretadas.txt").read().split("familias expandidas con")[0])
