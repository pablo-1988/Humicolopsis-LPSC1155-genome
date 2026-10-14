"""A4 paso 4 (A7.6b): conteos del repertorio de P con los proteomas reanotados de A7, ANCLADOS en los clados ya curados.

Problema: al rehacer los árboles (y OrthoFinder) con los proteomas nuevos, la regla automática de 02_clados_ortologos.py
puede recuperar otro clado aunque los genes no cambien (p. ej. PHO84: UFBoot 100 / 59 secuencias en el árbol original,
UFBoot 82 / 39 en el nuevo), y los ortogrupos cambian con la versión de OrthoFinder. Para que la diferencia refleje SOLO
la reanotación de 7 genomas, cada familia se ancla en su definición curada:
  - árbol: clado = MRCA (árbol nuevo, raíz en punto medio) de los genes curados antes en las 17 especies cuyos proteomas
    no cambiaron; se cuentan ahí los genes de los 7 genomas reanotados. CexA y PHO89 usan los nodos curados a mano
    (CexA: nodo de 28 secuencias, UFBoot 100; PHO89: nodo de 20 secuencias con ambas referencias, UFBoot 82).
  - ortogrupo: familias = ortogrupos de A7 que contienen los genes de la familia en las especies sin cambios.
Las 17 especies sin cambios conservan sus conteos curados (mismas proteínas). Control: genes de especies sin cambios
que caen en el clado nuevo sin estar en el curado ("extra") -> si son muchos, el clado nuevo es más amplio.
Salida: resultados_A7prot/tabla_final_repertorio_anclada.tsv y resultados_A7prot/control_anclaje.tsv
Uso: python 04_reanclar_A7prot.py   (env `comparativa`, Bio.Phylo)"""
import csv, collections, importlib.util
from pathlib import Path
from Bio import Phylo

HERE = Path(__file__).parent
W = Path.home() / "hongo_comparativa"
A3R = HERE.parent / "A3_repertorio_fosforo_comparativo"
RE = {"HumcCH11", "Allacyti", "Allalyco", "Xenopyri", "Leotlubr", "Cyclminu", "Chaeraph"}   # reanotados en A7.1

# filas y fuentes de la tabla final (mismas que 03_tabla_final.py)
src = open(HERE / "03_tabla_final.py").read()
ROWS = eval(src.split("ROWS = ", 1)[1].split("\n]\n", 1)[0] + "\n]")
old = {r["family"]: r for r in csv.DictReader(open(HERE / "resultados/tabla_final_repertorio.tsv"), delimiter="\t")}
SPECIES = next(csv.reader(open(HERE / "resultados/tabla_final_repertorio.tsv"), delimiter="\t"))[3:]

# genes curados antes (árbol) y genes de familias por ortogrupo (A3)
cur = collections.defaultdict(set)
for r in csv.DictReader(open(HERE / "resultados/genes_curados.tsv"), delimiter="\t"):
    cur[(r["arbol"], r["subfamilia"])].add(f'{r["especie"]}|{r["gen"]}')
ogfam = collections.defaultdict(set)
for r in csv.DictReader(open(A3R / "resultados/familias_genes.tsv"), delimiter="\t"):
    ogfam[r["familia"]].add(f'{r["especie"]}|{r["gen"]}')

def nodo_manual(fam, ref, n_leaves, need_refs):
    t = Phylo.read(HERE / f"arboles/{fam}.contree", "newick"); t.root_at_midpoint()
    r = next(l for l in t.get_terminals() if l.name.startswith("REF_") and f"_{ref}_" in l.name)
    for anc in reversed(t.get_path(r)[:-1]):
        refs = {l.name.split("_")[1] for l in anc.get_terminals() if l.name.startswith("REF_")}
        if anc.count_terminals() == n_leaves and need_refs <= refs:
            return {l.name for l in anc.get_terminals() if "|" in l.name}
    raise SystemExit(f"no encontré el nodo manual de {fam}")
cur[("CexA", "CexA")] = nodo_manual("CexA", "cexA", 28, {"cexA"})
cur[("PHO89", "PHO89")] = nodo_manual("PHO89", "PHO89", 20, {"PHO89", "pho-4"})

# ortogrupos de A7
og_of, og_mem = {}, collections.defaultdict(list)
OFD = W / "A7/prot/OrthoFinder/Results_A7/Orthogroups"
for f in ("Orthogroups.tsv", "Orthogroups_UnassignedGenes.tsv"):
    rows = list(csv.reader(open(OFD / f), delimiter="\t"))
    for r in rows[1:]:
        for c in r[1:]:
            for g in filter(None, (x.strip() for x in c.split(","))): og_of[g] = r[0]; og_mem[r[0]].append(g)

out, qc = [], []
for cat, lab, fuente, key in ROWS:
    vals = {s: old[lab][s] for s in SPECIES}
    if fuente == "árbol":
        fam, sub = key
        anchors = {g for g in cur[(fam, sub)] if g.split("|")[0] not in RE}
        t = Phylo.read(HERE / f"arboles_A7prot/{fam}.contree", "newick"); t.root_at_midpoint()
        leaves = {l.name: l for l in t.get_terminals()}
        present = [leaves[g] for g in anchors if g in leaves]
        node = t.common_ancestor(present) if len(present) > 1 else present[0]
        mem = [l.name for l in node.get_terminals() if "|" in l.name]
        extra = [g for g in mem if g.split("|")[0] not in RE and g not in cur[(fam, sub)]]
        sup = node.confidence if node.confidence is not None else ""
        qc.append([lab, "árbol", len(anchors), len(present), len(mem), sup, len(extra), " ".join(sorted(extra)[:8])])
    else:
        anchors = {g for g in ogfam[key] if g.split("|")[0] not in RE}
        ogs = {og_of[g] for g in anchors if g in og_of}
        mem = [g for o in ogs for g in og_mem[o]]
        extra = [g for g in mem if g.split("|")[0] not in RE and g not in ogfam[key]]
        qc.append([lab, "OG", len(anchors), len([g for g in anchors if g in og_of]), len(mem), ",".join(sorted(ogs)), len(extra), " ".join(sorted(extra)[:8])])
    cnt = collections.Counter(g.split("|")[0] for g in mem)
    for s in RE: vals[s] = str(cnt.get(s, 0))
    out.append([cat, lab, fuente + "+anclado"] + [vals[s] for s in SPECIES])

# Correcciones tras revisar los hits contra las referencias (A3_A7prot/refs_vs_all.tsv), 2026-10-06:
#  PHO89: Allalyco FUN_007583 y Xenopyri FUN_007841 tienen 56–58% de identidad con N. crassa pho-4 (cobertura >= 0,9), igual que
#         CH-1.1; quedan fuera del clado anclado solo porque su soporte es bajo (UFBoot 29). Allacyti: sin hit -> 0.
#  SYG1:  en Allacyti/Allalyco los modelos de Funannotate están fragmentados (cobertura 0,23 y 0,08 de SYG1) y Xenopyri FUN_000734
#         cubre 0,63 como en otras especies: el gen está presente (1/1/1, igual que con miniprot).
#  PPX1:  cada genoma reanotado tiene un solo PPX1 de largo completo (cobertura ~1,0); los ortogrupos de A7 unen PPX1 con otra
#         familia (21 genes "extra" de especies sin cambios), lo que inflaba el conteo -> 1.
FIX = {"PHO89 (Na+/Pi)": {"Allalyco": "1", "Xenopyri": "1"},
       "SYG1 (Pi exporter, SPX)": {"Allacyti": "1", "Allalyco": "1", "Xenopyri": "1"},
       "PPX1 (exopolyphosphatase)": {s: "1" for s in RE}}
for r in out:
    for s, v in FIX.get(r[1], {}).items():
        r[3 + SPECIES.index(s)] = v; r[2] = r[2].replace("+manual", "") + "+manual"

R = HERE / "resultados_A7prot"
with open(R / "tabla_final_repertorio_anclada.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["category", "family", "source"] + SPECIES); w.writerows(out)
with open(R / "control_anclaje.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["family", "fuente", "anclas", "anclas_presentes", "miembros_clado", "UFBoot_u_OGs", "extra_especies_sin_cambio", "ejemplos_extra"]); w.writerows(qc)
print("cambios en los 7 genomas reanotados (antes -> ahora):")
for r in out:
    d = [f"{s} {old[r[1]][s]}->{r[3 + SPECIES.index(s)]}" for s in SPECIES if s in RE and old[r[1]][s] != r[3 + SPECIES.index(s)]]
    if d: print(f"  {r[1][:38]:38s} " + "; ".join(d))
print("\ncontrol (familias con genes 'extra' de especies sin cambios en el clado nuevo):")
for q in qc:
    if q[6]: print(f"  {q[0][:38]:38s} {q[1]:5s} anclas={q[2]} miembros={q[4]} soporte/OG={str(q[5])[:30]} extra={q[6]}")
