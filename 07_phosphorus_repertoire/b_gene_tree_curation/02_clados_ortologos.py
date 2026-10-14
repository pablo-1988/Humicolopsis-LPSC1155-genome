"""A4 paso 2: asignación de clados ortólogos anclados en referencias caracterizadas.

Para cada subfamilia (grupo de referencias con la misma función):
  - si el árbol tiene referencias de otras subfamilias: clado = nodo más alto (desde la referencia hacia la raíz,
    enraizado en punto medio) que NO contiene referencias ajenas;
  - si es la única subfamilia y tiene >= 2 referencias: clado = MRCA de sus referencias;
  - si es la única subfamilia con 1 referencia: todo el árbol (familia definida por dominio/homología).
Se reportan: soporte UFBoot del clado y número de copias por especie.
Salida: resultados/copias_curadas.tsv, resultados/genes_curados.tsv
Ejecutar con el python del env `filogenia` (Bio.Phylo).
"""
import csv
from collections import Counter
from pathlib import Path
import os
SUF = os.environ.get("PROT_SET", "")   # "" = proteomas originales de A3; "_A7prot" = proteomas reanotados de A7 (carpetas nuevas)
from Bio import Phylo

HERE = Path(__file__).parent
ARB = HERE / f"arboles{SUF}"; OUT = HERE / f"resultados{SUF}"; OUT.mkdir(exist_ok=True)

# árbol -> {subfamilia: [genes de referencia]}
SUB = {
    "PHO84": {"PHO84": ["PHO84", "pho-5"]},
    "PHO89": {"PHO89": ["PHO89", "pho-4"]},
    "PHO87_90_91": {"PHO87/90/91": ["PHO87", "PHO90", "PHO91"]},
    "SPX": {"PHO81": ["PHO81", "nuc-2"], "VTC2/3": ["VTC2", "VTC3"], "VTC4": ["VTC4"],
            "PHO87/90/91": ["PHO87", "PHO90", "PHO91"], "SYG1": ["SYG1"], "GDE1": ["GDE1"]},
    "PHO81": {"PHO81": ["PHO81", "nuc-2"]},
    "VTC": {"VTC1": ["VTC1"], "VTC2/3": ["VTC2", "VTC3"], "VTC4": ["VTC4"]},
    "PHO4_NUC1": {"PHO4/NUC-1": ["PHO4", "nuc-1"]},
    "PHO80": {"PHO80": ["PHO80", "preg"]},
    "PHO85": {"PHO85": ["PHO85", "mdk-1"]},
    "HAP_fosfatasas_fitasas": {"PHO5-like (ác. fosfatasa)": ["PHO5", "PHO3", "PHO11"], "Fitasa PhyA": ["phyA"], "Fitasa PhyB": ["phyB"]},
    "Fosfatasa_alcalina": {"PHO8": ["PHO8", "pho-2"]},
    "PAP": {"PAP": ["pho-3"]},
    "OAH_ICL_MCL": {"OAH": ["oahA", "Ss-oah1", "Bcoah"], "ICL": ["ICL1"], "MCL": ["ICL2"]},
    "Oxalato_descarboxilasa": {"Oxalato descarboxilasa": ["oxdC"]},
    "CexA": {"CexA": ["cexA"]},
    "Mae1_TDT": {"Mae1": ["mae1"]},
    "Citrato_sintasa": {"CIT mitocondrial": ["CIT1"], "CIT peroxisomal": ["CIT2"]},
    "Malato_DH": {"MDH1 (mit.)": ["MDH1"], "MDH2 (cit.)": ["MDH2"], "MDH3 (perox.)": ["MDH3"]},
}
SPECIES = ["HumcLPSC", "HumcCH11", "Allacyti", "Allalyco", "Xenopyri", "Leotlubr", "Clausp", "Pseudest", "Cyclminu",
           "Chaeraph", "Oidimaiu", "Hyalhepa", "Phiasuba", "Mollscop", "Cadosp", "Glarlozo", "Hymefrax", "Erysneca",
           "Sclescle", "Botrcine", "Chloaeru", "Neurcras", "Aspenige", "Sacccere"]

def refname(leaf):   # REF_<gen>_<acc>_<fam>
    return leaf.name.split("_")[1] if leaf.name.startswith("REF_") else None

rows, genes = [], []
for fam, subs in SUB.items():
    f = ARB / f"{fam}.contree"
    if not f.exists(): continue
    t = Phylo.read(f, "newick"); t.root_at_midpoint()
    reflea = {refname(l): l for l in t.get_terminals() if refname(l)}
    for sub, rg in subs.items():
        mine = [reflea[r] for r in rg if r in reflea]
        if not mine: continue
        other = {l for r, l in reflea.items() if r not in rg}
        node = t.common_ancestor(mine) if len(mine) > 1 else mine[0]
        if other & set(node.get_terminals()):
            rows.append([fam, sub, "CONFLICTO: referencias de otra subfamilia dentro del MRCA", ""] + [""] * len(SPECIES)); continue
        best = node
        if other:   # con otras subfamilias en el árbol: subir hasta toparse con una referencia ajena
            path = t.get_path(node)
            for anc in reversed([t.root] + path[:-1]):
                if other & set(anc.get_terminals()): break
                best = anc
        elif len(mine) == 1:
            best = t.root   # familia de una sola referencia, definida por dominio/homología: todo el árbol
        leaves = [l for l in best.get_terminals() if "|" in l.name]
        cnt = Counter(l.name.split("|")[0] for l in leaves)
        sup = best.confidence if best.confidence is not None else ""
        rows.append([fam, sub, len(leaves), sup] + [cnt.get(s, 0) for s in SPECIES])
        for l in leaves:
            sp, g = l.name.split("|", 1); genes.append([fam, sub, sp, g])

with open(OUT / "copias_curadas.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["arbol", "subfamilia", "n", "UFBoot"] + SPECIES); w.writerows(rows)
with open(OUT / "genes_curados.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["arbol", "subfamilia", "especie", "gen"]); w.writerows(genes)
print(f"{'subfamilia':28s} {'n':>4s} {'UFB':>4s}" + "".join(f"{s[:5]:>6s}" for s in SPECIES))
for r in rows:
    print(f"{(r[0][:10]+':'+str(r[1]))[:28]:28s} {str(r[2])[:4]:>4s} {str(r[3]):>4s}" + "".join(f"{str(x):>6s}" for x in r[4:]))
