"""A4 paso 3: tabla final del repertorio de P y ácidos orgánicos (copias por especie).

Fuente por familia:
 - "árbol": clado ortólogo curado en árbol génico (02_clados_ortologos.py -> copias_curadas.tsv)
 - "OG": ortogrupo(s) de OrthoFinder anclados en referencia (A3 04_repertorio_familias.py) — solo para
   genes de metabolismo central de copia baja donde las referencias de levadura son parálogos propios
   de levadura (citrato sintasa, malato DH, etc.) o cuando el árbol no resolvió el anclaje.
Salida: resultados/tabla_final_repertorio.tsv
"""
import csv
from pathlib import Path
import os
SUF = os.environ.get("PROT_SET", "")   # "" = proteomas originales de A3; "_A7prot" = proteomas reanotados de A7 (carpetas nuevas)
HERE = Path(__file__).parent
A3 = HERE.parent / "A3_repertorio_fosforo_comparativo" / f"resultados{SUF}" / "familias_copias.tsv"
CUR = HERE / f"resultados{SUF}" / "copias_curadas.tsv"

# (categoría, etiqueta en inglés, fuente, clave)
ROWS = [
    ("P sensing & regulation", "PHO85 (CDK)", "árbol", ("PHO85", "PHO85")),
    ("P sensing & regulation", "PHO80 (cyclin)", "OG", "PHO80"),
    ("P sensing & regulation", "PHO81 (CDK inhibitor, SPX)", "árbol", ("SPX", "PHO81")),
    ("P sensing & regulation", "PHO4/NUC-1 (bHLH TF)", "árbol", ("PHO4_NUC1", "PHO4/NUC-1")),
    ("P sensing & regulation", "KCS1 (IP6K)", "OG", "KCS1"),
    ("P sensing & regulation", "VIP1 (PP-IP5K)", "OG", "VIP1"),
    ("Pi transport", "PHO84 (H+/Pi, high affinity)", "árbol", ("PHO84", "PHO84")),
    ("Pi transport", "PHO89 (Na+/Pi)", "árbol", ("PHO89", "PHO89")),
    ("Pi transport", "PHO87/90/91 (low affinity, SPX)", "árbol", ("SPX", "PHO87/90/91")),
    ("Pi transport", "SYG1 (Pi exporter, SPX)", "árbol", ("SPX", "SYG1")),
    ("PolyP & vacuolar homeostasis", "VTC1", "árbol", ("VTC", "VTC1")),
    ("PolyP & vacuolar homeostasis", "VTC2/3 (SPX)", "árbol", ("SPX", "VTC2/3")),
    ("PolyP & vacuolar homeostasis", "VTC4 (SPX)", "árbol", ("SPX", "VTC4")),
    ("PolyP & vacuolar homeostasis", "PPX1 (exopolyphosphatase)", "OG", "PPX1"),
    ("PolyP & vacuolar homeostasis", "PPN1 (endopolyphosphatase)", "OG", "PPN1"),
    ("Phosphatases & phytases", "HAP phytase (PhyA-like)", "árbol", ("HAP_fosfatasas_fitasas", "Fitasa PhyA")),
    ("Phosphatases & phytases", "HAP phytase (PhyB-like)", "árbol", ("HAP_fosfatasas_fitasas", "Fitasa PhyB")),
    ("Phosphatases & phytases", "Alkaline phosphatase (PHO8)", "árbol", ("Fosfatasa_alcalina", "PHO8")),
    ("Phosphatases & phytases", "Purple acid phosphatase", "árbol", ("PAP", "PAP")),
    ("Phosphatases & phytases", "Glycerophosphodiesterase (GDE1)", "árbol", ("SPX", "GDE1")),
    ("Organic acids", "Oxaloacetate acetylhydrolase (OAH)", "árbol", ("OAH_ICL_MCL", "OAH")),
    ("Organic acids", "Oxalate decarboxylase", "árbol", ("Oxalato_descarboxilasa", "Oxalato descarboxilasa")),
    ("Organic acids", "Citrate synthase", "OG", "Citrato sintasa"),
    ("Organic acids", "ATP-citrate lyase", "OG", "ATP-citrato liasa"),
    ("Organic acids", "Malate dehydrogenase", "OG", "Malato DH"),
    ("Organic acids", "Isocitrate lyase", "árbol", ("OAH_ICL_MCL", "ICL")),
    ("Organic acid transport", "Citrate exporter (CexA-like)", "árbol", ("CexA", "CexA")),
    ("Organic acid transport", "Malate transporter (Mae1/TDT)", "árbol", ("Mae1_TDT", "Mae1")),
    ("Organic acid transport", "Mitochondrial citrate carrier (CTP1)", "OG", "Transp. citrato mitocondrial"),
]
# Correcciones manuales tras inspección del árbol (ver RESULTADOS_A4.md):
# PHO89: el MRCA levadura+N. crassa tiene soporte bajo (UFBoot 63) y excluye al ortólogo verdadero de Humicolopsis;
# se usa el clado con UFBoot 82 que contiene ambas referencias.
# CexA: árbol de 483 homólogos (superfamilia MFS); clado ortólogo CexA = nodo de 28 secuencias con UFBoot 100 alrededor de A. niger cexA.
# En CH-1.1 las 2 copias mapean al mismo locus de LPSC (FUN_005710, >99%): posible duplicación reciente o del ensamble de CH-1.1.
MANUAL_CEXA = {'HumcLPSC': 1, 'HumcCH11': 2, 'Allacyti': 1, 'Allalyco': 1, 'Xenopyri': 1, 'Leotlubr': 0, 'Clausp': 1, 'Pseudest': 2, 'Cyclminu': 1, 'Chaeraph': 2, 'Oidimaiu': 2, 'Hyalhepa': 2, 'Phiasuba': 2, 'Mollscop': 1, 'Cadosp': 1, 'Glarlozo': 1, 'Hymefrax': 2, 'Erysneca': 0, 'Sclescle': 1, 'Botrcine': 1, 'Chloaeru': 0, 'Neurcras': 1, 'Aspenige': 1, 'Sacccere': 0}
MANUAL = {"Citrate exporter (CexA-like)": MANUAL_CEXA, "PHO89 (Na+/Pi)": {"HumcLPSC": 2, "HumcCH11": 1, "Allacyti": 0, "Allalyco": 1, "Xenopyri": 1}}
# Las correcciones manuales se hicieron sobre los árboles de los proteomas originales: con otro set de proteomas NO se
# aplican; esas familias quedan marcadas "(revisar)" para inspeccionar sus árboles nuevos.
REVISAR = set(MANUAL) if SUF else set()
if SUF: MANUAL = {}
cur = {(r["arbol"], r["subfamilia"]): r for r in csv.DictReader(open(CUR), delimiter="\t")}
og = {r["familia"]: r for r in csv.DictReader(open(A3), delimiter="\t")}
SPECIES = [c for c in next(csv.reader(open(CUR), delimiter="\t"))[4:]]
out = []
for cat, lab, src, key in ROWS:
    r = cur.get(key) if src == "árbol" else og.get(key)
    if lab in MANUAL and (r is None or not str(r.get("n", "")).isdigit()):
        r = {sp: 0 for sp in SPECIES}
    if r is None or (src == "árbol" and not str(r.get("n", "")).isdigit()):
        out.append([cat, lab, src + " (pendiente)"] + [""] * len(SPECIES)); continue
    vals = [r[s] for s in SPECIES]
    if lab in MANUAL:
        for sp, v in MANUAL[lab].items(): vals[SPECIES.index(sp)] = str(v)
        src = src + "+manual"
    out.append([cat, lab, src + (" (revisar)" if lab in REVISAR else "")] + vals)
with open(HERE / f"resultados{SUF}" / "tabla_final_repertorio.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["category", "family", "source"] + SPECIES); w.writerows(out)
for r in out: print(f"{r[1][:38]:38s} {r[2]:8s} " + " ".join(f"{x:>2s}" for x in r[3:]))
