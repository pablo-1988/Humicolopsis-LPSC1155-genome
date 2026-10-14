"""A4 paso 1: árboles génicos por familia para curar ortología y número de copias.

Para cada familia clave:
 1. Homólogos = proteínas de los 24 proteomas con hit a las referencias de la familia
    (diamond, e <= 1e-10, cobertura de la referencia >= 50%) y/o con el dominio Pfam diagnóstico.
 2. + referencias caracterizadas (UniProt) como marcadores de función.
 3. MAFFT (--auto) -> trimAl (-gappyout) -> IQ-TREE 3 (MFP restringido a LG/WAG/JTT/Q.pfam/Q.yeast, UFBoot 1000).
Salida: arboles/<familia>.{faa,aln,treefile,contree}
La asignación de clados ortólogos se hace en el paso 2 (inspección + reglas por familia).
Uso: python 01_arboles_familia.py [familia ...]   (sin argumentos: todas)
"""
import csv, subprocess, sys
from collections import defaultdict
from pathlib import Path
import os
SUF = os.environ.get("PROT_SET", "")   # "" = proteomas originales de A3; "_A7prot" = proteomas reanotados de A7 (carpetas nuevas)

W = Path.home() / "hongo_comparativa"
HERE = Path(__file__).parent
A3 = HERE.parent / "A3_repertorio_fosforo_comparativo"
A = os.environ.get("A4_BIN", "/home/ubuntu/miniforge3/envs/comparativa/bin")   # Mac: .../comparativa_arm/bin
OUT = HERE / f"arboles{SUF}"; OUT.mkdir(exist_ok=True)
TMP = W / f"A4{SUF}"; TMP.mkdir(exist_ok=True)   # trabajo pesado fuera de Drive
EXCL = set()   # especies a excluir de los árboles (ninguna)
import os
THREADS = os.environ.get("A4_THREADS", "4")

# familia -> (referencias, dominios Pfam diagnósticos)
FAM = {
    "PHO84":   (["PHO84", "pho-5"], []),
    "PHO89":   (["PHO89", "pho-4"], ["PHO4"]),
    "PHO87_90_91": (["PHO87", "PHO90", "PHO91"], ["CitMHS"]),
    "SPX":     (["PHO81", "nuc-2", "VTC2", "VTC3", "VTC4", "PHO87", "PHO90", "PHO91", "SYG1", "GDE1"], ["SPX"]),
    "PHO81":   (["PHO81", "nuc-2"], []),
    "VTC":     (["VTC1", "VTC2", "VTC3", "VTC4"], ["VTC", "DUF202"]),
    "PHO4_NUC1": (["PHO4", "nuc-1"], []),
    "PHO80":   (["PHO80", "preg"], []),
    "PHO85":   (["PHO85", "mdk-1"], []),
    "HAP_fosfatasas_fitasas": (["PHO5", "PHO3", "PHO11", "phyA", "phyB"], ["His_Phos_2"]),
    "Fosfatasa_alcalina": (["PHO8", "pho-2"], ["Alk_phosphatase"]),
    "PhoD":    ([], ["PhoD"]),
    "PAP":     (["pho-3"], ["Pur_ac_phosph_N"]),
    "Fitasa_BPP_PTP": ([], ["Phytase", "PTPlike_phytase"]),
    "OAH_ICL_MCL": (["oahA", "Ss-oah1", "Bcoah", "ICL1", "ICL2"], ["PEP_mutase", "ICL"]),
    "Oxalato_descarboxilasa": (["oxdC"], []),
    "CexA":    (["cexA"], []),
    "Mae1_TDT": (["mae1"], ["SLAC1"]),
    "Citrato_sintasa": (["CIT1", "CIT2"], ["Citrate_synt"]),
    "Malato_DH": (["MDH1", "MDH2", "MDH3"], ["Ldh_1_N"]),
}

# proteomas y referencias
seqs, n = {}, None
for f in sorted((W / f"A3{SUF}/prot").glob("*.faa")):
    for l in open(f):
        if l.startswith(">"): n = l[1:].split()[0]; seqs[n] = []
        else: seqs[n].append(l.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
refseq, n = {}, None
for l in open(A3 / "refs/referencias.faa"):
    if l.startswith(">"): n = l[1:].strip(); refseq[n] = ""
    else: refseq[n] += l.strip()
ref_by_gene = {k.split("|")[1]: k for k in refseq}

hits = defaultdict(set)
for l in open(W / f"A3{SUF}/refs_vs_all.tsv"):
    q, s, pid, ln, ql, sl, ev, bs = l.rstrip().split("\t")
    if float(ev) <= 1e-10 and int(ln) / int(ql) >= 0.5:
        hits[q.split("|")[1]].add(s)
dom = defaultdict(set)
for l in open(W / f"A3{SUF}/pfam_familias.domtbl"):
    if l.startswith("#"): continue
    f = l.split(); dom[f[3]].add(f[0])

todo = sys.argv[1:] or list(FAM)
for fam in todo:
    refs, doms = FAM[fam]
    members = set().union(*[hits[r] for r in refs]) if refs else set()
    for d in doms: members |= dom[d]
    members = {m for m in members if m.split("|")[0] not in EXCL and len(seqs[m]) >= 80}
    faa = TMP / f"{fam}.faa"
    with open(faa, "w") as fh:
        for m in sorted(members): fh.write(f">{m}\n{seqs[m]}\n")
        for r in refs:
            k = ref_by_gene.get(r)
            if k: fh.write(f">{k.replace('|', '_')}\n{refseq[k]}\n")
    aln, trm = TMP / f"{fam}.aln", TMP / f"{fam}.trim.aln"
    with open(aln, "w") as fh:
        subprocess.run([f"{A}/mafft", "--auto", "--thread", "2", faa], stdout=fh, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([f"{A}/trimal", "-in", aln, "-out", trm, "-gappyout"], check=True)
    subprocess.run([f"{A}/iqtree3", "-s", trm, "-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam,Q.yeast", "-B", "1000", "-T", THREADS, "--prefix", TMP / fam,
                    "-redo", "--quiet"], check=True)
    for ext in (".faa", ".trim.aln", ".treefile", ".contree", ".iqtree"):
        src = TMP / f"{fam}{ext}"
        if src.exists(): (OUT / f"{fam}{ext}").write_bytes(src.read_bytes())
    print(f"{fam}: {len(members)} homólogos + {len(refs)} refs")
