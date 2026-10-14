"""A7.5 — Familias expandidas/contraídas por rama (CAFE5) y enriquecimiento funcional (eggNOG-mapper).

Ramas evaluadas (rama que lleva al MRCA de los taxa indicados):
  Humicolopsis_tronco  = HumcLPSC + HumcCH11   (contraste principal: el linaje)
  LPSC_terminal        = HumcLPSC              (detalle propio de la cepa)
  Phacidiales_tronco   = Allacyti + Allalyco + Xenopyri (contraste)
Familia expandida (contraída) en una rama = familia significativa en CAFE5 (p < 0,05) + p de la rama < 0,05
+ cambio > 0 (< 0) en esa rama.

Anotación de cada ortogrupo: un término (GO, ruta KEGG, KO, Pfam) se asigna si lo tiene >= 50% de los genes
del ortogrupo que tienen alguna anotación de ese tipo (eggNOG-mapper sobre todos los proteomas).
Prueba: Fisher exacta unilateral (sobrerrepresentación) + Benjamini-Hochberg; fondo = todos los ortogrupos
analizados por CAFE5 con anotación de ese tipo. Solo biblioteca estándar de Python.

Uso: python3 07_enriquecimiento.py <dir_cafe> <Orthogroups.tsv> <dir_emapper> <go.obo> <salida_dir>
     dir_emapper contiene <código>.emapper.annotations (encabezados de proteína ">código|id").
"""
import sys, re, csv, math
from pathlib import Path
from collections import defaultdict, Counter

cafe, ogtsv, emdir, obo, out = map(Path, sys.argv[1:6])
out.mkdir(parents=True, exist_ok=True)
RAMAS = {"Humicolopsis_tronco": ["HumcLPSC", "HumcCH11"],
         "LPSC_terminal": ["HumcLPSC"],
         "Phacidiales_tronco": ["Allacyti", "Allalyco", "Xenopyri"]}
ALPHA, MIN_FRAC, MIN_TEST = 0.05, 0.5, 2

# ---------- CAFE5 ----------
def tabla(p):
    rows = list(csv.reader(open(p), delimiter="\t"))
    return [h.lstrip("#").strip() for h in rows[0]], {r[0]: r[1:] for r in rows[1:] if r}

def parse_newick(s):
    """Devuelve nodos (nombre, padre) de un newick; nombres de CAFE: 'Taxon<id>_n' o '<id>_n'."""
    nodes, stack, i, cur = [], [], 0, None
    while i < len(s):
        ch = s[i]
        if ch == "(":
            nodes.append({"name": None, "parent": stack[-1] if stack else None}); stack.append(len(nodes) - 1); i += 1
        elif ch in ",)":
            if ch == ")": cur = stack.pop()
            i += 1
            if ch == ")":
                m = re.match(r"[^,():;]*", s[i:]); nodes[cur]["name"] = m.group(0); i += len(m.group(0))
        elif ch == ";": break
        elif ch == ":":
            m = re.match(r":[0-9.eE+-]+", s[i:]); i += len(m.group(0))
        else:
            m = re.match(r"[^,():;]+", s[i:])
            nodes.append({"name": m.group(0), "parent": stack[-1]}); i += len(m.group(0))
    return nodes

def nodo_id(name):
    return int(re.search(r"<(\d+)>", name).group(1))

tree_line = next(l for l in open(next(cafe.glob("*_asr.tre"))) if l.strip().upper().startswith("TREE"))
nodes = parse_newick(tree_line.split("=", 1)[1].strip())
def tips_bajo(k):
    hijos = defaultdict(list)
    for j, n in enumerate(nodes):
        if n["parent"] is not None: hijos[n["parent"]].append(j)
    res, pila = set(), [k]
    while pila:
        x = pila.pop()
        if not hijos[x]: res.add(re.sub(r"<.*", "", nodes[x]["name"]))
        pila += hijos[x]
    return res
def mrca(taxa):
    cands = [j for j in range(len(nodes)) if set(taxa) <= tips_bajo(j)]
    return min(cands, key=lambda j: len(tips_bajo(j)))

pref = next(cafe.glob("*_change.tab")).name.replace("_change.tab", "")
fam_hdr, fam = tabla(cafe / f"{pref}_family_results.txt")
ch_hdr, change = tabla(cafe / f"{pref}_change.tab")
bp_hdr, bprob = tabla(cafe / f"{pref}_branch_probabilities.tab")
fam_p = {k: float(v[0]) for k, v in fam.items()}
col = lambda hdr, nid: next(i for i, h in enumerate(hdr[1:]) if re.search(rf"<{nid}>", h))

sets, resumen = {}, []
for rama, taxa in RAMAS.items():
    nid = nodo_id(nodes[mrca(taxa)]["name"])
    ic, ip = col(ch_hdr, nid), col(bp_hdr, nid)
    exp, con = [], []
    for og, v in change.items():
        try: pb = float(bprob[og][ip])
        except (KeyError, ValueError): continue
        if fam_p.get(og, 1) < ALPHA and pb < ALPHA:
            d = int(v[ic])
            (exp if d > 0 else con if d < 0 else []).append(og)
    sets[(rama, "expandidas")], sets[(rama, "contraidas")] = exp, con
    tot_ch = Counter(int(v[ic]) > 0 for v in change.values() if int(v[ic]) != 0)
    resumen.append([rama, nid, tot_ch[True], tot_ch[False], len(exp), len(con)])
with open(out / "resumen_ramas.tsv", "w") as f:
    f.write("rama\tnodo_cafe\tfamilias_con_ganancia\tfamilias_con_perdida\texpandidas_signif\tcontraidas_signif\n")
    for r in resumen: f.write("\t".join(map(str, r)) + "\n")
universo = set(change)

# ---------- anotación por ortogrupo ----------
go_name, go_ns = {}, {}
cur = None
for l in open(obo):
    if l.startswith("id: GO:"): cur = l[4:].strip()
    elif cur and l.startswith("name: "): go_name[cur] = l[6:].strip()
    elif cur and l.startswith("alt_id: "): go_name[l[8:].strip()] = go_name.get(cur, ""); go_ns[l[8:].strip()] = go_ns.get(cur, "?")
    elif cur and l.startswith("namespace: "): go_ns[cur] = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}.get(l[11:].strip(), "?")
    elif l.startswith("[Term]"): cur = None

gene_ann = {}
for f in emdir.glob("*.emapper.annotations"):
    for l in open(f):
        if l.startswith("#"): continue
        r = l.rstrip("\n").split("\t")
        if len(r) < 21: continue
        t = {"GO": [g for g in r[9].split(",") if g.startswith("GO:")],
             "KEGG_pathway": [p for p in r[12].split(",") if p.startswith("map")],
             "KEGG_ko": [k for k in r[11].split(",") if k.startswith("ko:")],
             "Pfam": [p for p in r[20].split(",") if p not in ("-", "")]}
        gene_ann[r[0]] = t
og_genes = {}
for r in csv.reader(open(ogtsv), delimiter="\t"):
    if r[0] == "Orthogroup": continue
    og_genes[r[0]] = [g.strip() for c in r[1:] for g in c.split(",") if g.strip()]

og_terms = defaultdict(dict)  # tipo -> og -> set(términos)
for og in universo:
    for tipo in ("GO", "KEGG_pathway", "KEGG_ko", "Pfam"):
        anot = [gene_ann[g][tipo] for g in og_genes.get(og, []) if g in gene_ann and gene_ann[g][tipo]]
        if not anot: continue
        cnt = Counter(t for a in anot for t in set(a))
        og_terms[tipo][og] = {t for t, c in cnt.items() if c / len(anot) >= MIN_FRAC}

# ---------- Fisher + BH ----------
def lchoose(n, k): return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
def fisher_mayor(a, K, n, N):  # P(X >= a), X ~ hipergeométrica(N, K, n)
    return min(1.0, sum(math.exp(lchoose(K, x) + lchoose(N - K, n - x) - lchoose(N, n)) for x in range(a, min(K, n) + 1)))
def bh(ps):
    m, orden = len(ps), sorted(range(len(ps)), key=lambda i: ps[i])
    q, prev = [0] * m, 1.0
    for rank, i in reversed(list(enumerate(orden, 1))):
        prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q

def nombre(tipo, t):
    if tipo == "GO": return f"{go_ns.get(t, '?')}: {go_name.get(t, '')}"
    return ""

todas = []
for (rama, sentido), ogs in sets.items():
    for tipo, ann in og_terms.items():
        fondo = set(ann)
        prueba = set(ogs) & fondo
        if len(prueba) < MIN_TEST: continue
        N, n = len(fondo), len(prueba)
        K = Counter(t for og in fondo for t in ann[og])
        a = Counter(t for og in prueba for t in ann[og])
        res = [(t, a[t], n, K[t], N, fisher_mayor(a[t], K[t], n, N)) for t in a if a[t] >= MIN_TEST]
        qs = bh([r[-1] for r in res]) if res else []
        for r, q in zip(res, qs):
            ogs_t = ",".join(sorted(og for og in prueba if r[0] in ann[og]))
            todas.append([rama, sentido, tipo, r[0], nombre(tipo, r[0]), r[1], r[2], r[3], r[4],
                          f"{(r[1]/r[2])/(r[3]/r[4]):.2f}", f"{r[5]:.3g}", f"{q:.3g}", ogs_t])
todas.sort(key=lambda r: (r[0], r[1], r[2], float(r[11])))
with open(out / "enriquecimiento.tsv", "w") as f:
    f.write("rama\tsentido\ttipo\ttermino\tnombre\tn_prueba_con_termino\tn_prueba\tn_fondo_con_termino\tn_fondo\tenriquecimiento\tp\tq_BH\tortogrupos\n")
    for r in todas: f.write("\t".join(map(str, r)) + "\n")
for (rama, sentido), ogs in sets.items():
    (out / f"ogs_{rama}_{sentido}.txt").write_text("\n".join(sorted(ogs)) + "\n")
print(open(out / "resumen_ramas.tsv").read())
sig = [r for r in todas if float(r[11]) < ALPHA]
print(f"términos con q < {ALPHA}: {len(sig)}")
for r in sig[:40]: print("  ", r[0], r[1], r[2], r[3], r[4], r[5], r[11])
