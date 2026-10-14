"""A7.4a — Entrada de CAFE5 desde Orthogroups.GeneCount.tsv.

Quita las referencias funcionales (S. cerevisiae, A. niger) y aplica los filtros del tutorial de CAFE5:
  - fuera los ortogrupos presentes en < 2 especies (sin información sobre cambios en el árbol);
  - los ortogrupos con >= 100 copias en alguna especie van a un archivo aparte (se corren con lambda fijo).
Uso: python3 05a_preparar_cafe.py <GeneCount.tsv> <salida_dir> "Sacccere Aspenige"
"""
import sys, csv
from pathlib import Path

gc, out, excl = sys.argv[1], Path(sys.argv[2]), set(sys.argv[3].split())
out.mkdir(parents=True, exist_ok=True)
rows = list(csv.reader(open(gc), delimiter="\t"))
head = rows[0]
cols = [i for i, h in enumerate(head) if h not in excl and h not in ("Orthogroup", "Total")]
n = {"main": 0, "large": 0, "<2 especies": 0}
with open(out / "cafe_input.tsv", "w") as fm, open(out / "cafe_large.tsv", "w") as fl:
    hdr = "Desc\tFamily ID\t" + "\t".join(head[i] for i in cols) + "\n"
    fm.write(hdr); fl.write(hdr)
    for r in rows[1:]:
        c = [int(r[i]) for i in cols]
        if sum(x > 0 for x in c) < 2:
            n["<2 especies"] += 1; continue
        line = f"(null)\t{r[0]}\t" + "\t".join(map(str, c)) + "\n"
        if max(c) >= 100:
            fl.write(line); n["large"] += 1
        else:
            fm.write(line); n["main"] += 1
print("especies:", len(cols), "|", n)
