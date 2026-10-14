"""Inventario de genomas de Leotiomycetes en NCBI con orden/familia y métricas de calidad.

Entrada: salida de `datasets summary genome taxon Leotiomycetes --as-json-lines`
Salida:  inventario_leotiomycetes_ncbi.tsv (un ensamble por fila)
"""
import json, subprocess, sys, csv
from pathlib import Path

DATASETS = "/Users/pablo/.local/share/mamba/envs/filogenia/bin/datasets"
jsonl = Path(sys.argv[1])
out = Path(__file__).parent / "inventario_leotiomycetes_ncbi.tsv"

recs = [json.loads(l) for l in jsonl.open()]
taxids = sorted({str(r["organism"]["tax_id"]) for r in recs})

# Clasificación taxonómica (orden, familia) por lotes
lineage = {}
for i in range(0, len(taxids), 200):
    res = subprocess.run([DATASETS, "summary", "taxonomy", "taxon", *taxids[i:i + 200], "--as-json-lines"],
                         capture_output=True, text=True, check=True)
    for l in res.stdout.splitlines():
        t = json.loads(l)["taxonomy"]
        c = t.get("classification", {})
        lineage[str(t["tax_id"])] = (c.get("order", {}).get("name", "NA"), c.get("family", {}).get("name", "NA"))

cols = ["accession", "organism", "strain", "order", "family", "level", "size_Mb", "contigs",
        "contig_N50_kb", "scaffold_N50_kb", "annotated", "protein_coding", "annot_provider",
        "refseq_category", "busco_lineage", "busco_complete", "release_date"]
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in recs:
        a, s, o = r["assembly_info"], r["assembly_stats"], r["organism"]
        ann = r.get("annotation_info") or {}
        busco = ann.get("busco", {})
        order, fam = lineage.get(str(o["tax_id"]), ("NA", "NA"))
        w.writerow([
            r["accession"], o["organism_name"],
            (o.get("infraspecific_names") or {}).get("strain", ""),
            order, fam, a.get("assembly_level", ""),
            round(int(s["total_sequence_length"]) / 1e6, 1), s.get("number_of_contigs", ""),
            round(s.get("contig_n50", 0) / 1e3, 1), round(s.get("scaffold_n50", 0) / 1e3, 1),
            "yes" if ann else "no", ann.get("stats", {}).get("gene_counts", {}).get("protein_coding", ""),
            ann.get("provider", ""), r.get("assembly_info", {}).get("refseq_category", ""),
            busco.get("busco_lineage", ""), busco.get("complete", ""), a.get("release_date", ""),
        ])
print(f"{len(recs)} ensambles -> {out}")
