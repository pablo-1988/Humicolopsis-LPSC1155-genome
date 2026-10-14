#!/usr/bin/env bash
# A3 paso 3: ortología (OrthoFinder) + dominios Pfam (hmmsearch --cut_ga) + ubicación de referencias.
# Ejecutar con: caffeinate -i bash 03_orthofinder_dominios.sh
set -euo pipefail
W=~/hongo_comparativa
ARM=/Users/pablo/.local/share/mamba/envs/comparativa_arm/bin
OF=/Users/pablo/.local/share/mamba/envs/comparativa/bin/orthofinder
D="$(cd "$(dirname "$0")" && pwd)"   # carpeta de este script
export PATH=$ARM:/Users/pablo/.local/share/mamba/envs/comparativa/bin:$PATH   # diamond/mafft arm64 primero
cd $W

# 1. proteomas limpios (1 isoforma no requerida: NCBI protein.faa puede traer isoformas -> se conserva la más larga por gen)
mkdir -p A3/prot
python3 - <<'EOF'
import re, pathlib
W = pathlib.Path.home() / "hongo_comparativa"
for f in sorted((W / "faa").glob("*.faa")):
    code = f.stem
    gff = W / "gff" / f"{code}.gff"
    prot2gene = {}
    if gff.exists() and code != "HumcLPSC":
        for l in open(gff):
            if "\tCDS\t" in l:
                a = l.rstrip().split("\t")[8]
                p = re.search(r"protein_id=([^;]+)", a); g = re.search(r"(?:locus_tag|gene)=([^;]+)", a)
                if p: prot2gene[p.group(1)] = g.group(1) if g else p.group(1)
    seqs, n = {}, None
    for l in open(f):
        if l.startswith(">"): n = l[1:].split()[0]; seqs[n] = []
        else: seqs[n].append(l.strip())
    best = {}
    for n, s in seqs.items():
        s = "".join(s).rstrip("*")
        g = prot2gene.get(n, n)
        if g not in best or len(s) > len(best[g][1]): best[g] = (n, s)
    with open(W / "A3" / "prot" / f"{code}.faa", "w") as out:
        for g, (n, s) in best.items():
            out.write(f">{code}|{n}\n{s}\n")
    print(code, len(seqs), "->", len(best))
EOF

# 2. OrthoFinder
[ -d A3/prot/OrthoFinder ] || $OF -f A3/prot -t 11 -a 4 -S diamond -M msa -A mafft -T fasttree > A3/orthofinder.log 2>&1

# 3. Dominios Pfam
cat A3/prot/*.faa > A3/all_prot.faa
$ARM/hmmsearch --cut_ga --cpu 11 --noali --domtblout A3/pfam_familias.domtbl "$D/hmm/familias.hmm" A3/all_prot.faa > /dev/null

# 4. Referencias -> proteínas del set (diamond)
$ARM/diamond makedb --in A3/all_prot.faa -d A3/all_prot --quiet
$ARM/diamond blastp -q "$D/refs/referencias.faa" -d A3/all_prot --very-sensitive -k 500 -e 1e-5 \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore -o A3/refs_vs_all.tsv --quiet
echo "listo"
