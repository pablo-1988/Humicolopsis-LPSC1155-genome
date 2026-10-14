"""Proteomas por homología (miniprot) para los ensambles sin anotación.

Evidencia proteica: proteoma de LPSC 1155 + Glarea lozoyensis + Pseudogymnoascus destructans
+ Botrytis cinerea + referencias Swiss-Prot de P/ácidos orgánicos.
Por cada locus del genoma se conserva el alineamiento de mayor puntaje (sin solapamiento en la misma hebra).
Salida: ~/hongo_comparativa/faa/<code>.faa  (marcado en taxa_map como 'miniprot')
"""
import subprocess, re, sys
from pathlib import Path

W = Path.home() / "hongo_comparativa"
E = "/Users/pablo/.local/share/mamba/envs/comparativa_arm/bin"
REFS = Path(__file__).parent / "refs" / "referencias.faa"
EVID = ["HumcLPSC", "Glarlozo", "Pseudest", "Botrcine"]
SIN_ANOT = ["Allacyti", "Allalyco", "Xenopyri", "Leotlubr", "Cyclminu", "Chaeraph", "HumcCH11"]
THREADS = sys.argv[1] if len(sys.argv) > 1 else "10"

mp = W / "miniprot"; mp.mkdir(exist_ok=True)
evid = mp / "evidencia.faa"
if not evid.exists():
    with evid.open("w") as out:
        for code in EVID:
            for line in (W / "faa" / f"{code}.faa").open():
                out.write(line.replace(">", f">{code}|", 1) if line.startswith(">") else line)
        out.write(REFS.read_text())

for code in SIN_ANOT:
    faa = W / "faa" / f"{code}.faa"
    if faa.exists():
        print("skip", code); continue
    idx = mp / f"{code}.mpi"
    subprocess.run([f"{E}/miniprot", "-t", THREADS, "-d", idx, W / "fna" / f"{code}.fna"], check=True,
                   stderr=subprocess.DEVNULL)
    gff = mp / f"{code}.gff"
    with gff.open("w") as fh:
        subprocess.run([f"{E}/miniprot", "-t", THREADS, "--gff", "--trans", "-I", "--outs=0.9", idx, evid],
                       stdout=fh, stderr=subprocess.DEVNULL, check=True)

    # parsear mRNA (coordenadas, score) + traducción (##STA)
    hits, seq = [], None
    for line in gff.open():
        if line.startswith("##STA"):
            seq = line.split("\t")[1].strip()
        elif "\tmRNA\t" in line:
            f = line.rstrip().split("\t")
            hits.append((f[0], int(f[3]), int(f[4]), f[6], float(f[5]), seq.replace("*", "")))
    # mejor alineamiento por locus: ordenar por score y rechazar solapados
    hits.sort(key=lambda h: -h[4])
    kept, occ = [], {}
    for h in hits:
        key = (h[0], h[3])
        if any(not (h[2] < s or h[1] > e) for s, e in occ.get(key, [])):
            continue
        occ.setdefault(key, []).append((h[1], h[2]))
        kept.append(h)
    kept.sort(key=lambda h: (h[0], h[1]))
    with faa.open("w") as out:
        for i, h in enumerate(kept, 1):
            out.write(f">{code}_mp{i:05d} {h[0]}:{h[1]}-{h[2]}({h[3]})\n{h[5]}\n")
    print(f"{code}: {len(kept)} loci proteicos")
