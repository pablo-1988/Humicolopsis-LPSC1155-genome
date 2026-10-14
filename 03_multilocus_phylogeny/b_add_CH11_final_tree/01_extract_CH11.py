"""Extrae los 9 marcadores de H. cephalosporioides CH-1.1 (GCA_050941815.1) desde su genoma,
usando como query los marcadores de LPSC 1155 (filogenia_v2_revision/02_ncbi_raw/_T1)."""
import subprocess, os
from pathlib import Path
from Bio import SeqIO
P = Path(__file__).resolve().parents[3]
T1 = P / "filogenia_v2_revision/02_ncbi_raw/_T1"
OUT = Path(__file__).parent / "_CH11"; OUT.mkdir(exist_ok=True)
W = Path.home() / "hongo_comparativa"
G = W / "fna/HumcCH11.fna"; DB = W / "blastdb/CH11"
B = "/Users/pablo/.local/share/mamba/envs/filogenia/bin"
if not Path(f"{DB}.nsq").exists() and not Path(f"{DB}.00.nsq").exists():
    DB.parent.mkdir(exist_ok=True)
    subprocess.run([f"{B}/makeblastdb", "-in", G, "-dbtype", "nucl", "-out", DB], check=True, capture_output=True)
genome = SeqIO.to_dict(SeqIO.parse(G, "fasta"))
for mk in ["SSU", "LSU", "ITS", "mtSSU", "TUB2", "EF1a", "MCM7", "RPB1", "RPB2"]:
    q = T1 / f"T1__{mk}.fasta"
    qlen = len(next(SeqIO.parse(q, "fasta")).seq)
    r = subprocess.run([f"{B}/blastn", "-query", q, "-db", DB, "-evalue", "1e-20", "-max_target_seqs", "5",
                        "-outfmt", "6 sseqid sstart send pident length qstart qend bitscore"], capture_output=True, text=True)
    hits = [l.split("\t") for l in r.stdout.strip().split("\n") if l]
    if not hits: print(mk, "sin hit"); continue
    sid, s, e, pid, L, qs, qe, bs = max(hits, key=lambda h: float(h[7]))
    s, e = int(s), int(e)
    sub = genome[sid].seq[s - 1:e] if s <= e else genome[sid].seq[e - 1:s].reverse_complement()
    (OUT / f"CH11__{mk}.fasta").write_text(f">Humicolopsis_cephalosporioides_CH-1.1\n{sub}\n")
    print(f"{mk}: {sid}:{s}-{e} id={pid}% cobertura_query={(int(qe)-int(qs)+1)/qlen:.0%} ({len(sub)} pb)")
