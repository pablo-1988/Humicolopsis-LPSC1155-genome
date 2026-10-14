"""CDS (sin intrones) y proteína de los marcadores codificantes de H. cephalosporioides CBS 645.76 (ex-tipo;
sin genoma). Alineamiento proteína->nt con miniprot usando el proteoma de LPSC 1155 como molde.
Salida: CBS645_cds.fna, CBS645_prot.faa (id = marcador)."""
import subprocess
from pathlib import Path
from Bio.Seq import Seq
from Bio import SeqIO
P = Path(__file__).resolve().parents[3]
HERE = Path(__file__).parent
MP = "/Users/pablo/.local/share/mamba/envs/comparativa_arm/bin/miniprot"
PROT = Path.home() / "hongo_comparativa/faa/HumcLPSC.faa"
SRC = {"MCM7": P / "filogenia_v2_revision/02_ncbi_raw/_cache/Humicolopsis_cephalosporioides__MCM7.fasta",
       "RPB1": P / "filogenia_v2_revision/02_ncbi_raw/_cache/Humicolopsis_cephalosporioides__RPB1.fasta",
       "RPB2": P / "filogenia_v2_revision/02_ncbi_raw/_cache/Humicolopsis_cephalosporioides__RPB2.fasta",
       "TUB2": P / "filogenia_v2_revision/02_ncbi_raw/_cache/Humicolopsis_cephalosporioides__TUB2.fasta",
       "EF1a": P / "filogenia_v2_revision/01_query_T1/TEF1_CBS645.fasta"}
cds_out, prot_out = [], []
for mk, f in SRC.items():
    rec = next(SeqIO.parse(f, "fasta")); nt = str(rec.seq).upper()
    tgt = HERE / f"_{mk}.fna"; tgt.write_text(f">{mk}\n{nt}\n")
    gff = subprocess.run([MP, "--gff", "-I", "-t", "4", tgt, PROT], capture_output=True, text=True).stdout
    best, segs, cur = None, {}, None
    for l in gff.splitlines():
        f_ = l.split("\t")
        if len(f_) < 9: continue
        if f_[2] == "mRNA":
            cur = f_[8].split("ID=")[1].split(";")[0]
            sc = float(f_[5]); tgtp = f_[8].split("Target=")[1].split()[0]
            if best is None or sc > best[1]: best = (cur, sc, f_[6], tgtp)
        elif f_[2] == "CDS":
            segs.setdefault(f_[8].split("Parent=")[1].split(";")[0], []).append((int(f_[3]), int(f_[4])))
    if best is None:
        print(f"{mk}: sin alineamiento proteico (fragmento corto con intrones) -> no se usa a nivel proteína")
        tgt.unlink(); continue
    mid, sc, strand, tp = best
    s = sorted(segs[mid])
    cds = "".join(nt[a - 1:b] for a, b in s)
    if strand == "-": cds = str(Seq(cds).reverse_complement())
    aa = str(Seq(cds[: len(cds) // 3 * 3]).translate())
    cds_out.append(f">{mk} {rec.id} template={tp}\n{cds}\n"); prot_out.append(f">{mk}\n{aa}\n")
    print(f"{mk}: {rec.id} {len(nt)} nt -> CDS {len(cds)} nt ({len(s)} exones), {aa.count('*')} stops, molde {tp}")
    tgt.unlink()
(HERE / "CBS645_cds.fna").write_text("".join(cds_out)); (HERE / "CBS645_prot.faa").write_text("".join(prot_out))
