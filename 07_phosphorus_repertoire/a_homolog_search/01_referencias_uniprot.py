"""Descarga proteínas de referencia (Swiss-Prot/UniProt) con función caracterizada
para las familias de P y ácidos orgánicos. Salida: refs/referencias.faa + refs/referencias.tsv
Organismos: S. cerevisiae (559292), N. crassa (367110), A. niger (5061/425011), S. sclerotiorum (665079),
B. cinerea (332648)."""
import csv, time, urllib.request, urllib.parse
from pathlib import Path

OUT = Path(__file__).parent / "refs"
OUT.mkdir(exist_ok=True)

# familia, categoría, gen, query UniProt
REFS = [
    # Regulación / señalización PHO
    ("PHO85", "PHO_regulacion", "PHO85", "gene_exact:PHO85 AND organism_id:559292"),
    ("PHO80", "PHO_regulacion", "PHO80", "gene_exact:PHO80 AND organism_id:559292"),
    ("PHO81", "PHO_regulacion", "PHO81", "gene_exact:PHO81 AND organism_id:559292"),
    ("PHO4", "PHO_regulacion", "PHO4", "gene_exact:PHO4 AND organism_id:559292"),
    ("PHO4", "PHO_regulacion", "nuc-1", "gene_exact:nuc-1 AND organism_id:367110"),
    ("PHO2", "PHO_regulacion", "PHO2", "accession:P07269"),
    ("PHO80", "PHO_regulacion", "preg", "gene_exact:preg AND organism_id:367110"),
    ("PHO81", "PHO_regulacion", "nuc-2", "gene_exact:nuc-2 AND organism_id:367110"),
    ("PHO85", "PHO_regulacion", "mdk-1", "accession:Q7RVR4"),
    ("KCS1", "PHO_regulacion", "KCS1", "gene_exact:KCS1 AND organism_id:559292"),
    ("VIP1", "PHO_regulacion", "VIP1", "gene_exact:VIP1 AND organism_id:559292"),
    ("IPK1", "PHO_regulacion", "IPK1", "gene_exact:IPK1 AND organism_id:559292"),
    ("IPK2", "PHO_regulacion", "ARG82", "gene_exact:ARG82 AND organism_id:559292"),
    # Transporte de Pi
    ("PHO84", "Transporte_Pi", "PHO84", "gene_exact:PHO84 AND organism_id:559292"),
    ("PHO84", "Transporte_Pi", "pho-5", "gene_exact:pho-5 AND organism_id:367110"),
    ("PHO89", "Transporte_Pi", "PHO89", "gene_exact:PHO89 AND organism_id:559292"),
    ("PHO89", "Transporte_Pi", "pho-4", "gene_exact:pho-4 AND organism_id:367110"),
    ("PHO87/90/91", "Transporte_Pi", "PHO87", "gene_exact:PHO87 AND organism_id:559292"),
    ("PHO87/90/91", "Transporte_Pi", "PHO90", "gene_exact:PHO90 AND organism_id:559292"),
    ("PHO87/90/91", "Transporte_Pi", "PHO91", "gene_exact:PHO91 AND organism_id:559292"),
    ("PHO88", "Transporte_Pi", "PHO88", "gene_exact:PHO88 AND organism_id:559292"),
    ("PHO86", "Transporte_Pi", "PHO86", "gene_exact:PHO86 AND organism_id:559292"),
    ("GIT1", "Transporte_Pi", "GIT1", "gene_exact:GIT1 AND organism_id:559292"),
    ("PIC2/MIR1", "Transporte_Pi", "MIR1", "gene_exact:MIR1 AND organism_id:559292"),
    # Poli-P / homeostasis vacuolar
    ("VTC", "PolyP", "VTC1", "gene_exact:VTC1 AND organism_id:559292"),
    ("VTC", "PolyP", "VTC2", "gene_exact:VTC2 AND organism_id:559292"),
    ("VTC", "PolyP", "VTC3", "gene_exact:VTC3 AND organism_id:559292"),
    ("VTC", "PolyP", "VTC4", "gene_exact:VTC4 AND organism_id:559292"),
    ("PPX1", "PolyP", "PPX1", "gene_exact:PPX1 AND organism_id:559292"),
    ("PPN1", "PolyP", "PPN1", "gene_exact:PPN1 AND organism_id:559292"),
    ("PPN2", "PolyP", "PPN2", "gene_exact:PPN2 AND organism_id:559292"),
    ("SYG1", "PolyP", "SYG1", "gene_exact:SYG1 AND organism_id:559292"),
    # Fosfatasas y fitasas
    ("PHO5-like HAP", "Fosfatasas", "PHO5", "gene_exact:PHO5 AND organism_id:559292"),
    ("PHO5-like HAP", "Fosfatasas", "PHO3", "gene_exact:PHO3 AND organism_id:559292"),
    ("PHO5-like HAP", "Fosfatasas", "PHO11", "gene_exact:PHO11 AND organism_id:559292"),
    ("Fitasa HAP", "Fitasas", "phyA", "gene_exact:phyA AND taxonomy_id:5061 AND reviewed:true"),
    ("Fitasa HAP", "Fitasas", "phyB", "gene_exact:phyB AND taxonomy_id:5061 AND reviewed:true"),
    ("Fosfatasa alcalina", "Fosfatasas", "PHO8", "gene_exact:PHO8 AND organism_id:559292"),
    ("PHO13 HAD", "Fosfatasas", "PHO13", "gene_exact:PHO13 AND organism_id:559292"),
    ("Fosfatasa alcalina", "Fosfatasas", "pho-2", "gene_exact:pho-2 AND organism_id:367110"),
    ("PAP", "Fosfatasas", "pho-3", "gene_exact:pho-3 AND organism_id:367110"),
    ("GDE1", "Fosfatasas", "GDE1", "gene_exact:GDE1 AND organism_id:559292"),
    ("Nucleasa NUC1", "P_organico", "NUC1", "gene_exact:NUC1 AND organism_id:559292"),
    # Ácidos orgánicos
    ("OAH", "Acidos_organicos", "oahA", "gene_exact:oahA AND taxonomy_id:5061"),
    ("OAH", "Acidos_organicos", "Ss-oah1", "accession:A7ESB3"),
    ("OAH", "Acidos_organicos", "Bcoah", "accession:A0A384JYA5"),
    ("ICL", "Acidos_organicos", "ICL1", "gene_exact:ICL1 AND organism_id:559292"),
    ("MCL", "Acidos_organicos", "ICL2", "gene_exact:ICL2 AND organism_id:559292"),
    ("Oxalato descarboxilasa", "Acidos_organicos", "oxdC", "protein_name:oxalate decarboxylase AND reviewed:true AND taxonomy_id:4751"),
    ("Oxalil-CoA/oxalato", "Acidos_organicos", "PCS60", "gene_exact:PCS60 AND organism_id:559292"),
    ("Citrato sintasa", "Acidos_organicos", "CIT1", "gene_exact:CIT1 AND organism_id:559292"),
    ("Citrato sintasa", "Acidos_organicos", "CIT2", "gene_exact:CIT2 AND organism_id:559292"),
    ("ATP-citrato liasa", "Acidos_organicos", "ACL1", "accession:O93988"),
    ("ATP-citrato liasa", "Acidos_organicos", "ACL2", "accession:O13907"),
    ("Malato DH", "Acidos_organicos", "MDH1", "gene_exact:MDH1 AND organism_id:559292"),
    ("Malato DH", "Acidos_organicos", "MDH2", "gene_exact:MDH2 AND organism_id:559292"),
    ("Malato DH", "Acidos_organicos", "MDH3", "gene_exact:MDH3 AND organism_id:559292"),
    ("Fumarasa", "Acidos_organicos", "FUM1", "gene_exact:FUM1 AND organism_id:559292"),
    ("Piruvato carboxilasa", "Acidos_organicos", "PYC1", "gene_exact:PYC1 AND organism_id:559292"),
    ("Enzima málica", "Acidos_organicos", "MAE1", "gene_exact:MAE1 AND organism_id:559292"),
    ("Transp. citrato mitocondrial", "Transporte_acidos", "CTP1", "gene_exact:CTP1 AND organism_id:559292"),
    ("Transp. dicarboxilato mitocondrial", "Transporte_acidos", "DIC1", "gene_exact:DIC1 AND organism_id:559292"),
    ("Transp. oxaloacetato mitocondrial", "Transporte_acidos", "OAC1", "gene_exact:OAC1 AND organism_id:559292"),
    ("Exportador citrato CexA", "Transporte_acidos", "cexA", "gene_exact:cexA AND taxonomy_id:5061"),
    ("Transp. malato TDT (Mae1)", "Transporte_acidos", "mae1", "gene_exact:mae1 AND organism_id:284812"),
    ("Transp. carboxilato JEN1", "Transporte_acidos", "JEN1", "gene_exact:JEN1 AND organism_id:559292"),
    ("Transp. dicarboxilato ADY2", "Transporte_acidos", "ADY2", "gene_exact:ADY2 AND organism_id:559292"),
]

URL = "https://rest.uniprot.org/uniprotkb/search?format=fasta&size=3&query="
faa, tsv = [], [["familia", "categoria", "gen", "accession", "reviewed", "header"]]
for fam, cat, gene, q in REFS:
    for attempt in range(3):
        try:
            # preferir reviewed; si no hay, aceptar TrEMBL
            txt = ""
            for extra in (" AND reviewed:true", ""):
                txt = urllib.request.urlopen(URL + urllib.parse.quote(q + extra), timeout=60).read().decode()
                if txt.strip():
                    break
            break
        except Exception as e:
            time.sleep(3)
    recs = [r for r in txt.split(">") if r.strip()]
    if not recs:
        print(f"NO ENCONTRADO: {gene} ({q})")
        tsv.append([fam, cat, gene, "NA", "", ""])
        continue
    r = recs[0]
    head, seq = r.split("\n", 1)
    acc = head.split("|")[1]
    faa.append(f">REF|{gene}|{acc}|{fam.replace(' ', '_')}\n{seq.strip()}\n")
    tsv.append([fam, cat, gene, acc, head.startswith("sp"), head])
    time.sleep(0.3)

(OUT / "referencias.faa").write_text("".join(faa))
with (OUT / "referencias.tsv").open("w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(tsv)
print(f"{len(faa)} referencias descargadas")
