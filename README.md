# Genome of *Humicolopsis cephalosporioides* LPSC 1155: analysis scripts

Code and data (alignments, trees and CAFE5 files) for the analyses in the manuscript on the genome, evolutionary position and phosphorus-related gene repertoire of *Humicolopsis cephalosporioides* LPSC 1155. Script folders follow the order of the Materials and Methods; figure scripts are not included.

Script comments are in Spanish.

## Structure

| Folder | Methods section | What it does |
|---|---|---|
| `01_assembly_annotation/` | Genome assembly, polishing and scaffolding; Genome annotation | Two rounds of Pilon polishing (BWA-MEM, MAPQ ≥ 20) and the Funannotate v1.8.15 protocol (container `nextgenusfs/funannotate:v1.8.15`: scaffold renaming, tantan masking, `predict` with BUSCO sordariomycetes_odb10 / *N. crassa* seed). `00_qc_assembly_pipeline_template.txt` lists the fastp → SPAdes (`--careful`) → QUAST → BUSCO commands (generic template; its Funannotate block was not the one used). |
| `02_genome_selection/` | Table S3 | Inventory of Leotiomycetes genomes in NCBI used to select the comparative dataset. |
| `03_multilocus_phylogeny/` | Multilocus phylogeny | (a) Retrieval of the nine loci from GenBank (type material preferred) and extraction of the LPSC 1155 loci from its genome (BLASTn/tBLASTn); (b) extraction of the CH-1.1 loci and MAFFT → trimAl → concatenation → IQ-TREE (MFP+MERGE, 1,000 UFBoot). |
| `04_phylogenomics_14genes/` | Genome-scale phylogenies, 14-gene dataset | Genome download (NCBI datasets), CH-1.1 preparation, CBS 645.76 CDS delimitation with miniprot, gene extraction, codon-aware alignment and IQ-TREE. |
| `05_phylogenomics_BUSCO/` | Genome-scale phylogenies, BUSCO dataset | compleasm (leotiomycetes_odb10) on all assemblies, locus selection and alignments, partitioned IQ-TREE 3 (`MFP+MERGE -rcluster 10`), gene trees, ASTRAL-IV and gCF/sCF. |
| `06_gene_family_evolution/` | Gene family evolution and functional enrichment | Re-annotation of seven genomes with Funannotate (Apptainer), OrthoFinder, ultrametric tree (ape `chronos`), CAFE5 (error model + gamma), eggNOG-mapper and Fisher/BH enrichment, robustness run on Funannotate genomes only. `run_A7.sh` runs the whole pipeline; paths are set in `00_config.sh`; `entorno_a7.yml` is the conda environment; `LEEME_A7.md` describes each step. |
| `07_phosphorus_repertoire/` | Phosphorus and organic acid gene repertoire | (a) Reference proteins from UniProt, DIAMOND + hmmsearch (Pfam `--cut_ga`) against the proteomes; (b) gene trees per family (MAFFT, trimAl, IQ-TREE), orthologous-clade assignment, final copy-number table and placement of the re-annotated proteomes. |
| `08_secretome/` | Secretome prediction | DeepSig + TMbed predictions and secretome filtering (no TM outside the signal peptide, no HDEL/KDEL). |
| `09_humicolopsis_records_Chile/` | Supplementary Fig./Table SX | Re-analysis of *Humicolopsis* ASVs in the *Nothofagus pumilio* rhizosphere dataset of Saona et al. 2025 (*ISME Commun* 5:ycaf010); writes the supplementary table and figure. |

## Data (`data/`)

Alignments, trees and CAFE5 files reported in the manuscript.

| Folder | Contents |
|---|---|
| `01_multilocus_9loci/` | Per-locus MAFFT alignments (`alignments/`) and trimAl-trimmed alignments (`trimmed/`) of the nine loci; concatenated matrix, partition file and taxon × locus matrix (`concatenated/`); IQ-TREE trees (`trees/`): unconstrained tree (`Humicolopsis_v3.*`, IQ-TREE 2.4.0) and final tree with the outgroup constrained to be monophyletic (`Humicolopsis_v3_constrained.*`, IQ-TREE 3.1.3; constraint in `constraint_outgroup.nwk`). `.treefile` = ML tree; `.contree` = UFBoot consensus; `.iqtree` = full report; `.best_scheme.nex` = selected partitioning scheme. |
| `02_phylogenomics_14genes/` | Trimmed alignments of the 14 loci (codon-aware for protein-coding genes), concatenated matrix (28,618 bp, 19 taxa) with partitions, and IQ-TREE 2.4.0 tree. |
| `03_phylogenomics_BUSCO/` | `matrix/`: concatenated amino-acid matrix of 2,742 leotiomycetes_odb10 loci for 23 taxa (gzip), per-locus partitions, partial taxa (CBS 645.76). `partitioned_tree/`: partitioned IQ-TREE 3.1.3 tree (2,346 merged partitions; best scheme and models), with gCF/sCF (`part_concord.cf.*`). `single_model_tree/`: Q.YEAST+F+I+R5 tree, model selection on 300 loci, concordance factors. `gene_trees_ASTRAL/`: 2,742 gene trees (gzip) and the ASTRAL-IV species tree. |
| `04_gene_family_trees/` | For each of the 19 P-related families: input proteins (`.faa`), trimmed alignment (`.trim.aln`), ML tree (`.treefile`), UFBoot consensus (`.contree`) and report (`.iqtree`). `trees_original_proteomes/` = trees used for curation (24 proteomes); `trees_reannotated_proteomes/` = trees rebuilt with the seven Funannotate re-annotated proteomes, used to place their genes in the curated clades. `curated_tables/`: curated genes, copy numbers and final repertoire tables. |
| `05_CAFE5/` | `input/`: OrthoFinder gene counts (24 proteomes), CAFE5 input (15,346 orthogroups, 22 genomes), large families analysed separately (≥ 100 copies, 17 orthogroups) and ultrametric tree. `output/`: results of the selected gamma model (`Gamma_family_results.txt`, `Gamma_change.tab`, `Gamma_asr.tre`) and likelihoods of all base and gamma runs (`cafe_modelos.txt`). `robustness_funannotate_only/`: same analysis restricted to the eight genomes annotated with Funannotate. |

Taxon codes used in alignments and trees (e.g. `HumcLPSC`, `Allacyti`) are listed with species, strain and assembly accession in `data/taxon_codes.tsv`.

## Main software

fastp, SPAdes, Pilon 1.24, BWA-MEM 0.7.19, samtools 1.22.1, RagTag 2.1.0, minimap2 2.30, QUAST 5.2.0, BUSCO 6.0.0, compleasm 0.2.9, miniprot 0.18, Funannotate 1.8.15, eggNOG-mapper 2.1.13, HMMER 3, antiSMASH 8.0.4, BLAST+ 2.16, MAFFT 7.525/7.526, trimAl 1.4.1/1.5, IQ-TREE 2.4.0 and 3.1.3, ASTRAL-IV 1.25.4.8, OrthoFinder 2.5.5, DIAMOND 2.2.8, CAFE5 5.1.0, R/ape 5.8.1, DeepSig 0.9, TMbed 1.0.2, Python 3 with Biopython and pandas.

## Before running

- Scripts were run on macOS (Apple Silicon) and on a Linux server. Several scripts contain absolute paths to conda environments (`/Users/pablo/.local/share/mamba/envs/...`, `/home/ubuntu/miniforge3/envs/...`) and to a working directory outside the repository (`~/hongo_comparativa`, holding genomes `fna/`, proteomes `faa/` and intermediate files). Adjust these paths before re-running.
- Genome assemblies and proteomes are not included. Public genomes can be downloaded with `04_phylogenomics_14genes/10_fetch_genomes.py`. The accessions are listed in Table S3.
- `03_multilocus_phylogeny/.../01_collect_ncbi.py` needs a contact e-mail for NCBI Entrez (`Entrez.email`).
