---
output:
  word_document: default
  pdf_document: default
---
# Shotgun Metagenomics Workflow

Record of commands used at each step. Sample IDs are written generically as `SAMPLE01`–`SAMPLE24`.
Samples are split into two groups: facility A and facility B. 

Every step writes to its own sequentially numbered output directory.

**Directory map**

| Directory | Contents |
|---|---|
| `1_raw_data/` | Input FASTQ files |
| `2_qc_reads/` | Quality-filtered reads |
| `3_no_mouse_reads/` | Host-depleted reads |
| `4_assembly/` | Per-sample assemblies (metaSPAdes) |
| `5_coassembly/` | Per-facility co-assemblies (MEGAHIT) |
| `6_prodigal/` | Gene and protein predictions |
| `7_gene_catalog/` | Non-redundant mix-assembly gene catalog |
| `8_mapped_bam/` | Reads mapped back to gene catalog |
| `9_coverm/` | Gene count and TPM tables |
| `10_annotations/` | AsgeneDB, CAZyme, eggNOG annotations |
| `11_initial_binning/` | MetaBAT2 / MaxBin2 / CONCOCT bins |
| `12_bin_refinement/` | Refined bin sets |
| `13_quant_bins/` | Bin abundance tables |
| `15_bin_classification/` | Bin taxonomy |
| `16_MAG_function/` | MAG gene calls and annotation |
| `17_AsDB_bins/` | Per-MAG arsenic gene annotation |


---

## Part 1. Quality filtering and host read removal


### 1.1 Quality filter

```bash

for F in 1_raw_data/*_1.fastq; do
    R=${F%_*}_2.fastq
    BASE=${F##*/}
    SAMPLE=${BASE%_*}
    metawrap read_qc -1 $F -2 $R -t 24 --skip-bmtagger -o 2_qc_reads/$SAMPLE &
done
wait


```

### 1.2 Build host reference index

```bash


wget ftp://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

bowtie2-build \
    Mus_musculus.GRCm39.dna.primary_assembly.fa \
    mouse_index


```

### 1.3 Remove host reads

Reads are aligned to the host reference and only unmapped pairs are retained.

```bash
for i in $(seq -w 1 24); do
    S=SAMPLE${i}
    echo "Processing $S ..."

    bowtie2 -p 24 \
        -x ref/mouse_index \
        -1 2_qc_reads/${S}/final_pure_reads_1.fastq \
        -2 2_qc_reads/${S}/final_pure_reads_2.fastq \
        --very-sensitive-local \
        --un-conc-gz 3_no_mouse_reads/${S}_no_mouse.fastq.gz \
        -S /dev/null

    echo "Finished $S"
done
```

---

## Part 2. Assembly

### 2.1 Per-sample assembly (metaSPAdes)

```bash

#note file names need to be changed prior to input to avoid error from metaspades
#bowtie2 --un-conc-gz writes _no_mouse.fastq.1.gz / .2.gz, which must be renamed to _no_mouse_1.fastq.gz / _2.fastq.gz

for i in $(seq -w 1 24); do
    S=SAMPLE${i}
    echo "Starting $S at $(date)..."

    metawrap assembly \
        -1 3_no_mouse_reads/${S}_no_mouse_1.fastq.gz \
        -2 3_no_mouse_reads/${S}_no_mouse_2.fastq.gz \
        -m 200 -t 24 -l 500 \
        --metaspades \
        -o 4_assembly/${S}

    echo "Finished $S at $(date)"
done
```

### 2.2 Per-facility co-assembly (MEGAHIT)

Reads from all samples within a facility are pooled and co-assembled.

```bash
#command run twice, once for facility A, once for facility B

# Pool reads for this facility
cat 3_no_mouse_reads/SAMPLE{01,02,03,04,09,10,11,12,17,18,19,20}_no_mouse_1.fastq.gz \
    > 5_coassembly/A/all_R1.fastq.gz

metawrap assembly \
    -1 5_coassembly/A/all_R1.fastq.gz \
    -2 5_coassembly/A/all_R2.fastq.gz \
    -m 400 -t 24 -l 300 \
    --megahit \
    -o 5_coassembly/A/megahit_coassembly
```

---

## Part 3. Gene prediction and gene catalog construction

### 3.1 Gene prediction (Prodigal)

```bash
for i in $(seq -w 1 24); do
    S=SAMPLE${i}
    prodigal \
        -i 4_assembly/${S}/final_assembly.fasta \
        -p meta \
        -o 6_prodigal/${S}_genes.gff \
        -a 6_prodigal/${S}_proteins.faa \
        -d 6_prodigal/${S}_genes.fna
done
```

```bash
#command run twice, once for facility A, once for facility B

prodigal \
    -i 5_coassembly/A/megahit_coassembly/final_assembly.fasta \
    -p meta \
    -o 6_prodigal/A_coassembly_genes.gff \
    -a 6_prodigal/A_coassembly_proteins.faa \
    -d 6_prodigal/A_coassembly_genes.fna
```

### 3.2 Two-stage clustering (CD-HIT)

Individual-assembly proteins are clustered first, then co-assembly proteins are added and the combined set is clustered again (mix-assembly approach). 
Both stages use 95% identity with 90% alignment coverage of the shorter sequence.

```bash
# Stage 1: cluster individual assembly proteins
cat 6_prodigal/SAMPLE*_proteins.faa > 7_gene_catalog/all_individual_proteins.faa

cd-hit \
    -i 7_gene_catalog/all_individual_proteins.faa \
    -o 7_gene_catalog/individual_nr_proteins.faa \
    -c 0.95 -aS 0.9 -n 5 \
    -T 32 -M 64000 -d 0

# Stage 2: add co-assembly proteins and re-cluster
cat 7_gene_catalog/individual_nr_proteins.faa \
    6_prodigal/A_coassembly_proteins.faa \
    6_prodigal/B_coassembly_proteins.faa \
    > 7_gene_catalog/nr_individual_plus_coassemblies.faa

cd-hit \
    -i 7_gene_catalog/nr_individual_plus_coassemblies.faa \
    -o 7_gene_catalog/mix_catalog_proteins.faa \
    -c 0.95 -aS 0.9 -n 5 \
    -T 32 -M 64000 -d 0
```

#### Notes

- MEGAHIT reuses internal contig IDs (e.g., `k141_XXXXX`) across independent assemblies; Prodigal inherits these IDs, producing duplicate protein IDs for unrelated sequences.
- Duplicate IDs in the CD-HIT protein catalog were confirmed to represent distinct sequences and were renamed with `_dup1` suffixes rather than removed.
- Original → renamed IDs were saved in a lookup table for nucleotide extraction.
- Protein and nucleotide sequences were matched by exact ID using BioPython, with renamed IDs applied to corresponding nucleotide sequences.
- BioPython set-difference checks confirmed identical protein/nucleotide ID sets and sequence counts, with no sequences lost.

---

## Part 4. Read mapping and gene quantification

### 4.1 Build catalog index and map reads

```bash
bowtie2-build \
    7_gene_catalog/mix_catalog_genes_dedup.fna \
    7_gene_catalog/mix_catalog_index \
    --threads 32

for i in $(seq -w 1 24); do
    S=SAMPLE${i}
    echo "Mapping ${S} at $(date)..."

    bowtie2 \
        -x 7_gene_catalog/mix_catalog_index \
        -1 3_no_mouse_reads/${S}_no_mouse_1.fastq.gz \
        -2 3_no_mouse_reads/${S}_no_mouse_2.fastq.gz \
        --very-sensitive \
        --no-unal \
        -p 24 \
    | samtools view -bS - \
    | samtools sort -@ 8 -o 8_mapped_bam/${S}.bam

    samtools index 8_mapped_bam/${S}.bam
done
```

### 4.2 Count with CoverM

```bash
coverm contig \
    --bam-files 8_mapped_bam/*.bam \
    --methods count tpm length \
    --contig-end-exclusion 0 \
    --proper-pairs-only \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    -t 32 \
    -o 9_coverm/all_samples_coverm.tsv
```

---

## Part 5. Functional annotation of the gene catalog

### 5.1 Arsenic metabolism genes (AsgeneDB + DIAMOND)

```bash

wget https://data.cyverse.org/dav-anon/iplant/home/xinwei/AsgeneDB/AsgeneDB.zip
unzip AsgeneDB.zip

diamond makedb \
    --in AsgeneDB/AsgeneDB.fa \
    --db AsgeneDB


diamond blastp \
    --db db/AsgeneDB.dmnd \
    --query 7_gene_catalog/mix_catalog_proteins_dedup.faa \
    --out 10_annotations/AsDB/mix_catalog_AsDB.diamond.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
    --evalue 1e-10 \
    --max-target-seqs 10 \
    --threads 24 \
    --query-cover 80 \
    --subject-cover 50
```

### 5.2 CAZymes (run_dbcan)

```bash
run_dbcan CAZyme_annotation \
    --input_raw_data 7_gene_catalog/mix_catalog_proteins_dedup.faa \
    --output_dir 10_annotations/CAZyme \
    --db_dir db \
    --mode protein
```

### 5.3 COG/KEGG annotation (eggNOG-mapper)

```bash
#eggNOG database emapperdb-5.0.2 (eggnog.db, eggnog_proteins.dmnd, eggnog.taxa) downloaded to $EGGNOG_DATA_DIR prior to running

export EGGNOG_DATA_DIR=$(pwd)/db/eggnog-mapper-data

emapper.py \
    -i 7_gene_catalog/mix_catalog_proteins_dedup.faa \
    --itype proteins \
    -o 10_annotations/eggNOG/mix_catalog_eggnog \
    --data_dir $EGGNOG_DATA_DIR \
    --cpu 32 \
    --dbmem \
    -m diamond
```

---

## Part 6. Binning and MAG recovery

Binning was run separately on each facility co-assembly, using only the reads from that facility.

### 6.1 Stage reads

#reads were decompressed from `3_no_mouse_reads/` into `11_initial_binning/<facility>/reads/` with names ending in `_1.fastq` / `_2.fastq`, as required by the metawrap binning module

### 6.2 Initial binning

```bash
#command run twice, once for facility A, once for facility B

metawrap binning \
    -o 11_initial_binning/A \
    -t 32 \
    -a 5_coassembly/A/megahit_coassembly/final_assembly.fasta \
    --metabat2 --maxbin2 --concoct \
    11_initial_binning/A/reads/*_1.fastq \
    11_initial_binning/A/reads/*_2.fastq
```

### 6.3 Bin refinement

Bins retained at ≥90% completeness and ≤5% contamination.

```bash
#command run twice, once for facility A, once for facility B

metawrap bin_refinement \
    -o 12_bin_refinement/A \
    -t 32 \
    -A 11_initial_binning/A/metabat2_bins/ \
    -B 11_initial_binning/A/maxbin2_bins/ \
    -C 11_initial_binning/A/concoct_bins/ \
    -c 90 -x 5
```

### 6.4 Bin quantification

```bash
#command run twice, once for facility A, once for facility B

metawrap quant_bins \
    -b 12_bin_refinement/A/metawrap_90_5_bins/ \
    -o 13_quant_bins/A/ \
    -a 5_coassembly/A/megahit_coassembly/final_assembly.fasta \
    -t 32 \
    11_initial_binning/A/reads/*_1.fastq \
    11_initial_binning/A/reads/*_2.fastq
```

### 6.5 Bin classification

```bash
#command run twice, once for facility A, once for facility B

metawrap classify_bins \
    -b 12_bin_refinement/A/metawrap_90_5_bins/ \
    -o 15_bin_classification/A \
    -t 48
```

### 6.6 MAG gene calling and annotation

#Prokka (called by metawrap annotate_bins) rejects long contig names, so contigs in each refined bin were renamed to `<bin>_<n>` with awk before annotation

```bash
#command run twice, once for facility A, once for facility B

metaWRAP annotate_bins \
    -o 16_MAG_function/A \
    -t 96 \
    -b 12_bin_refinement/A/metawrap_90_5_bins_renamed/
```

---

## Part 7. Per-MAG functional annotation

### 7.1 Arsenic-cycling genes in MAGs


```bash
#command run twice, once for facility A, once for facility B

for faa in 16_MAG_function/A/bin_translated_genes/*.faa; do
    bin_name=$(basename "$faa" .faa)
    diamond blastp \
        --db db/AsgeneDB.dmnd \
        --query "$faa" \
        --out "17_AsDB_bins/A/${bin_name}.diamond.tsv" \
        --outfmt 6 \
        --evalue 1e-10 \
        --max-target-seqs 1 \
        --threads 24 \
        --query-cover 80 \
        --subject-cover 80
done
```

