# Mpox Sequence analysis from Illumina raw data (fastq files)
## The analysis pipeline follows raw sequence reads analysis from fastq to consensus fasta for amplicon based sequencing.

Hostile or bowtie or Kraken2 or BBMap (BBduk/BBsplit): Human read removal

FastQC: Quality check of raw reads

 Fastp or Trimmomatic: Adapter trimming and quality filtering

BWA mem or Minimap2: Mapping reads to references

Samtools: BAM processing and stats

Qualimap or Mosdepth: Mapping quality statistics

Freebayes or snippy for calling variants

bcftools for consensus building

## Step 1: Modules needed

Load modules if in hpc

```
module load hostile/2.0.0
module load fastqc/0.11.9
module load fastp/0.24.1
module load seqtk/1.3
module load bwa/0.7.17
module load freebayes/1.3.4
module unload samtools/1.17
module unload bcftools/1.13
```
Install modules if in the local terminal

```
conda install bioconda::fastp
conda install bioconda::fastqc
conda install bioconda::ivar
conda install bioconda::samtools
conda install bioconda::bfctools
conda install bioconda/label/broken::bcftools
conda install bioconda::minimap2
conda install bioconda::hostile
conda install bioconda::seqtk
conda install bioconda::freebayes
conda install bioconda::bwa
conda install bioconda::bwa-mem2
conda install bioconda::snakemake
conda install bioconda::mamba
conda install -n base mamba
conda create -c bioconda -c conda-forge -n squirrel -y squirrel
```
For geeks

```
#!/bin/bash

# List of required packages with their channels
declare -A packages=(
  [fastp]="bioconda"
  [fastqc]="bioconda"
  [ivar]="bioconda"
  [samtools]="bioconda"
  [samtools]="bioconda"
  [bcftools]="bioconda/label/broken"
  [minimap2]="bioconda"
  [hostile]="bioconda"
  [seqtk]="bioconda"
  [freebayes]="bioconda"
  [bwa]="bioconda"
  [bwa-mem2]="bioconda"
  [snakemake]="bioconda"
  [mamba]="bioconda"
)

# Function to check and install packages
for pkg in "${!packages[@]}"; do
  if conda list "$pkg" | grep -q "^$pkg"; then
    echo "$pkg is already installed."
  else
    echo "Installing $pkg from ${packages[$pkg]}..."
    conda install -y -c "${packages[$pkg]}" "$pkg"
  fi
done

# Install mamba in base environment if not already
if ! conda list -n base | grep -q "^mamba"; then
  echo "Installing mamba in base environment..."
  conda install -n base -y mamba
else
  echo "mamba is already installed in base."
fi

# Create squirrel environment if not exists
if ! conda info --envs | grep -q "^squirrel"; then
  echo "Creating 'squirrel' environment..."
  conda create -y -n squirrel -c bioconda -c conda-forge squirrel
else
  echo "'squirrel' environment already exists."
fi
```
Activate installed modules/make them available outisde mynev envt/run inside myenv or outside

```
echo 'export PATH=/home/woguta/anaconda3/envs/myenv/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```
Make eg squirrel executable in base

```
chmod +x /home/woguta/anaconda3/envs/squirrel/bin/squirrel
```
Confirm

```
head -n 1 /home/woguta/anaconda3/envs/squirrel/bin/squirrel
```
## Step 2: Define and create directories

Set directory paths

```
WORK_DIR="./mpox_files/mpox_sierra"
FASTQ_DIR="${WORK_DIR}/fastq_files"
FASTQC_DIR="${WORK_DIR}/fastqc_files"
FASTP_DIR="${WORK_DIR}/fastp_trimmed"
REF_DIR="${WORK_DIR}/refseqs"
DATABASE_DIR="${WORK_DIR}/databases"
OUT_DIR="${WORK_DIR}/results"
HOST_FILTERED_DIR="${WORK_DIR}/host_filtered"
BAM_DIR="${WORK_DIR}/mapped_bam"
VCF_DIR="${WORK_DIR}/vcf"
FASTA_DIR="${WORK_DIR}/fasta_files"
MPOX_REF1="${REF_DIR}/Mpox_ref_NC_063383.1.fasta"
MPOX_REF2="${REF_DIR}/mpox_ref_NC_003310.1.fasta"
```
Create directories

```
mkdir -p "$FASTQC_DIR" "$FASTA_DIR" "$REF_DIR" "$DATABASE_DIR" "$OUT_DIR" "$FASTP_DIR" "$TRIMMED_DIR" "$HOST_FILTERED_DIR" "$BAM_DIR" "$VCF_DIR" "$CONSENSUS_DIR" "$TMP_DIR"
```
Group the directories into array

```
DIRS=(
  "$FASTQC_DIR" "$FASTP_DIR" "$REF_DIR" "$DATABASE_DIR" "$OUT_DIR" "$TRIMMED_DIR"
  "$HOST_FILTERED_DIR" "$BAM_DIR" "$VCF_DIR" "$FASTA_DIR" "$TMP_DIR"
)

# Create directories if they don't exist
for dir in "${DIRS[@]}"; do
  [ -d "$dir" ] || mkdir -p "$dir"
done
```
## Step 3: Clean out human reads using default human-t2t-hla genome

```
hostile clean \
    --fastq1 "$FASTQ_DIR/515_S13_L001_R1_001.fastq.gz" \
    --fastq2 "$FASTQ_DIR/515_S13_L001_R2_001.fastq.gz" \
    --out-dir "$HOST_FILTERED_DIR" \
    --force 
```
## Step 4: First Quality control

```
fastqc \
    -f fastq "$HOST_FILTERED_DIR/515_S13_L001_R1_001.clean_1.fastq.gz" \
            "$HOST_FILTERED_DIR/515_S13_L001_R2_001.clean_2.fastq.gz" \
    -o "$FASTQC_DIR"
```
## Step 5: Adapter and low quality reads trimming

Using fastp simple/all default

```
fastp \
    -i "$HOST_FILTERED_DIR/515_S13_L001_R1_001.clean_1.fastq.gz" \
    -I "$HOST_FILTERED_DIR/515_S13_L001_R2_001.clean_2.fastq.gz" \
    -o "$FASTP_DIR/515_S13_trim_R1.fastq.gz" \
    -O "$FASTP_DIR/515_S13_trim_R2.fastq.gz"
```
Being stringent: preferred!

```
fastp \
  --in1 "$HOST_FILTERED_DIR/515_S13_L001_R1_001.clean_1.fastq.gz" \
  --in2 "$HOST_FILTERED_DIR/515_S13_L001_R2_001.clean_2.fastq.gz" \
  --out1 "$FASTP_DIR/515_S13_trim_R1.fastq.gz" \
  --out2 "$FASTP_DIR/515_S13_trim_R2.fastq.gz" \
  --detect_adapter_for_pe \
  --json "$FASTP_DIR/515_S13.fastp.json" \
  --html "$FASTP_DIR/515_S13.fastp.html" \
  --cut_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --length_required 20 \
  2> "$FASTP_DIR/515_S13.fastp.log"
```
## Step 6: Second Quality control/post trim fastQC

```
fastqc "$FASTP_DIR/515_S13_trim_R1.fastq.gz" "$FASTP_DIR/515_S13_trim_R2.fastq.gz" -o "$FASTQC_DIR"
```
## Step 7: Index the reference genome using bwa

```
mkdir -p "$REF_DIR/index" # create index directory
cp -rf "$REF_DIR/Mpox_ref_NC_063383.1.fasta" "$REF_DIR/index" #Copy recursively ref to index directory
bwa index -p "$REF_DIR/index/Mpox_ref_NC_063383.1" \
        "$REF_DIR/index/Mpox_ref_NC_063383.1.fasta"
```
Index the ref in refseqs directory using samtools

```
samtools faidx "$MPOX_REF1"
```
## Step 8: Map to MPXV reference using bwa and sort using samtools

```
bwa mem \
    "$REF_DIR/index/Mpox_ref_NC_063383.1" \
    "$FASTP_DIR/515_S13_trim_R1.fastq.gz" \
    "$FASTP_DIR/515_S13_trim_R1.fastq.gz" | \
    samtools sort -o "$BAM_DIR/515_S13.sorted.bam" 
```
Index the sorted bam alignment

```
samtools index -f "$BAM_DIR/515_S13.sorted.bam"
```
## Step 9: Variant calling using freebayes

Save into genomic variant call format for concensus building

```
freebayes \
    -p 1 \
    -f "$MPOX_REF1" \
    -F 0.2 \
    -C 1 \
    --pooled-continuous \
    --min-coverage 10 \
    --gvcf \
    --gvcf-dont-use-chunk true \
    "$BAM_DIR/515_S13.sorted.bam" > "$VCF_DIR/515_S13.gvcf"
```
# Save into vcf file for variant further studies

``
freebayes \
    -p 1 \
    -f "$MPOX_REF1" \
    -F 0.2 \
    -C 1 \
    --pooled-continuous \
    --min-coverage 10 \
    --vcf \
    --variant-input \
    "$BAM_DIR/515_S13.sorted.bam" > "$VCF_DIR/515_S13.vcf"
```
