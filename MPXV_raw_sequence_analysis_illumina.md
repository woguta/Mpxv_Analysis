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
conda install bioconda::bcftools
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
FASTQ_DIR="${WORK_DIR}/fastq"
FASTQC_DIR="${WORK_DIR}/fastqc"
FASTP_DIR="${WORK_DIR}/fastp"
REF_DIR="${WORK_DIR}/refseqs"
DATABASE_DIR="${WORK_DIR}/databases"
OUT_DIR="${WORK_DIR}/results"
HOST_FILTERED_DIR="${WORK_DIR}/hostile"
BAM_DIR="${WORK_DIR}/bam"
VCF_DIR="${WORK_DIR}/vcf"
FASTA_DIR="${WORK_DIR}/fasta"
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

Save into genomic variant call format for genomic concensus fasta building

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
Save into vcf file for variant further studies or use the *.variants.vcf file

```
freebayes \
    -p 1 \
    -f "$MPOX_REF1" \
    -F 0.2 \
    -C 1 \
    --pooled-continuous \
    --min-coverage 10 \
    "$BAM_DIR/515_S13.sorted.bam" > "$VCF_DIR/515_S13.vcf.gz"
```

## Step 10: Create consensus variants, low-frequency variants and a coverage mask

Create a function to process the gvcf and save in your scripts directory as "process_gvcf.py"

```
#!/usr/bin/env python3

import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Process gVCF for mask and variant filtering")
    parser.add_argument("gvcf", help="Input gVCF (can be .gz)")
    parser.add_argument("-d", "--min_depth", type=int, default=10, help="Minimum depth to keep")
    parser.add_argument("-l", "--min_af", type=float, default=0.25, help="Minimum allele frequency")
    parser.add_argument("-u", "--max_af", type=float, default=0.75, help="Maximum allele frequency")
    parser.add_argument("-m", "--mask_file", required=True, help="Output BED file with masked regions")
    parser.add_argument("-v", "--variant_file", required=True, help="Filtered variants VCF output")
    parser.add_argument("-c", "--consensus_file", required=True, help="Filtered consensus VCF output")
    return parser.parse_args()

def open_file(filename):
    return gzip.open(filename, 'rt') if filename.endswith(".gz") else open(filename, 'r')

def process_gvcf(args):
    with open_file(args.gvcf) as infile, \
         open(args.mask_file, 'w') as mask_out, \
         open(args.variant_file, 'w') as var_out, \
         open(args.consensus_file, 'w') as cons_out:

        for line in infile:
            if line.startswith("#"):
                var_out.write(line)
                cons_out.write(line)
                continue

            fields = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info, format_, sample = fields

            info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in info.split(";") if "=" in kv}

            depth = int(info_dict.get("DP", 0))
            af = float(info_dict.get("AF", 0.0)) if "AF" in info_dict else None
            end = int(info_dict.get("END", pos))  # gVCF block END

            # mask regions with low coverage
            if depth < args.min_depth:
                mask_out.write(f"{chrom}\t{int(pos)-1}\t{end}\n")
                continue

            if alt == "." or alt == "<NON_REF>":
                # likely reference or uninformative block, write to consensus only
                cons_out.write(line)
                continue

            if af is not None and (af < args.min_af or af > args.max_af):
                # ambiguous frequency, mask
                mask_out.write(f"{chrom}\t{int(pos)-1}\t{end}\n")
                continue

            # Otherwise, write variant
            var_out.write(line)
            cons_out.write(line)

if __name__ == "__main__":
    args = parse_args()
    process_gvcf(args)
```

Process the gvcf file

```
python ./myscripts/process_gvcf.py \
  -d 10 \
  -l 0.25 \
  -u 0.75 \
  -m "$VCF_DIR/515_S13.mask.txt" \
  -v "$VCF_DIR/515_S13.variants.vcf" \
  -c "$VCF_DIR/515_S13.consensus.vcf" \
  "$VCF_DIR/515_S13.gvcf.gz"
```
Compress and index the processd and compressed gVCF file

```
bgzip -c "$VCF_DIR/515_S13.variants.vcf" > "$VCF_DIR/515_S13.variants.vcf.gz"
bcftools index -f "$VCF_DIR/515_S13.variants.vcf.gz"
bcftools index -f "$VCF_DIR/515_S13.gvcf.gz"
bcftools index -f "$VCF_DIR/515_S13.vcf.gz"
```

## Step 11: Normalize variant records into canonical VCF representation

```
for v in "variants" "consensus"; do
    echo -e "normalising variants in: $v"
    bcftools norm \
        -f "$MPOX_REF1" \
        "$VCF_DIR/515_S13.$v.vcf" > "$VCF_DIR/515_S13.$v.norm.vcf"
done
```

## Step 12: Split consensus VCF file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF

```
for vt in "ambiguous" "fixed"; do
    echo "Splitting on ConsensusTag: $vt"
    awk -v vartag="ConsensusTag=$vt" \
        '$0 ~ /^#/ || $0 ~ vartag' \
        "$VCF_DIR/515_S13.consensus.norm.vcf" > "$VCF_DIR/515_S13.$vt.norm.vcf"

    bgzip -f "$VCF_DIR/515_S13.$vt.norm.vcf"
    tabix -f -p vcf "$VCF_DIR/515_S13.$vt.norm.vcf.gz"
done
```
## Step 13: Apply ambiguous variants first using IUPAC codes, has no indels

```
bcftools consensus \
    -f "$MPOX_REF1" \
    -I "$VCF_DIR/515_S13.ambiguous.norm.vcf.gz" > "$VCF_DIR/515_S13.ambiguous.fa"
```

Get viral contig name from reference

```
CTG_NAME=$(head -n1 "$MPOX_REF1" | sed 's/>//')
```

Make sure bed file is in correct coordinate formats

```
awk '{if ($2 == 0) $2 = 1; else $2=$2+1}1' OFS="\t" "$VCF_DIR/515_S13.mask.txt" > "$VCF_DIR/515_S13.mask.1based.txt"
```

## Step 14: Build consensus variants, low-frequency variants and a coverage mask generation into final genomic sequence in fasta using bcftools

```
bcftools consensus \
    -f "$VCF_DIR/515_S13.ambiguous.fa" \
    -m "$VCF_DIR/515_S13.mask.1based.txt" \
    "$VCF_DIR/515_S13.fixed.norm.vcf.gz" | \
    sed "s|$CTG_NAME|515_S13|" > "$FASTA_DIR/515_S13.fasta"
```
## Step 15: Create a loop/script to process all the fasta files
A. Save as mpox_fastq2fasta_sierra.sh, stringent script

```
#!/bin/bash

# Define paths
FASTQ_DIR="./mpox_files/mpox_sierra/fastq"
HOST_FILTERED_DIR="./mpox_files/mpox_sierra/hostile"
FASTQC_DIR="./mpox_files/mpox_sierra/fastqc"
FASTP_DIR="./mpox_files/mpox_sierra/fastp"
REF_DIR="./mpox_files/mpox_sierra/refseqs"
VCF_DIR="./mpox_files/mpox_sierra/vcf"
BAM_DIR="./mpox_files/mpox_sierra/bam"
FASTA_DIR="./mpox_files/mpox_sierra/fasta"
MPOX_REF1="$REF_DIR/Mpox_ref_NC_063383.1.fasta"

# Index the reference (only once)
mkdir -p "$REF_DIR/index"
cp -f "$MPOX_REF1" "$REF_DIR/index"
bwa index -p "$REF_DIR/index/Mpox_ref_NC_063383.1" "$REF_DIR/index/Mpox_ref_NC_063383.1.fasta"
samtools faidx "$MPOX_REF1"

# Loop over FASTQ R1 files
for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" | cut -d'_' -f1,2)
    R2="${R1/_R1_/_R2_}"
    echo "Processing sample: $SAMPLE"

    # Step 1: Remove host reads
    echo "Removing human reads for for $SAMPLE..."
    hostile clean \
        --fastq1 "$R1" \
        --fastq2 "$R2" \
        --out-dir "$HOST_FILTERED_DIR" \
        --force

    CLEAN_R1="$HOST_FILTERED_DIR/${SAMPLE}_L001_R1_001.clean_1.fastq.gz"
    CLEAN_R2="$HOST_FILTERED_DIR/${SAMPLE}_L001_R2_001.clean_2.fastq.gz"

    # Step 2: FastQC
    echo "Running FASTQC for $SAMPLE..."
    fastqc -f fastq "$CLEAN_R1" "$CLEAN_R2" -o "$FASTQC_DIR"

    # Step 3: Trim adapters and low-quality reads
    echo "Trimming adapters & low reads for $SAMPLE..."
    fastp \
        --in1 "$CLEAN_R1" \
        --in2 "$CLEAN_R2" \
        --out1 "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        --out2 "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --json "$FASTP_DIR/${SAMPLE}.fastp.json" \
        --html "$FASTP_DIR/${SAMPLE}.fastp.html" \
        --cut_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 20 \
        2> "$FASTP_DIR/${SAMPLE}.fastp.log"

    # Step 4: Post-trim FastQC
    echo "Post-trim FASTQC for $SAMPLE..."
    fastqc "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" -o "$FASTQC_DIR"

    # Step 5: Map to reference
    echo "Mapping to ref for $SAMPLE..."
    bwa mem "$REF_DIR/index/Mpox_ref_NC_063383.1" \
        "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" | \
        samtools sort -o "$BAM_DIR/${SAMPLE}.sorted.bam"

    samtools index -f "$BAM_DIR/${SAMPLE}.sorted.bam"

    # Step 6: Variant calling
    echo "Calling variants for $SAMPLE..."
    freebayes \
        -p 1 \
        -f "$MPOX_REF1" \
        -F 0.2 \
        -C 1 \
        --pooled-continuous \
        --min-coverage 10 \
        --gvcf \
        --gvcf-dont-use-chunk true \
        "$BAM_DIR/${SAMPLE}.sorted.bam" > "$VCF_DIR/${SAMPLE}.gvcf"

    bgzip -f "$VCF_DIR/${SAMPLE}.gvcf"
    bcftools index -f "$VCF_DIR/${SAMPLE}.gvcf.gz"

    # Step 7: Process gVCF
    echo "Processing gVCF for $SAMPLE..."
    python ./myscripts/process_gvcf.py \
        -d 10 \
        -l 0.25 \
        -u 0.75 \
        -m "$VCF_DIR/${SAMPLE}.mask.txt" \
        -v "$VCF_DIR/${SAMPLE}.variants.vcf" \
        -c "$VCF_DIR/${SAMPLE}.consensus.vcf" \
        "$VCF_DIR/${SAMPLE}.gvcf.gz"

    # Step 8: Normalize VCFs
    echo "Normalizing VCF for $SAMPLE..."
    for v in "variants" "consensus"; do
        bcftools norm \
            -f "$MPOX_REF1" \
            "$VCF_DIR/${SAMPLE}.${v}.vcf" > "$VCF_DIR/${SAMPLE}.${v}.norm.vcf"
    done

    # Step 9: Split ambiguous vs fixed
    echo "Splitting ambiguous vs fixed for $SAMPLE..."
    for vt in "ambiguous" "fixed"; do
        awk -v vartag="ConsensusTag=$vt" \
            '$0 ~ /^#/ || $0 ~ vartag' \
            "$VCF_DIR/${SAMPLE}.consensus.norm.vcf" > "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf"
        bgzip -f "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf"
        tabix -f -p vcf "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf.gz"
    done

    # Step 10: Apply ambiguous variants using IUPAC codes
    echo "Applying IUPAC codes for $SAMPLE..."
    bcftools consensus \
        -f "$MPOX_REF1" \
        -I "$VCF_DIR/${SAMPLE}.ambiguous.norm.vcf.gz" > "$VCF_DIR/${SAMPLE}.ambiguous.fa"

    # Step 11: Get contig name from FASTA
    CTG_NAME=$(head -n1 "$MPOX_REF1" | sed 's/>//')

    # Step 12: Correct mask file to 1-based
    echo "Correcting mask to 1-based for $SAMPLE..."
    awk '{if ($2 == 0) $2 = 1; else $2 = $2 + 1}1' OFS="\t" \
        "$VCF_DIR/${SAMPLE}.mask.txt" > "$VCF_DIR/${SAMPLE}.mask.1based.txt"

    # Step 13: Final consensus
    echo "Final concesus fasta for $SAMPLE..."
    bcftools consensus \
        -f "$VCF_DIR/${SAMPLE}.ambiguous.fa" \
        -m "$VCF_DIR/${SAMPLE}.mask.1based.txt" \
        "$VCF_DIR/${SAMPLE}.fixed.norm.vcf.gz" | \
        sed "s|$CTG_NAME|${SAMPLE}|" > "$FASTA_DIR/${SAMPLE}.fasta"
done
```
B. Save as mpox_fastq2fasta1_sierra.sh, dafault script  but with  gvcf
```
#!/bin/bash

# Define paths
FASTQ_DIR="./mpox_files/mpox_sierra/fastq"
FASTQC_DIR="./mpox_files/mpox_sierra/fastqc"
FASTP_DIR="./mpox_files/mpox_sierra/fastp"
REF_DIR="./mpox_files/mpox_sierra/refseqs"
VCF_DIR="./mpox_files/mpox_sierra/vcf"
BAM_DIR="./mpox_files/mpox_sierra/bam"
FASTA_DIR="./mpox_files/mpox_sierra/fasta"
MPOX_REF1="$REF_DIR/Mpox_ref_NC_063383.1.fasta"

# Index the reference (only once)
mkdir -p "$REF_DIR/index"
cp -f "$MPOX_REF1" "$REF_DIR/index"
bwa index -p "$REF_DIR/index/Mpox_ref_NC_063383.1" "$REF_DIR/index/Mpox_ref_NC_063383.1.fasta"
samtools faidx "$MPOX_REF1"

# Loop over FASTQ R1 files
for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" | cut -d'_' -f1,2)
    R2="${R1/_R1_/_R2_}"
    echo "Processing sample: $SAMPLE"

    # Step 1: Initial FASTQC
    fastqc -f fastq "$R1" "$R2" -o "$FASTQC_DIR"

    # Step 2: Adapter/quality trimming
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        -O "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --json "$FASTP_DIR/${SAMPLE}.fastp.json" \
        --html "$FASTP_DIR/${SAMPLE}.fastp.html" \
        2> "$FASTP_DIR/${SAMPLE}.fastp.log"

    # Step 3: Post-trim FastQC
    fastqc "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" -o "$FASTQC_DIR"

    # Step 4: Mapping to reference
    bwa mem "$REF_DIR/index/Mpox_ref_NC_063383.1" \
        "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" | \
        samtools sort -o "$BAM_DIR/${SAMPLE}.sorted.bam"

    samtools index -f "$BAM_DIR/${SAMPLE}.sorted.bam"

    # Step 5: Variant calling
    freebayes \
        -p 1 \
        -f "$MPOX_REF1" \
        "$BAM_DIR/${SAMPLE}.sorted.bam" > "$VCF_DIR/${SAMPLE}.gvcf"

    bgzip -f "$VCF_DIR/${SAMPLE}.gvcf"
    bcftools index -f "$VCF_DIR/${SAMPLE}.gvcf.gz"

    # Step 6: Process gVCF
    python ./myscripts/process_gvcf.py \
        -m "$VCF_DIR/${SAMPLE}.mask.txt" \
        -v "$VCF_DIR/${SAMPLE}.variants.vcf" \
        -c "$VCF_DIR/${SAMPLE}.consensus.vcf" \
        "$VCF_DIR/${SAMPLE}.gvcf.gz"

    # Step 7: Normalize VCFs
    for v in "variants" "consensus"; do
        bcftools norm \
            -f "$MPOX_REF1" \
            "$VCF_DIR/${SAMPLE}.${v}.vcf" > "$VCF_DIR/${SAMPLE}.${v}.norm.vcf"
    done

    # Step 8: Split ambiguous vs fixed
    for vt in "ambiguous" "fixed"; do
        awk -v vartag="ConsensusTag=$vt" \
            '$0 ~ /^#/ || $0 ~ vartag' \
            "$VCF_DIR/${SAMPLE}.consensus.norm.vcf" > "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf"
        bgzip -f "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf"
        tabix -f -p vcf "$VCF_DIR/${SAMPLE}.${vt}.norm.vcf.gz"
    done

    # Step 9: Apply ambiguous variants
    bcftools consensus \
        -f "$MPOX_REF1" \
        -I "$VCF_DIR/${SAMPLE}.ambiguous.norm.vcf.gz" > "$VCF_DIR/${SAMPLE}.ambiguous.fa"

    # Step 10: Contig name
    CTG_NAME=$(head -n1 "$MPOX_REF1" | sed 's/>//')

    # Step 11: Convert BED mask to 1-based
    awk '{if ($2 == 0) $2 = 1; else $2 = $2 + 1}1' OFS="\t" \
        "$VCF_DIR/${SAMPLE}.mask.txt" > "$VCF_DIR/${SAMPLE}.mask.1based.txt"

    # Step 12: Final consensus FASTA
    bcftools consensus \
        -f "$VCF_DIR/${SAMPLE}.ambiguous.fa" \
        -m "$VCF_DIR/${SAMPLE}.mask.1based.txt" \
        "$VCF_DIR/${SAMPLE}.fixed.norm.vcf.gz" | \
        sed "s|$CTG_NAME|${SAMPLE}|" > "$FASTA_DIR/${SAMPLE}.fasta"
done
```

C. Save the script as mpox_fastq2fasta2_sierra.sh. Use the default pipeline, but prioritize using VCF-based consensus generation. Since FastQC results showed good sequence coverage and the data comes from targeted amplicon-based sequencing, there's no need for host read removal. To capture low-frequency variants, VCF output is preferred—using overly stringent filtering may miss these important mutations.

```
#!/bin/bash

# Define paths
FASTQ_DIR="./mpox_files/mpox_sierra/fastq"
FASTQC_DIR="./mpox_files/mpox_sierra/fastqc"
FASTP_DIR="./mpox_files/mpox_sierra/fastp"
REF_DIR="./mpox_files/mpox_sierra/refseqs"
VCF_DIR="./mpox_files/mpox_sierra/vcf"
BAM_DIR="./mpox_files/mpox_sierra/bam"
FASTA_DIR="./mpox_files/mpox_sierra/fasta"
MPOX_REF1="$REF_DIR/Mpox_ref_NC_063383.1.fasta"

# Index the reference (only once)
mkdir -p "$REF_DIR/index"
cp -f "$MPOX_REF1" "$REF_DIR/index"
bwa index -p "$REF_DIR/index/Mpox_ref_NC_063383.1" "$REF_DIR/index/Mpox_ref_NC_063383.1.fasta"
samtools faidx "$MPOX_REF1"

# Loop over FASTQ R1 files
for R1 in "$FASTQ_DIR"/*_R1_001.fastq.gz; do
    SAMPLE=$(basename "$R1" | cut -d'_' -f1,2)
    R2="${R1/_R1_/_R2_}"
    FASTA_FILE="$FASTA_DIR/${SAMPLE}.fa"

    # Skip sample if FASTA already exists
    if [[ -f "$FASTA_FILE" ]]; then
        echo -e "Fasta already exists for $SAMPLE at $FASTA_FILE - Skipping.\n"
    continue
    else
        echo -e "Proceeding processing $SAMPLE"
    fi

    # Step 1: Initial FASTQC
    fastqc -f fastq "$R1" "$R2" -o "$FASTQC_DIR"

    # Step 2: Adapter/quality trimming
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        -O "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --json "$FASTP_DIR/${SAMPLE}.fastp.json" \
        --html "$FASTP_DIR/${SAMPLE}.fastp.html" \
        2> "$FASTP_DIR/${SAMPLE}.fastp.log"

    # Step 3: Post-trim FastQC
    fastqc "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" -o "$FASTQC_DIR"

    # Step 4: Mapping to reference
    bwa mem "$REF_DIR/index/Mpox_ref_NC_063383.1" \
        "$FASTP_DIR/${SAMPLE}_trim_R1.fastq.gz" \
        "$FASTP_DIR/${SAMPLE}_trim_R2.fastq.gz" | \
        samtools sort -o "$BAM_DIR/${SAMPLE}.sorted.bam"

    samtools index -f "$BAM_DIR/${SAMPLE}.sorted.bam"

    # Step 5: Variant calling
    freebayes \
        -p 1 \
        -f "$MPOX_REF1" \
        "$BAM_DIR/${SAMPLE}.sorted.bam" > "$VCF_DIR/${SAMPLE}.vcf"

    bgzip -f "$VCF_DIR/${SAMPLE}.vcf"
    bcftools index -f "$VCF_DIR/${SAMPLE}.vcf.gz"

    # Step 5: Contig name
    CTG_NAME=$(head -n1 "$MPOX_REF1" | sed 's/>//')

    # Step 6: Final consensus FASTA
    bcftools consensus \
        -f "$MPOX_REF1" \
        -I "$VCF_DIR/${SAMPLE}.vcf.gz" | \
       sed "s|$CTG_NAME|${SAMPLE}|" > "$FASTA_DIR/${SAMPLE}.fa"
done
```

Run for apobec3 signatures for sustained human-to-human transmission

```
squirrel \
#/home/woguta/anaconda3/envs/squirrel/bin/squirrel \
    "$OUT_DIR/algn/mpox_aligned-nuc_2025-06-13T0911.fasta" \
    --no-mask \
    --seq-qc \
    --outdir "$OUT_DIR/squirrel" \
    --outfile apobec3_test.fasta \
    --run-phylo \
    --run-apobec3-phylo \
    --interactive-tree \
    --clade cladeiib
    --clade cladeiib
```
