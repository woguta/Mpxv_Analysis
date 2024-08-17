# Sequence Annotations

## A. Prokka
1. Predict genes based on ref genome gff
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Emb_Prk_Mpox
#SBATCH -n 3
#SBATCH --mem=8G

# Exit with error reporting
set -e
set -x

# Load modules
module purge
module load clustalw/2.1
module load emboss/6.6.0
module load samtools/1.15.1
module load prokka/1.14.6

# Define directories
WORK_DIR="/var/scratch/${USER}/viral_genomes/mpox_data"
FASTA_DIR="${WORK_DIR}/mpox_fasta"
OUT_DIR="${WORK_DIR}/mpox_results/mpox_prokka"

# Make the output directory if it doesn't exist
if [ ! -e "${OUT_DIR}" ]; then
  mkdir -p "${OUT_DIR}"
fi

# Reference protein fasta in gb format
protein_db="${WORK_DIR}/mpox_refseqs/Mpox_ref_NC_063383.1.gb"

echo "Running Prokka"
# Loop through all fasta.gz files in the input directory
for gz_file in "${FASTA_DIR}"/*.fasta.gz; do
    if [ -f "$gz_file" ]; then
        sample=$(basename "$gz_file" .fasta.gz)
        temp_fasta="${FASTA_DIR}/${sample}.temp.fasta"
        final_fasta="${FASTA_DIR}/${sample}.fasta"
        compressed_fasta="${OUT_DIR}/${sample}.fasta.gz"
        output_path="${OUT_DIR}/${sample}"

        # Decompress the gz file to a temporary fasta file
        echo "Decompressing ${gz_file} to ${temp_fasta}"
        gunzip -c "$gz_file" > "$temp_fasta"

        # Check if the temporary file is in FASTA format
        if grep -q "^>" "$temp_fasta"; then
            echo "File ${temp_fasta} is in FASTA format"
            fasta_file="$temp_fasta"
        else
            # Convert to FASTA format if it's not
            fasta_file="${FASTA_DIR}/${sample}.converted.fasta"
            echo "Converting ${temp_fasta} to FASTA format"
            seqret -sequence "$temp_fasta" -outseq "$fasta_file"
        fi

        # Clean the FASTA headers and save to the final fasta file
        echo "Cleaning headers in ${fasta_file}"
        awk '/^>/ {header=$0; gsub(/[|\/\-]/, "_", header); print header} !/^>/ {print}' "$fasta_file" > "$final_fasta"

        # Ensure the final fasta file has a proper FASTA header
        if ! grep -q "^>" "$final_fasta"; then
        echo "ERROR!!! The file ${final_fasta} does not have a proper FASTA header. Adding default headers."

        # Add default headers to each sequence
        awk 'BEGIN {header_num=1} !/^>/ {print ">sequence_" header_num "\n" $0; header_num++} /^>/ {print}' "$final_fasta" > "${final_fasta}.tmp"
        mv "${final_fasta}.tmp" "$final_fasta"
        echo "Default headers added to ${final_fasta}"
        else
            echo "The file ${final_fasta} has been cleaned and is in proper FASTA format"

        fi

        # Compress the final fasta file
        echo "Compressing ${final_fasta} to ${compressed_fasta}"
        bgzip -c "$final_fasta" > "$compressed_fasta"

        # Index the compressed fasta file
        echo "Indexing ${compressed_fasta}"
        samtools faidx "$compressed_fasta"

        # Echo fasta file path for the current sample
        echo "contigs_fasta for ${sample}: ${final_fasta}"

        # Run the Prokka command on the current file
        echo "Processing Prokka for ${sample}: ${final_fasta}"
        prokka "$final_fasta" \
            --outdir "$output_path" \
            --cpus 3 \
            --mincontiglen 200 \
            --kingdom Viruses \
            --centre WHO \
            --addgenes \
            --addmrna \
            --locustag MPXV \
            --genus "Mpox virus" \
            --proteins "$protein_db" \
            --usegenus \
            --compliant \
            --rfam \
            --force \
            --debug

        # Clean up all temp files
        rm -f "$temp_fasta"

    else
        echo "ERROR!!! fasta.gz file not found: ${gz_file}"
    fi
done
```
2. Concatenate & group prokka data of interest (gffs, gbks, gbfs & faas) per sample in each samples own folder

```
#!/bin/bash

# Define directories
WORK_DIR="/var/scratch/${USER}/viral_genomes/mpox_data"
PROKKA_DIR="${WORK_DIR}/mpox_results/mpox_prokka"
OUTPUT_DIR="${WORK_DIR}/mpox_results/prokka_data"

# Make the output directory if it doesn't exist
if [ ! -e "${OUTPUT_DIR}" ]; then
  mkdir -p "${OUTPUT_DIR}"
fi

# Add a debug statement to check PLASMIDS_DIR
echo "PROKKA_DIR: ${PROKKA_DIR}"

# Iterate over subdirectories in GFF_DIR
for sample_dir in "${PROKKA_DIR}"/*; do
  if [ ! -d "${sample_dir}" ]; then
    continue
  fi

  # Extract the sample name from the directory name
  sample_name=$(basename "${sample_dir}") 
  gff_file="${sample_dir}/PROKKA_07132024.gff"
  gbk_file="${sample_dir}/PROKKA_07132024.gbk"
  gbf_file="${sample_dir}/PROKKA_07132024.gbf"
  faa_file="${sample_dir}/PROKKA_07132024.faa"

  #Check the  files of interests #gff #gbk #gbf #faa
  if [ ! -f "${gff_file}" ]; then
    echo "Error: ${gff_file} not found!"
  elif [ ! -s "${gff_file}" ]; then
    echo "Error: ${gff_file} is empty!"
  elif [ ! -r "${gff_file}" ]; then
    echo "Error: ${gff_file} is unreadable!"
  else
    echo "${gff_file} is present and readable."
  fi

  if [ ! -f "${gbk_file}" ]; then
    echo "Error: ${gbk_file} not found!"
  elif [ ! -s "${gbk_file}" ]; then
    echo "Error: ${gbk_file} is empty!"
  elif [ ! -r "${gbk_file}" ]; then
    echo "Error: ${gbk_file} is unreadable!"
  else
    echo "${gbk_file} is present and readable."
  fi

  if [ ! -f "${gbf_file}" ]; then
    echo "Error: ${gbf_file} not found!"
  elif [ ! -s "${gbf_file}" ]; then
    echo "Error: ${gbf_file} is empty!"
  elif [ ! -r "${gbf_file}" ]; then
    echo "Error: ${gbf_file} is unreadable!"
  else
    echo "${gbf_file} is present and readable."
  fi

  if [ ! -f "${faa_file}" ]; then
    echo "Error: ${faa_file} not found!"
  elif [ ! -s "${faa_file}" ]; then
    echo "Error: ${faa_file} is empty!"
  elif [ ! -r "${faa_file}" ]; then
    echo "Error: ${faa_file} is unreadable!"
  else
    echo "${faa_file} is present and readable."
  fi

  OUT="${OUTPUT_DIR}/${sample_name}"
  mkdir -p "${OUT}"

  cat_gff="${OUT}/${sample_name}.gff"
  cat_gbk="${OUT}/${sample_name}.gbk"
  cat_gbf="${OUT}/${sample_name}.gbf"
  cat_faa="${OUT}/${sample_name}.faa"

  echo "Concatenating files for ${sample_name}:"
  echo "Gff file: ${gff_file}"
  echo "Gbk file: ${gbk_file}"
  echo "Gbf file: ${gbf_file}"
  echo "Faa file: ${faa_file}"

  cat "${gff_file}" > "${cat_gff}"
  cat "${gbk_file}" > "${cat_gbk}"
  cat "${gbf_file}" > "${cat_gbf}"
  cat "${faa_file}" > "${cat_faa}"
  echo "Concatenating is complete for ${sample_name}"
  echo "Concatenated assembly: ${cat_gff}"
  echo "Concatenated assembly: ${cat_gbk}"
  echo "Concatenated assembly: ${cat_gbf}"
  echo "Concatenated assembly: ${cat_faa}"
done

echo "Concatenation for all samples is complete!"
```
3. Concatenate & group prokka data gffs, gbks, gbfs & faas in their own folders eg. groupind all gbks for all samples in gbk_data folder

```
#!/bin/bash

# Define directories
WORK_DIR="/var/scratch/woguta/viral_genomes/mpox_data"
PROKKA_DIR="${WORK_DIR}/mpox_results/prokka_data"

# Add a debug statement to check PPROKKA_DIR
echo "PROKKA_DIR: ${PROKKA_DIR}"

# Iterate over subdirectories in PROKKA_DIR
for sample_dir in "${PROKKA_DIR}"/*; do
  if [ ! -d "${sample_dir}" ]; then
    continue
  fi

  # Extract the sample name from the directory name
  sample_name=$(basename "${sample_dir}") 
  gff_file="${sample_dir}/${sample_name}.gff"
  gbk_file="${sample_dir}/${sample_name}.gbk"
  gbf_file="${sample_dir}/${sample_name}.gbf"
  faa_file="${sample_dir}/${sample_name}.faa"

  #Check the  files of interests #gff #gbk #gbf #faa
  if [ ! -f "${gff_file}" ]; then
    echo "Error: ${gff_file} not found!"
  elif [ ! -s "${gff_file}" ]; then
    echo "Error: ${gff_file} is empty!"
  elif [ ! -r "${gff_file}" ]; then
    echo "Error: ${gff_file} is unreadable!"
  else
    echo "${gff_file} is present and readable."
  fi

  if [ ! -f "${gbk_file}" ]; then
    echo "Error: ${gbk_file} not found!"
  elif [ ! -s "${gbk_file}" ]; then
    echo "Error: ${gbk_file} is empty!"
  elif [ ! -r "${gbk_file}" ]; then
    echo "Error: ${gbk_file} is unreadable!"
  else
    echo "${gbk_file} is present and readable."
  fi

  if [ ! -f "${gbf_file}" ]; then
    echo "Error: ${gbf_file} not found!"
  elif [ ! -s "${gbf_file}" ]; then
    echo "Error: ${gbf_file} is empty!"
  elif [ ! -r "${gbf_file}" ]; then
    echo "Error: ${gbf_file} is unreadable!"
  else
    echo "${gbf_file} is present and readable."
  fi

  if [ ! -f "${faa_file}" ]; then
    echo "Error: ${faa_file} not found!"
  elif [ ! -s "${faa_file}" ]; then
    echo "Error: ${faa_file} is empty!"
  elif [ ! -r "${faa_file}" ]; then
    echo "Error: ${faa_file} is unreadable!"
  else
    echo "${faa_file} is present and readable."
  fi

  OUT="${PROKKA_DIR}/prokka_gff"
  OUT1="${PROKKA_DIR}/prokka_gbk"
  OUT2="${PROKKA_DIR}/prokka_gbf"
  OUT3="${PROKKA_DIR}/prokka_faa"
  mkdir -p "${OUT}"
  mkdir -p "${OUT1}"
  mkdir -p "${OUT2}"
  mkdir -p "${OUT3}"

  cat_gff="${OUT}/${sample_name}.gff"
  cat_gbk="${OUT1}/${sample_name}.gbk"
  cat_gbf="${OUT2}/${sample_name}.gbf"
  cat_faa="${OUT3}/${sample_name}.faa"

  echo "Concatenating files for ${sample_name}:"
  echo "Gff file: ${gff_file}"
  echo "Gbk file: ${gbk_file}"
  echo "Gbf file: ${gbf_file}"
  echo "Faa file: ${faa_file}"

  cat "${gff_file}" > "${cat_gff}"
  cat "${gbk_file}" > "${cat_gbk}"
  cat "${gbf_file}" > "${cat_gbf}"
  cat "${faa_file}" > "${cat_faa}"
  echo "Concatenating is complete for ${sample_name}"
  echo "Concatenated assembly: ${cat_gff}"
  echo "Concatenated assembly: ${cat_gbk}"
  echo "Concatenated assembly: ${cat_gbf}"
  echo "Concatenated assembly: ${cat_faa}"
done

echo "Concatenation for all samples is complete!"
```
4. Renaming gff, gbk and gbf files with corrected header inside out

```
#!/bin/bash

# Define directories
WORK_DIR="/var/scratch/${USER}/viral_genomes/mpox_data"
GFF_DIR="${WORK_DIR}/mpox_results/prokka_data/prokka_gff"
GBK_DIR="${WORK_DIR}/mpox_results/prokka_data/prokka_gbk"
GBF_DIR="${WORK_DIR}/mpox_results/prokka_data/prokka_gbf"

# Loop through each sample in the GFF_DIR
for gff_file in "${GFF_DIR}"/*.gff; do
  # Extract the sample name from the file name
  sample_name=$(basename "${gff_file}" .gff) 
  gbk_file="${GBK_DIR}/${sample_name}.gbk"
  gbf_file="${GBF_DIR}/${sample_name}.gbf"
  
  # Function to check files
  check_file() {
    local file=$1
    if [ ! -f "${file}" ]; then
      echo "Error: ${file} not found!"
    elif [ ! -s "${file}" ]; then
      echo "Error: ${file} is empty!"
    elif [ ! -r "${file}" ]; then
      echo "Error: ${file} is unreadable!"
    else
      echo "${file} is present and readable."
      return 0
    fi
    return 1
  }

  # Function to process GFF file
  process_gff() {
    local file=$1
    local tmp_file=$(mktemp)
    while IFS= read -r line; do
      if [[ $line == "##sequence-region"* ]]; then
        echo "##sequence-region ${sample_name}" >> "${tmp_file}"
      elif [[ $line == *"gnl|WHO|MPXV_1"* ]]; then
        echo "${line//gnl|WHO|MPXV_1/gnl|WHO|${sample_name}}" >> "${tmp_file}"
      else
        echo "${line}" >> "${tmp_file}"
      fi
    done < "${file}"
    mv "${tmp_file}" "${file}"
    echo "Updated ${file} with sample name ${sample_name}."
  }

  # Function to process GBK and GBF files
  process_gbk_gbf() {
    local file=$1
    local tmp_file=$(mktemp)
    while IFS= read -r line; do
      if [[ $line == LOCUS* ]]; then
        echo "${line/MPXV_1/${sample_name}}" >> "${tmp_file}"
      else
        echo "${line}" >> "${tmp_file}"
      fi
    done < "${file}"
    mv "${tmp_file}" "${file}"
    echo "Updated ${file} with sample name ${sample_name}."
  }

  # Check and process GFF file
  if check_file "${gff_file}"; then
    process_gff "${gff_file}"
  fi

  # Check and process GBK file
  if check_file "${gbk_file}"; then
    process_gbk_gbf "${gbk_file}"
  fi

  # Check and process GBF file
  if check_file "${gbf_file}"; then
    process_gbk_gbf "${gbf_file}"
  fi
done
```

## B. Variants calling
Involves using snippy and variants annotation/prediction using snpeff, extraction using snpsift

1. Variants calling using snippy

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Snippy_Mpox
#SBATCH -n 4
#SBATCH --mem=8G

# Exit with error reporting
set -e
set -x

# load required modules
module purge
module load perl/5.22.3
module load snippy/4.6.0
module load python/3.9

# Directories path
WORK_DIR="/var/scratch/${USER}/viral_genomes/mpox_data"
DATA_DIR="${WORK_DIR}/mpox_fasta"
REF_DIR="${WORK_DIR}/mpox_refseqs"
OUT_DIR="${WORK_DIR}/mpox_results/mpox_snippy"

# Define ref
reference="${REF_DIR}/Mpox_ref_NC_063383.1.fasta"

# Make directories
if [ ! -e "${OUT_DIR}" ]; then
  mkdir -p "${OUT_DIR}"
fi

# Run snippy
for fasta_file in "${DATA_DIR}"/*.fasta; do
    sample=$(basename "${fasta_file}" .fasta)

    # create output directory for every sample
    OUT="${OUT_DIR}/${sample}"

    if [ ! -e "${OUT}" ]; then
        mkdir -p "${OUT}"
    fi

    vcf="${OUT}/${sample}.vcf"
    bed="${OUT}/${sample}.bed"
    gff="${OUT}/${sample}.gff"
    csv="${OUT}/${sample}.csv"
    html="${OUT}/${sample}.html"

    if [ ! -f "${vcf}" ] && [ ! -f "${bed}" ] && [ ! -f "${gff}" ] && [ ! -f "${csv}" ] && [ ! -f "${html}" ]; then
        echo -e "Calling variants on sample:\t${sample}; Fasta file: ${fasta_file}"
        snippy \
            --cpus 4 \
            --ram 8 \
            --prefix "${sample}" \
            --cleanup \
            --mapqual 60 \
            --basequal 15 \
            --mincov 10 \
            --force \
            --outdir "${OUT}" \
            --ref "${reference}" \
            --ctgs "${fasta_file}"
    fi
done

# run snippy-core
dirs=()
for d in "${OUT_DIR}"/*; do
    if [ -d "${d}" ] && [ "${d}" != 'reference' ]; then
        dirs+=("${d}")
    fi
done

DIRS=$(printf '%s ' "${dirs[@]}")
echo "${DIRS}"

CORE_DIR="${OUT_DIR}/snippy-core"

if [ ! -e "${CORE_DIR}" ]; then
    mkdir -p "${CORE_DIR}"
fi

cd "${CORE_DIR}"
echo -e "Running snippy-core"
snippy-core --ref "${reference}" --prefix core ${DIRS}
```
2. Variants annotation, prediction and extraction using snpeff and snpsift

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J mpox_snp_eff_sift
#SBATCH -n 4

# Exit immediately if a command exits with a non-zero status
set -e

# Load necessary modules
module purge
module load samtools/1.15.1
module load bcftools/1.15.1
module load snpeff/4.1g
module load java/17

# Directories path
WORK_DIR="/var/scratch/${USER}/viral_genomes/mpox_data"
REF_DIR="${WORK_DIR}/mpox_refseqs"
SNIPPY_DIR="${WORK_DIR}/mpox_results/mpox_snippy"

#snpEff config for Mpox virus
echo "Processing snpEff config for mpox"
mkdir -p "./database/snpEff/data/NC063383/"
cp -rf $REF_DIR/Mpox_ref_NC_063383.1.gff3 ./database/snpEff/data/NC063383/genes.gff
cp -rf $REF_DIR/Mpox_ref_NC_063383.1.fasta ./database/snpEff/data/NC063383/sequences.fa
echo -e "# Mpox viral genome, version MpoxNC063383\nNC063383.genome: NC063383" > ./database/snpEff/data/NC063383/snpEff.config

#Build the database
java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar build \
        -config ./database/snpEff/data/NC063383/snpEff.config \
        -dataDir ./../ \
        -gff3 \
        -v NC063383

# Define the output prefix for the index files and its directory
reference_prefix="${output_dir}/NC063383/NC063383"
mkdir -p "${reference_prefix}"

#Step 5: Compress the filtered VCF files #bgzip -c file.vcf > file.vcf.gz
echo "Compressing filtered VCF files..."
for sample_dir in "${SNIPPY_DIR}"/*; do
    sample_name=$(basename "${sample_dir}")
    filtered_vcf="${sample_dir}/${sample_name}.filt.vcf"
    echo "Compressing: ${filtered_vcf}"
    bgzip "${filtered_vcf}"
done

# Step 6: Index the filtered.vcf.gz VCF files
echo "Indexing filtered VCF files..."
for sample_dir in "${SNIPPY_DIR}"/*; do
    sample_name=$(basename "${sample_dir}")
    filtered_vcf_gz="${sample_dir}/${sample_name}.filt.vcf.gz"
    echo "Indexing: ${filtered_vcf_gz}"
    bcftools index -c --threads 4 "${filtered_vcf_gz}"
done

# Step 7: Variant Annotation with GFF file
# Define input files and directories
vcf_dir="${SNIPPY_DIR}"
annotated_dir="./mpox_results/annotated_mpox_variants_snippy"
reference_file="./database/snpEff/data/NC063383/sequences.fa"
gff_file="../database/snpEff/data/NC063383/genes.gff"
snpeff_jar="/export/apps/snpeff/4.1g/snpEff.jar"

# Create output directory for annotated variants
mkdir -p "$annotated_dir"

# Run SnpEff for variant annotation 
echo "Performing variant annotation..."
for sample_dir in "${SNIPPY_DIR}"/*; do
    sample_name=$(basename "${sample_dir}")
    vcf_file="${sample_dir}/${sample_name}.filt.vcf.gz"
    echo "Processing ${sample_name}: ${vcf_file}"
    bcftools view --threads 4 "$vcf_file" |
    java -Xmx4g -jar "$snpeff_jar" \
           -config ./database/snpEff/data/NC063383/snpEff.config \
           -dataDir ./../ \
           -v NC063383 ${vcf_file} > "${annotated_dir}/${sample_name}.snpEff.vcf"

echo "Variant annotation complete."

# Step 8: Rename summary.html and genes.txt and zip vcf files
mv ./snpEff_summary.html "${annotated_dir}/${sample_name}.snpEff.summary.html"
mv ./snpEff_genes.txt "${annotated_dir}/${sample_name}.snpEff.genes.txt"

# Compress vcf
bgzip -c "${annotated_dir}/${sample_name}.snpEff.vcf" > "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Create tabix index - Samtools
tabix -p vcf -f "${annotated_dir}/${sample_name}.snpEff.vcf.gz"

# Generate VCF files
bcftools stats "${annotated_dir}/${sample_name}.snpEff.vcf.gz" > "${annotated_dir}/${sample_name}.snpEff.stats.txt"

echo "Variant summary renaming, compressing complete."

done

# Step 9: Variant Extraction with SnpSift
# Define input files and directories
snpeff_dir="${annotated_dir}"
extracted_dir="./mpox_results/extracted_variants_snippy"

# Create output directory for extracted variants
mkdir -p "$extracted_dir"

# Run SnpSift for variant extraction #${sample}.filtered.snpEff.genes.txt
echo "Performing variant extraction..."
for snpeff_file in "$snpeff_dir"/*.snpEff.vcf.gz; do
    sample_name=$(basename "${snpeff_file}" .snpEff.vcf.gz)
    echo "Processing ${sample_name}: ${snpeff_file}"
    bcftools view --threads 4 "$snpeff_file" |
    java -Xmx4g -jar "/export/apps/snpeff/4.1g/SnpSift.jar" \
        extractFields \
        -s "," \
        -e "." \
        /dev/stdin \
        "ANN[*].GENE" "ANN[*].GENEID" \
        "ANN[*].IMPACT" "ANN[*].EFFECT" \
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
        "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" \
        "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" \
        > "${extracted_dir}/${sample_name}.snpsift.txt"
done

echo "SnpSift variant extraction complete"
```
3. Extract the variants for analysis and visualization in R

```
#!/usr/bin/Rscript

# Load necessary modules
module purge
module load R/4.3

# Open R console/studio in the terminal
R

# Remove all objects saved in R workspace
rm(list = ls())

# Set the input directory where the sample directories are located
input_dir <- "./mpox_results/extracted_variants_snippy"

# Get a list of all sample directories in the input directory
file_list <- list.files(input_dir, pattern = "\\.snpsift\\.txt", full.names = TRUE)

# Create an empty data frame to store the aggregated data
agg_df <- data.frame()

# Loop over each annotated result file
for (file in file_list) {
  # Read the annotated data from the file
  annotated_data <- read.delim(file, sep = "\t", header = TRUE)

  # Skip sample if no annotated data was found
  if (nrow(annotated_data) == 0) {
    next
  }

  # Extract the sample name from the file path
  sample_name <- sub("\\.snpsift\\.txt$", "", basename(file))

  # Add the sample name as a column in the data frame
  annotated_data$sample_name <- sample_name

  # Append the data to the aggregated data frame
  agg_df <- rbind(agg_df, annotated_data)
}

# Set the output file name
output_file <- "combined_annotated_data.csv"

# Check if the aggregated data frame is empty
if (nrow(agg_df) > 0) {
  # Write the aggregated data to a CSV file
  write.csv(agg_df, file = output_file, row.names = FALSE)
  cat("Variant extraction complete. Data saved to:", output_file, "\n")
} else {
  cat("No data extracted. Output file not created.\n")
}
```
