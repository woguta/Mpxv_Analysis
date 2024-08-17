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
OUTPUT_DIR="${WORK_DIR}/mpox_results/prokka_data"
PROKKA_DIR="${WORK_DIR}/mpox_results/mpox_prokka"

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
WORK_DIR="/home/woguta/mpox_files"
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
