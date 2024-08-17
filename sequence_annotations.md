# Sequence Annotations

## A. Prokka
```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J Em_prk_mpox
#SBATCH -n 3

# Exit with error reporting
set -e

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
        output_path="${OUT_DIR}/${sample}"
        clean_fasta="${FASTA_DIR}/${sample}.clean.fasta"
        compressed_fasta="${OUT_DIR}/${sample}.clean.fasta.gz"
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

        # Clean the FASTA headers
        echo "Cleaning headers in ${fasta_file}"
        awk '/^>/ {gsub(/[|\\/\\-]/, "_", $0)} {print}' "$fasta_file" > "$clean_fasta"

        # Compress the cleaned fasta file
        echo "Compressing ${clean_fasta} to ${compressed_fasta}"
        gzip -c "$clean_fasta" > "$compressed_fasta"

        # Index the compressed fasta file
        echo "Indexing ${compressed_fasta}"
        samtools faidx "$compressed_fasta"

        # Echo fasta file path for the current sample
        echo "contigs_fasta for ${sample}: ${clean_fasta}"

         # Run the Prokka command on the current file
        echo "Processing Prokka for ${sample}: ${clean_fasta}"
        prokka "$clean_fasta" \
            --outdir "$output_path" \
            --cpus 3 \
            --mincontiglen 200 \
            --kingdom Viruses \
            --centre WHO \
            --addgenes \
            --addmrna \
            --locustag MPXV \
            --genus "Mpox" \
            --proteins "$protein_db" \
            --usegenus \
            --compliant \
            --rfam \
            --force \
            --debug

        # Clean up all files except the cleaned, compressed, and indexed files
        rm -f "$temp_fasta"
        rm -f "$fasta_file"
        rm -f "$gz_file"

    else
        echo "ERROR!!! fasta.gz file not found: ${gz_file}"
    fi
done

# Remove any remaining old files in the input directory
rm -f "${FASTA_DIR}"/*.fasta
rm -f "${FASTA_DIR}"/*.gzi
rm -f "${FASTA_DIR}"/*.fai
```
