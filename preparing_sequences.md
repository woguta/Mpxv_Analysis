# Prepare the sequences

1. Clean the seqs for easy use
```
#!/bin/bash

# Input and output file paths
work_dir=~/mpox_files
input_dir="$work_dir/mpox_sept/mpox_fasta"
input_file="$input_dir/mpox_all_480.fasta"
output_file="$input_dir/mpox_all_480_clean.fasta"

# Remove sequence length from headers and write to output file
sed 's/\/[0-9]\+-[0-9]\+//' "$input_file" > "${output_file}.tmp1"

# Replace special characters in headers with underscores and ensure proper FASTA format
awk '
    BEGIN {header_num=1}
    /^>/ {
        header=$0
        gsub(/[|\/\-]/, "_", header)
        print header
    }
    !/^>/ {print}
' "${output_file}.tmp1" > "${output_file}.tmp2"

mv "${output_file}.tmp2" "$output_file"
rm "${output_file}.tmp1"

# Ensure the final fasta file has a proper FASTA header
if ! grep -q "^>" "$output_file"; then
    echo "ERROR!!! The file ${output_file} does not have a proper FASTA header. Adding default headers."
    awk '
        BEGIN {header_num=1}
        !/^>/ {
            print ">sequence_" header_num
            header_num++
        }
        {print}
    ' "$output_file" > "${output_file}.tmp"
    mv "${output_file}.tmp" "$output_file"
    echo "Default headers added to ${output_file}"
else
    echo "The file ${output_file} has been cleaned and is in proper FASTA format"
fi
```
2. Split the combined fasta files
```
#!/bin/bash

# Define variables
WORK_DIR=~/mpox_files
INPUT_DIR="$WORK_DIR/mpox_sept/mpox_fasta"
FASTA_DIR="$WORK_DIR/mpox_sept/fasta_files"
INPUT_FILE="$INPUT_DIR/mpox_all_480_clean.fasta"

# Create directories
echo "Creating directories..."
mkdir -p "$FASTA_DIR"

# Display a message indicating splitting of FASTA files
echo "Splitting the combined FASTA file..."

# Change directory to the input directory
cd "$INPUT_DIR" || exit

# Extract original filenames without renaming
awk '/^>/{filename = substr($0, 2); gsub(/[\[\]\/\\&:\(\)\{\}<>!@#\$%\^*+=`|~'"'"';,?"'"'"']/, "_", filename); print > (filename ".fasta"); next} {print >> (filename ".fas>
# Compress each FASTA file using bgzip, create index using samtools faidx & unzip again
for fasta_file in *.fasta; do
    bgzip "$fasta_file"  # Compress FASTA file
    samtools faidx "${fasta_file}.gz"  # Create index file
    gunzip -k "${fasta_file}.gz"  # Unzip the fasta files
done

# Move specific files to fasta dir
mv -f *pxV_*.fasta "$FASTA_DIR"
mv -f *.fasta.* "$FASTA_DIR"

# Clean folder moving/copying
rm -f *pxV_*.fasta
rm -f *.fasta.*
```

3. Copying files to and from local comp/laptop and hpc 

i. From local comp to remote server (outside hpc)
```
scp -r ./mpox_files/mpox_august/fasta_files/* woguta@hpc.server.org:./mpox_data/mpox_fasta
```
ii. From remote server hpc to var/scartch for analysis (inside hpc)
```
cp -r ./mpox_data/mpox_fasta/* /var/scratch/woguta/viral_genomes/mpox_data/mpox_fasta
```
iii. From var/scratch to hpc homepage (inside hpc) /var/scratch/woguta/viral_genomes/mpox_data
```
cp -r /var/scratch/woguta/viral_genomes/mpox_data/mpox_results/prokka_data/ ~/mpox_data/
```
```
cp -r /var/scratch/woguta/viral_genomes/mpox_data/mpox_fasta/ ~/mpox_data/mpox_fasta
```
```
cp -r /var/scratch/woguta/viral_genomes/mpox_data/combined_annotated_data.csv ~/mpox_data/
```
iv. From remote server to homepage comp/laptop (outside hpc)
```
scp -r woguta@hpc.server.org:~/mpox_data/prokka_data/ ~/mpox_files/mpox_results/
```
```
scp -r woguta@hpc.server.org:~/mpox_data/mpox_fasta/ ~/mpox_files/mpoxv_pangenome_analysis/fasta_files/
```
```
scp -r woguta@hpc.server.org:~/mpox_data/combined_annotated_data.csv ~/mpox_files/mpox_results/
```
