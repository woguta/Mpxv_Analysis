## Mpox Sequence analysis from Illumina raw data (fastq files)
#The analysis pipeline follows raw sequence reads analysis from fastq to consensus fasta. Examples of tools

a. FastQC: Quality check of raw reads
b. Trimmomatic or fastp: Adapter trimming and quality filtering
c. hostile or bowtie or Kraken2 or BBMap (BBduk/BBsplit): Human read removal
d. Minimap2 or BWA mem: Mapping reads to references
e. Samtools: BAM processing and stats
f. Qualimap or Mosdepth: Mapping quality statistics
g. freebayes or snippy for calling variants
h. bcftools for consensus building

## 1. Modules needed
```
module load hostile/2.0.0
module load fastqc/0.11.9
module load fastp/0.24.1
module load seqtk/1.3
module load bwa/0.7.17
module load freebayes/1.3.4
module unload bcftools/1.17
module unload bcftools/1.13
```

