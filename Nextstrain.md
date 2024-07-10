# Phylogenetic and evolutionary analysis

## A. Nextstrain
Installing nextstrain
1. Install Nextstrain CLI
```
curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash
```
2. Set up a Nextstrain runtime
```
nextstrain setup --set-default conda
```
3. Enter an interactive Nextstrain shell in the current directory (.).
```
nextstrain shell .
```
4. Run Augur.
```
augur --help
```
5. Run Auspice.
```
auspice --help
```
6. Exit the Nextstrain shell.
```
exit
```
7. Index the sequences
```
augur index \
  --sequences ./mpox.fasta \
  --output results/mpox_sequences_index.tsv
```
8. Filter the Sequences
```
augur filter \
  --sequences ./mpox.fasta \
  --sequence-index ./results/mpox_sequences_index.tsv \
  --metadata ./mpox_212_tree_labels.csv \
  --output ./results/mpox_filtered.fasta \
  --group-by country year \
  --sequences-per-group 20 \
  --min-date 1970
```
9. Filter the Sequences to drop some already analyzed seqs
```
augur filter \
  --sequences ./mpox.fasta \
  --sequence-index results/mpox_sequences_index.tsv \
  --metadata ./mpox_212_tree_labels.csv \
  --exclude config/dropped_strains.txt \
  --output ./results/mpox_filtered.fasta \
  --group-by country year month \
  --sequences-per-group 20 \
  --min-date 1970
```
10. Align the Sequences
```
augur align \
  --sequences ./results/mpox_filtered.fasta \
  --reference-sequence ./ref_mpox.fasta \
  --output ./results/mpox_aligned.fasta \
  --fill-gaps
```

11. Get a Time-Resolved Tree
    
NB Remove fasta seq length in the header

Shell
```
#!/bin/bash

# Input and output file paths
input_file="./mpox_data/mpox_270_aligned.fasta"
output_file="./mpox_data/mpox_270_cleaned_aligned.fasta"

# Remove sequence length from headers and write to output file
sed 's/\/[0-9]\+-[0-9]\+//' "$input_file" > "$output_file"
```
R script
```
# Read the FASTA file
fasta_file <- "./mpox_data/mpox_270_aligned.fasta"
fasta_content <- readLines(fasta_file)

# Remove sequence length from headers
cleaned_content <- gsub("/[0-9]+-[0-9]+", "", fasta_content)

# Write the cleaned content back to the file
writeLines(cleaned_content, "mpox_270_cleaned_aligned.fasta")
```
Then use

```
augur refine \
  --tree ./results/mpox_fasttree_data_58.nwk \
  --alignment ./results/mpox_aligned_data_58.fasta \
  --metadata ./mpox_212_tree_labels.csv \
  --output-tree ./results/mpox_tree.nwk \
  --output-node-data ./results/mpox_branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4
```

OR

```
augur refine \
  --tree ./results/mpox_fasttree_data_58.nwk \
  --alignment ./results/mpox_aligned_data_58.fasta \
  --metadata ./mpox_212_tree_labels.csv \
  --output-tree ./results/mpox_tree.nwk \
  --output-node-data ./results/mpox_branch_lengths.json \
  --timetree \
  --year-bounds 1970 2024 \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4
```
OR - this is the best option

```
augur refine \
  --tree ./mpox_data/mpox_272_tree.nwk \
  --alignment ./mpox_data/mpox_270_cleaned_aligned.fasta \
  --metadata ./mpox_data/mpox_data_270.csv \
  --output-tree ./results/mpox_tree.nwk \
  --output-node-data ./results/mpox_branch_lengths.json \
  --timetree \
  --stochastic-resolve \
  --clock-rate  0.001 \
  --clock-std-dev  0.0001 \
  --year-bounds 1970 2024 \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4
```
12. Annotate the Phylogeny

a. Reconstruct Ancestral Traits
```
augur traits \
  --tree ./results/mpox_tree.nwk \
  --metadata ./mpox_212_tree_labels.csv \
  --output-node-data ./results/mpox_traits.json \
  --columns clade date country \
  --confidence
```
b. Infer Ancestral Sequences
```
augur ancestral \
  --tree ./results/mpox_tree.nwk \
  --alignment ./results/mpox_aligned_data_58.fasta \
  --output-node-data ./results/mpox_nt_muts.json \
  --inference joint
```
c. Identify Amino-Acid Mutations ref nust be in gb file format
```
augur translate \
  --tree ./results/mpox_tree.nwk \
  --ancestral-sequences ./results/mpox_nt_muts.json \
  --reference-sequence ./mpox_ref.gb \
  --output-node-data ./results/mpox_aa_muts.json
```

d. Infer tip frequencies of mutations or clades

```
augur frequencies \
	--method kde \
	--metadata ./mpox_212_tree_labels.csv \
	--pivot-interval-units months \
        --pivot-interval 12 \
	--min-date 1970-01-01 \
	--max-date 2030-12-31 \
	--tree ./results/mpox_tree.nwk \
	--include-internal-nodes \
	--alignment ./results/mpox_aligned_data_58.fasta \
	--output ./results/mpox_tip-frequencies.json
```
 
13. Export the Results
Create the Config Directory
```
mkdir config
```
```
cd config
```
Upgrade Nextstrain
```
nextstrain update conda
```
Check how to use augur
```
augur export v2 --help
```
Create a json file and validate in json file validator
```
 {
  "title": "Phylodynamics of Monkeypox virus in Africa",
  "colorings": [
    {
      "key": "year",
      "title": "Year",
      "type": "continuous"
    },
    {
      "key": "host",
      "title": "Host",
      "type": "categorical"
    },
    {
      "key": "country",
      "title": "Country",
      "type": "categorical"
    },
    {
      "key": "clade",
      "title": "Clade",
      "type": "categorical"
    },
    {
      "key": "country",
      "title": "Country",
      "type": "categorical"
    }
  ],
  "geo_resolutions": [
      "Country"
    ],
  "maintainers": [
    {"name": "WHO/EPR/HIR"},
    {"name": "Walter Oguta"},
    {"name": "Rachel Achilla"}
    ],

  "filters": [
    "Year", "Clade", "Country", "Host"
  ],
  "display_defaults": {
    "layout": "rect",
    "color_by": "Clade",
    "branch_label": "Clade",
    "geo_resolution": "Country"

  }
}
```
Create a colours.tsv file
```
attribute   value       color
country     Nigeria     #FFFF00
country     Central_African_Republic #1CE6FF
country     Benin       #FF34FF
clade       Clade I     #809693
clade       Clade IIa   #D16100
clade       Clade IIb A #300018
date        1970-XX-XX  #FF913F
date        1971-XX-XX  #938A81
```
Create lat_long.tsv file
```
country Benin   9.961027        2.327362
country Central_African_Republic        7.0323598       19.9981227
country Cameroon        4.6125522       13.1535811
country Congo   -0.516362       15.877811
country DRC     -4.0383 21.7587
country Egypt   26.8206 30.8025
country Gabon   -0.8999695      11.6899699
country Liberia 5.7499721       -9.3658524
country Nigeria 9.6000359       7.9999721
country Sierra_Leone    8.6400349       -11.8400269
country South_Africa    -28.8166235     24.991639
country Sudan   14.5844444      29.4917691
```
Export the results into json file for view in auspice ```https://auspice.us```

```
augur export v2 \
  --tree ./results/mpox_tree.nwk \
  --metadata ./mpox_212_tree_labels.csv \
  --node-data ./results/mpox_branch_lengths.json \
  --node-data ./results/mpox_traits.json \
  --node-data ./results/mpox_nt_muts.json \
  --node-data ./results/mpox_aa_muts.json \
  --node-data ./results/mpox_tip-frequencies.json \
  --colors ./config/colors.tsv \
  --lat-longs ./config/lat_longs.tsv \
  --auspice-config ./config/auspice_config.json \
  --output ./results/auspice/mpox.json
```
##B Pan-Genome Analysis using Panx
Clone to conda terminal
git clone https://github.com/neherlab/pan-genome-analysis.git
