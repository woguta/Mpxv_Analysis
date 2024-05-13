# Phylogenetic And Evlutionary Analysis

## A. Nextstrain
## Installing nextstrain

## 1. Install Nextstrain CLI
```
curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash
```

## 2. Set up a Nextstrain runtime
```
nextstrain setup --set-default conda
```

## 3. Enter an interactive Nextstrain shell in the current directory (.).
```
nextstrain shell .
```
## 4. Run Augur.
```
augur --help
```
## 5. Run Auspice.
```
auspice --help
```

## 6. Exit the Nextstrain shell.
```
exit
```
## 6. Index the sequences
```
augur index \
  --sequences ./mpox.fasta \
  --output results/mpox_sequences_index.tsv
```
## 7. Filter the Sequences
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
## 8. Filter the Sequences to drop some already analyzed seqs
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
## 9. Align the Sequences
```
augur align \
  --sequences ./results/mpox_filtered.fasta \
  --reference-sequence ./ref_mpox.fasta \
  --output ./results/mpox_aligned.fasta \
  --fill-gaps
```

## 10. Get a Time-Resolved Tree
```
augur refine \
  --tree ./results/mpox_fasttree_data_38.nwk \
  --alignment ./results/Mpox_alignment_mafft3.fasta \
  --metadata ./mpox_212_tree_labels.csv \
  --output-tree ./results/mpox_tree.nwk \
  --output-node-data ./results/mpox_branch_lengths.json \
  --timetree \
  --coalescent opt \
  --date-confidence \
  --date-inference marginal \
  --clock-filter-iqd 4
```

## 11. Annotate the Phylogeny

## a. Reconstruct Ancestral Traits
```
augur traits \
  --tree ./results/mpox_tree.nwk \
  --metadata ./mpox_212_tree_labels.csv \
  --output-node-data ./results/mpox_traits.json \
  --columns region country \
  --confidence
```
## b. Infer Ancestral Sequences
```
augur ancestral \
  --tree ./results/mpox_tree.nwk \
  --alignment ./results/Mpox_alignment_mafft3.fasta \
  --output-node-data ./results/mpox_nt_muts.json \
  --inference joint
```
## c. Identify Amino-Acid Mutations ref nust be in gb file format
```
augur translate \
  --tree ./results/mpox_tree.nwk \
  --ancestral-sequences ./results/mpox_nt_muts.json \
  --reference-sequence ./mpox_ref.gb \
  --output-node-data results/mpox_aa_muts.json
```
## 12. Export the Results
## Create the Config Directory
```
mkdir config
```
```
cd config
```
```
nano auspice_config.json
```
## Verify
```
ls config
```
```
#augur export v2 --help
```
## Export results with
```
augur export v2 \
  --tree ./results/mpox_tree.nwk \
  --metadata ./mpox_212_tree_labels.csv \
  --node-data ./results/mpox_branch_lengths.json \
              ./results/mpox_traits.json \
              ./results/mpox_nt_muts.json \
              ./results/mpox_aa_muts.json \
  --colors config/colors.tsv \
  --lat-longs config/lat_longs.tsv \ 
  --output ./results/mpox.json
```
0R
```
augur export v2 \
  --tree ./results/mpox_tree.nwk \
  --metadata ./mpox_212_tree_labels.csv \
  --node-data ./results/mpox_branch_lengths.json \
  --node-data ./results/mpox_traits.json \
  --node-data ./results/mpox_nt_muts.json \
  --node-data ./results/mpox_aa_muts.json \
  --colors config/colors.tsv \
  --lat-longs config/lat_longs.tsv \
  --auspice-config config/auspice_config.json \
  --output auspice/mpox.json
```
