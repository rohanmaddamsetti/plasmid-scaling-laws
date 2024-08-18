#!/usr/bin/env bash

## concatenate-and-filter-plasmid-GhostKOALA-results.sh

## concatenate the results of the ten batches.
cat ../results/plasmid-GhostKOALA-output-files/GhostKOALA_batch1_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch2_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch3_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch4_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch5_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch6_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch7_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch8_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch9_ko_results.tsv ../results/GhostKOALA-output-files/GhostKOALA_batch10_ko_results.tsv > ../results/plasmid_GhostKOALA_concatenated_ko_results.tsv
## filter for rows with successful KEGG mappings
cat ../results/plasmid_GhostKOALA_concatenated_ko_results.tsv | grep "\tK" > ../results/successful_plasmid_GhostKOALA_concatenated_ko_results.tsv
