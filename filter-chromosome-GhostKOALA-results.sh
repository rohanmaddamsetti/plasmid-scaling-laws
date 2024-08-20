#!/usr/bin/env bash

## concatenate-and-filter-plasmid-GhostKOALA-results.sh

## filter for rows with successful KEGG mappings
cat ../results/chromosome_GhostKOALA_ko_results.tsv | grep "\tK" > ../results/successful_chromosome_GhostKOALA_ko_results.tsv
