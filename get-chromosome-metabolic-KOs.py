#!/usr/bin/env python

"""
get-plasmid-metabolic-KOs.py by Rohan Maddamsetti.

Usage: python get-chromosome-metabolic-KOs.py

I submitted ../results/unique_chromosome_KEGG_IDs.tsv to the KEGG Mapper
Reconstruct webserver, and saved a list of plasmid 2086 KO IDs that map to
KEGG Pathway 01100 "Metabolic pathways", in the file ../results/unique_chromosome_metabolic_KEGG_IDs.txt.

This script filters ../results/successful_chromosome_GhostKOALA_ko_results.tsv
for chromosome proteins with KOs in the list in
../results/unique_plasmid_metabolic_KEGG_IDs.txt.

"""


def main():

    """ IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   """
    replicon_length_file = "../results/replicon-lengths-and-protein-counts.csv"

    ## for debugging other scripts: what is the set of SeqID prefixes?
    NCBI_prefixes = set()
    
    metabolic_KEGG_IDs = set()
    with open("../results/unique_chromosome_metabolic_KEGG_IDs.txt", "r") as metabolic_KEGG_fh:
        for line in metabolic_KEGG_fh:
            metabolic_KEGG_ID = line.strip()
            metabolic_KEGG_IDs.add(metabolic_KEGG_ID)

    line_buffer = [] ## This is to save the lines to write to file.
    with open("../results/successful_chromosome_GhostKOALA_ko_results.tsv") as chromosome_GhostKOALA_fh:
        for line in chromosome_GhostKOALA_fh:
            line = line.strip()
            gene_string, KO_ID = line.split("\t")
            AnnotationAccession, SeqID, SeqType, locus_tag, product, = [ x.split("=")[-1] for x in gene_string.split("|") ]
            NCBI_prefix = SeqID.split("_")[0]
            NCBI_prefixes.add(NCBI_prefix)

            if KO_ID in metabolic_KEGG_IDs:
                line_buffer.append("\t".join([AnnotationAccession, SeqID, SeqType, locus_tag, product]))

    ##  write line_buffer to the output file.
    with open("../results/chromosome-proteins-in-KEGG-metabolism.tsv", "w") as outfh:
        outfh.write("AnnotationAccession\tSeqID\tSeqType\tlocus_tag\tproduct\n")
        for line in line_buffer:
            outfh.write(line + "\n")

    ## print the set of NCBI prefixes.
    print("THE SET OF NCBI PREFIXES:")
    print(NCBI_prefixes)
    print("*******")

    return


if __name__ == "__main__":
    main()
