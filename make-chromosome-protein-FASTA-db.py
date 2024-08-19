#!/usr/bin/env python

"""
make-chromosome-protein-FASTA-db.py by Rohan Maddamsetti.

This script goes through all the gbk annotation files in ../results/gbk-annotation,
and makes a FASTA database of all proteins found on plasmids.

IMPORTANT TODO: pick a better, more representative set of genomes spanning 

"""

import os
from tqdm import tqdm
from Bio import SeqIO
import polars as pl


CHROMOSOME_MIN_LENGTH = 500000 ## a small endosymbiont has a 515,000 bp major chromosome

def select_genomes_for_chromosome_control(replicon_length_file, ecological_annotation_file,NUM_GENOMES=100):
    """ choose a set of 100 genomes, distributed equally over the rank
    distribution of chromosome lengths, over the set of genomes
    with ecological annotation.

    NOTE: This is approximate, since we don't focus on the major chromosome in each genome
    (each row is one chromosome in a genome, including secondary chromosomes), but should be good enough.
    """

    replicon_length_df = pl.read_csv(replicon_length_file)
    ecological_annotation_df = pl.read_csv(ecological_annotation_file)

    annotated_replicon_df = replicon_length_df.join(
        ecological_annotation_df,
        on="AnnotationAccession", how="left", coalesce=True).filter(
            ## filter out genomes with NA or blank Annotation
            ## then filger for chromosomes
            ~pl.col("Annotation").is_in(["NA", "blank"])).filter(
                pl.col("SeqType") == "chromosome").filter(
                    ## focus on large major chromosomes
                    pl.col("replicon_length") > CHROMOSOME_MIN_LENGTH).sort(
                    ## sort by replicon_length
                    "replicon_length")

    annotation_accession_list = annotated_replicon_df["AnnotationAccession"].to_list()
    annotation_accession_list_length = len(annotation_accession_list)
    my_modulus = int(annotation_accession_list_length/NUM_GENOMES)

    selected_genomes_list = list()
    ## select NUM_GENOMES that are roughly equally apart in terms of size range.
    for i, annotation_accession in enumerate(annotation_accession_list):
        if i % my_modulus == 0:
            selected_genomes_list.append(annotation_accession)

    return selected_genomes_list


def find_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file):
    bad_genomes_list  = []
    cur_annotation_accession = None
    cur_chromosome_length = -1
    with open(replicon_length_file, "r") as replicon_length_fh:
        for i, line in enumerate(replicon_length_fh):
            if i == 0: continue ## skip the header
            line = line.strip() ## remove leading and lagging whitespace.
            AnnotationAccession, SeqID, SeqType, replicon_length, protein_count = line.split(",")
            replicon_length = int(replicon_length)
            protein_count = int(protein_count)
            if AnnotationAccession != cur_annotation_accession and SeqType == "chromosome":
                cur_annotation_accession = AnnotationAccession
                cur_chromosome_length = replicon_length
            else: ## AnnotationAccession == cur_annotation_accession
                if (SeqType == "plasmid") and (replicon_length > cur_chromosome_length) and (cur_annotation_accession not in bad_genomes_list):
                    bad_genomes_list.append(cur_annotation_accession)
    return bad_genomes_list


def generate_chromosome_protein_fasta_db(refgenomes_dir, fasta_outfile, selected_genomes_list, bad_genomes_list):
    """
    IMPORTANT: This function only writes out proteins that are found on plasmids.
    """
    
    with open(fasta_outfile, "w") as outfh:
        gbk_filelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff")]
        good_gbk_filelist = [x for x in gbk_filelist if x not in bad_genomes_list]
        for gbk_file in tqdm(good_gbk_filelist):
            gbk_path = os.path.join(refgenomes_dir, gbk_file)
            genome_id = gbk_file.split(".gbff")[0]

            AnnotationAccession = genome_id.split("_genomic")[0]
            """ skip genomes with a chromosome smaller than a plasmid in the genome.
            This is probably a RefSeq genome with a 'plasmid' that is actually the chromosome.
            """
            if AnnotationAccession in bad_genomes_list: continue
            """ skip genomes that are not in the list of selected genomes."""
            if AnnotationAccession not in selected_genomes_list: continue
            
            with open(gbk_path,'r') as gbk_fh:
                SeqID = None
                SeqType = None
                for record in SeqIO.parse(gbk_fh, "genbank"):
                    SeqID = record.id
                    if "chromosome" in record.description:
                        SeqType = "chromosome"
                    else:
                        continue
                    for feature in record.features:
                        ## only analyze protein-coding genes.
                        if feature.type != "CDS": continue
                        ##print(feature.qualifiers)
                        locus_tag = feature.qualifiers["locus_tag"][0]
                        ## for downstream parsing, replace spaces with underscores in the product annotation field.
                        product = feature.qualifiers["product"][0].replace(" ","_")
                        ## skip over CDS that don't have an annotated translation.
                        if "translation" not in feature.qualifiers: continue
                        protein_seq = feature.qualifiers["translation"][0]
                        header = ">" + "|".join(["AnnotationAccession="+AnnotationAccession,"SeqID="+SeqID,"SeqType="+SeqType,"locus_tag="+locus_tag,"product="+product])
                        outfh.write(header + "\n")
                        outfh.write(str(protein_seq) + "\n")
    return


def main():
    refgenomes_dir = "../results/gbk-annotation/"
    replicon_length_file = "../results/replicon-lengths-and-protein-counts.csv"

    """ make a list of genomes in which the chromosome is smaller than some plasmid.
    These are often misannotated genomes in which the real chromosome is switched up with a plasmid."""
    bad_genomes_list = find_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file)

    print("*******")
    print("BAD GENOMES WITH A PLASMID APPARENTLY LARGER THAN THE CHROMOSOME:")
    for x in bad_genomes_list:
        print(x)
    print("*******")

    ecological_annotation_file = "../results/computationally-annotated-gbk-annotation-table.csv"
    selected_genomes_list = select_genomes_for_chromosome_control(replicon_length_file, ecological_annotation_file)
    
    """ make a database of proteins found on chromosomes in the selected genomes
    as input to GhostKOALA to get KEGG IDs.  """
    fasta_outfile = "../results/chromosome_GhostKOALA_input.fasta"
    generate_chromosome_protein_fasta_db(refgenomes_dir, fasta_outfile, selected_genomes_list, bad_genomes_list)
    return


## run the script.
if __name__ == "__main__":
    main()
