#!/usr/bin/env python

"""
 calculate-CDS-rRNA-fractions.py by Rohan Maddamsetti.

This script answers a basic question that Lingchong asked:
for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?

Repeat the calculation for 5S, 16S, 23S, rRNAs.
CORE ASSUMPTION: I am assuming these rRNAs are always found in archaea as well--
seems to be the case in known archaeal species at least.

Usage: python calculate-CDS-rRNA-fractions.py > ../results/CDS-rRNA-fractions.csv
"""

import os
import re
import gzip
from Bio import SeqIO


def main():
    ## print the header of the file
    print("AnnotationAccession,SeqID,SeqLength,CodingLength,coding_fraction,CDS_count,summed_CDS_length,summed_CDS_fraction,rRNA_5S_count,rRNA_5S_length,rRNA_5S_fraction,rRNA_16S_count,rRNA_16S_length,rRNA_16S_fraction,rRNA_23S_count,rRNA_23S_length,rRNA_23S_fraction,total_rRNA_count,total_rRNA_length,total_rRNA_fraction")

    gbk_annotation_path = "../results/gbk-annotation/"
    gbk_gz_files = [x for x in os.listdir(gbk_annotation_path) if x.endswith(".gbff.gz")]
    
    for gbk_gz in gbk_gz_files:        
        AnnotationAccession = gbk_gz.split("_genomic.gbff.gz")[0]
        gbk_gz_path = os.path.join(gbk_annotation_path, gbk_gz)
        with gzip.open(gbk_gz_path, "rt") as gbk_handle:
            records = SeqIO.parse(gbk_handle, 'genbank')
            ## Iterate through each SeqRecord in the GenBank file
            for record in records:
                total_sequence_length = len(record.seq)
                cds_count = 0
                summed_cds_length = 0
                ## count the number of distinct sites in protein-coding regions to
                ## account for overlapping genes.
                coding_positions_set = set()
                rRNA_5S_count = 0
                rRNA_5S_length = 0
                rRNA_16S_count = 0
                rRNA_16S_length = 0
                rRNA_23S_count = 0
                rRNA_23S_length = 0
                ## Iterate through features and calculate CDS length
                for feature in record.features:
                    my_seq_length = len(feature.location)
                    if feature.type == 'CDS':
                        cds_count += 1
                        for part in feature.location.parts:
                            ## just in case genes are discontinuous
                            ## (there shouldn't be introns but won't hurt)
                            summed_cds_length += len(part)
                            coding_positions_set.update(set(part))
                    elif feature.type == 'rRNA':
                        try: 
                            my_gene_product = feature.qualifiers['product'][0]
                        except:
                            raise AssertionError("ERROR: unknown rRNA!")
                        if my_gene_product == "5S ribosomal RNA":
                            rRNA_5S_count += 1
                            rRNA_5S_length += my_seq_length
                        elif my_gene_product == "16S ribosomal RNA":
                            rRNA_16S_count += 1
                            rRNA_16S_length += my_seq_length
                        elif my_gene_product == "23S ribosomal RNA":
                            rRNA_23S_count += 1
                            rRNA_23S_length += my_seq_length
                ## Calculate the fraction of sequence covered by CDS and rRNA.
                coding_length = len(coding_positions_set)
                ## calculate the sum of all cds lengths (ignoring overlaps)
                summed_cds_length
                coding_fraction = coding_length / total_sequence_length
                summed_CDS_fraction = summed_cds_length / total_sequence_length
                rRNA_5S_fraction = rRNA_5S_length / total_sequence_length
                rRNA_16S_fraction = rRNA_16S_length / total_sequence_length
                rRNA_23S_fraction = rRNA_23S_length / total_sequence_length

                total_rRNA_count = rRNA_5S_count + rRNA_16S_count + rRNA_23S_count
                total_rRNA_length = rRNA_5S_length + rRNA_16S_length + rRNA_23S_length
                total_rRNA_fraction = total_rRNA_length / total_sequence_length

                print(",".join([AnnotationAccession, record.id, str(total_sequence_length),
                                str(coding_length), str(coding_fraction),
                                str(cds_count), str(summed_cds_length), str(summed_CDS_fraction),
                                str(rRNA_5S_count), str(rRNA_5S_length), str(rRNA_5S_fraction),
                                str(rRNA_16S_count), str(rRNA_16S_length), str(rRNA_16S_fraction),
                                str(rRNA_23S_count), str(rRNA_23S_length), str(rRNA_23S_fraction),
                                str(total_rRNA_count), str(total_rRNA_length), str(total_rRNA_fraction)]))
    return


if __name__ == "__main__":
    main()
