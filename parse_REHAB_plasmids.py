#!/usr/bin/env python

"""
parse_REHAB_plasmids.py by Rohan Maddamsetti.

This script makes a csv file called "NCBI_REHAB_plasmid_metadata.csv"
with the data needed to compare my PCN estimates against the PCN estimates
published in Supplementary Table 2 of L. Shaw et al. (2021) in Science Advances.
These PCN data are in ../data/Shaw2021_PCN_data.csv.

We filter complete-prokaryotes-with-plasmids.txt to get all rows with the phrase "The REHAB Consortium",
and parse each row for RefSeq AnnotationAccession, Strain, plasmid name, and plasmid SeqID.

IMPORTANT: We assume all of these plasmids have a SeqID starting with an "NZ_" prefix.

Usage: python parse_REHAB_plasmids.py > ../results/NCBI_REHAB_plasmid_metadata.csv
"""


def main():
    print("AnnotationAccession,Strain,SeqType,SeqName,SeqID")
    complete_prokaryotes_file = "../results/complete-prokaryotes-with-plasmids.txt"
    with open(complete_prokaryotes_file, "r") as input_fh:
        for i,line in enumerate(input_fh):
            line = line.strip() ## remove leading/lagging whitespace
            if "The REHAB Consortium" in line:
                ## these assertions should be true based on previous filtering of prokaryotes.txt
                assert "Complete Genome" in line
                assert "plasmid" in line
                fields = line.split("\t")
                
                ftp_path = fields[-3]
                annotation_accession = ftp_path.split("/")[-1]
                ## This assertion should follow from upstream processing.
                assert annotation_accession.startswith("GCF_")  
                strain = fields[-1]

                replicon_list = fields[8].split(";")                
                for replicon_string in replicon_list:
                    ## only consider plasmids.
                    if replicon_string.startswith("chromosome"): continue
                    SeqType, SeqName_ID_pair = replicon_string.split()
                    SeqName, raw_SeqID = SeqName_ID_pair.split(":")
                    ## some processing to make the SeqID nice.
                    SeqID = "NZ_" + raw_SeqID.split("/")[-1]
                    metadata_row = ",".join([annotation_accession, strain, SeqType, SeqName, SeqID])
                    print(metadata_row)

if __name__ == "__main__":
    main()
