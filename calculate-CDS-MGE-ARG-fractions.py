#!/usr/bin/env python

"""
 calculate-CDS-MGE-ARG-fractions.py by Rohan Maddamsetti.

This script answers a basic question that Lingchong asked:
for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?

repeat the calculation for MGE-associated proteins and ARGs.

CRITICAL TODO: figure out how to add duplicated ARGs as well, accounting for duplicated ARGs
distributed over different replicons.

Usage: python calculate-CDS-MGE-ARG-fractions.py > ../results/CDS-MGE-ARG-fractions.csv
"""

import os
import re
from Bio import SeqIO


def isARG(product_annotation):
    chloramphenicol_keywords = "chloramphenicol|Chloramphenicol"
    tetracycline_keywords = "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
    MLS_keywords = "macrolide|lincosamide|streptogramin"
    multidrug_keywords = "Multidrug resistance|multidrug resistance|antibiotic resistance"
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*"
    glycopeptide_keywords = "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
    polypeptide_keywords = "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
    diaminopyrimidine_keywords = "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
    sulfonamide_keywords = "sulfonamide|Sul1|sul1|sulphonamide"
    quinolone_keywords = "quinolone|Quinolone|oxacin|qnr|Qnr"
    aminoglycoside_keywords = "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
    macrolide_keywords = "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
    antimicrobial_keywords = "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\S*"
    ARG_regex = "|".join([chloramphenicol_keywords, tetracycline_keywords,
                          MLS_keywords, multidrug_keywords, beta_lactam_keywords,
                          glycopeptide_keywords, polypeptide_keywords, diaminopyrimidine_keywords,
                          sulfonamide_keywords, quinolone_keywords, aminoglycoside_keywords,
                          macrolide_keywords, antimicrobial_keywords])
    if re.search(ARG_regex, product_annotation): return True
    return False


def isMGE(product_annotation):
    transposon_keywords = "IS|transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|tra[A-Z]|conjugate transposon|Transpos\S*|Tn[0-9]|tranposase|Tnp|Ins|ins"
    plasmid_keywords = "relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation|Mob\S*|Plasmid|Rep|Conjug\S*"
    phage_keywords = "capsid|phage|Tail|tail|head|tape measure|antiterminatio|Phage|virus|Baseplate|baseplate|coat|entry exclusion"
    other_HGT_keywords = "Integrase|integrase|excision\S*|exonuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip|intron|antitoxin|toxin|Toxin|Reverse transcriptase|hok|Hok|competence|addiction"
    MGE_regex = "|".join([transposon_keywords, plasmid_keywords, phage_keywords, other_HGT_keywords])
    if re.search(MGE_regex, product_annotation): return True
    return False


def get_prot_data(feature):
    try:
        prot_seq = feature.qualifiers['translation'][0]
    except:
        prot_seq = "NA"
        
    try: ## replace all commas with semicolons for csv formatting.
        prot_product = feature.qualifiers['product'][0].replace(',',';')
    except:
        prot_product = "NA"
        
    prot_location = str(feature.location)
                        
    cur_prot = { "seq" : prot_seq,
                 "product" : prot_product,
                 "location" : prot_location }
    return cur_prot


def main():
    gbk_annotation_path = "../results/gbk-annotation/"
    gbk_files = [x for x in os.listdir(gbk_annotation_path) if x.endswith(".gbff")]

    ## print the header of the file
    print("AnnotationAccession,SeqID,SeqLength,CDS_count,CDS_length,CDS_fraction,MGE_count,MGE_length,MGE_fraction,ARG_count,ARG_length,ARG_fraction")
    for gbk in gbk_files:
        AnnotationAccession = gbk.split("_genomic.gbff")[0]
        gbk_path = os.path.join(gbk_annotation_path, gbk)
        records = SeqIO.parse(gbk_path, 'genbank')

        ## Iterate through each SeqRecord in the GenBank file
        for record in records:
            total_sequence_length = len(record.seq)
            cds_count = 0
            cds_length = 0
            mge_count = 0
            mge_length = 0
            arg_count = 0
            arg_length = 0

            ## Iterate through features and calculate CDS length
            for feature in record.features:
                if feature.type == 'CDS':
                    my_seq_length = len(feature.location)
                    cds_count += 1
                    cds_length += my_seq_length
                    cur_prot_dict = get_prot_data(feature)
                    if isMGE(cur_prot_dict['product']):
                        mge_count += 1
                        mge_length += my_seq_length
                    elif isARG(cur_prot_dict['product']):
                        arg_count += 1
                        arg_length += my_seq_length
            
            ## Calculate the fraction of sequence covered by CDS
            CDS_fraction = cds_length / total_sequence_length
            ## Calculate the fraction of sequence covered by MGE-associated proteins
            MGE_fraction = mge_length / total_sequence_length
            ## Calculate the percentage of sequence covered by ARGs
            ARG_fraction = arg_length / total_sequence_length

            ## TODO: calculate the fraction of sequence covered by duplicated ARGs
            ##dupARG_fraction = dupARG_length / total_sequence_length
            
            
            print(",".join([AnnotationAccession, record.id, str(total_sequence_length), str(cds_count), str(cds_length), str(CDS_fraction), str(mge_count), str(mge_length), str(MGE_fraction), str(arg_count), str(arg_length), str(ARG_fraction)]))


if __name__ == "__main__":
    main()
