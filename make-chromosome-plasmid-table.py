#!/usr/bin/env python

'''
make-chromosome-plasmid-table.py by Rohan Maddamsetti.

This script reads in ../results/complete-prokaryotes-with-plasmids.txt.
'''
import os

warning_count = 0
found_genomes = 0
accessions_not_found_in_refseq = []

with open("../results/chromosome-plasmid-table.csv",'w') as out_fh:
    header = "AnnotationAccession,SeqID,SeqType,Organism,Strain\n"
    out_fh.write(header)
    ## open the genome report file, and parse line by line.
    with open("../results/complete-prokaryotes-with-plasmids.txt", "r") as genome_report_fh:
        for i, line in enumerate(genome_report_fh):
            line = line.strip()
            if i == 0: ## get the names of the columns from the header.
                column_names_list = line.split('\t')
                continue ## don't process the header any further.
            fields = line.split('\t')
            organism = fields[0]
            strain = fields[-1]
            replicons = fields[8]
            ftp_path = fields[20]
            annotation_accession = os.path.basename(ftp_path)
            my_annotation_file = "../results/gbk-annotation/" + annotation_accession + "_genomic.gbff"
            ''' make sure that this file exists in the annotation directory--
            skip if this was not the case.
            this is important; we don't want to include genomes that were
            not in the search database in the first place. '''
            if not os.path.exists(my_annotation_file):
                warning_count += 1
                accessions_not_found_in_refseq.append(annotation_accession)
                continue
            ## if we got here, then the RefSeq genbank annotation file was found.
            found_genomes += 1
            ## get the SeqID and SeqType from the Genbank annotation file.


            replicon_list = replicons.split(';')
            for s in replicon_list:
                s = s.strip()
                seq_name, raw_seq_id_string = s.split(":")
                seq_type = seq_name.split(" ")[0]
                seq_id = raw_seq_id_string.split("/")[0]
                if seq_type not in ("chromosome", "plasmid"):
                    print("PROBABLE PHAGE (not chromosome or plasmid):")
                    print(s)
                    continue




                row_data = [annotation_accession, seq_id, seq_type, organism, strain]
                ## replace all commas with semicolon to respect the csv output format.
                row_string = ','.join([x.replace(',',';') for x in row_data]) + '\n'
                out_fh.write(row_string)

## now print out warnings.
print("Warning: the following accessions were not found in RefSeq. ")
for annotation_accession in accessions_not_found_in_refseq:
    print(annotation_accession)
print("A total of " + str(warning_count) + " accessions were not found in RefSeq.")
print("The most likely explanation is that these missing accessions did not pass RefSeq Quality Control.")
print("A total of " + str(found_genomes) + " genomes were found in RefSeq and processed.")
