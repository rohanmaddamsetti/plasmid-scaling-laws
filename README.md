# plasmid-scaling-laws by Rohan Maddamsetti


## Software requirements
## Python 3.12+

## COMPUTATIONAL PROTOCOL

Make a top-level directory with three directories inside, named "data", "results", and "src".  
Now copy all source code files in this repository into "src".

Now, download prokaryotes.txt into ../data/:  

wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt  

Then, filter the prokaryotes.txt genome data for those that have complete genomes,
and replace "GCA" with "GCF" throughout this file, so that RefSeq data and not Genbank data
is accessed in all downstream steps:  

python filter-genome-reports.py > ../results/complete-prokaryotes-with-plasmids.txt.  

To preprocess metadata for PCN comparisons with the PCN data in Shaw et al. (2021) in
Science Advances, run the following command:  

python parse_REHAB_plasmids.py > ../results/NCBI_REHAB_plasmid_metadata.csv  

Then, fetch genome annotation for each row in complete-prokaryotes-with-plasmids.txt,
fetch the protein-coding genes for all chromosomes and plasmids for
each row in best-prokaryotes.txt,

These steps can be done at the same time on the Duke Compute Cluster (DCC).
And make sure these scripts are called from the src directory.
fetch-gbk-annotation runs for several hours.  

First, make some output directories.  
mkdir ../results/gbk-annotation/

Then:  
sbatch --mem=16G -t 36:00:00 --wrap="python fetch-gbk-annotation.py"  

Now run the following scripts on DCC. Some run
quite quickly, so no need to submit them to a partition on DCC--
just run them in an interactive session on DCC.  

python make-chromosome-plasmid-table.py  
python make-gbk-annotation-table.py  
python count-proteins-and-replicon-lengths.py  

Then, copy the following files from the results/  
directory onto my local machine (same directory name and file structure).  

gbk-annotation-table.csv  
chromosome-plasmid-table.csv  
replicon-lengths-and-protein-counts.csv
../results/gbk-annotation

## POTENTIALLY IMPORTANT TODO FOR GhostKOALA plasmid and chromosome RESULTS:
## Annotate the LENGTH of each gene in the fasta headers, and parse
## and keep these data downstream, so that I can analyze the fraction of
## protein-coding sequence dedicated to metabolic pathways.


## Running GhostKOALA on plasmids.

We use the GhostKOALA algorithm associated with the KEGG database to annotate KEGG numbers:  
https://www.kegg.jp/ghostkoala/  
Then, to make input files for GhostKOALA, make a plasmid protein db as follows, and since
GhostKOALA requires a FASTA list of up to 500,000 sequences, we split the database into smaller input files.  

IMPORTANT NOTE: make-plasmid-protein-FASTA-db.py passes over genomes in which there is
a plasmid that is apparently larger than the chromosome. I found several genomes in
these complete genomes from RefSeq in which a 'plasmid' is actually the chromosome,
and the chromosome is most likely a plasmid. This is a small number of cases,
but this has a disproportionate effect on examining metabolic genes on plasmids.  

python make-plasmid-protein-FASTA-db.py  
mkdir ../results/plasmid-GhostKOALA-input-files/  
python make-plasmid-GhostKOALA-input-files.py  

Then, submit each of the input files in ../results/plasmid-GhostKOALA-input-files/ to the
GhostKOALA webserver at https://www.kegg.jp/ghostkoala.  

Save the output files generated by GhostKOALA in ../results/plasmid-GhostKOALA-output-files/,
with names like "plasmid_GhostKOALA_batch10_ko_results.tsv".  

Then, run concatenate-and-filter-plasmid-GhostKOALA-results.sh to generate
../results/plasmid_GhostKOALA_concatenated_ko_results.tsv, and 
successful_plasmid_GhostKOALA_concatenated_ko_results.tsv, which only contains rows that were
successfully mapped to a KEGG KO ID.  

Then run:  
python remove-chromosomes-from-plasmid-Ghost-KOALA-results.py > ../results/filtered_successful_plasmid_GhostKOALA_concatenated_ko_results.tsv  
python get_unique_KEGG_IDs.py > ../results/unique_plasmid_KEGG_IDs.tsv  

to get the union of all KEGG IDs found among the plasmid genes.  

Then, upload this file to the KEGG Mapper Reconstruct webserver at:  
https://www.kegg.jp/kegg/mapper/reconstruct.html,  

and copy-paste the set of KO IDs that map to 01100 metabolic pathways into the file:  
unique_plasmid_metabolic_KEGG_IDs.txt. This file should just have KO IDs in there.  

Then, run:  
python get-plasmid-metabolic-KOs.py
To generate ../results/plasmid-proteins-in-KEGG-metabolism.tsv.
This is a table of all KEGG metabolic genes in the plasmids.  

IMPORTANT NOTE: get-plasmid-metabolic-KOs.py passes over genomes in which there is
a plasmid that is apparently larger than the chromosome. I found several genomes in
these complete genomes from RefSeq in which a 'plasmid' is actually the chromosome,
and the chromosome is most likely a plasmid. This is a small number of cases,
but this has a disproportionate effect on examining metabolic genes on plasmids.  


## Running GhostKOALA on a set of 100 genomes (just chromosomes) for comparison.

Again, we use the GhostKOALA algorithm associated with the KEGG database to annotate KEGG numbers:  
https://www.kegg.jp/ghostkoala/  
Then, to make input files for GhostKOALA, make a small chromosome protein db as follows, and since
GhostKOALA requires a FASTA list of up to 500,000 sequences, we choose a set of 100 genomes,
distributed roughly equally over the rank-distribution of chromosome lengths, over the set of genomes
with ecological annotation.

python make-chromosome-protein-FASTA-db.py  

This makes one input file for GhostKOALA, ../results/chromosome_GhostKOALA_input.fasta
Submit this file to the GhostKOALA webserver at https://www.kegg.jp/ghostkoala.  

Save the output file generated by GhostKOALA as ../results/chromosome_GhostKOALA_ko_results.tsv.  

Then, run filter-chromosome-GhostKOALA-results.sh to generate
successful_chromosome_GhostKOALA_ko_results.tsv, which only contains rows that were
successfully mapped to a KEGG KO ID.  

Then run:
python get_unique_chromosome_KEGG_IDs.py > ../results/unique_chromosome_KEGG_IDs.tsv 

to get the union of all KEGG IDs found among the chromosome genes.  

Then, upload this file to the KEGG Mapper Reconstruct webserver at:  
https://www.kegg.jp/kegg/mapper/reconstruct.html,  

and copy-paste the set of KO IDs that map to 01100 metabolic pathways into the file:  
unique_chromosome_metabolic_KEGG_IDs.txt. This file should just have KO IDs in there.  

Then, run:  
python get-chromosome-metabolic-KOs.py

to generate ../results/chromosome-proteins-in-KEGG_metabolism.tsv.
This is a table of all KEGG metabolic genes in the sampled chromosomes.


## MOB-Typer plasmid annotation

We use MOB-SUITE to get important plasmid data on mobility: https://github.com/phac-nml/mob-suite  
To make input files for MOB-typer, do the following:  

mkdir ../results/plasmid-FASTA-references  
mkdir ../results/plasmid-MOB-typer-results  
python write-plasmid-seqs-for-MOB-typer.py  

Then install MOB-SUITE on your local supercomputer (Duke Compute Cluster)  
and run the python wrapper script run-MOB-typer.py using the following 
sbatch shell script:  

sbatch run-MOB-typer.sh  

Then, combine all the MOB-typer results into a single CSV file:
python comb_mob_genome_data.py  

Finally, Hye-in removed duplicated rows in combined_mob_typer_results.txt to generate
the final file of MOB-typer results, which is called "unique_mob_results.txt", as follows.  

These notes come from Hye-in's notebook:  

Some problem occurred:  
When I counted the number of lines of “combined_mob_typer_results.txt”, the number of lines was 60914.  
The number of data files is 41849.  
I noticed that there were many duplicate items in the “combined_mob_typer_results.txt” file.  

First identified the duplicate items from the file by running:  
$ sort combined_mob_typer_results.txt | sort | uniq -d > duplicated.txt  
Then, counted the # of rows from “duplicated.txt”:  
$ cat duplicated.txt | wc -l  
There were 19064 items in the “duplicated.txt”.  
Created a new .txt file that does not contain duplicated items  
$ sort -u combined_mob_typer_results.txt -o unique_mob_results.txt  

unique_mob_results.txt is the final file of MOB-typer results used in the Data analysis.  

## Data analysis.

Then run the R script "PCN-analysis.R" to generate figures and analyses.



