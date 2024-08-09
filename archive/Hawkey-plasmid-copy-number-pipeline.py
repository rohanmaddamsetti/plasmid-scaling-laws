#!/usr/bin/env python

"""
Hawkey-plasmid-copy-number-pipeline.py by Rohan Maddamsetti.

For each genome in the Hawkey et al. 2022 dataset:

1) Download sequencing Illumina reads.
2) generate a fasta file of gene sequences from the genbank annotation.
     mark each gene by chromosome or plasmid, and give the replicon an ID.
     Also annotate the product field, so that ARGs and MGE-associated genes
     can be scored.
3) Run Kallisto Index on the fasta file of gene sequences.
4) Run kallisto quant on the index and sequencing reads.
5) estimate copy number of each replicon by averaging over genes,
    and make a table of copy number for all plasmids and chromosomes.
6) estimate copy number of all ARGs in the genome, relative to chromosome.
     Specifically, divide est_counts by length for each ARG to get coverage per bp,
     and divide by chromosome coverage per bp to get 
     copy number relative to chromosome.
    Then make a table of ARG copy number relative to chromosome.
7) Make a table of chromosome and plasmid lengths, based on
    the reference genomes.

NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
I checked, and I did download the right SRA reads file. But it states at the following URL
that this RefSeq ID is suppressed. So maybe something wrong with this one?
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026154285.1/

8) For benchmarking, run breseq on each of the Hawkey genomes.
9) parse breseq output to estimate plasmid copy number.

"""

import subprocess
import os
import gzip
import re
from Bio import SeqIO
from bs4 import BeautifulSoup


def download_fastq_reads(SRA_data_download_dir, SRA_accession_list_file):
    """
    use sra-toolkit to download the SRA accession data in
    SRA_accession_list_file into SRA_data_download_dir.
    """
    sra_numbers = []
    with open(SRA_accession_list_file, "r") as SRA_acc_fh:
        for SRA_acc in SRA_acc_fh:
            sra_numbers.append(SRA_acc.strip())
            
    for sra_id in sra_numbers:
        """
        the sra_id has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        sra_dir_path_for_downloading = os.path.join(SRA_data_download_dir, sra_id)
        if os.path.exists(sra_dir_path_for_downloading): continue
        prefetch_args = ["prefetch", sra_id, "-O", sra_dir_path_for_downloading]
        print (" ".join(prefetch_args))
        subprocess.run(prefetch_args)

    ## have to change working directory for fasterq_dump
    my_cwd = os.getcwd()
    os.chdir(SRA_data_download_dir)
    for sra_id in sra_numbers:
        ## paired-end files for kallisto.
        sra_fastq_file_1 = sra_id + "_1.fastq"
        sra_fastq_file_2 = sra_id + "_2.fastq"
        if os.path.exists(sra_fastq_file_1) and os.path.exists(sra_fastq_file_2):
            continue
        else:
            print ("Generating fastq for: " + sra_id)
            ## run with 10+ threads  (default is 6 threads).
            fasterq_dump_args = ["fasterq-dump", "--threads", "10", sra_id]
            print(" ".join(fasterq_dump_args))
            subprocess.run(fasterq_dump_args)
    ## now change back to original working directory.
    os.chdir(my_cwd)
    return


def generate_fasta_reference_for_kallisto(gbk_gz_path, outfile):
    """ NOTE: this code depends on some assumptions that are true
    of the Hawkey et al. 2022 (by manual verification) that may not
    generalize for datasets uploaded by other people.
    """
    with open(outfile, "w") as outfh:
        with gzip.open(gbk_gz_path,'rt') as gbk_fh:
            SeqID = None
            SeqType = None
            for record in SeqIO.parse(gbk_fh, "genbank"):
                SeqID = record.id
                if "complete" in record.description:
                    if "plasmid" in record.description:
                        SeqType = "plasmid"
                    elif "chromosome" in record.description:
                        SeqType = "chromosome"
                    else:
                        continue
                else:
                    continue
                for feature in record.features:
                    ## only analyze protein-coding genes.
                    if feature.type != "CDS": continue
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    ## Important: for kallisto, we need to replace spaces with underscores in the product annotation field.
                    ## In addition, replace commas with semicolons for compatibility with *.csv format.
                    product = feature.qualifiers["product"][0].replace(" ","_").replace(",", ";")
                    DNAseq = feature.extract(record.seq)
                    header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"locus_tag="+locus_tag,"product="+product])
                    outfh.write(header + "\n")
                    outfh.write(str(DNAseq) + "\n")
    return


def make_Hawkey_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir):    
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        print("making: ", fasta_outfile)
        generate_fasta_reference_for_kallisto(gzpath, fasta_outfile)
    return


def make_Hawkey_kallisto_indices(kallisto_ref_dir, kallisto_index_dir):
    ref_fasta_filelist = [x for x in os.listdir(kallisto_ref_dir) if x.endswith(".fna")]
    for ref_fasta_file in ref_fasta_filelist:
        ref_fasta_path = os.path.join(kallisto_ref_dir, ref_fasta_file)
        genome_id = ref_fasta_file.split(".fna")[0]
        index_file = genome_id + ".idx"
        index_path = os.path.join(kallisto_index_dir, index_file)
        kallisto_index_args = ["kallisto", "index", "-i", index_path, ref_fasta_path]
        subprocess.run(kallisto_index_args)
    return


def run_kallisto_quant(Hawkey2022_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, results_dir):
    index_list = [x for x in os.listdir(kallisto_index_dir) if x.endswith(".idx")]
    for index_file in index_list:
        index_path = os.path.join(kallisto_index_dir, index_file)
        genome_id = index_file.split(".idx")[0]
        SRA_id = Hawkey2022_genomeID_to_SRA_ID_dict[genome_id]
        if SRA_id == "NA":
            continue
        else:
            read_file1 = SRA_id + "_1.fastq"
            read_file2 = SRA_id + "_2.fastq"
            read_path1 = os.path.join(SRA_data_dir, read_file1)
            read_path2 = os.path.join(SRA_data_dir, read_file2)
            output_path = os.path.join(results_dir, genome_id)
            ## run with 10 threads by default.
            kallisto_quant_args = ["kallisto", "quant", "-t", "10", "-i", index_path, "-o", output_path, "-b", "100", read_path1, read_path2]
            subprocess.run(kallisto_quant_args)
    return


def make_genome_to_SRA_dict(Hawkey2022_metadata_csv):
    genome_to_SRA_dict = dict()
    with open(Hawkey2022_metadata_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header.
            line = line.strip()
            NCBI_bioproject, ReferenceGenome, BioSample, SRA_Data = line.split(',')
            genome_id = ReferenceGenome.split(".gbff.gz")[0]
            genome_to_SRA_dict[genome_id] = SRA_Data
    return genome_to_SRA_dict


def parse_metadata_in_kallisto_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    locus_tag = fields[2].split("=")[-1]
    ## convert underscores back into spaces, based on the conversion made in
    ## generate_fasta_reference_for_kallisto() for compatibility with kallisto.
    ## This conversion back to spaces is important for regular expression matching
    ## with this product annotation to classify proteins.
    ## In addition, replace commas with semicolons for compatibility with *.csv format.
    product = fields[3].split("=")[-1].replace("_", " ").replace(",", ";")
    metadata_tuple = (SeqID, SeqType, locus_tag, product)
    return(metadata_tuple)


def estimate_chr_plasmid_copy_numbers(genecount_tsv_path):
    """
    Input: genecount_tsv_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
    output: copy_number_dict: keys are SeqIDs and values are a tuple:
    (seqtype, seq_est_counts, seq_gene_length, chromosome_est_counts, chromosome_gene_length, seq_coverage, chromosome_coverage, seq_copy_number)
    """

    ## check to see if genecount_tsv_path exists.
    if not os.path.exists(genecount_tsv_path):
        raise ValueError(f"The input file {genecount_tsv_path} was not found.")

    ## Sum the est_counts and the length of all measured genes for all replicons in the genome.
    genome_dict = dict()
    ## keys are SeqIDs.
    ## values are a dict: {SeqType: "chromosome", total_length: 10000, total_est_counts: 100}
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_kallisto_header(target_id)
            if SeqID in genome_dict:
                genome_dict[SeqID]["total_length"] += float(length)
                genome_dict[SeqID]["total_est_counts"] += float(est_counts)
            else: ## Initialize the dictionary.
                genome_dict[SeqID] = {"SeqType" : SeqType, "total_length" : float(length), "total_est_counts": float(est_counts)}

    """ Now, calculate the coverage for each replicon. At the same time, find the chromosome
     for normalization, and save chromosome_est_counts, chromosome_gene_length,
     chromosome_coverage, to associate these values with every replicon.

    NOTE: we initialize these values to -1 to easily catch errors downstream,
    since these values should be positive if everything goes right.
    """
    chromosome_est_counts = -1
    chromosome_gene_length = -1
    chromosome_coverage = -1

    ##keys are seq_ids, value is a tuple of (SeqType, seq_est_counts, total_length, seq_coverage).
    coverage_dict = dict()
    for SeqID, replicon_dict in genome_dict.items():
        seq_coverage = replicon_dict["total_est_counts"]/replicon_dict["total_length"]
        coverage_dict[SeqID] = (replicon_dict["SeqType"], replicon_dict["total_est_counts"], replicon_dict["total_length"], seq_coverage)
        if replicon_dict["SeqType"] == "chromosome":
            chromosome_est_counts = replicon_dict["total_est_counts"]
            chromosome_gene_length = replicon_dict["total_length"]
            chromosome_coverage = seq_coverage

    ## Raise an exception when nothing aligns to the chromosome.
    if chromosome_est_counts < 1 or chromosome_gene_length < 1:
        raise ValueError(f"WARNING: no reads pseudoaligned to chromosome in file {genecount_tsv_path}")
            
    """ now normalize by chromosome coverage to get copy number estimates,
     and associate chromosome_est_counts, chromosome_gene_length, chromosome_coverage
    with every replicon, so that the data going into the copy number calculations are transparent
    in downstream data processing."""
    copy_number_dict = dict()
    for SeqID, value_tuple in coverage_dict.items():
        seqtype, seq_est_counts, seq_gene_length, seq_coverage = value_tuple
        seq_copy_number = seq_coverage/chromosome_coverage
        copy_number_dict[SeqID] = (seqtype, seq_est_counts, seq_gene_length, chromosome_est_counts, chromosome_gene_length, seq_coverage, chromosome_coverage, seq_copy_number)
    return(copy_number_dict)


def measure_Hawkey2022_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, SeqEstCounts, SeqGeneLength, ChromosomeEstCounts, ChromosomeGeneLength, SeqCoverage, ChromosomeCoverage, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = []
    SeqTypeVec = []
    SeqEstCountsVec = []
    SeqGeneLengthVec = []
    ChromosomeEstCountsVec = []
    ChromosomeGeneLengthVec = []
    SeqCoverageVec = []
    ChromosomeCoverageVec = []
    CopyNumberVec = []
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genecount_tsv_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        ## If kallisto failed on this sample, skip.
        if not os.path.exists(genecount_tsv_path):
            print(f"The file '{genecount_tsv_path}' does not exist. Skipping...")
            continue
        try: ## Handle the case that no reads pseudoalign to the chromosome.
            copy_number_dict = estimate_chr_plasmid_copy_numbers(genecount_tsv_path)
        except ValueError as kallisto_error:
            print(kallisto_error)
            continue
        
        for SeqID, value_tuple in copy_number_dict.items():
            seqtype, seq_est_counts, seq_gene_length, chromosome_est_counts, chromosome_gene_length, seq_coverage, chromosome_coverage, seq_copy_number = value_tuple

            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            
            SeqEstCountsVec.append(seq_est_counts)
            SeqGeneLengthVec.append(seq_gene_length)

            ChromosomeEstCountsVec.append(chromosome_est_counts)
            ChromosomeGeneLengthVec.append(chromosome_gene_length)
            
            SeqCoverageVec.append(seq_coverage)
            ChromosomeCoverageVec.append(chromosome_coverage)
            CopyNumberVec.append(seq_copy_number)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(SeqEstCountsVec) == len(SeqGeneLengthVec) == len(ChromosomeEstCountsVec) == len(ChromosomeGeneLengthVec) == len(SeqCoverageVec) == len(ChromosomeCoverageVec) == len(CopyNumberVec)
    
    ## now write the copy number data to file.
    with open(copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,SeqEstCounts,SeqGeneLength,ChromosomeEstCounts,ChromosomeGeneLength,SeqCoverage,ChromosomeCoverage,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            row = ",".join([AnnotationAccessionVec[i],
                            SeqIDVec[i],
                            SeqTypeVec[i],
                            str(SeqEstCountsVec[i]),
                            str(SeqGeneLengthVec[i]),
                            str(ChromosomeEstCountsVec[i]),
                            str(ChromosomeGeneLengthVec[i]),
                            str(SeqCoverageVec[i]),
                            str(ChromosomeCoverageVec[i]),
                            str(CopyNumberVec[i])])
            outfh.write(row + "\n")
    return


################################################################################################
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


def estimate_ARG_copy_numbers(genecount_tsv_path):
    """
    Input: genecount_tsv_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
    output: coverage_dict: keys are locus_tags and values are a tuple:
    (SeqID, SeqType, product, ARG_est_counts, ARG_length, chromosome_est_counts, chromosome_gene_length, ARG_coverage, chromosome_coverage, ARG_copy_number)
    """
    
    ## check to see if genecount_tsv_path exists.
    if not os.path.exists(genecount_tsv_path):
        raise ValueError(f"The input file {genecount_tsv_path} was not found.")

    ## Sum the est_counts and the length of all measured genes on the chromosome.
    chromosome_est_counts = 0.0
    chromosomal_gene_length = 0.0

    ARG_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all ARGs.
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_kallisto_header(target_id)
            if SeqType == "chromosome":
                chromosomal_gene_length += float(length)
                chromosome_est_counts += float(est_counts)
            if isARG(product):
                ARG_coverage = float(est_counts) / float(length)
                ARG_coverage_dict[locus_tag] = (SeqID, SeqType, product, est_counts, length, ARG_coverage)
    ## NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
    ## Raise an exception when nothing aligns to the chromosome.
    if chromosome_est_counts < 1 or chromosomal_gene_length < 1:
        raise ValueError(f"WARNING: no reads pseudoaligned to chromosome in file {genecount_tsv_path}")

    chromosome_coverage = chromosome_est_counts/chromosomal_gene_length
    ## now normalize by chromosome coverage to get copy number estimates.
    ARG_copy_number_dict = dict()
    for locus_tag, value_tuple in ARG_coverage_dict.items():
        SeqID, SeqType, product, ARG_est_counts, ARG_length, ARG_coverage = value_tuple
        ARG_copy_number = ARG_coverage/chromosome_coverage
        ARG_copy_number_dict[locus_tag] = (SeqID, SeqType, product, ARG_est_counts, ARG_length, chromosome_est_counts, chromosomal_gene_length, ARG_coverage, chromosome_coverage, ARG_copy_number)
    return(ARG_copy_number_dict)


def measure_Hawkey2022_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, locus_tag, product, CopyNumber

    AnnotationAccession, SeqID, SeqType, locus_tag, product, GeneEstCounts, GeneLength, ChromosomeEstCounts, ChromosomeGeneLength, GeneCoverage, ChromosomeCoverage, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    LocusTagVec = []
    ProductVec = []
    GeneEstCountsVec = []
    GeneLengthVec = []
    ChromosomeEstCountsVec = []
    ChromosomeGeneLengthVec = []
    GeneCoverageVec = []
    ChromosomeCoverageVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genecount_tsv_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        ## If kallisto failed on this sample, skip.
        if not os.path.exists(genecount_tsv_path):
            print(f"The file '{genecount_tsv_path}' does not exist. Skipping...")
            continue
        try: ## Handle the case that no reads pseudoalign to the chromosome.
            ARG_copy_number_dict = estimate_ARG_copy_numbers(genecount_tsv_path)
        except ValueError as kallisto_error:
            print(kallisto_error)
            continue
    
        for locus_tag, value_tuple in ARG_copy_number_dict.items():
            SeqID, seqtype, product, ARG_est_counts, ARG_length, chromosome_est_counts, chromosome_gene_length, ARG_coverage, chromosome_coverage, ARG_copy_number = value_tuple
            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            LocusTagVec.append(locus_tag)
            ProductVec.append(product)
            GeneEstCountsVec.append(ARG_est_counts)
            GeneLengthVec.append(ARG_length)
            ChromosomeEstCountsVec.append(chromosome_est_counts)
            ChromosomeGeneLengthVec.append(chromosome_gene_length)
            GeneCoverageVec.append(ARG_coverage)
            ChromosomeCoverageVec.append(chromosome_coverage)
            CopyNumberVec.append(ARG_copy_number)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(LocusTagVec) == len(ProductVec) == len(GeneEstCountsVec) == len(GeneLengthVec) == len(ChromosomeEstCountsVec) == len(ChromosomeGeneLengthVec) == len(GeneCoverageVec) == len(ChromosomeCoverageVec) == len(CopyNumberVec)
    
    ## now write the ARG copy number data to file.
    with open(ARG_copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,locus_tag,product,ARGEstCounts,ARGLength,ChromosomeEstCounts,ChromosomeGeneLength,ARGCoverage,ChromosomeCoverage,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            row = ",".join([AnnotationAccessionVec[i],
                            SeqIDVec[i],
                            SeqTypeVec[i],
                            LocusTagVec[i],
                            ProductVec[i],
                            str(GeneEstCountsVec[i]),
                            str(GeneLengthVec[i]),
                            str(ChromosomeEstCountsVec[i]),
                            str(ChromosomeGeneLengthVec[i]),
                            str(GeneCoverageVec[i]),
                            str(ChromosomeCoverageVec[i]),
                            str(CopyNumberVec[i])])
            outfh.write(row + "\n")
    return


def tabulate_Hawkey2022_replicon_lengths(refgenomes_dir, replicon_length_csv_file):
    with open(replicon_length_csv_file, 'w') as outfh:
        header = "AnnotationAccession,SeqID,replicon_length\n"
        outfh.write(header)
        for gbk_gz in os.listdir(refgenomes_dir):
            if not gbk_gz.endswith(".gbff.gz"): continue
            annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]
            infile = os.path.join(refgenomes_dir, gbk_gz)
            with gzip.open(infile, "rt") as genome_fh:
                for replicon in SeqIO.parse(genome_fh, "gb"):
                    SeqID = replicon.id
                    replicon_length = str(len(replicon))
                    ## now write out the data for the replicon.
                    row = ','.join([annotation_accession, SeqID, replicon_length])
                    outfh.write(row + "\n")
    return


def run_breseq_on_DCC(Hawkey2022_metadata_csv, unzipped_refgenomes_dir, SRA_data_dir, breseq_results_dir):
    with open(Hawkey2022_metadata_csv, "r") as metadata_fh:
        for i, line in enumerate(metadata_fh):
            line = line.strip()
            if i == 0: continue ## skip the header.
            bioproject, gbk_gz, biosample, sra_id = line.split(',')
            unzipped_gbk = gbk_gz.split(".gz")[0]
            annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]

            breseq_outdir = os.path.join(breseq_results_dir, annotation_accession)
            cur_reference_gbk = os.path.join(unzipped_refgenomes_dir, unzipped_gbk)
            fastq1 = os.path.join(SRA_data_dir, sra_id+"_1.fastq")
            fastq2 = os.path.join(SRA_data_dir, sra_id+"_2.fastq")
    
            breseq_args = ["breseq", "-o", breseq_outdir, "-r", cur_reference_gbk, fastq1, fastq2]
            breseq_string = " ".join(breseq_args)
            sbatch_string = "sbatch -p scavenger -t 2:30:00 --mem=8G --wrap=\"" + breseq_string + "\""
            print(sbatch_string)
            subprocess.run(sbatch_string, shell=True)
    return


def parse_breseq_results(Hawkey2022_breseq_outdir, results_csv_path):
    with open(results_csv_path, "w") as csv_fh:
        ## print a header string.
        print("AnnotationAccession,SeqID,replicon_length,mean_coverage,description,SeqType,CopyNumber", file=csv_fh)
        ## filter on IDs starting with "GCF_"
        for genome_accession in [x for x in os.listdir(Hawkey2022_breseq_outdir) if x.startswith("GCF_")]:
            breseq_summary_path = os.path.join(Hawkey2022_breseq_outdir, genome_accession, "output", "summary.html")
            ## Read the HTML file if it exists.
            if not os.path.exists(breseq_summary_path): continue
            with open(breseq_summary_path, 'r') as summary_fh:
                html_content = summary_fh.read()

            ## Create a BeautifulSoup object
            soup = BeautifulSoup(html_content, 'html.parser')

            ## Find the table by its section heading
            table_section = soup.find('h2', string='Reference Sequence Information')
            reference_table = table_section.find_next('table')

            ## Extract the table data
            table_data = []
            for row in reference_table.find_all('tr'):
                row_data = [cell.get_text(strip=True) for cell in row.find_all('td')]
                if row_data and row_data[0] == "coverage":
                    table_data.append(row_data)

            ## Print the extracted table data
            for i, row in enumerate(table_data):
                SeqID = row[2]
                replicon_length = row[3].replace(",", "") ## remove commas in the length.
                mean_coverage = row[4]
                ## for CSV compatibility, replace commas with semicolons.
                description = row[-1].replace(",",";")
                if i == 0:
                    SeqType = "chromosome"
                    chromosome_coverage = float(mean_coverage)
                else:
                    SeqType = "plasmid"
                ## handle cases where coverage fit failed.
                if mean_coverage == "NA" or mean_coverage.startswith('['):
                    copy_number = "NA"
                else:
                    copy_number = str(float(mean_coverage)/chromosome_coverage)
                csv_string = ",".join([genome_accession, SeqID, replicon_length, mean_coverage, description, SeqType, copy_number])
                print(csv_string, file=csv_fh)
    return


def main():

    RUN_STAGE = 0

    SRA_data_dir = "../data/Hawkey2022-SRA-data/"
    SRA_accession_list_file = "../data/Hawkey2022-SRA-accessions.txt"
    
    refgenomes_dir = "../data/Hawkey2022-Hybrid-Assemblies-NCBI-BioProject-PRJNA646837/"
    kallisto_ref_dir = "../results/Hawkey2022_kallisto_references/"
    kallisto_index_dir = "../results/Hawkey2022_kallisto_indices/"
    kallisto_quant_results_dir = "../results/Hawkey2022_kallisto_quantification"

    Hawkey2022_metadata_csv = "../data/Hawkey2022-genome-metadata.csv"
    Hawkey2022_genomeID_to_SRA_ID_dict = make_genome_to_SRA_dict(Hawkey2022_metadata_csv)

    copy_number_csv_file = "../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv"
    ARG_copy_number_csv_file = "../results/Hawkey2022_ARG_copy_numbers.csv"

    replicon_length_csv_file = "../results/Hawkey2022_replicon_lengths.csv"

    ## I made the following directory by hand, by copying refgenomes_dir and unzipping the files
    ## with "gunzip *.gz". I could automate this in the future.
    unzipped_refgenomes_dir = "../results/unzipped-NCBI-BioProject-PRJNA646837/"
    Hawkey2022_breseq_outdir = "../results/Hawkey2022-breseq-runs/"
    Hawkey2022_breseq_csv = "../results/Hawkey2022-breseq-copy-number.csv"
    
    if RUN_STAGE == 0:
        pass ## do nothing. this is for debugging.
    elif RUN_STAGE == 1:
        download_fastq_reads(SRA_data_dir, SRA_accession_list_file)
    elif RUN_STAGE == 2:
        make_Hawkey_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_dir)
    elif RUN_STAGE == 3:
        make_Hawkey_kallisto_indices(kallisto_ref_dir, kallisto_index_dir)
    elif RUN_STAGE == 4:
        run_kallisto_quant(Hawkey2022_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, kallisto_quant_results_dir)
    elif RUN_STAGE == 5:
        measure_Hawkey2022_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file)
    elif RUN_STAGE == 6:
        measure_Hawkey2022_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file)
    elif RUN_STAGE == 7:
        tabulate_Hawkey2022_replicon_lengths(refgenomes_dir, replicon_length_csv_file)
    elif RUN_STAGE == 8:
        run_breseq_on_DCC(Hawkey2022_metadata_csv, unzipped_refgenomes_dir, SRA_data_dir, Hawkey2022_breseq_outdir)
    elif RUN_STAGE == 9:
        parse_breseq_results(Hawkey2022_breseq_outdir, Hawkey2022_breseq_csv)
        
    return


main()
