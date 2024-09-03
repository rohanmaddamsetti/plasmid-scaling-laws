## PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.

## TODO: use https://github.com/tidymodels/tidyclust to do 2-means clustering on plasmid size and copy number
## to see if we can recover two clusters corresponding to large and small plasmids.

library(tidyverse)
library(cowplot)
library(tidyclust)
library(ggExtra)
library(ggrepel)


################################################################################
## Functions and global variables.

ggplotRegression <- function(dat, xvar, yvar){
    ## code from:
    ## https://community.rstudio.com/t/annotate-ggplot2-with-regression-equation-and-r-squared/6112/7  
  fml <- paste(yvar, "~", xvar)

  fit <- lm(fml, dat)
  
  ggplot(fit$model, aes(x = .data[[xvar]], y = .data[[yvar]])) + 
      geom_point() +
      theme_classic() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


normalize.plasmid.lengths <- function(df.with.AnnotationAccessions) {
    ## Make a column representing plasmid length normalized by chromosome length in each AnnotationAccession
    ## Group the data frame by AnnotationAccession and calculate the maximum replicon length in each group
    max_replicon_lengths <- df.with.AnnotationAccessions %>%
        group_by(AnnotationAccession) %>%
        summarise(max_replicon_length = max(replicon_length))
    
    ## update the df with the maximum replicon lengths and normalized
    updated.df <- df.with.AnnotationAccessions %>%
        left_join(max_replicon_lengths, by = "AnnotationAccession") %>%
        ## Creating a new column 'normalized_replicon_length' by dividing 'replicon_length' by 'max_replicon_length'
        mutate(normalized_replicon_length = replicon_length / max_replicon_length)

    return(updated.df)
}


## calculate the number of plasmids in each ecological category.
make.plasmid.totals.col <- function(PIRA.PCN.estimates) {
    plasmid.totals <- PIRA.PCN.estimates %>%
        ## extra check to make sure no chromosomes here.
        filter(SeqType == "plasmid") %>%
        group_by(Annotation) %>%
        summarize(total_plasmids = n()) %>%
        arrange(desc(total_plasmids))
    return(plasmid.totals)
}


calc.plasmid.confints <- function(df) {
    ## See Wikipedia reference:
    ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    df %>%
        ## use the normal approximation for binomial proportion conf.ints
        mutate(se = sqrt(p*(1-p)/total_plasmids)) %>%
        ## See Wikipedia reference:
        ## https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
        mutate(Left = p - 1.96*se) %>%
        mutate(Right = p + 1.96*se) %>%
        ## truncate confidence limits to interval [0,1].
        rowwise() %>% mutate(Left = max(0, Left)) %>%
        rowwise() %>% mutate(Right = min(1, Right)) %>%
        ## Sort every table by the total number of plasmids.
        arrange(desc(total_plasmids))
}


make.highPCN.table <- function(PIRA.PCN.estimates, highPCNthreshold = 50) {
    
    ## count the number of isolates with high PCNs in each category.
    high.PCN.category.counts <- PIRA.PCN.estimates %>%
        ## extra check to make sure no chromosomes here.
        filter(SeqType == "plasmid") %>%
        filter(PIRACopyNumber > highPCNthreshold) %>%
        ## next two lines is to count isolates rather than genes
        select(AnnotationAccession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(high_PCN_plasmids = n)
    
    ## join columns to make the Table.
    Table <- make.plasmid.totals.col(PIRA.PCN.estimates) %>%
        left_join(high.PCN.category.counts) %>%
        mutate(high_PCN_plasmids =
                   replace_na(high_PCN_plasmids,0)) %>%
        mutate(p = high_PCN_plasmids/total_plasmids) %>%
        calc.plasmid.confints()
    
    return(Table)
}


make.lowPCN.table <- function(PIRA.PCN.estimates, lowPCNthreshold = 1) {
    
    ## count the number of isolates with high PCNs in each category.
    low.PCN.category.counts <- PIRA.PCN.estimates %>%
        ## extra check to make sure no chromosomes here.
        filter(SeqType == "plasmid") %>%
        filter(PIRACopyNumber < lowPCNthreshold) %>%
        ## next two lines is to count isolates rather than genes
        select(AnnotationAccession, Annotation) %>%
        distinct() %>%
        count(Annotation, sort = TRUE) %>%
        rename(low_PCN_plasmids = n)
    
    ## join columns to make the Table.
    Table <- make.plasmid.totals.col(PIRA.PCN.estimates) %>%
        left_join(low.PCN.category.counts) %>%
        mutate(low_PCN_plasmids =
                   replace_na(low_PCN_plasmids,0)) %>%
        mutate(p = low_PCN_plasmids/total_plasmids) %>%
        calc.plasmid.confints()
    
    return(Table)
}


make.confint.figure.panel <- function(Table, order.by.total.plasmids, title,
                                      no.category.label = FALSE) {    
    Fig.panel <- Table %>%
        mutate(Annotation = factor(
                   Annotation,
                   levels = rev(order.by.total.plasmids))) %>%
        ggplot(aes(y = Annotation, x = p)) +
        geom_point(size=1) +
        ylab("") +
        xlab("Proportion of Plasmids") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2)
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


## antibiotic-specific keywords.
chloramphenicol.keywords <- "chloramphenicol|Chloramphenicol"
tetracycline.keywords <- "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
MLS.keywords <- "macrolide|lincosamide|streptogramin"
multidrug.keywords <- "Multidrug resistance|multidrug resistance|antibiotic resistance"
beta.lactam.keywords <- "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
glycopeptide.keywords <- "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
polypeptide.keywords <- "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
diaminopyrimidine.keywords <- "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
sulfonamide.keywords <- "sulfonamide|Sul1|sul1|sulphonamide"
quinolone.keywords <- "quinolone|Quinolone|oxacin|qnr|Qnr"
aminoglycoside.keywords <- "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
macrolide.keywords <- "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
antimicrobial.keywords <- "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"


antibiotic.keywords <- paste(chloramphenicol.keywords, tetracycline.keywords, MLS.keywords, multidrug.keywords,
    beta.lactam.keywords, glycopeptide.keywords, polypeptide.keywords, diaminopyrimidine.keywords,
    sulfonamide.keywords, quinolone.keywords, aminoglycoside.keywords, macrolide.keywords, antimicrobial.keywords, sep="|")

## require that PCN estimates are supported by a minimum of MIN_READ_COUNT reads per replicon.
MIN_READ_COUNT <- 10000

## K-means clustering (2 clusters) for large and small plasmids.
kmeans_spec <- k_means(num_clusters = 2) %>%
  set_engine("stats")


################################################################################
## Import files for the superset of all complete genomes with plasmids.
## This is a superset of the genomes for which we have PCN information.

## annotate the genomes and plasmids.
gbk.annotation <- read.csv("../results/computationally-annotated-gbk-annotation-table.csv") %>%
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## collapse Annotations into a smaller number of categories as follows:
    ## Marine, Freshwater --> Water
    ## Sediment, Soil, Terrestrial --> Earth
    ## Plants, Agriculture --> Plants
    ## Animals --> Animals
    mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Animals", "Animals")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "blank", "Unannotated"))

## get the metadata for chromosomes and plasmids,
## generated from the table of complete genomes with plasmids.
replicon.annotation.data <- read.csv("../results/chromosome-plasmid-table.csv") %>%
    left_join(gbk.annotation) %>%
    ## Annotate the genera.
    mutate(Genus = stringr::word(Organism, 1))

## filter for the metadata for plasmids.
plasmid.annotation.data <- replicon.annotation.data %>%
    filter(SeqType == "plasmid")


## get chromosome and plasmid metadata for the genomes with linked sequencing reads,
## with PCN data. This file was generated by PCN_pipeline.py.
## This should be a subset of the data in replicon.annotation.data.
PCN.replicon.metadata <- read.csv("../results/NCBI-replicon_lengths.csv")


## This is the full file of MOBTyper results
## CRITICAL TODO: analyze all the additional correlates that may be in these data.
MOBTyper.results <- read.csv("../data/Maddamsetti2024_FileS5-MOBTyper-plasmid-annotations.csv") %>%
    select(
        sample_id, gc, rep_type.s., rep_type_accession.s., relaxase_type.s.,
        relaxase_type_accession.s., mpf_type, mpf_type_accession.s., orit_type.s.,
        orit_accession.s., predicted_mobility, primary_cluster_id, secondary_cluster_id,
        predicted_host_range_overall_name, observed_host_range_ncbi_name) %>%
    mutate(
        ## everything before the third-to-last underscore, then strip "_genomic" suffix
        AnnotationAccession = str_replace(str_extract(sample_id, "^.*(?=_[^_]*_[^_]*$)"), "_genomic", ""),
        ## everything after and including the second-to-last underscore
        SeqID = str_extract(sample_id, "[^_]+_[^_]+$")) 

## Import MOBTyper mobility results for all plasmids, for merging with the plasmid copy number data.
## This is just the mobility data from the above big table, generated by parse-MOBtyper-results.py.
plasmid.mobility.data <- read.csv("../results/mobility-results.csv")

## Get the length of each plasmid and chromosome, and number of proteins on them.
replicon.length.data <- read.csv("../results/replicon-lengths-and-protein-counts.csv")


## Get the length of each plasmid, and number of proteins on them.
plasmid.length.data <- replicon.length.data %>%
    filter(SeqType == "plasmid") %>%
    ## replace prefixes, for merging in the mobility results.
    mutate(Plasmid = str_replace(SeqID, "N[CZ]_","")) %>%
    ## join the mobility results.
    left_join(plasmid.mobility.data)

################################################################################
## KEGG metabolic pathways results.

######### POTENTIALLY IMPORTANT TODO: in the upstream pipeline for the GhostKOALA
######### annotation of metabolic genes-- INCLUDE THE LENGTH OF EACH GENE!
######### This will also me to analyze the FRACTION of sequence dedicated to metabolic proteins
######### on each chromosome and plasmid, to be consistent with the results in Figure 3.

## IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   
plasmid.proteins.in.KEGG.metabolism <- read.table("../results/plasmid-proteins-in-KEGG-metabolism.tsv", header = TRUE)

metabolic.genes.in.plasmids <- plasmid.proteins.in.KEGG.metabolism %>%
    group_by(SeqID, SeqType) %>%
    summarize(metabolic_protein_count = n()) %>%
    mutate(metabolic_protein_count = replace_na(metabolic_protein_count, 0))

## IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   
chromosome.proteins.in.KEGG.metabolism <- read.table("../results/chromosome-proteins-in-KEGG-metabolism.tsv", header = TRUE)

metabolic.genes.in.chromosomes <- chromosome.proteins.in.KEGG.metabolism %>%
    group_by(SeqID, SeqType) %>%
    summarize(metabolic_protein_count = n()) %>%
    mutate(metabolic_protein_count = replace_na(metabolic_protein_count, 0))

################################################################################
## CDS fraction results

##  get output from calculate-CDS-MGE-ARG-fractions.py.
## We want to answer a basic question that Lingchong asked:
## for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?
CDS.MGE.ARG.fraction.data <- read.csv("../results/CDS-MGE-ARG-fractions.csv") %>%
    ## make the dataframe compatible with replicon.annotation.data,
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(replicon.annotation.data) %>%
    ## add a column for nomalized plasmid lengths.
    normalize.plasmid.lengths() %>%
    ## add a column for non_MGE_count.
    mutate(non_MGE_protein_count = CDS_count - MGE_count) %>%
    ## calculate non-coding regions.
    mutate(non_CDS_length = SeqLength - CDS_length) %>%
    ## set nice units for linear-scale graphs.
    mutate(CDS_length_in_Mbp = CDS_length / 1000000) %>%
    mutate(non_CDS_length_in_Mbp = non_CDS_length / 1000000) %>%
    mutate(SeqLength_in_Mbp = SeqLength / 1000000)



## TODO WHEN RERUNNING FROM SCRATCH:
## Make sure all the genomes in ../results/gbk-annotation are consistent
## with the genomes annotated in computationally-annotated-genomes etc.
## (that is, there are no suppressed RefSeq genomes in this folder, and everything is annotated.)
## in the Annotation and SeqType columns, and rewrite upstream code to
## solve this problem.
## FOR DEBUGGING
bad.annotations.vec <- unique(filter(CDS.MGE.ARG.fraction.data, is.na(SeqType) | is.na(Annotation))$AnnotationAccession)
bad.annotations.df <- data.frame(BadAnnotationAccessions = bad.annotations.vec)
write.csv(bad.annotations.df, "../results/BAD-ANNOTATIONS-IN-CDS--MGE-ARG-FRACTIONS.csv", row.names=F, quote=F)

## get ARG copy number data-- this is only used for annotating ARGs.
kallisto.ARG.copy.number.data <- read.csv("../results/kallisto-ARG_copy_numbers.csv")


################################################################################
## Get the PIRA PCN estimates. These are the main data for this paper.
## IMPORTANT NOTE: We only have PCN estimates for ~10,000 plasmids in ~4,500 genomes,
## for which we can find linked short-read data in the NCBI Sequencing Read Archive (SRA).
## This is a subset of the genomes and plasmids considered in this paper.

## SMALL TODO: In PCN_pipeline.py, make sure the ThemistoID_right column is dropped
## and not written out to this CSV file.

## IMPORTANT: do NOT filter this for just plasmids just yet--
## we need to include chromosomes for proper comparison with breseq results.

## CRITICAL TODO: fix upstream annotation so that we don't have any
## "NA" or "blank" Annotation genomes in PIRA.estimates.
PIRA.estimates <- read.csv("../results/PIRA-PCN-estimates.csv") %>%
    rename(
        PIRACopyNumber = PIRA_CopyNumberEstimate,
        PIRAReadCount = ReadCount,
        PIRASequencingCoverage = SequencingCoverage,
        PIRALongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(PIRAReadCount > MIN_READ_COUNT) %>%
    mutate(RepliconDNAContent = replicon_length * PIRACopyNumber) %>%
    ## add taxonomic and ecological annotation.
    left_join(replicon.annotation.data) %>%
    normalize.plasmid.lengths()
## write the normalized data to disk.
write.csv(PIRA.estimates, "../results/PIRA-PCN-estimates-with-normalization.csv")

################################################################################
## Genomes for Supplementary Table 1.
## get concrete statistics for how many genome papers report plasmid copy number:

## sample 50 genomes with plasmids at random, and examine how many papers report PCN.
## pick genomes with some plasmid with PCN > 10, so that we are focusing on genomes where people
## might actually want to report PCN.
## for reproducibility, use random seed = 60, which I generated using random.org (see screenshot in ../data).

genomes.to.sample.from.PIRA.PCN.estimates <- PIRA.estimates %>%
    filter(PIRACopyNumber > 10) %>%
    select(AnnotationAccession) %>%
    distinct()

## set the random seed for reproducibility.
set.seed(60)

## make a vector of the chosen genomes
chosen.genomes <- sample(
    x = genomes.to.sample.from.PIRA.PCN.estimates$AnnotationAccession,
    replace = FALSE,
    size = 50,
    )

## turn into a data frame
chosen.genomes.df <- data.frame(AnnotationAccession = chosen.genomes)
## and write to file.
write.csv(chosen.genomes.df, "../results/Fifty-random-genomes-with-multicopy-plasmids.csv",
          row.names = FALSE, quote=FALSE)

################################################################################
## Supplementary Figure S3: Effect of PIRA compared to naive themisto estimates.

## naive calculation with themisto read counts, ignoring multireplicon reads.
naive.themisto.PCN.estimates <- read.csv("../results/naive-themisto-PCN-estimates.csv") %>%
    rename(
        ThemistoNaiveCopyNumber = CopyNumber,
        ThemistoNaiveReadCount = ReadCount,
        ThemistoNaiveSequencingCoverage = SequencingCoverage,
        ThemistoNaiveLongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(SeqType == "plasmid") %>%
    ## IMPORTANT: label plasmids with reads below the threshold.
    ## in the figure, put a different shape for points with insufficient reads.
    mutate(InsufficientReads = ifelse(ThemistoNaiveReadCount < MIN_READ_COUNT, TRUE, FALSE)) %>%
    select(AnnotationAccession, SeqID, SeqType, ThemistoNaiveCopyNumber, ThemistoNaiveReadCount, InsufficientReads)

PIRA.vs.naive.themisto.df <- naive.themisto.PCN.estimates %>%
    left_join(PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))

## 1,561 plasmids have sufficient reads with PIRA, but not with the naive themisto read mapping.
sum(PIRA.vs.naive.themisto.df$InsufficientReads == TRUE)
    

## Now make S3 Figure panel A.
S3FigA <- PIRA.vs.naive.themisto.df %>%
    ggplot(aes(
        x = log10(ThemistoNaiveCopyNumber),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN,
        shape = InsufficientReads)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic() +
    xlab("log10(Naive Themisto PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none', shape = 'none') +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(ThemistoNaiveCopyNumber), y=log10(PIRACopyNumber)),
        color="light blue",
        formula=y~x)


## S3 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S3FigB <- PIRA.vs.naive.themisto.df %>%
    filter(PIRA_low_PCN == TRUE) %>%
    ggplot(aes(
        x = log10(ThemistoNaiveCopyNumber),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN,
        shape = InsufficientReads)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("red")) +
    theme_classic() +
    xlab("log10(Naive Themisto PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none', shape = 'none')

## S3 Figure panels C and D remove points with insufficient reads.

## Now make S3 Figure panel C
S3FigC <- PIRA.vs.naive.themisto.df %>%
    filter(InsufficientReads == FALSE) %>%
    ggplot(aes(
        x = log10(ThemistoNaiveCopyNumber),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic() +
    xlab("log10(Naive Themisto PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none') +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(ThemistoNaiveCopyNumber), y=log10(PIRACopyNumber)),
        color="light blue",
        formula=y~x)


## S3 Figure panel D zooms in on the plasmids with PIRA PCN < 0.8.
S3FigD <- PIRA.vs.naive.themisto.df %>%
    filter(InsufficientReads == FALSE) %>%
    filter(PIRA_low_PCN == TRUE) %>%
    ggplot(aes(
        x = log10(ThemistoNaiveCopyNumber),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("red")) +
    theme_classic() +
    xlab("log10(Naive Themisto PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none')


## make Supplementary Figure S3.
S3Fig <- plot_grid(S3FigA, S3FigB, S3FigC, S3FigD, labels=c("A", "B", "C", "D"))
ggsave("../results/S3Fig.pdf", S3Fig, height=6, width=8)


################################################################################
## Supplementary Figure S4.
## Benchmarking of PIRA with themisto against PIRA with minimap2 alignments on 100 random genomes
## with low copy number plasmids (PCN < 0.8).

low.PCN.minimap2.estimates.df <- read.csv("../results/minimap2-PIRA-low-PCN-benchmark-estimates.csv") %>%
    rename("minimap2ReadCount" = ReadCount) %>%
    select("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length",
           "minimap2ReadCount", "minimap2_PIRA_CopyNumberEstimate") %>%
        ## IMPORTANT: only examine plasmids.
    filter(SeqType == "plasmid")

## merge with PIRA estimates
PIRA.vs.minimap2.df <- low.PCN.minimap2.estimates.df %>%
    left_join(PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))

## make S4 Figure panel A
S4FigA <- PIRA.vs.minimap2.df %>%
    ggplot(aes(
        x = log10(minimap2_PIRA_CopyNumberEstimate),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic() +
    xlab("log10(minimap2 PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none') +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(minimap2_PIRA_CopyNumberEstimate), y=log10(PIRACopyNumber)),
        color="light blue",
        formula=y~x)

## S4 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S4FigB <- PIRA.vs.minimap2.df %>%
    filter(PIRA_low_PCN == TRUE) %>%
    ggplot(aes(
        x = log10(minimap2_PIRA_CopyNumberEstimate),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("red")) +
    theme_classic() +
    xlab("log10(minimap2 PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none')

## Now make Supplementary Figure S4.
S4Fig <- plot_grid(S4FigA, S4FigB, labels=c("A", "B"))
ggsave("../results/S4Fig.pdf", S4Fig, height=4, width=8)

## make a linear model and examine it.
minimap2.PIRA.PCN.lm.model <- lm(
    formula = log10(PIRACopyNumber) ~ log10(minimap2_PIRA_CopyNumberEstimate),
    data = PIRA.vs.minimap2.df)
## look at the linear regression.
summary(minimap2.PIRA.PCN.lm.model)

## confidence intervals for parameters
minimap2.PIRA.conf.intervals <- confint(minimap2.PIRA.PCN.lm.model)
minimap2.PIRA.conf.intervals


###################################################################################
## Supplementary Figure S5.
## Benchmarking of these 100 random genomes with breseq as another gold standard control,
## This additional test makes sure these estimates are accurate,
## and not artifactual due to low sequencing coverage with minimap2 compared to breseq.
## probably due to stringent minimap2 parameters by default, for now will not explore or discuss in text.

## first get the metadata we need from the PIRA estimates.
PCN.benchmark.metadata.df <- PIRA.estimates %>%
    select("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length") %>%
    ## IMPORTANT: trim the ".1" suffixes of the SeqIDs so that we can properly merge
    ## with the TrimmedSeqID in the breseq low PCN benchmark summary data (see below).
    mutate(TrimmedSeqID = str_remove(SeqID, "\\..*$"))

## now get the breseq coverage results for the benchmarking genomes.
low.PCN.breseq.summary.df <- read.csv("../results/breseq-low-PCN-benchmark-estimates.csv") %>%
    left_join(PCN.benchmark.metadata.df) %>%
    ## remove rows with missing coverage
    filter(!is.na(mean_coverage))

## make a separate df with a column for the coverage for the longest replicon in each genome
low.PCN.breseq.chromosomal.coverage.df <- low.PCN.breseq.summary.df %>%
    group_by(AnnotationAccession) %>%
    filter(replicon_length == max(replicon_length)) %>%
    ungroup() %>%
    rename("chromosomal_coverage" = mean_coverage) %>%
    select("AnnotationAccession", "chromosomal_coverage")

## now join these data back to the original breseq summary data
low.PCN.breseq.estimate.df <- low.PCN.breseq.summary.df %>%
    left_join(low.PCN.breseq.chromosomal.coverage.df) %>%
    group_by(AnnotationAccession) %>%
    ## calculate breseq PCN estimates here.
    mutate(BreseqCopyNumberEstimate = mean_coverage/chromosomal_coverage) %>%
    ## and only consider plasmids.
    filter(SeqType == "plasmid")


## merge with the PIRA estimates, and benchmark copy number estimates.
PIRA.vs.breseq.df <- low.PCN.breseq.estimate.df %>%
    left_join(PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))


## Now make S5 Figure panel A.
S5FigA <- PIRA.vs.breseq.df %>%
    ggplot(aes(
        x = log10(BreseqCopyNumberEstimate),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("black", "red")) +
    theme_classic() +
    xlab("log10(breseq PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none') +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(BreseqCopyNumberEstimate), y=log10(PIRACopyNumber)),
        color="light blue",
        formula=y~x)


## S5 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S5FigB <- PIRA.vs.breseq.df %>%
    filter(PIRA_low_PCN == TRUE) %>%
    ggplot(aes(
        x = log10(BreseqCopyNumberEstimate),
        y = log10(PIRACopyNumber),
        color = PIRA_low_PCN)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values=c("red")) +
    theme_classic() +
    xlab("log10(breseq PCN)")  +
    ylab("log10(PIRA PCN)") +
    guides(color = 'none') 

## Now make Supplementary Figure S5.
S5Fig <- plot_grid(S5FigA, S5FigB, labels=c("A", "B"))
ggsave("../results/S5Fig.pdf", S5Fig, height=4, width=8)


## make a linear model and examine it.
breseq.PIRA.PCN.lm.model <- lm(
    formula = log10(PIRACopyNumber) ~ log10(BreseqCopyNumberEstimate),
    data = PIRA.vs.breseq.df)
## look at the linear regression.
summary(breseq.PIRA.PCN.lm.model)

## confidence intervals for parameters
breseq.PIRA.conf.intervals <- confint(breseq.PIRA.PCN.lm.model)
breseq.PIRA.conf.intervals


################################################################################
## Supplementary Figure S6:
## compare naive kallisto to naive themisto PCN estimates to show that PCN numbers by pseudoalignment
## are reproducible irrespective of the specific software implementation.

kallisto.replicon.PCN.estimates <- read.csv("../results/kallisto-replicon_copy_numbers.csv") %>%
    rename(KallistoNaiveCopyNumber = CopyNumber) %>%
    filter(SeqType == "plasmid") %>%
    left_join(PCN.replicon.metadata)

## compare naive themisto to naive kallisto estimates.
naive.themisto.vs.naive.kallisto.df <- naive.themisto.PCN.estimates %>%
    ## IMPORTANT: remove plasmids with insufficient reads.
    filter(InsufficientReads == FALSE) %>%
    left_join(kallisto.replicon.PCN.estimates) %>%
    filter(SeqType == "plasmid")

## make Supplementary Figure S6.
S6Fig <- naive.themisto.vs.naive.kallisto.df %>%
    ggplot(
        aes(x=log10(KallistoNaiveCopyNumber), y=log10(ThemistoNaiveCopyNumber))) +
    geom_point(size=1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_classic() +
    xlab("log10(Naive Kallisto PCN)")  +
    ylab("log10(Naive Themisto PCN)") +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(KallistoNaiveCopyNumber), y=log10(ThemistoNaiveCopyNumber)),
        color="light blue",
        formula=y~x)
## save Supplementary Figure S6.
S6Fig <- ggsave("../results/S6Fig.pdf", S6Fig, height=4, width=4)

################################################################################
## PLASMID BIOLOGY ANALYSIS
################################################################################
## Figure 1. Patterns in Plasmid copy number, and analysis of ecological and phylogenetic
## associations with PCN.

## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering to exclude NA Annotations.

## take the PIRA estimates, filter for plasmids,
## and annotate plasmids with ARGs, and MOB type.
PIRA.PCN.estimates <- PIRA.estimates %>%
    filter(SeqType == "plasmid") %>%
    ## add Plasmid column to merge plasmid.mobility.data,
    ## by splitting the SeqID string on the underscore and taking the second part.
    mutate(Plasmid = sapply(strsplit(SeqID, "_"), function(x) x[2])) %>%
    left_join(plasmid.mobility.data) %>%
    ## annotate by presence of ARGs
    mutate(has.ARG = ifelse(SeqID %in% kallisto.ARG.copy.number.data$SeqID, TRUE, FALSE))

## Run K-means clustering with K = 2 on the PIRA.PCN.estimates, based solely on replicon_length.
kmeans_fit <- kmeans_spec %>%
    fit(~ log10(replicon_length), data=PIRA.PCN.estimates)

cluster_assignment = extract_cluster_assignment(kmeans_fit)

## add the cluster assignments to PIRA.PCN.estimates.
PIRA.PCN.estimates$Cluster <- cluster_assignment$.cluster


## Make the basic plot for Fig1A, then add the marginal histograms.
Fig1A_base <- PIRA.PCN.estimates %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    guides(color="none") +
    geom_smooth(method = "lm", se = FALSE)

## Add the marginal histograms
Fig1A <- ggExtra::ggMarginal(Fig1A_base, margins="x")

## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.
unannotated.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    filter(is.na(Annotation) | (Annotation == "blank") | Annotation == "NA")

## plot the PIRA PCN estimates.
## Break down this result by predicted plasmid mobility.
Fig1B <- PIRA.PCN.estimates %>%
    filter(!is.na(PredictedMobility)) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    guides(color="none") +
    facet_grid(PredictedMobility ~ .)

## make Figure 1.
Fig1 <- plot_grid(Fig1A, Fig1B, labels=c('A', 'B'), ncol=1, rel_heights = c(1, 2.5))
ggsave("../results/Fig1.pdf", Fig1, height=8, width=4)


################################################################################
## Make a Supplementary Figure S7 that is the same as Figure 1,
## but plotting normalized plasmid length relative to the length of the longest
## chromosome.

## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering to exclude NA Annotations.

## CRITICAL TODO: there are still a point or two at normalized plasmid length == 1
## that look like bugs! Investigate and fix or verify!

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
S7FigA <- PIRA.PCN.estimates %>%
    ggplot(aes(
        x = log2(normalized_replicon_length),
        y = log2(PIRACopyNumber))) +
    geom_point(size=0.2, alpha=0.5) +
    scale_x_continuous(breaks = c(-1, -2, -5, -10, -12)) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    xlab("log2(Normalized length)") +
    ylab("log2(Copy number)")

## Break down this result by predicted plasmid mobility.
S7FigB <- PIRA.PCN.estimates %>%
    filter(!is.na(PredictedMobility)) %>%
    ggplot(aes(
        x = log2(normalized_replicon_length),
        y = log2(PIRACopyNumber))) +
    geom_point(size=0.2, alpha=0.5) +
    scale_x_continuous(breaks = c(-1, -2, -5, -10, -12)) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    xlab("log2(Normalized length)") +
    ylab("log2(Copy number)") +
    facet_grid(PredictedMobility ~ .)

## make Supplementary Figure S7.
S7Fig <- plot_grid(S7FigA, S7FigB, labels=c('A', 'B'),ncol=1, rel_heights = c(1, 2.5))
ggsave("../results/S7Fig.pdf", S7Fig, height=8, width=4)

################################################################################
## Supplementary Figures S8 through S11. Break down the result in Figure 1 by taxonomy
## and ecological category to show universality of the PCN vs. length anticorrelation.

## Supplementary Figure S8. 
## Break down by taxonomic group.
S8Fig <- Fig1A + facet_wrap(. ~ TaxonomicGroup) +
    geom_smooth(method = "lm", se = FALSE)
## save the plot.
ggsave("../results/S8Fig.pdf", S8Fig, height=8,width=8)


## Supplementary Figure S9
## Break down by taxonomic subgroup
S9Fig <- Fig1A + facet_wrap(. ~ TaxonomicSubgroup) +
    geom_smooth(method = "lm", se = FALSE)
## save the plot.
ggsave("../results/S9Fig.pdf", S9Fig, height=9, width=8)


## Supplementary FIgure S10
## Break down by genus.
S10Fig <- Fig1A + facet_wrap(. ~ Genus, ncol=12) +
    geom_smooth(method = "lm", se = FALSE)
## save the plot.
ggsave("../results/S10Fig.pdf", S10Fig, height=24, width=16)


## Supplementary Figure S11:
## show PCN distribution over ecology.
S11Fig <- Fig1A +
    geom_hline(yintercept=2,linetype="dashed",color="gray") +
    facet_wrap(.~Annotation) +
    geom_smooth(method = "lm", se = FALSE)
## save the plot.
ggsave("../results/S11Fig.pdf", S11Fig, height=8, width=8)


################################################################################
## Supplementary Figure S12. Inverse relationship between plasmid size and copy number
## in data from Yao et al. (2022) Supplementary Table S4
## and Bethke et al. (2023) Supplementary Table S1.

YouLab.PCN.data <- read.csv("../data/YouLab-PCN-data.csv") %>%
    mutate(replicon_length = Length_in_kbp*1000)

## Now combine PIRA PCN estimates with You lab data.
PIRA.Bethke.Yao.data <- PIRA.PCN.estimates %>%
    select(SeqID, SeqType, replicon_length, PIRACopyNumber) %>%
    mutate(Reference="NCBI") %>%
    ## for joining with Bethke and Yao estimates.
    mutate(CopyNumber = PIRACopyNumber) %>%
    mutate(Plasmid = SeqID) %>%
    full_join(YouLab.PCN.data)


## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
S12Fig <- ggplot(PIRA.Bethke.Yao.data,
                        aes(x=log10(replicon_length),
                            y=log10(CopyNumber),
                            color=Reference)) +
    geom_point(size=0.2, alpha=0.5) +
    scale_color_manual(values=c("purple", "light gray", "red")) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Length)") +
    theme(legend.position="top")
  
## save the plot.
ggsave("../results/S12Fig.pdf", S12Fig, height=4,width=5)

################################################################################
## Supplementary Figure S13.
## The PCN vs. plasmid length anticorrelation largely holds within individual genomes.

within.genome.correlation.data.df <- PIRA.PCN.estimates %>%
    group_by(AnnotationAccession) %>%
    summarize(
        replicons_within_genome = n(),
        correlation = cor(log10(replicon_length), log10(PIRACopyNumber)))

S13Fig <- within.genome.correlation.data.df %>%
    ## turn replicons_within_genome into a discrete factor for plotting.
    mutate(replicons_within_genome = as.factor(replicons_within_genome)) %>%
    ggplot(aes(x = replicons_within_genome, y = correlation)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "light gray") +
    geom_boxplot(outlier.size = 0.1) +
    facet_wrap(. ~ replicons_within_genome, scales = "free_x") +
    theme_classic() +
    theme(
    axis.text.x = element_blank(),   # Remove x-axis text labels
    axis.ticks.x = element_blank()   # Remove x-axis ticks
    ) +
    xlab("replicons within genome")

ggsave("../results/S13Fig.pdf", S13Fig, height=8, width=8)

################################################################################
## Supplementary Figure S14. let's make a histogram of PCN in these data.

S14Fig <- PIRA.PCN.estimates %>%
    ## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank") %>%
    ggplot(aes(x = log10(PIRACopyNumber))) +
    geom_histogram(bins=1000) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "light gray") +
    geom_vline(xintercept = 2, linetype = "dashed", color = "light gray") +
    facet_wrap(. ~ Annotation) +
    theme_classic()

ggsave("../results/S14Fig.pdf", S14Fig, height = 6, width = 6)

################################################################################
## Supplementary Figure S15.
## Examine the tails of the PCN distribution.
## are low PCN (PCN < 1) and high PCN (PCN > 50) plasmids associated with any ecology?
## there is an enrichment of high PCN plasmids in human-impacted environments.

## This vector is used for ordering axes in this figure and the next figure.
order.by.total.plasmids <- make.plasmid.totals.col(PIRA.PCN.estimates)$Annotation

## calculate the fraction of low PCN plasmids in each category.
## Make Z-distributed confidence intervals for the fraction of isolates with
## PCN < 1.

low.PCN.plasmids.table <- make.lowPCN.table(PIRA.PCN.estimates)
## plot the confidence intervals to see if there is any enrichment of low PCN plasmids in any ecological category.
S15FigA <- make.confint.figure.panel(low.PCN.plasmids.table, order.by.total.plasmids, "proportion of plasmids with PCN < 1")

## calculate the fraction of high PCN plasmids in each category.
## there is an enrichment of  PCN > 50 plasmids in human-impacted environments.
## Make Z-distributed confidence intervals for the fraction of isolates with
## PCN > 50.

high.PCN.plasmids.table <- make.highPCN.table(PIRA.PCN.estimates)
## plot the confidence intervals to see if there is any enrichment of high PCN plasmids in any ecological category.
S15FigB <- make.confint.figure.panel(high.PCN.plasmids.table, order.by.total.plasmids, "proportion of plasmids with PCN > 50")

S15Fig <- plot_grid(S15FigA, S15FigB, labels=c('A','B'), ncol=1)
ggsave("../results/S15Fig.pdf", S15Fig, height = 8, width = 6)

## calculate the total number of plasmids,
## and the number of plasmids with PCN < 1 and PCN > 50.
PCN.count <- nrow(PIRA.PCN.estimates) ## 10,261 plasmids here
PCN.count

## There are 2238 plasmids with PCN < 1 in these data.
## this is 22% of plasmids
low.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber < 1) %>% nrow()
low.PCN.count
low.PCN.count/PCN.count

## There are 473 plasmids with PCN > 50 in these data.
## This is 4.6% of plasmids.
high.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber > 50) %>% nrow()
high.PCN.count
high.PCN.count/PCN.count

## There are 2166 plasmids with PCN > 10 in these data
## This is 21% of plasmids.
multicopy10.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber > 10) %>% nrow()
multicopy10.PCN.count
multicopy10.PCN.count/PCN.count


################################################################################
## Supplementary Figure S16. Chromosome DNA content does not constrain Plasmid DNA content.

total.DNA.content.data <- PIRA.estimates %>%
    group_by(AnnotationAccession) %>%
    summarize(TotalDNAContent = sum(RepliconDNAContent))

chromosome.DNA.content.data <- PIRA.estimates %>%
    filter(SeqType == "chromosome") %>%
    group_by(AnnotationAccession) %>%
    summarize(ChromosomeDNAContent = sum(RepliconDNAContent))

plasmid.DNA.content.data <- PIRA.estimates %>%
    filter(SeqType == "plasmid") %>%
    group_by(AnnotationAccession) %>%
    summarize(PlasmidDNAContent = sum(RepliconDNAContent))

DNA.content.data <- total.DNA.content.data %>%
    left_join(chromosome.DNA.content.data) %>%
    left_join(plasmid.DNA.content.data) %>%
    left_join(PIRA.estimates) ## basically to get the Annotation column.


S16Fig <- DNA.content.data %>%
    ggplot(aes(
        x = log10(ChromosomeDNAContent),
        y = log10(PlasmidDNAContent))) +
    geom_point(size=0.2) +
    facet_wrap(. ~ Annotation) +
    theme_classic()
## save the plot
ggsave("../results/S16Fig.pdf", S16Fig)


################################################################################

## WORKING HERE ON FIGURE 2 and associated analyses.




#######################################################################
## Cliques of plasmids in the Acman et al. (2020) plasmid similarity network
## have similar sizes and copy numbers.

## Importantly, In part, the clique structure is determined by plasmid size, since large plasmids
## and small plasmids will have a lot of sequence that is not shared, by definition.

## Get the plasmid data with nice taxonomy and clique annotations from Acman et al. (2020).
Acman.metadata <- read.csv("../data/Acman2020-SupplementaryData.csv") %>%
    ## rename columns as need to merge with PIRA.PCN.estimates,
    ## and to resolve naming conflicts.
    rename(SeqID = Version) %>%
    rename(Acman_Accession = Accession) %>%
    rename(Acman_Plasmid = Plasmid) %>%
    rename(Acman_Organism = Organism) %>%
    rename(Acman_Species = Species) %>%
    rename(Acman_Genus = Genus) %>%
    rename(Acman_Family = Family) %>%
    rename(Acman_Order = Order) %>%
    rename(Acman_Class = Class) %>%
    rename(Acman_Phylum = Phylum) %>%
    rename(Acman_Domain = Domain) %>%
    rename(Acman_Length = Length)


## inner_join the Acman et al. (2020) metadata to the PIRA PCN estimates.
## we have 1,505 plasmids in this dataset.
Acman.typed.PIRA.plasmid.estimates <- PIRA.PCN.estimates %>%
    inner_join(Acman.metadata)

Acman.cliques.with.PIRA.plasmid.estimates <- Acman.typed.PIRA.plasmid.estimates %>%
    ## remove NA Cliques.
    filter(!is.na(Clique))

## sort plasmid cliques by replicon size.
Acman.mean.clique.sizes <- Acman.cliques.with.PIRA.plasmid.estimates %>%
    group_by(Clique) %>%
    summarize(mean_replicon_length = mean(replicon_length)) %>%
    mutate(Clique_replicon_length_rank = row_number(mean_replicon_length))

## now merge the ranks back to the original data.
Acman.cliques.with.PIRA.plasmid.estimates <- Acman.cliques.with.PIRA.plasmid.estimates %>%
    left_join(Acman.mean.clique.sizes) %>%
    ## set nice units for linear-scale graphs.
    mutate(replicon_length_in_Mbp = replicon_length / 1000000)


## Acman cliques show a very limited size distribution.
Acman.clique.size.plot <- Acman.cliques.with.PIRA.plasmid.estimates %>%
    ggplot(aes(
        x = Clique_replicon_length_rank,
        y = log10(replicon_length))) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic()

Acman.clique.PCN.plot <- Acman.cliques.with.PIRA.plasmid.estimates %>%
    ggplot(aes(
        x = Clique_replicon_length_rank,
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic()

Acman.plot <- plot_grid(Acman.clique.size.plot, Acman.clique.PCN.plot, labels=c('A','B'), nrow=2)


#######################################################################
## The same result holds for the PTUs in the Redondo-Salvo et al. (2020) paper.

## That is, PTUs  in the plasmid similarity network
## have similar sizes and copy numbers.

## Importantly, PTU structure is in part determined by plasmid size, since large plasmids
## and small plasmids will not have a lot of shared sequence, by definition.


## when I reformatted these data as CSV, I renamed the AccessionVersion column
## as "SeqID" to merge with my PIRA estimates.
RedondoSalvo.data <- read.csv("../data/RedondoSalvo2020-SupplementaryData/reformatted-SupplementaryData2.csv")

## 415 plasmids have copy numbers and host ranges.
RedondoSalvo.host.range.PIRA.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.data) %>%
    filter(Host_range != "-") %>%
    filter(!is.na(Host_range))

## This figure shows that copy number / plasmid size does not predict host range.
## there are narrow and broad host range plasmids both large and small.
## compare with the MOB-Typer result in this vein.
RedondoSalvo.host.range.plot <- RedondoSalvo.host.range.PIRA.estimates %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Host_range)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    facet_wrap(.~Host_range)

## 453 plasmids have copy numbers and assigned PTUs.
PTU.PIRA.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.data) %>%
    filter(PTU != "-") %>%
    filter(!is.na(PTU))

## sort PTUs by replicon size.
PTU.mean.sizes <- PTU.PIRA.estimates %>%
    group_by(PTU) %>%
    summarize(mean_replicon_length = mean(replicon_length)) %>%
    mutate(PTU_replicon_length_rank = row_number(mean_replicon_length))

## now merge the ranks back to the original data.
PTU.PIRA.estimates <- PTU.PIRA.estimates %>%
    left_join(PTU.mean.sizes) %>%
    ## set nice units for linear-scale graphs.
    mutate(replicon_length_in_Mbp = replicon_length / 1000000)


## Acman cliques show a very limited size distribution.
Redondo.Salvo.PTU.size.plot <- PTU.PIRA.estimates %>%
    ggplot(aes(
        x = PTU_replicon_length_rank,
        y = log10(replicon_length))) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic()

Redondo.Salvo.PTU.PCN.plot <- PTU.PIRA.estimates %>%
    ggplot(aes(
        x = PTU_replicon_length_rank,
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic()

Redondo.Salvo.plot <- plot_grid(Redondo.Salvo.PTU.size.plot, Redondo.Salvo.PTU.PCN.plot, labels=c('A','B'), nrow=2)

################################################################################
## examine PCN distribution over all the potential correlates in the MOB-typer results.
## let's plot PCN estimates across different correlates, as long as there are more than 10 data points in that group.

## This has 10,261 annotated plasmids with PCN data.
MOB.typed.PIRA.plasmid.estimates <- PIRA.PCN.estimates %>%
    left_join(MOBTyper.results)

## plot over rep_type groups.
rep_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(rep_type.s.)) %>%
    filter(rep_type.s. != '-') %>%
    count(rep_type.s.) %>%
    filter(n > 10)

SX1Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(rep_type.s. %in% rep_type_groups$rep_type.s.) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(rep_type.s. ~ .)
## save the plot
ggsave("../results/SX1Fig.pdf", SX1Fig,height=12,width=12)


## plot over relaxase_type.
relaxase_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(relaxase_type.s.)) %>%
    filter(relaxase_type.s. != '-') %>%
    count(relaxase_type.s.) %>%
    filter(n > 10)

SX2Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(relaxase_type.s. %in% relaxase_type_groups$relaxase_type.s.) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(relaxase_type.s. ~ .)
## save the plot
ggsave("../results/SX2Fig.pdf", SX2Fig,height=12,width=12)


## plot over mating-pair-formation type.
mpf_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(mpf_type)) %>%
    filter(mpf_type != '-') %>%
    count(mpf_type) %>%
    filter(n > 10)

SX3Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(mpf_type %in% mpf_type_groups$mpf_type) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(mpf_type ~ .)
## save the plot
ggsave("../results/SX3Fig.pdf", SX3Fig)


## plot over oriT types.
oriT_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(orit_type.s.)) %>%
    filter(orit_type.s. != '-') %>%
    count(orit_type.s.) %>%
    filter(n > 10)

SX4Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(orit_type.s. %in% oriT_type_groups$orit_type.s.) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(orit_type.s. ~ .)
## save the plot
ggsave("../results/SX4Fig.pdf", SX4Fig)


## plot over primary cluster type.
primary_cluster_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(primary_cluster_id)) %>%
    filter(primary_cluster_id != '-') %>%
    count(primary_cluster_id) %>%
    filter(n > 10)

SX5Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(primary_cluster_id %in% primary_cluster_type_groups$primary_cluster_id) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(primary_cluster_id ~ .)
## save the plot
ggsave("../results/SX5Fig.pdf", SX5Fig, height=12,width=12,limitsize=FALSE)


## plot over secondary cluster type.
secondary_cluster_type_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(secondary_cluster_id)) %>%
    filter(secondary_cluster_id != '-') %>%
    count(secondary_cluster_id) %>%
    filter(n > 10)

SX6Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(secondary_cluster_id %in% secondary_cluster_type_groups$secondary_cluster_id) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(secondary_cluster_id ~ .)
## save the plot
ggsave("../results/SX6Fig.pdf", SX6Fig, height=12,width=12,limitsize=FALSE)


## plot over predicted host range.
predicted_host_range_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(predicted_host_range_overall_name)) %>%
    filter(predicted_host_range_overall_name != '-') %>%
    count(predicted_host_range_overall_name) %>%
    filter(n > 10)

SX7Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(predicted_host_range_overall_name %in% predicted_host_range_groups$predicted_host_range_overall_name) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(predicted_host_range_overall_name ~ .)
## save the plot
ggsave("../results/SX7Fig.pdf", SX7Fig, height=12,width=12,limitsize=FALSE)


## plot over observed host range.
observed_host_range_groups <- MOB.typed.PIRA.plasmid.estimates %>%
    filter(!is.na(observed_host_range_ncbi_name)) %>%
    filter(observed_host_range_ncbi_name != '-') %>%
    count(observed_host_range_ncbi_name) %>%
    filter(n > 10)

SX8Fig <- MOB.typed.PIRA.plasmid.estimates %>%
    ## only plot groups with more than 10 data points.
    filter(observed_host_range_ncbi_name %in% observed_host_range_groups$observed_host_range_ncbi_name) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    facet_wrap(observed_host_range_ncbi_name ~ .)
## save the plot
ggsave("../results/SX8Fig.pdf", SX8Fig, height=12,width=12,limitsize=FALSE)

################################################################################
## IMPORTANT TODO: validate these MOB-Typer results using the Plasmid Finder results
## reported in Supplementary Table S5 of the Redondo-Salvo paper.





###################################################################################
## Figure 3. Coding Sequences (CDS) on plasmids follow an empirical scaling law.

## Main figure, all the points together.
## supplementary figure: same figure, separated by Annotation category.

Fig3A <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = SeqLength_in_Mbp,
            y = CDS_length_in_Mbp,
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("replicon length (Mbp)") +
    ylab("coding sequence length (Mbp)") +
    theme_classic() + guides(color = "none")


Fig3B <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(CDS_length),
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("log10(replicon length)") +
    ylab("log10(coding sequence length)") +
    theme_classic() + guides(color = "none")

Fig3 <- plot_grid(Fig3A, Fig3B, labels = c("A", "B"), nrow=1)

## save the plot.
ggsave("../results/Fig3.pdf", Fig3, height=3.5, width=9)


################################################################################
## Supplementary Figures S17 through S20. Break down the result in Figure 3 by taxonomy
## and ecological category to show universality of the CDS scaling relationship.

## Supplementary Figure S17.
## Break down by taxonomic group.
S17Fig <- Fig3B + facet_wrap(. ~ TaxonomicGroup)
## save the plot.
ggsave("../results/S17Fig.pdf", S17Fig, height=8,width=8)


## Supplementary Figure S18
## Break down by taxonomic subgroup
S18Fig <- Fig3B + facet_wrap(. ~ TaxonomicSubgroup)
## save the plot.
ggsave("../results/S18Fig.pdf", S18Fig, height=9, width=8)


## Supplementary FIgure S19
## Break down by genus.
S19Fig <- Fig3B + facet_wrap(. ~ Genus, ncol=50)
## save the plot.
ggsave("../results/S19Fig.pdf", S19Fig, height=50, width=50, limitsize = FALSE)


## Supplementary Figure S20:
## show generality over ecology.
S20Fig <- Fig3B + facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/S20Fig.pdf", S20Fig, height=8, width=8)


####################################
## look at non-coding fractions.
## CRITICAL TODO: WHY is this so much messier than
## the coding fraction??
## Is this real, or does this represent a bug or artifact??

CDS.fraction.data <- CDS.MGE.ARG.fraction.data <- 


S21FigA <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = SeqLength_in_Mbp,
            y = non_CDS_length_in_Mbp,
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("replicon length") +
    ylab("noncoding sequence length") +
    theme_classic() + guides(color = "none")


S21FigB <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(non_CDS_length),
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("log10(replicon length)") +
    ylab("log10(noncoding sequence length)") +
    theme_classic() + guides(color = "none")

S21Fig <- plot_grid(S21FigA, S21FigB, labels=c("A", "B"))
## save the plot.
ggsave("../results/S21Fig.pdf", S21Fig, height=3.5)


########################################################################
## Figure 4: plasmid size scales with metabolic capacity.
## These data suggest that metabolic capacity constrains plasmid size.

## CRITICAL TODO: we will need to carefully distinguish between true zeros and NA
## values created by a plasmid/chromosome not being included in these data.

metabolic.gene.plasmid.data <- plasmid.length.data %>%
    left_join(metabolic.genes.in.plasmids) %>%
    ## make the dataframe compatible with plasmid.annotation.data
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join with plasmid.annotation.data.
    left_join(plasmid.annotation.data) %>%
    ## get CDS data for each genome.
    left_join(CDS.MGE.ARG.fraction.data) %>%
    ## set NA values of metabolic_protein_count to zeros.
    mutate(metabolic_protein_count = ifelse(is.na(metabolic_protein_count), 0, metabolic_protein_count)) %>%
    mutate(metabolic_protein_fraction = metabolic_protein_count / protein_count)

## join chromosome.length.data to metabolic.genes.in.chromosomes,
## since only 100 genomes were examined.
## KEY ASSUMPTION: each chromosome has at least one metabolic gene.
metabolic.gene.chromosome.data <- metabolic.genes.in.chromosomes %>%
    left_join(replicon.length.data) %>%
    ## get CDS data for each genome.
    left_join(CDS.MGE.ARG.fraction.data) %>%
    mutate(metabolic_protein_fraction = metabolic_protein_count / protein_count)

## combine plasmid and chromosome metabolic gene data for Figure 4.
metabolic.gene.plasmid.and.chromosome.data <- full_join(
    metabolic.gene.plasmid.data, metabolic.gene.chromosome.data)

Fig4A <- metabolic.gene.plasmid.and.chromosome.data %>%
    ggplot(
        aes(
            x = SeqLength_in_Mbp,
            y = metabolic_protein_count,
            color = SeqType)) +
    geom_point(size=0.5,alpha=0.5) +
    xlab("replicon_length (Mbp)") +
    ylab("metabolic genes") +
    theme_classic() + guides(color = "none")

Fig4B <- metabolic.gene.plasmid.and.chromosome.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(metabolic_protein_count),
            color = SeqType)) +
    geom_point(size=0.5,alpha=0.5) +
    xlab("log10(replicon length)") +
    ylab("log10(metabolic genes)") +
    theme_classic() + guides(color = "none")


Fig4 <- plot_grid(Fig4A, Fig4B, labels = c("A", "B"), nrow=1)
## save the plot.
ggsave("../results/Fig4.pdf", Fig4, height=3.5, width=9)

################################################################################
## Supplementary Figures S22 through S25. Break down the result in Figure 4 by taxonomy
## and ecological category to show universality of the CDS scaling relationship.

## Supplementary Figure S22
## Break down by taxonomic group.
S22Fig <- Fig4B + facet_wrap(. ~ TaxonomicGroup)
## save the plot.
ggsave("../results/S22Fig.pdf", S22Fig, height=8,width=8)


## Supplementary Figure S23
## Break down by taxonomic subgroup
S23Fig <- Fig4B + facet_wrap(. ~ TaxonomicSubgroup)
## save the plot.
ggsave("../results/S23Fig.pdf", S23Fig, height=9, width=8)


## Supplementary FIgure S19
## Break down by genus.
S24Fig <- Fig4B + facet_wrap(. ~ Genus, ncol=50)
## save the plot.
ggsave("../results/S24Fig.pdf", S24Fig, height=50, width=50, limitsize = FALSE)


## Supplementary Figure S25:
## show generality over ecology.
S25Fig <- Fig4B + facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/S25Fig.pdf", S25Fig)


########################################################################
## examine the proportion of ARG-associated genes
## found on these plasmids across different ecological categories.

S26FigA <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = SeqLength_in_Mbp,
            y = ARG_count,
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("replicon length (Mbp)") +
    ylab("ARG count") +
    theme_classic() +
    guides(color = "none") +
    ggtitle("Antibiotic resistance genes")

S26FigB <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(ARG_count),
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.5) +
    xlab("log10(replicon length") +
    ylab("log10(ARG count)") +
    theme_classic() +
    guides(color = "none")

S26Fig <- plot_grid(S26FigA, S26FigB, labels=c('A','B'))
## save the plot.
ggsave("../results/S26Fig.pdf", S26Fig, height=3.5)

########################################################################
## examine MGE-associated genes.

S27FigA <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = SeqLength_in_Mbp,
            y = MGE_count,
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.2) +
     xlab("replicon length (Mbp)") +
    ylab("MGE count") +
    theme_classic() +
    guides(color = "none") +
    ggtitle("MGE-associated genes")

S27FigB <- CDS.MGE.ARG.fraction.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(MGE_count),
            color = SeqType)) +
    geom_point(size=0.05,alpha=0.2) +
     xlab("log10(replicon length") +
    ylab("log10(MGE count)") +
    theme_classic() +
    guides(color = "none")

S27Fig <- plot_grid(S27FigA, S27FigB, labels=c('A','B'))
## save the plot.
ggsave("../results/S27Fig.pdf", S27Fig, height=3.5)


###################################################################################
## The following are analyses which are not in the paper, but that I want to keep for now,
## in case they make it in later, or for future follow-up papers.
###################################################################################
###################################################################################
## let's compare long read PCN estimates from minimap2 and pseudoalignment for high PCN plasmids--
## This shows basically no correlation between PCN estimates! This, along with the paper
## "Recovery of small plasmid sequences via Oxford Nanopore sequencing" I found
## I found about inferring PCN with Oxford nanopore suggests that long-read sequencing
## is not appropriate for inferring PCN.

high.PCN.plasmids <- PIRA.PCN.estimates %>%
    filter(PIRACopyNumber > 100)

## merge high.PCN.plasmid.RunID.df with high.PCN.plasmids.
high.PCN.plasmid.RunID.df <- read.csv("../results/high-PCN-plasmid-RunID_table.csv")

longread.high.PCN.estimates <- read.csv("../results/longread-alignment-PCN-estimates.csv") %>%
    rename(longread.minimap2.CopyNumber = CopyNumber)

## 115 plasmids in this comparison.
long.read.alignment.versus.short.read.pseudoalignment.high.PCN.estimates <- high.PCN.plasmids %>%
    ## add RefSeq_ID column to merge RunID.df,
    ## by splitting the AnnotationAccession string on the second underscore in the string.
    mutate(RefSeq_ID = sapply(strsplit(AnnotationAccession, "_"), function(x) paste(x[1], x[2], sep="_"))) %>%
    inner_join(high.PCN.plasmid.RunID.df) %>%
    ## there are duplicate rows??? remove these.
    distinct() %>%
    inner_join(longread.high.PCN.estimates) %>%
    distinct()


high.PCN.estimate.comparison.plot <- long.read.alignment.versus.short.read.pseudoalignment.high.PCN.estimates %>%
    ggplot(aes(
        x = log10(PIRACopyNumber),
        y = log10(longread.minimap2.CopyNumber),
        color = LongReadDataType)) +
    geom_point() +
    theme_classic() +
    xlab("log10(PIRA Copy Number)")  +
    ylab("log10(Long Read Minimap2 Copy Number)") +
    theme(legend.position="top") +
    xlim(-2,3) +
    ylim(-2,3)

ggsave(
    "../results/long-read-alignment-versus-short-read-pseudoalignment-high-PCN-estimates.pdf",
    high.PCN.estimate.comparison.plot, height=5, width=5)


###################################################################################
## POTENTIAL TODO: compare a mixture model fit (two lines, with a breakpoint as an extra parameter),
## to a second-order polynomial fit.

## make a linear regression model:
## log10(Plasmid copy number) vs. log10(Plasmid length).
plasmid.lm.model <- lm(
    formula=log10(PIRACopyNumber)~log10(replicon_length),
    data=PIRA.PCN.estimates)
## look at the linear regression.
summary(plasmid.lm.model)

second.order.plasmid.lm.model <- lm(
    formula=log10(PIRACopyNumber)~poly(log10(replicon_length),2,raw=TRUE),
    data=PIRA.PCN.estimates)
## look at the second order regression.
summary(second.order.plasmid.lm.model)

## let's compare these models. The model with lower AIC is better.
AIC(plasmid.lm.model)
AIC(second.order.plasmid.lm.model)
## also compare the models using an ANOVA.
print(anova(
    plasmid.lm.model,
    second.order.plasmid.lm.model))

