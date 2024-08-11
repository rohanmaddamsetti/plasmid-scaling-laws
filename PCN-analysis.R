## NCBI-PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.

## CRITICAL TODO: recalculate CDS fraction data,
## and Investigate why CDS.fraction.data has so many NA values
## in the Annotation and SeqType columns, and rewrite upstream code to
## solve this problem.

## CRITICAL TODO: repeat the metabolic gene analysis with MGE-associated genes.
## examine the proportion of MGE-associated genes
## found on these plasmids across different ecological categories,
## as replicon length varies.


## POTENTIAL TODO: make a figure comparing the fit between PIRA and Naive Themisto to the alignment methods
## (minimap2 and breseq).

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


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
MOBTyper.results <- read.csv("../data/Maddamsetti2024_FileS5-MOBTyper-plasmid-annotations.csv")

## Import MOBTyper mobility results for all plasmids, for merging with the plasmid copy number data.
## This is just the mobility data from the above big table, generated by parse-MOBtyper-results.py.
plasmid.mobility.data <- read.csv("../results/mobility-results.csv")

## Get the length of each plasmid, and number of proteins on them.
plasmid.length.data <- read.csv("../results/replicon-lengths-and-protein-counts.csv") %>%
    filter(SeqType == "plasmid") %>%
    ## replace prefixes, for merging in the mobility results.
    mutate(Plasmid = str_replace(SeqID, "N[CZ]_","")) %>%
    ## join the mobility results.
    left_join(plasmid.mobility.data)


## IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   
plasmid.proteins.in.KEGG.metabolism <- read.table("../results/plasmid-proteins-in-KEGG-metabolism.tsv", header = TRUE)

metabolic.genes.per.plasmid <- plasmid.proteins.in.KEGG.metabolism %>%
    group_by(SeqID, SeqType) %>%
    summarize(metabolic_protein_count = n()) %>%
    mutate(metabolic_protein_count = replace_na(metabolic_protein_count, 0))


##  get output from calculate-CDS-fractions.py.
## We want to answer a basic question that Lingchong asked:
## for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?
CDS.fraction.data <- read.csv("../results/CDS-fractions.csv") %>%
    ## make the dataframe compatible with replicon.annotation.data,
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(replicon.annotation.data) %>%
    ## add a column for nomalized plasmid lengths.
    normalize.plasmid.lengths()

################################################################################
## Get the PIRA PCN estimates. These are the main data for this paper.
## IMPORTANT NOTE: We only have PCN estimates for ~10,000 plasmids in ~4,500 genomes,
## for which we can find linked short-read data in the NCBI Sequencing Read Archive (SRA).
## This is a subset of the genomes and plasmids considered in this paper.

## SMALL TODO: In PCN_pipeline.py, make sure the ThemistoID_right column is dropped
## and not written out to this CSV file.

## IMPORTANT: do NOT filter this for just plasmids just yet--
## we need to include chromosomes for proper comparison with breseq results.

PIRA.estimates <- read.csv("../results/PIRA-PCN-estimates.csv") %>%
    rename(
        PIRACopyNumber = PIRA_CopyNumberEstimate,
        PIRAReadCount = ReadCount,
        PIRASequencingCoverage = SequencingCoverage,
        PIRALongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(PIRAReadCount > MIN_READ_COUNT) %>%
    normalize.plasmid.lengths()

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

## get ARG copy number data-- this is only used for annotating ARGs.
kallisto.ARG.copy.number.data <- read.csv("../results/kallisto-ARG_copy_numbers.csv")


## take the PIRA estimates, filter for plasmids,
## and annotate plasmids with ARGs, and MOB type.
PIRA.PCN.estimates <- PIRA.estimates %>%
    filter(SeqType == "plasmid") %>%
    ## add Plasmid column to merge plasmid.mobility.data,
    ## by splitting the SeqID string on the underscore and taking the second part.
    mutate(Plasmid = sapply(strsplit(SeqID, "_"), function(x) x[2])) %>%
    left_join(plasmid.mobility.data) %>%
    ## add taxonomic and ecological annotation.
    left_join(replicon.annotation.data) %>%
    ## annotate by presence of ARGs
    mutate(has.ARG = ifelse(SeqID %in% kallisto.ARG.copy.number.data$SeqID, TRUE, FALSE))



## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.
unannotated.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    filter(is.na(Annotation) | (Annotation == "blank") | Annotation == "NA")

    

## CRITICAL TODO: POLISH THIS FIGURE.
## keep consistent color/symbol coding throughout the paper.

## plot the PIRA PCN estimates.
Fig1A <- PIRA.PCN.estimates %>%
    ## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank") %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Annotation)) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    theme(legend.position="bottom")

## Break down this result by predicted plasmid mobility.
Fig1B <- Fig1A + facet_wrap(. ~ PredictedMobility)

## Break down by taxonomic group.
Fig1C <- Fig1A + facet_wrap(. ~ TaxonomicGroup)

## Break down by taxonomic subgroup
Fig1D <- Fig1A + facet_wrap(. ~ TaxonomicSubgroup)


## CRITICAL TODO: examine PCN distribution over INC groups and MOB groups.

## Figure 1EF. PCN distribution over ecology.
## clear association between high PCN plasmids and particular ecological annotations.

Fig1E <- PIRA.PCN.estimates %>%
    ## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank") %>%
    ggplot(
        aes(
            x = log10(replicon_length),
            y = log10(PIRACopyNumber),
            color = Annotation)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_hline(yintercept=2,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    guides(color = "none") +
    facet_wrap(.~Annotation)

## let's examine the ecology of high copy number plasmids.
## very cool! we see association with high PCN with humans, human-impacted environments, and livestock.
Fig1F <- PIRA.PCN.estimates %>%
    ## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank") %>%
    ggplot(aes(
        x = log10(PIRACopyNumber),
        y = Annotation,
        color = Annotation)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Ecological Annotation") +
    xlab("log10(Plasmid copy number)") 


## make Figure 1.
#Fig1 <- plot_grid(Fig1A, Fig1B, labels=c('A', 'B', 'C','D'),nrow=2)
ggsave("../results/Fig1.pdf", Fig1A, height=4, width=4)

################################################################################
## Supplementary Figure S7. Inverse relationship between plasmid size and copy number
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
S7Fig <- ggplot(PIRA.Bethke.Yao.data,
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
ggsave("../results/S7Fig.pdf", S7Fig, height=3.75,width=5)


################################################################################
## Make a Supplementary Figure S8 that is the same as Figure 1,
## but plotting normalized plasmid length relative to the length of the longest
## chromosome.

## CRITICAL TODO: there are still a few points at normalized plasmid length == 1
## that look like bugs! Investigate and fix or verify!

## POTENTIAL TODO: move this result into the main figures.
## This is a really nice result showing how the variance in plasmid copy number decreases
## as length increases.

## CRITICAL TODO: kick the tires on this result to make sure it is true, and not an artifact of doing the math wrong.

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
S8Fig <- PIRA.PCN.estimates %>%
    filter(!is.na(Annotation)) %>%
    ggplot(
        aes(
            x = normalized_replicon_length,
            y = log10(PIRACopyNumber),
            color = Annotation)) +
    geom_point(size=0.1, alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("Normalized Plasmid length") +
    guides(color = "none")

## Break down this result by predicted plasmid mobility.
S8FigB <- S8FigA + facet_grid(. ~ PredictedMobility)

## make Supplementary Figure S8.
S8Fig <- plot_grid(S8FigA, S8FigB, labels=c('A', 'B'),nrow=2)

ggsave("../results/S8Fig.pdf", S8Fig, height=6, width=7)

################################################################################
## Analysis of ARGs and plasmid copy number.

## TODO: make a supplementary figure.

## CRITICAL TODO: specifically test whether plasmids with duplicated ARGs
## tend to have higher PCN than plasmids with singleton ARGs.

## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.

ARG.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.ARG == FALSE)

## plasmids with ARGs have lower PCN than those without ARGs.
mean(ARG.plasmid.data$PIRACopyNumber)
mean(no.ARG.plasmid.data$PIRACopyNumber)

median(ARG.plasmid.data$PIRACopyNumber)
median(no.ARG.plasmid.data$PIRACopyNumber)

################################################################################
## calculate the total number of plasmids and the number of plasmids with PCN > 100.
nrow(PIRA.PCN.estimates) ## 10,261 plasmids here
## There are 176 plasmids with PCN > 100 in these data.
PIRA.PCN.estimates %>% filter(PIRACopyNumber > 100) %>% nrow()


################################################################################
## Supplementary Figure S9. let's make a histogram of PCN in these data.

S9Fig <- PIRA.PCN.estimates %>%
    ggplot(aes(x = log10(PIRACopyNumber))) +
    geom_histogram(bins=1000) +
    theme_classic()

ggsave("../results/S9Fig.pdf", S9Fig, height = 6, width = 6)

################################################################################
## Supplementary Figure S10. examine total DNA (chromosomes and plasmids) per genome.

DNA.content.data <- PIRA.estimates %>%
    mutate(DNAContent = replicon_length * PIRACopyNumber) %>%
    group_by(AnnotationAccession) %>%
    summarize(TotalDNAContent = sum(DNAContent))

S10Fig <- DNA.content.data %>%
    ggplot(aes(x=TotalDNAContent)) +
    geom_histogram(bins=1000) +
    theme_classic() +
## IMPORTANT TODO: EXAMINE OUTLIERS IN THIS PLOT.
    xlim(0,15000000)


ggsave("../results/S10Fig.pdf", S10Fig, height = 6, width = 6)


###################################################################################
## Figure 2. Coding Sequence (CDS) analysis.

## Main figure, all the points together.
## supplementary figure: same figure, separated by Annotation category.

## CRITICAL TODO: Investigate why CDS.fraction.data has so many NA values
## in the Annotation and SeqType columns, and rewrite upstream code to
## solve this problem.

## TEMPORARY HACK:
## remove NA values. CRITICAL TODO: find and fix the causes for this.
clean.CDS.fraction.data <- CDS.fraction.data %>%
    filter(!is.na(Annotation)) %>%
    filter(!is.na(SeqType))

## FOR DEBUGGING
bad.annotations.vec <- unique(filter(CDS.fraction.data, is.na(SeqType) | is.na(Annotation))$AnnotationAccession)
bad.annotations.df <- data.frame(BadAnnotationAccessions = bad.annotations.vec)
write.csv(bad.annotations.df, "../results/BAD-ANNOTATIONS-IN-CDS-FRACTIONS.csv", row.names=F, quote=F)


Fig2A <- clean.CDS.fraction.data %>%
    ggplot(
        aes(
            x = log10(SeqLength),
            y = log10(CDSLength),
            color = SeqType)) +
    geom_point(size=0.2,alpha=0.5) +
    xlab("log10(replicon length)") +
    ylab("log10(coding sequence length)") +
    theme_classic() + guides(color = "none")

Fig2B <- clean.CDS.fraction.data %>%
    ggplot(    
        aes(
            x = log10(SeqLength),
            y = CDSFraction,
            color = SeqType)) +
    geom_point(size=0.2,alpha=0.5) +
    xlab("log10(replicon length)") +
    ylab("log10(coding sequence fraction)") +
    theme_classic() +
    guides(color = "none")

Fig2C <- clean.CDS.fraction.data %>%
    ggplot(
        aes(
            x = CDSFraction,
            fill = SeqType)) +
    geom_histogram(position = 'identity', bins=100,alpha=0.5) +
    coord_flip() +
    xlab("log10(coding sequence fraction)") +
    theme_classic() +
    guides(fill = "none")

Fig2 <- plot_grid(Fig2A, Fig2B, Fig2C, labels = c("A", "B", "C"), nrow=1)

## save the plot.
ggsave("../results/Fig2.pdf", Fig2, height=4, width=7.5)

## examine the same thing for plasmids, but normalized by chromosome length.
S11Fig <- clean.CDS.fraction.data %>%
    ggplot(
        aes(
            x = log10(normalized_replicon_length),
            y = CDSFraction,
            color = SeqType)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    guides(color = "none")

## save the plot.
ggsave("../results/S11Fig.pdf", S11Fig)


########################################################################
##  analysis of metabolic genes on plasmids.
## data structures for analysis of metabolic genes on plasmids.

## make the dataframe for the exploratory data analysis.
metabolic.gene.scatterplot.data <- plasmid.length.data %>%
    left_join(metabolic.genes.per.plasmid) %>%
    ## make the dataframe compatible with plasmid.annotation.data
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(plasmid.annotation.data) %>%
    ## set NA values of metabolic_protein_count to zeros.
    mutate(metabolic_protein_count = ifelse(is.na(metabolic_protein_count), 0, metabolic_protein_count))

## annotate big.plasmids as plasmids > 500 kB.
big.plasmid.data <- metabolic.gene.scatterplot.data %>%
    filter(replicon_length > 500000)

metabolic.gene.scatterplot.data <- metabolic.gene.scatterplot.data %>%
    mutate(big_plasmids = ifelse(SeqID %in% big.plasmid.data$SeqID, TRUE, FALSE))

metabolic.gene.log.scatterplot <- ggplot(
    data = metabolic.gene.scatterplot.data,
    aes(x = log2(replicon_length), y = log2(metabolic_protein_count), color = big_plasmids)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot.pdf", metabolic.gene.log.scatterplot)


metabolic.gene.scatterplot <- ggplot(
    data = metabolic.gene.scatterplot.data,
    aes(x = replicon_length, y = metabolic_protein_count, color = big_plasmids)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot.pdf", metabolic.gene.scatterplot)


metabolic.gene.log.scatterplot2 <- ggplotRegression(
    metabolic.gene.scatterplot.data,
    "log2(replicon_length)", "log2(metabolic_protein_count)")

## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot2.pdf", metabolic.gene.log.scatterplot2)


metabolic.gene.log.scatterplot3 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    ggplot(
    aes(x = log2(replicon_length), y = log2(metabolic_protein_count), color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot3.pdf", metabolic.gene.log.scatterplot3)

metabolic.gene.scatterplot3 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    ggplot(
    aes(x = replicon_length, y = metabolic_protein_count, color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot3.pdf", metabolic.gene.scatterplot3)


metabolic.gene.log.scatterplot4 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    ggplot(
    aes(x = log2(replicon_length), y = log2(metabolic_protein_count), color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() +
    facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot4.pdf", metabolic.gene.log.scatterplot4)

metabolic.gene.scatterplot4 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    ggplot(
    aes(x = replicon_length, y = metabolic_protein_count, color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() +
    facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot4.pdf", metabolic.gene.scatterplot4)


## Super interesting. the big plasmids basically all come from nitrogen-fixing bacteria and plant pathogens!
big.plasmid.data
write.csv(x=big.plasmid.data, file="../results/big-plasmids-threshold750proteins.csv")


## calculate slope for the correlation.
metabolic.gene.log.scatterplot3 <- ggplotRegression(
    metabolic.gene.scatterplot.data,
    "replicon_length", "metabolic_protein_count")
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot3.pdf", metabolic.gene.log.scatterplot3)

## TODO:
## calculate the slope of the regression across the different ecological categories.
## look at this distribution of slope parameters, and ask whether the slope
## use a binned average, so that each range of the data gives equal contribution,
## so that this calculation does not overly favor small plasmids with few metabolic genes.

## ALSO: repeat this analysis with MGE-associated genes. examine the proportion of MGE-associated genes
## found on these plasmids across different ecological categories.
    









###################################################################################
###################################################################################
## The following are analyses which are not in the paper, but that I want to keep for now,
## in case they make it in later.
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


################################################################################
################################################################################
################################################################################
## Let's play with Zhengqing's notion of plasmid capacity.

## summarize average plasmid length and copy number (geometric mean)
geom.mean.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    group_by(AnnotationAccession) %>%
    summarize(
        geomean_PIRACopyNumber = exp(mean(log(PIRACopyNumber))),
        geomean_replicon_length = exp(mean(log(replicon_length))))

geomean.plasmid.copy.number.plot <- ggplot(
    geom.mean.PIRA.PCN.estimates,
    aes(
        x=log10(geomean_replicon_length),
        y=log10(geomean_PIRACopyNumber))) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(geometric mean Plasmid copy number)")  +
    xlab("log10(geometric mean Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=geom.mean.PIRA.PCN.estimates,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(geomean_replicon_length),y=log10(geomean_PIRACopyNumber)),
        color="light blue",
        formula=y~x)
## save the plot.
ggsave("../results/NCBI-geometric-mean-plasmid-copy-number.pdf",
       geomean.plasmid.copy.number.plot,height=5.75,width=5.75)

geomean.plasmid.lm.model <- lm(
    formula=log10(geomean_PIRACopyNumber)~log10(geomean_replicon_length),
    data=geom.mean.PIRA.PCN.estimates)
summary(geomean.plasmid.lm.model)

testplot2 <- ggplotRegression(
    geom.mean.PIRA.PCN.estimates,
    "log10(geomean_replicon_length)",
    "log10(geomean_PIRACopyNumber)")
testplot2

## Notes from Zhengqing on Slack about his methods:
## Say one cell has 20 copies of 5kB plasmid, and 5 copies of 20 KB.
## Altogether there are 200K bps.
## One method (constant copy number) to get a coherent representation as a dot
## on the scatterplot is by defining a total copy number of 25, and weighted size
## as 200K/25=8KB — approximating the cell has 25 copies of identical
## hypothetical 8KB plasmids.
## The other method (average size) is assuming all the plasmids have a geometrical
## averaged size of 10KB, and with a hypothetical copy number of 200K/10K = 10 copies.

## This first method preserves total plasmid copy number and total plasmid DNA in bp,
## and partitions DNA equally among the "average" plasmid.

test3 <- PIRA.PCN.estimates %>%
    mutate(total.plasmid.bp = replicon_length*PIRACopyNumber) %>%
    group_by(AnnotationAccession) %>%
    summarize(
        totalCN = sum(PIRACopyNumber),
        totalBP = sum(total.plasmid.bp)) %>%
    mutate(weightedSize = totalBP/totalCN)

## write to disk
write.csv(test3, file="../results/data-for-Zhengqing-style-plot1.csv")

testplot3 <- ggplotRegression(
    test3,
    "log10(weightedSize)", "log10(totalCN)")
testplot3
## save the plot.
ggsave("../results/Zhengqing-style-plot1.pdf",
       testplot3,height=5.75,width=5.75)

testplot3.5 <- ggplotRegression(
    test3,
    "log10(weightedSize)", "log10(totalBP)")
testplot3.5

## This second method preserves total plasmid bp,
## and assigns a copy number based on the geometric mean of plasmid size.
test4 <- PIRA.PCN.estimates %>%
    mutate(total.plasmid.bp = replicon_length*PIRACopyNumber) %>%
    group_by(AnnotationAccession) %>%
    summarize(
        geomean_length = exp(mean(log(replicon_length))),
        totalPlasmidLength = sum(total.plasmid.bp)) %>%
    mutate(normalizedCN = totalPlasmidLength/geomean_length)

testplot4 <- ggplotRegression(
    test4,
    "log10(geomean_length)", "log10(normalizedCN)")
testplot4
## save the plot.
ggsave("../results/Zhengqing-style-plot2.pdf",
       testplot4,height=5.75,width=5.75)

testplot4.5 <- ggplotRegression(
    test4,
    "log10(geomean_length)", "log10(totalPlasmidLength)")
testplot4.5

################################################################################
