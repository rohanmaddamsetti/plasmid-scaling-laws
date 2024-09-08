## PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.

## CRITICAL POINT: Annotate megaplasmid/chromids to differentiate them from plasmids, based on size.
## The metabolic scaling law emerges in chromids.
## From Coluzzi et al. (2022): To avoid the misidentification of ICEs as
## conjugative plasmids in chromids or secondary chromosomes,
## we excluded from further study the 419 plasmids larger than 500 kb.

## CRTITICAL POINT: my analysis does not distinguish between chromids and megaplasmids,
## and the definition may be intrinsically blurry.
## See the very nice review on megaplasmids by James Hall and David Baltrus.

## The wikipedia page on chromids is very informative:
## https://en.wikipedia.org/wiki/Secondary_chromosome


library(segmented) ## put this first so that dplyr::select is not masked.
library(tidyverse)
library(cowplot)
library(tidyclust)
library(ggExtra)
library(ggrepel)
library(viridis)


################################################################################
## Functions and global variables.

filter.correlate.column <- function(df, correlate_column_name_string, min_group_size = 10) {
    ## this function filters data frames for groups with more than 10 data points in the column
    ## named in the string correlate_column_name_string, using tidy evaluation.
    ## I figured this out using ChatGPT and Chapter 20 of Advanced R by Hadley Wickham
    ## as a reference: https://adv-r.hadley.nz/evaluation.html#tidy-evaluation
    
    correlate_column_name <- sym(correlate_column_name_string)
    
    correlate.groups <- df %>%
        filter(!is.na(!!correlate_column_name)) %>%
        filter(!!correlate_column_name != '-') %>%
        count(!!correlate_column_name) %>%
        filter(n > min_group_size)

    df %>%
        filter(!!correlate_column_name %in% correlate.groups[[correlate_column_name_string]])
}


rank.correlate.column <- function(df, correlate_column_name_string) {
    ## this function ranks values in the correlate column by mean(replicon_length).
    ## the relevant correlate column is named in the string correlate_column_name_string,
    ## and parsed using tidy evaluation.
    correlate_column_name <- sym(correlate_column_name_string)

    correlate.column.ranks <- df %>%
        ## remove missing values (NA or '-).
        filter(!is.na(!!correlate_column_name)) %>%
        filter(!!correlate_column_name != "-") %>%
        group_by(!!correlate_column_name) %>%
        summarize(mean_replicon_length = mean(replicon_length)) %>%
        mutate(rank = row_number(mean_replicon_length))
    
    ## join the ranks with the original dataframe.
    df %>%
        inner_join(correlate.column.ranks) %>%
        ## and order the column by replicon length.
        mutate(!!correlate_column_name := fct_reorder(!!correlate_column_name, rank))
}


cluster_PIRA.PCN.estimates_by_plasmid_length <- function(PIRA.PCN.estimates) {
    ## Run K-means clustering with K = 2 on the PIRA.PCN.estimates, based solely on replicon_length.

    ## K-means clustering (2 clusters) for large and small plasmids.
    kmeans_spec <- k_means(num_clusters = 2) %>%
        set_engine("stats")
    
    kmeans_fit <- kmeans_spec %>%
        fit(~ log10(replicon_length), data=PIRA.PCN.estimates)
    
    cluster_assignment = extract_cluster_assignment(kmeans_fit)

    ## make a copy of the input dataframe
    clustered.PIRA.PCN.estimates <- PIRA.PCN.estimates
    ## add a column for the plasmid clusters
    clustered.PIRA.PCN.estimates$Cluster <- cluster_assignment$.cluster
    ## return the clustered dataframe.
    return(clustered.PIRA.PCN.estimates)
}


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


make_PCN_base_plot <- function(my.PCN.data) {
    ## Make the basic plot for Fig1A, before adding the marginal histograms,
    ## or facetting by column
    my.PCN.data %>%
        ggplot(aes(
            x = log10_replicon_length,
            y = log10_PIRACopyNumber,
            color = Cluster)) +
        geom_point(size=0.2,alpha=0.5) +
        geom_hline(yintercept=0,linetype="dashed",color="gray") +
        theme_classic() +
        scale_color_manual(values=c("#d95f02","#7570b3")) +
        xlab("log10(Length)")  +
        ylab("log10(Copy Number)") +
        guides(color="none") +
        theme(strip.background = element_blank())
}


make_normalized_PCN_base_plot <- function(my.PCN.data) {
    ## Make the basic plot for S7FigA, before adding the marginal histograms,
    ## or facetting by column
    my.PCN.data %>%
        ggplot(aes(
            x = log10_normalized_replicon_length,
            y = log10_PIRACopyNumber,
            color = Cluster)) +
        geom_point(size=0.2,alpha=0.5) +
        geom_hline(yintercept=0,linetype="dashed",color="gray") +
        theme_classic() +
        scale_color_manual(values=c("#d95f02","#7570b3")) +
        xlab("log10(Normalized length)")  +
        ylab("log10(Copy Number)") +
        guides(color="none") +
        theme(strip.background = element_blank())
}


make_CDS_scaling_base_plot <- function(CDS.MGE.ARG.fraction.data) {
    CDS.MGE.ARG.fraction.data %>%
        ggplot(
            aes(
                x = log10(SeqLength),
                y = log10(CDS_length),
                color = SeqType)) +
        geom_point(size=0.05,alpha=0.5) +
        xlab("log10(replicon length)") +
        ylab("log10(coding sequence length)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank())
}


make_metabolic_scaling_base_plot <- function(metabolic.gene.and.chromosome.data) {
    metabolic.gene.plasmid.and.chromosome.data %>%
        ggplot(
            aes(
                x = log10(SeqLength),
                y = log10(metabolic_protein_count),
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.5) +
        xlab("log10(replicon length)") +
        ylab("log10(metabolic genes)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank())
}


## require that PCN estimates are supported by a minimum of MIN_READ_COUNT reads per replicon.
MIN_READ_COUNT <- 10000

## To avoid the misidentification of ICEs as conjugative plasmids in chromids or secondary chromosomes,
## specifically label megaplasmids or chromids larger than 500 kb (Cite Coluzzi et al. 2022 for this practice).
PLASMID_LENGTH_THRESHOLD <- 500000


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
    normalize.plasmid.lengths()


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


## make a linear model comparing breseq to PIRA estimates and examine it.
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
## take the PIRA estimates, filter for plasmids,
## and annotate plasmids with ARGs, and MOB type.
PIRA.PCN.estimates <- PIRA.estimates %>%
    ## IMPORTANT: this removes megaplasmids and chromids from this analysis!! BE AWARE!
    filter(SeqType == "plasmid") %>%
    ## add Plasmid column to merge plasmid.mobility.data,
    ## by splitting the SeqID string on the underscore and taking the second part.
    mutate(Plasmid = sapply(strsplit(SeqID, "_"), function(x) x[2])) %>%
    left_join(plasmid.mobility.data) %>%
    ## The next two lines are to get segmented regression working, using library(segmented).
    mutate(log10_replicon_length = log10(replicon_length)) %>%
    mutate(log10_PIRACopyNumber = log10(PIRACopyNumber)) %>%
    ## this next line is to check out the segmented regression
    ## with the normalized replicon lengths.
    mutate(log10_normalized_replicon_length = log10(normalized_replicon_length)) %>%
    ## create a column indicating how plasmids cluster by length.
    cluster_PIRA.PCN.estimates_by_plasmid_length()

    
## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering to exclude NA Annotations.
## CRITICAL TODO: FIGURE OUT WHY ~1000 GENOMES ARE NOT ANNOTATED RIGHT.
unannotated.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    filter(is.na(Annotation) | (Annotation == "blank") | Annotation == "NA")

## CRITICAL TODO: there are still a point or two at normalized plasmid length == 1
## that look like bugs! Investigate and fix or verify!
## CHECK WHAT IS GOING ON!
potential.buggy.normalized.plasmids <- PIRA.PCN.estimates %>%
    filter(normalized_replicon_length == 1)


################################################################################
## The data is best fit by piecewise regression, revealing a scaling law between length and copy number.
## we use library(segmented) for the piecewise regression (segmented regression).

## first make a linear fit model.
PCN.lm.model <- lm(log10_PIRACopyNumber ~ log10_replicon_length, data=PIRA.PCN.estimates)
## look at the linear regression.
summary(PCN.lm.model)

#fit piecewise regression model to original model, estimating a breakpoint.
segmented.PCN.model <- segmented(
    PCN.lm.model,
    seg.Z = ~log10_replicon_length,
    psi = list(log10_replicon_length = 5))

## the breakpoint is at 4.774. 10^4.774 = 59,429bp.
## So plasmids above ~60Kbp have a flatter slope.
summary(segmented.PCN.model)

## save the segmented regression fit as a dataframe.
segmented.fit.df = data.frame(
    log10_replicon_length = PIRA.PCN.estimates$log10_replicon_length,
    log10_PIRACopyNumber = broken.line(segmented.PCN.model)$fit)

## Compare the segmented PCN model fit (two lines, with a breakpoint as an extra parameter),
## to a second-order polynomial fit. This analysis shows that the segmented PCN model is a better fit
## than the second-order polynomial.

## second-order polynomial fit.
second.order.PCN.lm.model <- lm(
    formula=log10_PIRACopyNumber ~ poly(log10_replicon_length,2,raw=TRUE),
    data=PIRA.PCN.estimates)

## let's compare these models. The model with lower AIC is better.
## The segmented model has the best fit.
AIC(PCN.lm.model)
AIC(second.order.PCN.lm.model)
AIC(segmented.PCN.model)

## We will add a segmented regression line to Fig1CD.
## first make a linear fit model with the normalized replicon length.
normalized.PCN.lm.model <- lm(log10_PIRACopyNumber ~ log10_normalized_replicon_length, data=PIRA.PCN.estimates)

#fit piecewise regression model based on the normalized.PCN.lm.model.
segmented.normalized.PCN.model <- segmented(
    normalized.PCN.lm.model,
    seg.Z = ~log10_normalized_replicon_length,
    psi = list(log10_normalized_replicon_length = -1.5))

## the breakpoint is at -1.762. 10^-1.762 = 1.73% of the length of the chromosome.
summary(segmented.normalized.PCN.model)

## save the segmented regression fit as a dataframe.
normalized.segmented.fit.df = data.frame(
    log10_normalized_replicon_length = PIRA.PCN.estimates$log10_normalized_replicon_length,
    log10_PIRACopyNumber = broken.line(segmented.normalized.PCN.model)$fit)

## compare to a second-order polynomial fit.
second.order.normalized.PCN.lm.model <- lm(
    formula=log10_PIRACopyNumber ~ poly(log10_normalized_replicon_length,2,raw=TRUE),
    data=PIRA.PCN.estimates)


## let's compare these models. The model with lower AIC is better.
## Again, the normalized segmented model has the best fit.
AIC(normalized.PCN.lm.model)
AIC(second.order.normalized.PCN.lm.model)
AIC(segmented.normalized.PCN.model)


################################################################################
## Figure 1. Anticorrelation between plasmid length and copy number.

## For panels A and B, don't draw the segmented regression.
Fig1A_base <- make_PCN_base_plot(PIRA.PCN.estimates) 

## Add the marginal histograms
Fig1A <- ggExtra::ggMarginal(Fig1A_base, margins="x") 

## Break down this result by predicted plasmid mobility.
Fig1B <- PIRA.PCN.estimates %>%
    filter.correlate.column("PredictedMobility") %>%
    make_PCN_base_plot() +
    facet_wrap(.~PredictedMobility, nrow=3) +
    theme(strip.background = element_blank())

## make Figure 1AB
Fig1AB_base <- plot_grid(Fig1A, Fig1B, labels=c('A', 'B'), ncol=1, rel_heights = c(1, 2.5))
Fig1AB_title <- ggdraw() + draw_label("Plasmid length", fontface='bold')
Fig1AB <- plot_grid(Fig1AB_title, Fig1AB_base, ncol=1, rel_heights=c(0.1, 1))

## Figure 1CD.
##  plot normalized plasmid length relative to the length of the longest
## chromosome in the cell.
## This figure show that the scaling law holds well, until plasmids reach ~1.73% of the main chromosome in length.
## then, copy number is roughly flat.

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
Fig1C_base <- PIRA.PCN.estimates %>%
    make_normalized_PCN_base_plot() +
    ## draw the segmented regression.
    geom_line(data = normalized.segmented.fit.df, color = 'blue')

## Add the marginal histogram
Fig1C <- ggExtra::ggMarginal(Fig1C_base, margins="x")

## Break down this result by predicted plasmid mobility.
Fig1D <- PIRA.PCN.estimates %>%
    filter.correlate.column("PredictedMobility") %>%
    make_normalized_PCN_base_plot() +
    ## draw the segmented regression.
    geom_line(data = normalized.segmented.fit.df, color = 'blue') +
    facet_wrap(.~PredictedMobility, nrow=3)

## make Fig1CD
Fig1CD_base <- plot_grid(Fig1C, Fig1D, labels=c('C', 'D'),ncol=1, rel_heights = c(1, 2.5))
Fig1CD_title <- ggdraw() + draw_label("Plasmid length normalized by chromosome", fontface='bold')
Fig1CD <- plot_grid(Fig1CD_title, Fig1CD_base, ncol=1, rel_heights=c(0.1, 1))

## make Figure 1.
Fig1 <- plot_grid(Fig1AB, Fig1CD, ncol=2)
ggsave("../results/Fig1.pdf", Fig1, height=9, width=8)

## calculate basic statistics about the clusters of small and large plasmids.

small.plasmids <- PIRA.PCN.estimates %>%
    filter(Cluster == "Cluster_2")
## mean length of small plasmids is 6212 bp.
mean(small.plasmids$replicon_length)
## mean PCN of small plasmids is 28.9.
mean(small.plasmids$PIRACopyNumber)


large.plasmids <- PIRA.PCN.estimates %>%
    filter(Cluster == "Cluster_1")
## mean length of large plasmids is 144070 bp.
mean(large.plasmids$replicon_length)
## mean PCN of large plasmids is 4.02.
mean(large.plasmids$PIRACopyNumber)

## examine the tail of very large plasmids that are longer than 500Kbp.
very.large.plasmids <- large.plasmids %>%
    filter(replicon_length > 500000)

very.large.plasmids.by.mobility <- very.large.plasmids %>%
    count(PredictedMobility)
## Very large plasmids: these are chromids.
##       conjugative  29
##       mobilizable  29
##   non-mobilizable 101
##              <NA>  32


################################################################################
## Supplementary Figures S7 through S10. Break down the result in Figure 1 by taxonomy
## and ecological category to show universality of the PCN vs. length anticorrelation.

## Supplementary Figure S7. 
## Break down by taxonomic group.
S7Fig <- PIRA.PCN.estimates %>%
    filter.correlate.column("TaxonomicGroup") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ TaxonomicGroup) +
    theme(strip.background = element_blank())
## save the plot.
ggsave("../results/S7Fig.pdf", S7Fig, height=5,width=5)


## Supplementary Figure S8
## Break down by taxonomic subgroup
S8Fig <- PIRA.PCN.estimates %>%
    filter.correlate.column("TaxonomicSubgroup") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup, ncol=3) +
    theme(strip.background = element_blank())
## save the plot.
ggsave("../results/S8Fig.pdf", S8Fig, width=8)


## Supplementary FIgure S9
## Break down by genus.
S9Fig <- PIRA.PCN.estimates %>%
    filter.correlate.column("Genus") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ Genus, ncol = 6) +
    theme(strip.background = element_blank())
## save the plot.
ggsave("../results/S9Fig.pdf", S9Fig, height=11, width=8)


## Supplementary Figure S10:
## show PCN distribution over ecology.
S10Fig <- PIRA.PCN.estimates %>%
    filter.correlate.column("Annotation") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ Annotation) +
    theme(strip.background = element_blank()) +
    geom_hline(yintercept=2,linetype="dashed",color="gray")
## save the plot.
ggsave("../results/S10Fig.pdf", S10Fig, height=8, width=8)


################################################################################
## Supplementary Figure S11. Inverse relationship between plasmid size and copy number
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
S11Fig <- ggplot(PIRA.Bethke.Yao.data,
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
ggsave("../results/S11Fig.pdf", S11Fig, height=4,width=5)

################################################################################
## Supplementary Figure S12.
## The PCN vs. plasmid length anticorrelation holds within individual genomes.

within.genome.correlation.data.df <- PIRA.PCN.estimates %>%
    group_by(AnnotationAccession) %>%
    summarize(
        replicons_within_genome = n(),
        correlation = cor(log10(replicon_length), log10(PIRACopyNumber)))

S12Fig <- within.genome.correlation.data.df %>%
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

ggsave("../results/S12Fig.pdf", S12Fig, height=8, width=8)


################################################################################
## Figure 2. Small multicopy plasmid almost always coexist with large low-copy plasmids.
## look at the joint distribution of largest plasmid per cell and smallest plasmid per cell.
## Does this give any insight into how small plasmids persist?
## there is data suggesting that plasmids tend to live as "packs" in cells.
## This is really quite an interesting result.
## relevant recent paper: "A cryptic plasmid is among the most numerous genetic elements in the human gut"
## in Cell.


## Group the data frame by AnnotationAccession and calculate the maximum replicon length in each group
max.plasmid.lengths <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession) %>%
    summarise(max_replicon_length = max(replicon_length))

min.plasmid.lengths <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession) %>%
    summarise(min_replicon_length = min(replicon_length))

max.PCNs <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession) %>%
    summarise(max_PCN = max(PIRACopyNumber))

max.and.min.plasmid.lengths <- max.plasmid.lengths %>%
    inner_join(min.plasmid.lengths) %>%
    inner_join(max.PCNs)

## 10261 plasmids
PIRA.PCN.estimates %>% nrow()

## these 10,261 plasmids are found in 3,458 genomes.
PIRA.PCN.estimates %>%
    count(AnnotationAccession) %>%
    nrow()

##3458 genomes containing plasmids
nrow(max.and.min.plasmid.lengths)

## only 1022 of these have plasmids found by themselves: 1022/3458 = 29.5%
max.and.min.plasmid.lengths %>%
    filter(max_replicon_length == min_replicon_length) %>%
    nrow()


max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids <- max.and.min.plasmid.lengths %>%
    filter(max_PCN > 10)

##1081 genomes containing multicopy plasmids here
nrow(max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids)

## only 33 are found by themselves: 33/1081 = 2.2%
max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids %>%
    filter(max_replicon_length == min_replicon_length) %>%
    nrow()

## This is obviously statistically significant.
binom.test(x=33,n=1081, p = (1022/3458))

Fig2 <- max.and.min.plasmid.lengths %>%
    ## CRITICAL TODO: FIGURE OUT WHY THIS OUTLIER IS IN THESE DATA!
    filter(max_PCN < 1000) %>%
        ggplot(aes(
        x = log10(max_replicon_length),
        y = log10(min_replicon_length),
        color = log10(max_PCN))) +
    geom_point(size=0.5) +
    theme_classic() +
    scale_color_viridis(option="magma") +
    geom_abline(linetype="dashed", color="light gray",size=0.2) +
    xlab("log10(length of largest plasmid)") +
    ylab("log10(length of smallest plasmid)") +
    labs("log10(maximum plasmid copy number)") +
    guides(color = "none") +
    ggtitle("All genomes containing plasmids")

## save the plot.
ggsave("../results/Fig2.pdf", Fig2, height=3.5, width=3.5)


################################################################################
## Examination of Plasmid length and copy number across genetic correlates
## from plasmid typing metadata.

## Get the MOB-Typer results that Hye-in generated.
## This has 10,261 annotated plasmids with PCN data.
MOB.typed.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    left_join(MOBTyper.results)

## Now import all the plasmid metadata from other papers.

## Get the plasmid data with nice taxonomy and clique annotations from Acman et al. (2020).
Acman.PTU.data <- read.csv("../data/Acman2020-SupplementaryData.csv") %>%
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


## Get the plasmid taxonomy unit (PTU) and host range annotations from Redondo-Salvo et al. (2020).
## when I reformatted these data as CSV, I renamed the AccessionVersion column
## as "SeqID" to merge with my PIRA estimates.
RedondoSalvo.PTU.data <- read.csv(
    "../data/RedondoSalvo2020-SupplementaryData/reformatted-SupplementaryData2.csv")


## get the Plasmid Finder results reported in Supplementary Table S5 of Redondo-Salvo et al. (2020).
## I renamed the AccessionVersion column to SeqID by hand when reformatting, for the join to work.
PIRA.PCN.for.RedondoSalvo2020.plasmid.metadata <- read.csv(
    "../data/RedondoSalvo2020-SupplementaryData/reformatted_SupplementaryData5.csv") %>%
    select(SeqID, MOB, PFinder_80, PFinder_95) %>%
    right_join(PIRA.PCN.estimates)


## Analyze PCN in context of the correlates in the Coluzzi et al. (2022) paper.
## "Evolution of plasmid mobility: origin and fate of conjugative and nonconjugative plasmids"
## When reformatting, I manually renamed the Acc_No_NCBI to SeqID for the join,
## and renamed Plasmid to PlasmidType to prevent a column name collision in the join.
PIRA.PCN.for.Coluzzi2022.plasmid.metadata <- read.csv(
    "../data/Coluzzi2022-SupplementaryData/reformatted_SupplementaryTable.csv") %>%
    select(SeqID, Type, MOB_Class, MPF_prots, Mobilisation_type, PlasmidType, PTU) %>%
    right_join(PIRA.PCN.estimates)


## Analyze PCN in context of the correlates in the Ares-Arroyo et al. (2023) paper.
## "Origins of transfer establish networks of functional dependencies for plasmid transfer by conjugation"
## I renamed the Plasmid column to SeqID by hand when reformatting, for the join to work.
## I have PCN for 917 plasmids in these data.
PIRA.PCN.for.AresArroyo2023.data <- read.csv(
    "../data/Ares-Arroyo2023-SupplementaryData/reformatted_Table_S1.csv") %>%
    inner_join(PIRA.PCN.estimates)


########################################
## Plasmid length and copy number are conserved within plasmid taxonomic groups.

## Supplementary Figure S13AB
## Cliques of plasmids in the Acman et al. (2020) plasmid similarity network
## have similar sizes and copy numbers.

## Importantly, In part, the clique structure is determined by plasmid size, since large plasmids
## and small plasmids will have a lot of sequence that is not shared, by definition.


## inner_join the Acman et al. (2020) metadata to the PIRA PCN estimates.
## we have 1,505 plasmids in this dataset.
Acman.cliques.with.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    inner_join(Acman.PTU.data) %>%
    rank.correlate.column("Clique")

## Acman cliques show a very limited size distribution.
Acman.clique.size.plot <- Acman.cliques.with.PIRA.PCN.estimates %>%
    ggplot(aes(
        x = rank,
        y = log10(replicon_length),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("Cliques ranked by length")  +
    ylab("log10(Length)") +
    guides(color = "none") +
    ggtitle("Plasmid cliques in Acman et al. (2020)")


Acman.clique.PCN.plot <- Acman.cliques.with.PIRA.PCN.estimates %>%
    ggplot(aes(
        x = rank,
        y = log10(PIRACopyNumber),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("Cliques ranked by length")  +
    ylab("log10(Copy Number)") +
    guides(color = "none")

S13FigAB <- plot_grid(Acman.clique.size.plot, Acman.clique.PCN.plot, labels=c('A','B'), nrow=2)


## Supplementary Figure S13CD.
## The same result holds for the PTUs in the Redondo-Salvo et al. (2020) paper.

## That is, PTUs  in the plasmid similarity network
## have similar sizes and copy numbers.

## Importantly, PTU structure is in part determined by plasmid size, since large plasmids
## and small plasmids will not have a lot of shared sequence, by definition.

## 453 plasmids have copy numbers and assigned PTUs.
PTU.PIRA.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.PTU.data) %>%
    rank.correlate.column("PTU")

## Redondo-Salvo cliques show a very limited size distribution.
Redondo.Salvo.PTU.size.plot <- PTU.PIRA.estimates %>%
    ggplot(aes(
        x = rank,
        y = log10(replicon_length),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("PTUs ranked by length")  +
    ylab("log10(Length)") +
    guides(color = "none") +
    ggtitle("PTUs in Redondo-Salvo et al. (2020)")


Redondo.Salvo.PTU.PCN.plot <- PTU.PIRA.estimates %>%
    ggplot(aes(
        x = rank,
        y = log10(PIRACopyNumber),
    color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("PTUs ranked by length")  +
    ylab("log10(Copy Number)") +
    guides(color = "none")

S13FigCD <- plot_grid(Redondo.Salvo.PTU.size.plot, Redondo.Salvo.PTU.PCN.plot, labels=c('C','D'), nrow=2)

S13Fig <- plot_grid(S13FigAB, S13FigCD, nrow=1)
ggsave("../results/S13Fig.pdf", S13Fig, height=5, width=7)

#########################################################
## Supplementary Figure S14. MOB-Typer PTU and Rep type analysis.

## First plot over primary cluster type.
## from $ mob_cluster -h:
## the mash distance for assigning the primary cluster ID is 0.06 by default.

MOB.typed.PIRA.clusters <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("primary_cluster_id")
    
## Make the plots.
MOB.Typer.PTU.size.plot <- MOB.typed.PIRA.clusters %>%
    ggplot(aes(
        x = rank,
        y = log10(replicon_length),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("PTUs ranked by length")  +
    ylab("log10(Length)") +
    guides(color = "none") +
    ggtitle("MOB-Cluster Mash Distance < 0.06")


MOB.Typer.PTU.PCN.plot <- MOB.typed.PIRA.clusters %>%
    ggplot(aes(
        x = rank,
        y = log10(PIRACopyNumber),
    color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("PTUs ranked by length")  +
    ylab("log10(Copy Number)") +
    guides(color = "none")

S14FigAB <- plot_grid(MOB.Typer.PTU.size.plot, MOB.Typer.PTU.PCN.plot, labels=c('A','B'), nrow=2)

## The same thing holds for rep protein typing.
MOB.typed.PIRA.reptypes <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("rep_type.s.")

## MOB-typer Rep protein type classes by length
MOB.Typer.reptype.size.plot <- MOB.typed.PIRA.reptypes %>%
    ggplot(aes(
        x = rank,
        y = log10(replicon_length),
    color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("Rep types ranked by length")  +
    ylab("log10(Length)") +
    guides(color = "none") +
    ggtitle("MOB-Typer Rep types")


## MOB-typer Rep protein type classes by PCN
MOB.Typer.reptype.PCN.plot <- MOB.typed.PIRA.reptypes %>%
    ggplot(aes(
        x = rank,
        y = log10(PIRACopyNumber),
        color = Cluster)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    scale_color_manual(values=c("#d95f02","#7570b3")) +
    xlab("Rep types ranked by length")  +
    ylab("log10(Copy Number)") +
    guides(color = "none")

S14FigCD <- plot_grid(MOB.Typer.reptype.size.plot, MOB.Typer.reptype.PCN.plot, labels=c('C','D'), nrow=2)

S14Fig <- plot_grid(S14FigAB, S14FigCD, ncol=2)
## save the plot
ggsave("../results/S14Fig.pdf", S14Fig)


############################
## Supplementary Figure S15.
## plot over Rep protein types in Ares-Arroyo et al. (2023).

S15Fig <- PIRA.PCN.for.AresArroyo2023.data %>%
    rank.correlate.column("Replicase") %>%
    ## only plot groups with more than 10 data points.
    filter.correlate.column("Replicase") %>%
    make_PCN_base_plot() +
    facet_wrap(Replicase ~ .)
## save the plot
ggsave("../results/S15Fig.pdf", S15Fig,height=8,width=8)

########################################
## Supplementary Figure S16. mobility group (relaxase) type analysis.

## plot over relaxase_type.
S16Fig <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("relaxase_type.s.") %>%
    filter.correlate.column("relaxase_type.s.") %>%
    make_PCN_base_plot() +
    facet_wrap(relaxase_type.s. ~ .) +
    ggtitle("MOB-Typer relaxase types")
## save the plot
ggsave("../results/S16Fig.pdf", S16Fig)


########################################
## oriT analysis.

## plot over oriT
S17Fig <- PIRA.PCN.for.AresArroyo2023.data %>%
    rank.correlate.column("oriT") %>%
    ## only plot groups with more than 10 data points.
    filter.correlate.column("oriT") %>%
    make_PCN_base_plot() +
    facet_wrap(oriT ~ .) +
    ggtitle("oriT annotated by Ares-Arroyo et al. (2023)")
## save the plot
ggsave("../results/S17Fig.pdf", S17Fig)


## plot over oriT_Family
S18Fig <- PIRA.PCN.for.AresArroyo2023.data %>%
    rank.correlate.column("oriT_Family") %>%
    ## only plot groups with more than 10 data points.
    filter.correlate.column("oriT_Family") %>%
    make_PCN_base_plot() +
    facet_wrap(oriT_Family ~ .) +
    ggtitle("oriT family annotated by Ares-Arroyo et al. (2023)")
## save the plot
ggsave("../results/S18Fig.pdf", S18Fig,height=8,width=8)

##############################
## Host range analysis.

## plot over observed host range.
## This is: Taxon name of convergence of plasmids in MOB-suite plasmid DB
## See documentation here: https://github.com/phac-nml/mob-suite 
S19Fig <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("observed_host_range_ncbi_name") %>%
    filter.correlate.column("observed_host_range_ncbi_name") %>%
    make_PCN_base_plot() +
    facet_wrap(observed_host_range_ncbi_name ~ ., ncol=3) +
    ggtitle("Host range annotated by MOB-Typer")
## save the plot
ggsave("../results/S19Fig.pdf", S19Fig, height=12,width=10,limitsize=FALSE)


## 415 plasmids have copy numbers and host ranges.
RedondoSalvo.host.range.PIRA.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.PTU.data) %>%
    filter(Host_range != "-") %>%
    filter(!is.na(Host_range))

## This figure shows that copy number / plasmid size does not predict host range.
## there are narrow and broad host range plasmids both large and small.
## compare with the MOB-Typer result in this vein.
S20Fig <- RedondoSalvo.host.range.PIRA.estimates %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Host_range)) +
    geom_point(size=1,alpha=0.5) +
    theme_classic()
## save the plot
ggsave("../results/S20Fig.pdf", S20Fig, height=3.5, width=6)


################################################################################
## Supplementary Figure S21. let's make a histogram of PCN in these data.

S21Fig <- PIRA.PCN.estimates %>%
    ## CRITICAL TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank") %>%
    ggplot(aes(x = log10(PIRACopyNumber))) +
    geom_histogram(bins=1000) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "light gray") +
    geom_vline(xintercept = 2, linetype = "dashed", color = "light gray") +
    facet_wrap(. ~ Annotation) +
    theme_classic()

ggsave("../results/S21Fig.pdf", S21Fig, height = 6, width = 6)

################################################################################
## Supplementary Figure S22.
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
S22FigA <- make.confint.figure.panel(low.PCN.plasmids.table, order.by.total.plasmids, "proportion of plasmids with PCN < 1")

## calculate the fraction of high PCN plasmids in each category.
## there is an enrichment of  PCN > 50 plasmids in human-impacted environments.
## Make Z-distributed confidence intervals for the fraction of isolates with
## PCN > 50.

high.PCN.plasmids.table <- make.highPCN.table(PIRA.PCN.estimates)
## plot the confidence intervals to see if there is any enrichment of high PCN plasmids in any ecological category.
S22FigB <- make.confint.figure.panel(high.PCN.plasmids.table, order.by.total.plasmids, "proportion of plasmids with PCN > 50")

S22Fig <- plot_grid(S22FigA, S22FigB, labels=c('A','B'), ncol=1)
ggsave("../results/S22Fig.pdf", S22Fig, height = 8, width = 6)

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
## Supplementary Figure S23. Chromosome DNA content does not constrain Plasmid DNA content.

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


S23Fig <- DNA.content.data %>%
    ggplot(aes(
        x = log10(ChromosomeDNAContent),
        y = log10(PlasmidDNAContent))) +
    geom_point(size=0.2) +
    facet_wrap(. ~ Annotation) +
    theme_classic() +
    theme(strip.background = element_blank())
## save the plot
ggsave("../results/S23Fig.pdf", S23Fig, height=5)


###################################################################################
## Figure 3. Coding Sequences (CDS) on plasmids follow an empirical scaling law.

## Main figure, all the points together.
## supplementary figure: same figure, separated by Annotation category.

## Fig 3A: show the combined plot
Fig3A <- CDS.MGE.ARG.fraction.data %>%
    make_CDS_scaling_base_plot()

## Fig 3B: show generality over ecology.
Fig3B <- CDS.MGE.ARG.fraction.data %>%    
    make_CDS_scaling_base_plot() +
    facet_wrap(.~Annotation)

Fig3 <- plot_grid(Fig3A, Fig3B, labels = c("A", "B"), nrow=1)
## save the plot.
ggsave("../results/Fig3.pdf", Fig3, height=4, width=7.5)


################################################################################
## Supplementary Figures  S24 through S6. Break down the result in Figure 3 by taxonomy
## and ecological category to show universality of the CDS scaling relationship.

## Supplementary Figure S24.
## Break down by taxonomic group.
S24Fig <- CDS.MGE.ARG.fraction.data %>%
    filter.correlate.column("TaxonomicGroup") %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicGroup)
## save the plot.
ggsave("../results/S24Fig.pdf", S24Fig, height=5)


## Supplementary Figure S25
## Break down by taxonomic subgroup
S25Fig <- CDS.MGE.ARG.fraction.data %>%
    filter.correlate.column("TaxonomicSubgroup") %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup)
## save the plot.
ggsave("../results/S25Fig.pdf", S25Fig, height=6, width=10)


## Supplementary FIgure S26
## Break down by genus.
S26Fig <- CDS.MGE.ARG.fraction.data %>%
    filter.correlate.column("Genus") %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ Genus, ncol=12)

## save the plot.
ggsave("../results/S26Fig.pdf", S26Fig, height=15, width = 14)


########################################################################
## Figure 4: plasmid size scales with metabolic capacity.

## CRITICAL TODO: we will need to carefully distinguish between true zeros and NA
## values created by a plasmid/chromosome not being included in these data.

## IMPORTANT: removing plasmids longer than 500kB removes the scaling law.
## Therefore, the metabolic.scaling law largely holds for chromids and chromosomes.

metabolic.gene.plasmid.data <- plasmid.length.data %>%
    left_join(metabolic.genes.in.plasmids) %>%
    ## make the dataframe compatible with plasmid.annotation.data
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join with plasmid.annotation.data.
    left_join(plasmid.annotation.data) %>%
    ## get CDS data for each genome.
    left_join(CDS.MGE.ARG.fraction.data) %>%
    ## set NA values of metabolic_protein_count to zeros.
    mutate(metabolic_protein_count = ifelse(is.na(metabolic_protein_count), 0, metabolic_protein_count))

## join chromosome.length.data to metabolic.genes.in.chromosomes,
## since only 100 genomes were examined.
## KEY ASSUMPTION: each chromosome has at least one metabolic gene.
metabolic.gene.chromosome.data <- metabolic.genes.in.chromosomes %>%
    left_join(replicon.length.data) %>%
    ## get CDS data for each genome.
    left_join(CDS.MGE.ARG.fraction.data)

## combine plasmid and chromosome metabolic gene data for Figure 4.
metabolic.gene.plasmid.and.chromosome.data <- metabolic.gene.plasmid.data %>%
    full_join(metabolic.gene.chromosome.data) %>%
    ## IMPORTANT TODO: filter upstream Annotation so I don't need to do this ad hoc filtering here.
    filter(!is.na(Annotation)) %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType))
    

## Fig4A: show the combined plot
Fig4A <- metabolic.gene.plasmid.and.chromosome.data %>%
    make_metabolic_scaling_base_plot()

## Fig4B: show generality over ecology.
Fig4B <- metabolic.gene.plasmid.and.chromosome.data %>%    
    make_metabolic_scaling_base_plot() +
    facet_wrap(.~Annotation, nrow=3)

Fig4 <- plot_grid(Fig4A, Fig4B, labels = c('A', 'B'))
## save the plot.
ggsave("../results/Fig4.pdf", Fig4, height=4, width=7.5)


################################################################################
## Supplementary Figures S27 through S29. Break down the result in Figure 4 by taxonomy
## and ecological category to show universality of the CDS scaling relationship.

## Supplementary Figure S27
## Break down by taxonomic group.
S27Fig <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter.correlate.column("TaxonomicGroup") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicGroup)
## save the plot.
ggsave("../results/S27Fig.pdf", S27Fig, height=6,width=8)


## Supplementary Figure S28
## Break down by taxonomic subgroup
S28Fig <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter.correlate.column("TaxonomicSubgroup") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup, ncol=3)
## save the plot.
ggsave("../results/S28Fig.pdf", S28Fig, height=12, width=8)


## Supplementary FIgure S29
## Break down by genus.
S29Fig <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter.correlate.column("Genus") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ Genus, ncol=50)
## save the plot.
ggsave("../results/S29Fig.pdf", S29Fig, height=50, width=50, limitsize = FALSE)


