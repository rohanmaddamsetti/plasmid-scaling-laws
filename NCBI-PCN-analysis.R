## NCBI-PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## CRITICAL TODO:
## THIS Supplementary Figure does not look right!!
## is this a bug in the data???

## TODO: make a figure comparing the fit between PIRA and Naive Themisto to the alignment methods
## (minimap2 and breseq).

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


################################################################################
## functions and global variables.

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
## get the PIRA PCN estimates. These are the main data for this paper.

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
    filter(PIRAReadCount > MIN_READ_COUNT)


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
## Supplementary Figure S2: Effect of PIRA compared to naive themisto estimates.

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
    

## Now make S2 Figure panel A.
S2FigA <- PIRA.vs.naive.themisto.df %>%
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


## S2 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S2FigB <- PIRA.vs.naive.themisto.df %>%
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

## S2 Figure panels C and D remove points with insufficient reads.

## Now make S2 Figure panel C
S2FigC <- PIRA.vs.naive.themisto.df %>%
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


## S2 Figure panel D zooms in on the plasmids with PIRA PCN < 0.8.
S2FigD <- PIRA.vs.naive.themisto.df %>%
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


## make Supplementary Figure S2.
S2Fig <- plot_grid(S2FigA, S2FigB, S2FigC, S2FigD, labels=c("A", "B", "C", "D"))
ggsave("../results/S2Fig.pdf", S2Fig, height=6, width=8)



################################################################################
## Supplementary Figure S3.
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

## make S3 Figure panel A
S3FigA <- PIRA.vs.minimap2.df %>%
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

## S3 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S3FigB <- PIRA.vs.minimap2.df %>%
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

## Now make Supplementary Figure S3.
S3Fig <- plot_grid(S3FigA, S3FigB, labels=c("A", "B"))
ggsave("../results/S3Fig.pdf", S3Fig, height=4, width=8)

## make a linear model and examine it.
minimap2.PIRA.PCN.lm.model <- lm(
    formula = log10(PIRACopyNumber) ~ log10(minimap2_PIRA_CopyNumberEstimate),
    data = PIRA.vs.minimap2.df)
## look at the linear regression.
summary(minimap2.PIRA.PCN.lm.model)


###################################################################################
## Supplementary Figure S4.
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


## Now make S4 Figure panel A.
S4FigA <- PIRA.vs.breseq.df %>%
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


## S4 Figure panel B zooms in on the plasmids with PIRA PCN < 0.8.
S4FigB <- PIRA.vs.breseq.df %>%
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

## Now make Supplementary Figure S4.
S4Fig <- plot_grid(S4FigA, S4FigB, labels=c("A", "B"))
ggsave("../results/S4Fig.pdf", S4Fig, height=4, width=8)


## make a linear model and examine it.
breseq.PIRA.PCN.lm.model <- lm(
    formula = log10(PIRACopyNumber) ~ log10(BreseqCopyNumberEstimate),
    data = PIRA.vs.breseq.df)
## look at the linear regression.
summary(breseq.PIRA.PCN.lm.model)

################################################################################
## Supplementary Figure S5:
## compare naive kallisto to naive themisto PCN estimates to show that PCN numbers by pseudoalignment
## are reproducible irrespective of the specific software implementation.

## kallisto replicon mapping PCN estimates.
kallisto.replicon.metadata <- read.csv("../results/NCBI-replicon_lengths.csv")

kallisto.replicon.PCN.estimates <- read.csv("../results/kallisto-replicon_copy_numbers.csv") %>%
    rename(KallistoNaiveCopyNumber = CopyNumber) %>%
    filter(SeqType == "plasmid") %>%
    left_join(kallisto.replicon.metadata)

## compare naive themisto to naive kallisto estimates.
naive.themisto.vs.naive.kallisto.df <- naive.themisto.PCN.estimates %>%
    ## IMPORTANT: remove plasmids with insufficient reads.
    filter(InsufficientReads == FALSE) %>%
    left_join(kallisto.replicon.PCN.estimates) %>%
    filter(SeqType == "plasmid")

## make Supplementary Figure S5.   
S5Fig <- naive.themisto.vs.naive.kallisto.df %>%
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
## save Supplementary Figure S5.
S5Fig <- ggsave("../results/S5Fig.pdf", S5Fig, height=4, width=4)

################################################################################
## BIOLOGY ANALYSES STARTING HERE

## This is the full file of MOBTyper results
## CRITICAL TODO: analyze all the additional correlates that may be in these data.
MOBTyper.results <- read.csv("../data/Maddamsetti2024_FileS5-MOBTyper-plasmid-annotations.csv")


## This is just the mobiliity data from the above big table.
## Import MOBTyper mobility results for merging with the plasmid copy number data.
mobility.results <- read.csv("../results/mobility-results.csv")

## get ARG copy number data-- this is only used for annotating ARGs.
kallisto.ARG.copy.number.data <- read.csv("../results/kallisto-ARG_copy_numbers.csv") %>%
    mutate(beta.lactam.resistance = ifelse(str_detect(product,beta.lactam.keywords), TRUE, FALSE))
## this is used for annotating plasmids containing beta-lactamases.
beta.lactam.ARGs <- filter(kallisto.ARG.copy.number.data, beta.lactam.resistance==TRUE)


## take the PIRA estimates, filter for plasmids,
## and annotate plasmids with ARGs, and MOB type.
PIRA.PCN.estimates <- PIRA.estimates %>%
    filter(SeqType == "plasmid") %>%
    ## add Plasmid column to merge mobility.results,
    ## by splitting the SeqID string on the underscore and taking the second part.
    mutate(Plasmid = sapply(strsplit(SeqID, "_"), function(x) x[2])) %>%
    left_join(mobility.results) %>%
    mutate(has.ARG = ifelse(SeqID %in% kallisto.ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification,
                                    `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases"))

## Make a column representing plasmid length normalized by chromosome length in each AnnotationAccession
# Group the data frame by AnnotationAccession and calculate the maximum replicon length in each group
max_replicon_lengths <- PIRA.PCN.estimates %>%
  group_by(`AnnotationAccession`) %>%
  summarise(max_replicon_length = max(replicon_length))

# Merging the maximum replicon lengths back to the original data frame
PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    left_join(max_replicon_lengths, by = "AnnotationAccession") %>%
    ## Creating a new column 'normalized_replicon_length' by dividing 'replicon_length' by 'max_replicon_length'
    mutate(normalized_replicon_length = replicon_length / max_replicon_length)


## CRITICAL TODO: POLISH THIS CODE
## shape := ARG- or ARG+ (don't worry about beta-lactamases)
## color := conjugative, mobilizable, etc.
## keep consistent coding throughout the paper.

## plot the PIRA PCN estimates. 10,261 plasmids in these data.
Fig3A <- PIRA.PCN.estimates %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        shape = has.ARG,
        color = PredictedMobility)) +
    geom_point(size=0.2,alpha=0.5) +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    theme_classic() +
    xlab("log10(Length)")  +
    ylab("log10(Copy Number)") +
    theme(legend.position="bottom")

## Break down this result by predicted plasmid mobility.
Fig3B <- Fig3A + facet_grid(`Plasmid class`~ PredictedMobility)

## make Figure 3.
#Fig3 <- plot_grid(Fig3A, Fig3B, labels=c('A', 'B'),nrow=2)
ggsave("../results/Fig3.pdf", Fig3A, height=4, width=4)

## Make a Supplementary Figure S6 that is the same as Figure 3,
## but plotting normalized plasmid length relative to the length of the longest
## chromosome.

## CRITICAL TODO:
## THIS Supplementary Figure does not look right!!
## is this a bug in the data???

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
S6FigA <- ggplot(PIRA.PCN.estimates,
                 aes(
                     x = normalized_replicon_length,
                     y = log10(PIRACopyNumber),
                     color = `Plasmid class`)) +
    geom_point(size=0.1, alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("Normalized Plasmid length") +
    theme(legend.position="right")

## Break down this result by predicted plasmid mobility.
S6FigB <- S6FigA + facet_grid(`Plasmid class`~ PredictedMobility)

## make Supplementary Figure S6.
S6Fig <- plot_grid(S6FigA, S6FigB, labels=c('A', 'B'),nrow=2)
ggsave("../results/S6Fig.pdf", S6Fig, height=12, width=12)


################################################################################
## Figure 4. PCN distribution over INC groups, MOB groups, and ecology.
## clear association between high PCN plasmids and particular ecological annotations.

## First examine over ecology.
ecological.annotation <- read.csv("../results/computationally-annotated-gbk-annotation-table.csv")

ecologically.annotated.PCN.data <- PIRA.PCN.estimates %>%
    left_join(ecological.annotation) %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank")

Fig4A <- ggplot(
    ecologically.annotated.PCN.data,
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
Fig4B <- ecologically.annotated.PCN.data %>%
    ggplot(aes(
        x = log10(PIRACopyNumber),
        y = Annotation,
        color = Annotation)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Ecological Annotation") +
    xlab("log10(Plasmid copy number)") 

Fig4 <- plot_grid(Fig4A, Fig4B, labels=c("A","B"), nrow=2)
## save the plot
ggsave("../results/Fig4.pdf", Fig4, height=12, width=12)


################################################################################
## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.

beta.lactamase.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.beta.lactamase==TRUE)

no.beta.lactamase.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.beta.lactamase == FALSE)

ARG.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- PIRA.PCN.estimates %>%
    filter(has.ARG == FALSE)

## plasmids with beta-lactamases have lower PCN than those without beta-lactamases.
mean(beta.lactamase.plasmid.data$PIRACopyNumber)
mean(no.beta.lactamase.plasmid.data$PIRACopyNumber)

median(beta.lactamase.plasmid.data$PIRACopyNumber)
median(no.beta.lactamase.plasmid.data$PIRACopyNumber)

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
## Supplementary Figure S7. let's make a histogram of PCN in these data.

S7Fig <- PIRA.PCN.estimates %>%
    ggplot(aes(x = log10(PIRACopyNumber))) +
    geom_histogram(bins=1000) +
    theme_classic()

ggsave("../results/S7Fig.pdf", S7Fig, height = 6, width = 6)


################################################################################
## Supplementary Figure S8. examine total DNA (chromosomes and plasmids) per genome.

DNA.content.data <- PIRA.estimates %>%
    mutate(DNAContent = replicon_length * PIRACopyNumber) %>%
    group_by(AnnotationAccession) %>%
    summarize(TotalDNAContent = sum(DNAContent))

S8Fig <- DNA.content.data %>%
    ggplot(aes(x=TotalDNAContent)) +
    geom_histogram(bins=1000) +
    theme_classic() +
## IMPORTANT TODO: EXAMINE OUTLIERS IN THIS PLOT.
    xlim(0,15000000)


ggsave("../results/S8Fig.pdf", S8Fig, height = 6, width = 6)

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


