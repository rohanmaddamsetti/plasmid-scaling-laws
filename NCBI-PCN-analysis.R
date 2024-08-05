## NCBI-PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## Make a scatterplot of plasmid copy numbers against plasmid length,
## and color dots by presence of ARGs.

## redo the same analysis, focusing on conjugative, mobilizable, and non-mobilizable plasmids.

library(tidyverse)
library(cowplot)
library(ggrepel)
library(data.table)


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

################################################################################










## get lengths of all the replicons.
NCBI.replicon.length.data <- read.csv("../results/NCBI-replicon_lengths.csv")

## get ARG copy number data.
kallisto.ARG.copy.number.data <- read.csv("../results/NCBI-ARG_copy_numbers.csv") %>%
    mutate(beta.lactam.resistance = ifelse(str_detect(product,beta.lactam.keywords), TRUE, FALSE))

beta.lactam.ARGs <- filter(NCBI.ARG.copy.number.data, beta.lactam.resistance==TRUE)
non.beta.lactam.ARGs <- filter(NCBI.ARG.copy.number.data, beta.lactam.resistance==FALSE)

kallisto.gene.averaged.copy.number.data <- read.csv("../results/NCBI-replicon_copy_numbers_from_genes.csv") %>%
    full_join(NCBI.replicon.length.data) %>%
    mutate(has.ARG = ifelse(SeqID %in% NCBI.ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification, `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases")) ##%>%
    ## remove outlier points with very low coverage.
##    filter(CopyNumber > 0.5)

## Make a column representing plasmid length normalized by chromosome length in each AnnotationAccession

# Group the data frame by AnnotationAccession and calculate the maximum replicon length in each group
max_replicon_lengths <- NCBI.chromosome.plasmid.copy.number.data %>%
  group_by(`AnnotationAccession`) %>%
  summarise(max_replicon_length = max(replicon_length))

# Merging the maximum replicon lengths back to the original data frame
NCBI.chromosome.plasmid.copy.number.data <- NCBI.chromosome.plasmid.copy.number.data %>%
  left_join(max_replicon_lengths, by = "AnnotationAccession") %>%
## Creating a new column 'normalized_replicon_length' by dividing 'replicon_length' by 'max_replicon_length'
    mutate(normalized_replicon_length = replicon_length / max_replicon_length)

## Import MOBTyper mobility results for merging with the plasmid copy number data.
mobility.results <- read.csv("../results/mobility-results.csv")


NCBI.plasmid.copy.number.data <- NCBI.chromosome.plasmid.copy.number.data %>%
    filter(SeqType == "plasmid") %>%
    arrange(CopyNumber) %>%
    ## add Plasmid column to merge mobility.results,
    ## by splitting the SeqID string on the underscore and taking the second part.
    mutate(Plasmid = sapply(strsplit(SeqID, "_"), function(x) x[2])) %>%
    left_join(mobility.results) %>%
    ## for now, remove rows without predicted mobility.
    filter(!is.na(PredictedMobility))

################################################################################
## get concrete statistics for how many genome papers report plasmid copy number:

## sample 50 genomes with plasmids at random, and examine how many papers report PCN.
## pick genomes with some plasmid with PCN > 10, so that we are focusing on genomes where people
## might actually want to report PCN.
## for reproducibility, use random seed = 60, which I generated using random.org (see screenshot in ../data).

head(NCBI.plasmid.copy.number.data)

genomes.to.sample.from.NCBI.plasmid.copy.number.data <- NCBI.plasmid.copy.number.data %>%
    filter(CopyNumber > 10) %>%
    select(AnnotationAccession) %>%
    distinct()

## set the random seed for reproducibility.
set.seed(60)

chosen.genomes <- sample(
    x = genomes.to.sample.from.NCBI.plasmid.copy.number.data$AnnotationAccession,
    replace = FALSE,
    size = 50,
    )

write.csv(chosen.genomes, "../results/Fifty-random-genomes-with-multicopy-plasmids.csv", row.names = FALSE)

################################################################################
## compare results with using each replicon as a 'gene' for read mapping with kallisto.
kallisto.replicon.data <- read.csv("../results/NCBI-replicon_copy_numbers.csv") %>%
    full_join(NCBI.replicon.length.data) %>%
    mutate(has.ARG = ifelse(SeqID %in% NCBI.ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification, `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases"))

NCBI.replicon.data.for.comparison <- NCBI.replicon.data %>%
    mutate(KallistoRepliconCopyNumber = CopyNumber) %>%
    select(RefSeqID, SeqID, SeqType, AnnotationAccession, KallistoRepliconCopyNumber)

## merge the gene-level and replicon-level estimates.
NCBI.plasmid.estimate.data <- NCBI.replicon.data.for.comparison %>%
    filter(SeqType == "plasmid") %>%
    full_join(NCBI.plasmid.copy.number.data) %>%
    mutate(KallistoGeneAveragedCopyNumber = CopyNumber)

kallisto.comparison.plot <- NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoGeneAveragedCopyNumber, y=KallistoRepliconCopyNumber)) +
    geom_point() + theme_classic()

## There are 1944 big outliers, on the ends of both axes!
kallisto.comparison.plot

## let's remove the outliers and see what the plot looks like.
cropped.kallisto.comparison.plot <- kallisto.comparison.plot +
    xlim(0,1000) + ylim(0,1000) +
    ## Add diagonal line to compare estimates for the whole replicon,
    ## compared to the average over all genes.
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoGeneAveragedCopyNumber,y=KallistoRepliconCopyNumber),
        color="light blue",
        formula=y~x)

cropped.kallisto.comparison.plot
ggsave("../results/NCBI-kallisto-PCN-method-comparison.pdf",
       cropped.kallisto.comparison.plot,height=5.75,width=5.75)


## make a linear regression model:
## Kallisto gene-level PCN estimates vs. replicon-level PCN estimates.
## but first, trim the data.
trimmed.plasmid.estimate.data <- NCBI.plasmid.estimate.data %>%
    filter(KallistoGeneAveragedCopyNumber < 1000) %>%
    filter(KallistoRepliconCopyNumber < 1000)

kallisto.PCN.lm.model <- lm(
    formula = KallistoRepliconCopyNumber ~ KallistoGeneAveragedCopyNumber,
    data=trimmed.plasmid.estimate.data)
## look at the linear regression.
summary(kallisto.PCN.lm.model)

################################################################################

## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.

beta.lactamase.plasmid.data <- NCBI.plasmid.copy.number.data %>%
    filter(has.beta.lactamase==TRUE)

no.beta.lactamase.plasmid.data <- NCBI.plasmid.copy.number.data %>%
    filter(has.beta.lactamase == FALSE)

ARG.plasmid.data <- NCBI.plasmid.copy.number.data %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- NCBI.plasmid.copy.number.data %>%
    filter(has.ARG == FALSE)

## plasmids with beta-lactamases have lower PCN than those without beta-lactamases.
mean(beta.lactamase.plasmid.data$CopyNumber)
mean(no.beta.lactamase.plasmid.data$CopyNumber)

median(beta.lactamase.plasmid.data$CopyNumber)
median(no.beta.lactamase.plasmid.data$CopyNumber)

## plasmids with ARGs have lower PCN than those without ARGs.
mean(ARG.plasmid.data$CopyNumber)
mean(no.ARG.plasmid.data$CopyNumber)

median(ARG.plasmid.data$CopyNumber)
median(no.ARG.plasmid.data$CopyNumber)

## make a linear regression model:
## log10(Plasmid copy number) vs. log10(Plasmid length).
#plasmid.lm.model <- lm(
#    formula=log10(CopyNumber)~log10(replicon_length),
 #   data=NCBI.plasmid.copy.number.data)
## look at the linear regression.
#summary(plasmid.lm.model)

#second.order.plasmid.lm.model <- lm(
#    formula=log10(CopyNumber)~poly(log10(replicon_length),2,raw=TRUE),
#    data=NCBI.plasmid.copy.number.data)
## look at the second order regression.
#summary(second.order.plasmid.lm.model)

## let's compare these models. The model with lower AIC is better.
#AIC(plasmid.lm.model)
#AIC(second.order.plasmid.lm.model)
## also compare the models using an ANOVA.
#print(anova(
 #   plasmid.lm.model,
 #   second.order.plasmid.lm.model))

## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
NCBI.plasmid.copy.number.plot <- ggplot(NCBI.plasmid.copy.number.data,
                                   aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    facet_grid(`Plasmid class`~PredictedMobility)
    ## add the linear regression.
    #geom_smooth(
     #   data=NCBI.plasmid.copy.number.data,
      #  inherit.aes=FALSE,
       # method='lm',
      #  aes(x=log10(replicon_length),y=log10(CopyNumber)),
      #  color="light blue",
      #  formula=y~x) +
    ## let's look at second-order polynomial fit.
    #geom_smooth(
     #   data=NCBI.plasmid.copy.number.data,
      #  inherit.aes=FALSE,
       # method='lm',
       # aes(x=log10(replicon_length),y=log10(CopyNumber)),
       # color="light green",
        #formula=y~poly(x, 2, raw=TRUE))
## save the plot.
ggsave("../results/NCBI2024-plasmid-copy-number.pdf",
       NCBI.plasmid.copy.number.plot,height=5.75,width=5.75)

## repeat this plot, using kallisto replicon-level estimates.
NCBI.plasmid.replicon.copy.number.plot <- NCBI.replicon.data %>%
    filter(SeqType == "plasmid") %>%
    ggplot(
        aes(x=log10(replicon_length),y=log10(CopyNumber),
            color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top")
## save the plot.
ggsave("../results/NCBI2024-replicon-plasmid-copy-number.pdf",
       NCBI.plasmid.replicon.copy.number.plot,height=5.75,width=5.75)

## plot for Duke retreat
Duke.Retreat.plasmid.replicon.copy.number.plot <- NCBI.replicon.data %>%
    filter(SeqType == "plasmid") %>%
    filter(CopyNumber > 0.1) %>%
    filter(CopyNumber < 5000) %>%
    ggplot(
        aes(x=log10(replicon_length),y=log10(CopyNumber),
            color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    guides(color = 'none')
## save the plot.
ggsave("../results/Duke-Retreat-2024-replicon-plasmid-copy-number.pdf",
       Duke.Retreat.plasmid.replicon.copy.number.plot,height=4,width=4)



## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
normalized.NCBI.plasmid.copy.number.plot1 <- ggplot(NCBI.plasmid.copy.number.data,
                                   aes(x=log10(normalized_replicon_length),y=log10(CopyNumber),
                                       color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    facet_grid(`Plasmid class`~PredictedMobility)

ggsave("../results/NCBI2024-plasmid-copy-number-normalized-length-facets.pdf",
       normalized.NCBI.plasmid.copy.number.plot1,height=5.75,width=5.75)

## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
normalized.NCBI.plasmid.copy.number.plot2 <- ggplot(NCBI.plasmid.copy.number.data,
                                   aes(x=log10(normalized_replicon_length),y=log10(CopyNumber),
                                       color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_vline(xintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top")

ggsave("../results/NCBI2024-plasmid-copy-number-normalized-length.pdf",
       normalized.NCBI.plasmid.copy.number.plot2,height=5.75,width=5.75)

## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
NCBI.plasmid.copy.number.plot2 <- ggplot(NCBI.plasmid.copy.number.data,
                                   aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="bottom")

ggsave("../results/NCBI2024-PCN.pdf",
       NCBI.plasmid.copy.number.plot2,height=5.75,width=5.75)

## just plot plasmids with PCN > 100.
NCBI.plasmid.copy.number.plot3 <- NCBI.plasmid.copy.number.data %>%
    mutate(highPCN =  ifelse(CopyNumber > 100, TRUE, FALSE)) %>%
    ggplot(aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=highPCN)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="bottom")

ggsave("../results/NCBI2024-high-PCN.pdf",
       NCBI.plasmid.copy.number.plot3,height=5.75,width=5.75)

## calculate the total number of plasmids and the number of plasmids with PCN > 100.
nrow(NCBI.plasmid.copy.number.data) ## 10,729 plasmids here

NCBI.plasmid.copy.number.data %>% filter(CopyNumber > 100) %>% nrow() ## 265 plasmids here.

################################################################################
## let's make a histogram of PCN in these data.

NCBI.plasmid.histogram <- NCBI.plasmid.copy.number.data %>%
    ggplot(aes(x = log10(CopyNumber))) +
    geom_histogram(bins=1000) +
    theme_classic()

NCBI.plasmid.histogram


################################################################################
## let's examine total plasmid DNA per genome.

plasmid.DNA.content.data <- NCBI.plasmid.copy.number.data %>%
    filter(SeqType == "plasmid") %>%
    mutate(DNAContent = replicon_length * CopyNumber) %>%
    group_by(AnnotationAccession, max_replicon_length) %>%
    summarize(TotalPlasmidDNAContent = sum(DNAContent))

plasmid.DNA.content.plot <- plasmid.DNA.content.data %>%
    ggplot(aes(x=max_replicon_length, y=TotalPlasmidDNAContent)) +
    geom_point() +
    theme_classic() +
    geom_smooth()

plasmid.DNA.content.plot

plasmid.DNA.content.plot2 <- plasmid.DNA.content.plot + ylim(0,2000000)

plasmid.DNA.content.plot2


#################################################################################
## compare PCN estimates using themisto versus kallisto.

## naive calculation with themisto read counts, ignoring multireplicon reads.
naive.themisto.PCN.estimates <- read.csv("../results/naive-themisto-PCN-estimates.csv") %>%
    rename(
        ThemistoNaiveCopyNumber = CopyNumber,
        ThemistoReadCount = ReadCount,
        ThemistoNaiveSequencingCoverage = SequencingCoverage,
        ThemistoNaiveLongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(SeqType == "plasmid") 

## now compare to kallisto estimates.
kallisto.themisto.NCBI.plasmid.estimate.data <- NCBI.plasmid.estimate.data %>%
    inner_join(naive.themisto.PCN.estimates)

## compare naive themisto estimates to gene averaged kallisto estimates.
themisto.kallisto.comparison.plot1 <- kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoGeneAveragedCopyNumber, y=ThemistoNaiveCopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")

themisto.kallisto.comparison.plot1
## again, there are some big outliers.

## let's remove the outliers and see what the plot looks like.
cropped.themisto.kallisto.comparison.plot1 <- themisto.kallisto.comparison.plot1 +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoGeneAveragedCopyNumber, y=ThemistoNaiveCopyNumber),
        color="light blue",
        formula=y~x)

cropped.themisto.kallisto.comparison.plot1

## compare naive themisto estimates to kallisto replicon estimates.
    
themisto.kallisto.comparison.plot2 <- kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoRepliconCopyNumber, y=ThemistoNaiveCopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")

themisto.kallisto.comparison.plot2
## again, there are some big outliers.

## let's remove the outliers and see what the plot looks like.
cropped.themisto.kallisto.comparison.plot2 <- themisto.kallisto.comparison.plot2 +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoRepliconCopyNumber, y=ThemistoNaiveCopyNumber),
        color="light blue",
        formula=y~x)

cropped.themisto.kallisto.comparison.plot2
ggsave("../results/themisto-kallisto-naive-replicon-PCN-estimate-comparison.pdf",
       cropped.themisto.kallisto.comparison.plot2,
       height=5.75,width=5.75)


##################################################################
## now let's compare the simple themisto PCN estimates to the other methods.

## write out a list of all plasmids with PCN > 100, and SequencingCoverage > 0.1 and LongestRepliconCoverage > 0.1
## for structural variation analysis for the darwinian circuit project.
high.PCN.plasmids <- read.csv("../results/simple-themisto-PCN-estimates.csv") %>%
    filter(SeqType == "plasmid") %>%
    filter(SequencingCoverage > 0.1) %>%
    filter(LongestRepliconCoverage > 0.1) %>%
    filter(CopyNumber > 100) %>%
    ## remove polars cruft from the dataframe.
    select(-SeqType_right, -replicon_length_right) %>%
    ## and sort the rows.
    arrange(AnnotationAccession, desc(CopyNumber))
write.csv(high.PCN.plasmids, "../results/high-PCN-plasmids.csv",row.names=FALSE,quote=FALSE)


## simple calculation with themisto read counts, partitioning multireplicon reads equally among replicons.
simple.themisto.PCN.estimates <- read.csv("../results/simple-themisto-PCN-estimates.csv") %>%
    rename(
        ThemistoSimpleCopyNumber = CopyNumber,
        ThemistoSimpleCount = ReadCount,
        ThemistoSimpleSequencingCoverage = SequencingCoverage,
        ThemistoSimpleLongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(SeqType == "plasmid") 

## make the scatterplot.
NCBI.simple.themisto.plasmid.copy.number.plot <- ggplot(simple.themisto.PCN.estimates,
                                   aes(x=log10(replicon_length),y=log10(ThemistoSimpleCopyNumber))) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_hline(yintercept=2,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top")

NCBI.simple.themisto.plasmid.copy.number.plot


full.kallisto.themisto.NCBI.plasmid.estimate.data <- kallisto.themisto.NCBI.plasmid.estimate.data %>%
    full_join(simple.themisto.PCN.estimates)

cropped.themisto.kallisto.comparison.plot3 <- full.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoGeneAveragedCopyNumber, y=ThemistoSimpleCopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoGeneAveragedCopyNumber, y=ThemistoSimpleCopyNumber),
        color="light blue",
        formula=y~x)

cropped.themisto.kallisto.comparison.plot3

cropped.themisto.kallisto.comparison.plot4 <- full.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=ThemistoNaiveCopyNumber, y=ThemistoSimpleCopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=ThemistoNaiveCopyNumber, y=ThemistoSimpleCopyNumber),
        color="light blue",
        formula=y~x)

cropped.themisto.kallisto.comparison.plot4

ggsave("../results/simple-themisto-naive-themisto-PCN-estimate-comparison.pdf",
       cropped.themisto.kallisto.comparison.plot4,
       height=5.75,width=5.75)

###################################################################################
## now let's compare the PIRA PCN estimates to the simple and naive themisto PCN estimates, and to the
## kallisto estimates.

## This analysis shows that PIRA essentially gives the same estimate as the themisto naive method.
## the themisto naive method is the fastest method, and correlates well with the kallisto replicon method.
## from a computational point of view, PIRA is not so good.
## however, PIRA does yield PCN estimates for 2212 plasmids for which the other methods don't get any
## answer, due to multireads being excluded from the analysis.

PIRA.PCN.estimates <- read.csv("../results/PIRA-PCN-estimates.csv") %>%
    rename(
        PIRACopyNumber = PIRA_CopyNumberEstimate,
        PIRAReadCount = ReadCount,
        PIRASequencingCoverage = SequencingCoverage,
        PIRALongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(SeqType == "plasmid") 

## now compare to simple and naive themisto estimates, and to kallisto estimates.
PIRA.kallisto.themisto.NCBI.plasmid.estimate.data <- PIRA.PCN.estimates %>%
    inner_join(full.kallisto.themisto.NCBI.plasmid.estimate.data) %>%
    ## only examine plasmids with at least 10,000 mapped reads.
    filter(PIRAReadCount > 10000)


PIRA.naive.themisto.comparison.plot <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=ThemistoNaiveCopyNumber, y=PIRACopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=ThemistoNaiveCopyNumber, y=PIRACopyNumber),
        color="light blue",
        formula=y~x)

ggsave("../results/PIRA-naive-themisto-PCN-estimate-comparison.pdf",
       PIRA.naive.themisto.comparison.plot,
       height=5.75,width=5.75)

## take a look at the statistics of the fit
PIRA.naive.PCN.lm.model <- lm(
    formula=PIRACopyNumber~ThemistoNaiveCopyNumber,
    data=PIRA.kallisto.themisto.NCBI.plasmid.estimate.data)
## look at the linear regression.
summary(PIRA.naive.PCN.lm.model)


PIRA.simple.themisto.comparison.plot <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=ThemistoSimpleCopyNumber, y=PIRACopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=ThemistoSimpleCopyNumber, y=PIRACopyNumber),
        color="light blue",
        formula=y~x)

## now compare to kallisto estimates.
PIRA.kallisto.comparison.plot1 <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoGeneAveragedCopyNumber, y=PIRACopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoGeneAveragedCopyNumber, y=PIRACopyNumber),
        color="light blue",
        formula=y~x)

PIRA.kallisto.comparison.plot2 <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(aes(x=KallistoRepliconCopyNumber, y=PIRACopyNumber)) +
    geom_point() + theme_classic() + geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    xlim(0,1000) + ylim(0,1000) +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=KallistoRepliconCopyNumber, y=PIRACopyNumber),
        color="light blue",
        formula=y~x)


## do we get more plasmid estimates (AKA plasmids that only have multiread coverage) when we use PIRA?
## YES!!!
PIRA.kallisto.themisto.NCBI.plasmid.estimate.data$PIRACopyNumber
PIRA.kallisto.themisto.NCBI.plasmid.estimate.data$ThemistoNaiveCopyNumber


## PIRA gets estimates for 12,829 plasmids.
PIRA.PCN.plot <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    ggplot(
        aes(x=log10(replicon_length),y=log10(PIRACopyNumber),
            color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") ##+
    #facet_grid(`Plasmid class`~PredictedMobility)


## filter based on PIRAReadCount
PIRA.PCN.plot2 <- PIRA.kallisto.themisto.NCBI.plasmid.estimate.data %>%
    filter(PIRAReadCount > 10000) %>%
    ggplot(
        aes(x=log10(replicon_length),y=log10(PIRACopyNumber),
            color=`Plasmid class`)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top")

PIRA.PCN.plot


##################################################################
## let's examine PCN distribution over INC groups, MOB groups, and ecology.
## clear association between high PCN plasmids and particular ecological annotations.

## First examine over ecology.
ecological.annotation <- read.csv("../results/computationally-annotated-gbk-annotation-table.csv")

ecologically.annotated.PCN.data <- simple.themisto.PCN.estimates %>%
    left_join(ecological.annotation) %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank")

ecological.NCBI.simple.themisto.plasmid.copy.number.plot <- ggplot(
    ecologically.annotated.PCN.data,
    aes(
        x = log10(replicon_length),
        y = log10(ThemistoSimpleCopyNumber),
        color = Annotation)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    geom_hline(yintercept=2,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    guides(color = "none") +
    facet_wrap(.~Annotation)

ecological.NCBI.simple.themisto.plasmid.copy.number.plot
## save the plot
ggsave(
    "../results/ecological-associations-with-PCN.pdf",
    ecological.NCBI.simple.themisto.plasmid.copy.number.plot, height=5, width=5)



## let's examine the ecology of high copy number plasmids.
## very cool! we see association with high PCN with humans, human-impacted environments, and livestock.

ecological.high.plasmid.copy.number.plot <- ecologically.annotated.PCN.data %>%
    filter(ThemistoSimpleCopyNumber > 100) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(ThemistoSimpleCopyNumber),
        color = Annotation)) +
    geom_point(size=0.2,alpha=0.5) +
    theme_classic() +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)")

ecological.high.plasmid.copy.number.plot



MOBTyper.results <- read.csv("../data/Maddamsetti2024_FileS5-MOBTyper-plasmid-annotations.csv")



###################################################################################
## let's compare long read PCN estimates from minimap2 and pseudoalignment for high PCN plasmids--
## This shows basically no correlation between PCN estimates! This, along with the paper
## "Recovery of small plasmid sequences via Oxford Nanopore sequencing" I found
## I found about inferring PCN with Oxford nanopore suggests that long-read sequencing
## is not appropriate for inferring PCN.

## merge high.PCN.plasmid.RunID.df with high.PCN.plasmids.
high.PCN.plasmid.RunID.df <- read.csv("../results/high-PCN-plasmid-RunID_table.csv")

longread.high.PCN.estimates <- read.csv("../results/longread-alignment-PCN-estimates.csv") %>%
    rename(longread.minimap2.CopyNumber = CopyNumber)

## 126 plasmids in this comparison.
long.read.alignment.versus.short.read.pseudoalignment.high.PCN.estimates <- high.PCN.plasmids %>%
    ## add RefSeq_ID column to merge RunID.df,
    ## by splitting the AnnotationAccession string on the second underscore in the string.
    mutate(RefSeq_ID = sapply(strsplit(AnnotationAccession, "_"), function(x) paste(x[1], x[2], sep="_"))) %>%
    inner_join(high.PCN.plasmid.RunID.df) %>%
    ## there are duplicate rows??? remove these.
    distinct() %>%
    inner_join(longread.high.PCN.estimates) %>%
    rename(shortread.themisto.CopyNumber = CopyNumber) %>%
    distinct()


high.PCN.estimate.comparison.plot <- long.read.alignment.versus.short.read.pseudoalignment.high.PCN.estimates %>%
    ggplot(aes(
        x = log10(shortread.themisto.CopyNumber),
        y = log10(longread.minimap2.CopyNumber),
        color = LongReadDataType)) +
    geom_point() +
    theme_classic() +
    xlab("log10(Short Read Themisto Copy Number)")  +
    ylab("log10(Long Read Minimap2 Copy Number)") +
    theme(legend.position="top") +
    xlim(-2,3) +
    ylim(-2,3)

ggsave(
    "../results/long-read-alignment-versus-short-read-pseudoalignment-high-PCN-estimates.pdf",
    high.PCN.estimate.comparison.plot, height=5, width=5)

###################################################################################
## Benchmarking of PIRA with themisto against PIRA with minimap2 alignments on 100 random genomes
## with low copy number plasmids (PCN < 0.8) to see whether these results are trustworthy.

## BIG QUESTION: why is the sequencing coverage so low? Might this be artifactual?
## I should probably re-run these 100 random genomes with breseq as another "gold standard control"
## to make sure these estimates are accurate, and not artifactual due to low sequencing coverage.

low.PCN.PIRA.themisto.estimates.df <- read.csv("../results/PIRA-low-PCN-benchmark-estimates.csv") %>%
    rename("ThemistoReadCount" = ReadCount) %>%
    select("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length", "ThemistoReadCount", "PIRA_CopyNumberEstimate")

low.PCN.PIRA.minimap2.estimates.df <- read.csv("../results/minimap2-PIRA-low-PCN-benchmark-estimates.csv") %>%
    rename("minimap2ReadCount" = ReadCount) %>%
    select("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length", "minimap2ReadCount", "minimap2_PIRA_CopyNumberEstimate")

PIRA.alignment.versus.pseudoalignment.comparison.df <- full_join(
    low.PCN.PIRA.themisto.estimates.df, low.PCN.PIRA.minimap2.estimates.df) %>%
    ## IMPORTANT: only examine plasmids.
    filter(SeqType == "plasmid")


PIRA.alignment.versus.pseudoalignment.comparison.plot <- PIRA.alignment.versus.pseudoalignment.comparison.df %>%
    ggplot(aes(
        x = log10(minimap2_PIRA_CopyNumberEstimate),
        y = log10(PIRA_CopyNumberEstimate))) +
    geom_point() +
    theme_classic() +
    xlab("log10(Short Read pseudoalignment Copy Number)")  +
    ylab("log10(Short Read alignment Copy Number)")
## Nice result, a few cases of underestimation for low PCN plasmids with themisto, but by and large highly accurate.
ggsave(
    "../results/short-read-alignment-versus-short-read-pseudoalignment-low-PCN-estimates.pdf",
    PIRA.alignment.versus.pseudoalignment.comparison.plot, height=5, width=5)

## now zoom in on the plasmids with minimap2 alignment PCN < 0.8.
## SUPER IMPORTANT: These data show that there are more aligned reads with themisto than with minimap2,
## given the parameters that I am using. Perhaps the default parameters I am using for alignment is very stringent?
low.PCN.PIRA.alignment.versus.pseudoalignment.comparison.df <- PIRA.alignment.versus.pseudoalignment.comparison.df %>% filter(minimap2_PIRA_CopyNumberEstimate < 0.8)

low.PCN.PIRA.alignment.versus.pseudoalignment.comparison.plot <- low.PCN.PIRA.alignment.versus.pseudoalignment.comparison.df %>%
        ggplot(aes(
        x = log10(minimap2_PIRA_CopyNumberEstimate),
        y = log10(PIRA_CopyNumberEstimate))) +
    geom_point() +
    theme_classic() +
    xlab("log10(Short Read pseudoalignment Copy Number)")  +
    ylab("log10(Short Read alignment Copy Number)")

ggsave(
    "../results/only-low-PCN-short-read-alignment-versus-short-read-pseudoalignment-low-PCN-estimates.pdf",
    low.PCN.PIRA.alignment.versus.pseudoalignment.comparison.plot, height=5, width=5)

###################################################################################
## Benchmarking of these 100 random genomes with breseq as another gold standard control,
## to make sure these estimates are accurate, and not artifactual due to low sequencing coverage.

## first get the metadata we need from the PIRA estimates.
low.PCN.benchmark.metadata.df <- low.PCN.PIRA.themisto.estimates.df %>%
    select("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length") %>%
    ## IMPORTANT: trim the ".1" suffixes of the SeqIDs so that we can properly merge
    ## with the TrimmedSeqID in the breseq low PCN benchmark summary data (see below).
    mutate(TrimmedSeqID = str_remove(SeqID, "\\..*$"))

## now get the breseq coverage results for the benchmarking genomes.
low.PCN.breseq.summary.df <- read.csv("../results/breseq-low-PCN-benchmark-estimates.csv") %>%
    full_join(low.PCN.benchmark.metadata.df) %>%
    filter(!is.na(mean_coverage))

## make a separate df with a column for the coverage for the longest replicon in each genome
low.PCN.breseq.chromosomal.coverage.df <- low.PCN.breseq.summary.df %>%
    group_by(AnnotationAccession) %>%
    filter(replicon_length == max(replicon_length)) %>%
    ungroup() %>%
    rename("chromosomal_coverage" = mean_coverage) %>%
    select("AnnotationAccession", "chromosomal_coverage")

## now join these data back to the original breseq summary data
low.PCN.breseq.PCN.estimate.df <- low.PCN.breseq.summary.df %>%
    left_join(low.PCN.breseq.chromosomal.coverage.df) %>%
    group_by(AnnotationAccession) %>%
    ## calculate breseq PCN estimates here.
    mutate(BreseqCopyNumberEstimate = mean_coverage/chromosomal_coverage) %>%
    ## and only consider plasmids.
    filter(SeqType == "plasmid")


## merge with the PIRA estimates, and benchmark copy number estimates.
low.PCN.breseq.PIRA.comparison.df <- low.PCN.breseq.PCN.estimate.df %>%
    left_join(PIRA.alignment.versus.pseudoalignment.comparison.df)

breseq.versus.PIRA.alignment.comparison.plot <- low.PCN.breseq.PIRA.comparison.df %>%
    ggplot(aes(
        x = log10(BreseqCopyNumberEstimate),
        y = log10(minimap2_PIRA_CopyNumberEstimate))) +
    geom_point() +
    theme_classic() +
    xlab("log10(Short Read Breseq Copy Number)")  +
    ylab("log10(Short Read minimap2 alignment Copy Number)") +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(BreseqCopyNumberEstimate), y=log10(minimap2_PIRA_CopyNumberEstimate)),
        color="light blue",
        formula=y~x)

ggsave(
    "../results/breseq-versus-short-read-alignment-low-PCN-estimates.pdf",
    breseq.versus.PIRA.alignment.comparison.plot, height=5, width=5)


breseq.PIRA.alignment.PCN.lm.model <- lm(
    formula = log10(minimap2_PIRA_CopyNumberEstimate) ~ log10(BreseqCopyNumberEstimate),
    data=low.PCN.breseq.PIRA.comparison.df)
## look at the linear regression.
summary(breseq.PIRA.alignment.PCN.lm.model)


breseq.versus.PIRA.pseudoalignment.comparison.plot <- low.PCN.breseq.PIRA.comparison.df %>%
    ggplot(aes(
        x = log10(BreseqCopyNumberEstimate),
        y = log10(PIRA_CopyNumberEstimate))) +
    geom_point() +
    theme_classic() +
    xlab("log10(Short Read Breseq Copy Number)")  +
    ylab("log10(Short Read pseudoalignment Copy Number)") +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(BreseqCopyNumberEstimate), y=log10(PIRA_CopyNumberEstimate)),
        color="light blue",
        formula=y~x)

ggsave(
    "../results/breseq-versus-short-read-pseudoalignment-low-PCN-estimates.pdf",
    breseq.versus.PIRA.pseudoalignment.comparison.plot, height=5, width=5)

breseq.PIRA.pseudoalignment.PCN.lm.model <- lm(
    formula = log10(PIRA_CopyNumberEstimate) ~ log10(BreseqCopyNumberEstimate),
    data=low.PCN.breseq.PIRA.comparison.df)
## look at the linear regression.
summary(breseq.PIRA.pseudoalignment.PCN.lm.model)


#### QUICK CHECK OF PIRA estimates.
log.PIRA.PCN.plot <- PIRA.PCN.estimates %>%
    filter(PIRAReadCount > 10000) %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber))) +
    geom_point(size=0.1,alpha=0.5) +
    theme_classic() +
    xlab("log10(Replicon Length)")  +
    ylab("log10(PIRA Plasmid Copy Number)")
ggsave("../results/log-PIRA-PCN-plot.pdf", log.PIRA.PCN.plot, height=5, width=5)
