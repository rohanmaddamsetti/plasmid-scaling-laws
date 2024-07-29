## NCBI-Hawkey-Bethke-Yao-plasmid-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## plasmid-copy-number-pipeline.py.

##     Make a scatterplot of plasmid copy numbers against plasmid length,
##     and color dots by presence of ARGs.
##    
##     based on these results, plasmids with ARGs do NOT have high copy number.

##    The Hawkey et al. 2022 paper specifically focuses on ESBL resistance,
##    so let's focus on beta-lactamases, and can compare to other kinds of resistances in these data.

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
Hawkey.replicon.length.data <- read.csv("../results/Hawkey2022_replicon_lengths.csv")

## get ARG copy number data.
Hawkey.ARG.copy.number.data <- read.csv("../results/Hawkey2022_ARG_copy_numbers.csv") %>%
    mutate(beta.lactam.resistance = ifelse(str_detect(product,beta.lactam.keywords), TRUE, FALSE))

beta.lactam.ARGs <- filter(Hawkey.ARG.copy.number.data, beta.lactam.resistance==TRUE)
non.beta.lactam.ARGs <- filter(Hawkey.ARG.copy.number.data, beta.lactam.resistance==FALSE)

Hawkey.chromosome.plasmid.copy.number.data <- read.csv("../results/Hawkey2022_chromosome_plasmid_copy_numbers.csv") %>%
    full_join(Hawkey.replicon.length.data) %>%
    mutate(has.ARG = ifelse(SeqID %in% Hawkey.ARG.copy.number.data$SeqID, TRUE, FALSE)) %>%
    mutate(has.beta.lactamase = ifelse(SeqID %in% beta.lactam.ARGs$SeqID, TRUE, FALSE)) %>%
    ## 0 == no ARG, 1 == has ARG, 2 == has beta-lactamase.
    mutate(ARG.classification = has.ARG + has.beta.lactamase) %>%
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification, `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases")) %>%
    ## remove outlier points with very low coverage.
    filter(CopyNumber > 0.5)



## beta-lactamases have higher copy number compared to other ARGs in these strains.
wilcox.test(beta.lactam.ARGs$CopyNumber, non.beta.lactam.ARGs$CopyNumber,alternative="greater")$p.value

mean(beta.lactam.ARGs$CopyNumber)
mean(non.beta.lactam.ARGs$CopyNumber)

median(beta.lactam.ARGs$CopyNumber)
median(non.beta.lactam.ARGs$CopyNumber)

## Plasmids with ARGs actually have lower copy numbers than
## plasmids without ARGs.

Hawkey.plasmid.copy.number.data <- Hawkey.chromosome.plasmid.copy.number.data %>%
    filter(SeqType == "plasmid") %>%
    arrange(CopyNumber)

## write to file to cross-check with ChatGPT generated code.
Hawkey.plasmid.copy.number.data %>%
    rename(Strain = AnnotationAccession) %>%
write.csv(file="../results/Combined.csv")


beta.lactamase.plasmid.data <- Hawkey.plasmid.copy.number.data %>%
    filter(has.beta.lactamase==TRUE)

no.beta.lactamase.plasmid.data <- Hawkey.plasmid.copy.number.data %>%
    filter(has.beta.lactamase == FALSE)

ARG.plasmid.data <- Hawkey.plasmid.copy.number.data %>%
    filter(has.ARG==TRUE)

no.ARG.plasmid.data <- Hawkey.plasmid.copy.number.data %>%
    filter(has.ARG == FALSE)

mean(beta.lactamase.plasmid.data$CopyNumber)
mean(no.beta.lactamase.plasmid.data$CopyNumber)

mean(ARG.plasmid.data$CopyNumber)
mean(no.ARG.plasmid.data$CopyNumber)

median(beta.lactamase.plasmid.data$CopyNumber)
median(no.beta.lactamase.plasmid.data$CopyNumber)

median(ARG.plasmid.data$CopyNumber)
median(no.ARG.plasmid.data$CopyNumber)

## make a linear regression model:
## log10(Plasmid copy number) vs. log10(Plasmid length).
plasmid.lm.model <- lm(
    formula=log10(CopyNumber)~log10(replicon_length),
    data=Hawkey.plasmid.copy.number.data)
## look at the linear regression.
summary(plasmid.lm.model)

second.order.plasmid.lm.model <- lm(
    formula=log10(CopyNumber)~poly(log10(replicon_length),2,raw=TRUE),
    data=Hawkey.plasmid.copy.number.data)
## look at the second order regression.
summary(second.order.plasmid.lm.model)

## let's compare these models. The model with lower AIC is better.
AIC(plasmid.lm.model)
AIC(second.order.plasmid.lm.model)
## also compare the models using an ANOVA.
print(anova(
    plasmid.lm.model,
    second.order.plasmid.lm.model))

## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
Hawkey.plasmid.copy.number.plot <- ggplot(Hawkey.plasmid.copy.number.data,
                                   aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=`Plasmid class`)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=Hawkey.plasmid.copy.number.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(replicon_length),y=log10(CopyNumber)),
        color="light blue",
        formula=y~x) +
    ## let's look at second-order polynomial fit.
    geom_smooth(
        data=Hawkey.plasmid.copy.number.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(replicon_length),y=log10(CopyNumber)),
        color="light green",
        formula=y~poly(x, 2, raw=TRUE))
## save the plot.
ggsave("../results/Hawkey2022-plasmid-copy-number.pdf",
       Hawkey.plasmid.copy.number.plot,height=5.75,width=5.75)

################################################################################
## analyze the plasmid copy number data in the Supporting Information
## of Bethke et al. (2022) and Yao et al. (2022) from the Lingchong You lab.

YouLab.PCN.data <- read.csv("../data/YouLab-PCN-data.csv") %>%
    mutate(replicon_length = Length_in_kbp*1000) %>%
    ## For merging with plasmid.copy.data in the next section.
    mutate(SeqID = Plasmid)

## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
YouLab.PCN.plot <- ggplot(YouLab.PCN.data,
                                   aes(x=log10(replicon_length),y=log10(CopyNumber),
                                       color=Reference)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=YouLab.PCN.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(replicon_length),y=log10(CopyNumber)),
        color="light blue",
        formula=y~x)
## save the plot.
ggsave("../results/YouLab-PCN.pdf",
       YouLab.PCN.plot,height=5.75,width=5.75)

YouLab.PCN.plot2 <- ggplotRegression(
    YouLab.PCN.data,
    "log10(replicon_length)",
    "log10(CopyNumber)")
ggsave("../results/YouLab-PCN-2.pdf",
       YouLab.PCN.plot2,height=5.75,width=7)

################################################################################
## Import the NCBI plasmid copy number data that Maggie generated.
NCBI.plasmid.data <- read.csv("../results/NCBI_plasmid_copy_number_2200.csv") %>%
    mutate(Reference="NCBI") %>%
    ## do these things for compatibility with HBY data.
    mutate(ARG.classification = as.factor(ARG.classification)) %>%
    mutate(`Plasmid class` = recode(ARG.classification, `0` = "No ARGs",
                                    `1` = "Non-beta-lactamase ARGs",
                                    `2` = "Beta-lactamases")) %>%
    ## remove outlier points with very low coverage.
    filter(CopyNumber > 0.5)

## let's make a histogram of PCN in these data.

NCBI.plasmid.histogram <- ggplot(NCBI.plasmid.data, aes(x = log10(CopyNumber))) +
    geom_histogram() +
    theme_classic()

NCBI.plasmid.histogram


################################################################################
## Now combine Hawkey et al. data with You lab data.
Hawkey.Bethke.Yao.data <- Hawkey.plasmid.copy.number.data %>%
    mutate(Reference="Hawkey2022") %>%
    full_join(YouLab.PCN.data)

## Now combine the NCBI data with the HBY data.
NCBI.Hawkey.Bethke.Yao.data <- Hawkey.Bethke.Yao.data %>%
    full_join(NCBI.plasmid.data)


## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
NHBY.PCN.plot <- ggplot(NCBI.Hawkey.Bethke.Yao.data,
                        aes(x=log10(replicon_length),
                            y=log10(CopyNumber),
                            color=Reference)) +
    geom_point(alpha=0.5) +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=NCBI.Hawkey.Bethke.Yao.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(replicon_length),y=log10(CopyNumber)),
        color="light blue",
        formula=y~x)
## save the plot.
ggsave("../results/NHBY-PCN.pdf",
       NHBY.PCN.plot,height=5.75,width=5.75)


## make a linear regression model:
## log10(Plasmid copy number) vs. log10(Plasmid length).
NHBY.lm.model <- lm(
    formula=log10(CopyNumber)~log10(replicon_length),
    data=NCBI.Hawkey.Bethke.Yao.data)
## look at the linear regression.
summary(NHBY.lm.model)

################################################################################
## Let's play with Zhengqing's his notion of plasmid capacity.
NCBI.Hawkey.Bethke.Yao.data2 <- NCBI.Hawkey.Bethke.Yao.data %>%
    mutate(PlasmidCapacity1 = replicon_length*CopyNumber) %>%
    mutate(PlasmidCapacity2 = log10(replicon_length)*log10(CopyNumber))

plasmid.copy.number.plot2 <- ggplot(
    NCBI.Hawkey.Bethke.Yao.data2,
    aes(x=PlasmidCapacity1,
        fill=`Plasmid class`)) +
    geom_histogram() +
    theme_classic() +
    theme(legend.position="top")

plasmid.copy.number.plot3 <- ggplot(
    NCBI.Hawkey.Bethke.Yao.data2,
    aes(x=PlasmidCapacity2,
        fill=`Plasmid class`)) +
    geom_histogram() +
    theme_classic() +
    theme(legend.position="top")

## summarize average plasmid length and copy number (geometric mean)
geom.mean.NCBI.Hawkey.Bethke.Yao.data <- NCBI.Hawkey.Bethke.Yao.data %>%
    group_by(AnnotationAccession) %>%
    summarize(
        geomean_CopyNumber = exp(mean(log(CopyNumber))),
        geomean_replicon_length = exp(mean(log(replicon_length))))
              

geomean.plasmid.copy.number.plot <- ggplot(geom.mean.NCBI.Hawkey.Bethke.Yao.data,
                                   aes(x=log10(geomean_replicon_length),y=log10(geomean_CopyNumber))) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(geometric mean Plasmid copy number)")  +
    xlab("log10(geometric mean Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=geom.mean.NCBI.Hawkey.Bethke.Yao.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(geomean_replicon_length),y=log10(geomean_CopyNumber)),
        color="light blue",
        formula=y~x)
## save the plot.
ggsave("../results/Hawkey2022-geometric-mean-plasmid-copy-number.pdf",
       geomean.plasmid.copy.number.plot,height=5.75,width=5.75)

geomean.plasmid.lm.model <- lm(
    formula=log10(geomean_CopyNumber)~log10(geomean_replicon_length),
    data=geom.mean.NCBI.Hawkey.Bethke.Yao.data)

summary(geomean.plasmid.lm.model)


testplot <- ggplotRegression(
    NCBI.Hawkey.Bethke.Yao.data,
    "log10(replicon_length)", "log10(CopyNumber)")
testplot

testplot2 <- ggplotRegression(
    geom.mean.NCBI.Hawkey.Bethke.Yao.data,
    "log10(geomean_replicon_length)",
    "log10(geomean_CopyNumber)")
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

test3 <- NCBI.Hawkey.Bethke.Yao.data %>%
    mutate(total.plasmid.bp = replicon_length*CopyNumber) %>%
    group_by(AnnotationAccession) %>%
    summarize(
        totalCN = sum(CopyNumber),
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
test4 <- NCBI.Hawkey.Bethke.Yao.data %>%
    mutate(total.plasmid.bp = replicon_length*CopyNumber) %>%
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
## benchmark plasmid copy number estimates with kallisto against  estimates with breseq.
## This benchmarking only runs on the Hawkey data.

## IMPORTANT TODO: figure out what's going on with the NA missing values.

breseq.Hawkey.data <- read.csv("../results/Hawkey2022-breseq-copy-number.csv") %>%
    filter(SeqType == "plasmid")

## let's first plot the trend with these data.
## scatterplot of log10(Plasmid copy number) vs. log10(Plasmid length).
Hawkey.breseq.PCN.plot <- ggplot(breseq.Hawkey.data,
                                 aes(x=log10(replicon_length),
                                     y=log10(CopyNumber))) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept=0,linetype="dashed",color="gray") +
    ylab("log10(Plasmid copy number)")  +
    xlab("log10(Plasmid length in bp)") +
    theme(legend.position="top") +
    ## add the linear regression.
    geom_smooth(
        data=breseq.Hawkey.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(replicon_length),y=log10(CopyNumber)),
        color="light blue",
        formula=y~x)
## save the plot.
ggsave("../results/breseq-Hawkey-PCN.pdf",
       Hawkey.breseq.PCN.plot,height=5.75,width=5.75)

## now let's compare kallisto and breseq PCN estimates.
breseq.plasmid.comp.data <- breseq.Hawkey.data %>%
    rename(BreseqCopyNumber = CopyNumber) %>%
    select(AnnotationAccession,SeqID,SeqType,BreseqCopyNumber)


benchmarking.PCN.estimate.data <- Hawkey.plasmid.copy.number.data %>%
    select(AnnotationAccession,SeqID,SeqType,CopyNumber) %>%
    ## hack to make the SeqIDs compatible with breseq output.
    mutate(SeqID = substr(SeqID, 1,11)) %>%
    full_join(breseq.plasmid.comp.data)


benchmarking.scatterplot <- ggplot(
    benchmarking.PCN.estimate.data,
    aes(x=log10(CopyNumber),y=log10(BreseqCopyNumber))) +
    geom_point() +
    theme_classic() +
    ## add the linear regression.
    geom_smooth(
        data=benchmarking.PCN.estimate.data,
        inherit.aes=FALSE,
        method='lm',
        aes(x=log10(CopyNumber),y=log10(BreseqCopyNumber)),
        color="light blue",
        formula=y~x)
ggsave("../results/Hawkey2022-kallisto-breseq-benchmarking.pdf", benchmarking.scatterplot)

## make a linear regression model:
## log10(kallisto plasmid copy number) vs. log10(Breseq plasmid copy number).
benchmarking.lm.model <- lm(
    formula=log10(BreseqCopyNumber)~log10(CopyNumber),
    data=benchmarking.PCN.estimate.data)
## look at the linear regression.
summary(benchmarking.lm.model)
