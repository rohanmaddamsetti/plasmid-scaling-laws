## PCN-analysis.R by Rohan Maddamsetti.
## analyze the plasmid copy number results made by
## PCN-pipeline.py.

## CRITICAL TODO DURING REVISION:
## REMOVE plasmids that are clearly contigs: those with "unlocalized" or "unplaced"
## in the Definition in the Genbank annotation. Right now I remove all sequences < 1000bp in length,
## which is probably good enough.

## CRITICAL TODO DURING REVISION:
## There are 483 genomes that DO NOT have ecological annotation-- these
## are genomes with "Chromosome" level assembly in the PCN pipeline,
## which are not found in the "Complete Genome" annotations used in the main
## ecological annotation pipeline. This inconsistency between the sets of genomes
## used across the two analyses is causing this issue. This does not affect any fundamental
## analyses here, but this inconsistency is extremely annoying and needs to be fixed.

## One simple choice would be to just remove all the chromosome-level assemblies,
## and stick to 'Complete Genome' level assemblies throughout. This is probably the cleanest solution.
## This might also take care of small plasmid contigs,
## but I should add an extra filtering step to remove these regardless.

## TODO WHEN RERUNNING FROM SCRATCH:
## Make sure all the genomes in ../results/gbk-annotation are consistent
## with the genomes annotated in computationally-annotated-genomes etc.
## (that is, there are no suppressed RefSeq genomes in this folder, and everything is annotated.)
## in the Annotation and SeqType columns, and rewrite upstream code to solve this problem.


## IMPORTANT: In Figures 2 and 3,
## I annotate megaplasmid/chromids to differentiate them from plasmids, based on size.
## The metabolic scaling law emerges in chromids.
## From Coluzzi et al. (2022): To avoid the misidentification of ICEs as
## conjugative plasmids in chromids or secondary chromosomes,
## we excluded from further study the 419 plasmids larger than 500 kb.

## IMPORTANT: my analysis does not distinguish between chromids and megaplasmids,
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
## Global variables and functions.

## require that PCN estimates are supported by a minimum of MIN_READ_COUNT reads per replicon.
MIN_READ_COUNT <- 10000

## To avoid the misidentification of ICEs as conjugative plasmids in chromids or secondary chromosomes,
## specifically label megaplasmids or chromids larger than 500 kb (Cite Coluzzi et al. 2022 for this practice).
PLASMID_LENGTH_THRESHOLD <- 500000

## Threshold for showing the emergent metabolic scaling relation.
METABOLIC_GENE_THRESHOLD <- 100


make.PIRA.vs.naive.themisto.plot <- function(PIRA.vs.naive.themisto.df, show.linear.regression=TRUE) {
    figure.panel <- PIRA.vs.naive.themisto.df %>%
        ggplot(aes(
            x = log10(ThemistoNaiveCopyNumber),
            y = log10(PIRACopyNumber),
            color = PIRA_low_PCN,
            shape = InsufficientReads)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        scale_color_manual(values=c("black", "red")) +
        theme_classic() +
        xlab("log10(direct themisto PCN)")  +
        ylab("log10(PIRA PCN)") +
        guides(color = 'none', shape = 'none')

    if (show.linear.regression) {
        figure.panel <- figure.panel +
            ## add the linear regression.
            geom_smooth(
                method='lm',
                aes(x=log10(ThemistoNaiveCopyNumber), y=log10(PIRACopyNumber)),
                color="light blue",
                formula=y~x)
    }
    
    return(figure.panel)
}


filter.and.group.together.smaller.groups.in.the.correlate.column <- function(df, correlate_column_name_string, lumped_group_name, min_group_size = 50) {    
    ## this function filters data frames for groups with more than min_group_size data points in the column
    ## named in the string correlate_column_name_string, using tidy evaluation,
    ## and groups together data points in groups that fall below the min_group_size into lumped_group_name.
    
    correlate_column_name <- sym(correlate_column_name_string)
    
    correlate.groups <- df %>%
        filter(!is.na(!!correlate_column_name)) %>%
        filter(!!correlate_column_name != '-') %>%
        count(!!correlate_column_name) %>%
        filter(n >= min_group_size)
    
    data_in_nice_correlate_groups_df <- df %>%
        filter(!!correlate_column_name %in% correlate.groups[[correlate_column_name_string]])

    remaining_data_df <- df %>%
        ## negate the previous filter,
        filter(!(!!correlate_column_name %in% correlate.groups[[correlate_column_name_string]])) %>%
        ## but be sure to remove any NA or unannotated values in the correlate column.
        filter(!is.na(!!correlate_column_name)) %>%
        filter(!!correlate_column_name != '-') %>%
        ## and lump these point together into the whatever the lumped_group_name is.
        mutate(!!correlate_column_name := lumped_group_name)

    final_df <- bind_rows(data_in_nice_correlate_groups_df, remaining_data_df)

    ## put the group named by lumped_group_name as the final group.
    group_values_without_lumped_group <- sort(setdiff(unique(final_df[[correlate_column_name_string]]), lumped_group_name))
    group_levels <- c(group_values_without_lumped_group, lumped_group_name)

    final_df <- final_df %>%
        mutate(!!correlate_column_name := factor(!!correlate_column_name, levels = group_levels))

    return(final_df)
}


filter.correlate.column <- function(df, correlate_column_name_string, min_group_size = 50) {
    ## this function filters data frames for groups with more than min_group_size data points in the column
    ## named in the string correlate_column_name_string, using tidy evaluation.
    ## I figured this out using ChatGPT and Chapter 20 of Advanced R by Hadley Wickham
    ## as a reference: https://adv-r.hadley.nz/evaluation.html#tidy-evaluation
    
    correlate_column_name <- sym(correlate_column_name_string)
    
    correlate.groups <- df %>%
        filter(!is.na(!!correlate_column_name)) %>%
        filter(!!correlate_column_name != '-') %>%
        count(!!correlate_column_name) %>%
        filter(n >= min_group_size)

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


add_mean_PTU_length_column <- function(df, PTU_column_name_string) {
    ## the relevant PTU column is named in the string correlate_column_name_string,
    ## and parsed using tidy evaluation.
    PTU_column_name <- sym(PTU_column_name_string)

    mean.PTU.length.df <- df %>%
        group_by(!!PTU_column_name) %>%
        summarize(mean_PTU_length = mean(replicon_length))

    ## join the mean PTU length column to the original dataframe.
    df %>%
        inner_join(mean.PTU.length.df)
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
    ## add a column for the plasmid clusters by size
    clustered.PIRA.PCN.estimates$Size_Cluster <- cluster_assignment$.cluster
    ## return the clustered dataframe.
    return(clustered.PIRA.PCN.estimates)
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


remove.genomes.with.bad.chromosomes <- function(full.PIRA.estimates) {
    ## remove genomes where the longest replicon is not the chromosome.
    ## this can be caused by either misannotation of the chromosome/plasmid,
    ## or by genomes that didn't get any CopyNumber estimate for the chromosome.

    genomes.with.plasmids.that.have.max.replicon.lengths.df <- full.PIRA.estimates %>%
        filter(replicon_length == max_replicon_length) %>%
        filter(SeqType == "plasmid")

    bad.genomes.with.plasmids.that.have.max.replicon.lengths <- unique(
        genomes.with.plasmids.that.have.max.replicon.lengths.df$AnnotationAccession
    )
    ## remove these bad genomes from the analysis.
    full.PIRA.estimates %>%
        filter(!(AnnotationAccession %in% bad.genomes.with.plasmids.that.have.max.replicon.lengths))
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
        xlab("proportion of plasmids") +
        theme_classic() +
        ggtitle(title) +
        ## plot CIs.
        geom_errorbarh(aes(xmin=Left,xmax=Right), height=0.2, size=0.2) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
    
    if (no.category.label)
        Fig.panel <- Fig.panel +
            theme(axis.text.y=element_blank())
    
    return(Fig.panel)
}


make_PCN_base_plot <- function(my.PCN.data) {
    ## Make the basic plot for S2Fig, before adding the marginal histograms,
    ## or facetting by column
    my.PCN.data %>%
        ggplot(aes(
            x = log10_replicon_length,
            y = log10_PIRACopyNumber,
            color = PredictedMobility)) +
        geom_point(size=0.5,alpha=0.8) +
        geom_hline(yintercept=0,linetype="dashed",color="gray") +
        theme_classic() +
        scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
        ## make the points in the legend larger.
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("log10(length)")  +
        ylab("log10(copy number)") +
        theme(legend.position = "bottom") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_normalized_PCN_base_plot <- function(my.PCN.data) {
    ## Make the basic plot for Fig1B, before adding the marginal histograms,
    ## or facetting by column
    my.PCN.data %>%
        ggplot(aes(
            x = log10_normalized_replicon_length,
            y = log10_PIRACopyNumber,
            color = PredictedMobility)) +
        geom_point(size=0.5,alpha=0.8) +
        geom_hline(yintercept=0,linetype="dashed",color="gray") +
        theme_classic() +
        scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
        ## make the points in the legend larger.
        guides(color = guide_legend(override.aes = list(size = 5))) +
        xlab("log10(normalized length)")  +
        ylab("log10(copy number)") +
        theme(legend.position = "bottom") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_PTU_length_rank_plot <- function(PTUs.with.PIRA.PCN.estimates) {
    PTUs.with.PIRA.PCN.estimates %>%
        ggplot(aes(
            x = rank,
            y = log10(replicon_length),
            color = PredictedMobility)) +
        geom_point(size=0.5,alpha=0.8) +
        theme_cowplot() + ## for 7pt margins
        scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
        xlab("PTUs")  +
        ylab("log10(length)") +
        guides(color = "none") +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_PTU_PCN_rank_plot <- function(PTUs.with.PIRA.PCN.estimates) {
    PTUs.with.PIRA.PCN.estimates %>%
        ggplot(aes(
            x = rank,
            y = log10(PIRACopyNumber),
            color = PredictedMobility)) +
        geom_point(size=0.5,alpha=0.8) +
        theme_cowplot() + ## for 7pt margins
        scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
        xlab("PTUs")  +
        ylab("log10(copy number)") +
        guides(color = "none") +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_PTU_mean_length_PCN_base_plot <- function(df_with_PCN_mean_lengths) {
    ## Make the basic plot for S4Fig panels, before adding the marginal histograms or anything else.
    my_base_plot <- df_with_PCN_mean_lengths %>%
        ggplot(aes(
            x = log10(mean_PTU_length),
            y = log10(PIRACopyNumber),
            color = PredictedMobility)) +
        geom_point(size=0.5,alpha=0.8) +
        geom_hline(yintercept=0,linetype="dashed",color="gray") +
        theme_cowplot() + ## for 7pt margins
        scale_color_manual(values=c("#fc8d62","#66c2a5","#8da0cb"), name="plasmid mobility") +
        xlab("log10(PTU mean length)")  +
        ylab("log10(copy number)") +
        guides(color = "none") +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))

    plot_with_marginals <- ggExtra::ggMarginal(my_base_plot, groupColour = TRUE, groupFill = TRUE, margins="both")

    return (plot_with_marginals)
}


make_CDS_scaling_base_plot <- function(CDS.fraction.data) {
    CDS.fraction.data %>%
        ggplot(
            aes(
                x = log10(SeqLength),
                y = log10(CDS_length),
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.8) +
        xlab("log10(length)") +
        ylab("log10(coding sequence length)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_normalized_CDS_base_plot <- function(CDS.fraction.data) {
    CDS.fraction.data %>%
        ggplot(
            aes(
                x = log10(normalized_replicon_length),
                y = log10(CDS_length),
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.8) +
        xlab("log10(length normalized by chromosome)") +
        ylab("log10(coding sequence length)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


make_metabolic_scaling_base_plot <- function(metabolic.gene.plasmid.and.chromosome.data) {

    ## fit a linear regression to the chromosome data, and save the fit,
    ## to show the emergent scaling trend on the figures.
    ## Using chromosome here to preserve independence between data
    ## and the predicted scaling.
    chromosome.metabolic.scaling.model <- lm(
        formula = log10_metabolic_protein_count ~ log10_replicon_length,
        data = filter(
            metabolic.gene.plasmid.and.chromosome.data,
            SeqType == "chromosome"))
        
    ## save the model fit as a dataframe for plotting.
    chromosome.metabolic.scaling.fit.df <- data.frame(
        log10_replicon_length =  filter(
            metabolic.gene.plasmid.and.chromosome.data,
            SeqType== "chromosome")$log10_replicon_length,
        log10_metabolic_protein_count = chromosome.metabolic.scaling.model$fit
    )

   
    metabolic.gene.plasmid.and.chromosome.data %>%
        ggplot(
            aes(
                x = log10_replicon_length,
                y = log10_metabolic_protein_count,
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.8) +
        xlab("log10(length)") +
        ylab("log10(metabolic genes)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11)) +
        geom_line(data = chromosome.metabolic.scaling.fit.df,
                  color = 'black', linetype = "solid", size=0.5)
}


make_normalized_metabolic_scaling_base_plot <- function(metabolic.gene.plasmid.and.chromosome.data) {
    metabolic.gene.plasmid.and.chromosome.data %>%
        ggplot(
            aes(
                x = log10(normalized_replicon_length),
                y = log10_metabolic_protein_count,
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.8) +
        xlab("log10(length normalized by chromosome)") +
        ylab("log10(metabolic genes)") +
        theme_classic() +
        guides(color = "none") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))
}


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
    mutate(Annotation = replace(Annotation, Annotation == "Animals", "Animals"))

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

######### TODO: in the upstream pipeline for the GhostKOALA
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

##  get output from calculate-CDS-rRNA-fractions.py.
## We want to answer a basic question that Lingchong asked:
## for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?
CDS.rRNA.fraction.data <- read.csv("../results/CDS-rRNA-fractions.csv") %>%
    ## make the dataframe compatible with replicon.annotation.data,
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(replicon.annotation.data) %>%
    ## add a column for nomalized plasmid lengths.
    normalize.plasmid.lengths()

################################################################################
## Get the PCN estimates. These are the main data for this paper.
## IMPORTANT NOTE: We only have PCN estimates for ~10,000 plasmids in ~4,500 genomes,
## for which we can find linked short-read data in the NCBI Sequencing Read Archive (SRA).
## This is a subset of the genomes and plasmids considered in this paper.


## IMPORTANT: do NOT filter this for just plasmids just yet--
## we need to include chromosomes for proper comparison with breseq results.

## IMPORTANT: PIRA was only run on the genomes that had multireads.
## So, we need to combine the PIRA estimates with the subset of
## naive themisto estimates for which no multireads were called.
PIRA.estimates <- read.csv("../results/PIRA-PCN-estimates.csv") %>%
    ## SMALL TODO: In PCN_pipeline.py, make sure the ThemistoID_right column is dropped
## and not written out to this CSV file.
    select(-ThemistoID_right) %>%
    ## get rid of unneeded columns.
    select(-ThemistoID) %>%
    rename(CopyNumber = PIRA_CopyNumberEstimate)


## now get the data for the genomes without multireads.
no.multiread.naive.themisto.estimates <- read.csv("../results/naive-themisto-PCN-estimates.csv") %>%
    ## SMALL TODO: In PCN_pipeline.py, make sure the ThemistoID_right column is dropped
    ## and not written out to this CSV file.
    select(-SeqType_right) %>%
    ## filter for the genomes without multireads.
    filter(!(AnnotationAccession %in% PIRA.estimates$AnnotationAccession)) %>%
    ## no additional read counts in these data, add these columns for compatibility with PIRA.estimates.
    mutate(InitialReadCount = ReadCount) %>%
    mutate(AdditionalReadCount = 0) %>%
    mutate(InitialCopyNumberEstimate =  CopyNumber)

  
## now merge the datasets together.
## CRITICAL TODO: fix upstream annotation so that we don't have any
## "NA" or "blank" Annotation genomes in full.PIRA.estimates.
full.PIRA.estimates <- full_join(no.multiread.naive.themisto.estimates, PIRA.estimates) %>%
    ## rename columns to so that we can compare estimates for benchmarking.
    rename(
        PIRACopyNumber = CopyNumber,
        PIRAReadCount = ReadCount,
        PIRASequencingCoverage = SequencingCoverage,
        PIRALongestRepliconCoverage = LongestRepliconCoverage) %>%
    filter(PIRAReadCount > MIN_READ_COUNT) %>%
    mutate(RepliconDNAContent = replicon_length * PIRACopyNumber) %>%
    ## add taxonomic and ecological annotation.
    left_join(replicon.annotation.data) %>%
    normalize.plasmid.lengths() %>%
    ## remove any genomes that either don't have a chromosome with a copy number
    ## estimate, or that have a 'plasmid' that has the max_replicon_length in the genome.
    remove.genomes.with.bad.chromosomes() %>%
    ## Remove outlier 'plasmids' < 1000bp, these are likely small, unassembled contigs.
    ## CRITICAL TODO: REMOVE plasmids that are clearly contigs: those with "unlocalized" or "unplaced"
    ## in the Definition in the Genbank annotation. Right now I remove all sequences < 1000bp in length,
    ## which is probably good enough.
    filter(replicon_length > 1000)

## 4,512 genomes in the final PCN dataset.
length(unique(full.PIRA.estimates$AnnotationAccession))
## 11,338 plasmids in the final PCN dataset.
nrow(filter(full.PIRA.estimates, SeqType=="plasmid"))

## write the normalized data to disk.
write.csv(full.PIRA.estimates, "../results/S2Data-PIRA-PCN-estimates-with-normalization.csv")


## TODO: fix upstream annotation so I don't have to do this filtering to exclude NA Annotations.
## TODO: FIGURE OUT WHY 1946 GENOMES ARE NOT ANNOTATED RIGHT.
## This is actually not critical for the PCN analysis, but would be nice to have this fixed.
unannotated.full.PIRA.estimates <- full.PIRA.estimates %>%
    filter(is.na(Annotation) | (Annotation == "blank") | Annotation == "NA")


################################################################################
## Genomes for Supplementary Table 1.
## get concrete statistics for how many genome papers report plasmid copy number:

## sample 50 genomes with plasmids at random, and examine how many papers report PCN.
## pick genomes with some plasmid with PCN > 10, so that we are focusing on genomes where people
## might actually want to report PCN.
## for reproducibility, use random seed = 60, which I generated using random.org (see screenshot in ../data).

genomes.to.sample.from.PIRA.PCN.estimates <- full.PIRA.estimates %>%
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
## Supplementary Figure S1: Plasmid copy number benchmarking.
################################################################################
## S1AFig and S1BFig: Effect of PIRA compared to naive themisto estimates.

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
    left_join(full.PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))

## 1,563 plasmids have sufficient reads with PIRA, but not with the naive themisto read mapping.
sum(PIRA.vs.naive.themisto.df$InsufficientReads == TRUE)

## In S1 Figure panel A, show points with insufficient reads, and remove the linear regression.
## ## but let's include them to show how PIRA recovers PCN for more plasmids.
S1FigA <- make.PIRA.vs.naive.themisto.plot(PIRA.vs.naive.themisto.df, FALSE) +
    ggtitle("PIRA recovers more plasmids\nby including multiread data")

## In S1 Figure panel B, remove points with insufficient reads,
S1FigB <- PIRA.vs.naive.themisto.df %>%
    filter(InsufficientReads == FALSE) %>%
    make.PIRA.vs.naive.themisto.plot()

###################################################################################
## Supplementary Figure S1C
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
    left_join(full.PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))

## make S1 Figure panel C
S1FigC <- PIRA.vs.minimap2.df %>%
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
## Supplementary Figure S1D.
## Benchmarking of these 100 random genomes with breseq as another gold standard control,
## This additional test makes sure these estimates are accurate,
## and not artifactual due to low sequencing coverage with minimap2 compared to breseq.
## probably due to stringent minimap2 parameters by default, for now will not explore or discuss in text.

## first get the metadata we need from the PIRA estimates.
PCN.benchmark.metadata.df <- full.PIRA.estimates %>%
    select("AnnotationAccession", "SeqID", "SeqType", "replicon_length") %>%
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
    left_join(full.PIRA.estimates) %>%
    ## color points with PIRA PCN < 0.8
    mutate(PIRA_low_PCN = ifelse(PIRACopyNumber< 0.8, TRUE, FALSE))


## Now make S1 Figure panel D.
S1FigD <- PIRA.vs.breseq.df %>%
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
## Supplementary Figure S1E:
## compare kallisto to naive themisto PCN estimates to show that PCN numbers by pseudoalignment
## are reproducible irrespective of the specific software implementation.

kallisto.replicon.PCN.estimates <- read.csv("../results/kallisto-replicon_copy_numbers.csv") %>%
    rename(KallistoNaiveCopyNumber = CopyNumber) %>%
    filter(SeqType == "plasmid") %>%
    left_join(PCN.replicon.metadata)

## compare naive themisto to naive kallisto estimates.
kallisto.vs.naive.themisto.df <- naive.themisto.PCN.estimates %>%
    ## IMPORTANT: remove plasmids with insufficient reads.
    filter(InsufficientReads == FALSE) %>%
    left_join(kallisto.replicon.PCN.estimates) %>%
    filter(SeqType == "plasmid")

## make Supplementary Figure S1E.
S1FigE <- kallisto.vs.naive.themisto.df %>%
    ggplot(aes(x = log10(ThemistoNaiveCopyNumber), y = log10(KallistoNaiveCopyNumber))) +
    geom_point(size=1) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_classic() +
    xlab("log10(direct themisto PCN)") +
    ylab("log10(kallisto PCN)")  +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(ThemistoNaiveCopyNumber), y=log10(KallistoNaiveCopyNumber)),
        color="light blue",
        formula=y~x)

################################################################################
## Supplementary Figure 1F. Benchmarking against external PCN estimates from Shaw et al. (2021)
## in Sciences Advances.


Shaw2021_NCBI_metadata <- read.csv("../results/NCBI_REHAB_plasmid_metadata.csv")

## Get PIRA estimates for plasmids in Shaw et al. (2021).
Shaw2021.PIRA.estimates <- Shaw2021_NCBI_metadata %>%
    ## trim the initial "p" in SeqName column for the merge.
    mutate(SeqName = str_replace(SeqName, "^p", "")) %>%
    left_join(full.PIRA.estimates)

## Get the PCN estimates from Shaw et al. (2021).
Shaw2021.PCN.estimates <- read.csv("../data/Shaw2021_PCN_data.csv") %>%
    select(Plasmid.name, Parent.isolate, Total.length..bp., Depth..relative.to.chromosomal.coverage.) %>%
    rename(SeqName = Plasmid.name) %>%
    rename(Strain = Parent.isolate) %>%
    rename(replicon_length = Total.length..bp.) %>%
    rename(Shaw2021_PCN = Depth..relative.to.chromosomal.coverage.)


## merge PCN estimates to compare.
PIRA.vs.ShawPCN.df <- Shaw2021.PIRA.estimates %>%
    inner_join(Shaw2021.PCN.estimates)

S1FigF <- PIRA.vs.ShawPCN.df %>%
    ggplot(aes(
        x = log10(Shaw2021_PCN),
        y = log10(PIRACopyNumber)
    )) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    xlab("log10(Shaw et al. (2021) PCN)")  +
    ylab("log10(PIRA PCN)") +
    ## add the linear regression.
    geom_smooth(
        method='lm',
        aes(x=log10(Shaw2021_PCN), y=log10(PIRACopyNumber)),
        color="light blue",
        formula=y~x) +
    ggtitle("PIRA recapitulates PCN estimates\nin Shaw et al. (2021)\nSupplementary Table S2")


S1Fig <- plot_grid(
    S1FigA, S1FigB, S1FigC, S1FigD, S1FigE, S1FigF,
    labels = c("A", "B", "C", "D", "E", "F"),
    nrow=3)

## save Supplementary Figure S1.
ggsave("../results/S1Fig.pdf", S1Fig, height=10, width=7)

################################################################################
## PLASMID BIOLOGY ANALYSIS
################################################################################
## take the PIRA estimates, filter for plasmids,
## and annotate plasmids with ARGs, and MOB type.
PIRA.PCN.estimates <- full.PIRA.estimates %>%
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

################################################################################
## make Supplementary Data File 3 with means, confidence intervals,
## and quantile statistics (including minima and maxima PCN per bin)
## for plasmids binned every 1000bp in length.

binned.PIRA.PCN.estimate.summary <- PIRA.PCN.estimates %>%
    arrange(replicon_length) %>%
    mutate(replicon_length_percentile = ntile(replicon_length, 100)) %>%
    group_by(replicon_length_percentile) %>%
    summarize(
        mean_replicon_length = mean(replicon_length),
        mean_normalized_replicon_length = mean(normalized_replicon_length),
        mean_log10_replicon_length = mean(log10(replicon_length)),
        mean_log10_normalized_replicon_length = mean(log10(normalized_replicon_length)),
        mean_PCN = mean(PIRACopyNumber),
        mean_log10_PCN = mean(log10(PIRACopyNumber)),
        min_PCN = min(PIRACopyNumber),
        max_PCN = max(PIRACopyNumber),
        std_PCN = sd(PIRACopyNumber),
        std_log10_PCN = sd(log10(PIRACopyNumber)),
        count = n(),
        ## calculate quantiles
        q25_PCN = quantile(PIRACopyNumber, 0.25),
        q50_PCN = quantile(PIRACopyNumber, 0.50),  ## Median
        q75_PCN = quantile(PIRACopyNumber, 0.75),
        q25_log10_PCN = quantile(log10(PIRACopyNumber), 0.25),
        q50_log10_PCN = quantile(log10(PIRACopyNumber), 0.50),  ## Median
        q75_log10_PCN = quantile(log10(PIRACopyNumber), 0.75),
        ## normal approximation
        CI_Lower_PCN = mean_PCN - 1.96 * (std_PCN / sqrt(count)),
        CI_Upper_PCN = mean_PCN - 1.96 * (std_PCN / sqrt(count)),
        CI_Lower_log10_PCN = mean_log10_PCN - 1.96 * (std_log10_PCN / sqrt(count)),
        CI_Upper_log10_PCN = mean_log10_PCN + 1.96 * (std_log10_PCN / sqrt(count))
    )

## write this summary to disk.
write.csv(binned.PIRA.PCN.estimate.summary, "../results/S3Data-PIRA-PCN-Kbp-bin-summary.csv")

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

## the breakpoint is at 4.75. 10^4.75 = 56,234bp.
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

## the breakpoint is at -1.761. 10^-1.761 = 1.73% of the length of the chromosome.
summary(segmented.normalized.PCN.model)

## save the segmented regression fit as a dataframe.
normalized.segmented.fit.df = data.frame(
    log10_normalized_replicon_length = PIRA.PCN.estimates$log10_normalized_replicon_length,
    log10_PIRACopyNumber = broken.line(segmented.normalized.PCN.model)$fit)

## compare to a second-order polynomial fit.
second.order.normalized.PCN.lm.model <- lm(
    formula=log10_PIRACopyNumber ~ poly(log10_normalized_replicon_length,2, raw=TRUE),
    data=PIRA.PCN.estimates)


## let's compare these models. The model with lower AIC is better.
## Again, the normalized segmented model has the best fit.
AIC(normalized.PCN.lm.model)
AIC(second.order.normalized.PCN.lm.model)
AIC(segmented.normalized.PCN.model)

## get R^2 values for these models
print("LINEAR MODEL WITH NORMALIZED LENGTH")
summary(normalized.PCN.lm.model)
print("QUADRATIC MODEL WITH NORMALIZED LENGTH")
summary(second.order.normalized.PCN.lm.model)
print("SEGMENTED MODEL WITH NORMALIZED LENGTH")
summary(segmented.normalized.PCN.model)


################################################################################
## Figure 1BC and Supplementary Figure S2.
## Plasmid copy number pipeline and inverse correlation between plasmid length and copy number.

## Supplementary Figure S3 shows the unnormalized anticorrelation.
## Break down this result by predicted plasmid mobility.
S2FigA_base <- PIRA.PCN.estimates %>%
    filter(!is.na(Annotation)) %>%
    filter.correlate.column("PredictedMobility") %>%
    make_PCN_base_plot() +
    theme(strip.background = element_blank())

## Get the legend.
S2Fig_legend <- get_legend(S2FigA_base)

## draw the segmented regression,
## and remove the legend.
S2FigA_without_marginals <- S2FigA_base +
    geom_line(data = segmented.fit.df, color = 'maroon') +
    guides(color = "none")

## Add the marginal histograms
S2FigA <- ggExtra::ggMarginal(S2FigA_without_marginals, groupColour = TRUE, groupFill = TRUE, margins="both") 

## S2Fig panel B: facet by ecological annotation.
S2FigB <- S2FigA_base + guides(color = "none") + facet_wrap(.~Annotation)

## make the panels without the legend
S2FigAB <- plot_grid(S2FigA, S2FigB, labels=c('A', 'B'), ncol=2, rel_widths = c(1, 1))

## make S2 Figure with the legend and save to file.
S2Fig <- plot_grid(S2FigAB, S2Fig_legend, ncol=1, rel_heights = c(1,0.1))
ggsave("../results/S2Fig.pdf", S2Fig, height=4.25, width=7.25)


## Figure 1BC.
##  plot normalized plasmid length relative to the length of the longest
## chromosome in the cell.
## This figure show that the scaling law holds well, until plasmids reach ~1.73% of the main chromosome in length.
## then, copy number is roughly flat.

## scatterplot of log10(Normalized plasmid copy number) vs. log10(plasmid length).
Fig1B_base <- PIRA.PCN.estimates %>%
    filter(!is.na(Annotation)) %>%
    filter.correlate.column("PredictedMobility") %>%
    make_normalized_PCN_base_plot() +
     theme(strip.background = element_blank())

## Get the legend.
Fig1BC_legend <- get_legend(Fig1B_base)

## draw the segmented regression,
## and remove the legend.
Fig1B_without_marginals <- Fig1B_base +
    geom_line(data = normalized.segmented.fit.df, color = 'maroon') +
    guides(color = "none")

## add the marginal histogram to Figure 1B.
Fig1B <- ggExtra::ggMarginal(Fig1B_without_marginals, groupColour = TRUE, groupFill = TRUE, margins="both")

## Figure 1C: facet by ecological annotation.
Fig1C <- Fig1B_base + guides(color = "none") + facet_wrap(.~Annotation)

## make Fig1BC and save to file.
Fig1BC_title <- ggdraw() + draw_label("Plasmid length normalized by chromosome", fontface='bold')
Fig1BC <- plot_grid(Fig1B, Fig1C, labels=c('B', 'C'), ncol=2, rel_widths = c(1, 1))
Fig1BC_with_title_and_legend <- plot_grid(Fig1BC_title, Fig1BC, Fig1BC_legend, ncol=1, rel_heights = c(0.1,1,0.1))

ggsave("../results/Fig1BC.pdf", Fig1BC_with_title_and_legend, height=4.5, width=7.25)


## Show summary statistics on a version of Figure 1B for Reviewer 1.
Fig1B.with.confint <- Fig1B_without_marginals +
    geom_line(data = binned.PIRA.PCN.estimate.summary,
        aes(
            x = mean_log10_normalized_replicon_length,
            y = mean_log10_PCN),
        color='red') +
        geom_line(data = binned.PIRA.PCN.estimate.summary,
        aes(
            x = mean_log10_normalized_replicon_length,
            y = CI_Upper_log10_PCN),
        color='black') +
    geom_line(data = binned.PIRA.PCN.estimate.summary,
        aes(
            x = mean_log10_normalized_replicon_length,
            y = CI_Lower_log10_PCN),
        color='black') +
    geom_line(data = binned.PIRA.PCN.estimate.summary,
        aes(
            x = mean_log10_normalized_replicon_length,
            y = q25_log10_PCN),
        color='blue') +
    geom_line(data = binned.PIRA.PCN.estimate.summary,
        aes(
            x = mean_log10_normalized_replicon_length,
            y = q75_log10_PCN),
        color='blue')

## make the figure for Reviewer 1.
ggsave("../results/Fig1B-for-reviewer-1.pdf", Fig1B.with.confint, height=4.5, width=4.5)

################################################################################
## calculate basic statistics about the clusters of small and large plasmids.

small.plasmids <- PIRA.PCN.estimates %>%
    filter(Size_Cluster == "Cluster_2")
## mean length of small plasmids is 6433 bp.
mean(small.plasmids$replicon_length)
## mean PCN of small plasmids is 28.4.
mean(small.plasmids$PIRACopyNumber)


large.plasmids <- PIRA.PCN.estimates %>%
    filter(Size_Cluster == "Cluster_1")
## mean length of large plasmids is 137704 bp.
mean(large.plasmids$replicon_length)
## mean PCN of large plasmids is 1.79.
mean(large.plasmids$PIRACopyNumber)

## examine the tail of very large plasmids that are longer than 500Kbp.
very.large.plasmids <- large.plasmids %>%
    filter(replicon_length > 500000)

very.large.plasmids.by.mobility <- very.large.plasmids %>%
    count(PredictedMobility)
## Very large plasmids: these are chromids.
##       conjugative  29
##       mobilizable  28
##   non-mobilizable 101
##              <NA>  31

## get the percentages / numbers of predicted conjugative plasmids in the large.plasmids cluster.
conjugative.plasmids <- PIRA.PCN.estimates %>%
    filter(PredictedMobility == "conjugative")

## 3126 conjugative plasmids
nrow(conjugative.plasmids)
## 3103 conjugative plasmids out of 3126 (99.3% are in the large cluster).
nrow(filter(conjugative.plasmids, Size_Cluster == "Cluster_1"))

################################################################################
## Supplementary Figures S3. Break down the result in Figure 1 
## by genus to show universality of the PCN vs. length anticorrelation.

## Supplementary Figure S3 
S3FigA <- PIRA.PCN.estimates %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column(
        "TaxonomicGroup", "All other\ntaxonomic groups") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ TaxonomicGroup, nrow = 1) +
    theme(strip.background = element_blank()) +
    guides(color = "none") +
    ggtitle("Taxonomic groups")


S3FigB <- PIRA.PCN.estimates %>%
    ## put new lines after slasks to improve the aspect ratio of the subpanels.
    mutate(TaxonomicSubgroup = str_replace_all(TaxonomicSubgroup, "/", "/\n")) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column(
        "TaxonomicSubgroup", "All other\ntaxonomic subgroups") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup, ncol = 5) +
    theme(strip.background = element_blank()) +
    guides(color = "none") +
    ggtitle("Taxonomic subgroups")

## Break down by genus.
S3FigC <- PIRA.PCN.estimates %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("Genus", "All other genera") %>%
    make_PCN_base_plot() +
    facet_wrap(. ~ Genus, ncol = 7) +
    theme(strip.background = element_blank()) +
    guides(color = "none") +
    ggtitle("Microbial genera")

S3Fig <- plot_grid(S3FigA, S3FigB, S3FigC, labels=c('A','B','C'), ncol = 1, rel_heights=c(1.5,2.5,3.5))

## save the plot.
ggsave("../results/S3Fig.pdf", S3Fig, height=11, width=7.5, limitsize=FALSE)


################################################################################
## The PCN vs. plasmid length anticorrelation holds within individual genomes.

analyze.within.genome.correlations <- function(MIN_PLASMIDS_PER_GENOME=2, df=PIRA.PCN.estimates) {
    within.genome.correlation.data.df <- df %>%
        group_by(AnnotationAccession) %>%
        summarize(
            plasmids_within_genome = n(),
            correlation = cor(log10(replicon_length), log10(PIRACopyNumber))) %>%
        ## make sure there is at least MIN_PLASMIDS_PER_GENOME
        filter(plasmids_within_genome >= MIN_PLASMIDS_PER_GENOME)

    negative.within.genome.correlation.data.df <- within.genome.correlation.data.df %>%
        filter(correlation < 0)

    print("the number of genomes with an inverse correlation")
    print(nrow(negative.within.genome.correlation.data.df))
    print("the average inverse correlation for these genomes:")
    print(mean(negative.within.genome.correlation.data.df$correlation))

    positive.within.genome.correlation.data.df <- within.genome.correlation.data.df %>%
        filter(correlation > 0)
    print("the number of genomes with an positive correlation")
    print(nrow(positive.within.genome.correlation.data.df))
    print("the average correlation for these genomes:")
    print(mean(positive.within.genome.correlation.data.df$correlation))

    print("")
}

analyze.within.genome.correlations(2)

analyze.within.genome.correlations(3)

################################################################################
## Small multicopy plasmid almost always coexist with large low-copy plasmids.
## look at the joint distribution of largest plasmid per cell and smallest plasmid per cell.
## Does this give any insight into how small plasmids persist?
## there is data suggesting that plasmids tend to live as "packs" in cells.
## This is really quite an interesting result.
## relevant recent paper: "A cryptic plasmid is among the most numerous genetic elements in the human gut"
## in Cell.


## Group the data frame by AnnotationAccession and calculate the maximum replicon length in each group
max.plasmid.lengths <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession, Annotation) %>%
    summarise(max_replicon_length = max(replicon_length))

min.plasmid.lengths <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession, Annotation) %>%
    summarise(min_replicon_length = min(replicon_length))

max.PCNs <- PIRA.PCN.estimates %>%
group_by(AnnotationAccession, Annotation) %>%
    summarise(max_PCN = max(PIRACopyNumber))

max.and.min.plasmid.lengths <- max.plasmid.lengths %>%
    inner_join(min.plasmid.lengths) %>%
    inner_join(max.PCNs)

## 11388 plasmids
PIRA.PCN.estimates %>% nrow()

## these 11,388 plasmids are found in 4,317 genomes.
PIRA.PCN.estimates %>%
    count(AnnotationAccession) %>%
    nrow()

##4317 genomes containing plasmids
nrow(max.and.min.plasmid.lengths)

## only 1688 of these have plasmids found by themselves: 1688/4317 = 39.1%
max.and.min.plasmid.lengths %>%
    filter(max_replicon_length == min_replicon_length) %>%
    nrow()


max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids <- max.and.min.plasmid.lengths %>%
    filter(max_PCN > 10)

##1314 genomes containing multicopy plasmids here
nrow(max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids)

## only 159 are found by themselves: 159/1314 = 12%
max.and.min.plasmid.lengths.filtered.for.multicopy.plasmids %>%
    filter(max_replicon_length == min_replicon_length) %>%
    nrow()

## This is obviously statistically significant.
binom.test(x=159,n=1314, p = (1688/4317))

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
## Supplementary Figure S4.
## Plasmid length and copy number are conserved within plasmid taxonomic groups.

## IMPORTANT TODO: there is a bug, in which there are plasmids with NA PredictedMobility
## in the Acman and Redondo-Salvo panels, but not in any of the other panels.
## This doesn't change any message or anything, but needs to be sorted out.

## Supplementary Figure S4ABC
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
    make_PTU_length_rank_plot()

## plot copy numbers for all the PTUs
Acman.clique.PCN.plot <- Acman.cliques.with.PIRA.PCN.estimates %>%
    make_PTU_PCN_rank_plot()


## plot PCN against mean PTU length.
Acman.clique.PCN.vs.mean.size.plot <- Acman.cliques.with.PIRA.PCN.estimates %>%
    ## calculate mean PTU length and add as a column.
    add_mean_PTU_length_column("Clique") %>%
    make_PTU_mean_length_PCN_base_plot()

S4FigABC_title <- ggdraw() + draw_label("PTUs defined by Acman et al. (2020)", size=10, fontface="bold")

S4FigABC <- plot_grid(
    S4FigABC_title,
    plot_grid(Acman.clique.size.plot,
              Acman.clique.PCN.plot,
              Acman.clique.PCN.vs.mean.size.plot,
              labels=c('A','B','C'), nrow=1),
    nrow=2,rel_heights=c(0.1,1))


## Supplementary Figure S4DEF.
## The same result holds for the PTUs in the Redondo-Salvo et al. (2020) paper.

## That is, PTUs  in the plasmid similarity network
## have similar sizes and copy numbers.

## Importantly, PTU structure is in part determined by plasmid size, since large plasmids
## and small plasmids will not have a lot of shared sequence, by definition.

## 453 plasmids have copy numbers and assigned PTUs.
PTU.full.PIRA.PCN.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.PTU.data) %>%
    rank.correlate.column("PTU")

## Redondo-Salvo cliques show a very limited size distribution.
Redondo.Salvo.PTU.size.plot <- PTU.full.PIRA.PCN.estimates %>%
    make_PTU_length_rank_plot()

Redondo.Salvo.PTU.PCN.plot <- PTU.full.PIRA.PCN.estimates %>%
    make_PTU_PCN_rank_plot()

## plot PCN against mean PTU length.
Redondo.Salvo.PCN.vs.mean.size.plot <- PTU.full.PIRA.PCN.estimates %>%
    ## calculate mean PTU length and add as a column.
    add_mean_PTU_length_column("PTU") %>%
    make_PTU_mean_length_PCN_base_plot()

S4FigDEF_title <- ggdraw() + draw_label("PTUs defined by Redondo-Salvo et al. (2020)", size=10, fontface="bold")

S4FigDEF <- plot_grid(
    S4FigDEF_title,
    plot_grid(Redondo.Salvo.PTU.size.plot,
              Redondo.Salvo.PTU.PCN.plot,
              Redondo.Salvo.PCN.vs.mean.size.plot,
              labels=c('D','E','F'), nrow=1),
    nrow=2,
    rel_heights=c(0.1,1))


#########################################################
## Supplementary Figure S4GHI. MOB-Cluster PTU analysis

## First plot over primary cluster type.
## from $ mob_cluster -h:
## the mash distance for assigning the primary cluster ID is 0.06 by default.

MOB.typed.PIRA.clusters <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("primary_cluster_id")
    
## Make the plots.
MOB.Cluster.PTU.size.plot <- MOB.typed.PIRA.clusters %>%
    make_PTU_length_rank_plot()

MOB.Cluster.PTU.PCN.plot <- MOB.typed.PIRA.clusters %>%
    make_PTU_PCN_rank_plot()

## plot PCN against mean PTU length.
MOB.Cluster.PCN.vs.mean.size.plot <- MOB.typed.PIRA.clusters %>%
    ## calculate mean PTU length and add as a column.
    add_mean_PTU_length_column("primary_cluster_id") %>%
    make_PTU_mean_length_PCN_base_plot()

S4FigGHI_title <- ggdraw() + draw_label("PTUs defined by MOB-Cluster (Mash distance < 0.06)", size=10, fontface="bold")

S4FigGHI <- plot_grid(
    S4FigGHI_title,
    plot_grid(MOB.Cluster.PTU.size.plot,
              MOB.Cluster.PTU.PCN.plot,
              MOB.Cluster.PCN.vs.mean.size.plot,
              labels=c('G','H', 'I'), nrow=1),
    nrow=2,rel_heights=c(0.1,1))


#########################################################
## Supplementary Figure S4JKL. MOB-Typer Rep type analysis.
## The same thing holds for rep protein typing.
MOB.typed.PIRA.reptypes <- MOB.typed.PIRA.PCN.estimates %>%
    rank.correlate.column("rep_type.s.")

## MOB-typer Rep protein type classes by length
MOB.Typer.reptype.size.plot <- MOB.typed.PIRA.reptypes %>%
        make_PTU_length_rank_plot()

## MOB-typer Rep protein type classes by PCN
MOB.Typer.reptype.PCN.plot <- MOB.typed.PIRA.reptypes %>%
    make_PTU_PCN_rank_plot()

MOB.Typer.PCN.vs.mean.size.plot <- MOB.typed.PIRA.reptypes %>%
    ## calculate mean PTU length and add as a column.
    add_mean_PTU_length_column("rep_type.s.") %>%
    make_PTU_mean_length_PCN_base_plot()

S4FigJKL_title <- ggdraw() + draw_label("PTUs defined by MOB-Typer (Rep protein typing)", size=10, fontface="bold")

S4FigJKL <- plot_grid(
    S4FigJKL_title,
    plot_grid(
        MOB.Typer.reptype.size.plot,
        MOB.Typer.reptype.PCN.plot,
        MOB.Typer.PCN.vs.mean.size.plot,
        labels=c('J','K','L'), nrow=1),
    nrow=2, rel_heights=c(0.1,1))

############################
## Supplementary Figure S4MNO.
## the same thing holds over Rep protein types in Ares-Arroyo et al. (2023).
AresArroyo2023.typed.PIRA.reptypes <- PIRA.PCN.for.AresArroyo2023.data %>%
    rank.correlate.column("Replicase")

## Ares-Arroyo Rep protein type classes by length
AresArroyo2023.reptype.size.plot <- AresArroyo2023.typed.PIRA.reptypes %>%
    make_PTU_length_rank_plot()

## Ares-Arroyo Rep protein type classes by PCN
AresArroyo2023.reptype.PCN.plot <- AresArroyo2023.typed.PIRA.reptypes %>%
    make_PTU_PCN_rank_plot()

AresArroyo2023.PCN.vs.mean.size.plot <- AresArroyo2023.typed.PIRA.reptypes %>%
    ## calculate mean PTU length and add as a column.
    add_mean_PTU_length_column("Replicase") %>%
    make_PTU_mean_length_PCN_base_plot()

S4FigMNO_title <- ggdraw() + draw_label("PTUs defined by Rep types in Ares-Arroyo et al. (2023)", size=10, fontface="bold")

S4FigMNO <- plot_grid(
    S4FigMNO_title,
    plot_grid(
        AresArroyo2023.reptype.size.plot,
        AresArroyo2023.reptype.PCN.plot,
        AresArroyo2023.PCN.vs.mean.size.plot,
        labels=c('M','N','O'), nrow=1),
    nrow=2, rel_heights=c(0.1,1))

S4Fig_title <- ggdraw() + draw_label("Length is more conserved than copy numbers within plasmid taxonomy units (PTUs)", fontface='bold')

## save the full plot.
S4Fig <- plot_grid(S4Fig_title, S4FigABC, S4FigDEF, S4FigGHI, S4FigJKL, S4FigMNO, rel_heights = c(0.15,1,1,1,1,1), ncol=1)
ggsave("../results/S4Fig.pdf", S4Fig, width=8, height=13)


########################################
## Supplementary Figure S5. mobility group (relaxase) type analysis.

## plot over relaxase_type.
S5Fig <- MOB.typed.PIRA.PCN.estimates %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("relaxase_type.s.", "All other relaxase types") %>%
    make_PCN_base_plot() +
    facet_wrap(relaxase_type.s. ~ ., ncol=4) +
    ggtitle("MOB-Typer relaxase types") +
    guides(color = "none")
## save the plot
ggsave("../results/S5Fig.pdf", S5Fig)


##############################
## Supplementary Figures S6 Host range analysis.

## Supplementary Figure S6A
## plot over observed host range.
## This is "Taxon name of convergence of plasmids in MOB-suite plasmid DB",
## following the documentation here: https://github.com/phac-nml/mob-suite 
S6FigA <- MOB.typed.PIRA.PCN.estimates %>%
    ## put new lines in the host ranges to improve the aspect ratio of the subpanels.
    mutate(observed_host_range_ncbi_name = str_replace_all(observed_host_range_ncbi_name, ",", ",\n")) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("observed_host_range_ncbi_name", "All other host ranges") %>%
    make_PCN_base_plot() +
    facet_wrap(observed_host_range_ncbi_name ~ ., nrow=3) +
    ggtitle("Host range annotated by MOB-Typer") +
    guides(color = "none")

## Supplementary Figure S7B.
## 415 plasmids have copy numbers and host ranges.
RedondoSalvo.host.range.full.PIRA.estimates <- PIRA.PCN.estimates %>%
    inner_join(RedondoSalvo.PTU.data) %>%
    filter(Host_range != "-") %>%
    filter(!is.na(Host_range))

## This figure shows that copy number / plasmid size does not predict host range.
## there are narrow and broad host range plasmids both large and small.
## compare with the MOB-Typer result in this vein.
S6FigB <- RedondoSalvo.host.range.full.PIRA.estimates %>%
    ggplot(aes(
        x = log10(replicon_length),
        y = log10(PIRACopyNumber),
        color = Host_range)) +
    xlab("log10(length)") +
    ylab("log10(copy number)") +
    ## rename the legend
    labs(color = "host range") +
    geom_point(size=1,alpha=0.5) +
    theme_classic() +
    ggtitle("Host range annotated by\nRedondo-Salvo et al. (2020)") +
    theme(legend.position = "right") +
    ## make the points in the legend larger.
    guides(color = guide_legend(override.aes = list(size = 5)))

S6Fig <- plot_grid(
    S6FigA,
    plot_grid(S6FigB,"", rel_widths=c(1.5,1)),
    labels=c('A','B'), ncol=1,
    rel_heights=c(2.5,1))
## save the plot
ggsave("../results/S6Fig.pdf", S6Fig, height=8.5,width=7.1)


################################################################################
## Supplementary Figure S7.

## This vector is used for ordering Annotation levels in this figure.
order.by.total.plasmids <- make.plasmid.totals.col(PIRA.PCN.estimates)$Annotation

## Supplementary Figure S7A. let's make a histogram of PCN in these data.
S7FigA <- PIRA.PCN.estimates %>%
    ## order the Annotation levels.
    mutate(Annotation = factor(Annotation, levels = order.by.total.plasmids)) %>%
    ## TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    ggplot(aes(x = log10(PIRACopyNumber))) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "light gray") +
    ## place the high PCN at log10(50).
    geom_vline(xintercept = log10(50), linetype = "dashed", color = "light gray") +
    xlab("log10(copy number)") +
    ylab("count") +
    facet_wrap(. ~ Annotation) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    theme(
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x  = element_text(size=11),
        axis.text.y  = element_text(size=11))

## Supplementary Figure S7BC.
## Examine the tails of the PCN distribution.
## are low PCN (PCN < 1) and high PCN (PCN > 50) plasmids associated with any ecology?
## there is an enrichment of high PCN plasmids in human-impacted environments.

## How many total plasmids? 11,338.
PIRA.PCN.estimates %>%
    nrow()

## How many plasmids have PCN < 1? 2,376.
PIRA.PCN.estimates %>%
    filter(PIRACopyNumber < 1) %>%
    nrow()

## How many plasmids have PCN >= 1? 8,962
PIRA.PCN.estimates %>%
    filter(PIRACopyNumber >= 1) %>%
    nrow()


## calculate the fraction of low PCN plasmids in each category.
## Make Z-distributed confidence intervals for the fraction of isolates with
## PCN < 1.

low.PCN.plasmids.table <- make.lowPCN.table(PIRA.PCN.estimates) %>%
    ## TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA")

## plot the confidence intervals to see if there is any enrichment of low PCN plasmids in any ecological category.
S7FigB <- make.confint.figure.panel(
                                     low.PCN.plasmids.table,
                                     order.by.total.plasmids,
                                     "Proportion of plasmids with PCN < 1")

## calculate the fraction of high PCN plasmids in each category.
## there is an enrichment of  PCN > 50 plasmids in human-impacted environments.
## Make Z-distributed confidence intervals for the fraction of isolates with
## PCN > 50.

high.PCN.plasmids.table <- make.highPCN.table(PIRA.PCN.estimates) %>%
    ## TODO: fix upstream annotation so I don't have to do this filtering.
    filter(Annotation != "NA") %>%
    filter(Annotation != "blank")

## plot the confidence intervals to see if there is any enrichment of high PCN plasmids in any ecological category.
S7FigC <- make.confint.figure.panel(
    high.PCN.plasmids.table,
    order.by.total.plasmids,
    "Proportion of plasmids with PCN > 50")

## save the figure.
S7Fig <- plot_grid(S7FigA, S7FigB, S7FigC, labels=c('A','B','C'), ncol=1, rel_heights=c(2,1,1))
ggsave("../results/S7Fig.pdf", S7Fig, height = 8, width = 5)

## calculate the total number of plasmids,
## and the number of plasmids with PCN < 1 and PCN > 50.
PCN.count <- nrow(PIRA.PCN.estimates) ## 11,338 plasmids here
PCN.count

## There are 2376 plasmids with PCN < 1 in these data.
## this is 21% of plasmids
low.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber < 1) %>% nrow()
low.PCN.count
low.PCN.count/PCN.count

## There are 527 plasmids with PCN > 50 in these data.
## This is 4.6% of plasmids.
high.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber > 50) %>% nrow()
high.PCN.count
high.PCN.count/PCN.count

## There are 2434 plasmids with PCN > 10 in these data
## This is 21.5% of plasmids.
multicopy10.PCN.count <- PIRA.PCN.estimates %>% filter(PIRACopyNumber > 10) %>% nrow()
multicopy10.PCN.count
multicopy10.PCN.count/PCN.count

###################################################################################
## Figure 2. Coding Sequences (CDS) on plasmids follow an empirical scaling law.

## Main figure, all the points together.
## supplementary figure: same figure, separated by Annotation category.

## Fig 2A: show the combined plot
Fig2A <- CDS.rRNA.fraction.data %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB,
    ## but we have to make these annotations right before we make the figure, so
    ## that we don't accidentally filter out chromids when selecting plasmids.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    make_CDS_scaling_base_plot()

## Fig 2B: noncoding fraction analysis.
noncoding.fraction.data <- CDS.rRNA.fraction.data %>%
    mutate(noncoding_length = replicon_length - CDS_length) %>%
    mutate(noncoding_fraction = noncoding_length / replicon_length) %>%
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType))

## make a dataframe to add the trace line for Fig 2B.
Fig2B.mean.noncoding.fraction.per.length <- noncoding.fraction.data %>%
    filter(SeqType != "chromosome") %>%
    ## assign data into 100 bins by length
    mutate(bin = ntile(replicon_length, 100)) %>%
    group_by(bin) %>%
    summarize(
        noncoding_fraction = mean(noncoding_fraction),
        replicon_length = mean(replicon_length))


## now make Figure 2B.
Fig2B <- noncoding.fraction.data %>%
    ggplot(
        aes(
            x = log10(replicon_length),
            y = noncoding_fraction,
            color = SeqType)) +
    geom_point(size=0.5,alpha=0.8) +
    xlab("log10(length)") +
    ylab("noncoding fraction") +
    theme_classic() +
    guides(color = "none") +
    theme(strip.background = element_blank()) +
    theme(
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text.x  = element_text(size=11),
        axis.text.y  = element_text(size=11)) +
        geom_smooth(
        data = Fig2B.mean.noncoding.fraction.per.length,
        size = 0.8, alpha = 0.2, color = "dark gray", se=FALSE)


## Fig 2C: show generality over ecology.
Fig2C <- CDS.rRNA.fraction.data %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB,
    ## but we have to make these annotations right before we make the figure, so
    ## that we don't accidentally filter out chromids when selecting plasmids.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(.~Annotation)


## Now put together the complete Figure 2.
Fig2 <- plot_grid(
    plot_grid(Fig2A, Fig2B, rel_heights=c(3,2), nrow=2, labels=c('A','B')),
Fig2C, labels=c("",'C'),nrow=1, rel_widths=c(1,1.5))
## save the plot.
ggsave("../results/Fig2.pdf", Fig2, height=5.325, width=7.1)


## make plot for Yasa
Fig2X <- make_normalized_CDS_base_plot(CDS.rRNA.fraction.data)
## save the plot.
ggsave("../results/normalizedFig2X.pdf", Fig2X, height=4, width=7.1)


################################################################################
## Supplementary Figure S8. Break down the result in Figure 2 by genus
## to show universality of the CDS scaling relationship.

## Break down by taxonomic group.
S8FigA <- CDS.rRNA.fraction.data %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("TaxonomicGroup", "All other\ntaxonomic groups") %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB,
    ## but we have to make these annotations right before we make the figure, so
    ## that we don't accidentally filter out chromids when selecting plasmids.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicGroup, ncol=6) +
    ## improve the y-axis labels.
    scale_y_continuous(breaks=c(2, 7)) +
    ggtitle("Taxonomic groups")

## Break down by taxonomic subgroup
S8FigB <- CDS.rRNA.fraction.data %>%
    ## put new lines after slasks to improve the aspect ratio of the subpanels.
    mutate(TaxonomicSubgroup = str_replace_all(TaxonomicSubgroup, "/", "/\n")) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("TaxonomicSubgroup", "All other\ntaxonomic subgroups") %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB,
    ## but we have to make these annotations right before we make the figure, so
    ## that we don't accidentally filter out chromids when selecting plasmids.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup, ncol=6) +
    ## improve the y-axis labels.
    scale_y_continuous(breaks=c(2, 7)) +
    ggtitle("Taxonomic subgroups")

## Break down by genus.
S8FigC <- CDS.rRNA.fraction.data %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("Genus", "All other genera") %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB,
    ## but we have to make these annotations right before we make the figure, so
    ## that we don't accidentally filter out chromids when selecting plasmids.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    make_CDS_scaling_base_plot() +
    facet_wrap(. ~ Genus, ncol=10) +
    ## improve the y-axis labels.
    scale_y_continuous(breaks=c(2, 7)) +
    ggtitle("Microbial genera")


## save the plot over two pages.
S8FigAB_page <- plot_grid(S8FigA, S8FigB, labels=c("A","B"), ncol = 1, rel_heights=c(1.5,2))
S8FigC_page <- plot_grid(S8FigC, labels=c("C"))
ggsave("../results/S8FigAB.pdf", S8FigAB_page, width=9.5,height=11, limitsize=FALSE)
ggsave("../results/S8FigC.pdf", S8FigC_page, width=13,height=15, limitsize=FALSE)


########################################################################
## Figure 3: plasmid size scales with metabolic capacity.

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
    left_join(CDS.rRNA.fraction.data) %>%
    ## set NA values of metabolic_protein_count to zeros.
    mutate(metabolic_protein_count = ifelse(is.na(metabolic_protein_count), 0, metabolic_protein_count))

## join chromosome.length.data to metabolic.genes.in.chromosomes,
## since only 100 genomes were examined.
## KEY ASSUMPTION: each chromosome has at least one metabolic gene.
metabolic.gene.chromosome.data <- metabolic.genes.in.chromosomes %>%
    left_join(replicon.length.data) %>%
    ## get CDS data for each genome.
    left_join(CDS.rRNA.fraction.data)

## combine plasmid and chromosome metabolic gene data for Figure 3.
metabolic.gene.plasmid.and.chromosome.data <- metabolic.gene.plasmid.data %>%
    full_join(metabolic.gene.chromosome.data) %>%
    ## TODO: filter upstream Annotation so I don't need to do this ad hoc filtering here.
    filter(!is.na(Annotation)) %>%
    ## IMPORTANT: annotate chromids as plasmids that are longer than 500kB.
    mutate(SeqType = ifelse(
               SeqType == "plasmid" & replicon_length > PLASMID_LENGTH_THRESHOLD,
               "chromid", SeqType)) %>%
    ## The next two lines is for convenience when plotting the emergent metabolic scaling fit.
    mutate(log10_replicon_length = log10(replicon_length)) %>%
    mutate(log10_metabolic_protein_count = log10(metabolic_protein_count))

## fit a linear regression to the points with >= 100 metabolic genes, and save the fit,
## to show the emergent scaling trend on the figures.
metabolic.scaling.model <- lm(
    formula = log10_metabolic_protein_count ~ log10_replicon_length,
    data = filter(
        metabolic.gene.plasmid.and.chromosome.data,
        metabolic_protein_count >= METABOLIC_GENE_THRESHOLD))
## the slope estimate = 1.03, and the Adjusted R-squared = 0.8774.
summary(metabolic.scaling.model)

## look at the chromosome metabolic scaling regression model.
chromosome.metabolic.scaling.model <- lm(
    formula = log10_metabolic_protein_count ~ log10_replicon_length,
    data = filter(
        metabolic.gene.plasmid.and.chromosome.data,
        SeqType == "chromosome"))
## slope estimate = 0.987, Adjusted R-squared = 0.8903.
summary(chromosome.metabolic.scaling.model)


## make a table of the mean lengths per metabolic gene for plasmids & chromids for Figure 3.
Fig3.mean.length.per.metabolic.gene.table <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter(SeqType != "chromosome") %>%
    group_by(Annotation, log10_metabolic_protein_count) %>%
    summarize(log10_replicon_length = mean(log10_replicon_length))

## make a table of the mean number of metabolic genes for plasmids and chromosomes of a given length for Fig 3.
## IMPORTANT POINT: The smoothed curve omits all points with zero metabolic genes.
Fig3.mean.metabolic.genes.per.length <-  metabolic.gene.plasmid.and.chromosome.data %>%
    filter(SeqType != "chromosome") %>% 
    group_by(Annotation, log10_replicon_length) %>%
    summarize(metabolic_protein_count = mean(metabolic_protein_count)) %>%
    mutate(log10_metabolic_protein_count = log10(metabolic_protein_count))

## Fig3A: show the combined plot
Fig3A <- metabolic.gene.plasmid.and.chromosome.data %>%
    make_metabolic_scaling_base_plot() +
    geom_smooth(
        ##data = Fig3.mean.metabolic.genes.per.length,
        data = Fig3.mean.length.per.metabolic.gene.table,
        linewidth = 0.8, alpha = 0.2, color = "dark gray", se=FALSE)
   

## Fig3B: show generality over ecology.
Fig3B <- metabolic.gene.plasmid.and.chromosome.data %>%
    make_metabolic_scaling_base_plot() +
    geom_smooth(
        ##data = Fig3.mean.metabolic.genes.per.length,
        data = Fig3.mean.length.per.metabolic.gene.table,
        linewidth = 0.8, alpha = 0.2, color = "dark gray", se=FALSE) +
    facet_wrap(.~Annotation, nrow=3)

Fig3 <- plot_grid(Fig3A, Fig3B, labels = c('A', 'B'))
## save the plot.
ggsave("../results/Fig3.pdf", Fig3, height=4, width=7.1)

## make plot for Yasa
Fig3X <- make_normalized_metabolic_scaling_base_plot(metabolic.gene.plasmid.and.chromosome.data)
## save the plot.
ggsave("../results/normalizedFig3X.pdf", Fig3X, height=4, width=7.1)


################################################################################
## PLAYING AROUND. Make a similar figure to Figure 3, looking at numbers of ribosomal RNAs.
## (this is not on log-scale).
ribosomal.RNA.plot <- metabolic.gene.plasmid.and.chromosome.data %>%
    mutate(log10_total_rRNA_count = log10(total_rRNA_count)) %>%
    ggplot(
        aes(
            x = replicon_length,
            y = total_rRNA_count,
                color = SeqType)) +
        geom_point(size=0.5,alpha=0.8) +
        xlab("log10(length)") +
        ylab("log10(total ribosomal RNAs)") +
        theme_classic() +
##        guides(color = "none") +
        theme(strip.background = element_blank()) +
        theme(
            axis.title.x = element_text(size=11),
            axis.title.y = element_text(size=11),
            axis.text.x  = element_text(size=11),
            axis.text.y  = element_text(size=11))

ribosomal.RNA.plot

################################################################################
## Supplementary Figure S9. Break down the result in Figure 4 by genus
## to show universality of the CDS scaling relationship among genera containing chromids.
genera.containing.chromids <- filter(metabolic.gene.plasmid.and.chromosome.data, SeqType == "chromid")$Genus


## examining taxonomic groups 
S9FigA <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter(Genus %in% genera.containing.chromids) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("TaxonomicGroup", "All other\ntaxonomic groups") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicGroup, nrow = 1) +
    ggtitle("Taxonomic groups")

## examining taxonomic subgroups
S9FigB <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter(Genus %in% genera.containing.chromids) %>%
    ## put new lines after slasks to improve the aspect ratio of the subpanels.
    mutate(TaxonomicSubgroup = str_replace_all(TaxonomicSubgroup, "/", "/\n")) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column(
        "TaxonomicSubgroup", "All other\ntaxonomic subgroups") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ TaxonomicSubgroup, ncol=6) +
    ggtitle("Taxonomic subgroups")

S9FigAB <- plot_grid(S9FigA, S9FigB, labels=c("A","B"),ncol =1, rel_heights=c(1.25,2))

## Break down by genus, only showing genera containing megaplasmids.
S9FigC <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter(Genus %in% genera.containing.chromids) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("Genus", "All other genera") %>%
    make_metabolic_scaling_base_plot() +
    facet_wrap(. ~ Genus, ncol=5) +
    ggtitle("Microbial genera")

S9Fig <- plot_grid(S9FigAB, S9FigC, labels = c("","C"), ncol=1, rel_heights=c(2,3.5))

## save the plot.
ggsave("../results/S9Fig.pdf", S9Fig, height=14, width=8.5,limitsize=FALSE)

## ChatGPT says these are the default 3 colors in ggplot2.
## "#F8766D" (a red shade)
## "#00BA38" (a green shade)
## "#619CFF" (a blue shade)

## as a control, examine those genera that don't contain megaplasmids.
## we see the same trend, but it is less obvious.
S11Fig <- metabolic.gene.plasmid.and.chromosome.data %>%
    filter(!(Genus %in% genera.containing.chromids)) %>%
    filter.and.group.together.smaller.groups.in.the.correlate.column("Genus", "All other genera") %>%
    make_metabolic_scaling_base_plot() +
    scale_color_manual(values=c("#00BA38","#619CFF")) +
    facet_wrap(. ~ Genus, ncol=6)
## save the plot.
ggsave("../results/S11Fig.pdf", S11Fig, height=10, width=8)

################################################################################
## Supplementary Figure S12. Plot plasmid DNA content normalized by chromosome DNA content.
## plasmid DNA content is almost always < 0.5 of chromosome DNA content in the genome,
## indicated by the vertical dashed line.

plasmid.DNA.content.data <- full.PIRA.estimates %>%
    filter(SeqType != "chromosome") %>%
    group_by(AnnotationAccession, Annotation) %>%
    summarize(plasmid.DNA.content = sum(RepliconDNAContent))

chromosome.DNA.content.data <- full.PIRA.estimates %>%
    filter(SeqType == "chromosome") %>%
    group_by(AnnotationAccession, Annotation) %>%
    summarize(chromosome.DNA.content = sum(RepliconDNAContent))

normalized.plasmid.DNA.content.data <- full_join(plasmid.DNA.content.data, chromosome.DNA.content.data) %>%
    group_by(AnnotationAccession, Annotation) %>%
    summarize(normalized.plasmid.DNA.content = plasmid.DNA.content / chromosome.DNA.content)

S12Fig <- normalized.plasmid.DNA.content.data %>%
    ggplot(aes(x = log10(normalized.plasmid.DNA.content))) +
    geom_histogram(bins=100) +
    facet_grid(Annotation~.) +
    theme_classic() +
    geom_vline(xintercept = log10(0.5), linetype = "dashed", color="light gray") +
    theme(strip.background = element_blank())
## save the plot.
ggsave("../results/S12Fig.pdf", S12Fig, height=10, width=8)
