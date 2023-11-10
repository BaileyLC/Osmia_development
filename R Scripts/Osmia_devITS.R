##### Project: Osmia developmental microbiome

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of ITS rRNA gene data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(phyloseq) # Version 1.44.0
  library(vegan) # Version 2.6-4
  library(magrittr) # Version 2.0.3
  library(decontam) # Version 1.20.0
  library(nlme) # Version 3.1-163
  library(grDevices) # Version 4.3.1
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(ShortRead) # Version 1.58.0
  library(dplyr) # Version 1.1.3
  library(DESeq2) # Version 1.40.2

# Import data
  seqtab.nochim <- readRDS("Osmia_dev_seqsITS.rds")
  taxa <- readRDS("Osmia_dev_taxaITS.rds")
  metaITS_dev <- read.csv("Osmia_dev_master - ITS_all.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples <- data.frame(metaITS_dev)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  nesting_tube <- samples$nesting_tube
  sample_or_control <- samples$sample_or_control
  sampleinfo <- data.frame(extractionID = extractionID, sample_type = sample_type, sampleID = sampleID, nesting_tube = nesting_tube, sample_or_control = sample_or_control)
  rownames(sampleinfo) <- samples.out

# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                  sample_data(sampleinfo), 
                  tax_table(taxa))
  ps1

## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Format data into a ggplot-friendly df
  df <- as.data.frame(sample_data(ps1)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps1)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))

# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_type)) + 
    geom_point()

# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps1)$is.neg <- sample_data(ps1)$sample_or_control == "control"
  contamdf.prev <- isContaminant(ps1, method = "prevalence", neg = "is.neg")

# How many contaminants are there?
  table(contamdf.prev$contaminant)

# Which ASVs are contaminants?
  head(which(contamdf.prev$contaminant))

# Determine which ASVs are contaminants based on prevalence (presence/absence) higher than 0.5 in negative controls
  contamdf.prev05 <- isContaminant(ps1, method = "prevalence", neg = "is.neg", threshold = 0.5)

# How many contaminants are there?
  table(contamdf.prev05$contaminant)

# Which ASVs are contaminants?
  head(which(contamdf.prev05$contaminant))

# Make phyloseq object of presence-absence in negative controls
  ps.neg <- prune_samples(sample_data(ps1)$sample_or_control == "control", ps1)

# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))

# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- prune_samples(sample_data(ps1)$sample_or_control == "sample", ps1)

# Calculate taxa abundance in samples from sample counts
  ps.pos.presence <- transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))

# Make data.frame of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contam.prev=contamdf.prev$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.prev)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")

# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps1)
  ps.noncontam

# Remove control samples used for identifying contaminants
  ps_sub <- subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub

# Remove samples without any reads  
  ps2 <- prune_samples(sample_sums(ps_sub) != 0, ps_sub)

# How many reads are in each sample? 
  sample_sums(ps2)
  
# What is the mean number of reads in all samples?
  mean(sample_sums(ps2))  
  
# Add Seq to each taxa name
  taxa_names(ps2) <- paste0("Seq", seq(ntaxa(ps2)))

# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps2), TRUE), 
                           sorted = 1:ntaxa(ps2),
                           type = "OTUs")

# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps2), TRUE), 
                                             sorted = 1:nsamples(ps2), 
                                             type = "Samples"))

# Plot number of reads per OTU & sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") +
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")

## Species richness ----

# Calculate species richness
  fungrich <- estimate_richness(ps2, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))

# Build df with metadata
  fungrich$sample_type <- sample_data(ps2)$sample_type
  fungrich$sampleID <- sample_data(ps2)$sampleID
  fungrich$nesting_tube <- sample_data(ps2)$nesting_tube

# Plot species richness  
  plot_richness(ps2, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
    theme_bw()

# Remove samples with 0 species richness 
  fungrich[fungrich == 0] <- NA
  fungrich <- fungrich[complete.cases(fungrich), ]

# Examine the effects of sample_type on Shannon richness
  mod4 <- lme(Shannon ~ sample_type, random = ~1|sampleID, data = fungrich)
  anova(mod4)

# Examine the effects of sample_type on Simpson richness
  mod5 <- lme(Simpson ~ sample_type, random = ~1|sampleID, data = fungrich)
  anova(mod5)

# Examine the effects of sample_type on observed richness
  mod6 <- lme(Observed ~ sample_type, random = ~1|sampleID, data = fungrich)
  anova(mod6)

# Order samples on x-axis
  fungrich$sample_type <- factor(fungrich$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "emerged", "dead"))

# Boxplot of Shannon richness
  ggplot(fungrich, aes(x = sample_type, y = Shannon, fill = sample_type)) + 
    geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
    theme_bw() +
    scale_fill_manual(name = "Developmental Stage",
                      values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
    labs(title = "") +
    xlab("Developmental Stage") +
    ylab("Shannon richness")

# Boxplot of Simpson richness
  ggplot(fungrich, aes(x = sample_type, y = Simpson, fill = sample_type)) + 
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
    labs(title = "Simpson richness") +
    xlab("Developmental Stage")

# Boxplot of Observed richness
  ggplot(fungrich, aes(x = sample_type, y = Observed, fill = sample_type)) + 
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
    labs(title = "Observed richness") +
    xlab("Developmental Stage")

## Rarefaction ----

# Produce rarefaction curves
  tab <- otu_table(ps2)
  class(tab) <- "matrix"
  tab <- t(tab)
  rare <- rarecurve(tab, step = 20, label = FALSE)

# Rarefy
  set.seed(1234)
  rareps <- rarefy_even_depth(ps2, sample.size = 20)

# Perform PERMANOVA to test effects of developmental stage on bacterial community composition
  fung_bray <- phyloseq::distance(rareps, method = "bray")
  samplefung <- data.frame(sample_data(rareps))
  adonis2(fung_bray ~ sample_type, data = samplefung)

## Ordination ----

# Calculate the relative abundance of each otu
  ps.prop <- transform_sample_counts(rareps, function(otu) otu/sum(otu))
  ps.prop

# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- ordinate(ps.prop, method = "PCoA", distance = "bray")

# Plot ordination
  plot_ordination(rareps, ord.pcoa.bray, color = "sample_type") + 
    theme_bw() +
    theme(text = element_text(size = 16)) +
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    geom_point(size = 3) +
    scale_color_manual(values = c("#616161", "#9575CD", "#E4511E", "#FDD835", "#43A047", "#0288D1")) +
    labs(title = "PCoA: Bacteria", color = "Developmental Stage")

## Stacked community plot ----

# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Stretch palette (define more intermediate color options)
  okabe_ext <- usecol(Okabe_Ito, n = 52)
  colors <- sample(okabe_ext)

# Remove patterns in tax_table   
  tax_table(rareps)[, colnames(tax_table(rareps))] <- gsub(tax_table(rareps)[, colnames(tax_table(rareps))], pattern = "[a-z]__", replacement = "")

# Sort data by Family
  y1 <- tax_glom(rareps, taxrank = 'Family') # agglomerate taxa
  y2 <- transform_sample_counts(y1, function(x) x/sum(x))
  y3 <- psmelt(y2)
  y3$Family <- as.character(y3$Family)
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  y3$Family <- as.factor(y3$Family)
  head(y3)

# Reorder x-axis 
  y3$sample_type <- factor(y3$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "emerged", "dead"))

# Plot Family by sample type
  ggplot(data = y3, aes(x = sample_type, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Fungi")

# Plot Family for each sample
  ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) +
    facet_grid(~ sample_type, scale = "free", space = "free") +
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Fungi")

# Sort data by Genus
  y4 <- tax_glom(rareps, taxrank = 'Genus') # agglomerate taxa
  y5 <- transform_sample_counts(y4, function(x) x/sum(x))
  y6 <- psmelt(y5)
  y6$Genus <- as.character(y6$Genus)
  y6$Genus[y6$Abundance < 0.01] <- "Genera < 1% abund."
  y6$Genus <- as.factor(y6$Genus)
  head(y6)

# Order samples on x-axis
  y6$sample_type <- factor(y6$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "dead"))

# Plot Genus by sample type
  ggplot(data = y6, aes(x = sample_type, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Fungi")

# Plot Genus for each sample
  ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) +
    facet_grid(~ sample_type, scale = "free", space = "free") +
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Treatment") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Fungi")

## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html

# Remove patterns in tax_table   
  tax_table(ps2)[, colnames(tax_table(ps2))] <- gsub(tax_table(ps2)[, colnames(tax_table(ps2))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq_to_deseq2(ps2, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }

# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq_obj), 1, gm_mean)
  
# Estimate size factors
  desq_dds <- estimateSizeFactors(desq_obj, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds <- DESeq(desq_dds, fitType = "local")

# Set significance factor  
  alpha <- 0.05
  
# Initial vs final provision
  
# Extract results from differential abundance table for initial vs final provision
  init_final <- results(desq_dds, contrast = c("sample_type", "initial provision", "final provision"))
  
# Order differential abundances by their padj value
  init_final <- init_final[order(init_final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_final_p05 <- init_final[(init_final$padj < alpha & !is.na(init_final$padj)), ]
  
# Check to see if any padj is below alpha
  init_final_p05
  
# Final provision vs larva
  
# Extract results from differential abundance table for final provision vs larva
  final_larva <- results(desq_dds, contrast = c("sample_type", "final provision", "larva"))
  
# Order differential abundances by their padj value
  final_larva <- final_larva[order(final_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_larva_p05 <- final_larva[(final_larva$padj < alpha & !is.na(final_larva$padj)), ]
  
# Check to see if any padj is below alpha
  final_larva_p05
  
# Larva vs pre-wintering adult
  
# Extract results from differential abundance table for larva vs pre-wintering adult
  larva_pre <- results(desq_dds, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva_pre <- larva_pre[order(larva_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_pre_p05 <- larva_pre[(larva_pre$padj < alpha & !is.na(larva_pre$padj)), ]
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  larva_pre_p05 <- cbind(as(larva_pre_p05, "data.frame"),
                         as(tax_table(ps2)[rownames(larva_pre_p05), ], "matrix"))
  
# Plot
  ggplot(larva_pre_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw()
  
# Pre-wintering adult vs emerged
  
# Extract results from differential abundance table for pre-wintering adult vs emerged
  pre_emerg <- results(desq_dds, contrast = c("sample_type", "pre.wintering.adult", "emerged"))
  
# Order differential abundances by their padj value
  pre_emerg<- pre_emerg[order(pre_emerg$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_emerg_p05 <- pre_emerg[(pre_emerg$padj < alpha & !is.na(pre_emerg$padj)), ]
  
# Check to see if any padj is below alpha
  pre_emerg_p05
  
# Emerged vs dead
  
# Extract results from differential abundance table for emerged vs dead
  emerg_dead <- results(desq_dds, contrast = c("sample_type", "emerged", "dead"))
  
# Order differential abundances by their padj value
  emerg_dead <- emerg_dead[order(emerg_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_dead_p05 <- emerg_dead[(emerg_dead$padj < alpha & !is.na(emerg_dead$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_dead_p05
  