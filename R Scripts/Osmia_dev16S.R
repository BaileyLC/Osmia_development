##### Project: Osmia developmental microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of 16S rRNA gene data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(Biostrings) # Version 2.68.1
  library(ggplot2) # Version 3.4.3
  library(phyloseq) # Version 1.44.0
  library(vegan) # Version 2.6-4
  library(microbiome) # Version 1.22.0
  library(knitr) # Version 1.45
  library(magrittr) # Version 2.0.3
  library(decontam) # Version 1.20.0
  library(nlme) # Version 3.1-163
  library(grDevices) # Version 4.3.1
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(RVAideMemoire) # Version 0.9-83-7
  library(dplyr) # Version 1.1.3
  library(DESeq2) # Version 1.40.2

# Import data
  seqtab.nochim <- readRDS("Osmia_dev_seqs16S.rds")
  taxa <- readRDS("Osmia_dev_taxa16S.rds")
  meta16S_dev <- read.csv("Osmia_dev_master - 16S_all.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples <- data.frame(meta16S_dev)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  nesting_tube <- samples$nesting_tube
  sample_or_control <- samples$sample_or_control
  sampleinfo <- data.frame(extractionID = extractionID,
                           sample_type = sample_type,
                           sampleID = sampleID,
                           nesting_tube = nesting_tube,
                           sample_or_control = sample_or_control)
  rownames(sampleinfo) <- samples.out

# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                  sample_data(sampleinfo), 
                  tax_table(taxa))
  ps1
  
# Summarize the phyloseq obj contents before processing
  summarize_phyloseq(ps1)

## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Format data into a ggplot-friendly df
  df <- as.data.frame(sample_data(ps1))
  df$LibrarySize <- sample_sums(ps1)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))

# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_or_control)) + 
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
                        contam.prev = contamdf.prev$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.prev)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")

# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps1)
  ps.noncontam

# Remove control samples used for identifying contaminants
  ps_sub <- subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub

# Remove DNA from mitochondria & chloroplast
  ps2 <- ps_sub %>%
    subset_taxa(
      Kingdom == "Bacteria" &
        Family  != "mitochondria" &
        Class   != "Chloroplast"
    )

# Remove DNA from Eukarya, Eukaryota & Streptophyta
  ps2 <- ps2 %>%
    subset_taxa(Kingdom != "Eukarya" &
                Kingdom != "Eukaryota" &
                  Family != "Streptophyta")

# Remove DNA from Archaea
  ps2 <- ps2 %>%
    subset_taxa(Kingdom != "Archaea")

# Remove DNA from cyanobacteria & chlorplasts
  ps2 <- ps2 %>%
    subset_taxa(Phylum != "Cyanobacteria/Chloroplast")

# What remains in the phyloseq object?
  ps2

# Remove samples without any reads  
  ps3 <- prune_samples(sample_sums(ps2) != 0, ps2)
  ps3

# How many reads are in each sample? 
  sample_sums(ps3)

# What is the mean number of reads in all samples?
  mean(sample_sums(ps3))

# Summarize the phyloseq obj contents before processing
  summarize_phyloseq(ps3)
  
# Add Seq to each taxa name
  taxa_names(ps3) <- paste0("Seq", seq(ntaxa(ps3)))

# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps3), TRUE), 
                           sorted = 1:ntaxa(ps3), 
                           type = "OTUs")

# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps3), TRUE), 
                                  sorted = 1:nsamples(ps3), 
                                  type = "Samples"))

# Plot number of reads per OTU & sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") + 
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
## Species richness ----  
  
# Estimate Shannon, Simpson & observed richness
  bactrich <- estimate_richness(ps3, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bactrich$sample_type <- sample_data(ps3)$sample_type
  bactrich$sampleID <- sample_data(ps3)$sampleID
  bactrich$nesting_tube <- sample_data(ps3)$nesting_tube
  
# Plot Shannon, Simpson & observed richness
  plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
    theme_bw()
  
# Remove samples with 0 species richness 
  bactrich[bactrich == 0] <- NA
  bactrich <- bactrich[complete.cases(bactrich), ]
  
# Examine the effects of sample_type on Shannon index
  mod1 <- lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  anova(mod1)
  
# Examine the effects of sample_type on Simpson index
  mod2 <- lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  anova(mod2)
  
# Examine the effects of sample_type on observed richness
  mod3 <- lme(Observed ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  anova(mod3)
  
# Order samples on x-axis
  bactrich$sample_type <- factor(bactrich$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "dead"))
  
# Boxplot of Shannon index
  Osmia_dev_Shannon_bact <- ggplot(bactrich, aes(x = sample_type, y = Shannon, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = c("#FDD835", "#E4511E", "#43A047", "#0288D1", "#616161")) +
                                labs(title = "A") +
                                xlab("Sample type") +
                                ylab("Shannon index")
  Osmia_dev_Shannon_bact
  
# Boxplot of Simpson index
    Osmia_dev_Simpson_bact <- ggplot(bactrich, aes(x = sample_type, y = Simpson, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = c("#FDD835", "#E4511E", "#43A047", "#0288D1", "#616161")) +
                                labs(title = "A") +
                                xlab("Sample type") +
                                ylab("Simpson index")
    Osmia_dev_Simpson_bact

# Boxplot of Observed richness
  Osmia_dev_Observed_bact <- ggplot(bactrich, aes(x = sample_type, y = Observed, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = c("#FDD835", "#E4511E", "#43A047", "#0288D1", "#616161")) +
                                labs(title = "A") +
                                xlab("Sample type") +
                                ylab("Observed richness")
  Osmia_dev_Observed_bact
  
## Rarefaction ----
  
# Produce rarefaction curves
  tab <- otu_table(ps3)
  class(tab) <- "matrix"
  tab <- t(tab)
  rare <- rarecurve(tab, label = FALSE, tidy = FALSE)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_bact <- rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia_dev_rare_bact <- ggplot(rare_tidy_bact, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            geom_vline(xintercept = 20) +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "A") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  Osmia_dev_rare_bact

# Set seed and rarefy  
  set.seed(1234)
  rareps_bact <- rarefy_even_depth(ps3, sample.size = 20)
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray <- phyloseq::distance(rareps_bact, method = "bray")
  
# Convert to data frame
  samplebact <- data.frame(sample_data(rareps_bact))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact_perm <- adonis2(bact_bray ~ sample_type, data = samplebact)
  bact_perm
  
# Follow up with pairwise comparisons - which sample types differ?
  bact_perm_BH <- pairwise.perm.manova(bact_bray, samplebact$sample_type, p.method = "BH")
  bact_perm_BH
  
## Test for homogeneity of multivariate dispersion ----

# Calculate the average distance of group members to the group centroid
  disp_bact <- betadisper(bact_bray, samplebact$sample_type)
  disp_bact
  
# Do any of the group dispersions differ?
  disp_bact_an <- anova(disp_bact)
  disp_bact_an
  
# Which group dispersions differ?
  disp_bact_ttest <- permutest(disp_bact, 
                               control = permControl(nperm = 999),
                               pairwise = TRUE)
  disp_bact_ttest
  
# Which group dispersions differ?
  disp_bact_tHSD <- TukeyHSD(disp_bact)
  disp_bact_tHSD
  
## Ordination ----
  
# Calculate the relative abundance of each otu  
  ps.prop_bact <- transform_sample_counts(rareps_bact, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- ordinate(ps.prop_bact, method = "PCoA", distance = "bray")  
  
# Plot ordination
  Osmia_dev_PCoA_bact <- plot_ordination(ps.prop_bact, ord.pcoa.bray, color = "sample_type") + 
                            theme_bw() +
                            theme(text = element_text(size = 16)) +
                            theme(legend.justification = "left", 
                                  legend.title = element_text(size = 16, colour = "black"), 
                                  legend.text = element_text(size = 14, colour = "black")) +
                            theme(legend.position = "none") +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            geom_point(size = 3) +
                            scale_color_manual(values = c("#616161", "#E4511E", "#FDD835", "#43A047", "#0288D1")) + 
                            labs(color = "Developmental Stage") +
                            ggtitle("A")
  Osmia_dev_PCoA_bact
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe_ext <- usecol(Okabe_Ito, n = 80)
  colors <- sample(okabe_ext) 
  
# Sort data by Family 
  y1 <- tax_glom(rareps_bact, taxrank = 'Family') # agglomerate taxa
  y2 <- transform_sample_counts(y1, function(x) x/sum(x))
  y3 <- psmelt(y2)
  y3$Family <- as.character(y3$Family)
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  y3$Family <- as.factor(y3$Family)
  head(y3)
  
# Order samples on x-axis
  y3$sample_type <- factor(y3$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "dead"))
  
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
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Family for each sample
  Osmia_dev_fam_relabund_bact <- ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type,
                                               scale = "free", 
                                               space = "free") +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    xlab("Treatment") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 14)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 14, colour = "black"), 
                                          legend.text = element_text(size = 10, colour = "black")) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    ggtitle("A")
  Osmia_dev_fam_relabund_bact
  
# Sort data by Genus
  y4 <- tax_glom(rareps_bact, taxrank = 'Genus') # agglomerate taxa
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
    theme(text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 14, colour = "black"), 
          legend.text = element_text(size = 10, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Genus for each sample
  Osmia_dev_gen_relabund_bact <- ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free") +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    xlab("Treatment") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 14)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 14, colour = "black"), 
                                          legend.text = element_text(size = 10, colour = "black")) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    ggtitle("A")
  Osmia_dev_gen_relabund_bact
  
## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# When run again, change ps2 to rareps!   
  
# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq_to_deseq2(rareps_bact, ~ sample_type)
  
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
  
# Initial vs final provisions
  
# Extract results from differential abundance table for initial vs final provisions
  init_final <- results(desq_dds, contrast = c("sample_type", "initial provision", "final provision"))
  
# Order differential abundances by their padj value
  init_final <- init_final[order(init_final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_final_p05 <- init_final[(init_final$padj < alpha & !is.na(init_final$padj)), ]
  
# Check to see if any padj is below alpha
  init_final_p05
  
# Initial provisions vs larvae
  
# Extract results from differential abundance table for initial provisions vs larvae
  init_larva <- results(desq_dds, contrast = c("sample_type", "initial provision", "larva"))
  
# Order differential abundances by their padj value
  init_larva <- init_larva[order(init_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_larva_p05 <- init_larva[(init_larva$padj < alpha & !is.na(init_larva$padj)), ]
  
# Check to see if any padj is below alpha
  init_larva_p05
  
# Initial provisions vs pre-wintering adults
  
# Extract results from differential abundance table for initial provisions vs pre-wintering adults
  init_pre <- results(desq_dds, contrast = c("sample_type", "initial provision", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  init_pre <- init_pre[order(init_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_pre_p05 <- init_pre[(init_pre$padj < alpha & !is.na(init_pre$padj)), ]
  
# Check to see if any padj is below alpha
  init_pre_p05
  
# Initial provisions vs dead adults

# Extract results from differential abundance table for initial provisions vs dead adults
  init_dead <- results(desq_dds, contrast = c("sample_type", "initial provision", "dead"))
  
# Order differential abundances by their padj value
  init_dead <- init_dead[order(init_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_dead_p05 <- init_dead[(init_dead$padj < alpha & !is.na(init_dead$padj)), ]
  
# Check to see if any padj is below alpha
  init_dead_p05
  
# Final provisions vs larvae
  
# Extract results from differential abundance table for final provisions vs larvae
  final_larva <- results(desq_dds, contrast = c("sample_type", "final provision", "larva"))
  
# Order differential abundances by their padj value
  final_larva <- final_larva[order(final_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_larva_p05 <- final_larva[(final_larva$padj < alpha & !is.na(final_larva$padj)), ]
  
# Check to see if any padj is below alpha
  final_larva_p05
  
# Final provisions vs pre-wintering adults
  
# Extract results from differential abundance table for final provisions vs pre-wintering adults
  final_pre <- results(desq_dds, contrast = c("sample_type", "final provision", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  final_pre <- final_pre[order(final_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_pre_p05 <- final_pre[(final_pre$padj < alpha & !is.na(final_pre$padj)), ]
  
# Check to see if any padj is below alpha
  final_pre_p05
  
# Final provisions vs dead adults
  
# Extract results from differential abundance table for final provisions vs dead adults
  final_dead <- results(desq_dds, contrast = c("sample_type", "final provision", "dead"))
  
# Order differential abundances by their padj value
  final_dead <- final_dead[order(final_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_dead_p05 <- final_dead[(final_dead$padj < alpha & !is.na(final_dead$padj)), ]
  
# Check to see if any padj is below alpha
  final_dead_p05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva_pre <- results(desq_dds, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva_pre <- larva_pre[order(larva_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_pre_p05 <- larva_pre[(larva_pre$padj < alpha & !is.na(larva_pre$padj)), ]
  
# Check to see if any padj is below alpha
  larva_pre_p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva_dead <- results(desq_dds, contrast = c("sample_type", "larva", "dead"))
  
# Order differential abundances by their padj value
  larva_dead <- larva_dead[order(larva_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_dead_p05 <- larva_dead[(larva_dead$padj < alpha & !is.na(larva_dead$padj)), ]
  
# Check to see if any padj is below alpha
  larva_dead_p05
  
