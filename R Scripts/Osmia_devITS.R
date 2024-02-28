##### Project: Osmia Developmental Microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of ITS rRNA gene data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(stringr) # Version 1.5.1
  library(phyloseq) # Version 1.44.0
  library(plotrix) # Version 3.8-4
  library(decontam) # Version 1.20.0
  library(microbiome) # Version 1.22.0
  library(ggplot2) # Version 3.4.3
  library(vegan) # Version 2.6-4
  library(magrittr) # Version 2.0.3
  library(nlme) # Version 3.1-163
  library(RColorBrewer) # Version 1.1-3
  library(unikn) # Version 0.9.0
  library(RVAideMemoire) # Version 0.9-83-7
  library(DESeq2) # Version 1.40.2
  library(tidyverse) # Version 1.2.0

# Import data
  seqtab.nochim <- readRDS("Osmia_dev_seqsITS.rds")
  taxa <- readRDS("Osmia_dev_taxaITS.rds")
  metaITS.dev <- read.csv("Osmia_dev_master - ITS_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples.out <- stringr::str_sort(samples.out, numeric = TRUE)
  samples <- data.frame(metaITS.dev)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  nesting_tube <- samples$nesting_tube
  sample_or_control <- samples$sample_or_control
  DNA_conc <- samples$DNA_conc
  sample.info <- data.frame(extractionID = extractionID,
                            sample_type = sample_type,
                            sampleID = sampleID,
                            nesting_tube = nesting_tube,
                            sample_or_control = sample_or_control,
                            DNA_conc = DNA_conc)
  rownames(sample.info) <- samples.out

# Format your data to work with phyloseq
  ps5 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                  sample_data(sample.info), 
                  tax_table(taxa))
  ps5

# Display total number of reads, mean, and se in phyloseq obj before processing
  sum(sample_sums(ps5))
  mean(sample_sums(ps5))
  print(plotrix::std.error(sample_sums(ps5)))

## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Format data into a ggplot-friendly df
  df <- as.data.frame(sample_data(ps5)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps5)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))

# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_type)) + 
    geom_point()

# Determine which ASVs are contaminants based on frequency of DNA in negative controls
  contamdf.freq <- decontam::isContaminant(ps5, conc = DNA_conc, method = "frequency", threshold = 0.1)
  table(contamdf.freq$contaminant)
  head(which(contamdf.freq$contaminant))
  
# Determine which ASVs are contaminants based on frequency of DNA in negative controls with a higher threshold
  contamdf.freq05 <- decontam::isContaminant(ps5, conc = DNA_conc, method = "frequency", threshold = 0.5)
  table(contamdf.freq05$contaminant)
  head(which(contamdf.freq05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps5)$is.neg <- sample_data(ps5)$sample_or_control == "control"
  contamdf.prev <- decontam::isContaminant(ps5, method = "prevalence", neg = "is.neg", threshold = 0.1)
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))

# Determine which ASVs are contaminants based on prevalence (presence/absence) higher than 0.5 in negative controls
  contamdf.prev05 <- decontam::isContaminant(ps5, method = "prevalence", neg = "is.neg", threshold = 0.5)
  table(contamdf.prev05$contaminant)
  head(which(contamdf.prev05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls
  contamdf.comb <- decontam::isContaminant(ps5, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.1)
  table(contamdf.comb$contaminant)
  head(which(contamdf.comb$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls with a higher threshold  
  contamdf.comb05 <- decontam::isContaminant(ps5, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.5)
  table(contamdf.comb05$contaminant)
  head(which(contamdf.comb05$contaminant))

# Make phyloseq object of presence-absence in negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps5)$sample_or_control == "control", ps5)

# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))

# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps5)$sample_or_control == "sample", ps5)

# Calculate taxa abundance in samples from sample counts
  ps.pos.presence <- phyloseq::transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))

# Make data.frame of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contamdf.comb05 = contamdf.comb05$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contamdf.comb05)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")

# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.comb05$contaminant, ps5)
  ps.noncontam

# Remove control samples used for identifying contaminants
  ps.sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps.sub

# Remove DNA from Pseudogymnoascus, the causative agent of white-nose syndrome in bats
  ps6 <- ps.sub %>%
    phyloseq::subset_taxa(Genus != "g__Pseudogymnoascus")
  
# Remove samples without any reads  
  ps7 <- phyloseq::prune_samples(sample_sums(ps6) != 0, ps6)
  ps7
  
# Display total number of reads, mean, and se in phyloseq obj after processing
  sum(sample_sums(ps7))
  mean(sample_sums(ps7))
  print(plotrix::std.error(sample_sums(ps7)))

# Calculate the reads per sample
  reads.sample <- microbiome::readcount(ps7)
  head(reads.sample)
  
# Add reads per sample to meta data
  sample_data(ps7)$reads.sample <- reads.sample  
  
# Save sample metadata
  meta <- sample_data(ps7)
  
# How many samples for each developmental stage?  
  meta %>%
    group_by(sample_type) %>%
    summarise(N = n(),
              mean = mean(reads.sample),
              se = sd(reads.sample)/sqrt(N),
              max = max(reads.sample),
              min = min(reads.sample))
  
# Remove patterns in tax_table   
  tax_table(ps7)[, colnames(tax_table(ps7))] <- gsub(tax_table(ps7)[, colnames(tax_table(ps7))], pattern = "[a-z]__", replacement = "")
  
# Save taxonomic and ASV counts
  write.csv(tax_table(ps7), "Osmia_dev_ITStaxa.csv")
  write.csv(otu_table(ps7), "Osmia_dev_ITSotu.csv")
  
# Add Seq to each taxa name
  taxa_names(ps7) <- paste0("Seq", seq(ntaxa(ps7)))

# Create a df containing the number of reads per OTU
  readsumsdf <- data.frame(nreads = sort(taxa_sums(ps7), TRUE), 
                           sorted = 1:ntaxa(ps7),
                           type = "OTUs")

# Add a column containing the number of reads per sample
  readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps7), TRUE), 
                                             sorted = 1:nsamples(ps7), 
                                             type = "Samples"))

# Plot number of reads per OTU & sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") +
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")

## Richness and alpha diversity ----

# All pollen and bee samples  
  
# Estimate Shannon index, Simpson index & observed richness
  fung.rich <- phyloseq::estimate_richness(ps7, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))

# Build df with metadata
  fung.rich$sample_type <- sample_data(ps7)$sample_type
  fung.rich$nesting_tube <- sample_data(ps7)$nesting_tube

# Plot richness and diversity
  phyloseq::plot_richness(ps7, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
              theme_bw()

# Remove samples with 0 reads
  fung.rich[fung.rich == 0] <- NA
  fung.rich <- fung.rich[complete.cases(fung.rich), ]

# Examine the effects of sample_type on Shannon index
  mod7 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = fung.rich)
  anova(mod7)

# Examine the effects of sample_type on Simpson index
  mod8 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = fung.rich)
  anova(mod8)

# Examine the effects of sample_type on observed richness
  mod9 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = fung.rich)
  anova(mod9)

# Set color scheme  
  dev.colors <- c("fresh pollen egg" = "#FDD835",
                  "aged pollen" = "#E4511E",
                  "larva" = "#43A047",
                  "pre-wintering adult" = "#0288D1",
                  "dead adult" = "#616161")  
  
# Order samples on x-axis
  fung.rich$sample_type <- factor(fung.rich$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Boxplot of Shannon index
  Osmia.dev.Shannon.fungi <- ggplot(fung.rich, aes(x = sample_type, y = Shannon, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(axis.text.x = element_text(size = 12, colour = "black"),
                                      axis.text.y = element_text(size = 12, colour = "black"),
                                      axis.title.x = element_text(size = 16, colour = "black"),
                                      axis.title.y = element_text(size = 16, colour = "black")) +
                                scale_color_manual(values = dev.colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "B") +
                                xlab("Sample type") +
                                ylim(0, 5) +
                                ylab("Shannon index")
  Osmia.dev.Shannon.fungi

# Boxplot of Simpson index
  Osmia.dev.Simpson.fungi <- ggplot(fung.rich, aes(x = sample_type, y = Simpson, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(axis.text.x = element_text(size = 12, colour = "black"),
                                      axis.text.y = element_text(size = 12, colour = "black"),
                                      axis.title.x = element_text(size = 16, colour = "black"),
                                      axis.title.y = element_text(size = 16, colour = "black")) +
                                scale_color_manual(values = dev.colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "B") +
                                xlab("Sample Type") +
                                ylim(0, 1) +
                                ylab("Simpson Index")
  Osmia.dev.Simpson.fungi

# Boxplot of Observed richness
  Osmia.dev.Observed.fungi <- ggplot(fung.rich, aes(x = sample_type, y = Observed, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                theme(axis.text.x = element_text(size = 12, colour = "black"),
                                      axis.text.y = element_text(size = 12, colour = "black"),
                                      axis.title.x = element_text(size = 14, colour = "black"),
                                      axis.title.y = element_text(size = 14, colour = "black")) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = dev.colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "B") +
                                xlab("Sample Type") +
                                ylim(0, 40) +
                                ylab("Observed Richness")
  Osmia.dev.Observed.fungi

# Only bee samples
  
# Subset phyloseq object to only include samples from larvae, pre-wintering adults, emerged adults, and dead adults
  ps8 <- phyloseq::subset_samples(ps7, sample_type != "fresh pollen egg")
  ps8 <- phyloseq::subset_samples(ps8, sample_type != "aged pollen")
  ps8

# Estimate Shannon index, Simpson index & observed richness
  fung.rich.bee <- phyloseq::estimate_richness(ps8, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata  
  fung.rich.bee$sample_type <- sample_data(ps8)$sample_type
  fung.rich.bee$nesting_tube <- sample_data(ps8)$nesting_tube
  
# Plot richness and diversity
  phyloseq::plot_richness(ps8, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
    theme_bw()
  
# Remove samples with 0 reads 
  fung.rich.bee[fung.rich.bee == 0] <- NA
  fung.rich.bee <- fung.rich.bee[complete.cases(fung.rich.bee), ]
  
# Examine the effects of sample_type on Shannon index
  mod10 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = fung.rich.bee)
  anova(mod10)
  
# Examine the effects of sample_type on Simpson index
  mod11 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = fung.rich.bee)
  anova(mod11)
  
# Examine the effects of sample_type on observed richness
  mod12 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = fung.rich.bee)
  anova(mod12)
  
## Beta diversity with relative abundance data ----  
  
# All pollen and bee samples  
  
# Calculate the relative abundance of each otu
  ps.prop.fung <- phyloseq::transform_sample_counts(ps7, function(otu) otu/sum(otu))
  ps.prop.fung
  
# Save relative abundance data
  write.csv(otu_table(ps.prop.fung), "Osmia_dev_ITSotu_relabund.csv")
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray <- phyloseq::distance(ps.prop.fung, method = "bray")
  
# Convert to data frame
  sample.fung <- data.frame(sample_data(ps7))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm <- vegan::adonis2(fung.bray ~ sample_type, data = sample.fung)
  fung.perm
  
# Follow up with pairwise comparisons - which sample types differ?
  fungi.perm.BH <- RVAideMemoire::pairwise.perm.manova(fung.bray, sample.fung$sample_type, p.method = "BH")
  fungi.perm.BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.relabund <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = sample.fung$nesting_tube,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition, dealing with pseudoreplication
  fung.perm.pseudo <- vegan::adonis2(fung.bray ~ sample_type, permutations = perm.relabund, data = sample.fung)
  fung.perm.pseudo
  
# Only bee samples
  
# Calculate the relative abundance of each otu
  ps.prop.fung.bee <- phyloseq::transform_sample_counts(ps8, function(otu) otu/sum(otu))
  ps.prop.fung.bee
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.bee <- phyloseq::distance(ps.prop.fung.bee, method = "bray")
  
# Convert to data frame
  sample.fung.bee <- data.frame(sample_data(ps8))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm.bee <- vegan::adonis2(fung.bray.bee ~ sample_type, data = sample.fung.bee)
  fung.perm.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  fungi.perm.BH.bee <- RVAideMemoire::pairwise.perm.manova(fung.bray.bee, sample.fung.bee$sample_type, p.method = "BH")
  fungi.perm.BH.bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.relabund.bee <- how(within = Within(type = "free"),
                           plots = Plots(type = "none"),
                           blocks = sample.fung.bee$nesting_tube,
                           observed = FALSE,
                           complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition, dealing with pseudoreplication
  fung.perm.pseudo.bee <- vegan::adonis2(fung.bray.bee ~ sample_type, permutations = perm.relabund.bee, data = sample.fung.bee)
  fung.perm.pseudo.bee
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# All pollen and bee samples  
  
# Calculate the average distance of group members to the group centroid
  disp.fung <- vegan::betadisper(fung.bray, sample.fung$sample_type)
  disp.fung
  
# Do any of the group dispersions differ?
  disp.fung.an <- anova(disp.fung)
  disp.fung.an
  
# Which group dispersions differ?
  disp.fung.tHSD <- stats::TukeyHSD(disp.fung)
  disp.fung.tHSD
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.fung.bee <- vegan::betadisper(fung.bray.bee, sample.fung.bee$sample_type)
  disp.fung.bee
  
# Do any of the group dispersions differ?
  disp.fung.an.bee <- anova(disp.fung.bee)
  disp.fung.an.bee

# Which group dispersions differ?
  disp.fung.tHSD.bee <- stats::TukeyHSD(disp.fung.bee)
  disp.fung.tHSD.bee

## Ordination with relative abundance data ----
  
# All pollen and bee samples  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- phyloseq::ordinate(ps.prop.fung, method = "PCoA", distance = "bray")

# Order samples
  sample_data(ps.prop.fung)$sample_type <- factor(sample_data(ps.prop.fung)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.fungi <- plot_ordination(ps.prop.fung, ord.pcoa.bray, color = "sample_type") + 
                            theme_bw() +
                            theme(text = element_text(size = 16),
                                  legend.justification = "left") + 
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            theme(axis.text.x = element_text(size = 12, colour = "black"),
                                  axis.text.y = element_text(size = 12, colour = "black"),
                                  axis.title.x = element_text(size = 14, colour = "black"),
                                  axis.title.y = element_text(size = 14, colour = "black")) +
                            geom_point(size = 4) +
                            scale_color_manual(values = dev.colors,
                                               labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                            labs(title = "B",
                                 color = "Sample Type")
  Osmia.dev.PCoA.fungi
  
# Only bee samples
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bee <- phyloseq::ordinate(ps.prop.fung.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.fung.bee)$sample_type <- factor(sample_data(ps.prop.fung.bee)$sample_type, levels = c("larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.fungi.bee <- plot_ordination(ps.prop.fung.bee, ord.pcoa.bray.bee, color = "sample_type") + 
                                  theme_bw() +
                                  theme(legend.justification = "left", 
                                        text = element_text(size = 16)) + 
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  theme(axis.text.x = element_text(size = 12, colour = "black"),
                                        axis.text.y = element_text(size = 12, colour = "black"),
                                        axis.title.x = element_text(size = 14, colour = "black"),
                                        axis.title.y = element_text(size = 14, colour = "black")) +
                                  geom_point(size = 4) +
                                  scale_color_manual(values = dev.colors,
                                                     labels = c("larvae", "pre-wintering adults", "dead adults")) +
                                  labs(title = "B",
                                       color = "Sample Type")
  Osmia.dev.PCoA.fungi.bee
  
## Rarefaction ----

# All pollen and bee samples

# Produce rarefaction curves
  tab <- otu_table(ps7)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a "tidy" df
  rare.tidy.fungi <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia.dev.rare.fungi <- ggplot(rare.tidy.fungi, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "B") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  Osmia.dev.rare.fungi

# Set seed and rarefy
  set.seed(1234)
  fung.rareps <- phyloseq::rarefy_even_depth(ps7, sample.size = 20)

# Bee samples only
  
# Produce rarefaction curves
  tab.bee <- otu_table(ps8)
  class(tab.bee) <- "matrix"
  tab.bee <- t(tab.bee)
  
# Save rarefaction data as a "tidy" df
  rare.tidy.fungi.bee <- vegan::rarecurve(tab.bee, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia.dev.rare.fungi.bee <- ggplot(rare.tidy.fungi.bee, aes(x = Sample, y = Species, group = Site)) +
                                geom_line() +
                                theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                labs(title = "B") + 
                                xlab("Number of reads") +
                                ylab("Number of species")
  Osmia.dev.rare.fungi.bee
  
# Set seed and rarefy
  set.seed(1234)
  fung.rareps.bee <- phyloseq::rarefy_even_depth(ps8, sample.size = 20)

## Beta diversity with rarefied data ----  
  
# All pollen and bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.rare <- phyloseq::distance(fung.rareps, method = "bray")
  
# Convert to data frame
  sample.fung.rare <- data.frame(sample_data(fung.rareps))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm.rare <- vegan::adonis2(fung.bray.rare ~ sample_type, data = sample.fung.rare)
  fung.perm.rare
  
# Follow up with pairwise comparisons - which sample types differ?
  fungi.perm.BH.rare <- RVAideMemoire::pairwise.perm.manova(fung.bray.rare, sample.fung.rare$sample_type, p.method = "BH")
  fungi.perm.BH.rare  
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.rare <- how(within = Within(type = "free"),
                   plots = Plots(type = "none"),
                   blocks = sample.fung.rare$nesting_tube,
                   observed = FALSE,
                   complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm.rare.pseudo <- vegan::adonis2(fung.bray.rare ~ sample_type, permutations = perm.rare, data = sample.fung.rare)
  fung.perm.rare.pseudo
  
# Follow up with pairwise comparisons - which sample types differ?
  fungi.perm.BH.rare.pseudo <- RVAideMemoire::pairwise.perm.manova(fung.bray.rare, sample.fung.rare$sample_type, p.method = "BH")
  fungi.perm.BH.rare.pseudo
  
# Only bee samples
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung.bray.rare.bee <- phyloseq::distance(fung.rareps.bee, method = "bray")
  
# Convert to data frame
  sample.fung.rare.bee <- data.frame(sample_data(fung.rareps.bee))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm.rare.bee <- vegan::adonis2(fung.bray.rare.bee ~ sample_type, data = sample.fung.rare.bee)
  fung.perm.rare.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  #fungi.perm.BH.rare.bee <- RVAideMemoire::pairwise.perm.manova(fung.bray.rare.bee, sample.fung.rare.bee$sample_type, p.method = "BH")
  #fungi.perm.BH.rare.bee  
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.rare.bee <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = sample.fung.rare.bee$nesting_tube,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung.perm.rare.pseudo <- vegan::adonis2(fung.bray.rare.bee ~ sample_type, permutations = perm.rare.bee, data = sample.fung.rare.bee)
  fung.perm.rare.pseudo

## Test for homogeneity of multivariate dispersion with rarefied data ----
  
# All pollen and bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare <- vegan::betadisper(fung.bray.rare, sample.fung.rare$sample_type)
  disp.fung.rare
  
# Do any of the group dispersions differ?
  disp.fung.an.rare <- anova(disp.fung.rare)
  disp.fung.an.rare
  
# Which group dispersions differ?
  disp.fung.tHSD.rare <- stats::TukeyHSD(disp.fung.rare)
  disp.fung.tHSD.rare
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.fung.rare.bee <- vegan::betadisper(fung.bray.rare.bee, sample.fung.rare.bee$sample_type)
  disp.fung.rare.bee
  
# Do any of the group dispersions differ?
  disp.fung.an.rare.bee <- anova(disp.fung.rare.bee)
  disp.fung.an.rare.bee
  
# Which group dispersions differ?
  disp.fung.tHSD.rare.bee <- stats::TukeyHSD(disp.fung.rare.bee)
  disp.fung.tHSD.rare.bee

## Ordination with rarefied data ----

# All pollen and bee samples  
  
# Calculate the relative abundance of each otu
  ps.prop.rare <- phyloseq::transform_sample_counts(fung.rareps, function(otu) otu/sum(otu))
  ps.prop.rare

# PCoA using Bray-Curtis distance
  ord.pcoa.bray.rare <- phyloseq::ordinate(ps.prop.rare, method = "PCoA", distance = "bray")

# Order samples
  sample_data(ps.prop.rare)$sample_type <- factor(sample_data(ps.prop.rare)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.fungi.rare <- plot_ordination(ps.prop.rare, ord.pcoa.bray.rare, color = "sample_type") + 
                                  theme_bw() +
                                  theme(legend.justification = "left", 
                                        text = element_text(size = 16)) + 
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  theme(axis.text.x = element_text(size = 12, colour = "black"),
                                        axis.text.y = element_text(size = 12, colour = "black"),
                                        axis.title.x = element_text(size = 14, colour = "black"),
                                        axis.title.y = element_text(size = 14, colour = "black")) +
                                  geom_point(size = 4) +
                                  scale_color_manual(values = dev.colors,
                                                     labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                  labs(title = "B",
                                       color = "Sample Type")
  Osmia.dev.PCoA.fungi.rare
  
# Only bee samples
  
# Calculate the relative abundance of each otu
  ps.prop.rare.bee <- phyloseq::transform_sample_counts(fung.rareps.bee, function(otu) otu/sum(otu))
  ps.prop.rare.bee
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.rare.bee <- phyloseq::ordinate(ps.prop.rare.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.rare.bee)$sample_type <- factor(sample_data(ps.prop.rare.bee)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.fungi.rare.bee <- plot_ordination(ps.prop.rare.bee, ord.pcoa.bray.rare.bee, color = "sample_type") + 
                                      theme_bw() +
                                      theme(legend.justification = "left", 
                                            text = element_text(size = 16)) + 
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      theme(axis.text.x = element_text(size = 12, colour = "black"),
                                            axis.text.y = element_text(size = 12, colour = "black"),
                                            axis.title.x = element_text(size = 14, colour = "black"),
                                            axis.title.y = element_text(size = 14, colour = "black")) +
                                      geom_point(size = 4) +
                                      scale_color_manual(values = dev.colors,
                                                         labels = c('larvae', 'pre-wintering adults', 'dead adults')) +
                                      labs(title = "B",
                                           color = "Sample Type")
  Osmia.dev.PCoA.fungi.rare.bee

## Stacked community plots ----

# Generate colorblind friendly palette
  Okabe.Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Stretch palette (define more intermediate color options)
  #okabe.ext <- unikn::usecol(Okabe.Ito, n = 98)
  #colors <- sample(okabe.ext)

# Remove patterns in tax_table
  tax_table(ps7)[, colnames(tax_table(ps7))] <- gsub(tax_table(ps7)[, colnames(tax_table(ps7))], pattern = "[a-z]__", replacement = "")

# New labels for facet_wrap
  new.labs <- c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")
  names(new.labs) <- c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult")  
  
# Agglomerate taxa by Family
  y7 <- phyloseq::tax_glom(ps7, taxrank = 'Family') # agglomerate taxa
  
# Transform counts to relative abundances
  y8 <- phyloseq::transform_sample_counts(y7, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df  
  y9 <- phyloseq::psmelt(y8)
  
# Ensure Family is a chr 
  y9$Family <- as.character(y9$Family)
  
# Group Families with less that 1% abundance and rename 
  y9$Family[y9$Abundance < 0.01] <- "Family < 1% abund."
  
# Ensure Family is a factor 
  y9$Family <- as.factor(y9$Family)

# Order samples on x-axis
  y9$sample_type <- factor(y9$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot Family by sample type
  ggplot(data = y9, aes(x = sample_type, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample Type") +
    scale_x_discrete(labels = c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")) +
    theme_bw() + 
    theme(text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 14, colour = "black"), 
          legend.text = element_text(size = 12, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "Fungi")
  
# Plot Family for each sample
  Osmia.dev.fam.relabund.fung <- ggplot(data = y9, aes(x = sampleID, y = Abundance, fill = Family)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type, 
                                               scale = "free",
                                               space = "free",
                                               labeller = labeller(sample_type = new.labs)) +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    xlab("Sample") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 12, colour = "black"),
                                          strip.text = element_text(size = 14)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    theme(panel.spacing.x=unit(0.1, "lines")) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    labs(title = "B")
  Osmia.dev.fam.relabund.fung

# Agglomerate taxa by Genus
  y10 <- phyloseq::tax_glom(ps7, taxrank = 'Genus')
  
# Transform counts to relative abundances  
  y11 <- phyloseq::transform_sample_counts(y10, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y12 <- phyloseq::psmelt(y11)
  
# Ensure Genus is a chr   
  y12$Genus <- as.character(y12$Genus)
  
# Group Genera with less that 1% abundance and rename  
  y12$Genus[y12$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor 
  y12$Genus <- as.factor(y12$Genus)

# Order samples on x-axis
  y12$sample_type <- factor(y12$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "emerged adult", "dead adult"))
  
# Plot Genus by sample type
  ggplot(data = y12, aes(x = sample_type, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample Type") +
    scale_x_discrete(labels = c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "emerged adults", "dead adults")) +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 12, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    labs(title = "Fungi")

# New labels for facet wrap (because these include emerged adults)
  all.labs <- c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "emerged adults", "dead adults")
  names(all.labs) <- c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "emerged adult", "dead adult")
  
# Plot Genus for each sample
  Osmia.dev.gen.relabund.fung <- ggplot(data = y12, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    theme_bw() +
                                    theme(legend.position = "right") +
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free",
                                               labeller = as_labeller(all.labs, default = label_wrap_gen(multi_line = TRUE, width = 13))) +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    scale_x_discrete(expand = c(0, 1.5)) +
                                    xlab("Sample") +
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "top", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 12, colour = "black"),
                                          strip.text = element_text(size = 14)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    theme(panel.spacing.x = unit(0.1, "lines")) +
                                    guides(fill = guide_legend(ncol = 2)) +
                                    labs(fill = "Genera",
                                         title = "B")
  Osmia.dev.gen.relabund.fung

# Agglomerate taxa by Genus
  y10 <- phyloseq::tax_glom(ps7, taxrank = 'Genus')
  
# Identify the top 15 genera
  top15.fung.gen <- microbiome::top_taxa(y10, n = 15)
  
# Remove taxa that are not in the top 15  
  ps.top15.fung.gen <- phyloseq::prune_taxa(top15.fung.gen, y10)
  
# Remove samples with 0 reads from the top 15 genera
  ps.top15.fung.gen <- phyloseq::prune_samples(sample_sums(ps.top15.fung.gen) != 0, ps.top15.fung.gen)
  
# Transform counts to relative abundances
  ps.top15.fung.gen.trans <- phyloseq::transform_sample_counts(ps.top15.fung.gen, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df  
  ps.top15.fung.gen.trans.melt <- phyloseq::psmelt(ps.top15.fung.gen.trans)
  
# Order samples on x-axis
  ps.top15.fung.gen.trans.melt$sample_type <- factor(ps.top15.fung.gen.trans.melt$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot top 15 genera for each sample
  Osmia.dev.15gen.relabund.fung <- ggplot(data = ps.top15.fung.gen.trans.melt, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                      geom_bar(stat = "identity", position = "fill") + 
                                      scale_fill_manual(values = colors) +
                                      facet_grid(~ sample_type, 
                                                 scale = "free", 
                                                 space = "free",
                                                 labeller = as_labeller(new.labs, default = label_wrap_gen(multi_line = TRUE, width = 13))) +
                                      ylab("Relative abundance") + 
                                      ylim(0, 1.0) +
                                      scale_x_discrete(expand = c(0, 3)) +
                                      xlab("Sample") +
                                      theme_bw() + 
                                      theme(text = element_text(size = 16)) +
                                      theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank()) + 
                                      theme(legend.position = "right",
                                            legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 12, colour = "black"),
                                            strip.text = element_text(size = 14)) +
                                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                      theme(panel.spacing.x = unit(0.1, "lines")) +
                                      guides(fill = guide_legend(ncol = 1)) +
                                      labs(fill = "Genera",
                                           title = "B")
  Osmia.dev.15gen.relabund.fung
  
## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html  

# All pollen and bee samples  
  
# Remove patterns in tax_table   
  tax_table(fung.rareps)[, colnames(tax_table(fung.rareps))] <- gsub(tax_table(fung.rareps)[, colnames(tax_table(fung.rareps))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq.obj <- phyloseq::phyloseq_to_deseq2(fung.rareps, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm.mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq.obj), 1, gm.mean)
  
# Estimate size factors
  desq.dds <- DESeq2::estimateSizeFactors(desq.obj, geoMeans = geoMeans)
  
# Fit a local regression
  desq.dds <- DESeq2::DESeq(desq.dds, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# Fresh pollen + egg vs aged pollen
  
# Extract results from differential abundance table for fresh pollen + egg vs aged pollen
  init.final <- results(desq.dds, contrast = c("sample_type", "fresh pollen egg", "aged pollen"))
  
# Order differential abundances by their padj value
  init.final <- init.final[order(init.final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init.final.p05 <- init.final[(init.final$padj < alpha & !is.na(init.final$padj)), ]
  
# Check to see if any padj is below alpha
  init.final.p05
  
# Fresh pollen + egg vs larvae
  
# Extract results from differential abundance table for fresh pollen + egg vs larvae
  init.larva <- DESeq2::results(desq.dds, contrast = c("sample_type", "fresh pollen egg", "larva"))
  
# Order differential abundances by their padj value
  init.larva <- init.larva[order(init.larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init.larva.p05 <- init.larva[(init.larva$padj < alpha & !is.na(init.larva$padj)), ]
  
# Check to see if any padj is below alpha
  init.larva.p05
  
# Fresh pollen + egg vs pre-wintering adults
  
# Extract results from differential abundance table for fresh pollen + egg vs pre-wintering adults
  init.pre <- DESeq2::results(desq.dds, contrast = c("sample_type", "fresh pollen egg", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  init.pre <- init.pre[order(init.pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init.pre.p05 <- init.pre[(init.pre$padj < alpha & !is.na(init.pre$padj)), ]
  
# Check to see if any padj is below alpha
  init.pre.p05
  
# Fresh pollen + egg vs dead adults
  
# Extract results from differential abundance table for fresh pollen + egg vs dead adults
  init.dead <- DESeq2::results(desq.dds, contrast = c("sample_type", "fresh pollen egg", "dead adult"))
  
# Order differential abundances by their padj value
  init.dead <- init.dead[order(init.dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init.dead.p05 <- init.dead[(init.dead$padj < alpha & !is.na(init.dead$padj)), ]
  
# Check to see if any padj is below alpha
  init.dead.p05
  
# Aged pollen vs larvae
  
# Extract results from differential abundance table for aged pollen vs larvae
  final.larva <- DESeq2::results(desq.dds, contrast = c("sample_type", "aged pollen", "larva"))
  
# Order differential abundances by their padj value
  final.larva <- final.larva[order(final.larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final.larva.p05 <- final.larva[(final.larva$padj < alpha & !is.na(final.larva$padj)), ]
  
# Check to see if any padj is below alpha
  final.larva.p05
  
# Aged pollen vs pre-wintering adults
  
# Extract results from differential abundance table for aged pollen vs pre-wintering adults
  final.pre <- DESeq2::results(desq.dds, contrast = c("sample_type", "aged pollen", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  final.pre <- final.pre[order(final.pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final.pre.p05 <- final.pre[(final.pre$padj < alpha & !is.na(final.pre$padj)), ]
  
# Check to see if any padj is below alpha
  final.pre.p05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva.pre <- DESeq2::results(desq.dds, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva.pre <- larva.pre[order(larva.pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva.pre.p05 <- larva.pre[(larva.pre$padj < alpha & !is.na(larva.pre$padj)), ]
  
# Check to see if any padj is below alpha
  larva.pre.p05

# Only bee samples  
  
# Remove patterns in tax_table   
  tax_table(fung.rareps.bee)[, colnames(tax_table(fung.rareps.bee))] <- gsub(tax_table(fung.rareps.bee)[, colnames(tax_table(fung.rareps.bee))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq.obj.bee <- phyloseq::phyloseq_to_deseq2(fung.rareps.bee, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm.mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }

# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq.obj.bee), 1, gm.mean)
  
# Estimate size factors
  desq.dds.bee <- DESeq2::estimateSizeFactors(desq.obj.bee, geoMeans = geoMeans)
  
# Fit a local regression
  desq.dds.bee <- DESeq2::DESeq(desq.dds.bee, fitType = "local")

# Set significance factor  
  alpha <- 0.05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva.pre.bee <- DESeq2::results(desq.dds.bee, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva.pre.bee <- larva.pre.bee[order(larva.pre.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva.pre.bee.p05 <- larva.pre.bee[(larva.pre.bee$padj < alpha & !is.na(larva.pre.bee$padj)), ]
  
# Check to see if any padj is below alpha
  larva.pre.bee.p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva.dead.bee <- DESeq2::results(desq.dds.bee, contrast = c("sample_type", "larva", "dead adult"))
  
# Order differential abundances by their padj value
  larva.dead.bee <- larva.dead.bee[order(larva.dead.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva.dead.bee.p05 <- larva.dead.bee[(larva.dead.bee$padj < alpha & !is.na(larva.dead.bee$padj)), ]
  
# Check to see if any padj is below alpha
  larva.dead.bee.p05
  
# Pre-wintering vs dead adults
  
# Extract results from differential abundance table for pre-wintering vs dead adults
  pre.dead.bee <- DESeq2::results(desq.dds.bee, contrast = c("sample_type", "pre.wintering.adult", "dead adult"))
  
# Order differential abundances by their padj value
  pre.dead.bee <- pre.dead.bee[order(pre.dead.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre.dead.bee.p05 <- pre.dead.bee[(pre.dead.bee$padj < alpha & !is.na(pre.dead.bee$padj)), ]
  
# Check to see if any padj is below alpha
  pre.dead.bee.p05
  
