##### Project: Osmia Developmental Microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of 16S rRNA gene data

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
  seqtab.nochim <- readRDS("Osmia_dev_seqs16S.rds")
  taxa <- readRDS("Osmia_dev_taxa16S.rds")
  meta16S.dev <- read.csv("Osmia_dev_master - 16S_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples.out <- stringr::str_sort(samples.out, numeric = TRUE)
  samples <- data.frame(meta16S.dev)
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
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                  sample_data(sample.info), 
                  tax_table(taxa))
  ps1

# Display total number of reads, mean, and se in phyloseq obj before processing
  sum(sample_sums(ps1))
  mean(sample_sums(ps1))
  print(plotrix::std.error(sample_sums(ps1)))

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

# Determine which ASVs are contaminants based on frequency of DNA in negative controls
  contamdf.freq <- decontam::isContaminant(ps1, conc = DNA_conc, method = "frequency", threshold = 0.1)
  table(contamdf.freq$contaminant)
  head(which(contamdf.freq$contaminant))
  
# Determine which ASVs are contaminants based on frequency of DNA in negative controls with a higher threshold
  contamdf.freq05 <- decontam::isContaminant(ps1, conc = DNA_conc, method = "frequency", threshold = 0.5)
  table(contamdf.freq05$contaminant)
  head(which(contamdf.freq05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls
  sample_data(ps1)$is.neg <- sample_data(ps1)$sample_or_control == "control"
  contamdf.prev <- decontam::isContaminant(ps1, method = "prevalence", neg = "is.neg", threshold = 0.1)
  table(contamdf.prev$contaminant)
  head(which(contamdf.prev$contaminant))

# Determine which ASVs are contaminants based on prevalence (presence/absence) in negative controls with a higher threshold
  contamdf.prev05 <- decontam::isContaminant(ps1, method = "prevalence", neg = "is.neg", threshold = 0.5)
  table(contamdf.prev05$contaminant)
  head(which(contamdf.prev05$contaminant))
  
# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls
  contamdf.comb <- decontam::isContaminant(ps1, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.1)
  table(contamdf.comb$contaminant)
  head(which(contamdf.comb$contaminant))

# Determine which ASVs are contaminants based on prevalence (presence/absence) and frequency in negative controls with a higher threshold  
  contamdf.comb05 <- decontam::isContaminant(ps1, conc = DNA_conc, neg = "is.neg", method = "combined", threshold = 0.5)
  table(contamdf.comb05$contaminant)
  head(which(contamdf.comb05$contaminant))

# Make phyloseq object of negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "control", ps1)

# Transform abundances in negative controls
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))

# Make phyloseq object of samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "sample", ps1)

# Transform abundances in samples
  ps.pos.presence <- phyloseq::transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))

# Make df of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contamdf.comb = contamdf.comb$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contamdf.comb)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")

# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.comb05$contaminant, ps1)
  ps.noncontam

# Remove control samples
  ps.sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps.sub

# Remove DNA from mitochondria & chloroplast
  ps2 <- ps.sub %>%
    phyloseq::subset_taxa(Kingdom == "Bacteria" &
                            Family  != "mitochondria" &
                            Class   != "Chloroplast"
    )

# Remove DNA from Eukarya, Eukaryota & Streptophyta
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Kingdom != "Eukarya" &
                          Kingdom != "Eukaryota" &
                            Family != "Streptophyta")

# Remove DNA from Archaea
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Kingdom != "Archaea")

# Remove DNA from cyanobacteria & chloroplasts
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Phylum != "Cyanobacteria/Chloroplast")

# Remove DNA from Legionella
  ps2 <- ps2 %>%
    phyloseq::subset_taxa(Genus != "Legionella")
  
# Remove samples without any reads
  ps3 <- phyloseq::prune_samples(sample_sums(ps2) != 0, ps2)
  ps3

# Display total number of reads, mean, and se in phyloseq obj after processing
  sum(sample_sums(ps3))
  mean(sample_sums(ps3))
  print(plotrix::std.error(sample_sums(ps3)))
  
# Calculate the reads per sample
  reads.sample <- microbiome::readcount(ps3)
  head(reads.sample)
  
# Add reads per sample to meta data
  sample_data(ps3)$reads.sample <- reads.sample
  
# Save sample metadata
  meta <- sample_data(ps3)
  
# How many samples, mean, and se for each developmental stage?
  meta %>%
    group_by(sample_type) %>%
    summarise(N = n(),
              mean = mean(reads.sample),
              se = sd(reads.sample)/sqrt(N),
              max = max(reads.sample),
              min = min(reads.sample))
  
# Save taxonomy and ASV counts data
  write.csv(tax_table(ps3), "Osmia_dev_16Staxa.csv")
  write.csv(otu_table(ps3), "Osmia_dev_16Sotu.csv")
  
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
  
## Richness and alpha diversity ----
  
# All pollen and bee samples
  
# Estimate Shannon index, Simpson index & observed richness
  bact.rich <- phyloseq::estimate_richness(ps3, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bact.rich$sample_type <- sample_data(ps3)$sample_type
  bact.rich$nesting_tube <- sample_data(ps3)$nesting_tube
  
# Plot richness and diversity
  phyloseq::plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
              theme_bw()
  
# Remove samples with 0 reads 
  bact.rich[bact.rich == 0] <- NA
  bact.rich <- bact.rich[complete.cases(bact.rich), ]
  
# Examine the effects of sample_type on Shannon index
  mod1 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = bact.rich)
  stats::anova(mod1)
  
# Examine the effects of sample_type on Simpson index
  mod2 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = bact.rich)
  stats::anova(mod2)
  
# Examine the effects of sample_type on observed richness
  mod3 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = bact.rich)
  stats::anova(mod3)
  
# Set color scheme  
  dev.colors <- c("fresh pollen egg" = "#FDD835",
                  "aged pollen" = "#E4511E",
                  "larva" = "#43A047",
                  "pre-wintering adult" = "#0288D1",
                  "dead adult" = "#616161")  
  
# Order samples on x-axis
  bact.rich$sample_type <- factor(bact.rich$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Boxplot of Shannon index
  Osmia.dev.Shannon.bact <- ggplot(bact.rich, aes(x = sample_type, y = Shannon, color = sample_type)) + 
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
                                scale_color_manual(values = dev.colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "A") +
                                xlab("Sample Type") +
                                ylim(0, 5) +
                                ylab("Shannon Index")
  Osmia.dev.Shannon.bact
  
# Boxplot of Simpson index
  Osmia.dev.Simpson.bact <- ggplot(bact.rich, aes(x = sample_type, y = Simpson, color = sample_type)) + 
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
                                   scale_color_manual(values = dev.colors) +
                                   scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                   labs(title = "A") +
                                   xlab("Sample Type") +
                                   ylim(0, 1) +
                                   ylab("Simpson Index")
    Osmia.dev.Simpson.bact

# Boxplot of Observed richness
  Osmia.dev.Observed.bact <- ggplot(bact.rich, aes(x = sample_type, y = Observed, color = sample_type)) + 
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
                                scale_color_manual(values = dev.colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "A") +
                                xlab("Sample Type") +
                                ylim(0, 40) +
                                ylab("Observed Richness")
  Osmia.dev.Observed.bact
  
# Only bee samples
  
# Subset phyloseq object to only include samples from larvae, pre-wintering adults, emerged adults, and dead adults
  ps4 <- phyloseq::subset_samples(ps3, sample_type != "fresh pollen egg")
  ps4 <- phyloseq::subset_samples(ps4, sample_type != "aged pollen")
  ps4
  
# Estimate Shannon index, Simpson index & observed richness
  bee.rich <- phyloseq::estimate_richness(ps4, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  bee.rich$sample_type <- sample_data(ps4)$sample_type
  bee.rich$nesting_tube <- sample_data(ps4)$nesting_tube
  
# Plot richness and diversity
  phyloseq::plot_richness(ps4, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
      theme_bw()
  
# Remove samples with 0 species reads 
  bee.rich[bee.rich == 0] <- NA
  bee.rich <- bee.rich[complete.cases(bee.rich), ]
  
# Examine the effects of sample_type on Shannon index
  mod4 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = bee.rich)
  stats::anova(mod4)
  
# Examine the effects of sample_type on Simpson index
  mod5 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = bee.rich)
  stats::anova(mod5)
  
# Examine the effects of sample_type on observed richness
  mod6 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = bee.rich)
  stats::anova(mod6)
  
## Beta diversity with relative abundance data ----
  
# All pollen and bee samples
  
# Calculate the relative abundance of each otu  
  ps.prop.bact <- phyloseq::transform_sample_counts(ps3, function(otu) otu/sum(otu))
  
# Save relative abundance data
  write.csv(otu_table(ps.prop.bact), "Osmia_dev_16Sotu_relabund.csv")
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray <- phyloseq::distance(ps.prop.bact, method = "bray")
  
# Convert to data frame
  sample.bact <- data.frame(sample_data(ps3))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact.perm.relabund <- vegan::adonis2(bact.bray ~ sample_type, data = sample.bact)
  bact.perm.relabund
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.BH <- RVAideMemoire::pairwise.perm.manova(bact.bray, sample.bact$sample_type, p.method = "BH")
  bact.perm.BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.relabund <- permute::how(within = Within(type = "free"),
                            plots = Plots(type = "none"),
                            blocks = sample.bact$nesting_tube,
                            observed = FALSE,
                            complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact.perm.relabund.pseudo <- vegan::adonis2(bact.bray ~ sample_type, permutations = perm.relabund, data = sample.bact)
  bact.perm.relabund.pseudo
  
# Only bee samples
  
# Calculate the relative abundance of each otu  
  ps.prop.bact.bee <- phyloseq::transform_sample_counts(ps4, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.bee <- phyloseq::distance(ps.prop.bact.bee, method = "bray")
  
# Convert to data frame
  sample.bact.bee <- data.frame(sample_data(ps4))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact.perm.relabund.bee <- vegan::adonis2(bact.bray.bee ~ sample_type, data = sample.bact.bee)
  bact.perm.relabund.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.BH.bee <- RVAideMemoire::pairwise.perm.manova(bact.bray.bee, sample.bact.bee$sample_type, p.method = "BH")
  bact.perm.BH.bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.relabund.bee <- permute::how(within = Within(type = "free"),
                                    plots = Plots(type = "none"),
                                    blocks = sample.bact.bee$nesting_tube,
                                    observed = FALSE,
                                    complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact.perm.relabund.pseudo.bee <- vegan::adonis2(bact.bray.bee ~ sample_type, permutations = perm.relabund.bee, data = sample.bact.bee)
  bact.perm.relabund.pseudo.bee
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# All pollen and bee samples  
  
# Calculate the average distance of group members to the group centroid
  disp.bact <- vegan::betadisper(bact.bray, sample.bact$sample_type)
  disp.bact
  
# Do any of the group dispersions differ?
  disp.bact.an <- stats::anova(disp.bact)
  disp.bact.an
  
# Which group dispersions differ?
  disp.bact.tHSD <- stats::TukeyHSD(disp.bact)
  disp.bact.tHSD
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.bact.bee <- vegan::betadisper(bact.bray.bee, sample.bact.bee$sample_type)
  disp.bact.bee
  
# Do any of the group dispersions differ?
  disp.bact.an.bee <- stats::anova(disp.bact.bee)
  disp.bact.an.bee
  
# Which group dispersions differ?
  disp.bact.tHSD.bee <- stats::TukeyHSD(disp.bact.bee)
  disp.bact.tHSD.bee

## Ordination with relative abundance data ----
  
# All pollen and bee samples
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- phyloseq::ordinate(ps.prop.bact, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact)$sample_type <- factor(sample_data(ps.prop.bact)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.bact <- plot_ordination(ps.prop.bact, ord.pcoa.bray, color = "sample_type") + 
                            theme_bw() +
                            theme(legend.position = "none",
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
                            labs(title = "A",
                                 color = "Sample Type")
  Osmia.dev.PCoA.bact

# Only bee samples  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.bee <- phyloseq::ordinate(ps.prop.bact.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.bact.bee)$sample_type <- factor(sample_data(ps.prop.bact.bee)$sample_type, levels = c("larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.bact.bee <- plot_ordination(ps.prop.bact.bee, ord.pcoa.bray.bee, color = "sample_type") + 
                                  theme_bw() +
                                  theme(legend.position = "none",
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
                                  labs(title = "A",
                                       color = "Sample Type")
  Osmia.dev.PCoA.bact.bee
  
## Rarefaction ----
  
# All pollen and bee samples  
  
# Produce rarefaction curves
  tab <- otu_table(ps3)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a "tidy" df
  rare.tidy.bact <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia.dev.rare.bact <- ggplot(rare.tidy.bact, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "A") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  Osmia.dev.rare.bact

# Set seed and rarefy  
  set.seed(1234)
  rareps.bact <- phyloseq::rarefy_even_depth(ps3, sample.size = 16)

# Only bee samples
  
# Produce rarefaction curves
  tab.bee <- otu_table(ps4)
  class(tab.bee) <- "matrix"
  tab.bee <- t(tab.bee)
  
# Save rarefaction data as a "tidy" df
  rare.tidy.bact.bee <- vegan::rarecurve(tab.bee, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia.dev.rare.bact.bee <- ggplot(rare.tidy.bact.bee, aes(x = Sample, y = Species, group = Site)) +
                                  geom_line() +
                                  theme_bw() +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  labs(title = "A") + 
                                  xlab("Number of reads") +
                                  ylab("Number of species")
  Osmia.dev.rare.bact.bee
  
# Set seed and rarefy  
  set.seed(1234)
  rareps.bact.bee <- phyloseq::rarefy_even_depth(ps4, sample.size = 15)
  
## Beta diversity with rarefied data ----  
  
# All pollen and bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.rare <- phyloseq::distance(rareps.bact, method = "bray")
  
# Convert to data frame
  sample.bact.rare <- data.frame(sample_data(rareps.bact))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact.perm.rare <- vegan::adonis2(bact.bray.rare ~ sample_type, data = sample.bact.rare)
  bact.perm.rare
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.BH.rare <- RVAideMemoire::pairwise.perm.manova(bact.bray.rare, sample.bact.rare$sample_type, p.method = "BH")
  bact.perm.BH.rare
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.rare <- permute::how(within = Within(type = "free"),
                            plots = Plots(type = "none"),
                            blocks = sample.bact.rare$nesting_tube,
                            observed = FALSE,
                            complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact.perm.rare.pseudo <- vegan::adonis2(bact.bray.rare ~ sample_type, permutations = perm.rare, data = sample.bact.rare)
  bact.perm.rare.pseudo
  
# Only bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact.bray.rare.bee <- phyloseq::distance(rareps.bact.bee, method = "bray")
  
# Convert to data frame
  sample.bact.rare.bee <- data.frame(sample_data(rareps.bact.bee))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact.perm.rare.bee <- vegan::adonis2(bact.bray.rare.bee ~ sample_type, data = sample.bact.rare.bee)
  bact.perm.rare.bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact.perm.BH.rare.bee <- RVAideMemoire::pairwise.perm.manova(bact.bray.rare.bee, sample.bact.rare.bee$sample_type, p.method = "BH")
  bact.perm.BH.rare.bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm.rare.bee <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = sample.bact.rare.bee$nesting_tube,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact.perm.rare.pseudo.bee <- vegan::adonis2(bact.bray.rare.bee ~ sample_type, permutations = perm.rare.bee, data = sample.bact.rare.bee)
  bact.perm.rare.pseudo.bee
  
## Test for homogeneity of multivariate dispersion with rarefied data ----

# All pollen and bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.bact.rare <- vegan::betadisper(bact.bray.rare, sample.bact.rare$sample_type)
  disp.bact.rare
  
# Do any of the group dispersions differ?
  disp.bact.an.rare <- stats::anova(disp.bact.rare)
  disp.bact.an.rare
  
# Which group dispersions differ?
  disp.bact.tHSD.rare <- stats::TukeyHSD(disp.bact.rare)
  disp.bact.tHSD.rare
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp.bact.rare.bee <- vegan::betadisper(bact.bray.rare.bee, sample.bact.rare.bee$sample_type)
  disp.bact.rare.bee
  
# Do any of the group dispersions differ?
  disp.bact.an.rare.bee <- stats::anova(disp.bact.rare.bee)
  disp.bact.an.rare.bee
  
# Which group dispersions differ?
  disp.bact.tHSD.rare.bee <- stats::TukeyHSD(disp.bact.rare.bee)
  disp.bact.tHSD.rare.bee
  
## Ordination with rarefied data ----
  
# All pollen and bee samples
  
# Calculate the relative abundance of each otu
  ps.prop.rare <- phyloseq::transform_sample_counts(rareps.bact, function(otu) otu/sum(otu))
  ps.prop.rare
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.rare <- phyloseq::ordinate(ps.prop.rare, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.rare)$sample_type <- factor(sample_data(ps.prop.rare)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.bact.rare <- plot_ordination(ps.prop.rare, ord.pcoa.bray.rare, color = "sample_type") + 
                                  theme_bw() +
                                  theme(legend.position = "none",
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
                                  labs(color = "Developmental Stage") +
                                  ggtitle("A")
  Osmia.dev.PCoA.bact.rare
  
# Only bee samples
  
# Calculate the relative abundance of each otu
  ps.prop.rare.bee <- phyloseq::transform_sample_counts(rareps.bact.bee, function(otu) otu/sum(otu))
  ps.prop.rare.bee
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray.rare.bee <- phyloseq::ordinate(ps.prop.rare.bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop.rare.bee)$sample_type <- factor(sample_data(ps.prop.rare.bee)$sample_type, levels = c("larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia.dev.PCoA.bact.rare.bee <- plot_ordination(ps.prop.rare.bee, ord.pcoa.bray.rare.bee, color = "sample_type") + 
                                      theme_bw() +
                                      theme(legend.position = "none",
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
                                      labs(color = "Developmental Stage") +
                                      ggtitle("A")
  Osmia.dev.PCoA.bact.rare.bee
  
## Stacked community plots ----
  
# Generate colorblind friendly palette
  Okabe.Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe.ext <- unikn::usecol(Okabe.Ito, n = 111)
  colors <- sample(okabe.ext) 
  
# New labels for facet_wrap
  new.labs <- c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")
  names(new.labs) <- c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult")  
  
# Agglomerate taxa by Family
  y1 <- phyloseq::tax_glom(ps3, taxrank = 'Family')
  
# Transform counts to relative abundances
  y2 <- phyloseq::transform_sample_counts(y1, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y3 <- phyloseq::psmelt(y2)
  
# Ensure Family is a chr  
  y3$Family <- as.character(y3$Family)
  
# Group Family be less that 1% abundance and rename
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  
# Ensure Genus is a factor 
  y3$Family <- as.factor(y3$Family)
  
# Order samples on x-axis
  y3$sample_type <- factor(y3$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot Family by sample type
  ggplot(data = y3, aes(x = sample_type, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample Type") +
    scale_x_discrete(labels = new.labs) +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left",
          legend.title = element_text(size = 14, colour = "black"), 
          legend.text = element_text(size = 7, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Family for each sample
  Osmia.dev.fam.relabund.bact <- ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
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
                                    theme(text = element_text(size = 16),
                                          legend.justification = "left", 
                                          legend.title = element_text(size = 14, colour = "black"),
                                          legend.text = element_text(size = 8, colour = "black"),
                                          strip.text = element_text(size = 7)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    theme(panel.spacing.x = unit(0.1, "lines")) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    ggtitle("A")
  Osmia.dev.fam.relabund.bact
  
# Agglomerate taxa by Genus
  y4 <- phyloseq::tax_glom(ps3, taxrank = 'Genus')
  
# Transform counts to relative abundances
  y5 <- phyloseq::transform_sample_counts(y4, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  y6 <- phyloseq::psmelt(y5)
  
# Ensure Genus is a chr   
  y6$Genus <- as.character(y6$Genus)
  
# Group Family be less that 1% abundance and rename
  y6$Genus[y6$Abundance < 0.01] <- "Genera < 1% abund."
  
# Ensure Genus is a factor  
  y6$Genus <- as.factor(y6$Genus)

# Order samples on x-axis
  y6$sample_type <- factor(y6$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))

# Plot Genus by sample type 
  ggplot(data = y6, aes(x = sample_type, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample Type") +
    scale_x_discrete(labels = new.labs) +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 14, colour = "black"), 
          legend.text = element_text(size = 7, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Genus for each sample
  Osmia.dev.gen.relabund.bact <- ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free",
                                               labeller = labeller(sample_type = new.labs)) +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    scale_x_discrete(expand = c(0, 1.5)) +
                                    xlab("Sample") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 16)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 16, colour = "black"), 
                                          legend.text = element_text(size = 12, colour = "black"),
                                          strip.text = element_text(size = 10)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    theme(panel.spacing.x = unit(0.1, "lines")) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    labs(fill = "Genera") +
                                    ggtitle("A")
  Osmia.dev.gen.relabund.bact
  
# Save plot  
  ggsave("Osmia_dev_16Sgenera.png", plot = Osmia.dev.gen.relabund.bact, width = 30, height = 10, unit = "in")

# Top 15 Genera
  
# Agglomerate taxa by Genus
  y4 <- phyloseq::tax_glom(ps3, taxrank = 'Genus')
  
# Identify the top 15 genera
  top15.bact.gen <- microbiome::top_taxa(y4, n = 15)
  
# Remove taxa that are not in the top 15
  ps.top15.bact.gen <- phyloseq::prune_taxa(top15.bact.gen, y4)
  
# Remove samples with 0 reads from the top 15 genera
  ps.top15.bact.gen <- phyloseq::prune_samples(sample_sums(ps.top15.bact.gen) != 0, ps.top15.bact.gen)
  
# Transform counts to relative abundances
  ps.top15.bact.gen.trans <- phyloseq::transform_sample_counts(ps.top15.bact.gen, function(x) x/sum(x))
  
# Convert to a ggplot2-friendly df
  ps.top15.bact.gen.trans.melt <- phyloseq::psmelt(ps.top15.bact.gen.trans)
  
# Order samples on x-axis
  ps.top15.bact.gen.trans.melt$sample_type <- factor(ps.top15.bact.gen.trans.melt$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot top 15 genera for each sample
  Osmia.dev.15gen.relabund.bact <- ggplot(data = ps.top15.bact.gen.trans.melt, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                      geom_bar(stat = "identity", position = "fill") + 
                                      scale_fill_manual(values = colors) +
                                      facet_grid(~ sample_type, 
                                                 scale = "free", 
                                                 space = "free",
                                                 labeller = labeller(sample_type = new.labs)) +
                                      theme(legend.position = "right") +
                                      ylab("Relative abundance") + 
                                      ylim(0, 1.0) +
                                      scale_x_discrete(expand = c(0, 1.5)) +
                                      xlab("Sample") +
                                      theme_bw() + 
                                      theme(text = element_text(size = 16)) +
                                      theme(panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank()) + 
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 14, colour = "black"), 
                                            legend.text = element_text(size = 8, colour = "black"),
                                            strip.text = element_text(size = 7)) + 
                                      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                      theme(panel.spacing.x = unit(0.1, "lines")) +
                                      guides(fill = guide_legend(ncol = 1)) +
                                      labs(fill = "Genera") +
                                      ggtitle("A")
  Osmia.dev.15gen.relabund.bact

## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html  

# All pollen and bee samples
  
# Convert from a phyloseq to a deseq obj
  desq.obj <- phyloseq::phyloseq_to_deseq2(rareps.bact, ~ sample_type)
  
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
  init.final <- DESeq2::results(desq.dds, contrast = c("sample_type", "fresh pollen egg", "aged pollen"))
  
# Order differential abundances by their padj value
  init.final <- init.final[order(init.final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init.final.p05 <- init.final[(init.final$padj < alpha & !is.na(init.final$padj)), ]
  
# Check to see if any padj is below alpha
  init.final.p05
  
# Fresh pollen egg vs larvae
  
# Extract results from differential abundance table for fresh pollen egg vs larvae
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
  
# Aged pollen vs dead adults
  
# Extract results from differential abundance table for aged pollen vs dead adults
  final.dead <- DESeq2::results(desq.dds, contrast = c("sample_type", "aged pollen", "dead adult"))
  
# Order differential abundances by their padj value
  final.dead <- final.dead[order(final.dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final.dead.p05 <- final.dead[(final.dead$padj < alpha & !is.na(final.dead$padj)), ]
  
# Check to see if any padj is below alpha
  final.dead.p05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva.pre <- DESeq2::results(desq.dds, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva.pre <- larva.pre[order(larva.pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva.pre.p05 <- larva.pre[(larva.pre$padj < alpha & !is.na(larva.pre$padj)), ]
  
# Check to see if any padj is below alpha
  larva.pre.p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva.dead <- DESeq2::results(desq.dds, contrast = c("sample_type", "larva", "dead adult"))
  
# Order differential abundances by their padj value
  larva.dead <- larva.dead[order(larva.dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva.dead.p05 <- larva.dead[(larva.dead$padj < alpha & !is.na(larva.dead$padj)), ]
  
# Check to see if any padj is below alpha
  larva.dead.p05
  
# Only bee samples
  
# Convert from a phyloseq to a deseq obj
  desq.obj.bee <- phyloseq::phyloseq_to_deseq2(rareps.bact.bee, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm.mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq.obj.bee), 1, gm.mean)
  
# Estimate size factors
  desq.dds.bee <- DESeq2::estimateSizeFactors(desq.obj.bee, geoMeans = geoMeans)
  
# Fit a local regression
  desq.obj.bee <- DESeq2::DESeq(desq.obj.bee, fitType = "local")
  
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
  
# Extract results from differential abundance table for larvae vs dead adults
  pre.dead.bee <- DESeq2::results(desq.dds.bee, contrast = c("sample_type", "pre.wintering.adult", "dead adult"))
  
# Order differential abundances by their padj value
  pre.dead.bee <- pre.dead.bee[order(pre.dead.bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre.dead.bee.p05 <- pre.dead.bee[(pre.dead.bee$padj < alpha & !is.na(pre.dead.bee$padj)), ]
  
# Check to see if any padj is below alpha
  pre.dead.bee.p05
  
