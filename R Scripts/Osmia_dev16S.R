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
  meta16S_dev <- read.csv("Osmia_dev_master - 16S_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples.out <- stringr::str_sort(samples.out, numeric = TRUE)
  samples <- data.frame(meta16S_dev)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  nesting_tube <- samples$nesting_tube
  sample_or_control <- samples$sample_or_control
  DNA_conc <- samples$DNA_conc
  sampleinfo <- data.frame(extractionID = extractionID,
                           sample_type = sample_type,
                           sampleID = sampleID,
                           nesting_tube = nesting_tube,
                           sample_or_control = sample_or_control,
                           DNA_conc = DNA_conc)
  rownames(sampleinfo) <- samples.out

# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
                  sample_data(sampleinfo), 
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
  ps_sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub

# Remove DNA from mitochondria & chloroplast
  ps2 <- ps_sub %>%
    phyloseq::subset_taxa(
      Kingdom == "Bacteria" &
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
  reads_sample <- microbiome::readcount(ps3)
  head(reads_sample)
  
# Add reads per sample to meta data
  sample_data(ps3)$reads_sample <- reads_sample
  
# Save sample metadata
  meta <- sample_data(ps3)
  
# How many samples, mean, and se for each developmental stage?
  meta %>%
    group_by(sample_type) %>%
    summarise(N = n(),
              mean = mean(reads_sample),
              se = sd(reads_sample)/sqrt(N))
  
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
  
## Alpha diversity ----
  
# All pollen and bee samples
  
# Estimate Shannon, Simpson & observed richness
  bactrich <- phyloseq::estimate_richness(ps3, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bactrich$sample_type <- sample_data(ps3)$sample_type
  bactrich$nesting_tube <- sample_data(ps3)$nesting_tube
  
# Plot alpha diversity
  phyloseq::plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
              theme_bw()
  
# Remove samples with 0 reads 
  bactrich[bactrich == 0] <- NA
  bactrich <- bactrich[complete.cases(bactrich), ]
  
# Examine the effects of sample_type on Shannon index
  mod1 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  stats::anova(mod1)
  
# Examine the effects of sample_type on Simpson index
  mod2 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  stats::anova(mod2)
  
# Examine the effects of sample_type on observed richness
  mod3 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = bactrich)
  stats::anova(mod3)
  
# Set color scheme  
  dev_colors <- c("fresh pollen egg" = "#FDD835",
                  "aged pollen" = "#E4511E",
                  "larva" = "#43A047",
                  "pre-wintering adult" = "#0288D1",
                  "dead adult" = "#616161")  
  
# Order samples on x-axis
  bactrich$sample_type <- factor(bactrich$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Boxplot of Shannon index
  Osmia_dev_Shannon_bact <- ggplot(bactrich, aes(x = sample_type, y = Shannon, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(values = dev_colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "A") +
                                xlab("Sample Type") +
                                ylab("Shannon Index")
  Osmia_dev_Shannon_bact
  
# Boxplot of Simpson index
  Osmia_dev_Simpson_bact <- ggplot(bactrich, aes(x = sample_type, y = Simpson, color = sample_type)) + 
                                   geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                   geom_jitter(size = 1, alpha = 0.9) +
                                   theme_bw() +
                                   theme(legend.position = "none") +
                                   theme(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank()) +
                                   scale_color_manual(values = dev_colors) +
                                   scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                   labs(title = "A") +
                                   xlab("Sample Type") +
                                   ylab("Simpson Index")
    Osmia_dev_Simpson_bact

# Boxplot of Observed richness
  Osmia_dev_Observed_bact <- ggplot(bactrich, aes(x = sample_type, y = Observed, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(values = dev_colors) +
                                scale_x_discrete(labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) +
                                labs(title = "A") +
                                xlab("Sample Type") +
                                ylab("Observed Richness")
  Osmia_dev_Observed_bact
  
# Only bee samples
  
# Subset phyloseq object to only include samples from larvae, pre-wintering adults, emerged adults, and dead adults
  ps4 <- phyloseq::subset_samples(ps3, sample_type != "fresh pollen egg")
  ps4 <- phyloseq::subset_samples(ps4, sample_type != "aged pollen")
  ps4
  
# Estimate Shannon, Simpson & observed richness
  bee_rich <- phyloseq::estimate_richness(ps4, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata
  bee_rich$sample_type <- sample_data(ps4)$sample_type
  bee_rich$nesting_tube <- sample_data(ps4)$nesting_tube
  
# Plot alpha diversity
  phyloseq::plot_richness(ps4, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
      theme_bw()
  
# Remove samples with 0 species reads 
  bee_rich[bee_rich == 0] <- NA
  bee_rich <- bee_rich[complete.cases(bee_rich), ]
  
# Examine the effects of sample_type on Shannon index
  mod4 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = bee_rich)
  stats::anova(mod4)
  
# Examine the effects of sample_type on Simpson index
  mod5 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = bee_rich)
  stats::anova(mod5)
  
# Examine the effects of sample_type on observed richness
  mod6 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = bee_rich)
  stats::anova(mod6)
  
## Beta diversity with relative abundance data ----
  
# All pollen and bee samples
  
# Calculate the relative abundance of each otu  
  ps.prop_bact <- phyloseq::transform_sample_counts(ps3, function(otu) otu/sum(otu))
  
# Save relative abundance data
  write.csv(otu_table(ps.prop_bact), "Osmia_dev_16Sotu_relabund.csv")
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray <- phyloseq::distance(ps.prop_bact, method = "bray")
  
# Convert to data frame
  samplebact <- data.frame(sample_data(ps3))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact_perm_relabund <- vegan::adonis2(bact_bray ~ sample_type, data = samplebact)
  bact_perm_relabund
  
# Follow up with pairwise comparisons - which sample types differ?
  bact_perm_BH <- RVAideMemoire::pairwise.perm.manova(bact_bray, samplebact$sample_type, p.method = "BH")
  bact_perm_BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_relabund <- permute::how(within = Within(type = "free"),
                            plots = Plots(type = "none"),
                            blocks = samplebact$nesting_tube,
                            observed = FALSE,
                            complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact_perm_relabund_pseudo <- vegan::adonis2(bact_bray ~ sample_type, permutations = perm_relabund, data = samplebact)
  bact_perm_relabund_pseudo
  
# Only bee samples
  
# Calculate the relative abundance of each otu  
  ps.prop_bact_bee <- phyloseq::transform_sample_counts(ps4, function(otu) otu/sum(otu))
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray_bee <- phyloseq::distance(ps.prop_bact_bee, method = "bray")
  
# Convert to data frame
  samplebact_bee <- data.frame(sample_data(ps4))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact_perm_relabund_bee <- vegan::adonis2(bact_bray_bee ~ sample_type, data = samplebact_bee)
  bact_perm_relabund_bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact_perm_BH_bee <- RVAideMemoire::pairwise.perm.manova(bact_bray_bee, samplebact_bee$sample_type, p.method = "BH")
  bact_perm_BH_bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_relabund_bee <- permute::how(within = Within(type = "free"),
                                    plots = Plots(type = "none"),
                                    blocks = samplebact_bee$nesting_tube,
                                    observed = FALSE,
                                    complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact_perm_relabund_pseudo_bee <- vegan::adonis2(bact_bray_bee ~ sample_type, permutations = perm_relabund_bee, data = samplebact_bee)
  bact_perm_relabund_pseudo_bee
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# All pollen and bee samples  
  
# Calculate the average distance of group members to the group centroid
  disp_bact <- vegan::betadisper(bact_bray, samplebact$sample_type)
  disp_bact
  
# Do any of the group dispersions differ?
  disp_bact_an <- stats::anova(disp_bact)
  disp_bact_an
  
# Which group dispersions differ?
  disp_bact_ttest <- vegan::permutest(disp_bact, 
                                      control = permControl(nperm = 999),
                                      pairwise = TRUE)
  disp_bact_ttest
  
# Which group dispersions differ?
  disp_bact_tHSD <- stats::TukeyHSD(disp_bact)
  disp_bact_tHSD
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_bact_bee <- vegan::betadisper(bact_bray_bee, samplebact_bee$sample_type)
  disp_bact_bee
  
# Do any of the group dispersions differ?
  disp_bact_an_bee <- stats::anova(disp_bact_bee)
  disp_bact_an_bee
  
# Which group dispersions differ?
  disp_bact_ttest_bee <- vegan::permutest(disp_bact_bee, 
                                          control = permControl(nperm = 999),
                                          pairwise = TRUE)
  disp_bact_ttest_bee
  
# Which group dispersions differ?
  disp_bact_tHSD_bee <- stats::TukeyHSD(disp_bact_bee)
  disp_bact_tHSD_bee

## Ordination with relative abundance data ----
  
# All pollen and bee samples
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- phyloseq::ordinate(ps.prop_bact, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop_bact)$sample_type <- factor(sample_data(ps.prop_bact)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
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
                            scale_color_manual(values = dev_colors,
                                               labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) + 
                            labs(color = "Developmental Stage") +
                            ggtitle("A")
  Osmia_dev_PCoA_bact

# Only bee samples  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_bee <- phyloseq::ordinate(ps.prop_bact_bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop_bact_bee)$sample_type <- factor(sample_data(ps.prop_bact_bee)$sample_type, levels = c("larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia_dev_PCoA_bact_bee <- plot_ordination(ps.prop_bact_bee, ord.pcoa.bray_bee, color = "sample_type") + 
                                  theme_bw() +
                                  theme(text = element_text(size = 16)) +
                                  theme(legend.justification = "left", 
                                  legend.title = element_text(size = 16, colour = "black"), 
                                  legend.text = element_text(size = 14, colour = "black")) +
                                  theme(legend.position = "none") +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  geom_point(size = 3) +
                                  scale_color_manual(values = dev_colors,
                                                     labels = c("larvae", "pre-wintering adults", "dead adults")) + 
                                  labs(color = "Developmental Stage") +
                                  ggtitle("A")
  Osmia_dev_PCoA_bact_bee
  
## Rarefaction ----
  
# All pollen and bee samples  
  
# Produce rarefaction curves
  tab <- otu_table(ps3)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_bact <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia_dev_rare_bact <- ggplot(rare_tidy_bact, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "A") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  Osmia_dev_rare_bact

# Set seed and rarefy  
  set.seed(1234)
  rareps_bact <- phyloseq::rarefy_even_depth(ps3, sample.size = 16)

# Only bee samples
  
# Produce rarefaction curves
  tab_bee <- otu_table(ps4)
  class(tab_bee) <- "matrix"
  tab_bee <- t(tab_bee)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_bact_bee <- vegan::rarecurve(tab_bee, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia_dev_rare_bact_bee <- ggplot(rare_tidy_bact_bee, aes(x = Sample, y = Species, group = Site)) +
                                  geom_line() +
                                  theme_bw() +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  labs(title = "A") + 
                                  xlab("Number of reads") +
                                  ylab("Number of species")
  Osmia_dev_rare_bact_bee
  
# Set seed and rarefy  
  set.seed(1234)
  rareps_bact_bee <- phyloseq::rarefy_even_depth(ps4, sample.size = 15)
  
## Beta diversity with rarefied data ----  
  
# All pollen and bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray_rare <- phyloseq::distance(rareps_bact, method = "bray")
  
# Convert to data frame
  samplebact_rare <- data.frame(sample_data(rareps_bact))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact_perm_rare <- vegan::adonis2(bact_bray_rare ~ sample_type, data = samplebact_rare)
  bact_perm_rare
  
# Follow up with pairwise comparisons - which sample types differ?
  bact_perm_BH_rare <- RVAideMemoire::pairwise.perm.manova(bact_bray_rare, samplebact_rare$sample_type, p.method = "BH")
  bact_perm_BH_rare
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_rare <- permute::how(within = Within(type = "free"),
                            plots = Plots(type = "none"),
                            blocks = samplebact_rare$nesting_tube,
                            observed = FALSE,
                            complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact_perm_rare_pseudo <- vegan::adonis2(bact_bray_rare ~ sample_type, permutations = perm_rare, data = samplebact_rare)
  bact_perm_rare_pseudo
  
# Only bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray_rare_bee <- phyloseq::distance(rareps_bact_bee, method = "bray")
  
# Convert to data frame
  samplebact_rare_bee <- data.frame(sample_data(rareps_bact_bee))
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition
  bact_perm_rare_bee <- vegan::adonis2(bact_bray_rare_bee ~ sample_type, data = samplebact_rare_bee)
  bact_perm_rare_bee
  
# Follow up with pairwise comparisons - which sample types differ?
  bact_perm_BH_rare_bee <- RVAideMemoire::pairwise.perm.manova(bact_bray_rare_bee, samplebact_rare_bee$sample_type, p.method = "BH")
  bact_perm_BH_rare_bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_rare_bee <- permute::how(within = Within(type = "free"),
                                plots = Plots(type = "none"),
                                blocks = samplebact_rare_bee$nesting_tube,
                                observed = FALSE,
                                complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on bacterial community composition, dealing with pseudoreplication
  bact_perm_rare_pseudo_bee <- vegan::adonis2(bact_bray_rare_bee ~ sample_type, permutations = perm_rare_bee, data = samplebact_rare_bee)
  bact_perm_rare_pseudo_bee
  
## Test for homogeneity of multivariate dispersion with rarefied data ----

# All pollen and bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_bact_rare <- vegan::betadisper(bact_bray_rare, samplebact_rare$sample_type)
  disp_bact_rare
  
# Do any of the group dispersions differ?
  disp_bact_an_rare <- stats::anova(disp_bact_rare)
  disp_bact_an_rare
  
# Which group dispersions differ?
  disp_bact_ttest_rare <- vegan::permutest(disp_bact_rare, 
                                           control = permControl(nperm = 999),
                                           pairwise = TRUE)
  disp_bact_ttest_rare
  
# Which group dispersions differ?
  disp_bact_tHSD_rare <- stats::TukeyHSD(disp_bact_rare)
  disp_bact_tHSD_rare
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_bact_rare_bee <- vegan::betadisper(bact_bray_rare_bee, samplebact_rare_bee$sample_type)
  disp_bact_rare_bee
  
# Do any of the group dispersions differ?
  disp_bact_an_rare_bee <- stats::anova(disp_bact_rare_bee)
  disp_bact_an_rare_bee
  
# Which group dispersions differ?
  disp_bact_ttest_rare_bee <- vegan::permutest(disp_bact_rare_bee, 
                                               control = permControl(nperm = 999),
                                               pairwise = TRUE)
  disp_bact_ttest_rare_bee
  
# Which group dispersions differ?
  disp_bact_tHSD_rare_bee <- stats::TukeyHSD(disp_bact_rare_bee)
  disp_bact_tHSD_rare_bee
  
## Ordination with rarefied data ----
  
# All pollen and bee samples
  
# Calculate the relative abundance of each otu
  ps.prop_rare <- phyloseq::transform_sample_counts(rareps_bact, function(otu) otu/sum(otu))
  ps.prop_rare
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_rare <- phyloseq::ordinate(ps.prop_rare, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop_rare)$sample_type <- factor(sample_data(ps.prop_rare)$sample_type, levels = c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia_dev_PCoA_bact_rare <- plot_ordination(ps.prop_rare, ord.pcoa.bray_rare, color = "sample_type") + 
                                  theme_bw() +
                                  theme(text = element_text(size = 16)) +
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 16, colour = "black"), 
                                        legend.text = element_text(size = 14, colour = "black")) +
                                  theme(legend.position = "none") +
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  geom_point(size = 3) +
                                  scale_color_manual(values = dev_colors,
                                                     labels = c('fresh pollen + egg', 'aged pollen', 'larvae', 'pre-wintering adults', 'dead adults')) + 
                                  labs(color = "Developmental Stage") +
                                  ggtitle("A")
  Osmia_dev_PCoA_bact_rare
  
# Only bee samples  
  
# Calculate the relative abundance of each otu
  ps.prop_rare_bee <- phyloseq::transform_sample_counts(rareps_bact_bee, function(otu) otu/sum(otu))
  ps.prop_rare_bee
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_rare_bee <- phyloseq::ordinate(ps.prop_rare_bee, method = "PCoA", distance = "bray")
  
# Order samples
  sample_data(ps.prop_rare_bee)$sample_type <- factor(sample_data(ps.prop_rare_bee)$sample_type, levels = c("larva", "pre-wintering adult", "dead adult"))
  
# Plot ordination
  Osmia_dev_PCoA_bact_rare_bee <- plot_ordination(ps.prop_rare_bee, ord.pcoa.bray_rare_bee, color = "sample_type") + 
                                      theme_bw() +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) +
                                      theme(legend.position = "none") +
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      geom_point(size = 3) +
                                      scale_color_manual(values = dev_colors,
                                                         labels = c("larvae", "pre-wintering adults", "dead adults")) + 
                                      labs(color = "Developmental Stage") +
                                      ggtitle("A")
  Osmia_dev_PCoA_bact_rare_bee
  
## Stacked community plots ----
  
# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe_ext <- unikn::usecol(Okabe_Ito, n = 80)
  colors <- sample(okabe_ext) 
  
# New labels for facet_wrap
  new_labs <- c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")
  names(new_labs) <- c("fresh pollen egg", "aged pollen", "larva", "pre-wintering adult", "dead adult")  
  
# Sort data by Family 
  y1 <- phyloseq::tax_glom(rareps_bact, taxrank = 'Family') # agglomerate taxa
  y2 <- phyloseq::transform_sample_counts(y1, function(x) x/sum(x))
  y3 <- phyloseq::psmelt(y2)
  y3$Family <- as.character(y3$Family)
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  y3$Family <- as.factor(y3$Family)
  head(y3)
  
# Save relative abundance data
  write.csv(y3, "Osmia_dev_Fam_bact_relabund.csv")
  
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
    scale_x_discrete(labels = c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")) +
    theme_bw() +
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left",
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 8, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Family for each sample
  Osmia_dev_fam_relabund_bact <- ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type,
                                               scale = "free", 
                                               space = "free",
                                               labeller = labeller(sample_type = new_labs)) +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    xlab("Sample") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 14)) +
                                    theme(panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 14, colour = "black"), 
                                          legend.text = element_text(size = 8, colour = "black"),
                                          strip.text = element_text(size = 8)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    ggtitle("A")
  Osmia_dev_fam_relabund_bact
  
# Sort data by Genus
  y4 <- phyloseq::tax_glom(rareps_bact, taxrank = 'Genus') # agglomerate taxa
  y5 <- phyloseq::transform_sample_counts(y4, function(x) x/sum(x))
  y6 <- phyloseq::psmelt(y5)
  y6$Genus <- as.character(y6$Genus)
  y6$Genus[y6$Abundance < 0.01] <- "Genera < 1% abund."
  y6$Genus <- as.factor(y6$Genus)
  head(y6)
  
# Save relative abundance data
  write.csv(y6, "Osmia_dev_Gen_bact_relabund.csv")

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
    scale_x_discrete(labels = c("fresh pollen + egg", "aged pollen", "larvae", "pre-wintering adults", "dead adults")) +
    theme_bw() + 
    theme(text = element_text(size = 14)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 14, colour = "black"), 
          legend.text = element_text(size = 8, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Bacteria")
  
# Plot Genus for each sample
  Osmia_dev_gen_relabund_bact <- ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) +
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free",
                                               labeller = labeller(sample_type = new_labs)) +
                                    theme(legend.position = "right") +
                                    ylab("Relative abundance") + 
                                    ylim(0, 1.0) +
                                    xlab("Sample") +
                                    theme_bw() + 
                                    theme(text = element_text(size = 14)) +
                                    theme(panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank()) + 
                                    theme(legend.justification = "left", 
                                          legend.title = element_text(size = 14, colour = "black"), 
                                          legend.text = element_text(size = 8, colour = "black"),
                                          strip.text = element_text(size = 8)) + 
                                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                                    guides(fill = guide_legend(ncol = 3)) +
                                    ggtitle("A")
  Osmia_dev_gen_relabund_bact
  
## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html  

# All pollen and bee samples  
  
# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq::phyloseq_to_deseq2(rareps_bact, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq_obj), 1, gm_mean)

# Estimate size factors
  desq_dds <- DESeq2::estimateSizeFactors(desq_obj, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds <- DESeq2::DESeq(desq_dds, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# Fresh pollen + egg vs aged pollen
  
# Extract results from differential abundance table for fresh pollen + egg vs aged pollen
  init_final <- DESeq2::results(desq_dds, contrast = c("sample_type", "fresh pollen egg", "aged pollen"))
  
# Order differential abundances by their padj value
  init_final <- init_final[order(init_final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_final_p05 <- init_final[(init_final$padj < alpha & !is.na(init_final$padj)), ]
  
# Check to see if any padj is below alpha
  init_final_p05
  
# Fresh pollen egg vs larvae
  
# Extract results from differential abundance table for fresh pollen egg vs larvae
  init_larva <- DESeq2::results(desq_dds, contrast = c("sample_type", "fresh pollen egg", "larva"))
  
# Order differential abundances by their padj value
  init_larva <- init_larva[order(init_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_larva_p05 <- init_larva[(init_larva$padj < alpha & !is.na(init_larva$padj)), ]
  
# Check to see if any padj is below alpha
  init_larva_p05
  
# Fresh pollen + egg vs pre-wintering adults
  
# Extract results from differential abundance table for fresh pollen + egg vs pre-wintering adults
  init_pre <- DESeq2::results(desq_dds, contrast = c("sample_type", "fresh pollen egg", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  init_pre <- init_pre[order(init_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_pre_p05 <- init_pre[(init_pre$padj < alpha & !is.na(init_pre$padj)), ]
  
# Check to see if any padj is below alpha
  init_pre_p05
  
# Fresh pollen + egg vs dead adults
  
# Extract results from differential abundance table for fresh pollen + egg vs dead adults
  init_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "fresh pollen egg", "dead adult"))
  
# Order differential abundances by their padj value
  init_dead <- init_dead[order(init_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_dead_p05 <- init_dead[(init_dead$padj < alpha & !is.na(init_dead$padj)), ]
  
# Check to see if any padj is below alpha
  init_dead_p05
  
# Aged pollen vs larvae
  
# Extract results from differential abundance table for aged pollen vs larvae
  final_larva <- DESeq2::results(desq_dds, contrast = c("sample_type", "aged pollen", "larva"))
  
# Order differential abundances by their padj value
  final_larva <- final_larva[order(final_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_larva_p05 <- final_larva[(final_larva$padj < alpha & !is.na(final_larva$padj)), ]
  
# Check to see if any padj is below alpha
  final_larva_p05
  
# Aged pollen vs pre-wintering adults
  
# Extract results from differential abundance table for aged pollen vs pre-wintering adults
  final_pre <- DESeq2::results(desq_dds, contrast = c("sample_type", "aged pollen", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  final_pre <- final_pre[order(final_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_pre_p05 <- final_pre[(final_pre$padj < alpha & !is.na(final_pre$padj)), ]
  
# Check to see if any padj is below alpha
  final_pre_p05
  
# Aged pollen vs dead adults
  
# Extract results from differential abundance table for aged pollen vs dead adults
  final_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "aged pollen", "dead adult"))
  
# Order differential abundances by their padj value
  final_dead <- final_dead[order(final_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_dead_p05 <- final_dead[(final_dead$padj < alpha & !is.na(final_dead$padj)), ]
  
# Check to see if any padj is below alpha
  final_dead_p05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva_pre <- DESeq2::results(desq_dds, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva_pre <- larva_pre[order(larva_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_pre_p05 <- larva_pre[(larva_pre$padj < alpha & !is.na(larva_pre$padj)), ]
  
# Check to see if any padj is below alpha
  larva_pre_p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "larva", "dead adult"))
  
# Order differential abundances by their padj value
  larva_dead <- larva_dead[order(larva_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_dead_p05 <- larva_dead[(larva_dead$padj < alpha & !is.na(larva_dead$padj)), ]
  
# Check to see if any padj is below alpha
  larva_dead_p05
  
# Only bee samples
  
# Convert from a phyloseq to a deseq obj
  desq_obj_bee <- phyloseq::phyloseq_to_deseq2(rareps_bact_bee, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
  
# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq_obj_bee), 1, gm_mean)
  
# Estimate size factors
  desq_dds_bee <- DESeq2::estimateSizeFactors(desq_obj_bee, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds_bee <- DESeq2::DESeq(desq_dds_bee, fitType = "local")
  
# Set significance factor  
  alpha <- 0.05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva_pre_bee <- DESeq2::results(desq_dds_bee, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva_pre_bee <- larva_pre_bee[order(larva_pre_bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_pre_bee_p05 <- larva_pre_bee[(larva_pre_bee$padj < alpha & !is.na(larva_pre_bee$padj)), ]
  
# Check to see if any padj is below alpha
  larva_pre_bee_p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva_dead_bee <- DESeq2::results(desq_dds_bee, contrast = c("sample_type", "larva", "dead adult"))
  
# Order differential abundances by their padj value
  larva_dead_bee <- larva_dead_bee[order(larva_dead_bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_dead_bee_p05 <- larva_dead_bee[(larva_dead_bee$padj < alpha & !is.na(larva_dead_bee$padj)), ]
  
# Check to see if any padj is below alpha
  larva_dead_bee_p05

# Pre-wintering vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  pre_dead_bee <- DESeq2::results(desq_dds_bee, contrast = c("sample_type", "pre.wintering.adult", "dead adult"))
  
# Order differential abundances by their padj value
  pre_dead_bee <- pre_dead_bee[order(pre_dead_bee$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_dead_bee_p05 <- pre_dead_bee[(pre_dead_bee$padj < alpha & !is.na(pre_dead_bee$padj)), ]
  
# Check to see if any padj is below alpha
  pre_dead_bee_p05
  
