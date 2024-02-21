##### Project: Osmia Developmental Microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of ITS rRNA gene data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(phyloseq) # Version 1.44.0
  library(plotrix) # Version 3.8-4
  library(decontam) # Version 1.20.0
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
  metaITS_dev <- read.csv("Osmia_dev_master - ITS_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples <- data.frame(metaITS_dev)
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
  df <- as.data.frame(sample_data(ps1)) # Put sample_data into a ggplot-friendly data.frame
  df$LibrarySize <- sample_sums(ps1)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))

# Plot samples by library size
  ggplot(data = df, aes(x = Index, y = LibrarySize, color = sample_type)) + 
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

# Determine which ASVs are contaminants based on prevalence (presence/absence) higher than 0.5 in negative controls
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

# Make phyloseq object of presence-absence in negative controls
  ps.neg <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "control", ps1)

# Calculate taxa abundance in samples from sample counts
  ps.neg.presence <- phyloseq::transform_sample_counts(ps.neg, function(abund) 1*(abund > 0))

# Make phyloseq object of presence-absence in true positive samples
  ps.pos <- phyloseq::prune_samples(sample_data(ps1)$sample_or_control == "sample", ps1)

# Calculate taxa abundance in samples from sample counts
  ps.pos.presence <- phyloseq::transform_sample_counts(ps.pos, function(abund) 1*(abund > 0))

# Make data.frame of prevalence in positive and negative samples
  df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence), 
                        prevalence.neg = taxa_sums(ps.neg.presence),
                        contam.comb05 = contamdf.comb05$contaminant)

# Plot
  ggplot(data = df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.comb05)) + 
    geom_point() +
    xlab("Prevalence (Controls)") +
    ylab("Prevalence (Samples)")

# Make a new phyloseq object without contaminant taxa 
  ps.noncontam <- phyloseq::prune_taxa(!contamdf.comb05$contaminant, ps1)
  ps.noncontam

# Remove control samples used for identifying contaminants
  ps_sub <- phyloseq::subset_samples(ps.noncontam, sample_or_control != "control")
  ps_sub

# Remove DNA from Pseudogymnoascus, the causative agent of white-nose syndrome in bats
  ps2 <- ps_sub %>%
    phyloseq::subset_taxa(Genus != "g__Pseudogymnoascus")
  
# Remove samples without any reads  
  ps3 <- phyloseq::prune_samples(sample_sums(ps2) != 0, ps2)

# Display total number of reads, mean, and se in phyloseq obj after processing
  sum(sample_sums(ps3))
  mean(sample_sums(ps3))
  print(plotrix::std.error(sample_sums(ps3)))
  
# Save sample metadata
  meta <- sample_data(ps3)
  
# How many samples for each developmental stage?  
  meta %>%
    group_by(sample_type) %>%
    summarise(N = n())
  
# Remove patterns in tax_table   
  tax_table(ps3)[, colnames(tax_table(ps3))] <- gsub(tax_table(ps3)[, colnames(tax_table(ps3))], pattern = "[a-z]__", replacement = "")
  
# Save taxonomic and ASV counts
  write.csv(tax_table(ps3), "Osmia_dev_ITStaxa.csv")
  write.csv(otu_table(ps3), "Osmia_dev_ITSotu.csv")
  
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
  
# Calculate species richness
  fungrich <- phyloseq::estimate_richness(ps3, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))

# Build df with metadata
  fungrich$sample_type <- sample_data(ps3)$sample_type
  fungrich$nesting_tube <- sample_data(ps3)$nesting_tube

# Plot species richness  
  phyloseq::plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
              theme_bw()

# Remove samples with 0 species richness 
  fungrich[fungrich == 0] <- NA
  fungrich <- fungrich[complete.cases(fungrich), ]

# Examine the effects of sample_type on Shannon richness
  mod7 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = fungrich)
  anova(mod7)

# Examine the effects of sample_type on Simpson richness
  mod8 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = fungrich)
  anova(mod8)

# Examine the effects of sample_type on observed richness
  mod9 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = fungrich)
  anova(mod9)

# Order samples on x-axis
  fungrich$sample_type <- factor(fungrich$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "emerged", "dead"))

# Boxplot of Shannon index
  Osmia_dev_Shannon_fungi <- ggplot(fungrich, aes(x = sample_type, y = Shannon, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
                                labs(title = "B") +
                                xlab("Sample type") +
                                ylab("Shannon index")
  Osmia_dev_Shannon_fungi

# Boxplot of Simpson index
  Osmia_dev_Simpson_fungi <- ggplot(fungrich, aes(x = sample_type, y = Simpson, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
                                labs(title = "B") +
                                xlab("Sample type") +
                                ylab("Simpson index")
  Osmia_dev_Simpson_fungi

# Boxplot of Observed richness
  Osmia_dev_Observed_fungi <- ggplot(fungrich, aes(x = sample_type, y = Observed, color = sample_type)) + 
                                geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                                geom_jitter(size = 1, alpha = 0.9) +
                                theme_bw() +
                                theme(legend.position = "none") +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                scale_color_manual(name = "Developmental Stage",
                                                   values = c("#FDD835", "#E4511E", "#43A047", "#0288D1","#9575CD", "#616161")) +
                                labs(title = "B") +
                                xlab("Sample type") +
                                ylab("Observed richness")
  Osmia_dev_Observed_fungi

# Only bee samples
  
# Subset phyloseq object to only include samples from larvae, pre-wintering adults, emerged adults, and dead adults
  ps4 <- phyloseq::subset_samples(ps3, sample_type != "initial provision")
  ps4 <- phyloseq::subset_samples(ps4, sample_type != "final provision")
  ps4

# Build df with metadata
  fungrich_bee <- phyloseq::estimate_richness(ps4, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  fungrich_bee$sample_type <- sample_data(ps4)$sample_type
  fungrich_bee$nesting_tube <- sample_data(ps4)$nesting_tube
  
# Plot species richness  
  phyloseq::plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "nesting_tube") + 
    theme_bw()
  
# Remove samples with 0 species richness 
  fungrich_bee[fungrich_bee == 0] <- NA
  fungrich_bee <- fungrich_bee[complete.cases(fungrich_bee), ]
  
# Examine the effects of sample_type on Shannon richness
  mod10 <- nlme::lme(Shannon ~ sample_type, random = ~1|nesting_tube, data = fungrich_bee)
  anova(mod10)
  
# Examine the effects of sample_type on Simpson richness
  mod11 <- nlme::lme(Simpson ~ sample_type, random = ~1|nesting_tube, data = fungrich_bee)
  anova(mod11)
  
# Examine the effects of sample_type on observed richness
  mod12 <- nlme::lme(Observed ~ sample_type, random = ~1|nesting_tube, data = fungrich_bee)
  anova(mod12)
  
## Beta diversity with relative abundance data ----  
  
# All pollen and bee samples  
  
# Calculate the relative abundance of each otu
  ps.prop_fung <- phyloseq::transform_sample_counts(ps3, function(otu) otu/sum(otu))
  ps.prop_fung
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung_bray <- phyloseq::distance(ps.prop_fung, method = "bray")
  
# Convert to data frame
  samplefung <- data.frame(sample_data(ps3))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm <- vegan::adonis2(fung_bray ~ sample_type, data = samplefung)
  fung_perm
  
# Follow up with pairwise comparisons - which sample types differ?
  #fungi_perm_BH <- RVAideMemoire::pairwise.perm.manova(fung_bray, samplefung$sample_type, p.method = "BH")
  #fungi_perm_BH
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_relabund <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = samplefung$nesting_tube,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition, dealing with pseudoreplication
  fung_perm_pseudo <- vegan::adonis2(fung_bray ~ sample_type, permutations = perm_relabund, data = samplefung)
  fung_perm_pseudo
  
# Only bee samples
  
# Calculate the relative abundance of each otu
  ps.prop_fung_bee <- phyloseq::transform_sample_counts(ps4, function(otu) otu/sum(otu))
  ps.prop_fung_bee
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung_bray_bee <- phyloseq::distance(ps.prop_fung_bee, method = "bray")
  
# Convert to data frame
  samplefung_bee <- data.frame(sample_data(ps4))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm_bee <- vegan::adonis2(fung_bray_bee ~ sample_type, data = samplefung_bee)
  fung_perm_bee
  
# Follow up with pairwise comparisons - which sample types differ?
  #fungi_perm_BH_bee <- RVAideMemoire::pairwise.perm.manova(fung_bray_bee, samplefung_bee$sample_type, p.method = "BH")
  #fungi_perm_BH_bee
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_relabund_bee <- how(within = Within(type = "free"),
                           plots = Plots(type = "none"),
                           blocks = samplefung_bee$nesting_tube,
                           observed = FALSE,
                           complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition, dealing with pseudoreplication
  fung_perm_pseudo_bee <- vegan::adonis2(fung_bray_bee ~ sample_type, permutations = perm_relabund_bee, data = samplefung_bee)
  fung_perm_pseudo_bee
  
## Test for homogeneity of multivariate dispersion with relative abundance data ----
  
# All pollen and bee samples  
  
# Calculate the average distance of group members to the group centroid
  disp_fung <- vegan::betadisper(fung_bray, samplefung$sample_type)
  disp_fung
  
# Do any of the group dispersions differ?
  disp_fung_an <- anova(disp_fung)
  disp_fung_an
  
# Which group dispersions differ?
  #disp_fung_ttest <- vegan::permutest(disp_fung, 
                                      #control = permControl(nperm = 999),
                                      #pairwise = TRUE)
 #disp_fung_ttest
  
# Which group dispersions differ?
  #disp_fung_tHSD <- stats::TukeyHSD(disp_fung)
  #disp_fung_tHSD
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_fung_bee <- vegan::betadisper(fung_bray_bee, samplefung_bee$sample_type)
  disp_fung_bee
  
# Do any of the group dispersions differ?
  disp_fung_an_bee <- anova(disp_fung_bee)
  disp_fung_an_bee
  
# Which group dispersions differ?
  disp_fung_ttest_bee <- vegan::permutest(disp_fung_bee, 
                                          control = permControl(nperm = 999),
                                          pairwise = TRUE)
  disp_fung_ttest_bee
  
# Which group dispersions differ?
  disp_fung_tHSD_bee <- stats::TukeyHSD(disp_fung_bee)
  disp_fung_tHSD_bee
  
## Plot distance to centroid ----

# Create df with sample metadata
  #sam_dat <- as.data.frame(sample_data(ps2))
  
# Create df with distance to centroid measures
  #disp_fung_df <- as.data.frame(disp_fung$distances)
  #disp_fung_df$extractionID <- row.names(disp_fung_df)
  
# Merge dfs
  #disp_df <- merge(sam_dat, disp_fung_df, by = "extractionID")
  
# Plot
  #ggplot(disp_df, aes(x = sample_type, y = disp_fung$distance, color = sample_type)) + 
  #geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) + 
  #geom_jitter(size = 1, alpha = 0.9) +
  #theme_bw() +
  #theme(legend.position = "none") +
  #theme(panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank()) +
  #scale_color_manual(values = c("#FDD835", "#E4511E", "#43A047", "#0288D1", "#616161")) +
  #labs(title = "B") +
  #xlab("Sample type") +
  #ylab("Distance to centroid")
  
## Ordination with relative abundance data ----
  
# All pollen and bee samples  
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- phyloseq::ordinate(ps.prop_fung, method = "PCoA", distance = "bray")
  
# Plot ordination
  Osmia_dev_PCoA_fungi <- plot_ordination(ps.prop_fung, ord.pcoa.bray, color = "sample_type") + 
                            theme_bw() +
                            theme(text = element_text(size = 16)) +
                            theme(legend.justification = "left", 
                                  legend.title = element_text(size = 16, colour = "black"), 
                                  legend.text = element_text(size = 14, colour = "black")) + 
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            geom_point(size = 3) +
                            scale_color_manual(values = c("#616161", "#9575CD", "#E4511E", "#FDD835", "#43A047", "#0288D1")) +
                            labs(title = "B")
  Osmia_dev_PCoA_fungi  
  
# Only bee samples
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_bee <- phyloseq::ordinate(ps.prop_fung_bee, method = "PCoA", distance = "bray")
  
# Plot ordination
  Osmia_dev_PCoA_fungi_bee <- plot_ordination(ps.prop_fung_bee, ord.pcoa.bray_bee, color = "sample_type") + 
                                  theme_bw() +
                                  theme(text = element_text(size = 16)) +
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 16, colour = "black"), 
                                        legend.text = element_text(size = 14, colour = "black")) + 
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  geom_point(size = 3) +
                                  scale_color_manual(values = c("#616161", "#9575CD", "#E4511E", "#FDD835", "#43A047", "#0288D1")) +
                                  labs(title = "B")
  Osmia_dev_PCoA_fungi_bee
  
## Rarefaction ----

# All pollen and bee samples

# Produce rarefaction curves
  tab <- otu_table(ps3)
  class(tab) <- "matrix"
  tab <- t(tab)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_fungi <- vegan::rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia_dev_rare_fungi <- ggplot(rare_tidy_fungi, aes(x = Sample, y = Species, group = Site)) +
                            geom_line() +
                            theme_bw() +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            labs(title = "B") + 
                            xlab("Number of reads") +
                            ylab("Number of species")
  Osmia_dev_rare_fungi

# Set seed and rarefy
  set.seed(1234)
  fung_rareps <- phyloseq::rarefy_even_depth(ps3, sample.size = 20)

# Bee samples only
  
# Produce rarefaction curves
  tab_bee <- otu_table(ps4)
  class(tab_bee) <- "matrix"
  tab_bee <- t(tab_bee)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_fungi_bee <- vegan::rarecurve(tab_bee, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  Osmia_dev_rare_fungi_bee <- ggplot(rare_tidy_fungi_bee, aes(x = Sample, y = Species, group = Site)) +
                                geom_line() +
                                theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank()) +
                                labs(title = "B") + 
                                xlab("Number of reads") +
                                ylab("Number of species")
  Osmia_dev_rare_fungi_bee
  
# Set seed and rarefy
  set.seed(1234)
  fung_rareps_bee <- phyloseq::rarefy_even_depth(ps4, sample.size = 20)

## Beta diversity with rarefied data ----  
  
# All pollen and bee samples  
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung_bray_rare <- phyloseq::distance(fung_rareps, method = "bray")
  
# Convert to data frame
  samplefung_rare <- data.frame(sample_data(fung_rareps))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm_rare <- vegan::adonis2(fung_bray_rare ~ sample_type, data = samplefung_rare)
  fung_perm_rare
  
# Follow up with pairwise comparisons - which sample types differ?
  #fungi_perm_BH_rare <- RVAideMemoire::pairwise.perm.manova(fung_bray_rare, samplefung_rare$sample_type, p.method = "BH")
  #fungi_perm_BH_rare  
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_rare <- how(within = Within(type = "free"),
                   plots = Plots(type = "none"),
                   blocks = samplefung_rare$nesting_tube,
                   observed = FALSE,
                   complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm_rare_pseudo <- vegan::adonis2(fung_bray_rare ~ sample_type, permutations = perm_rare, data = samplefung_rare)
  fung_perm_rare_pseudo
  
# Only bee samples
  
# Create a distance matrix using Bray Curtis dissimilarity
  fung_bray_rare_bee <- phyloseq::distance(fung_rareps_bee, method = "bray")
  
# Convert to data frame
  samplefung_rare_bee <- data.frame(sample_data(fung_rareps_bee))
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm_rare_bee <- vegan::adonis2(fung_bray_rare_bee ~ sample_type, data = samplefung_rare_bee)
  fung_perm_rare_bee
  
# Follow up with pairwise comparisons - which sample types differ?
  #fungi_perm_BH_rare_bee <- RVAideMemoire::pairwise.perm.manova(fung_bray_rare_bee, samplefung_rare_bee$sample_type, p.method = "BH")
  #fungi_perm_BH_rare_bee  
  
# Set permutations to deal with pseudoreplication of bee nests
  perm_rare_bee <- how(within = Within(type = "free"),
                       plots = Plots(type = "none"),
                       blocks = samplefung_rare_bee$nesting_tube,
                       observed = FALSE,
                       complete = FALSE)
  
# Perform the PERMANOVA to test effects of developmental stage on fungal community composition
  fung_perm_rare_pseudo <- vegan::adonis2(fung_bray_rare_bee ~ sample_type, permutations = perm_rare_bee, data = samplefung_rare_bee)
  fung_perm_rare_pseudo

## Test for homogeneity of multivariate dispersion with rarefied data ----
  
# All pollen and bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_fung_rare <- vegan::betadisper(fung_bray_rare, samplefung_rare$sample_type)
  disp_fung_rare
  
# Do any of the group dispersions differ?
  disp_fung_an_rare <- anova(disp_fung_rare)
  disp_fung_an_rare
  
# Which group dispersions differ?
  disp_fung_ttest_rare <- vegan::permutest(disp_fung_rare, 
                                           control = permControl(nperm = 999),
                                           pairwise = TRUE)
  disp_fung_ttest_rare
  
# Which group dispersions differ?
  disp_fung_tHSD_rare <- stats::TukeyHSD(disp_fung_rare)
  disp_fung_tHSD_rare
  
# Only bee samples
  
# Calculate the average distance of group members to the group centroid
  disp_fung_rare_bee <- vegan::betadisper(fung_bray_rare_bee, samplefung_rare_bee$sample_type)
  disp_fung_rare_bee
  
# Do any of the group dispersions differ?
  disp_fung_an_rare_bee <- anova(disp_fung_rare_bee)
  disp_fung_an_rare_bee
  
# Which group dispersions differ?
  disp_fung_ttest_rare_bee <- vegan::permutest(disp_fung_rare_bee, 
                                               control = permControl(nperm = 999),
                                               pairwise = TRUE)
  disp_fung_ttest_rare_bee
  
# Which group dispersions differ?
  disp_fung_tHSD_rare_bee <- stats::TukeyHSD(disp_fung_rare_bee)
  disp_fung_tHSD_rare_bee

## Ordination with rarefied data ----

# All pollen and bee samples  
  
# Calculate the relative abundance of each otu
  ps.prop_rare <- phyloseq::transform_sample_counts(fung_rareps, function(otu) otu/sum(otu))
  ps.prop_rare

# PCoA using Bray-Curtis distance
  ord.pcoa.bray_rare <- phyloseq::ordinate(ps.prop_rare, method = "PCoA", distance = "bray")

# Plot ordination
  Osmia_dev_PCoA_fungi_rare <- plot_ordination(ps.prop_rare, ord.pcoa.bray_rare, color = "sample_type") + 
                                  theme_bw() +
                                  theme(text = element_text(size = 16)) +
                                  theme(legend.justification = "left", 
                                        legend.title = element_text(size = 16, colour = "black"), 
                                        legend.text = element_text(size = 14, colour = "black")) + 
                                  theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
                                  geom_point(size = 3) +
                                  scale_color_manual(values = c("#616161", "#9575CD", "#E4511E", "#FDD835", "#43A047", "#0288D1")) +
                                  labs(title = "B")
  Osmia_dev_PCoA_fungi_rare
  
# Only bee samples
  
# Calculate the relative abundance of each otu
  ps.prop_rare_bee <- phyloseq::transform_sample_counts(fung_rareps_bee, function(otu) otu/sum(otu))
  ps.prop_rare_bee
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray_rare_bee <- phyloseq::ordinate(ps.prop_rare_bee, method = "PCoA", distance = "bray")
  
# Plot ordination
  Osmia_dev_PCoA_fungi_rare_bee <- plot_ordination(ps.prop_rare_bee, ord.pcoa.bray_rare_bee, color = "sample_type") + 
                                      theme_bw() +
                                      theme(text = element_text(size = 16)) +
                                      theme(legend.justification = "left", 
                                            legend.title = element_text(size = 16, colour = "black"), 
                                            legend.text = element_text(size = 14, colour = "black")) + 
                                      theme(panel.grid.major = element_blank(),
                                            panel.grid.minor = element_blank()) +
                                      geom_point(size = 3) +
                                      scale_color_manual(values = c("#616161", "#9575CD", "#E4511E", "#FDD835", "#43A047", "#0288D1")) +
                                      labs(title = "B")
  Osmia_dev_PCoA_fungi_rare_bee

## Stacked community plot ----

# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Stretch palette (define more intermediate color options)
  okabe_ext <- unikn::usecol(Okabe_Ito, n = 78)
  colors <- sample(okabe_ext)

# Remove patterns in tax_table
  tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))] <- gsub(tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))], pattern = "[a-z]__", replacement = "")

# Sort data by Family
  y7 <- phyloseq::tax_glom(fung_rareps, taxrank = 'Family') # agglomerate taxa
  y8 <- phyloseq::transform_sample_counts(y7, function(x) x/sum(x))
  y9 <- phyloseq::psmelt(y8)
  y9$Family <- as.character(y9$Family)
  y9$Family[y9$Abundance < 0.01] <- "Family < 1% abund."
  y9$Family <- as.factor(y9$Family)
  head(y9)
  
# Save relative abundance data
  write.csv(y9, "Osmia_dev_Fam_fungi_relabund.csv")

# Reorder x-axis 
  y9$sample_type <- factor(y9$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "emerged", "dead"))

# Plot Family by sample type
  ggplot(data = y9, aes(x = sample_type, y = Abundance, fill = Family)) + 
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
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Fungi")

# Plot Family for each sample
  Osmia_dev_fam_relabund_fungi <- ggplot(data = y9, aes(x = sampleID, y = Abundance, fill = Family)) + 
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
                                    ggtitle("B")
  Osmia_dev_fam_relabund_fungi

# Sort data by Genus
  y10 <- phyloseq::tax_glom(fung_rareps, taxrank = 'Genus') # agglomerate taxa
  y11 <- phyloseq::transform_sample_counts(y10, function(x) x/sum(x))
  y12 <- phyloseq::psmelt(y11)
  y12$Genus <- as.character(y12$Genus)
  y12$Genus[y12$Abundance < 0.01] <- "Genera < 1% abund."
  y12$Genus <- as.factor(y12$Genus)
  head(y12)
  
# Save relative abundance data
  write.csv(y12, "Osmia_dev_Gen_fungi_relabund.csv")

# Order samples on x-axis
  y12$sample_type <- factor(y12$sample_type, levels = c("initial provision", "final provision", "larva", "pre-wintering adult", "dead"))

# Plot Genus by sample type
  ggplot(data = y12, aes(x = sample_type, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
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
    guides(fill = guide_legend(ncol = 3)) +
    ggtitle("Fungi")

# Plot Genus for each sample
  Osmia_dev_gen_relabund_fungi <- ggplot(data = y12, aes(x = sampleID, y = Abundance, fill = Genus)) + 
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
                                    ggtitle("B")
  Osmia_dev_gen_relabund_fungi

## Differential abundance with raw data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html  

# Remove patterns in tax_table   
  tax_table(ps3)[, colnames(tax_table(ps3))] <- gsub(tax_table(ps3)[, colnames(tax_table(ps3))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq::phyloseq_to_deseq2(ps3, ~ sample_type)
  
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
  
# Initial vs final provisions
  
# Extract results from differential abundance table for initial vs final provisions
  init_final <- results(desq_dds, contrast = c("sample_type", "initial provision", "final provision"))
  
# Order differential abundances by their padj value
  init_final <- init_final[order(init_final$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_final_p05 <- init_final[(init_final$padj < alpha & !is.na(init_final$padj)), ]
  
# Check to see if any padj is below alpha
  init_final_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  init_final_p05 <- cbind(as(init_final_p05, "data.frame"),
                          as(tax_table(ps3)[rownames(init_final_p05), ], "matrix"))
  
# Plot
  ggplot(init_final_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Initial vs. Final Provision")
  
# Initial provisions vs larvae
  
# Extract results from differential abundance table for initial provisions vs larvae
  init_larva <- DESeq2::results(desq_dds, contrast = c("sample_type", "initial provision", "larva"))
  
# Order differential abundances by their padj value
  init_larva <- init_larva[order(init_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_larva_p05 <- init_larva[(init_larva$padj < alpha & !is.na(init_larva$padj)), ]
  
# Check to see if any padj is below alpha
  init_larva_p05
  
# Initial provisions vs pre-wintering adults
  
# Extract results from differential abundance table for initial provisions vs pre-wintering adults
  init_pre <- DESeq2::results(desq_dds, contrast = c("sample_type", "initial provision", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  init_pre <- init_pre[order(init_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_pre_p05 <- init_pre[(init_pre$padj < alpha & !is.na(init_pre$padj)), ]
  
# Check to see if any padj is below alpha
  init_pre_p05
  
# Initial provisions vs dead adults
  
# Extract results from differential abundance table for initial provisions vs dead adults
  init_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "initial provision", "dead"))
  
# Order differential abundances by their padj value
  init_dead <- init_dead[order(init_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_dead_p05 <- init_dead[(init_dead$padj < alpha & !is.na(init_dead$padj)), ]
  
# Check to see if any padj is below alpha
  init_dead_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  init_dead_p05 <- cbind(as(init_dead_p05, "data.frame"),
                         as(tax_table(ps3)[rownames(init_dead_p05), ], "matrix"))
  
# Plot
  ggplot(init_dead_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Initial Provision vs. Dead Adults")  
  
# Final provisions vs larvae
  
# Extract results from differential abundance table for final provisions vs larvae
  final_larva <- DESeq2::results(desq_dds, contrast = c("sample_type", "final provision", "larva"))
  
# Order differential abundances by their padj value
  final_larva <- final_larva[order(final_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_larva_p05 <- final_larva[(final_larva$padj < alpha & !is.na(final_larva$padj)), ]
  
# Check to see if any padj is below alpha
  final_larva_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  final_larva_p05 <- cbind(as(final_larva_p05, "data.frame"),
                           as(tax_table(ps3)[rownames(final_larva_p05), ], "matrix"))
  
# Plot
  ggplot(final_larva_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Final Provision vs. Larvae")  
  
# Final provisions vs pre-wintering adults
  
# Extract results from differential abundance table for final provisions vs pre-wintering adults
  final_pre <- DESeq2::results(desq_dds, contrast = c("sample_type", "final provision", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  final_pre <- final_pre[order(final_pre$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  final_pre_p05 <- final_pre[(final_pre$padj < alpha & !is.na(final_pre$padj)), ]
  
# Check to see if any padj is below alpha
  final_pre_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  final_pre_p05 <- cbind(as(final_pre_p05, "data.frame"),
                         as(tax_table(ps3)[rownames(final_pre_p05), ], "matrix"))
  
# Plot
  ggplot(final_pre_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Final Provision vs. Pre-wintering Adults") 
  
# Final provisions vs dead adults
  
# Extract results from differential abundance table for final provisions vs dead adults
  final_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "final provision", "dead"))
  
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
  larva_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "larva", "dead"))
  
# Order differential abundances by their padj value
  larva_dead <- larva_dead[order(larva_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_dead_p05 <- larva_dead[(larva_dead$padj < alpha & !is.na(larva_dead$padj)), ]
  
# Check to see if any padj is below alpha
  larva_dead_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  larva_dead_p05 <- cbind(as(larva_dead_p05, "data.frame"),
                         as(tax_table(ps3)[rownames(larva_dead_p05), ], "matrix"))
  
  # Plot
  ggplot(larva_dead_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Final Provision vs. Dead Adults")
  
# Pre-wintering vs emerged adults
  
# Extract results from differential abundance table for pre-wintering vs emerged adults
  pre_emerg <- DESeq2::results(desq_dds, contrast = c("sample_type", "pre.wintering.adult", "emerged"))
  
# Order differential abundances by their padj value
  pre_emerg <- pre_emerg[order(pre_emerg$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_emerg_p05 <- pre_emerg[(pre_emerg$padj < alpha & !is.na(pre_emerg$padj)), ]
  
# Check to see if any padj is below alpha
  pre_emerg_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  pre_emerg_p05 <- cbind(as(pre_emerg_p05, "data.frame"),
                          as(tax_table(ps3)[rownames(pre_emerg_p05), ], "matrix"))
  
# Plot
  ggplot(pre_emerg_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Pre-wintering vs. Emerged Adults")
  
# Pre-wintering vs dead adults
  
# Extract results from differential abundance table for pre-wintering vs dead adults
  pre_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "pre.wintering.adult", "dead"))
  
# Order differential abundances by their padj value
  pre_dead <- pre_dead[order(pre_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_dead_p05 <- pre_dead[(pre_dead$padj < alpha & !is.na(pre_dead$padj)), ]
  
# Check to see if any padj is below alpha
  pre_dead_p05
  
# Emerged vs dead adults
  
# Extract results from differential abundance table for emerged vs dead
  emerg_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "emerged", "dead"))
  
# Order differential abundances by their padj value
  emerg_dead <- emerg_dead[order(emerg_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_dead_p05 <- emerg_dead[(emerg_dead$padj < alpha & !is.na(emerg_dead$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_dead_p05
  
# Emerged adults vs initial provisions
  
# Extract results from differential abundance table for emerged vs initial provisions
  emerg_init <- DESeq2::results(desq_dds, contrast = c("sample_type", "emerged", "initial provision"))
  
# Order differential abundances by their padj value
  emerg_init <- emerg_init[order(emerg_init$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_init_p05 <- emerg_init[(emerg_init$padj < alpha & !is.na(emerg_init$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_init_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  emerg_init_p05 <- cbind(as(emerg_init_p05, "data.frame"),
                         as(tax_table(ps3)[rownames(emerg_init_p05), ], "matrix"))
  
# Plot
  ggplot(emerg_init_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Emerged Adults vs. Initial Provisions")
  
# Emerged adults vs larvae
  
# Extract results from differential abundance table for emerged vs initial provisions
  emerg_larva <- DESeq2::results(desq_dds, contrast = c("sample_type", "emerged", "larva"))
  
# Order differential abundances by their padj value
  emerg_larva <- emerg_larva[order(emerg_larva$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_larva_p05 <- emerg_larva[(emerg_larva$padj < alpha & !is.na(emerg_larva$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_larva_p05
  
# Combine filtered differential abundance data with taxonomic names from phyloseq obj
  emerg_larva_p05 <- cbind(as(emerg_larva_p05, "data.frame"),
                          as(tax_table(ps3)[rownames(emerg_larva_p05), ], "matrix"))
  
# Plot
  ggplot(emerg_larva_p05, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
    geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    theme_bw() +
    ggtitle("Emerged Adults vs. Larvae")
  
## Differential abundance with rarefied data ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html
  
# Remove patterns in tax_table   
  tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))] <- gsub(tax_table(fung_rareps)[, colnames(tax_table(fung_rareps))], pattern = "[a-z]__", replacement = "")
  
# Convert from a phyloseq to a deseq obj
  desq_obj_rare <- phyloseq::phyloseq_to_deseq2(fung_rareps, ~ sample_type)
  
# Calculate the geometric mean and remove rows with NA
  gm_mean <- function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }

# Add a count of 1 to all geometric means
  geoMeans <- apply(counts(desq_obj_rare), 1, gm_mean)
  
# Estimate size factors
  desq_dds_rare <- DESeq2::estimateSizeFactors(desq_obj_rare, geoMeans = geoMeans)
  
# Fit a local regression
  desq_dds_rare <- DESeq2::DESeq(desq_dds_rare, fitType = "local")

# Set significance factor  
  alpha <- 0.05
  
# Initial provisions vs larvae
  
# Extract results from differential abundance table for initial provisions vs larvae
  init_larva_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "initial provision", "larva"))
  
# Order differential abundances by their padj value
  init_larva_rare <- init_larva_rare[order(init_larva_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_larva_rare_p05 <- init_larva_rare[(init_larva_rare$padj < alpha & !is.na(init_larva_rare$padj)), ]
  
# Check to see if any padj is below alpha
  init_larva_rare_p05
  
# Initial provisions vs pre-wintering adults
  
# Extract results from differential abundance table for initial provisions vs pre-wintering adults
  init_pre_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "initial provision", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  init_pre_rare <- init_pre_rare[order(init_pre_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_pre_rare_p05 <- init_pre_rare[(init_pre_rare$padj < alpha & !is.na(init_pre_rare$padj)), ]
  
# Check to see if any padj is below alpha
  init_pre_rare_p05
  
# Initial provisions vs dead adults
  
# Extract results from differential abundance table for initial provisions vs dead adults
  init_dead_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "initial provision", "dead"))
  
# Order differential abundances by their padj value
  init_dead_rare <- init_dead_rare[order(init_dead_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  init_dead_rare_p05 <- init_dead_rare[(init_dead_rare$padj < alpha & !is.na(init_dead_rare$padj)), ]
  
# Check to see if any padj is below alpha
  init_dead_rare_p05
  
# Larvae vs pre-wintering adults
  
# Extract results from differential abundance table for larvae vs pre-wintering adults
  larva_pre_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "larva", "pre.wintering.adult"))
  
# Order differential abundances by their padj value
  larva_pre_rare <- larva_pre_rare[order(larva_pre_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_pre_rare_p05 <- larva_pre[(larva_pre_rare$padj < alpha & !is.na(larva_pre_rare$padj)), ]
  
# Check to see if any padj is below alpha
  larva_pre_rare_p05
  
# Larvae vs dead adults
  
# Extract results from differential abundance table for larvae vs dead adults
  larva_dead_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "larva", "dead"))
  
# Order differential abundances by their padj value
  larva_dead_rare <- larva_dead_rare[order(larva_dead_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  larva_dead_rare_p05 <- larva_dead_rare[(larva_dead_rare$padj < alpha & !is.na(larva_dead_rare$padj)), ]
  
# Check to see if any padj is below alpha
  larva_dead_rare_p05
  
# Pre-wintering vs emerged adults
  
# Extract results from differential abundance table for pre-wintering vs emerged adults
  pre_emerg_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "pre.wintering.adult", "emerged"))
  
# Order differential abundances by their padj value
  pre_emerg_rare <- pre_emerg_rare[order(pre_emerg_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_emerg_rare_p05 <- pre_emerg_rare[(pre_emerg_rare$padj < alpha & !is.na(pre_emerg_rare$padj)), ]
  
# Check to see if any padj is below alpha
  pre_emerg_rare_p05
  
# Pre-wintering vs dead adults
  
# Extract results from differential abundance table for pre-wintering vs dead adults
  pre_dead_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "pre.wintering.adult", "dead"))
  
# Order differential abundances by their padj value
  pre_dead_rare <- pre_dead_rare[order(pre_dead_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  pre_dead_rare_p05 <- pre_dead_rare[(pre_dead_rare$padj < alpha & !is.na(pre_dead_rare$padj)), ]
  
# Check to see if any padj is below alpha
  pre_dead_rare_p05

# Emerged vs dead adults
  
# Extract results from differential abundance table for emerged vs dead
  emerg_dead <- DESeq2::results(desq_dds, contrast = c("sample_type", "emerged", "dead"))
  
# Order differential abundances by their padj value
  emerg_dead <- emerg_dead[order(emerg_dead$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_dead_p05 <- emerg_dead[(emerg_dead$padj < alpha & !is.na(emerg_dead$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_dead_p05
  
# Emerged adults vs initial provisions
  
# Extract results from differential abundance table for emerged vs initial provisions
  emerg_init_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "emerged", "initial provision"))
  
# Order differential abundances by their padj value
  emerg_init_rare <- emerg_init_rare[order(emerg_init_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_init_rare_p05 <- emerg_init_rare[(emerg_init_rare$padj < alpha & !is.na(emerg_init_rare$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_init_rare_p05
  
# Emerged adults vs larvae
  
# Extract results from differential abundance table for emerged vs initial provisions
  emerg_larva_rare <- DESeq2::results(desq_dds_rare, contrast = c("sample_type", "emerged", "larva"))
  
# Order differential abundances by their padj value
  emerg_larva_rare <- emerg_larva_rare[order(emerg_larva_rare$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  emerg_larva_rare_p05 <- emerg_larva_rare[(emerg_larva_rare$padj < alpha & !is.na(emerg_larva_rare$padj)), ]
  
# Check to see if any padj is below alpha
  emerg_larva_rare_p05
  
