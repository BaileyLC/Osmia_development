##### Project: Osmia climate change

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of 16S rRNA gene data

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
  seqtab.nochim <- readRDS("OsmiaCC_seqs16Sall.rds")
  taxa <- readRDS("OsmiaCC_taxa16Sall.rds")
  meta16S_CC <- read.csv("OsmiaCC_master - 16S_worked.csv")

## Create phyloseq object ----

# Re-create your df
  samples.out <- rownames(seqtab.nochim)
  samples <- data.frame(meta16S_CC)
  extractionID <- samples$extractionID
  sample_type <- samples$sample_type
  sampleID <- samples$sampleID
  temp_treat <- samples$temp_treat
  micro_treat <- samples$micro_treat
  combo_treat <- samples$combo_treat
  sample_or_control <- samples$sample_or_control
  sex <- samples$sex
  sampleinfo <- data.frame(extractionID = extractionID, 
                           sample_type = sample_type, 
                           sampleID = sampleID,  
                           temp_treat = temp_treat, 
                           micro_treat = micro_treat, 
                           combo_treat = combo_treat, 
                           sample_or_control = sample_or_control,
                           sex = sex)
  rownames(sampleinfo) <- samples.out
  
# Format your data to work with phyloseq
  ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), sample_data(sampleinfo), tax_table(taxa))
  ps1
  
## Inspect & remove contaminants ----
# Resource: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
  
# Create df with LibrarySize for each sample
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
  
# What remains in the phyloseq object?
  ps2
  
# Remove samples without any reads
  ps3 <- prune_samples(sample_sums(ps2) != 0, ps2)
  ps3
  
# Remove samples from females
  ps3 <- prune_samples(sample_data(ps3)$sex != "F", ps3)
  ps3
  
# How many reads are in each sample? 
  sample_sums(ps3)
  
# What is the mean number of reads in all samples?
  mean(sample_sums(ps3))
  
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
  
# Plot number of reads per sample
  ggplot(readsumsdf, aes(x = sorted, y = nreads)) + 
    geom_bar(stat = "identity") +
    ggtitle("Total number of reads") + 
    scale_y_log10() + 
    facet_wrap(~ type, 1, scales = "free")
  
## Species richness ----  
  
# Estimate Shannon, Simpson & observed richness
  bactrich <- estimate_richness(ps3, split = TRUE, measures = c("Shannon", "Simpson", "Observed"))
  
# Build df with metadata 
  bactrich$sampleID <- sample_data(ps3)$sampleID
  bactrich$sample_type <- sample_data(ps3)$sample_type
  bactrich$temp_treat <- sample_data(ps3)$temp_treat
  bactrich$micro_treat <- sample_data(ps3)$micro_treat
  bactrich$combo_treat <- sample_data(ps3)$combo_treat
  
# Plot Shannon, Simpson & observed richness  
  plot_richness(ps3, x = "sample_type", measures = c("Shannon", "Simpson", "Observed"), color = "combo_treat") + 
    theme_bw() +
    xlab("")
  
# Remove samples with 0 species richness 
  bactrich[bactrich == 0] <- NA
  bactrich <- bactrich[complete.cases(bactrich), ]
  
# Examine interactive effects of temperature and microbiome treatments on Shannon diversity
  mod1 <- lme(Shannon ~ temp_treat*micro_treat, random = ~1|sampleID, data = bactrich)
  anova(mod1)
  
# Examine interactive effects of temperature and microbiome treatments on Simpson diversity
  mod2 <- lme(Simpson ~ temp_treat*micro_treat, random = ~1|sampleID, data = bactrich)
  anova(mod2)
  
# Examine interactive effects of temperature and microbiome treatments on observed richness
  mod3 <- lme(Observed ~ temp_treat*micro_treat, random = ~1|sampleID, data = bactrich)
  anova(mod3)
  
# Reorder x-axis
  bactrich$combo_treat <- factor(bactrich$combo_treat, levels = c("CS", "CN", "AS", "AN", "WS", "WN"))

# New names for facet_grid
  type_names <- c('final provision' = "provisions with bee",
                  'provision w/o bee' = "provisions without bee", 
                  'initial provision' = "initial provisions")

# Boxplot of Shannon index
  OsmiaCC_Shannon_bact <- ggplot(bactrich, aes(x = combo_treat, y = Shannon, fill = combo_treat)) + 
                            geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                            geom_jitter(size = 1, alpha = 0.9) +
                            theme_bw() +
                            theme(legend.position = "none") +
                            theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank()) +
                            facet_grid(~ sample_type,
                                       scale = "free",
                                       space = "free",
                                       labeller = as_labeller(type_names)) +
                            scale_fill_manual(name = "Treatment", 
                                              values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                              labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                            labs(title = "A") +
                            xlab("Treatment") +
                            ylab("Shannon index")
  OsmiaCC_Shannon_bact
  
# Boxplot of Simpson index
  OsmiaCC_Simpson_bact <- ggplot(bactrich, aes(x = combo_treat, y = Simpson, fill = combo_treat)) + 
                            geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                            geom_jitter(size = 1, alpha = 0.9) +
                            theme_bw() +
                            theme(legend.position = "none") +
                            facet_grid(~ sample_type,
                                       scale = "free",
                                       space = "free",
                                       labeller = as_labeller(type_names)) +
                            scale_fill_manual(name = "Treatment", 
                                              values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                              labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                            labs(title = "A") + 
                            xlab("Treatment") +
                            ylab("Simpson index")
  OsmiaCC_Simpson_bact
  
# Boxplot of Observed richness
  OsmiaCC_Observed_bact <- ggplot(bactrich, aes(x = combo_treat, y = Observed, fill = combo_treat)) + 
                              geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(width = 0.1)) +
                              geom_jitter(size = 1, alpha = 0.9) +
                              theme_bw() +
                              theme(legend.position = "none") +
                              facet_grid(~ sample_type,
                                         scale = "free",
                                         space = "free",
                                         labeller = as_labeller(type_names)) +
                              scale_fill_manual(name = "Treatment",
                                                values = c("#64B5F6","#1565C0", "#9E9E9E", "#616161", "#E57373", "#C62828"),
                                                labels = c('Cool: Sterile', 'Cool: Natural', 'Ambient: Sterile', 'Ambient: Natural', 'Warm: Sterile', 'Warm: Natural')) +
                              xlab("Treatment") +
                              ylab("Observed richness") +
                              ggtitle("A")
  OsmiaCC_Observed_bact
  
## Rarefaction ----
  
# Produce rarefaction curves
  tab <- otu_table(ps3)
  class(tab) <- "matrix"
  tab <- t(tab)
  rare <- rarecurve(tab, step = 20, label = FALSE)
  
# Save rarefaction data as a "tidy" df
  rare_tidy_bact <- rarecurve(tab, label = FALSE, tidy = TRUE)
  
# Plot rarefaction curve
  OsmiaCC_rare_bact <- ggplot(rare_tidy_bact, aes(x = Sample, y = Species, group = Site)) +
                          geom_line() +
                          theme_bw() +
                          theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank()) +
                          labs(title = "A") + 
                          xlab("Number of reads") +
                          ylab("Number of species")
  OsmiaCC_rare_bact
  
# Set seed and rarefy
  set.seed(1234)
  rareps <- rarefy_even_depth(ps3, sample.size = 26)
  
# Create a distance matrix using Bray Curtis dissimilarity
  bact_bray <- phyloseq::distance(rareps, method = "bray")
  
# Convert to data frame
  samplebact <- data.frame(sample_data(rareps))
  
# Perform the PERMANOVA to test effects of treatment on bacterial community composition  
  bact_perm <- adonis2(bact_bray ~ temp_treat*micro_treat, data = samplebact)
  bact_perm
  
## Test for homogeneity of multivariate dispersion ----
  
# Calculate the average distance of group members to the group centroid
  disp_bact <- betadisper(bact_bray, samplebact$combo_treat)
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
  ps.prop <- transform_sample_counts(ps3, function(otu) otu/sum(otu))
  
# PCoA using Bray-Curtis distance
  ord.pcoa.bray <- ordinate(ps.prop, method = "PCoA", distance = "bray")
  
# Plot ordination
  OsmiaCC_PCoA_bact <- plot_ordination(ps.prop, ord.pcoa.bray, color = "combo_treat", shape = "sample_type") + 
                          theme_bw() +
                          theme(legend.position = "none") +
                          theme(text = element_text(size = 16)) +
                          theme(legend.justification = "left", 
                                legend.title = element_text(size = 16, colour = "black"), 
                                legend.text = element_text(size = 14, colour = "black")) + 
                          geom_point(size = 3) +
                          scale_color_manual(values = c("#616161", "#9E9E9E", "#1565C0", "#64B5F6", "#C62828", "#E57373")) +
                          labs(title = "A", color = "Treatment", shape = "Sample Type")
  OsmiaCC_PCoA_bact
  
## Stacked community plot ----
  
# Generate colorblind friendly palette
  Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  
# Stretch palette (define more intermediate color options)
  okabe_ext <- usecol(Okabe_Ito, n = 29)
  colors <- sample(okabe_ext)
  
# Sort data by Family
  y1 <- tax_glom(rareps, taxrank = 'Family') # agglomerate taxa
  y2 <- transform_sample_counts(y1, function(x) x/sum(x))
  y3 <- psmelt(y2)
  y3$Family <- as.character(y3$Family)
  y3$Family[y3$Abundance < 0.01] <- "Family < 1% abund."
  y3$Family <- as.factor(y3$Family)
  head(y3)
  
# Reorder x-axis  
  y3$combo_treat <- factor(y3$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# New names for facet_grid
  type_names <- c('final provision' = "provisions with bee",
                  'provision w/o bee' = "provisions without bee")
  
# Plot treatment by Family
  OsmiaCC_fam_relabund_bact <- ggplot(data = y3, aes(x = combo_treat, y = Abundance, fill = Family)) + 
                                  geom_bar(stat = "identity", position = "fill") + 
                                  scale_fill_manual(values = colors) + 
                                  facet_grid(~ sample_type,
                                            scale = "free",
                                            space = "free",
                                            labeller = as_labeller(type_names)) +
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
                                  guides(fill = guide_legend(ncol = 2)) +
                                  ggtitle("A")
  OsmiaCC_fam_relabund_bact
  
# Plot Family for each sample
  ggplot(data = y3, aes(x = sampleID, y = Abundance, fill = Family)) + 
    geom_bar(stat = "identity", position = "fill") + 
    scale_fill_manual(values = colors) + 
    facet_grid(~ combo_treat, 
               scale = "free", 
               space = "free") +
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample ID") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Bacteria")
  
# Sort data by Genus
  y4 <- tax_glom(rareps, taxrank = 'Genus') # agglomerate taxa
  y5 <- transform_sample_counts(y4, function(x) x/sum(x))
  y6 <- psmelt(y5)
  y6$Genus <- as.character(y6$Genus)
  y6$Genus[y6$Abundance < 0.01] <- "Genera < 1% abund."
  write.csv(y6, file = "y6.csv")
  y6$Genus <- as.factor(y6$Genus)
  head(y6)
  
# Reorder x-axis  
  y6$combo_treat <- factor(y6$combo_treat,levels = c("CS", "CN", "AS", "AN", "WS", "WN"))
  
# Plot Genus by treatment
  OsmiaCC_gen_relabund_bact <- ggplot(data = y6, aes(x = combo_treat, y = Abundance, fill = Genus)) + 
                                    geom_bar(stat = "identity", position = "fill") + 
                                    scale_fill_manual(values = colors) + 
                                    facet_grid(~ sample_type, 
                                               scale = "free", 
                                               space = "free",
                                               labeller = as_labeller(type_names)) +
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
                                    guides(fill = guide_legend(ncol = 2)) +
                                    ggtitle("A")
  OsmiaCC_gen_relabund_bact
  
# Plot Genus for each sample
  ggplot(data = y6, aes(x = sampleID, y = Abundance, fill = Genus)) + 
    geom_bar(stat = "identity",
             position = "fill") + 
    scale_fill_manual(values = colors) + 
    facet_grid(~ sample_type,
               scale = "free",
               space = "free",
               labeller = as_labeller(type_names)) +
    theme(legend.position = "right") +
    ylab("Relative abundance") + 
    ylim(0, 1.0) +
    xlab("Sample ID") +
    theme_bw() + 
    theme(text = element_text(size = 16)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(legend.justification = "left", 
          legend.title = element_text(size = 16, colour = "black"), 
          legend.text = element_text(size = 14, colour = "black")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    guides(fill = guide_legend(ncol = 2)) +
    ggtitle("Bacteria")

## Differential abundance ----
# Resource: https://joey711.github.io/phyloseq-extensions/DESeq2.html

# Convert from a phyloseq to a deseq obj
  desq_obj <- phyloseq_to_deseq2(ps3, ~ combo_treat)
  
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
  
# WN vs AN
  
# Extract results from differential abundance table for initial vs final provision
  WN_AN <- results(desq_dds, contrast = c("combo_treat", "WN", "AN"))
  
# Order differential abundances by their padj value
  WN_AN <- WN_AN[order(WN_AN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_AN_p05 <- WN_AN[(WN_AN$padj < alpha & !is.na(WN_AN$padj)), ]
  
# Check to see if any padj is below alpha
  WN_AN_p05

# AN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  AN_CN <- results(desq_dds, contrast = c("combo_treat", "AN", "CN"))
  
# Order differential abundances by their padj value
  AN_CN <- AN_CN[order(AN_CN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  AN_CN_p05 <- AN_CN[(AN_CN$padj < alpha & !is.na(AN_CN$padj)), ]
  
# Check to see if any padj is below alpha
  AN_CN_p05

# WN vs CN
  
# Extract results from differential abundance table for initial vs final provision
  WN_CN <- results(desq_dds, contrast = c("combo_treat", "WN", "CN"))
  
# Order differential abundances by their padj value
  WN_CN <- WN_CN[order(WN_CN$padj, na.last = NA), ]
  
# Filter data to only include padj < alpha and remove NAs
  WN_CN_p05 <- WN_CN[(WN_CN$padj < alpha & !is.na(WN_CN$padj)), ]
  
# Check to see if any padj is below alpha
  WN_CN_p05
  
