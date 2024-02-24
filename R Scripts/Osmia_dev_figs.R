##### Project: Osmia Developmental Microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Create nice figures

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(patchwork) # Version 1.1.3

## Shannon index ----

# Create plot
  Osmia_dev_Shannon <- Osmia_dev_Shannon_bact + Osmia_dev_Shannon_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Shannon
  
# Save plot  
  ggsave("Osmia_dev_Shannon.png", plot = Osmia_dev_Shannon)
  
## Simpson index ----

# Create plot
  Osmia_dev_Simpson <- Osmia_dev_Simpson_bact + Osmia_dev_Simpson_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Simpson
  
# Save plot  
  ggsave("Osmia_dev_Simpson.png", plot = Osmia_dev_Simpson)
  
## Observed richness ----
  
# Create plot
  Osmia_dev_Observed <- Osmia_dev_Observed_bact + Osmia_dev_Observed_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Observed
  
# Save plot  
  ggsave("Osmia_dev_Observed.png", plot = Osmia_dev_Observed)
  
## Rarefaction curves ----

# All pollen and bee samples  
  
# Create plot
  Osmia_dev_rare <- Osmia_dev_rare_bact + Osmia_dev_rare_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_rare
  
# Save plot
  ggsave("Osmia_dev_rare.png", plot = Osmia_dev_rare)
  
# Only bee samples
  
# Create plot
  Osmia_dev_rare_bee <- Osmia_dev_rare_bact_bee + Osmia_dev_rare_fungi_bee + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_rare_bee
  
# Save plot
  ggsave("Osmia_dev_rare_bee.png", plot = Osmia_dev_rare_bee)
  
## PCoA plots ----
  
# All pollen and bee samples  
  
# Relative abundance data  
  
# Create plot
  Osmia_dev_PCoA_relabund <- Osmia_dev_PCoA_bact + Osmia_dev_PCoA_fungi + plot_layout(ncol = 2, nrow = 1)
  Osmia_dev_PCoA_relabund
  
# Save plot
  ggsave("Osmia_dev_PCoA_relabund.png", plot = Osmia_dev_PCoA_relabund)
  
# Rarefied data
  
# Create plot
  Osmia_dev_PCoA_rare <- Osmia_dev_PCoA_bact_rare + Osmia_dev_PCoA_fungi_rare + plot_layout(ncol = 2, nrow = 1)
  Osmia_dev_PCoA_rare
  
# Save plot
  ggsave("Osmia_dev_PCoA_rare.png", plot = Osmia_dev_PCoA_rare)

# Only bee samples
  
# Relative abundance data  
  
# Create plot
  Osmia_dev_PCoA_relabund_bee <- Osmia_dev_PCoA_bact_bee + Osmia_dev_PCoA_fungi_bee + plot_layout(ncol = 2, nrow = 1)
  Osmia_dev_PCoA_relabund_bee
  
# Save plot
  ggsave("Osmia_dev_PCoA_relabund_bee.png", plot = Osmia_dev_PCoA_relabund_bee)
  
# Rarefied data
  
# Create plot
  Osmia_dev_PCoA_rare_bee <- Osmia_dev_PCoA_bact_rare_bee + Osmia_dev_PCoA_fungi_rare_bee + plot_layout(ncol = 2, nrow = 1)
  Osmia_dev_PCoA_rare_bee
  
# Save plot
  ggsave("Osmia_dev_PCoA_rare_bee.png", plot = Osmia_dev_PCoA_rare_bee)
  
## Stacked Plot by Family ----
  
# Create plot
  Osmia_dev_fam_relabund <- Osmia_dev_fam_relabund_bact + Osmia_dev_fam_relabund_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_fam_relabund
  
# Save plot
  ggsave("Osmia_dev_fam_relabund.png", plot = Osmia_dev_fam_relabund)
  
## Stacked Plot by Genus ----
  
# Create plot
  Osmia_dev_gen_relabund <- Osmia_dev_gen_relabund_bact + Osmia_dev_gen_relabund_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_gen_relabund
  
# Save plot
  ggsave("Osmia_dev_gen_relabund.png", plot = Osmia_dev_gen_relabund)
  
## Stacked Plot by Top 15 Genera ----
  
# Create plot
  Osmia_dev_gen_top15 <- Osmia_dev_15gen_relabund_bact + Osmia_dev_15gen_relabund_fung + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_gen_top15
  
# Save plot
  ggsave("Osmia_dev_gen_top15.png", plot = Osmia_dev_gen_top15, width = 15, height = 10, units = "in")
