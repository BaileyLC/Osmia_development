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
  Osmia.dev.Shannon <- Osmia.dev.Shannon.bact + Osmia.dev.Shannon.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.Shannon
  
# Save plot  
  ggsave("Osmia_dev_Shannon.png", plot = Osmia.dev.Shannon)
  
## Simpson index ----

# Create plot
  Osmia.dev.Simpson <- Osmia.dev.Simpson.bact + Osmia.dev.Simpson.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.Simpson
  
# Save plot  
  ggsave("Osmia_dev_Simpson.png", plot = Osmia.dev.Simpson)
  
## Observed richness ----
  
# Create plot
  Osmia.dev.Observed <- Osmia.dev.Observed.bact + Osmia.dev.Observed.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.Observed
  
# Save plot  
  ggsave("Osmia_dev_Observed.png", plot = Osmia.dev.Observed)
  
## Rarefaction curves ----

# All pollen and bee samples  
  
# Create plot
  Osmia.dev.rare <- Osmia.dev.rare.bact + Osmia.dev.rare.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.rare
  
# Save plot
  ggsave("Osmia_dev_rare.png", plot = Osmia.dev.rare)
  
# Only bee samples
  
# Create plot
  Osmia.dev.rare.bee <- Osmia.dev.rare.bact.bee + Osmia.dev.rare.fungi.bee + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.rare.bee
  
# Save plot
  ggsave("Osmia_dev_rare_bee.png", plot = Osmia.dev.rare.bee)
  
## PCoA plots ----
  
# All pollen and bee samples  
  
# Relative abundance data  
  
# Create plot
  Osmia.dev.PCoA.relabund <- Osmia.dev.PCoA.bact + Osmia.dev.PCoA.fungi + plot_layout(ncol = 2, nrow = 1)
  Osmia.dev.PCoA.relabund
  
# Save plot
  ggsave("Osmia_dev_PCoA_relabund.png", plot = Osmia.dev.PCoA.relabund)
  
# Rarefied data
  
# Create plot
  Osmia.dev.PCoA.rare <- Osmia.dev.PCoA.bact.rare + Osmia.dev.PCoA.fungi.rare + plot_layout(ncol = 2, nrow = 1)
  Osmia.dev.PCoA.rare
  
# Save plot
  ggsave("Osmia_dev_PCoA_rare.png", plot = Osmia.dev.PCoA.rare)

# Only bee samples
  
# Relative abundance data  
  
# Create plot
  Osmia.dev.PCoA.relabund.bee <- Osmia.dev.PCoA.bact.bee + Osmia.dev.PCoA.fungi.bee + plot_layout(ncol = 2, nrow = 1)
  Osmia.dev.PCoA.relabund.bee
  
# Save plot
  ggsave("Osmia_dev_PCoA_relabund_bee.png", plot = Osmia.dev.PCoA.relabund.bee)
  
# Rarefied data
  
# Create plot
  Osmia.dev.PCoA.rare.bee <- Osmia.dev.PCoA.bact.rare.bee + Osmia.dev.PCoA.fungi.rare.bee + plot_layout(ncol = 2, nrow = 1)
  Osmia.dev.PCoA.rare.bee
  
# Save plot
  ggsave("Osmia_dev_PCoA_rare_bee.png", plot = Osmia.dev.PCoA.rare.bee)
  
## Stacked Plot by Family ----
  
# Create plot
  Osmia.dev.fam.relabund <- Osmia.dev.fam.relabund.bact + Osmia.dev.fam.relabund.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.fam.relabund
  
# Save plot
  ggsave("Osmia_dev_fam_relabund.png", plot = Osmia.dev.fam.relabund)
  
## Stacked Plot by Genus ----
  
# Create plot
  Osmia.dev.gen.relabund <- Osmia.dev.gen.relabund.bact + Osmia.dev.gen.relabund.fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.gen.relabund
  
# Save plot
  ggsave("Osmia_dev_gen_relabund.png", plot =  Osmia.dev.gen.relabund)
  
## Stacked Plot by Top 15 Genera ----
  
# Create plot
  Osmia.dev.gen.top15 <- Osmia.dev.15gen.relabund.bact + Osmia.dev.15gen.relabund_fung + plot_layout(ncol = 1, nrow = 2)
  Osmia.dev.gen.top15
  
# Save plot
  ggsave("Osmia_dev_gen_top15.png", plot = Osmia.dev.gen.top15, width = 15, height = 10, units = "in")
