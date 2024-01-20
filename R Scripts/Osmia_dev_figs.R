##### Project: Osmia developmental microbiome

#### Ownership: Bailey Crowley & Robert N. Schaeffer

### Purpose: Create nice figures

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(patchwork) # Version 1.1.3

## Shannon index ----

# Create plot
  Osmia_dev_Shannon <- Shannon_bact + Shannon_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Shannon
  
# Save plot  
  ggsave("Osmia_dev_Shannon.png", plot = Osmia_dev_Shannon)
  
## Simpson index ----

# Create plot
  Osmia_dev_Simpson <- Simpson_bact + Simpson_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Simpson
  
# Save plot  
  ggsave("Osmia_dev_Simpson.png", plot = Osmia_dev_Simpson)
  
## Observed richness ----
  
# Create plot
  Osmia_dev_Observed <- Observed_bact + Observed_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_Observed
  
# Save plot  
  ggsave("Osmia_dev_Observed.png", plot = Osmia_dev_Observed)
  
## Rarefaction curves ----

# Create plot
  Osmia_dev_rare <- rare_bact + rare_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_rare
  
# Save plot
  ggsave("Osmia_dev_rare.png", plot = Osmia_dev_rare)
  
## PCoA plots ----
  
# Create plot
  Osmia_dev_PCoA <- PCoA_bact + PCoA_fungi + plot_layout(ncol = 2, nrow = 1)
  Osmia_dev_PCoA
  
# Save plot
  ggsave("Osmia_dev_PCoA.png", plot = Osmia_dev_PCoA)

## Relative abundance by Family ----
  
# Create plot
  Osmia_dev_fam_relabund <- fam_relabund_bact + fam_relabund_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_fam_relabund
  
# Save plot
  ggsave("Osmia_dev_fam_relabund.png", plot = Osmia_dev_fam_relabund)
  
## Relative abundance by Genus ----
  
# Create plot
  Osmia_dev_gen_relabund <- gen_relabund_bact + gen_relabund_fungi + plot_layout(ncol = 1, nrow = 2)
  Osmia_dev_gen_relabund
  
# Save plot
  ggsave("Osmia_dev_gen_relabund.png", plot = Osmia_dev_gen_relabund)
