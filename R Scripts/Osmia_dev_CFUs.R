##### Project: Osmia Developmental Microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Analysis of culture-dependent data

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")
  
# Load necessary packages
  library(ggplot2) # Version 3.4.3
  library(scales) # Version 1.3.0
  library(patchwork) # Version 1.1.3
  library(magrittr) # Version 2.0.3
  library(tidyverse) # Version 1.2.0
  library(stats) # Version 4.3.1
  library(ggpubr) # Version 0.6.0
  
# Import data
  CFUs <- read.csv("sampling_CFUs - Data.csv")

## Clean data ----
  
# Subset data to include pollen provision only  
  provision_CFUs <- subset(CFUs, development_stage %in% c("fresh", "aged"))

# Subset df to include just development stage, provision weight, and CFU counts
  provision_CFUs <- provision_CFUs %>%
    select(development_stage, weight_g, CFU_YM_aer, CFU_R2A_aer)
  
# Replace "TNTC" colony counts with 300
  provision_CFUs$CFU_YM_aer[provision_CFUs$CFU_YM_aer == "TNTC"] <- 300
  
# Replace "0" provision weights with 0.09 (scale was limited to measuring in grams)
  provision_CFUs$weight_g[provision_CFUs$weight_g == "0"] <- 0.09
  
# Format CFU counts as numeric
  provision_CFUs$CFU_YM_aer <- as.numeric(provision_CFUs$CFU_YM_aer)
  provision_CFUs$CFU_R2A_aer <- as.numeric(provision_CFUs$CFU_R2A_aer)

## Analyses ----

# Divide CFUs by provision weight, then times by 10 (dilution factor = 1/10)
  provision_CFUs$YM_new <- (provision_CFUs$CFU_YM_aer / provision_CFUs$weight_g)*10
  provision_CFUs$R2A_new <- (provision_CFUs$CFU_R2A_aer / provision_CFUs$weight_g)*10
  
# Remove rows with NAs and 0 CFUs
  provision_CFUs <- na.omit(provision_CFUs)
  provision_CFUs <- provision_CFUs[apply(provision_CFUs != 0, 1, all), ]
  
# Summary statistics
  provision_CFUs %>%
    group_by(development_stage) %>%
    summarise(YM_mean = mean(YM_new),
              YM_sd = sd(YM_new),
              R2A_mean = mean(R2A_new),
              R2A_sd = sd(R2A_new))

# Shapiro-Wilk tests for normality
  stats::shapiro.test(provision_CFUs$YM_new)
  stats::shapiro.test(provision_CFUs$R2A_new)

# Fungi

# Group YM CFU data for initial provisions
  YM_initial <- provision_CFUs %>%
    filter(development_stage == "fresh") %>%
    pull(YM_new)

# Group YM CFU data for final provisions  
  YM_final <- provision_CFUs %>%
    filter(development_stage == "aged") %>%
    pull(YM_new)
  
# Perform Mann Whitney rank sums test on fungi
  stats::wilcox.test(YM_initial, YM_final, exact = FALSE)
  
# Bacteria
  
# Group R2A CFU data for initial provisions
  R2A_initial <- provision_CFUs %>%
    filter(development_stage == "fresh") %>%
    pull(R2A_new)
  
# Group R2A CFU data for final provisions  
  R2A_final <- provision_CFUs %>%
    filter(development_stage == "aged") %>%
    pull(R2A_new)
  
# Perform Mann Whitney rank sums test on bacteria
  stats::wilcox.test(R2A_initial, R2A_final, exact = FALSE)

## Plotting ----
  
# Reorder the x-axis
  provision_CFUs$development_stage <- factor(provision_CFUs$development_stage, levels = c("fresh", "aged"))
  
# Plot bacteria
  R2A_plot <- ggplot(provision_CFUs, aes(x = development_stage, y = R2A_new, color = development_stage)) + 
                geom_boxplot(outlier.shape = NA,
                             width = 0.5,
                             position = position_dodge(width = 0.1)) + 
                geom_jitter(size = 1, 
                            alpha = 0.9) +
                theme_bw() +
                theme(legend.position = "none") + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                theme(text = element_text(size = 12)) +
                scale_color_manual(values = c("#D32F2F","#E57373")) +
                scale_y_continuous(trans = log10_trans(),
                                   breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)),
                                   limits = c(1, 1e6)) +
                ggtitle("A") +
                ylab("CFUs") +
                xlab("Provision") + 
                stat_compare_means(method = "wilcox.test",
                                   label.x = 1.25,
                                   label.y.npc = 'top')
  R2A_plot
  
# Plot fungi
  YM_plot <- ggplot(provision_CFUs, aes(x = development_stage, y = YM_new, color = development_stage)) + 
                geom_boxplot(outlier.shape = NA,
                             width = 0.5,
                             position = position_dodge(width = 0.1)) + 
                geom_jitter(size = 1, 
                            alpha = 0.9) +
                theme_bw() +
                theme(legend.position = "none") + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                theme(text = element_text(size = 12)) +
                scale_color_manual(values = c("#1565C0","#64B5F6")) +
                scale_y_continuous(trans = log10_trans(),
                                   breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)),
                                   limits = c(1, 1e6)) +
                ggtitle("B") +
                ylab("CFUs") +
                xlab("Provision") +
                stat_compare_means(method = "wilcox.test",
                                   label.x = 1.25,
                                   label.y.npc = 'top')
  YM_plot
    
# Arrange plots
  culture_dep <- R2A_plot + YM_plot + plot_layout(ncol = 2, nrow = 1)
  culture_dep
  