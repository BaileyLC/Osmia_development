##### Project: Osmia developmental microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Quality control of 16S rRNA data using the DADA2 pipeline
### Resource: https://benjjneb.github.io/dada2/tutorial_1_8.html

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(dada2) # Version: 1.28.0

## Prepare samples for quality control ----

# Specify the path where your FASTQ files are located
  path <- "seqs/Osmia_dev16S" # where to find the directory containing the unzipped fastq files
  list.files(path) # call files to ensure you are in the correct location

# Sort files to ensure forward and reverse reads are in the same order
  fnFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
  sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) 

## Check the quality of your clean reads ----  

# Plot average quality scores for forward reads
  plotQualityProfile(fnFs[1:2]) # Use Q30 to determine where to cut

# Plot average quality scores for reverse reads
  plotQualityProfile(fnRs[1:2]) # Use Q30 to determine where to cut

##  Prepare a path for your filtered reads ----  

# Create a new file path to store filtered and trimmed reads
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered subdirectory

# Rename filtered files and place them in the filtered subdirectory
  filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names

## Clean your reads ----   

# Quality filtering and trimming
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,  # specify object that contains your unfiltered reads, followed by the object that contains your unfiltered reads
                       truncLen = c(280, 200), # truncates reads a specified base (forward, reverse) 
                       maxN = 0, # set the maximum number of ambiguous bps allowed (must be zero for DADA2)
                       maxEE = c(2, 2), # set the maximum number of estimated errors allowed for an individual read
                       truncQ = 2, # truncates reads at the first instance of a Q score less than or equal to the specified value
                       rm.phix = TRUE, # removes reads identified as belonging to the phiX phage
                       compress = TRUE, # zips FASTQ files
                       multithread = TRUE) # allows FASTQ files to be processed in parallel

  head(out)

# Estimate the error model for DADA2 algorithm using reads
  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)  

# Visualize the estimated error rates
  plotErrors(errF, nominalQ = TRUE)

# Apply the core sample inference algorithmm to the filtered and trimmed sequence data
  dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
  dadaRs <- dada(filtRs, err = errR, multithread = TRUE)  

# Inspect the returned object
  dadaFs
  dadaRs

## Merge reads ----  

# Merge denoised reads
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
  head(mergers[[1]])

# Organize ribosomal sequence variants (RSVs) in to a sequence table
  seqtab <- makeSequenceTable(mergers)
  dim(seqtab)

# View the length of all total RSVs
  table(nchar(getSequences(seqtab))) # first row: size; second row: frequency

## Check for and remove chimeras ----  

# Perform de novo chimera sequence detection and removal
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
  dim(seqtab.nochim)

# Calculate the proportion of non-chimeric RSVs  
  sum(seqtab.nochim)/sum(seqtab)

## Save metadata ----  

# Check the number of reads that made it through each step in the pipeline
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
  rownames(track) <- sample.names
  head(track)

## Assign taxonomy ----  

# Assign taxonomy using a reference database  
  taxa <- assignTaxonomy(seqtab.nochim, "refseqs/rdp_train_set_18.fa.gz", multithread = TRUE) # What is the prob that train set and your data match
  taxa <- addSpecies(taxa, "refseqs/rdp_species_assignment_18.fa.gz")

# Removing sequence row names for display only
  taxa.print <- taxa
  rownames(taxa.print) <- NULL
  head(taxa.print)

## Save output ---- 

# Save files so you don't have to work through the computationally heavy work again
  saveRDS(seqtab.nochim, file = "Osmia_dev_seqs16S.rds")
  saveRDS(taxa, file = "Osmia_dev_taxa16S.rds")
