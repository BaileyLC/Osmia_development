##### Project: Osmia developmental microbiome

#### Owners: Bailey Crowley & Robert N. Schaeffer

### Purpose: Quality control of ITS rRNA data using the DADA2 pipeline
###: Resource: https://benjjneb.github.io/dada2/ITS_workflow.html

## Prepare work space ----

# Set working directory
  setwd("~/Downloads")

# Load necessary packages
  library(dada2) # Version: 1.28.0
  library(ShortRead) # Version 1.58.0
  library(Biostrings) # Version 2.68.1

## Prepare samples for quality control ----

# Specify the path where your FASTQ files are located
  path <- "Osmia_devITS"
  list.files(path)
  
# Sort files to ensure forward and reverse reads are in the same order
  fnFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
  
# Identify primers ----
  
# Define primers
  FWD <- "GTGAATCATCGAATCTTTGAA"
  REV <- "TCCTCCGCTTATTGATATGC"
  
# Create a function that will list all orientation of the primers
  allOrients <- function(primer) {
    require(Biostrings)
    dna <- DNAString(primer) # convert to DNAString object rather than character vectors
    orients <- c(Forward = dna,
                 Complement = Biostrings::complement(dna), 
                 Reverse = Biostrings::reverse(dna),
                 RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString)) # convert back to character vectors
  }
  
# Use function to list all orientations of the primers
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  
# Create a subdirectory to place sequences with ambiguous bases
  fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
  fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
  
# Pre-filter sequences to remove those with Ns
  filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Count the number of times primers (in any orientation) appear in reads  
  primerHits <- function(primer, fn) {
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  
# Display the number of times primers appear in reads
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
  
## Remove primers ----
 
# Specify the path where cutadapt is located 
  cutadapt <- "/Users/baileycrowley/miniconda3/bin/cutadapt"
  system2(cutadapt, args = "--version") # Run shell commands from R
  
# Create a new file path to store cutadapt-ed files
  path.cut <- file.path(path, "cutadapt")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  fnFs.cut <- file.path(path.cut, basename(fnFs))
  fnRs.cut <- file.path(path.cut, basename(fnRs))

# Define reverse complements of the FWD and REV primers
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)

# Flag the FWD and the reverse-complement of REV off the forward reads
  R1.flags <- paste("-g", FWD, "-a", REV.RC)
  
# Flag the REV and the reverse-complement of FWD off the reverse reads
  R2.flags <- paste("-G", REV, "-A", FWD.RC)
  
# Run cutadapt
  for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
  }

# Check to see if it worked by determining whether there are any primers left in the cutadapt-ed samples
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
  
# Sort files to ensure forward and reverse reads are in the same order
  cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))
  cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
  
# Extract sample names
  sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 1)
  head(sample.names)
  
## Check the quality of your clean reads ----  
  
# Plot average quality scores for forward reads
  plotQualityProfile(cutFs[1:2]) # Use Q30 to determine where to cut
  
# Plot average quality scores for reverse reads
  plotQualityProfile(cutRs[1:2]) # Use Q30 to determine where to cut
  
##  Prepare a path for your filtered reads ----  
  
# Rename filtered files and place them in the filtered subdirectory
  filtFs <- file.path(path.cut, "filtered", basename(cutFs))
  filtRs <- file.path(path.cut, "filtered", basename(cutRs))
  
  ## Clean your reads ----   
  
# Quality filtering and trimming
  out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                       maxN = 0, 
                       maxEE = c(2, 2), 
                       truncQ = 2,
                       minLen = 50, 
                       rm.phix = TRUE, 
                       compress = TRUE, 
                       multithread = TRUE)  # on windows, set multithread = FALSE
  
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
  unite.ref <- "sh_general_release_dynamic_s_25.07.2023.fasta"
  taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
  
# Removing sequence rownames for display only
  taxa.print <- taxa
  rownames(taxa.print) <- NULL
  head(taxa.print)
  
## Save output ----
  
# Save files so you don't have to work through the computationally heavy work again
  saveRDS(seqtab.nochim, file = "Osmia_dev_seqsITS.rds")
  saveRDS(taxa, file = "Osmia_dev_taxaITS.rds")
