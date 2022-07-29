# Set working directory
setwd("/media/edwin/Transcend/datasets/dada2/meta_analysis")
# Load packages -----------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")
library(ShortRead)
library(dada2)
library(tidyr)
library(Hmisc)
library(plotly)
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dunn.test)
library(vegan)
library(randomForest)
library(dplyr)
library(DESeq2)
library(ALDEx2)
library(metagenomeSeq)
library(microbiome)

# Load in data ------------------------------------------------------------

data_path <- ("/media/edwin/Transcend/datasets/dada2/meta_analysis") 
list.files(data_path)

# identify primers --------------------------------------------------------

FWD <- "GTGCCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Create new file path to Place filtered files in filtered/ subdirectory
filt_path <- file.path(data_path, "filtered_N") 
if(!file_test("-d", filt_path)) dir.create(filt_path)

fnFs <- sort(list.files(data_path, pattern="_1.fastq.gz"))
fnRs <- sort(list.files(data_path, pattern="_2.fastq.gz"))

sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)

fnFs <- file.path(data_path, fnFs)
fnRs <- file.path(data_path, fnRs)
fnFs

# Rename filtered files
fnFs.filtN <- file.path(filt_path, paste0(sampleNames, "_filtN_1.fastq.gz"))
fnRs.filtN <- file.path(filt_path, paste0(sampleNames, "_filtN_2.fastq.gz"))


# Filter and trim the forward and reverse reads:
filterAndTrim(fnFs, fnFs.filtN , fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#check for primer hits
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

cutadapt <- "/home/edwin/anaconda3/envs/qiime2-2021.2/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(data_path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 12, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
# Filter and Trim ---------------------------------------------------------
new_path <- "/home/edwin/dada2/PRJNA313391/cutadapt"

list.files(new_path)

# Sort to ensure forward/reverse reads are in same order
New_fnFs <- sort(list.files(data_path, pattern="fastq"))
head(New_fnFs)
str(data_path)
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(unlist(strsplit(New_fnFs, "[.]")), `[`, 1)
write.csv(sampleNames, "trying.csv")


sampleNames <- read.table("meta_analysis.txt")
y <- rownames(sampleNames)
z <- as.list(as.data.frame(sampleNames$V1))
x <- list()



# Specify the full path to the fnFs and fnRs
New_fnFs <- file.path(data_path, New_fnFs)
New_fnFs[1:3]

# Plot to check the quality of reads
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(New_fnFs) <= 20) {
  fwd_qual_plots<-plotQualityProfile(New_fnFs)
} else {
  rand_samples <- sample(size = 20, 1:length(New_fnFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(New_fnFs[rand_samples])
}

# Now that we have the filenames, we're going to inspect the quality profile of our data:
fwd_qual_plots <- plotQualityProfile(New_fnFs[1:4])
# write plots to disk
saveRDS(fwd_qual_plots, "fwd_qual_plots_1.rds")

ggsave(plot = fwd_qual_plots, filename = "fwd_qual_plots_4.png", 
       width = 10, height = 10, dpi = "retina")

# We define the filenames for the filtered fastq.gz files:
# Create new file path to Place filtered files in filtered/ subdirectory
filt_path <- file.path(data_path, "filtered") 
if(!file_test("-d", filt_path)) dir.create(filt_path)

# Rename filtered files
filtFs <- file.path(filt_path, paste0(z$`sampleNames$V1`, "_filt.fastq.gz"))


# Filter and trim the forward and reverse reads:

out <- filterAndTrim(New_fnFs, filtFs, maxN = 0, trimLeft = 20, minLen = 100, truncLen = 210,
                     maxEE=2, truncQ=2, rm.phix=TRUE, 
                     compress=TRUE,multithread=TRUE) # On Windows set multithread=FALSE
# look at how many reads were kept
tail(out)

#plot quality of trimdata
filt_plot <- plotQualityProfile(filtFs[1:4])
ggsave(plot = filt_plot, filename = "filt_plots_2.png", 
       width = 10, height = 10, dpi = "retina")

# summary of samples in filt_out by percentage
out %>% 
  data.frame() %>% 
  mutate(Samples= rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))

# Infer sequence variants and Dereplication -------------------------------------------------

derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- z$`sampleNames$V1`

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF)

# Infer sample composition 
dadaFS<-dada(derepFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)

dadaFS[[1]]

# Construct sequence table and remove chimeras ----------------------------

# Merge forward/reverse reads
#merge.reads <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#head(merge.reads[[1]])

# Construct sequence table
seqtabAll <- makeSequenceTable(dadaFS)
table(nchar(getSequences(seqtabAll)))

# Remmove chimeras
seqtabNoC <- removeBimeraDenovo(seqtabAll, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtabNoC)
sum(seqtabNoC)/sum(seqtabAll)
str(seqtabNoC)
# Track Reads through the DADA2 Pipeline ----------------------------------

getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dadaFS, getN), rowSums(seqtabNoC))

colnames(track.nbr.reads) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track.nbr.reads) <- z$`sampleNames$V1`
track.nbr.reads

write.csv(track.nbr.reads, file = "meta_analysis_prepro.csv")

uniquesToFasta(seqtabAll, "meta_analysis_prepro.csv.fasta")



# Assign taxonomy ---------------------------------------------------------

taxTab <- assignTaxonomy(seqtabNoC, "/media/edwin/Transcend/projectdata/RefSeq-RDP16S_v3_May2018.fa.gz", multithread=TRUE)

taxTab <- assignTaxonomy(seqtabNoC, "/home/edwin/gg_13_5.fasta.gz", multithread=TRUE)

unname(head(taxTab))

# Optionally assign species
taxTab <- addSpecies(taxTab, "/home/qcifrimpong/Microbiomics/dada2_tutorial_dog/silva_species_assignment_v132.fa.gz")

taxa.print <- taxTab # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
tail(taxa.print)

# Save the ASVs table in your working directory ---------------------------

write.csv(taxTab, file="ASVs_tax_meta_analysis.csv")
saveRDS(taxTab, "ASVs_tax_meta_analysis.rds")

asv_headers <- vector(dim(seqtabNoC)[2], mode="character")
count.asv.tab <- t(seqtabNoC)
#row.names(count.asv.tab) <- sub(">", "", asv_headers)

write.csv(count.asv.tab, file="ASVs_counts_meta_analysis.csv")
saveRDS(count.asv.tab, file="ASVs_counts_meta_analysis.rds")

# 10. Alignment -----------------------------------------------------------

library(DECIPHER)
seqs <- getSequences(seqtabNoC)
names(seqs) <- ASV.IDs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

# 11. Construct Phylogenetic Tree -----------------------------------------

#install.packages('phangorn')
library(phangorn)
library(Biostrings)
library(ggplot2)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA") 
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

# change the negative edges length to 0
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

saveRDS(fitGTR, "phangorn.tree.RDS")

# Load custom functions, set working directories, import data --------------------------
source("C:/Users/Numb3rs/Documents/Wk2/microbiome_custom_functions_tutorial.R")

ASVs <- readRDS("~/Downloads/ASVs_counts.RDS")#Default options with primer trimming & new refseq_rdp DB
str(ASVs)
dim(ASVs)

#Read in taxonomy table
taxa <- readRDS("~/Downloads/ASVs_taxonomy.RDS")
dim(taxa)

#Read in ASV counts#
identical(rownames(ASVs),rownames(taxTab))#TRUE
#Assign user-friendly ASV IDs to replace sequences
head(rownames(ASVs))
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
head(ASVs)

#Named vector:
rownames(ASVs) <- ASV.IDs
rownames(taxa) <- ASV.IDs
colnames(ASVs) <- z$`sampleNames$V1`
names(seqs) <- ASV.IDs
head(seqs)
seq_lens <- nchar(seqs)
seq_lens
plot(density(seq_lens))

head(ASVs)
head(taxa)
write.csv(ASVs, "ASV_count_meta.csv")
uniquesToFasta(seqtabNoC, fout="meta.fasta", ids = ASV.IDs)


#Merge ASV table and taxonomic table as a phyloseq object
phy <- phyloseq(otu_table(ASVs,taxa_are_rows=TRUE),tax_table(taxa))

identical(taxa_names(phy),rownames(ASVs))#TRUE
taxa_names(phy) <- names(seqs)
str(phy)
tax_df <- data.frame(tax_table(phy))
tax_df$seq_len <- seq_lens
head(tax_df)
dim(tax_df)
length(which(!is.na(tax_df[,"Species"])))# if GG ~ 32% assigned | If RefSeq-RDP ~ 70% assigned at species-level (3048/4415)

# Import metadata --------------------------------------------------
meta <-  read.delim("~/Downloads/practice.dataset1.metadata.tsv", header =TRUE, row.names=1)
head(meta)
##       Dog Treatment
## Dog1    B         2
## Dog2    G         3
## Dog3    K         3
## Dog8    B         4
## Dog9    G         0
## Dog10   K         4

rownames(meta)

##  [1] "Dog1"  "Dog2"  "Dog3"  "Dog8"  "Dog9"  "Dog10" "Dog15" "Dog16"
##  [9] "Dog17" "Dog22" "Dog23" "Dog24" "Dog29" "Dog30" "Dog31"

head(sample_names(phy))

## [1] "Dog10" "Dog15" "Dog16" "Dog17" "Dog1"  "Dog22"

length(sample_names(phy))#15

## [1] 15

length(rownames(meta))#15 (check if same number of samples in .biom file and metadatafile)
head(meta)
## [1] 15

length(intersect(rownames(meta),sample_names(phy)))#15 (check that the sample names match in all cases)

## [1] 15

#--------------------
# Assign the metadata to the phyloseq object 'phy' (phyloseq will put these in the right order)
sample_data(phy) <- meta 
identical(sample_names(phy),colnames(otu_table(phy)))#see
#--------------------
str(sample_data(phy)) # Need to change treatment column to factor variable

# 'data.frame':	15 obs. of  2 variables:
#   Formal class 'sample_data' [package "phyloseq"] with 4 slots
# ..@ .Data    :List of 2
# .. ..$ : Factor w/ 3 levels "B","G","K": 1 3 1 2 3 2 1 2 3 1 ...
# .. ..$ : int  2 4 1 4 0 3 3 1 2 0 ...
# ..@ names    : chr  "Dog" "Treatment"
# ..@ row.names: chr  "Dog1" "Dog10" "Dog15" "Dog16" ...
# ..@ .S3Class : chr "data.frame"

sample_data(phy)[,"Conditions"] <- as.factor(unlist(sample_data(phy)[,"Conditions"]))
saveRDS(phy, "meta_analysis_phyloseq.RDS") 

phy <- readRDS("/media/edwin/Transcend/datasets/dada2/meta_analysis/meta_analysis_phyloseq.RDS")
TREE = read_tree("/media/edwin/Transcend/datasets/dada2/meta_analysis/picrust/meta.tre")
phy_tree(phy) <- TREE
setwd("/media/edwin/Transcend/datasets/dada2/meta_analysis/")
# Check the number of reads per sample
reads <- sample_sums(phy) 
head(reads)
table(tax_table(phy)[, "Kingdom"], exclude = NULL)
ps <- subset_taxa(phy, !is.na(Species) & !Phylum %in% c("", "uncharacterized"))

filterKingdom = c("Archaea", "NA")
phy = subset_taxa(ps, !Kingdom %in% filterKingdom)


# Check the number of reads per sample
reads <- sample_sums(phy) 
reads

# Standardize sample read count to compare between samples
## Set (arbitrary) cutoff for number of acceptable reads/sample
length(which(reads<1000))
## Find median sample read count
total = median(sample_sums(phy))
## Function to standardize to median sample read count
stand_f = function(x, t=total) round(t * (x / sum(x)))

# Apply to phyloseq object
phy_std = transform_sample_counts(phy, stand_f)
ntaxa(phy_std)
tax_table(phy_std)[,"Species"] <- gsub("\\s*\\([^\\)]+\\)","",tax_table(phy_std)[,"Species"])
tax_table(phy_std)[,"Species"] <-gsub("'","",tax_table(phy_std)[,"Species"])

ps.taxa <- tax_glom(phy_std, taxrank = 'Species', NArm = TRUE)
ntaxa(ps.taxa)
M.f = filter_taxa(ps.taxa,function(x) sum(x > 100) > (0.2*length(x)) | sum(x) > 0.01*total, TRUE)
ntaxa(M.f)
tax_tab <- data.frame(tax_table(M.f))
dim(tax_tab)
x <- tax_tab$Species
y <- x[duplicated(x)]
length(y)

M.f <- subset_taxa(M.f, !Species %in% filterspecies)
taxa_names(M.f) <- tax_table(M.f)[,"Species"]

saveRDS(M.f, file = "filtered_meta_analysis_phyloseq.RDS")