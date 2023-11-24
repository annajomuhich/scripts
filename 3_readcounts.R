################# .bam to Read Counts (RSubread)
################# Anna Jo Muhich
################# August2023


###Install packages
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("Rsubread")
#install.packages("tidyverse")
#BiocManager::install("rtracklayer")



###Load packages
library(tidyverse)
library(rtracklayer)
library(Rsubread)
library(tools)


######## Prepare references (if needed)
# ###Convert gff to gtf
# setwd("~/UCDavis/Klieb_Lab/Projects/Cucurbit/Cuc_RNAseq_Pilot/reference")
# gff <- import.gff("ChineseLong_v3.gff3")
# #create a column "gene_id" that contains the gene name for every entry
# gff$gene_id <- ifelse(is.na(gff$ID),gff$Parent,gff$ID)
# #export
# export(gff,"ChineseLong_v3.gtf",format="gtf")
# 
# ###Convert gff to gtf
# setwd("~/UCDavis/Klieb_Lab/Projects/Cucurbit/Cuc_RNAseq_Pilot/reference")
# gff <- import.gff("Gy14_gene_gff_v2")
# #create a column "gene_id" that contains the gene name for every entry
# gff$gene_id <- ifelse(is.na(gff$ID),gff$Parent,gff$ID)
# #export
# export(gff,"Gy14_gene_gff_v2.gtf",format="gtf")
# 
# ###Convert gff to gtf
# setwd("~/UCDavis/Klieb_Lab/Projects/Cucurbit/Cuc_RNAseq_Pilot/reference")
# gff <- import.gff("Cpepo_v4.1.gff3")
# #create a column "gene_id" that contains the gene name for every entry
# gff$gene_id <- ifelse(is.na(gff$ID),gff$Parent,gff$ID)
# #export
# export(gff,"Cpepo_v4.1.gtf",format="gtf")

###Repeat for Bcin
# gff <- import.gff("Botrytis_cinerea.ASM83294v1.56.gff3")
# #create a column "gene_id" that contains the gene name for every entry
# gff$gene_id <- ifelse(is.na(gff$ID),gff$Parent,gff$ID)
# #export
# export(gff,"Botrytis_cinerea.ASM83294v1.5.gtf",format="gtf")


######## Loop for getting read counts that discerns between Host and Bcin bams.

#start where your sample directories containing .bams are
setwd("./bams")

# get the list of subdirectories in the main directory
subdirs <- setdiff(list.dirs(), ".")  # remove "." from subdirs
# Create an empty list to hold count matrices for each file
count_list <- list()

# loop through the subdirectories
for (subdir in subdirs) {
  # change the working directory to the subdirectory
  setwd(subdir)
  # get the list of files in the subdirectory
  file_list <- list.files()
  for (file_name in file_list) {
    # Check if file is a BAM file mapped to the plant genome (these files end in "Host.bam")
    if (grepl("Host.bam$", file_name)) {
      
      # Run featureCounts on the file with Host genes
      count_matrix <- featureCounts(file = file_name,
                                    annot.ext = "~/fastq2readcounts/reference/Vunguiculata_540_v1.2.gene_exons.gtf",
                                    isGTFAnnotationFile = TRUE,
                                    isPairedEnd = TRUE)
      # Assign the count matrix to an object with a name based on the file name
      count_list[[file_name]] <- count_matrix
    }
    # Check if the file is a BAM file mapped to the Bcin genome (these files end in "Bcin.bam")
    if (grepl("Bcin.bam$", file_name)) {
      
      # Run featureCounts on the file with Bcin genes
      count_matrix <- featureCounts(file = file_name,
                                    annot.ext = "~/fastq2readcounts/reference/Botrytis_cinerea.ASM83294v1.5.gtf",
                                    isGTFAnnotationFile = TRUE,
                                    isPairedEnd = TRUE)
      # Assign the count matrix to an object with a name based on the file name
      count_list[[file_name]] <- count_matrix
    }
  }
  setwd("..")
}

#make an empty vector for sample list output
Host_spl_list <- list()
Bcin_spl_list <- list()

######### Loop to get RSubread output into dataframes and reformatted
for (i in 1:length(count_list)) {
  df <- count_list[[i]]                                       #pulls out a df from the big object
  counts <- df$counts                                         #get the counts out of the df
  counts <- as.data.frame(counts)                             #make the counts a df
                                                              #convert transcripts from rownames to column
  counts <- rownames_to_column(counts,
                               var = "transcript")
  spl_id <- colnames(counts[2]) %>%                           #creates a simple sample ID w/ the plate ID and whether its Host/Bcin
    str_replace("(.*)_(.*).bam", "\\1_\\2")                   #MAKE CHANGES HERE for your naming convention if you want
  colnames(counts) <- c("transcript", paste0(spl_id))         #simplify the column name with the sample ID
  new_object_name <- paste0("counts_", spl_id)                #make a new object name with the sample ID
  assign(new_object_name, counts)                             #put the finished counts df in an object with its new name
  if (grepl("Host$", new_object_name)) {                       #check if sample name contains Host or Bcin.
    Host_spl_list <- append(Host_spl_list, list(counts))        #add the newly generated df to its corresponding list.
  }
  if (grepl("Bcin$", new_object_name)) {
    Bcin_spl_list <- append(Bcin_spl_list, list(counts))
  }           
}

######## MERGE the counts of diff samples.

#Put the first df in the list into the df
counts_Host_combined <- Host_spl_list[[1]]
counts_Bcin_combined <- Bcin_spl_list[[1]]

#Loop to join the remaining Host dfs together
for (i in 2:length(Host_spl_list)) {
  counts_Host_combined <- full_join(counts_Host_combined,
                                   Host_spl_list[[i]],
                                   by = "transcript")
}
#Loop to join the remaining Bcin dfs together
for (i in 2:length(Bcin_spl_list)) {
  counts_Bcin_combined <- full_join(counts_Bcin_combined,
                                   Bcin_spl_list[[i]],
                                   by = "transcript")
}

#clean up column names for samples
colnames(counts_Host_combined) <- gsub("_Host", "", colnames(counts_Host_combined))
colnames(counts_Bcin_combined) <- gsub("_Bcin", "", colnames(counts_Bcin_combined))

#save em!
write.csv(x = counts_Host_combined, file = "~/fastq2readcounts/readcounts/Host_readcounts.csv", row.names = F)
write.csv(x = counts_Bcin_combined, file = "~/fastq2readcounts/readcounts/Bcin_readcounts.csv", row.names = F)
