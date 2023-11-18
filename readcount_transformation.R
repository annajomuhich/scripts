########### Readcount transformation and TMM normalization
########### September 2023
########### Anna Jo Muhich

library(tidyverse)
library(ggplot2)
library(edgeR)

#read in data
host <- read.csv("ZE_1/readcounts/Host_readcounts.csv")
bcin <- read.csv("ZE_1/readcounts/Bcin_readcounts.csv")

#clean up transcript names
bcin$transcript <- sub("^transcript:", "", bcin$transcript)

#remove noncoding RNA for botrytis
bcin <- bcin %>% filter(!grepl("^E", transcript))

##### Collapse exon reads into total reads for each gene
#make gene column
# #Use this one for Csat
# host <- host %>%
# 	mutate(gene = sub("\\.[0-9].exon[0-9]+$", "", transcript)) %>%
# 	select(transcript, gene, everything())
#Use this one for Cpepo (e.g. Cp4.1_LG00_1000g00010.1:exon:001 becomes Cp4.1_LG00_1000g00010)
host <- host %>%
	mutate(gene = sub("\\.[0-9]:exon:[0-9]+$", "", transcript)) %>%
	select(transcript, gene, everything())
#This one for Botrytis
bcin <- bcin %>%
	mutate(gene = sub("\\.[0-9]+$", "", transcript)) %>%
	select(transcript, gene, everything())

#add together counts within each gene for each sample
#These for Csat:
# host <- host %>%
# 	group_by(gene) %>%
# 	summarise(across(starts_with("B"), sum)) #starts_with("B") because my sample IDs are e.g. "B27_1"
# bcin <- bcin %>%
# 	group_by(gene) %>%
# 	summarise(across(starts_with("B"), sum))
#These for Cpep:
host <- host %>%
	group_by(gene) %>%
	summarise(across(starts_with("Z"), sum)) #starts_with("Z") because my sample IDs are e.g. "ZE27_1"
bcin <- bcin %>%
	group_by(gene) %>%
	summarise(across(starts_with("Z"), sum))

#check colnames match
colnames(host) == colnames(bcin)

#bind readcounts together
df <- rbind(host, bcin)

#remove genes with no counts across any sample.
df <- df[rowSums(df[, -1] ==0) != ncol(df) - 1, ]

##### reformat as matrix where colnames are samples, rownames are genes, and cells contain numeric read counts
count_data <- df
#remove spaces from the readcounts
count_data[, -1] <- apply(count_data[, -1], 2, function(x) trimws(x))
#remove gene column
count_data <- count_data[,-1]
#convert readcounts to numeric values
count_data[!is.numeric(count_data)] <- lapply(count_data[!is.numeric(count_data)], as.numeric)
#convert to matrix
count_data <- as.matrix(count_data)
#put in gene rownames
row.names(count_data) <- df$gene

# Create a DGEList object
dge <- DGEList(counts = count_data)

# #Filter Low-expressed genes (not sure if I want to do this yet)
# keep <- rowSums(cpm(dge) >= 1) >= 2
# dge <- dge[keep, ]

#perform normalization using TMM
dge <- calcNormFactors(dge)

#access normalization factors
normalization_factors <- dge$samples$norm.factors

#extract normalized counts
normalized_counts <- cpm(dge, normalized.lib.sizes = T)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- rownames_to_column(normalized_counts,var = "gene")

#write csv
#For Csat
#write.csv(normalized_counts, "Bos_1/readcounts/norm_counts_all.csv", row.names = F, col.names = T)
#For Cpepo
write.csv(normalized_counts, "ZE_1/readcounts/norm_counts_all.csv", row.names = F, col.names = T)
