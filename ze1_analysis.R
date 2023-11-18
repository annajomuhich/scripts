library(tidyverse)
library(ggplot2)
library(ggrepel)
library(edgeR)

df <- read.csv("ZE_1/readcounts/norm_counts_all.csv")
sample_key <- read.csv("sample_info/cuc_rnaseq_sampleIDs.csv")

#pivot long
df <- df %>% pivot_longer(cols = !gene, names_to = "sample_ID", values_to = "count")

#join w sample key
df <- left_join(df, sample_key, by = "sample_ID")

####### % bcin reads of total reads
total_count <- df %>%
	group_by(sample_ID) %>%
	summarise(total_count = sum(count))

bcin_count <- df %>% 
	filter(grepl("^Bcin", gene)) %>%
	group_by(sample_ID) %>%
	summarise(bcin_count = sum(count))

bcin_pct <- left_join(total_count, bcin_count, by = "sample_ID") %>%
	mutate(bcin_pct = (bcin_count/total_count * 100))

bcin_pct <- left_join(bcin_pct, sample_key, by = "sample_ID")
bcin_pct <- arrange(bcin_pct, bcin_pct)

g <- ggplot(data = bcin_pct, aes(x = reorder(iso_name, bcin_pct), y = bcin_pct)) +
	geom_bar(stat = "identity", fill = "grey") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
				legend.position = "none") +
	xlab("Isolate") +
	ylab("% Bcin reads") +
	ggtitle("%Bcin reads of total in ZE_1")
		scale_x_discrete(labels = bcin_pct$iso_name)

ggsave("bcin_pct.pdf")

##### lesion size vs bcin pct
#quick mean lesion size calculation
lesion <- read.csv("../Cuc_Pheno/CucPhenoFull.csv")
colnames(lesion)[colnames(lesion) == "Isolate_number"] <- "iso_number"
lesion <- lesion %>%
	filter(Plant_geno_shorthand == "ZE")
means <- lesion %>%
	group_by(Plant_geno_shorthand, iso_number) %>%
	summarise(lesion_mean = mean(Lesion.Size.mm2),
						n = n_distinct(Lesion.Size.mm2))

bcin_pct$iso_number <- gsub("73", "Mock", bcin_pct$iso_number)
bcin_pct$iso_number <- gsub("74", "Mock", bcin_pct$iso_number)

bcin_pct <- left_join(bcin_pct, means, by = "iso_number")

ggplot(data = bcin_pct, aes(x = lesion_mean, y = bcin_pct)) +
	geom_point() +
	geom_smooth(method = "lm", se = F)

ggsave("lesionXbcin_pct.pdf",
			 width = 10,
			 height = 10)

##### PCA - Bcin reads

#filter so only bcin reads
bcin <- df %>%
	filter(grepl("^Bcin", gene))
means <- bcin_pct %>% select(sample_ID, lesion_mean)
bcin <- left_join(bcin, means, "sample_ID")

#Reformat for PCA
#pivot data wider again
bcin_forpca_cat <- bcin %>%
	pivot_wider(values_from = count, names_from = gene)
head(bcin_forpca_cat)

#remove categorical variables
bcin_forpca <- bcin_forpca_cat[,c(10:11173)]

#put in gene rownames
row.names(bcin_forpca) <- bcin_forpca_cat$sample_ID

# Create a DGEList object
dge <- DGEList(counts = bcin_forpca)

#run PCA
pca <- prcomp(dge)
#take a look
summary(pca)
#convert to dataframe
pca_df <- as.data.frame(pca$x)
#a few plots
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

ggplot(data = pca_df, aes(PC1,
													PC2,
													color = bcin_forpca_cat$lesion_mean)) +
	geom_point() +
	geom_text_repel(aes(label = bcin_forpca_cat$iso_name),
						position = position_jitter(width = 20, height = 20)) +
	theme_minimal() +
	labs(color = "Mean Lesion Size mm2",
			 title = "PCA of normalized Bcin reads from ZE_1")

ggsave("Bcin_PCA_lesions.pdf",
			 height = 19,
			 width = 28,
			 units = "cm")

#####try removing outliers
#1_04_12, B05_10, 01_03_22, 01_04_17, 01_05_22, KGB1
#pivot data wider again
bcin_forpca_cat <- bcin %>%
	pivot_wider(values_from = count, names_from = gene)
head(bcin_forpca_cat)
#remove outliers
bcin_forpca_cat <- bcin_forpca_cat %>%
	filter(sample_ID != "Z17_1") %>%
	filter(sample_ID != "Z39_1") %>%
	filter(sample_ID != "Z13_1") %>%
	filter(sample_ID != "Z19_1") %>%
	filter(sample_ID != "Z28_1") %>%
	filter(sample_ID != "Z54_1")
#remove categorical variables
bcin_forpca <- bcin_forpca_cat[,c(10:11173)]
#put in gene rownames
row.names(bcin_forpca) <- bcin_forpca_cat$sample_ID
# Create a DGEList object
dge <- DGEList(counts = bcin_forpca)
#run PCA
pca <- prcomp(dge)
#take a look
summary(pca)
#convert to dataframe
pca_df <- as.data.frame(pca$x)
#a few plots
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

ggplot(data = pca_df, aes(PC1,
													PC2,
													color = bcin_forpca_cat$lesion_mean)) +
	geom_point() +
	geom_text_repel(aes(label = bcin_forpca_cat$iso_name),
									position = position_jitter(width = 20, height = 20)) +
	theme_minimal() +
	labs(color = "Mean Lesion Size mm2",
			 title = "PCA of normalized Bcin reads from ZE_1 - no outliers")
ggsave("Bcin_PCA_lesions_noOutliers.pdf",
			 height = 19,
			 width = 28,
			 units = "cm")

##### PCA - Host reads
#filter so only bcin reads
host <- df %>%
	filter(grepl("^C", gene))
host <- left_join(host, means, "sample_ID")

#Reformat for PCA
#pivot data wider again
host_forpca_cat <- host %>%
	pivot_wider(values_from = count, names_from = gene)
head(host_forpca_cat)

#remove categorical variables
host_forpca <- host_forpca_cat[,c(10:11173)]

#put in gene rownames
row.names(host_forpca) <- host_forpca_cat$sample_ID

# Create a DGEList object
dge <- DGEList(counts = host_forpca)

#run PCA
pca <- prcomp(dge)
#take a look
summary(pca)
#convert to dataframe
pca_df <- as.data.frame(pca$x)
#a few plots
ggplot(data = pca_df, aes(PC1, PC2)) + geom_point()

ggplot(data = pca_df, aes(PC1,
													PC2,
													color = host_forpca_cat$lesion_mean)) +
	geom_point() +
	geom_text_repel(aes(label = host_forpca_cat$iso_name),
									position = position_jitter(width = 20, height = 20)) +
	theme_minimal() +
	labs(color = "Mean Lesion Size mm2",
			 title = "PCA of normalized Host reads from ZE_1")
ggsave("Host_PCA_Lesions.pdf",
			 height = 40,
			 width = 60,
			 units = "cm")
