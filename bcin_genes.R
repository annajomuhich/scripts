library(dplyr)
library(ggplot2)

#### Bos_1

df <- read.csv("Bos_1/readcounts/norm_counts_all.csv")
sample_key <- read.csv("sample_info/cuc_rnaseq_sampleIDs.csv")

#pivot long
df <- df %>% pivot_longer(cols = !gene, names_to = "sample_ID", values_to = "count")

#join w sample key
df <- left_join(df, sample_key, by = "sample_ID")

#filter to retain only Botrytis genes
df <- df %>% filter(grepl("^Bcin", gene))

#filter to remove genes with count of 0
df <- df %>% filter(count != 0)

#summarise count of bcin genes with any signal by isolate
genes <- df %>% group_by(iso_name) %>%
	summarise(bcin_gene_count = length(gene))

#plot
ggplot(data = genes, aes(x = reorder(x = iso_name, bcin_gene_count), y = bcin_gene_count)) +
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle=90, hjust=1)) +
	labs(
		x = "Isolate",
		title = "Bcin gene count in each sample - Bos_1")

#dev.size("cm")
ggsave("bcin_gene_count_Bos1.pdf", width = 33, height = 20, units = "cm")

### Gene by gene count for a few isolates
df_sel <- df %>% filter(iso_name == c("DavisNavel",
														"KernA2",
														"01_03_19",
														"01_04_17", 
														"Mock_0HAI",
														"Mock_16HAI",
														"Peachy",
														"Triple3"))
#compiling their gene count qualitative info based on the bcin_gene_count figure
gene_count_info <- data.frame(
	iso_name = c("DavisNavel",
							 "KernA2",
							 "01_03_19",
							 "01_04_17", 
							 "Mock_0HAI",
							 "Mock_16HAI",
							 "Peachy",
							 "Triple3"),
	gene_count = c("low",
								 "low",
								 "mid",
								 "mid",
								 "mock",
								 "mock",
								 "high",
								 "high")
)

df_sel <- left_join(df_sel, gene_count_info, by = "iso_name")

ggplot(df_sel, aes(x = gene, y = count, color = gene_count)) +
	geom_point(size = 0.5) +
	facet_wrap(~ iso_name, scales = "fixed") +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

ggsave("gene_counts_select1.pdf",
			 height = 20, width = 33, units = "cm")

df_sel2 <- df_sel %>% filter(gene_count == c("high", "mock"))

ggplot(df_sel2, aes(x = gene, y = count, color = gene_count)) +
	geom_point(size = 0.5) +
	facet_wrap(~ iso_name, scales = "fixed") +
	theme(axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

ggsave("gene_counts_select2.pdf",
			 height = 20, width = 33, units = "cm")

#### ZE_1

df <- read.csv("ZE_1/readcounts/norm_counts_all.csv")
sample_key <- read.csv("sample_info/cuc_rnaseq_sampleIDs.csv")

#pivot long
df <- df %>% pivot_longer(cols = !gene, names_to = "sample_ID", values_to = "count")

#join w sample key
df <- left_join(df, sample_key, by = "sample_ID")

#filter to retain only Botrytis genes
df <- df %>% filter(grepl("^Bcin", gene))

#filter to remove genes with count of 0
df <- df %>% filter(count != 0)

#summarise count of bcin genes with any signal by isolate
genes <- df %>% group_by(iso_name) %>%
	summarise(bcin_gene_count = length(gene))

#plot
ggplot(data = genes, aes(x = reorder(x = iso_name, bcin_gene_count), y = bcin_gene_count)) +
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle=90, hjust=1)) +
	labs(
		x = "Isolate",
		title = "Bcin gene count in each sample - ZE_1")

#dev.size("cm")
ggsave("bcin_gene_count_ZE1.pdf", width = 33, height = 20, units = "cm")
