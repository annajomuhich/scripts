#install.packages("gtools")
library(dplyr)
library(gtools)

df <- read.csv("sample_info/cuc_rnaseq_sampleIDs_batch.csv")

#remove instances where no value for batch
df <- df %>% filter(batch != "")

#read in list of well IDs
well_ID <- read.csv("sample_info/well_IDs.csv")

#get list of batch names
batch_list <- unique(df$batch)

#get a blank df
batches <- data.frame(batch = character(),
											order_for_plating = numeric(),
											well_ID = character(),
											spl = character())

#Loop to generate individual randomized batch lists
for(batch_no in batch_list) {
	#filter df by batch
	df_filt <- df %>% filter(batch == batch_no)
	#get list of samples in the batch
	spl <- df_filt$sample_ID
	#add a blank tot he list
	spl <- c(spl, "blank")
	#randomize the list
	spl <- sample(spl)
	#put in dataframe with well ID
	spl <- cbind(well_ID, as.data.frame(spl))
	#add batch column
	spl$batch <- batch_no
	#reorder columns
	spl <- spl %>% select(batch, everything())
	#append to batches df
	batches <- rbind(batches, spl)
	#sort
	sorted <- mixedsort(spl$spl)
	#save a sorted version of the sample list
	assign(batch_no, sorted)
}

#save sorted batches in case i want it for printing
sorted_batches <- data.frame(Bos_1, Bos_2, Bos_3,
					 ZE_1, ZE_2, ZE_3)
write.csv(sorted_batches, "sorted_batch_lists.csv", row.names = F)

#save final batch lists
write.csv(batches, "batch_lists.csv", row.names = F)
