library(tidyverse)

setwd("./sample_info/")
iso <- read.csv("72Iso_Assign_April2022.csv")

#add mock
mock <- data.frame(Isolate_number = "Mock", Isolate_name = "Mock")
iso <- rbind(iso, mock)

#make 6 replicates
iso <- rbind(iso, iso, iso, iso, iso, iso)

#randomize the isolate list
set.seed(123)
iso_rand <- sample(iso$Isolate_number) %>%
	as.data.frame()
colnames(iso_rand) <- "Isolate_number"

#join isolate names in
iso <- read.csv("72Iso_Assign_April2022.csv")
iso$Isolate_number <- as.character(iso$Isolate_number)
df <- left_join(iso_rand,
					iso,
					by = "Isolate_number")

write.csv(df, "cuc_rnaseq_sampleIDs.csv")

#Manually filled in tray info in Excel.

#reuploading:
df <- read.csv("cuc_rnaseq_sampleIDs.csv")

#I want to see if the same isolate is occurring multiple times within a tray on any of the trays with this present configuration
# Function to check for repetition
check_repetition <- function(df) {
	result <- sapply(split(df$iso_number, df$tray_number), function(x) {
		any(duplicated(x))
	})
	return(result)
}

repetition_check <- check_repetition(df)
print(repetition_check)

#there are many repetitions. I manually swapped out isolate positions in excel, then reuploaded to check again.
#once confirmed no repetitions, manually copied the iso layout for Bos IDs and saved.

#Manually added geno, rep number, and 16HAI plus 0HAI Mock samples

#Now want to generate a unique sample ID in the format "geno_isonumber_repnumber", e.g. Z1_1

df <- read.csv("./sample_info/cuc_rnaseq_sampleIDs.csv")
df <- df %>%
	mutate(sample_ID = paste0(str_extract(genotype, "^."),
														iso_number,
														"_",
														rep)) %>%
	select(sample_ID, genotype, iso_number, iso_name, rep, everything())
df %>%
	write.csv("./sample_info/cuc_rnaseq_sampleIDs.csv")

#Making printout for tray labels
df <- read.csv("sample_info/cuc_rnaseq_sampleIDs.csv") %>% select(!X)

df <- df %>%
	mutate(tray_section = paste0(tray_number,tray_section)) %>%
	filter(genotype == "Bos", iso_number != 73) %>%
	select(iso_number, tray_section, inoc_position) %>%
	pivot_wider(names_from = inoc_position,
							values_from = iso_number)

df <- rbind(df, df)
df %>%
	write.csv("tray_labels.csv")
