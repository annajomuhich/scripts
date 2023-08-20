library(tidyverse)

df <- read.csv("alignment_summary.csv", header = F)

#sep bcin and plant alignments
bcin <- df %>%
	filter(str_detect(V1, "Bcin"))
host <- df %>%
	filter(!str_detect(V1, "Bcin"))

#adjust colnames
colnames(host) <- c("ID", "host_aligment_pct")
colnames(bcin) <- c("ID", "bcin_alignment_pct")

#put host and bcin data together
df <- cbind(host, bcin[,2])

#adjust colnames
colnames(df) <- c("ID", "host_alignment_pct", "bcin_alignment_pct")

#calc bcin pct alignment of total reads
df <- df %>%
	mutate(bcin_alignment_oftotal =
				 	(100-host_alignment_pct)*(bcin_alignment_pct/100))

# #If you want to join the alignment summaries to a more detailed sample key:
# #get columns for plate location and genome
# df <- df %>%
# 	mutate(plate_location = str_extract(string = ID,
# 																			pattern = "^[^_]+"))
# 
# #read in sample key
# key <- read.csv("../sample_info/sample_key.csv") %>%
# 	select(!X)
# 
# #join key
# df <- df %>%
# 	left_join(x = df,
# 						y = key,
# 						by = "plate_location")

#overwrite alignment summary
df %>%
	write.csv("alignment_summary.csv")
