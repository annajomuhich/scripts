#!/bin/bash

##### Alignment Summary
##### Anna Jo Muhich
##### August 2023

#Loop to pull out alignment values from .txt files within the /bams outputs
#In order to evaluate alignment rates across many samples

cd bams

# Initialize an array to hold the data
data=()

# Loop through each directory in the current directory
for dir in */; do
  # Loop through each .txt file in the directory
  for file in "$dir"*.txt; do
    # Get the last line of the file
    last_line=$(tail -n 1 "$file")
    # Extract the numerical value from the last line
    value=$(echo "$last_line" | grep -Eo '[0-9]+(\.[0-9]+)?')
    # Add the filename and value to the data array
    data+=("$file,$value")
  done
done

# Save the table as a CSV file
printf '%s\n' "${data[@]}" > alignment_summary.csv

#Reformat the alignment_summary.csv in 2b_alignment_summary2.R