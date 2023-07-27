#!/bin/bash

# Set the path to the parent directory containing the subdirectories (directories with data)
data_dir="/home/garvena/projects/def-dmberman/garvena/singlecellseq/Data/InputData/Lai"

# Change to the data directory
cd "$data_dir"

# Loop through the subdirectories
for subdir in */; do
    echo "Processing files in $subdir"

    # Change to the subdirectory
    cd "$subdir"

    # Loop through the files and rename them
    for file in *.fastq.gz; do
        if [[ $file == *_S1_L*_R1_* ]]; then
            new_name="${file/_S1_L*_R1_/_R1_}"
        elif [[ $file == *_S1_L*_R2_* ]]; then
            new_name="${file/_S1_L*_R2_/_R2_}"
        else
            echo "Skipping file: $file (Unexpected filename format)"
            continue
        fi

        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed: $file -> $new_name"
    done

    # Change back to the parent directory for the next subdirectory
    cd "$data_dir"
done