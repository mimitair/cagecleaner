#!/bin/bash

echo Downloading genomes...

# Download the dehydrated genome package:
datasets download genome accession $(cat assemblies.txt) --dehydrated

# Unzip the file in a folder called "genomes". This folder is located two steps up in the data folder
unzip -d ../../data/genomes ncbi_dataset.zip

# Remove the zip file
rm ncbi_dataset.zip

# Rehydrate:
datasets rehydrate --directory ../../data/genomes

echo Dereplicating genomes...

# Run skDER on the downloaded genomes:
skder -g ../../data/genomes/ncbi_dataset/data/GC* -o ../../data/skder_out

echo All done! The results can be found in /data/skder_out
