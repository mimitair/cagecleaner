#!/bin/bash

#TODO If output folders exist, they need to be deleted:


echo Got $(cat assemblies.txt | wc -w) genomes...

# Download the dehydrated genome package:
datasets download genome accession $(cat assemblies.txt) --dehydrated

# Unzip the file in a folder called "genomes". This folder is located two steps up in the data folder
unzip -d ../../data/genomes ncbi_dataset.zip

# Remove the zip file
rm ncbi_dataset.zip

# Rehydrate:
datasets rehydrate --directory ../../data/genomes

# Put all genomes into one directory
mkdir ../../data/genomes/all
mv ../../data/genomes/ncbi_dataset/data/GC*/* ../../data/genomes/all


# Run skDER on the downloaded genomes:
#TODO play with threads and ani cutoffs
skder -g ../../data/genomes/all/* -o ../../data/skder_out

echo Dereplication done! $(cat assemblies.txt | wc -w) genomes were converted into $(ls ../../data/skder_out/Dereplicated_Representative_Genomes | wc -w) genomes. 

# Delete the assemblies.txt file
rm assemblies.txt

# Now we have to write the accession numbers of the dereplicated genomes to a file:
ls ../../data/skder_out/Dereplicated_Representative_Genomes > dereplicated_assemblies.txt

