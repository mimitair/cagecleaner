#!/bin/bash

echo The helper bash script has been initiated!

echo Preparing to download and dereplicate genomes...

#If intermediary/output files and folders exist, they need to be deleted:
if [ -f "../../data/output/cleaned_binary.csv" ]; then
    rm ../../data/output/cleaned_binary.csv
fi

if [ -f "../../data/output/clusters.txt" ]; then
    rm ../../data/output/clusters.txt
fi

if [ -d "../../data/genomes/all" ]; then
    rm -r ../../data/genomes/all
fi 

if [ -d "../../data/skder_out" ]; then
    rm -r ../../data/skder_out
fi

# Set some useful variables:
assembly_count=$(cat assemblies.txt | wc -w)
batch_count=$(cat assemblies.txt | wc -l)

# The first argument is the ani-cutoff for skDER:
pi_cutoff=$1

echo Got $assembly_count assembly IDs. Downloading in $batch_count batches...
 
# Create the directory to store all genomes
mkdir ../../data/genomes/all

while read line; do
    # Removing files that are autmoatically donwnloaded by NCBI datasets. 
    # This way datasets does not stop and ask the user to override the files.
    if [ -d "../../data/genomes/ncbi_dataset" ]; then
        rm -r ../../data/genomes/ncbi_dataset
    fi

    if [ -f "../../data/genomes/md5sum.txt" ]; then
        rm ../../data/genomes/md5sum.txt
    fi
    
    if [ -f "../../data/genomes/README.md" ]; then
        rm ../../data/genomes/README.md
    fi 

    # Download the dehydrated genome package. xargs removes the trailing whitespace:
    datasets download genome accession $(echo "$line" | xargs) --dehydrated

    # Unzip the file in a folder called "genomes". This folder is located two steps up in the data folder
    unzip -d ../../data/genomes ncbi_dataset.zip

    # Rehydrate:
    datasets rehydrate --directory ../../data/genomes

    # Put all genomes into one directory
    mv ../../data/genomes/ncbi_dataset/data/GC*/* ../../data/genomes/all
done < assemblies.txt

echo Downloading finished!
echo $(ls ../../data/genomes/all | wc -w) genomes in /data/genomes/all.
echo Directory size: $(du --human-readable -s ../../data/genomes/all).

echo Dereplicating genomes with percent identity cutoff of $pi_cutoff.

# Run skDER on the downloaded genomes:
#TODO play with threads and ani cutoffs
skder -g ../../data/genomes/all/* -o ../../data/skder_out -i $pi_cutoff

echo Dereplication done! $(ls ../../data/genomes/all | wc -w) genomes were converted into $(ls ../../data/skder_out/Dereplicated_Representative_Genomes | wc -w) genomes. 

# Now we have to write the accession numbers of the dereplicated genomes to a file:
ls ../../data/skder_out/Dereplicated_Representative_Genomes > dereplicated_assemblies.txt

# CLEAN UP:
rm ncbi_dataset.zip assemblies.txt 
