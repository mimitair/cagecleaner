#!/bin/bash

# Download the dehydrated genome package:
datasets download genome accession $(cat assemblies.txt) --dehydrated

# Unzip the file in a folder calles "genomes"
unzip -d genomes ncbi_dataset.zip

# Rehydrate:
datasets rehydrate --directory genomes

# Run skDER on the downloaded genomes:
#skder -g genomes/ncbi_dataset/data/GC*
