#!/home/mimitair/miniforge3/envs/ibp/bin/python

"""
This script takes a binary output file from the cblaster tool, 
extracts the scaffold accession numbers from each hit and couples them
to their respective genome assembly accessions.
Then each genome assembly is downloaded and clustered based on similarity.
"""

import polars as pl
import sys
import subprocess

# If direct link to assembly exists:
#COMMAND1 = "esearch -db nucleotide -query '{scaffold}' | elink -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"

# If esearch returns no result and no direct link to assembly is present:
#COMMAND2 = "esummary -db nucleotide -id {scaffold} | xtract -pattern DocumentSummary -element BioSample | elink -db biosample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"

def get_scaffolds(path_to_binary: str) -> list:
	"""
	This function extracts the scaffold IDs from the cblaster binary output file.
	
	:param str path_to_binary: Path to the cblaster binary output file.
	:rtype list: A list containing the scaffold IDs.
	"""
	# Read the file as comma separated, extract the first column (containing the scaffold IDs) as a polars series and convert to a python list
	print("Extracting scaffold accession IDs...")
	scaffolds = pl.read_csv(path_to_binary).to_series(1).to_list()
	print(f"Extracted {len(scaffolds)} scaffold IDs")
	
	return scaffolds
	

def get_assemblies(scaffolds: list) -> list:
	"""
	This function obtains the genome assembly ID for each scaffold ID in the cblaster binary output file.

	:param list scaffolds: A list containing scaffold IDs.
	:rtype list: A list containing all associated genome assembly IDs.
	"""
	print("Extracting genome assembly IDs...")
	for scaffold in scaffolds:
		
		# First we have to check if the scaffold accession ID exists in the NCBI database
		# TODO use esummary or esearch? some scaffold IDs are recognized by esummary but not by esearch
		# TODO check for metagenome assembled genomes?
		
		# For each scaffold ID, we extract the corresponding genome assembly through the NCBI BioProject database:
		# TODO check if a genome assembly even exists
		COMMAND = f"esummary -db nucleotide -id {scaffold} | xtract -pattern DocumentSummary -element BioSample | elink -db biosample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"
		subprocess.run(COMMAND, shell=True)	
		
		# TODO Check exit status

def main():
	# Path to the cblatser binary output file, given as command line argument	
	path_to_binary = sys.argv[1]
	
	get_assemblies(get_scaffolds(path_to_binary))	
	#subprocess.run("esearch -db nucleotide -query NC_013893.1", shell=True)

if __name__ == "__main__":
	main()
	
