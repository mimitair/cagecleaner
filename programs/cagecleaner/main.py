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

def get_scaffolds(path_to_binary: str) -> list:
	"""
	This function extracts the scaffold IDs from the cblaster binary output file.
	
	:param str path_to_binary: Path to the cblaster binary output file.
	:rtype list: A list containing the scaffold IDs.
	"""
	print("Extracting scaffold accession IDs from the input file...")

	# Read the file as comma separated, extract the first column (containing the scaffold IDs) as a polars series and convert to a python list
	scaffolds = pl.read_csv(path_to_binary, truncate_ragged_lines=True).to_series(1).to_list()

	print(f"Extracted {len(scaffolds)} scaffold IDs")

	return scaffolds
	

def get_assemblies(scaffolds: list) -> list:
	"""
	This function obtains the genome assembly ID for each scaffold ID in the cblaster binary output file.

	:param list scaffolds: A list containing scaffold IDs.
	:rtype list: A list containing all associated genome assembly IDs.
	"""
	print("Contacting the NCBI servers and retreiving genome assembly IDs for each hit. This may take a while...")
	
	# Create empty list to store genome assembly IDs:
	result = []
	
	for scaffold in scaffolds:
		
		# First we have to check if the scaffold accession ID exists in the NCBI database
		# TODO use esummary or esearch? some scaffold IDs are recognized by esummary but not by esearch
		# TODO check for metagenome assembled genomes?
		
		# For each scaffold ID, we extract the corresponding genome assembly through the NCBI BioProject database:
		# We use the NCBO E-direct utilities, available from the command line.
		# TODO check if a genome assembly even exists
		# TODO check if we can use bioproject id and get th genomes from there. We might not need to go all the way to accession level in this step
		COMMAND = f"esummary -db nucleotide -id {scaffold} | xtract -pattern DocumentSummary -element BioSample | elink -db biosample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"
		process = subprocess.run(COMMAND, shell=True, capture_output=True, text=True, check=True)	
		
		# Append the output of the above command to the result list:
		result.append(process.stdout.strip())

	# Afterwards, we write the genome accession assemblies to a file as a long string:
	with open('assemblies.txt', 'w') as file:
		file.write(' '.join(str(assembly) for assembly in result))
	
	# Run the bash script to download genomes:
	subprocess.run("./download_genomes.sh", shell=True, check=True)
	
	return result

def main():
	# Path to the cblatser binary output file, given as command line argument	
	path_to_binary = sys.argv[1]

		
	get_assemblies(get_scaffolds(path_to_binary))	

if __name__ == "__main__":
	main()
	
