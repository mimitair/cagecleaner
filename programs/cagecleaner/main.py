#!/home/mimitair/miniforge3/envs/ibp/bin/python

"""
This script takes a binary output file from the cblaster tool, 
extracts the scaffold accession numbers from each hit and couples them
to their respective genome assembly accessions.
Then each genome assembly is downloaded and clustered based on similarity.
As a final step, the representative genomes are coupled back to the original hits in the cblaster binary output file.
A "cleaned" binary output file is then generated.
"""

import polars as pl
import sys
import subprocess
import os

def validate_file(file_path: str) -> None:
	"""
	Check if the input file exists and is readable.
	os.path.isfile()
	True: if the path exists and is a regular file
	False: if the path doesn't exist or if it is a directory or any other type of non-regular file.
	"""
	if not os.path.isfile(file_path):
		raise FileNotFoundError(f"Input file '{file_path}' not found or is not a valid file.")


def get_scaffolds(path_to_binary: str) -> list:
	"""
	This function extracts the scaffold IDs from the cblaster binary output file.
	
	:param str path_to_binary: Path to the cblaster binary output file.
	:rtype list: A list containing the scaffold IDs.
	"""
	print("Extracting scaffold accession IDs from the input file...")

	# Read the file as comma separated, extract the first column (containing the scaffold IDs) as a polars series and convert to a python list
	try:
		scaffolds = pl.read_csv(path_to_binary, truncate_ragged_lines=True).to_series(1).to_list()
		if not scaffolds:
			raise ValueError("No scaffold IDs were found in the file.")
	except Exception as e: 
		raise ValueError(f"Error reading the binary output file: {e}")
	
	print(f"Extracted {len(scaffolds)} scaffold IDs")
	return scaffolds
	

def get_assemblies(scaffolds: list) -> dict:
	"""
	This function obtains the genome assembly ID for each scaffold ID in the cblaster binary output file.
	The assembly IDs are then written to a file 'assemblies.txt' as a long string to be passed as an argument to a helper bash script

	:param list scaffolds: A list containing scaffold IDs.
	:rtype dict: A dictionary with scaffold IDs as keys and assembly IDs as values.
	"""
	print("Contacting the NCBI servers and retreiving genome assembly IDs for each hit. This may take a while...")
	
	# Create empty dictionary to store genome assembly IDs with scaffold IDs as keys:
	result = {}
	
	for scaffold in scaffolds:
		
		# TODO First we have to check if the scaffold accession ID exists in the NCBI database
		# TODO use esummary or esearch? some scaffold IDs are recognized by esummary but not by esearch
		# TODO check for metagenome assembled genomes?
		
		# For each scaffold ID, we extract the corresponding genome assembly through the NCBI BioProject database:
		# We use the NCBI E-direct utilities, available from the command line.
		# TODO check if a genome assembly even exists
		# TODO check if we can use bioproject id and get the genomes from there. We might not need to go all the way to accession level in this step
		
		try:
			# Check scaffold accession ID validity
			COMMAND = f"esummary -db nucleotide -id {scaffold} | xtract -pattern DocumentSummary -element BioSample" 
			process = subprocess.run(COMMAND, shell=True, capture_output=True, text=True, check=True)	
			biosample = process.stdout.strip()
		
			if not biosample: # Raise error when BioSample for the scaffold ID doesn't exist
				raise ValueError(f"No BioSample found for scaffold ID: {scaffold}")
			
			# Retrive assembly ID
			command  = f"elink -db biosample -target assembly -id {biosample} | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession"
			process = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
			assembly = process.stdout.strip()

			if not assembly: # Check assembly ID exist
				raise ValueError(f"No assembly ID found for BioSample: {biosample}")
			
			#TODO CHECK FOR DUPLICATES?	

			# Append the output of the above command to the dictionary:
			result[scaffold] = assembly
			#print(process.stdout.strip())	

		# Check the shell commands
		except subprocess.CalledProcessError as e:
			print(f"Error executing NCBI command for scaffold {scaffold}: {e}")
		except Exception as e:
			print(f"Error retrieving assembly for scaffold {scaffold}: {e}")
	
	if not result: 
		raise RuntimeError("No assembly IDs were retrieved. Please check your input and try again.")
		
	# Print how many assembly IDs were extracted. This should equal the amount of scaffold IDs
	print(f"Extracted {len(result.values())} assembly IDs ")	
	print(result)
	return result
	

def cluster_genomes(assemblies: list) -> None:
	"""
	This function takes a list of genome assembly IDs and calls a helper bash script that downloads and dereplicates the genomes

	:param list assemblies: A list of assembly IDs
	"""
	if not assemblies:
		raise ValueError("No assemblies provided for clustering")
	
	# First we write the assembly IDs to a file that can be used by the helper bash script:	
	with open('assemblies.txt', 'w') as file:
		file.write(' '.join(str(assembly) for assembly in assemblies))

	try:	
		# Run the bash script to download and cluster genomes:
		subprocess.run("./helper.sh", shell=True, check=True)
		# The helper script writes the assembly IDs of the dereplicated genomes to a file called "dereplicated_assemblies.txt".	
	
	except subprocess.CalledProcessError as e:
		raise RuntimeError(f"Error in clustering genomes: {e}")
	if not os.path.exists("dereplicated_assemblies.txt"):
		raise FileNotFoundError("The clustering output file 'dereplicated_assemblies.txt' was not generated.")
	
	#TODO CHECK OUTPUT VALIDITY CATCH ERRORS
	

def main():
	# Check if the argument is provided
	if len(sys.argv) < 2: # at least one argument is given
		print("Error: Missing required argument: ")
		#print("Usage: python main.py <path_to_binary_output>")
		sys.exit(1)

	# Path to the cblatser binary output file, given as command line argument	
	path_to_binary = sys.argv[1]
	
	# Check validity of the input file with validate_file() function
	validate_file(path_to_binary)

	try: # Check if the  any error occurs during the execution of following steps

		# First we capture the dictionary with scaffold:assembly pairs:
		scaff_ass_pairs = get_assemblies(get_scaffolds(path_to_binary))

		# Then we download and cluster the genomes:
		cluster_genomes(scaff_ass_pairs.values())

		# The dereplicated genome assembly IDs are now stored in the file "dereplicated_assemblies.txt" on separate lines
		# We read this file and couple the assembly IDs back to their scaffold IDs through the dictionary.

		print("Generating cleaned binary output file...")

		with open('dereplicated_assemblies.txt', 'r') as file:
			for line in file:
				# Get the assembly accession out of the file name:
				# file names for the genomes are structured as follows: [GCF][_][nine digits][.][version][_][ASMblabla][.fna]
				# We only need the part before the second "_" character. Hence, we split the string two times based on the "_" character
				# and then join the first two parts back toegther. This results in the assembly ID as it is found in the dictionary
				assembly = "_".join(line.split("_", 2)[:2])
				
				# Now we extract the corresponding scaffold ID. The index of the assembly ID in the dictionary is calculated.
				# Then, this index is passed onto the list of keys to get the element at that specific index.
				scaffold = list(scaff_ass_pairs.keys())[list(scaff_ass_pairs.values()).index(assembly.strip())]
				
				# Now we have to match this scaffold ID with the corresponding hit in the original cblaster binary file.
				subprocess.run(f"grep '{scaffold}' {path_to_binary} > ../../data/output/cleaned_binary.csv", shell=True, check=True)
		print("All done! Results are written to '/data/output/cleaned_binary.csv'. Thank you for using cagecleaner :)")

	# If an exception is raised anywhere in the code within try block, raise error
	except Exception as e:
		print(f"Error in pipeline : {e}")
		sys.exit(1) # exit the script with a status code 1


if __name__ == "__main__":
	main()
	
