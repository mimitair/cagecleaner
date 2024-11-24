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
import re

def validate_file(file_path: str) -> bool:
	"""
    Validate the input file.
    - Check if the file exists and is readable.
    - Ensure the file ends with .csv.
    - Verify it has at least 5 columns, including essential columns like 'organism', 'scaffold', 'start'.
    - Check column names match expectations.
    - Ensure it contains at least a header row and two data rows.

    :param file_path: Path to the input file.
    :return: True if validation passes, False otherwise.
	"""
	# Check if file exists and readable
	if not os.path.isfile(file_path):
		print(f"Input file '{file_path}' not found or is not a valid file.")
		return False

	try: 
		# Check if file ends with .csv
		if not file_path.endswith(".csv"): 
			print(f"Error: File '{file_path}' does not have a .csv extension.")
			return False

		# Check the file has at least 6 columns
		data = pl.read_csv(file_path, truncate_ragged_lines=True)
		columns = data.columns
		if len(columns) < 6:
			print(f"Error: File '{file_path}' has fewer than 6 columns. Please check the file structure.")
			return False

		# Check column names (customization needed)
		required_columns = {"organism", "scaffold", "start", "end", "gene", "evalue"}
		if not required_columns.issubset(set(columns)):
			print(f"Error: File '{file_path}' is missing required columns: {required_columns - set(columns)}")
			return False

		# Check for header and minimum of 2 rows
		if data.height < 2: # need to be more than 3 rows
			print(f"Error: File '{file_path}' must have at least 2 data row")
	except Exception as e:
		print(f"Error reading file '{file_path}': {e}")
		return False

	return True

def get_scaffolds(path_to_binary: str) -> list:
	"""
	This function extracts the scaffold IDs from the cblaster binary output file.
	
	:param str path_to_binary: Path to the cblaster binary output file.
	:rtype list: A list containing the scaffold IDs.
	"""
	print("--- STEP 1 --- \nExtracting scaffold IDs from the input file...")

	# Read the file as comma separated, extract the first colum n (containing the scaffold IDs) as a polars series and convert to a python list
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
	print("--- STEP 2 --- \nContacting the NCBI servers and retreiving genome assembly IDs for each hit. This may take a while...")
    
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
			command = (
                f"esummary -db nucleotide -id {scaffold} | "
                "xtract -pattern DocumentSummary -element BioSample | "
                "elink -db biosample -target assembly | "
                "efetch -format docsum | "
                "xtract -pattern DocumentSummary -element AssemblyAccession"
            )
			process = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)    
			
			# Get the assembly ID from the subprocess output
			assembly = process.stdout.strip()

			# Validate assembly ID format
			if not assembly.startswith("GCF") or not re.match(r"GCF_\d{9}\.\d+", assembly): # re. regex
				raise ValueError(f"Invalid Assembly ID format: {assembly} for scaffold: {scaffold}")

			# This if statement prevents duplicate assembly IDs from populating the dictionary.     
			if process.stdout.strip() in result.values():
				continue
			else:
				result[scaffold] = assembly

		# Check the shell commands
		except subprocess.CalledProcessError as e:
			print(f"Error executing NCBI command for scaffold {scaffold}: {e.stderr}") # Error Reporting: captures and prints stderr output.
		except Exception as e:
			print(f"Error retrieving assembly for scaffold {scaffold}: {e}")
	
	if not result: 
		raise RuntimeError("No assembly IDs were retrieved. Please check your input and try again.")

	# Print how many assembly IDs were extracted. This should equal the amount of scaffold IDs
	print(f"Obtained {len(result.values())} assembly IDs ")	
	#print(result)
	
	return result
	

def dereplicate_genomes(assemblies: list) -> None:
	"""
	This function takes a list of genome assembly IDs and calls a helper bash script that downloads and dereplicates the genomes

	The assembly IDs are written to a file 'assemblies.txt' (space separated on one line) 
    to be passed as an argument to a helper bash script.
        
    A file called "dereplicated_assemblies.txt" is then generated by the helper script.

	:param list assemblies: A list of assembly IDs
	"""

	print("--- STEP 3 --- \nDownloading and dereplicating genomes.")

	if not assemblies:
		raise ValueError("No assemblies provided for dereplication")
	
	# First we write the assembly IDs to a file that can be used by the helper bash script:	
	with open('assemblies.txt', 'w') as file:
		file.write(' '.join(str(assembly) for assembly in assemblies))

	try:	
		# Run the bash script to download and cluster genomes:
		subprocess.run("./helper.sh", shell=True, check=True)
		# The helper script writes the assembly IDs of the dereplicated genomes to a file called "dereplicated_assemblies.txt".	
	
	except subprocess.CalledProcessError as e:
		print(f"Subprocess error during dereplication: {e.stderr.strip()}")
		raise RuntimeError("Error in genome dereplication.")
	
	if not os.path.exists("dereplicated_assemblies.txt"):
		raise FileNotFoundError("The clustering output file 'dereplicated_assemblies.txt' was not generated.")
	
	#TODO CHECK OUTPUT VALIDITY CATCH ERRORS

def extract_cleaned_hits(dereplicated_assemblies: str, path_to_binary:str, scaff_ass_pairs:dict) -> None:
        """
        This function reads the 'dereplicated_assemblies.txt' file generated by dereplicate_genomes().
        Each assembly ID is coupled to the corresponding scaffold ID via a dictionary.
        The scaffold ID is then used to grep the original cblaster binary output file and obtain the relevant hits.

        :param str dereplicated_assemblies: File path to the 'dereplicated_assemblies.txt' file.
        :param str path_to_binary: File path to the cblaster binary output file.
        :param dict scaff_ass_pairs: A dictionary containing scaffold IDs as keys and assembly IDs as values.
        """  
        print("--- STEP 4 --- \nGenerating cleaned binary output file...")
        
        # The dereplicated genome assembly IDs are now stored in the file "dereplicated_assemblies.txt" on separate lines
        # We read this file and couple the assembly IDs back to their scaffold IDs through the dictionary.
        with open(dereplicated_assemblies, 'r') as file:
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
                        subprocess.run(f"grep '{scaffold}' {path_to_binary} >> ../../data/output/cleaned_binary.csv", shell=True, check=True)


def main():
	# Check if the argument is provided
	if len(sys.argv) < 2: # at least one argument is given
		print("Error: Missing required argument: ")
		sys.exit(1)
	
	# Path to the cblatser binary output file, given as command line argument	
	path_to_binary = sys.argv[1]
	
	#TODO CHECK VALIDITY OF FILE (is it csv? is it cblaster output?)
	#Step0: Check validity of the input file with validate_file() function
	if not validate_file(path_to_binary):
		print("File validation failed. Exiting.")
		sys.exit(1)

	try: # Check if the  any error occurs during the execution of following steps

		# Step 1: Extract scaffold IDs
		scaffolds = get_scaffolds(path_to_binary)

		# Step 2: Retrieve assembly IDs
		scaff_ass_pairs = get_assemblies(scaffolds)

		# Step 3: Download and cluster the genomes
		# This generates a file called "dereplicated_assemblies.txt"
		dereplicate_genomes(scaff_ass_pairs.values())

        # Step 4: Extract BGC hits from the original input file based on dereplicated genomes:
		extract_cleaned_hits('dereplicated_assemblies.txt', path_to_binary, scaff_ass_pairs)

		print("All done! Results are written to '/data/output/cleaned_binary.csv'. Thank you for using cagecleaner :)")

	except Exception as e:
		# Handle any errors that occur during the execution of the pipeline
		print(f"Error in pipeline : {e}")
		sys.exit(1) # exit the script with a status code 1


if __name__ == "__main__":
	main()
	
