

"""
This script takes a binary and summary output file from the cblaster tool.
Secondly, each scaffold ID in the binary output is coupled to an assembly ID through edirect utilities.
Thirdly, each genome assembly is downloaded and dereplicated using the skDER tool.
Finally, the representative genomes are coupled back to the original hits in the cblaster binary/summary output file.
A "cleaned" binary output file and a file containing the corresponding clusters is then generated.
"""

import polars as pl
import sys
import subprocess
import re
import os  # For validating file paths
import time

def validate_input_files(path_to_binary: str, path_to_summary:str) -> bool:
        """
        This function takes the path to a binary and summary file as input and validates if the files are useable for downstream analysis.
        Does the file exist? Is it in .csv format? Does it have the necessary columns? Does the summary file match the binary file?

        :param:path_to_binary: Path to the binary file.
        :param:path_to_summary: Path to summary file.
        :rtype bool: True if all checks pass.
        """
        
        print("--- STEP 1: Validating input files. ---")
        print("Validating binary file.")
        
        # First we check if the file exists:
        if not os.path.isfile(path_to_binary): 
            raise FileNotFoundError(f"'{path_to_binary}' does not exist. Please check if the file path name is correct and try again.")
        
        # Then we check if it ends in .csv:
        assert path_to_binary.endswith(".csv"), "File must be in .csv format. A script 'to_csv.sh' is provided to convert CAGECAT binary output files to csv format."

        # Now we check if the file has at least 6 columns:
        data = pl.read_csv(path_to_binary, truncate_ragged_lines=True)
        assert len(data.columns) > 6, "The amount of columns does not correspond to a cblaster binary output file. Please check the file structure. At least 6 columns should be present."
        
        # Check if the required column names are present:
        assert {"Organism", "Scaffold", "Start", "End", "Score"} < set(data.columns), "The column names do not correspond to a cblaster binary output file. Please check the file structure.\n The columns 'Organims', 'Scaffold', 'Start', 'End' and 'Score' should be present."

        # Check if there are at least two hits:
        assert data.height > 3, "We are expecting at least two hits in the cblaster binary output file."
        
        print("CHECK") 
        
        print("Validating summary file.")

        if not os.path.isfile(path_to_summary):
            raise FileNotFoundError(f"'{path_to_summary}' does not exist. Please check if the file path is correct and try again.")
        
        # The amount of appearances of the word "Cluster" in the summary file should correspond to the amount of lines
        # in the binary file:
        #with open (path_to_summary, 'r') as file:
            #assert len(re.findall("Cluster", file.read())) == data.height, "Binary and summary file do no match." 
        
        print("CHECK")
        
        return True


def get_scaffolds(path_to_binary: str) -> list:
        """
        This function extracts the scaffold IDs from the cblaster binary output file.
        
        :param str path_to_binary: Path to the cblaster binary output file.
        :rtype list: A list containing the scaffold IDs.
        """
        
        print("--- STEP 2: Extracting scaffold IDs from the binary input file. ---")

        # Read the file as comma separated, extract the first column (containing the scaffold IDs) as a polars series and convert to a python list
        scaffolds = pl.read_csv(path_to_binary, truncate_ragged_lines=True).to_series(1).to_list()

        print(f"Extracted {len(scaffolds)} scaffold IDs")
        
        return scaffolds
        

def get_assemblies(scaffolds: list) -> dict:
        """
        This function obtains the genome assembly ID for each scaffold ID obtained by get_scaffolds().

        :param list scaffolds: A list containing scaffold IDs.
        :rtype dict: A dictionary with scaffold IDs as keys and assembly IDs as values.
        """
        
        # Set begin time:
        begin = time.time()
    
        print("--- STEP 3: Contacting the NCBI servers and retreiving genome assembly IDs for each scaffold. This may take a while. ---")
                
        # Create empty dictionary to store genome assembly IDs with scaffold IDs as keys:
        result = {}

        # A counter to count duplicate assembly IDs:
        counter = 0
        
        # We download in chunks of 100:
        #chunks = [scaffolds[i:i+2] for i in range(0, len(scaffolds), 2)]
        #print(chunks)
        for scaffold in scaffolds:
                
                # TODO First we have to check if the scaffold accession ID exists in the NCBI database
                
                # For each scaffold ID, we extract the corresponding genome assembly through the NCBI BioProject database:
                # We use the NCBI E-direct utilities, available from the command line.
                # TODO check if a genome assembly even exists
                process = subprocess.run(f"esummary -db nucleotide -id {scaffold} | xtract -pattern DocumentSummary -element BioSample | elink -db biosample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession", shell=True, capture_output=True, text=True, check=True)    
         
                # Check if the output is indeed something that starts with "GC[F,A]_":
                if process.stdout.strip().startswith("GC"):
                    if process.stdout.strip() in result.values():
                        counter += 1
                        continue
                    else:
                        result[scaffold] = process.stdout.strip()
                else:
                    print(f"An error occurred in retreiving the assembly ID for {scaffold}.")
                    #raise ValueError("s")
                    continue
        # Print how many assembly IDs were extracted. This should equal the amount of scaffold IDs
        print(f"Obtained {len(result.values())} assembly IDs.\nGot {counter} duplicate assembly IDs.") 

        # Set end time:
        end = time.time()
        
        print(f"Time elapsed: {round(end - begin, 2)} s")

        return result
        

def dereplicate_genomes(assemblies: list, pi_cutoff:float=99.0) -> None:
        """
        This function takes a list of genome assembly IDs and calls a helper bash script that downloads and dereplicates the genomes.
        
        The assembly IDs are written to a file 'assemblies.txt' in batches (space separated on a line) 
        to be passed as an argument to a helper bash script.
        
        A file called "dereplicated_assemblies.txt" is then generated by the helper script.

        :param list assemblies: A list of assembly IDs
        :param float pi_cutoff: The percent identity cutoff for dereplicating genomes (see skDER docs)
        """
        
        assert 0 < pi_cutoff <= 100, "Percent identity cutoff must be between 0 and 100. E.g.: 99.0"
    
        # Set begin time:
        begin = time.time()

        print("--- STEP 4: Downloading and dereplicating genomes. ---")

        # First we write the assembly IDs to a file that can be used by the helper bash script: 
        with open('assemblies.txt', 'a') as file:
                for i in range(0, len(assemblies), 300):
                    file.write(' '.join(str(assembly) for assembly in assemblies[i:i+300]))
                    file.write("\n")
 
        # Run the bash script to download and cluster genomes:
        subprocess.run(["bash", "helper.sh", str(pi_cutoff)], check=True)
        
        #TODO CHECK OUTPUT VALIDITY CATCH ERRORS
        
        # Set end time:
        end = time.time()
    
        print(f"Time elapsed: {round(end - begin, 2)} s")
        
        return None
 
def get_dereplicated_scaffolds(path_to_dereplicated_assemblies: str, scaff_ass_pairs:dict) -> list:
        """
        This function reads the 'dereplicated_assemblies.txt' file generated by dereplicate_genomes().
        Each assembly ID is coupled to the corresponding scaffold ID via a dictionary.

        :param str dereplicated_assemblies: File path to the 'dereplicated_assemblies.txt' file.
        :param dict scaff_ass_pairs: A dictionary containing scaffold IDs as keys and assembly IDs as values.
        :rtype list: A list containing the dereplicated scaffolds
        """  
        
        print("--- STEP 5: Extracting cleaned scaffold IDs ---")
        # Initiate empty list to store the result:
        dereplicated_scaffolds = []
        
        #TODO try-except clause here when opening the file?
        # The dereplicated genome assembly IDs are now stored in the file "dereplicated_assemblies.txt" on separate lines
        # We read this file and couple the assembly IDs back to their scaffold IDs through the dictionary.
        with open(path_to_dereplicated_assemblies, 'r') as file:
                for line in file:
                        # Get the assembly accession out of the file name:
                        # file names for the genomes are structured as follows: [GCF][_][nine digits][.][version][_][ASMblabla][.fna]
                        # We only need the part before the second "_" character. Hence, we split the string two times based on the "_" character
                        # and then join the first two parts back toegther. This results in the assembly ID as it is found in the dictionary
                        assembly = "_".join(line.split("_", 2)[:2])
                        
                        # Now we extract the corresponding scaffold ID. The index of the assembly ID in the dictionary is calculated.
                        # Then, this index is passed onto the list of keys to get the element at that specific index.
                        scaffold = list(scaff_ass_pairs.keys())[list(scaff_ass_pairs.values()).index(assembly.strip())]
                        dereplicated_scaffolds.append(scaffold)

                        # Now we have to match this scaffold ID with the corresponding hit in the original cblaster binary file.
                        #subprocess.run(f"grep '{scaffold}' {path_to_binary} >> ../../data/output/cleaned_binary.csv", shell=True, check=True)
        
        print(f"Got {len(dereplicated_scaffolds)} scaffold IDs")
        return dereplicated_scaffolds


def write_output(dereplicated_scaffolds:list, path_to_summary:str, path_to_binary:str) -> None:
        """
        This function takes a list of scaffold IDs and the cblaster summary and binary output files.
        It writes the corresponding Cluster IDs and cleaned hits to a file.

        :param list dereplicated_scaffolds: A list containing the cleaned scaffold IDs.
        :param str path_to_summary: Path to summary file.
        :param str path_to_binary: Path to binary file.
        """
       
        print("--- STEP 6: Generating output files. ---")
        
        # Create intermediary list to store the cleaned hits:
        cleaned_hits = []

        # First we generate the cleaned binary output file:
        with open(path_to_binary, 'r') as file:
            file_content = file.read()
            for scaffold in dereplicated_scaffolds:
                pattern = f".*{scaffold}.*"
                cleaned_hits.append(re.search(pattern, file_content).group(0))
        
        with open('../../data/output/cleaned_binary.csv', 'a') as file:
            for cleaned_hit in cleaned_hits:
                file.write(f"{cleaned_hit}\n")

        # Create an intermediary list to store the cluster IDs:
        clusters = []
        
        # Open the summary file in read mode:
        with open(path_to_summary, 'r') as file:
            
            # Read the file contents:
            file_content = file.read()
            
            # Loop over each scaffold ID
            for scaffold in dereplicated_scaffolds:
                
                # Construct the regex pattern to be matched in the file contents:
                # This regex has the cluster line as the second group. First the scaffold is matched, followed by a new line 
                # containing the "-" character 11 times. Th next line is the one containing the cluster ID
                pattern = f"({scaffold}\n[-]*\n)(Cluster \\d*)"
    
                # Search for the pattern and append the result to the list:
                clusters.append(re.search(pattern, file_content).group(2))
        
        # Now we open a file called clusters.txt and write the contents of the list to it:
        with open('../../data/output/clusters.txt', 'a') as file:
            for cluster in clusters:
                file.write(f"{cluster}\n")

        print("All done! Results are written to /data/output")    


def main():
        
        # Assert correct amount of arguments are given:
        assert len(sys.argv) == 4, "Incorrect amount of arguments.\nUsage: python3 main.py <path_to_binary> <path_to_summary> <percent_identity_cutoff>"        

        # Path to the cblatser binary output file, given as command line argument       
        path_to_binary = sys.argv[1]  # path to binary file
        path_to_summary = sys.argv[2]  # path to summary file
        pi_cutoff = float(sys.argv[3])  # percent identity cutoff for skDER, must be float
         
        if validate_input_files(path_to_binary, path_to_summary):
            # First we capture the dictionary with scaffold:assembly pairs:
            scaff_ass_pairs = get_assemblies(get_scaffolds(path_to_binary))

            # Then we download and cluster the genomes through the helper script:
            # This generates a file called "dereplicated_assemblies.txt"
            dereplicate_genomes(list(scaff_ass_pairs.values()), pi_cutoff)
            
            # Check if the file 'dereplicated_assemblies.txt' exists:
            if not os.path.isfile("dereplicated_assemblies.txt"):
                raise FileNotFoundError("Could not find the intermediary output file 'dereplicated_assemblies.txt'.")

            # Final step: Couple dereplicated assembly IDs back to their respective hits in the binary and summary file:
            write_output(get_dereplicated_scaffolds('dereplicated_assemblies.txt', scaff_ass_pairs), path_to_summary, path_to_binary)
        
        else:
            print("Validation failed. Exiting the program.")
            sys.exit()

if __name__ == "__main__":
        main()
        
