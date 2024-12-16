"""
Usage: `python3 main.py <path_to_binary.csv> <path_to_summary.txt> <percent identity cutoff for dereplication>`

This python script is designed to clean up redundant hits from the cblaster tool.
It takes the binary and summary output files as arguments, 
along with a percent identity cutoff to dereplicate the genomes using the skDER tool.

After dereplication of the genomes from which the original hits came, two files are generated as output:
    - cleaned_binary.csv: a file structured in the same way as the cblaster binary output, except that redundant hits have been removed. 
    - clusters.txt: the corresponding cluster IDs from the summary file for each cleaned hit.

The code is written in a chronological way. 
Readers can go through the functions top to bottom, the order in which they will be executed by main().
"""

import polars as pl  # To read csv file
import sys  # For cli argument parsing
import subprocess  # To execute shell commands
import re  # To capture patterns
import os  # For validating file paths
import time  # To measure the time
from progress.bar import Bar

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
        
        # Check if the required column names are present (Organism, Scaffold...):
        assert {"Organism", "Scaffold", "Start", "End", "Score"} < set(data.columns), "The column names do not correspond to a cblaster binary output file. Please check the file structure.\n The columns 'Organims', 'Scaffold', 'Start', 'End' and 'Score' should be present."

        # Check if there are at least two hits (3 because of header):
        assert data.height > 3, "We are expecting at least two hits in the cblaster binary output file."
        
        print("CHECK") 
        
        print("Validating summary file.")

        # Again, check if the summary file exists:
        if not os.path.isfile(path_to_summary):
            raise FileNotFoundError(f"'{path_to_summary}' does not exist. Please check if the file path is correct and try again.")
        
        # The amount of appearances of the word "Cluster" in the summary file should correspond to the amount of lines
        # in the binary file:
        with open (path_to_summary, 'r') as file:
            assert len(re.findall("Cluster", file.read())) == data.height, "Binary and summary file do no match." 
        
        print("CHECK")
        
        return True

def get_stats(path_to_summary: str) -> None:
    """
    This function reads the cblaster summary file and provides some descriptive statistics.
    More specifically: What species are present, and for each species; how many strains are there?

    :param str path_to_summary: path to the cblaster summary file
    """ 
    # Initiate empty dictionary to store 'species:amount of strains' key-value pairs. 
    result = {}
    
    # This is the regex pattern to be acptured in the summary file. Breakdown:
    # ([A_Z][a-z]+ ): GENUS, any string starting with capital letter and only containing lowercase alphabetical characters afterwards. Ends with space.
    # ([a-z.]+): SPECIES, any string with only lowercase alphabetical characters or a point. Sometimes names can be like Staphylococcus sp.
    # (.*\n): STRAIN, anything that comes afterwards and ends with a new line
    # (=+: The sequence of equal characters that's on the line beneath)
    pattern = "([A-Z][a-z]+ )([a-z.]+)(.*\n)(=+)"
    
    # Open the file and find all matches. Matches are returned in a list of tuples.
    with open(path_to_summary, 'r') as file:
        matches = re.findall(pattern, file.read())
    
    # Each match is now in a tuple of strings where the first string correpsonds to the first capture group in the pattern, 
    # second string to second capture group etc.
    # For example: [('Staphylococcus', 'aureus', 'HKU19'), (...), (...)]
    # We loop over the matches and extract the relevant info
    for hit in matches:
        organism_name = hit[0] + hit[1]  # genus + species
        strain = hit[2]
        if organism_name not in list(result.keys()):  # if 'genus species' is not in the dict yet, add it as a key
            result[organism_name] = []  # The value will be a list of strains
            result[organism_name].append(strain)  # Add the first strain
                
        else:  # If the species is already in the dictionary, check if the strain is already in the values. If not, append it to the value (which is a list).
            if strain not in list(result[organism_name]):
                result[organism_name].append(strain)
            else:
                continue
    
    # Print out the results:
    print("\nDescriptive statistics of the cblaster summary file:")
    print(f"Amount of species: {len(list(result.keys()))}")
    for species in result:
        print(f"{species}: {len(result[species])} strains")
    

def get_scaffolds(path_to_binary: str) -> list:
        """
        This function extracts the scaffold IDs from the cblaster binary output file using polars.
        
        :param str path_to_binary: Path to the cblaster binary output file.
        :rtype list: A list containing all the scaffold IDs.
        """
        print("\n--- STEP 2: Extracting scaffold IDs from the binary input file. ---")

        # Read the file as comma separated, extract the second column (containing the scaffold IDs) as a polars series 
        # and convert to a python list. truncate_ragged_lines took care of reading errors.
        # One could also do this without polars.
        scaffolds = pl.read_csv(path_to_binary, truncate_ragged_lines=True).to_series(1).to_list()

        print(f"Extracted {len(scaffolds)} scaffold IDs")
        
        # Return the list:
        return scaffolds
        

def get_assemblies(scaffolds: list) -> dict:
        """
        This function obtains the genome assembly ID for each scaffold ID obtained by get_scaffolds().
        The subprocess module is used to execute a piped command using NCBI e-utilities.

        :param list scaffolds: A list containing scaffold IDs.
        :rtype dict: A dictionary with scaffold IDs as keys and assembly IDs as values.
        """        
        # Set begin time:
        begin = time.time()
    
        print("\n--- STEP 3: Contacting the NCBI servers and retreiving genome assembly IDs for each scaffold. This may take a while. ---")
                
        # Create empty dictionary to store genome assembly IDs with scaffold IDs as keys:
        result = {}

        # A counter to count duplicate assembly IDs:
        counter = 0
        
        # Initiate a progress bar:
        bar = Bar('Fetching assembly IDs:', max=len(scaffolds))
    
        # Loop over each scaffold ID:
        for i in range(0, len(scaffolds)):
                # BREAKDOWN:
                # esummary: get the scaffold ID entry information in xml
                # xtract: extract the BioSample field in the xml. This can then be coupled to the genome database
                # elink: link the previously obtained BioSample to the 'assembly' database.
                # efetch: fetch the information in xml (DocumentSummary) format
                # xtract: extract the AssemblyAccession field
                # TODO Check if the scaffold accession ID exists in the NCBI database?
                # TODO Properly handle errors. Possibly using try-except. `process.stderr` is av for catching content of errors.
                # TODO Make it faster (suggestion: using multiprocess module and NCBI API key)

                process = subprocess.run(f"esummary -db nucleotide -id {scaffolds[i]} | xtract -pattern DocumentSummary -element BioSample | elink -db biosample -target assembly | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession", shell=True, capture_output=True, text=True, check=True)    
         
                # Check if the output is indeed something that starts with "GC[F,A]_". If not, it is not a valid assembly ID:
                # TODO Properly handle this scenario. What to do if it fails?
                if process.stdout.strip().startswith("GC"):
                    # Check if the assembly ID is a duplicate. If yes, duplicate couunter + 1. If not, add it to the dictionary.
                    if process.stdout.strip() in result.values():
                        counter += 1
                        continue
                    else:
                        result[scaffolds[i]] = process.stdout.strip()
                else:
                    print(f"\nAn error occurred in retreiving the assembly ID for {scaffolds[i]}.")
                    # raise ValueError("message"), but we don't want the program to stop completely for one wrong assembly ID?
                    # Perhaps write error messages to a file. Mostly it's just HTTP error 400: bad request
                    continue
                # Update the progress bar:
                bar.next()
        
        # Close the progress bar:
        bar.finish()
                
        # Print how many assembly IDs were extracted. This should equal the amount of scaffold IDs
        print(f"Obtained {len(result.values())} unique assembly IDs. {counter} were duplicates") 

        # Set end time:
        end = time.time()
        
        print(f"Time elapsed: {round(end - begin, 2)} s")

        return result
        

def dereplicate_genomes(assemblies: list, pi_cutoff:float=99.0) -> None:
        """
        This function takes a list of genome assembly IDs and calls a helper bash script that downloads and dereplicates the genomes.
        The standard ANI cutoff for dereplication is 99.
        
        First, the assembly IDs are written to a file 'assemblies.txt' in batches of 300
        to be passed as an argument to a helper bash script.
        
        A file called "dereplicated_assemblies.txt" is then generated by the helper script.

        :param list assemblies: A list of assembly IDs
        :param float pi_cutoff: The percent identity cutoff for dereplicating genomes (see skDER docs)
        """
        # Assert that pi cutoff is valid:
        assert 0 < pi_cutoff <= 100, "Percent identity cutoff must be between 0 and 100. E.g.: 99.0"
    
        # Set begin time:
        begin = time.time()

        print("\n--- STEP 4: Downloading and dereplicating genomes. ---")

        # First we write the assembly IDs to a file that can be used by the helper bash script (helper.sh):
        # Easiest is to write assembly IDs space separated on one line. 
        # F.e.: GCA_1 GCA_2 GCF_3 ... GCA_300 \n 
        with open('assemblies.txt', 'a') as file:
                for i in range(0, len(assemblies), 300):  # We download the genomes in batches of 300
                    # If this chunk of the list does not exceed a length of 300, we can safely write 300 IDs on one line:
                    if len(assemblies[i:]) > 300:
                        file.write(' '.join(assemblies[i:i+300]))
                    # If the length is smaller than 300, just go to the end of the list:
                    else:
                        file.write(' '.join(assemblies[i:]))
                    # Write a new line in any case:
                    file.write('\n') 

        # Run the bash script to download and cluster genomes:
        # First in the list is the command, the rest are considered as arguments. check=True checks for errors.
        subprocess.run(["bash", "helper.sh", str(pi_cutoff)], check=True)
         
        #TODO CHECK OUTPUT VALIDITY CATCH ERRORS
       
        # Check if the intermediary file 'dereplicated_assemblies.txt' exists:
        if not os.path.isfile("../../data/output/dereplicated_assemblies.txt"):
            raise FileNotFoundError("Could not find the intermediary output file 'dereplicated_assemblies.txt'.")

        # Set end time:
        end = time.time()
    
        print(f"Time elapsed: {round(end - begin, 2)} s")
        
        return None

 
def get_dereplicated_scaffolds(path_to_dereplicated_assemblies: str, scaff_ass_pairs:dict) -> list:
        """
        This function reads the 'dereplicated_assemblies.txt' file generated by dereplicate_genomes().
        Assembly IDs are written on seperate lines. One ID per line.
        Each assembly ID is coupled to the corresponding scaffold ID via a previously generated dictionary.

        :param str dereplicated_assemblies: File path to the 'dereplicated_assemblies.txt' file.
        :param dict scaff_ass_pairs: A dictionary containing scaffold IDs as keys and assembly IDs as values. (generated by get_assemblies())
        :rtype list: A list containing the dereplicated scaffolds
        """  
        print("\n--- STEP 5: Extracting cleaned scaffold IDs ---")
        # Initiate empty list to store the result:
        dereplicated_scaffolds = []
        
        # TODO try-except clause here when opening the file?
        # The dereplicated genome assembly IDs are now stored in the file "dereplicated_assemblies.txt" on separate lines
        # We read this file and couple the assembly IDs back to their scaffold IDs through the dictionary.
        with open(path_to_dereplicated_assemblies, 'r') as file:
                for line in file:
                        # Get the assembly accession out of the file name:
                        # file names for the genomes (assembly IDs) are structured as follows: [GCF][_][nine digits][.][version][_][ASMblabla][.fna]
                        # We only need the part before the second "_" character. Hence, we split the string two times based on the "_" character
                        # and then join the first two parts back together. This results in the assembly ID as it is found in the dictionary
                        assembly = "_".join(line.split("_", 2)[:2])
                        
                        # Now we extract the corresponding scaffold ID. The index of the assembly ID in the dictionary is calculated.
                        # Then, this index is passed onto the list of keys to get the element at that specific index.
                        scaffold = list(scaff_ass_pairs.keys())[list(scaff_ass_pairs.values()).index(assembly.strip())]
                        dereplicated_scaffolds.append(scaffold)

        print(f"Got {len(dereplicated_scaffolds)} cleaned scaffold IDs")
        return dereplicated_scaffolds


def write_output(dereplicated_scaffolds:list, path_to_summary:str, path_to_binary:str) -> None:
        """
        This function takes a list of (cleaned) scaffold IDs and the cblaster summary and binary output files.
        It writes the corresponding Cluster IDs and cleaned hits to a file.

        :param list dereplicated_scaffolds: A list containing the cleaned scaffold IDs.
        :param str path_to_summary: Path to summary file.
        :param str path_to_binary: Path to binary file.
        """
        print("\n--- STEP 6: Generating output files. ---")
        
        ## First we do the binary file:

        # Create intermediary list to store the cleaned hits:
        cleaned_hits = []

        # Open the binary file
        with open(path_to_binary, 'r',) as file:
            # Read the file contents
            file_content = file.read()
            # We also capture the header to write to our cleaned file later:
            header = file_content.split("\n")[0] + "\n" 

            # Loop over the cleaned scaffold IDs and match them in the binary file using regex.
            for scaffold in dereplicated_scaffolds:
                pattern = f".*{scaffold}.*"
                # Append to the list of cleaned hits:
                cleaned_hits.append(re.search(pattern, file_content).group(0))
            
        # Now we generate a file to write the cleaned hits to:
        with open('../../data/output/cleaned_binary.csv', 'a') as file:
            # The header line comes first:
            file.write(header)
            for cleaned_hit in cleaned_hits:
                file.write(f"{cleaned_hit}\n")

        ## Secondly, the cluster IDs. Same principle but then with the summary file:

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
        """
        This function executes all the previous functions in a logical manner.
        """
        # Set begin time:
        begin = time.time()

        # Assert correct amount of arguments are given:
        assert len(sys.argv) == 4, "Incorrect amount of arguments.\nUsage: python3 main.py <path_to_binary> <path_to_summary> <percent_identity_cutoff>"        
        
        # Read the arguments given at the command line:       
        path_to_binary = sys.argv[1]  # path to binary file
        path_to_summary = sys.argv[2]  # path to summary file
        pi_cutoff = float(sys.argv[3])  # percent identity cutoff for skDER, must be converted to float
         
        # Validate the input files. If it returns True, we execute the program.
        if validate_input_files(path_to_binary, path_to_summary):
            # Spit out some descriptive statistics of the summary file:
            # This function is commented out because cblaster already provides similar stats, and NCBI sometimes
            # gives us organism names such as: 'No organism staphylococcus'. The function cannot handle this properly.
            #get_stats(path_to_summary)
            
            # Capture the dictionary with scaffold:assembly pairs:
            scaff_ass_pairs = get_assemblies(get_scaffolds(path_to_binary))

            # Then we download and cluster the genomes through the helper script:
            # This generates a file called "dereplicated_assemblies.txt"
            # This function expects a list of assembly IDs to download, which can be found in the above dictionary as values
            dereplicate_genomes(list(scaff_ass_pairs.values()), pi_cutoff)
            
            # Final step: Couple dereplicated assembly IDs back to their respective hits in the binary and summary file:
            write_output(get_dereplicated_scaffolds('../../data/output/dereplicated_assemblies.txt', scaff_ass_pairs), path_to_summary, path_to_binary)
                
        # If the validation fails, we exit the program:
        else:
            print("Validation failed. Exiting the program.")
            sys.exit()

        # Set end time:
        end = time.time()
        
        print(f"Total runtime: {round(end - begin, 2)} s")

if __name__ == "__main__":
        main()
        
