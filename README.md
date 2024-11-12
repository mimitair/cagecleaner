# What is this repo?
For the Integrated Bioinformatics Project, we have been assigned the task of optimizing the CAGECat tool.
The INDEX.txt file explains the file hierarchy for this project.

# Problem definition
CompArative GEne Cluster Analysis Toolbox (CAGECAT) is an online web server used to find novel Biosynthetic Gene Clusters 
(BGCs). In the back, CAGECAT relies on [cblaster](https://github.com/gamcil/cblaster) and [clinker](https://github.com/gamcil/clinker).

The current issue in this workflow is that the output consists of many duplicate hits (see figure below). 

<img src="text/figures/example_output.png" alt="Example output of cblaster run." width=500> 

When we provide a query (a set of genes of which we believe are clustered and involved in seondary metabolite production), cblaster uses
the NCBI API to BLAST these sequences and find similar ones. cblaster then uses the Identical Protein Group database to
fetch genomic coordinates for our genes of interest. If the genes are close enough to each other in the genome, they are considered as a hit. 

However, due to the many duplicate sequence and genome entries in the NCBI database, the output of cblaster quickly becomes flooded with BGCs 
essentially coming from the "same" genome. Hence, the output of this tool becomes difficult to interpret and use in further analysis.

# Solution
Our solution consists of using genome dereplication tools to get rid of the NCBI duplicate entries. 
For each hit in the cblaster output, our tool fetches the corresponding genome assembly.
Once all genomes are fetched, skDER clusters the genomes based on Average Nucleotide Identity (ANI) 
and selects representative genomes for each cluster. These representative genomes are then used to extract the 
relevant BGCs.

# TODO
- [x] Get genome assemblies for each hit in the cblaster binary output
	- [ ] Ensure all assemblies exist and are obtained
- [x] Download genome assemblies using NCBI datasets
	- [ ] Test with large files (hpc)
- [ ] Dereplicate genomes
	- [x] skDER
	- [ ] Pangenome approach?
- [ ] Extract BGCs from representative genomes
- [ ] Create Docker engine for the whole workflow? (optional)
- [ ] Integrate it all into one tool (optional)
