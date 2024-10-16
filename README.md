This is the readme.

# What is this repo?
For the Integrated Bioinformatics Project, we have been assigned the task of optimizing the CAGECat tool.


# Problem definition
TODO
CompArative GEne Cluster Analysis Toolbox (CAGECAT) is an online web server used to find novel biosynthetic gene clusters 
(BGCs). In the back, CAGECAT relies on [cblaster](https://github.com/gamcil/cblaster) and [clinker](https://github.com/gamcil/clinker).
The current issue in this workflow is that the output consists of many duplicates that we want to filter out. When we provide
a query (a set of genes of which we believe are clustered and involved in seondary metabolite production), cblaster uses
the NCBI API to BLAST these sequences and find homologues in the database. However, often times, many different strains of
the same organism are present in the NCBI database, which all have the same BGC, and make their way to the output. 
This is inconvenient, since we only want to find novel BGCs, and not have to manually select one strain of a specific species
......... 


# Time schedule
IMPORT GANTT

# TODO
IMPORT TODO FILE
- [] test
	- [x] test
