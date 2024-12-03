#!/usr/bin/bash


# cblaster binary output file should be given as first argument.
file_path=$1

echo Converting $1 to .csv format.

# First convert multiple spaces to tab, then replace tab with comma, then write to file where .txt is replaced by .csv
sed 's/ \+ /\t/g' $1 | tr -s "\t" "," > ${file_path/.txt/.csv}

echo Done! Your new file is named ${file_path/.txt/.csv}. 
