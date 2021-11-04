#! /bin/bash
input_file=$1
output_file=$2

# if input and output file exists
if [[ -f "$input_file" && -f "$output_file" ]]
then
	# read it line by line
	while IFS= read -r line || [[ -n "$line" ]] 
	do
		url="https://www.rcsb.org/fasta/entry/$line"
		wget $url -nv -O - >> $output_file
	done < "$input_file"
fi
