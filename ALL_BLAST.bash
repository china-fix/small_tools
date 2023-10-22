#!/bin/bash

# Initialize variables
output_file=""
genome_files=""
query_file=""

# Function to display program usage
usage() {
    echo "Welcome to use Xiao's in silico robot"
    echo "Usage: $0 -o <output_file> -i <genome_files> -q <query_file>"
    echo "  -o <output_file>: Output file for BLAST results."
    echo "  -i <genome_files>: Path to genome files in .gz format (e.g., '/path/to/genome_folder/*.gz')."
    echo "  -q <query_file>: Path to query sequence file."
    exit 1
}

while getopts "ho:i:q:" opt; do
    case $opt in
        h)
            usage
            ;;
        o)
            output_file="$OPTARG"
            ;;
        i)
            genome_files="$OPTARG"
            ;;
        q)
            query_file="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done

# Check if required parameters are provided
if [ -z "$output_file" ] || [ -z "$genome_files" ] || [ -z "$query_file" ]; then
    echo "Missing required options. See usage with -h option."
    exit 1
fi

# Initialize the output file with a header
echo -e "Query_ID\tSubject_ID\t%_Identity\tAlignment_Length\tQuery_Length\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE-Value\tBit_Score\tSubject_File_Name\tCoverage(%)" > "$output_file"


# Loop through the genome files and append results to the output file
for genome_file in $genome_files; do
    subject_file_name="$(basename "$genome_file")"
    zcat "$genome_file" | blastn -query "$query_file" -subject - -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" | awk -v subject_file_name="$subject_file_name"  '{coverage = ($4 / $5) * 100; print $0 "\t" subject_file_name "\t" coverage}' >> "$output_file"
done

