#!/bin/bash

# Function to print usage
usage() {
    echo "Usage: $0 --input_folder INPUT_FOLDER --output_folder OUTPUT_FOLDER"
    echo
    echo "This script processes genomic data using Prokka and PEPPAN tools. It:"
    echo "  1. Renames contig names in the input FASTA files."
    echo "  2. Runs Prokka to annotate the FASTA files."
    echo "  3. Copies the generated GFF files to a working directory."
    echo "  4. Runs PEPPAN on the GFF files."
    echo "  5. Parses the output of PEPPAN."
    echo
    echo "Options:"
    echo "  --input_folder   The directory containing input FASTA files"
    echo "  --output_folder  The directory where output files will be stored"
    echo "  --help           Display this help message"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_folder) input_folder="$2"; shift ;;
        --output_folder) output_folder="$2"; shift ;;
        --help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if input_folder and output_folder are set
if [ -z "$input_folder" ] || [ -z "$output_folder" ]; then
    echo "Error: --input_folder and --output_folder are required."
    usage
fi

# Remove existing directories if they exist
rm -rf "$output_folder/prokka_ann"
rm -rf "$output_folder/PEPPAN_working"
rm -rf "$output_folder/wash_input"

# Get the starting directory
start_dir=$(pwd)

# Create necessary directories
mkdir -p "$output_folder/prokka_ann"
mkdir -p "$output_folder/PEPPAN_working/gffs"
mkdir -p "$output_folder/wash_input"

# Rename contig names in the fasta files in the specified input directory
find "$input_folder" -type f | parallel --verbose -j 40 "awk '/^>/ {print \">contig_\" ++i; next} {print}' {} > $output_folder/wash_input/{/}"

# Run Prokka on the fasta files in the "wash_input" directory
find "$output_folder/wash_input" -type f | parallel --verbose -j 40 "prokka --prefix {/} --compliant --rfam --cpus 6 --outdir $output_folder/prokka_ann/{/} {}"

# Copy GFF files to the "PEPPAN_working/gffs" directory
find "$output_folder/wash_input" -type f | parallel --verbose -j 40 "cp $output_folder/prokka_ann/{/}/{/.}.gff $output_folder/PEPPAN_working/gffs/"

# Run PEPPAN on the GFF files
PEPPAN -p OUT $output_folder/PEPPAN_working/gffs/*.gff -t 50

# Run PEPPAN_parser on the output GFF file
PEPPAN_parser -g $output_folder/PEPPAN_working/OUT.PEPPAN.gff -p $output_folder/PEPPAN_working/PAR_OUT -t -a 98 -c
