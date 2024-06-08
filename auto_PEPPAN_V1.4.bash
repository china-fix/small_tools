#!/bin/bash

# Remove existing directories if they exist
rm -rf prokka_ann
rm -rf PEPPAN_working
rm -rf wash_input

# Get the starting directory
start_dir=$(pwd)

# Create necessary directories
mkdir -p prokka_ann
mkdir -p PEPPAN_working/gffs
mkdir -p wash_input

# Rename contig names in the fasta files in the "input" directory
ls input/ | parallel --verbose -j 40 "awk '/^>/ {print \">contig_\" ++i; next} {print}' input/{} > wash_input/{}"

# Run Prokka on the fasta files in the "wash_input" directory
ls wash_input/ | parallel --verbose -j 40 "prokka --prefix {} --compliant --rfam --cpus 6 --outdir $start_dir/prokka_ann/{} wash_input/{}"

# Copy GFF files to the "PEPPAN_working/gffs" directory
ls wash_input/ | parallel --verbose -j 40 "cp $start_dir/prokka_ann/{}/{}.gff $start_dir/PEPPAN_working/gffs/"

# Change directory to "PEPPAN_working"
cd $start_dir/PEPPAN_working

# Run PEPPAN on the GFF files
PEPPAN -p OUT $start_dir/PEPPAN_working/gffs/*.gff -t 50

# Run PEPPAN_parser on the output GFF file
PEPPAN_parser -g $start_dir/PEPPAN_working/OUT.PEPPAN.gff -p PAR_OUT -t -a 98 -c
