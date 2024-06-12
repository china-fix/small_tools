#!/bin/bash

# Set starting directory and necessary paths
start_dir=$(pwd)
gffs_path="$start_dir/PEPPAN_working/gffs"
prokka_ann_out="$start_dir/prokka_ann"

# Function to rename contig names in input files
rename_contigs() {
    local file=$1
    local output_file=${file%.fasta}_renamed.fasta

    awk '/^>/ {print \">contig_\" ++i; next} {print}' "$file" > "$output_file"
}

# Function to run Prokka on input files
run_prokka() {
    local file=$1
    local out_prefix=${file%.*}
    local out_dir="$prokka_ann_out/$out_prefix"

    prokka --prefix $out_prefix --compliant --rfam --cpus 6 --outdir "$out_dir" "$start_dir/input/$file"
}

# Function to run PEPPAN on input files
run_peppan() {
    local file=$1

    # Run PEPPAN with appropriate parameters and time limit
    PEPPAN -p OUT "$gffs_path/$file" -t 50 &
    
    # Wait for the child process (PEPPAN) to complete before continuing
    wait $!
}

# Function to run PEPPAN_parser on input file
run_peppan_parser() {
    local gff_file=$1
    local out_gff="${gff_file%.gff}_parsed.gff"
    local out_path="$start_dir/PEPPAN_working/PAR_OUT"

    PEPPAN_parser -g "$gff_file" -p "$out_path" -t -a 98 -c > "$out_gff"
}

# Rename contig names in input files and store results
rename_contigs input/*
parallel --jobs 40 "cat {} | rename_contigs.sh" wash_input/*.fasta

# Run Prokka on renamed files
run_prokka $(ls wash_input/*_renamed.fasta)

# Move GFF files to PEPPAN_working/gffs directory and run PEPPAN
mv wash_input/*.gff "$gffs_path"
for file in $gffs_path/*.gff; do
    run_peppan "$file"
done

# Run PEPPAN_parser on output files from previous step
run_peppan_parser $(ls "$gffs_path"/*)
