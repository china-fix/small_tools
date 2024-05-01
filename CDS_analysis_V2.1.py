import os
import subprocess
import argparse
from Bio import SeqIO

"""
Xiao 2024-04-25 this robot use tblastn to compare the proteins to an assembly fata file, and check the proteins integrity
update 20240426 this robot use blastp to compare the reference proteins to an assembly faa file (suggested use prokka .ffn file and translate to faa by seqkit), and check the proteins integrity 
"""

def read_fasta(filename):
    """Read FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences[record.id] = str(record.seq)
    return sequences

def run_blastp(query_sequence, subject_sequences, output_file, qcovs_threshold=90):
    """Run BLASTn to search for query sequence in subject sequences."""
    # Define the required fields for BLASTp output format (-outfmt 6)
    blastp_fields = "qseqid sseqid length pident qcovs sseq"
    # Construct the BLASTn command with the specified output format and filtering parameters
    with open(output_file, 'w') as outfile:
        blastp_command = (
            f"blastp -query {query_sequence} -subject {subject_sequences} "
            f"-out {outfile.name} -outfmt '6 {blastp_fields}' "
            f"-qcov_hsp_perc {qcovs_threshold} "
            "-subject_besthit"
        )
        # Run the BLASTn command using subprocess
        subprocess.run(blastp_command, shell=True, check=True)

def is_pseudo_sequence(aligned_sequence):
    """Check if the aligned sequence is a pseudo-sequence."""
    internal_stop_codon = "*"
    start_codons = {"M", "V","L","I"}
    stop_codon = "*"

    # Check for internal stop codons
    if aligned_sequence[:-1].upper().find(internal_stop_codon) != -1:
        return 'internal_stop'

    # Check if the sequence starts with the correct start codon
    if not any(aligned_sequence.startswith(start_codon) for start_codon in start_codons):
        return 'incorrect_start_codon'

    # Check if the sequence ends with a stop codon
    if not aligned_sequence.endswith(stop_codon):
        return 'no_stop_codon'

    # If none of the above conditions are met, the sequence is likely valid
    return 'PASS'


    
    # Add more checks here as needed...

    return 'PASS'



def parse_blastp_results(blastp_output):
    """Parse BLASTp output and return a dictionary of hits."""
    hits = {}
    with open(blastp_output, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            query_id = fields[0]
            subject_id = fields[1]
            length = int(fields[2])
            identity = float(fields[3])
            qcovs = float(fields[4])
            aligned_sequence = fields[5]
            is_pseudo = is_pseudo_sequence(aligned_sequence)
            hits[query_id] = (subject_id, length, identity, qcovs, aligned_sequence, is_pseudo)
    return hits

def main():
    parser = argparse.ArgumentParser(description="Check if CDS aa sequences are present in assembly contigs's faa file using BLASTp.")
    parser.add_argument("--assembly", "-a", type=str, required=True, help="Path to assembly contigs's annotated faa file.")
    parser.add_argument("--cds", "-c", type=str, required=True, help="Path to CDS sequences Faa file.")
    # parser.add_argument("--perc_identity", type=float, default=80, help="Minimum percentage identity threshold (default: 80).")
    parser.add_argument("--qcovs", type=float, default=90, help="Minimum query coverage threshold (default: 90).")
    args = parser.parse_args()

    # Read assembly contigs and CDS sequences
    # assembly_contigs = read_fasta(args.assembly)
    cds_sequences = read_fasta(args.cds)

    # Write CDS sequences to a temporary file
    cds_temp_file = "cds_temp.fasta"
    with open(cds_temp_file, 'w') as f:
        for cds_id, cds_seq in cds_sequences.items():
            f.write(f">{cds_id}\n{cds_seq}\n")

    # Run tBLASTn to search for CDS sequences in assembly contigs
    blastp_output_file = "tblastn_results.txt"
    run_blastp(cds_temp_file, args.assembly, blastp_output_file, args.qcovs)

    # Parse BLASTn results
    blastp_results = parse_blastp_results(blastp_output_file)

    # Analyze BLASTn results and print table
    print("protein_ID\tSituation\tIdentity\tCoverage\tIntegrity\tSubject_id\tAlignment_length")
    for cds_id, cds_seq in cds_sequences.items():
        if cds_id in blastp_results:
            subject_id, length, identity, qcovs, aligned_sequence, is_pseudo = blastp_results[cds_id]
            print(f"{cds_id}\tPresent\t{identity:.2f}\t{qcovs:.2f}\t{is_pseudo}\t{subject_id}\t{length}")
        else:
            print(f"{cds_id}\tAbsent")

    # Remove temporary files
    os.remove(cds_temp_file)
    os.remove(blastp_output_file)

if __name__ == "__main__":
    main()
