import os
import subprocess
import argparse
from Bio import SeqIO

def read_fasta(filename):
    """Read FASTA file and return a dictionary of sequences."""
    sequences = {}
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences[record.id] = str(record.seq)
    return sequences

def run_blastn(query_sequence, subject_sequences, output_file, perc_identity_threshold=80, qcovs_threshold=80):
    """Run BLASTn to search for query sequence in subject sequences."""
    # Define the required fields for BLASTn output format (-outfmt 6)
    blastn_fields = "qseqid sseqid length pident qcovs sseq"
    # Construct the BLASTn command with the specified output format and filtering parameters
    with open(output_file, 'w') as outfile:
        blastn_command = (
            f"blastn -query {query_sequence} -subject {subject_sequences} "
            f"-out {outfile.name} -outfmt '6 {blastn_fields}' "
            f"-perc_identity {perc_identity_threshold} -qcov_hsp_perc {qcovs_threshold} "
            "-max_target_seqs 1"
        )
        # Run the BLASTn command using subprocess
        subprocess.run(blastn_command, shell=True, check=True)

def is_pseudo_sequence(aligned_sequence):
    """Check if the aligned sequence is a pseudo-sequence."""
    premature_stop_codons = {"TAA", "TAG", "TGA"}

    # Check for frameshifts (e.g., insertion or deletion of nucleotides)
    if len(aligned_sequence) % 3 != 0:
        return 'frameshifts'
    
    # Check the last three codons for stop codons
    last_three_codons = aligned_sequence[-3:].upper()
    if last_three_codons not in premature_stop_codons:
        return 'no_stop_condon'
    
    # Iterate through the aligned sequence in steps of 3 to extract codons
    for i in range(0, len(aligned_sequence) - 3, 3):
        codon = aligned_sequence[i:i+3].upper()
        if codon in premature_stop_codons:
            return 'premature'

    # Add more checks here as needed...

    return 'PASS'



def parse_blastn_results(blastn_output):
    """Parse BLASTn output and return a dictionary of hits."""
    hits = {}
    with open(blastn_output, 'r') as file:
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
    parser = argparse.ArgumentParser(description="Check if CDS sequences are present in assembly contigs using BLASTn.")
    parser.add_argument("--assembly", "-a", type=str, required=True, help="Path to assembly contigs FASTA file.")
    parser.add_argument("--cds", "-c", type=str, required=True, help="Path to CDS sequences FASTA file.")
    parser.add_argument("--perc_identity", type=float, default=80, help="Minimum percentage identity threshold (default: 80).")
    parser.add_argument("--qcovs", type=float, default=80, help="Minimum query coverage threshold (default: 80).")
    args = parser.parse_args()

    # Read assembly contigs and CDS sequences
    assembly_contigs = read_fasta(args.assembly)
    cds_sequences = read_fasta(args.cds)

    # Write CDS sequences to a temporary file
    cds_temp_file = "cds_temp.fasta"
    with open(cds_temp_file, 'w') as f:
        for cds_id, cds_seq in cds_sequences.items():
            f.write(f">{cds_id}\n{cds_seq}\n")

    # Run BLASTn to search for CDS sequences in assembly contigs
    blastn_output_file = "blastn_results.txt"
    run_blastn(cds_temp_file, args.assembly, blastn_output_file, args.perc_identity, args.qcovs)

    # Parse BLASTn results
    blastn_results = parse_blastn_results(blastn_output_file)

    # Analyze BLASTn results and print table
    print("CDS_ID\tSituation\tIdentity\tCoverage\tIntegrity")
    for cds_id, cds_seq in cds_sequences.items():
        if cds_id in blastn_results:
            subject_id, length, identity, qcovs, aligned_sequence, is_pseudo = blastn_results[cds_id]
            print(f"{cds_id}\tPresent\t{identity:.2f}\t{qcovs:.2f}\t{is_pseudo}")
        else:
            print(f"{cds_id}\tAbsent")

    # Remove temporary files
    os.remove(cds_temp_file)
    os.remove(blastn_output_file)

if __name__ == "__main__":
    main()
