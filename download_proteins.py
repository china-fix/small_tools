import argparse
import time
from Bio import Entrez, SeqIO
import urllib.error
import os

def download_and_save_protein_sequences(gene_ids, organism, output_filename, delay):
    # Provide your email address to the NCBI (required)
    Entrez.email = "robot_ncbi@edu.com"
    
    # Truncate the output file to ensure it's empty at the beginning
    with open(output_filename, 'w') as output_file:
        output_file.truncate(0)
        
    sequences = []

    for gene_id in gene_ids:
        # Search for the protein records using the Gene ID and organism
        term = "{}[Gene ID] AND {}[Organism]".format(gene_id, organism)
        while True:
            try:
                handle = Entrez.esearch(db="protein", term=term)
                break
            except  urllib.error.HTTPError:
                print("HTTPError encountered during esearch. Waiting for {} seconds before retrying.".format(delay))
                time.sleep(delay)

        record = Entrez.read(handle)

        if int(record["Count"]) == 0:
            print("Gene ID {} not found for organism {}.".format(gene_id, organism))
        else:
            print("hey xiao, i am extrating Gene ID {} from organism {} now!".format(gene_id, organism))
            protein_ids = record["IdList"]

            for protein_id in protein_ids:
                # Retrieve the protein sequence
                while True:
                    try:
                        # Retrieve the protein sequence
                        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                        record = SeqIO.read(handle, "fasta")
                        sequences.append(record)
                        break  # Exit the retry loop if successful
                    except urllib.error.HTTPError:
                        print("HTTPError encountered during efetch. Waiting for {} seconds before retrying.".format(delay))
                        time.sleep(delay)

                # Add a user-specified delay to control the download frequency
                time.sleep(delay)

        # Save the sequences to a FASTA file
        SeqIO.write(sequences, output_filename, "fasta")
        print("Protein sequences for Gene ID {} saved to {}".format(gene_id, output_filename))

def main():
    parser = argparse.ArgumentParser(description="Download protein sequences based on a list of gene IDs.")
    parser.add_argument("--list", help="File containing a list of gene IDs (one per line)", required=True)
    parser.add_argument("--organism", help="Organism name (e.g., 'Homo sapiens') for Entrez search", required=True)
    parser.add_argument("--output", help="Output filename for the FASTA file (default: protein_sequences.fasta)", default="protein_sequences.fasta")
    parser.add_argument("--delay", type=float, help="Delay between requests in seconds (default: 1)", default=1)
    
    args = parser.parse_args()

    # Read the list of gene IDs from the input file
    with open(args.list, "r") as file:
        gene_ids = [line.strip() for line in file]

    # Download and save protein sequences with user-specified delay
    download_and_save_protein_sequences(gene_ids, args.organism, args.output, args.delay)
    print("Protein sequences saved to {}".format(args.output))

if __name__ == "__main__":
    main()

