import argparse
import time
from Bio import Entrez, SeqIO
import urllib.error
import os

def download_and_save_protein_sequences(gene_ids, output_filename, delay):
    # Provide your email address to the NCBI (required)
    Entrez.email = "robot_ncbi@edu.com"
    
    # Truncate the output file to ensure it's empty at the beginning
    with open(output_filename, 'w') as output_file:
        output_file.truncate(0)
        
    sequences = []

    for gene_id in gene_ids:
        # Search for the protein records using the Gene ID and organism
        term = "{}[uid]".format(gene_id)
        while True:
            try:
                gene_handle = Entrez.esearch(db="gene", term=term)
                break
            except  urllib.error.HTTPError:
                print("HTTPError encountered during esearch. Waiting for {} seconds before retrying.".format(delay))
                time.sleep(delay)

        record = Entrez.read(gene_handle)

        if int(record["Count"]) == 0:
            print("Gene  {} not found.".format(gene_id))
        else:
            print("hey xiao, i am extrating related proteins from Gene {} now!".format(gene_id))
            gene_ids = record["IdList"]

            for e_gene_id in gene_ids:
                # Retrieve the related protein records
                while True:
                    try:
                        # Use efetch to retrieve the related protein IDs
                        protein_handle = Entrez.elink(dbfrom="gene", db="protein", id=e_gene_id)
                        protein_record = Entrez.read(protein_handle)
                        if len(protein_record) > 0 and "LinkSetDb" in protein_record[0] and len(protein_record[0]["LinkSetDb"]) > 0:
                            protein_ids = [link["Id"] for link in protein_record[0]["LinkSetDb"][0]["Link"]]
                        else:
                            print("No protein records found for Gene ID {}.".format(gene_id))
                            break  # Exit the loop if there are no protein records
                        break  # Exit the retry loop if successful
                        
                    except urllib.error.HTTPError:
                        print("HTTPError encountered during elink. Waiting for {} seconds before retrying.".format(delay))
                        time.sleep(delay)
                    except RuntimeError:
                        print("RuntimeError encountered during elink. Waiting for {} seconds before retrying.".format(delay))
                        time.sleep(delay)

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
    parser.add_argument("--list", help="File containing a list of gene ids (one per line)", required=True)
    parser.add_argument("--output", help="Output filename for the FASTA file (default: protein_sequences.fasta)", default="protein_sequences.fasta")
    parser.add_argument("--delay", type=float, help="Delay between requests in seconds (default: 1)", default=1)
    
    args = parser.parse_args()

    # Read the list of gene IDs from the input file
    with open(args.list, "r") as file:
        gene_ids = [line.strip() for line in file]

    # Download and save protein sequences with user-specified delay
    download_and_save_protein_sequences(gene_ids, args.output, args.delay)
    print("Protein sequences saved to {}".format(args.output))

if __name__ == "__main__":
    main()

