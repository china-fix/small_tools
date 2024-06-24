import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def extract_sequences(gbk_file, locus_tags_file, output_fasta):
    # Read locus tags into a set for quick lookup
    with open(locus_tags_file) as f:
        locus_tags = set(line.strip() for line in f)
    
    found_locus_tags = set()
    sequences = []
    
    # Parse the GBK file and extract sequences with matching locus tags
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
                if locus_tag in locus_tags:
                    seq_record = SeqRecord(
                        feature.extract(record).seq,
                        id=locus_tag,
                        description=""
                    )
                    sequences.append(seq_record)
                    found_locus_tags.add(locus_tag)
    
    # Write the extracted sequences to a FASTA file
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"Extracted {len(sequences)} sequences to {output_fasta}")
    
    # Find and report locus tags that were not matched
    unmatched_locus_tags = locus_tags - found_locus_tags
    if unmatched_locus_tags:
        print("The following locus tags were not found in the GBK file:")
        for locus_tag in unmatched_locus_tags:
            print(locus_tag)

def main():
    parser = argparse.ArgumentParser(description="Extract sequences from a GBK file based on locus tags and save them to a FASTA file.")
    parser.add_argument('--gbk_file', required=True, help='Path to the GBK file')
    parser.add_argument('--locus_tags_file', required=True, help='Path to the file containing locus tags')
    parser.add_argument('--output_fasta', required=True, help='Path to the output FASTA file')
    
    args = parser.parse_args()
    
    extract_sequences(args.gbk_file, args.locus_tags_file, args.output_fasta)

if __name__ == "__main__":
    main()
