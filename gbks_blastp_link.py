import os
import argparse
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Blast import NCBIXML

def extract_cds(gbk_file, output_fasta):
    with open(output_fasta, "w") as out_fasta:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if 'translation' in feature.qualifiers:
                        seq = feature.qualifiers['translation'][0]
                        header = f">{feature.qualifiers['locus_tag'][0]}"
                        out_fasta.write(f"{header}\n{seq}\n")

def run_blastp(query_fasta, db_fasta, output_xml):
    # Create BLAST database
    makeblastdb_cmd = ["makeblastdb", "-in", db_fasta, "-dbtype", "prot"]
    subprocess.run(makeblastdb_cmd, check=True)

    # Run BLASTp
    blastp_cmd = ["blastp", "-query", query_fasta, "-db", db_fasta, "-evalue", "0.001", "-outfmt", "5", "-out", output_xml]
    subprocess.run(blastp_cmd, check=True)

def parse_blastp_output(blastp_xml):
    matches = []
    with open(blastp_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            query_id = blast_record.query
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    match = {
                        "query_id": query_id,
                        "subject_id": alignment.hit_def,
                        "e_value": hsp.expect,
                        "score": hsp.score,
                        "identity": hsp.identities,
                        "alignment_length": hsp.align_length
                    }
                    matches.append(match)
    return matches

def main(args):
    with tempfile.TemporaryDirectory() as temp_dir:
        cds_fasta1 = os.path.join(temp_dir, "genome1_cds.fasta")
        cds_fasta2 = os.path.join(temp_dir, "genome2_cds.fasta")
        blastp_output = os.path.join(temp_dir, "blastp_output.xml")
        
        # Extract CDS sequences
        extract_cds(args.in_gbk1, cds_fasta1)
        extract_cds(args.in_gbk2, cds_fasta2)
        
        # Run BLASTp
        run_blastp(cds_fasta1, cds_fasta2, blastp_output)
        
        # Parse BLASTp output
        matches = parse_blastp_output(blastp_output)
        
        # Output matches to a table
        with open(args.out_file, "w") as out_table:
            out_table.write("Query ID,Subject ID,E-value,Score,Identity,Alignment Length\n")
            for match in matches:
                out_table.write(f"{match['query_id']},{match['subject_id']},{match['e_value']},"
                                f"{match['score']},{match['identity']},{match['alignment_length']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare CDSs between two bacterial genomes.")
    parser.add_argument("--in_gbk1", required=True, help="Input GenBank file 1")
    parser.add_argument("--in_gbk2", required=True, help="Input GenBank file 2")
    parser.add_argument("--out_file", default="match_table.csv", help="Output CSV file for matches")
    
    args = parser.parse_args()
    main(args)
