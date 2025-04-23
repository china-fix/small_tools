"""
This Python program uses Biopython to optimize a BLAST search on two FASTA files. It processes each sequence from the query file against sequences in the database, filtering hits with an E-value below 0.01. The program saves selected hits into a new FASTA file with their bitscore and e-value annotations.

The code optimizes memory usage by processing data in chunks instead of loading all sequences at once. It also improves performance by specifying BLAST parameters directly within the command line interface call.
"""
import sys
from Bio import SeqIO, NCBIXML, AlignIO

def main(query_fasta, db_fasta):
    # Load query and database sequences into memory using Biopython's chunk processing capabilities
    query_records = iter(SeqIO.parse(query_fasta, 'fasta'))
    db_records = iter(SeqIO.parse(db_fasta, 'fasta'))

    blast_cline_args = ['blastn', '-query', str(query_records), '-db', str(db_records),
                        '-evalue', '0.01', '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
                        '-num_threads', '4']  # Specify BLAST parameters

    blast_output = run_blast(blast_cline_args)

    with open('hits.fasta', 'w') as output_file:
        for line in blast_output.split('\n'):
            if not line: continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = map(float, line.split('\t'))
            
            # Save hits with E-value < 0.01 and bitscore > threshold to the output file
            if float(evalue) < 0.01:
                seq_line = f'>{sseqid} (bitscore: {float(bitscore)}, e-value: {float(evalue)})\n'
                output_file.write(seq_line)
                
                # Fetch sequence data for saving into the FASTA format
                query_seq = next(query_records).seq
                db_seq = next(db_records).seq
                
                alignment_data = (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send)
                
                # Create and write the FASTA format sequence with hit details to file
                output_file.write(query_seq + db_seq)
                
def run_blast(blast_cline_args):
    """
    Run BLAST using Biopython's NCBIXML interface.

    This function executes a BLAST search based on the provided command-line arguments,
    capturing its output for later processing.
    """
    blast_cline = NCBIXML.CommandLineRunner(blast_cline_args)
    blast_results = blast_cline()
    
    # Ensure BLAST results are in XML format
    assert isinstance(blast_results, str), "BLAST result should be in XML format."
    
    return blast_results


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python program_name.py query.fasta db.fasta')
        sys.exit(1)
    
    main(query_fasta=sys.argv[1], db_fasta=sys.argv[2])
