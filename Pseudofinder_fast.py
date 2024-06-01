'''
xiao.fei@sund.ku.dk 
This is designed for 'pseudofinder.py annotate', to run it in parallel for larger genomes to be caculate

'''

import os
import subprocess
import argparse
import concurrent.futures

# Function to run Pseudofinder on a single genome file
def run_pseudofinder(genome_file, output_dir, pseudofinder_path, blast_database, threads_per_run, extra_params):
    base_name = os.path.basename(genome_file)
    for ext in ['.gbk', '.gbff', '.gbf']:
        if base_name.endswith(ext):
            genome_name = base_name[:-len(ext)]
            break
    output_prefix = os.path.join(output_dir, genome_name)
    
    command = [
        'python', pseudofinder_path, 'annotate',
        '-g', genome_file,
        '-db', blast_database,
        '-op', output_prefix,
        '-t', str(threads_per_run)
    ] + extra_params
    
    try:
        subprocess.run(command, check=True)
        print(f"Successfully processed {genome_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {genome_file}: {e}")

# Function to get all .gbk files from the input directory
def get_genome_files(input_dir):
    return [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(('.gbk', '.gbff', '.gbf'))]


# Main function to run Pseudofinder in parallel
def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    genome_files = get_genome_files(args.input_dir)
    
    extra_params = []
    if args.intergenic_length:
        extra_params += ['-i', str(args.intergenic_length)]
    if args.length_pseudo:
        extra_params += ['-l', str(args.length_pseudo)]
    if args.shared_hits:
        extra_params += ['-s', str(args.shared_hits)]
    if args.evalue:
        extra_params += ['-e', str(args.evalue)]
    if args.hitcap:
        extra_params += ['-hc', str(args.hitcap)]
    if args.contig_ends:
        extra_params.append('--contig_ends')
    if args.intergenic_threshold:
        extra_params += ['-it', str(args.intergenic_threshold)]
    if args.reference:
        extra_params += ['-ref', args.reference]
    if args.diamond:
        extra_params.append('--diamond')
    if args.skip_checkdb:
        extra_params.append('--skip_checkdb')
    if args.no_bidirectional_length:
        extra_params.append('--no_bidirectional_length')
    if args.use_alignment:
        extra_params.append('--use_alignment')
    if args.perc_cov:
        extra_params += ['-cov', str(args.perc_cov)]
    if args.perc_id:
        extra_params += ['-id', str(args.perc_id)]
    if args.max_dnds:
        extra_params += ['-dnds', str(args.max_dnds)]
    if args.use_deviation:
        extra_params.append('--use_deviation')

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_parallel_runs) as executor:
        futures = [
            executor.submit(
                run_pseudofinder, 
                genome_file, 
                args.output_dir, 
                args.pseudofinder_path, 
                args.blast_database, 
                args.threads_per_run, 
                extra_params
            ) for genome_file in genome_files
        ]
        concurrent.futures.wait(futures)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Pseudofinder in parallel for multiple GenBank files.')
    parser.add_argument('--pseudofinder_path', required=True, help='Path to the Pseudofinder script')
    parser.add_argument('--blast_database', required=True, help='Path to the BLAST database')
    parser.add_argument('--input_dir', required=True, help='Directory containing input GenBank files')
    parser.add_argument('--output_dir', required=True, help='Directory to store the output')
    parser.add_argument('--threads_per_run', type=int, default=4, help='Number of threads to use per Pseudofinder run')
    parser.add_argument('--max_parallel_runs', type=int, default=4, help='Maximum number of parallel Pseudofinder runs')

    # Additional optional arguments for Pseudofinder
    parser.add_argument('-i', '--intergenic_length', type=int, help='Length of intergenic regions to check, default is 30 bp.')
    parser.add_argument('-l', '--length_pseudo', type=float, help='Percentage of length for pseudo candidates, default is 0.65.')
    parser.add_argument('-s', '--shared_hits', type=float, help='Percentage of blast hits that must be shared to join two nearby regions, default is 0.5.')
    parser.add_argument('-e', '--evalue', type=float, help='E-value for blast searches, default is 1e-4.')
    parser.add_argument('-hc', '--hitcap', type=int, help='Maximum number of allowed hits for BLAST, default is 15.')
    parser.add_argument('--contig_ends', action='store_true', help='Include intergenic regions at contig ends.')
    parser.add_argument('-it', '--intergenic_threshold', type=float, help='Number of BlastX hits needed to annotate an intergenic region as a pseudogene, default is 0.3.')
    parser.add_argument('-ref', '--reference', help='Reference genome for maximum-likelihood phylogenentic analysis using PAML.')
    parser.add_argument('--diamond', action='store_true', help='Use DIAMOND BLAST as the search engine.')
    parser.add_argument('--skip_checkdb', action='store_true', help='Skip the check if a dmnd or blast db already exists for the provided database.')
    parser.add_argument('--no_bidirectional_length', action='store_true', help='Do not check for run-on pseudogenes.')
    parser.add_argument('--use_alignment', action='store_true', help='Consider alignment quality between query and hits.')
    parser.add_argument('-cov', '--perc_cov', type=float, help='Proportion of reference gene length for frameshift and dN/dS analysis, default is 0.25.')
    parser.add_argument('-id', '--perc_id', type=float, help='Percent identical to reference gene for frameshift and dN/dS analysis, default is 0.25.')
    parser.add_argument('-dnds', '--max_dnds', type=float, help='Maximum dN/dS value for gene to be considered "intact", default is 0.3.')
    parser.add_argument('--use_deviation', action='store_true', help='Call pseudogenes if the gene length falls outside of two standard deviations from the mean blast hit length.')

    args = parser.parse_args()
    main(args)
