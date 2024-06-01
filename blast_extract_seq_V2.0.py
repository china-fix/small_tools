'''
xiao.fei@sund.ku.dk 
This robot can help to locally blast the query to the reference genomes and extract the matched subject seq to matched seq.fa
20220415---update remove the temp_blast.xml file, handle it in memory
-----------update tblastn as optional function for protein
20220508---update to output a .XIAO file record the summary of match or no-match things
-----------update to output DEFAULT_NAME.file_name which is the same as DEFAULT_NAME but the naming system is based on the file name of reference genome
20221027---update to add group annotation function to replace the Fixed "XIAO_ROBOT" annotation
20221031---upgrade V1.3 to add function to caculate the extracted seq types distribution

*********20240501---major upgrade, add more blast functions, and csv summary file for further data analaysis***********
20240501---update to add multiprocessing to hugely improve the running speed!
20240509---update the function relate to blastn up cut
'''


import glob
import os
import subprocess
import sys
import argparse
import io
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import logging
from itertools import product


def setup_logger(log_file):
    # Create logger for file logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Create file handler for file logging and set level to INFO
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add formatter to file handler
    file_handler.setFormatter(formatter)

    # Add file handler to logger
    logger.addHandler(file_handler)

    # Create logger for screen output
    screen_logger = logging.getLogger('screen_output')
    screen_logger.setLevel(logging.INFO)

    # Create console handler for screen output and set level to INFO
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    # Add formatter to console handler
    console_handler.setFormatter(formatter)

    # Add console handler to screen logger
    screen_logger.addHandler(console_handler)

    # Create file handler for screen output and set level to INFO
    screen_file_handler = logging.FileHandler(log_file)
    screen_file_handler.setLevel(logging.INFO)

    # Add formatter to screen file handler
    screen_file_handler.setFormatter(formatter)

    # Add screen file handler to screen logger
    screen_logger.addHandler(screen_file_handler)

    return logger, screen_logger


def parse_arguments():
    parser = argparse.ArgumentParser(description="Xiao_Fei_Robot: Local BLAST and Sequence Analysis Tool")
    parser.add_argument('--query', required=True, type=str, metavar='FILENAME',
                        help="Query FASTA filename to blast")
    parser.add_argument('--reference_folder', required=True, type=str, metavar='FOLDER',
                        help="Folder containing reference genomes for blasting")
    parser.add_argument('--cutoff_identity', default=0.8, type=float, metavar='DEFAULT 0.80',
                        help="Minimum similarity value to classify as a match")
    parser.add_argument('--cutoff_coverage', default=0.8, type=float, metavar='DEFAULT 0.80',
                        help="Minimum coverage value to classify as a match")
    parser.add_argument('--blast_type', default="tblastn", type=str, metavar='COMMAND',
                        help="Blast command name, either 'tblastn' (default) or 'blastn' (this is for up cut) or 'blastp'")
    parser.add_argument('--output', default="DEFAULT_NAME", type=str, metavar='OUTPUT_FOLDER',
                        help="Output folder name, 'DEFAULT_NAME' by default")
    parser.add_argument('--num_threads', default=4, type=int, metavar='4',
                        help="Number of threads (CPUs) to use in the BLAST search")
    parser.add_argument('--xiao_dev', default="none", type=str, metavar='advanced_data_analysis',
                        help="Advance data analysis only use by developer, 'none' by default, if you want to use, read the raw code and do the modification by yourself! quick_see; ")
    parser.add_argument('--xiao_dev_input', default=False, type=str, metavar='advanced_data_analysis_input',
                        help="extra data for advance data analysis, only use by developer, False by default, if you want to use, read the raw code and do the modification by yourself! ")
    parser.add_argument('--bp_tolerance', default=False, type=int, metavar='bp_num',
                        help="Number of bp tolerance for the aligned aa seq have no stop codon or incorrect start codon, suggest 5-10 aa")
    return parser.parse_args()

def perform_blast(query, reference, blast_type, num_threads):
    try:
        blast_result = subprocess.run([blast_type, "-query", query, "-subject", reference, "-outfmt", "5", "-num_threads", str(num_threads)],
                                      check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print("Error running BLAST command:", e)
        sys.exit(1)
    blast_xml = io.BytesIO(blast_result.stdout)
    return blast_xml

def filter_matches_combined(cutoff, blast_xml, reference_name, coverage_cutoff, blast_type, reference_file=None):
    matched_seqs = []
    blast_records = NCBIXML.parse(blast_xml)
    for blast_record in blast_records:
        query_length = blast_record.query_letters
        query = blast_record.query
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]
            alignment_coverage = hsp.align_length / query_length
            alignment_identity =hsp.identities / hsp.align_length
            if (alignment_identity  >= cutoff) and (alignment_coverage >= coverage_cutoff): #and hsp.expect <= 1e-4:
                if blast_type == 'blastp':
                    with open(reference_file, 'r') as ref_file:
                        for seq_record in SeqIO.parse(ref_file, 'fasta'):
                            if seq_record.id == alignment.hit_id:
                                matched_seq = (blast_record.query, SeqRecord(seq_record.seq, id=alignment.hit_def, description=(reference_name + " " + 
                                                                                                                                query + " " + 
                                                                                                                                str(alignment_coverage)+ " " + 
                                                                                                                                str(alignment_identity)+ " " +
                                                                                                                                str(hsp.align_length) + " " +
                                                                                                                                str(query_length))))
                                matched_seqs.append(matched_seq)
                                break
                elif blast_type == 'tblastn':
                    matched_seq = (blast_record.query, SeqRecord(Seq(hsp.sbjct), id=alignment.hit_def, description=(reference_name + " " + 
                                                                                                                                query + " " + 
                                                                                                                                str(alignment_coverage)+ " " + 
                                                                                                                                str(alignment_identity)+ " " +
                                                                                                                                str(hsp.align_length) + " " +
                                                                                                                                str(query_length))))
                    matched_seqs.append(matched_seq)
                ##else is for blastn runing in upstream cut mode
                else:
                    with open(reference_file, 'r') as ref_file:
                        for seq_record in SeqIO.parse(ref_file, 'fasta'):
                            if seq_record.id == alignment.hit_id:
                                # Calculate the start and end positions of the upstream region
                                # based on the strand direction of the subject sequence
                                upstream_length = 100
                                if hsp.sbjct_start - hsp.sbjct_end <=0:  # Positive strand
                                    upstream_start = max(0, hsp.sbjct_start - upstream_length) 
                                    upstream_end = hsp.sbjct_start
                                    upstream_seq = seq_record.seq[upstream_start-1:upstream_end-1]
                                else:  # Negative strand
                                    upstream_start = min(len(seq_record.seq), hsp.sbjct_start + upstream_length)
                                    upstream_end = hsp.sbjct_start
                                    upstream_seq = seq_record.seq[upstream_end-1:upstream_start-1].reverse_complement()
                                # Extract the upstream sequence from the subject sequence



                                
                                matched_seq = (blast_record.query, SeqRecord(Seq(upstream_seq), id=alignment.hit_def, description=(reference_name + " " + 
                                                                                                                                query + " " + 
                                                                                                                                str(alignment_coverage)+ " " + 
                                                                                                                                str(alignment_identity)+ " " +
                                                                                                                                str(hsp.align_length) + " " +
                                                                                                                                str(query_length))))
                                matched_seqs.append(matched_seq)
                                break
    return matched_seqs


def convert_to_original_format(modified_matched_seqs):
    """
    Convert the modified matched_seqs back to the original format.
    
    Args:
    - modified_matched_seqs (list): List containing tuples of (query, matched_seq_record).
    
    Returns:
    - original_matched_seqs (list): List containing SeqRecord objects representing the matched sequences.
    """
    original_matched_seqs = []
    
    for query, matched_seq_record in modified_matched_seqs:
        original_matched_seqs.append(matched_seq_record)
    
    return original_matched_seqs


# def blast_and_extract(query_file, reference_folder, cutoff, blast_type, output_file, coverage_cutoff,num_threads):
#     reference_files = glob.glob(os.path.join(reference_folder, '*'))
#     matched_seqs = []

#     for reference_file in reference_files:
#         reference_name = os.path.basename(reference_file)
#         blast_xml = perform_blast(query_file, reference_file, blast_type,num_threads)
#         is_protein = (blast_type == 'blastp')
#         matched_seqs.extend(filter_matches_combined(cutoff, blast_xml, reference_name, coverage_cutoff, is_protein, reference_file))
        
#     matched_seqs_seq = convert_to_original_format(matched_seqs)
#     SeqIO.write(matched_seqs_seq, output_file, "fasta")
#     return matched_seqs
import multiprocessing
from functools import partial  # Add this import statement
def process_reference(reference_file, query_file, cutoff, blast_type, coverage_cutoff, num_threads, matched_seqs):
    reference_name = os.path.basename(reference_file)
    blast_xml = perform_blast(query_file, reference_file, blast_type, num_threads)  
    matched_seqs_partial = filter_matches_combined(cutoff, blast_xml, reference_name, coverage_cutoff, blast_type, reference_file)
    matched_seqs.extend(matched_seqs_partial)

def blast_and_extract(query_file, reference_folder, cutoff, blast_type, output_file, coverage_cutoff, num_threads):
    reference_files = glob.glob(os.path.join(reference_folder, '*'))
    
    # Create a manager to share a list between processes
    manager = multiprocessing.Manager()
    matched_seqs = manager.list()

    # Create a multiprocessing Pool with the number of processes equal to the number of CPU cores
    with multiprocessing.Pool() as pool:
        # Define partial function with fixed arguments
        process_func = partial(process_reference, query_file=query_file, cutoff=cutoff, blast_type=blast_type, coverage_cutoff=coverage_cutoff, num_threads=num_threads, matched_seqs=matched_seqs)
        pool.map(process_func, reference_files)

    matched_seqs_seq = convert_to_original_format(matched_seqs)
    SeqIO.write(matched_seqs_seq, output_file, "fasta")
    return matched_seqs


def is_pseudo_sequence(aligned_sequence):
    """Check if the aligned sequence is a pseudo-sequence."""
    internal_stop_codon = "*"
    start_codons = {"M", "V", "L", "I"}
    stop_codon = "*"
    pseudo_status = []

    # Check if the sequence starts with the correct start codon
    if not any(aligned_sequence.startswith(start) for start in start_codons):
        pseudo_status.append('incorrect_start_codon')

    # Check for internal stop codons
    if len(aligned_sequence) > 1 and aligned_sequence[:-1].upper().find(internal_stop_codon) != -1:
        pseudo_status.append('internal_stop')

    # Check if the sequence ends with a stop codon
    if not aligned_sequence.endswith(stop_codon):
        pseudo_status.append('no_stop_codon')

    # If all checks pass, mark the sequence as 'Intact'
    if not pseudo_status:
        pseudo_status = 'Intact'
    else:
        pseudo_status = ','.join(pseudo_status)

    return pseudo_status


def generate_summary(output_file, blast_type, reference_folder,query_file, bp_tolerance):
    try:
        subprocess.run(["seqkit", "version"], check=True, stdout=subprocess.PIPE)
    except FileNotFoundError:
        print("Error: seqkit is not installed or not found in the system path.")
        print("Please install seqkit to perform sequence analysis.")
        sys.exit(1)

    cmd = 'seqkit fx2tab ' + output_file + ' > ' + output_file + '.tab'
    os.system(cmd)

    try:
        subprocess.run(["seqkit", "version"], check=True, stdout=subprocess.PIPE)
    except FileNotFoundError:
        print("Error: seqkit is not installed or not found in the system path.")
        print("Please install seqkit to perform sequence analysis.")
        sys.exit(1)

    cmd = 'seqkit fx2tab ' + query_file + ' > ' + output_file + '_query_file.tab'
    os.system(cmd)

    tab_df = pd.read_csv(output_file + '.tab', sep='\t', names=['assembly_name', 'seq'], index_col=False)
    split_names = tab_df['assembly_name'].str.rsplit(n=6, expand=True)
    tab_df[['subject_id', 'subject_file', 'query_id','alignment_coverage','alignment_identity','align_length','query_length']] = split_names
    if blast_type == 'blastp' or blast_type == 'tblastn':
        tab_df['is_pseudo'] = tab_df['seq'].apply(is_pseudo_sequence)
    else:
        tab_df['is_pseudo'] = 'NO_NEED'
        
        

    query_id_list=pd.read_csv(output_file + '_query_file.tab', sep='\t',names=['query_id','query_seq'],index_col=False)['query_id']
    subject_file_list_raw = glob.glob(os.path.join(reference_folder, '*'))
    subject_file_list = [os.path.basename(file_path) for file_path in subject_file_list_raw]
    # filter_list=tab_df['subject_file'].unique()
    link=list(product(query_id_list,subject_file_list))
    df_link =  pd.DataFrame(link, columns=['query_id','subject_file'])
    # print(df_link)
    tab_df_merge=pd.merge(df_link, tab_df, left_on=['query_id','subject_file'], right_on=['query_id','subject_file'], how='left')
    # print(tab_df_merge)
    tab_df_merge['is_pseudo'] = tab_df_merge['is_pseudo'].fillna('Absent')
    tab_df_merge[['alignment_coverage','alignment_identity']] = tab_df_merge[['alignment_coverage','alignment_identity']].fillna(0)
    tab_df_merge['seq'].fillna('NO_MATCH', inplace=True)

    # tab_df_merge.to_csv(output_file + "_SUM" + ".csv")
    # tab_df_merge.to_pickle(output_file + "_SUM" + ".pickle")

    ## steps added 2024-05-07
    if bp_tolerance:
        print("Hey xiao, bp tolerance mode is on, you set bp(aa) is "+str(bp_tolerance))
        tab_df_merge[['alignment_coverage','query_length','align_length']] = tab_df_merge[['alignment_coverage','query_length','align_length']].astype(float)
        # Add a new column using pd.query
        tab_df_merge['length_difference'] = tab_df_merge.query('query_length.notna() and align_length.notna()')['query_length'] - tab_df_merge.query('query_length.notna() and align_length.notna()')['align_length']
        # Filter rows where conditions are met
        filtered_rows = tab_df_merge.query('(length_difference < 20) and (is_pseudo.isin(["no_stop_codon", "incorrect_start_codon"]))')
        # Update 'is_pseudo' column for filtered rows to be 'Intact'
        tab_df_merge.loc[filtered_rows.index, 'is_pseudo'] = 'Intact'
        
    tab_df_merge.to_csv(output_file + "_SUM" + ".csv")
    tab_df_merge.to_pickle(output_file + "_SUM" + ".pickle")

    return tab_df_merge



### quick_see module is use to analysis the _sum.csv, summarized the matched gene distribution
def quick_see(input, output_file):
    df_one_match_each_queryfile = input.groupby(['query_id', 'subject_file'], as_index=False).apply(lambda x: x.sort_values('alignment_identity', ascending=False).iloc[0])
    result = df_one_match_each_queryfile.groupby(['query_id','seq'])['is_pseudo'].apply(lambda x: x.value_counts()).reset_index(name='count')
    
    
    # Add a new column 'mutation' with default value 'pseudo'
    result['mutation'] = 'Pseudo'

    # Update 'mutation' column where 'level_2' is 'Absent' or 'Intact'
    result.loc[result['level_2'] == 'Absent', 'mutation'] = 'Absent'
    result.loc[result['level_2'] == 'Intact', 'mutation'] = 'Intact'
    
    result['total_subject_count'] = result.groupby(['query_id'])['count'].transform('sum')
    result['percentage'] = (result['count'] / result['total_subject_count'])

    # Function to calculate mutation_score based on mutation type
    def calculate_mutation_score(row):
        if row['mutation'] == 'Pseudo':
            return row['percentage'] * 0.5
        elif row['mutation'] == 'Intact':
            return row['percentage'] * 1
        elif row['mutation'] == 'Absent':
            return row['percentage'] * 0
    # Apply the function to create the new column 'mutation_score'
    result['mutation_score'] = result.apply(calculate_mutation_score, axis=1)

    
    result.to_csv(output_file+"_quicksee.csv")
    result.to_pickle(output_file+"_quicksee.pickle")



    T=result.groupby(['query_id','level_2'])[['count','percentage']].sum().reset_index()
    T.to_csv(output_file+"_quicksee_1.csv")

    K=result.groupby(['query_id','mutation'])[['count','percentage']].sum().reset_index()
    query_id =K['query_id'].unique()
    mutation =['Absent','Pseudo','Intact']
    link=list(product(query_id,mutation))
    df_link =  pd.DataFrame(link, columns=['query_id','mutation'])
    new_K= pd.merge(df_link,K, on=['query_id','mutation'], how='left')
    new_K.fillna(0, inplace=True)
    new_K.to_csv(output_file+"_quicksee_2.csv")

    F=result.groupby('query_id')[['count','mutation_score']].sum().reset_index()
    F.to_csv(output_file+"_quicksee_3.csv")

    return result

 # get a subset of quicksee use the sublist   
def quick_see_plus(result_in,sublist,output_file):
    with open(sublist, 'r') as file:
        # Read lines into a list
        lines = [line.strip() for line in file.readlines()]
    result =result_in[result_in['query_id'].isin(lines)]
    result.to_csv(output_file+"_quicksee_plus.csv")
    result.to_pickle(output_file+"_quicksee_plus.pickle")


    T=result.groupby(['query_id','level_2'])[['count','percentage']].sum().reset_index()
    T.to_csv(output_file+"_quicksee_plus_1.csv")

    K=result.groupby(['query_id','mutation'])[['count','percentage']].sum().reset_index()
    query_id =K['query_id'].unique()
    mutation =['Absent','Pseudo','Intact']
    link=list(product(query_id,mutation))
    df_link =  pd.DataFrame(link, columns=['query_id','mutation'])
    new_K= pd.merge(df_link,K, on=['query_id','mutation'], how='left')
    new_K.fillna(0, inplace=True)
    new_K.to_csv(output_file+"_quicksee_plus_2.csv")

    F=result.groupby('query_id')[['count','mutation_score']].sum().reset_index()
    F.to_csv(output_file+"_quicksee_plus_3.csv")

    return result
    



def main():
    args = parse_arguments()
    
    # Create the output folder if it doesn't exist
    output_folder = args.output # Specify your desired output folder path
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Construct the full output file path
    output_file = os.path.join(output_folder, args.output)
    print("The result is located at " + output_file)

    # Create a log file path
    log_file = os.path.join(output_folder, "logfile.txt")
    
    # Setup logger to output to both console and file
    logger, screen_logger = setup_logger(log_file)

   # Log the command used
    command_used = " ".join(sys.argv)
    logger.info("Command used: %s", command_used)
    # Log screen information
    screen_logger.info("Starting script execution...")


    # Perform BLAST search and analysis
    blast_and_extract(args.query, args.reference_folder, args.cutoff_identity,
                      args.blast_type, output_file, args.cutoff_coverage, args.num_threads)
    
    # Generate summary and save it to the output folder
    tab_df = generate_summary(output_file, args.blast_type, args.reference_folder,args.query,bp_tolerance=args.bp_tolerance)
    
    print('Xiao_Fei_Robot: Job done!')

    if not args.xiao_dev == 'none':
        print("hey xiao, i am doing advance module now!")
        if args.xiao_dev == 'quick_see':
            print("run quick_see module now!")
            quick_see(tab_df, output_file)
        elif args.xiao_dev == 'quick_see_plus':
            print("run quick_see module now!")
            resultin=quick_see(tab_df, output_file)
            print("run quick_see_plus module now!")
            quick_see_plus(result_in=resultin, sublist=args.xiao_dev_input, output_file=output_file)
        else:
            print("hey xiao, no match to you advance function!")
    print('THE END')



if __name__ == '__main__':
    sys.exit(main())
