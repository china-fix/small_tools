"""
Hi there,
This is a conserved SNP (core-SNP) analysis workflow followed on the matrix table created by Enterobase SNP analysis pipeline.
1. reading into the matrix as dataframe
2. parse the dataframe, and adding groups as user defined
3. using for loop to go through all the SNP mutation sites, and extract the core-SNP sites between groups

Enjoy!
Xiao Fei (china-fixing@hotmail.com, xiao.fei@sund.ku.dk)
update20230301, rename output file if file already exsited
update20250416, read matrix updated
"""

import sys
# import subprocess
import argparse
import pandas as pd
import os
import gzip

def parse_args():
    parser=argparse.ArgumentParser(description='''
Welcome to use Xiao's robot
Hi there,
This is a conserved SNP (core-SNP) analysis workflow followed on the matrix table created by Enterobase SNP analysis pipeline.
1. reading into the matrix as dataframe
2. parse the dataframe, and adding groups as user defined
3. using for loop to go through all the SNP mutation sites, and extract the core-SNP sites between groups

downstream usage examples for the output of bed-file (e.g. core_snp.bed):
grep  'Conserved_SNP' core_snp.bed > c.bed
bedtools annotate -counts -i A.gff -files c.bed > core_snp_annotated.bed
grep -P '\\t[1-9]\d*$' core_snp_annotated.bed > core_snp_annotated_positive.bed


Enjoy!
Xiao Fei (china-fixing@hotmail.com, xiao.fei@sund.ku.dk)
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--IN_MATRIX', required=True, type=str, metavar='FILENAME', help="the matrix table created by Enterobase SNP analysis pipeline")
    parser.add_argument('--IN_GROUP', required=True, type=str, metavar='FILENAME', help="user defined grouping table for the strians in matrix file, can be modifed from the mapping file created by Enterobase SNP analysis pipeline，第一列文件名，第二列A或B")
    parser.add_argument('--OUTPUT', default="core_snp_ann.tab", type=str, metavar='FILENAME', help="Output bed-like file name, default is core_snp_ann.tab")
    # parser.add_argument('--REF_NAME', default="XF_ROBOT", type=str, metavar='FILENAME', help="the reference genome name used for the bed header output (inside)")
    parser.add_argument('--CUT_OFF', default=0.01, type=int, metavar='CUTOFF', help="cutoff for random error torlence for the SNP mutation, default is 0.01")
    return parser.parse_args()
    
### reading into the matrix as dataframe ### 
def read_matrix(input):
    # 第一步：定位标题行行号
    with gzip.open(input, 'rt') as f:
        header_line_num = 0
        for line in f:
            if line.startswith('#Seq'):
                break  # 找到标题行，停止搜索
            header_line_num += 1
        else:  # 若循环未找到标题行
            raise ValueError("ERROR: File missing '#Seq' header line")

    # 第二步：读取数据（跳过所有注释行）
    read_in = pd.read_csv(
        input,
        sep='\t',
        header=0,  # 关键修改：跳过后剩余的第一行作为标题
        skiprows=lambda x: x < header_line_num,  # 跳过标题行之前的行
        low_memory=False
    )
    
    # 第三步：处理数据
    snp_matrix = read_in.iloc[:,:-1]
    snp_matrix_ann = read_in.iloc[:,[1,-1]]
    snp_matrix_t = snp_matrix.set_index(['#Seq','#Site']).T
    
    return snp_matrix_t, snp_matrix_ann

### add grouping info ### 
def add_set(input, snp_matrix_t):
    map = pd.read_csv(input, sep='\t', index_col=0, names=['meta_name', 'set'])
    data = pd.concat([snp_matrix_t, map], axis=1)
    data = data.dropna()
    return data

### core SNP analysis ###
def analysis(data, cutoff):
    print(data.groupby('set').count())
    for sites in data.columns[:-2]:
        c_run = data.groupby('set')[[sites]].value_counts()
        
        c_run_A = c_run['A']/c_run['A'].sum()
        kk_A = c_run_A[c_run_A > cutoff]
        set_A = set(kk_A.index.values)

        c_run_B = c_run['B']/c_run['B'].sum()
        kk_B = c_run_B[c_run_B > cutoff]
        set_B = set(kk_B.index.values)

        result = len(set_A.intersection(set_B))
        tsites = str(sites[1])
        ref_name = str(sites[0])
        if result == 0:
            print(ref_name+'\t'+tsites+'\t'+ 'Conserved_SNP')
            print(ref_name+'\t'+tsites+'\t'+ tsites + '\t' + 'Conserved_SNP', file=open("core_snp.bed", "a"))
        else:
            print(ref_name+'\t'+tsites+'\t'+ tsites + '\t' + 'not_SNP_Conserved', file=open("core_snp.bed", "a"))

#main
def main():
    # print("main work")
    args=parse_args()
    print(args)
    original = "core_snp.bed"
    output = "core_snp_old.bed"
    try:
        os.rename(original, output)
    except FileNotFoundError:
        pass

    snp_matrix_t, snp_matrix_ann = read_matrix(args.IN_MATRIX)
    data = add_set(args.IN_GROUP, snp_matrix_t)
    analysis(data, args.CUT_OFF)

    bed_output = pd.read_csv('core_snp.bed', sep='\t', names=['reference', '#Site','pos2','ann'])
    bed_output_ann = bed_output.merge(snp_matrix_ann, on='#Site', how='left')
    bed_output_ann.to_csv(args.OUTPUT, sep='\t', index=False, header=False)

    print("Hi xiao, job done!")

if __name__ == '__main__':
    sys.exit(main())


