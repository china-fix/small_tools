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
update20250416, output directory management
"""

import sys
import argparse
import pandas as pd
import os
import gzip

def parse_args():
    parser = argparse.ArgumentParser(
        description=r'''
Welcome to use Xiao's robot
Hi there,
This is a conserved SNP (core-SNP) analysis workflow followed on the matrix table created by Enterobase SNP analysis pipeline.
1. reading into the matrix as dataframe
2. parse the dataframe, and adding groups as user defined
3. using for loop to go through all the SNP mutation sites, and extract the core-SNP sites between groups

Welcome to use Xiao's Core SNP Analysis Robot!

[主要功能]
1. 读取Enterobase生成的SNP矩阵文件
2. 根据用户分组信息进行核心SNP分析
3. 输出带注释的核心SNP结果

[输入输出]
输入文件要求：
- 矩阵文件：Enterobase SNP分析生成的压缩矩阵文件
- 分组文件：制表符分隔，第一列样本名，第二列样品别名，第三列分组（A/B）

输出文件结构：


[使用示例]
python core_SNP.py \
  --IN_MATRIX align.matrix.gz \
  --IN_GROUP sample_groups.list \
  --OUTPUT 20240524_results \
  --CUT_OFF 0.01

Xiao Fei (china-fixing@hotmail.com, xiao.fei@sund.ku.dk， xiaofei@fafu.edu.cn)

        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--IN_MATRIX', required=True, type=str, 
                      help="Input matrix file path (gzipped)")
    parser.add_argument('--IN_GROUP', required=True, type=str,
                      help="Grouping file path (TSV format)")
    parser.add_argument('--OUTPUT', required=True, type=str,
                      help="Output directory path")
    parser.add_argument('--CUT_OFF', default=0.01, type=float,
                      help="SNP filtering cutoff (default: 0.01)")
    return parser.parse_args()

def handle_output_dir(path):
    """处理输出目录"""
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise ValueError(f"Path exists but is not a directory: {path}")
        
        print(f"Output directory exists: {path}")
        choice = input("Continue? Existing files may be overwritten. [y/N]: ").lower()
        if choice not in ('y', 'yes'):
            print("Operation aborted by user")
            sys.exit(0)
    else:
        os.makedirs(path, exist_ok=True)
        print(f"Created output directory: {path}")

def read_matrix(input):
    # 保持原有读取逻辑不变
    with gzip.open(input, 'rt') as f:
        header_line_num = 0
        for line in f:
            if line.startswith('#Seq'):
                break
            header_line_num += 1
        else:
            raise ValueError("ERROR: File missing '#Seq' header line")

    read_in = pd.read_csv(
        input,
        sep='\t',
        header=0,
        skiprows=lambda x: x < header_line_num,
        low_memory=False
    )
    
    snp_matrix = read_in.iloc[:,:-1]
    snp_matrix_ann = read_in.iloc[:,[1,-1]]
    snp_matrix_t = snp_matrix.set_index(['#Seq','#Site']).T
    
    return snp_matrix_t, snp_matrix_ann

def add_set(input, snp_matrix_t):
    # 保持原有分组逻辑不变
    map = pd.read_csv(input, sep='\t', index_col=0, header=None, names=['meta_name', 'set'])
    data = pd.concat([snp_matrix_t, map], axis=1)
    return data.dropna()

def analysis(data, cutoff, output_dir):
    # 修改输出路径
    output_bed = os.path.join(output_dir, "core_snp.bed")
    
    
    # 保留原有分析逻辑
    print(data)
    print(data.groupby('set').count())
    
    with open(output_bed, 'w') as f:  # 使用with确保文件正确关闭
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
            
            line_content = f"{ref_name}\t{tsites}\t{tsites}\t"
            if result == 0:
                line_content += "Conserved_SNP"
                print(line_content)
            else:
                line_content += "not_SNP_Conserved"
            
            f.write(line_content + "\n")

def main():
    args = parse_args()
    
    # 处理输出目录
    output_dir = os.path.abspath(args.OUTPUT)
    handle_output_dir(output_dir)
    
    # 处理旧文件
    original = os.path.join(output_dir, "core_snp.bed")
    backup = os.path.join(output_dir, "core_snp_old.bed")
    try:
        if os.path.exists(original):
            os.rename(original, backup)
    except Exception as e:
        print(f"Notice: {str(e)}")

    # 保持原有处理流程
    snp_matrix_t, snp_matrix_ann = read_matrix(args.IN_MATRIX)
    data = add_set(args.IN_GROUP, snp_matrix_t)
    analysis(data, args.CUT_OFF, output_dir)

    # 生成最终结果
    final_output = os.path.join(output_dir, "core_snp_ann.tab")
    bed_output = pd.read_csv(os.path.join(output_dir, "core_snp.bed"), 
                           sep='\t', 
                           names=['reference', '#Site','pos2','ann'])
    bed_output_ann = bed_output.merge(snp_matrix_ann, on='#Site', how='left')
    bed_output_ann.to_csv(final_output, sep='\t', index=False, header=False)

    print(f"\nJob completed! Results saved in: {output_dir}")
    print("Files generated:")
    print(f"- {os.path.join(output_dir, 'core_snp.bed')}")
    print(f"- {final_output}")

if __name__ == '__main__':
    sys.exit(main())