import pandas as pd
import os
import re
import argparse

# 创建命令行参数解析器
parser = argparse.ArgumentParser(description='整合各样本的计数信息：生成总的表达矩阵')
parser.add_argument('-c', '--count_dir', type=str, default='./', help='directory containing *_union_counts.txt files')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='output_file_prefix')

# 解析命令行参数
args = parser.parse_args()

# 获取所有 *_union_counts.txt 文件
count_dir = args.count_dir
count_files = [f for f in os.listdir(count_dir) if f.endswith('_union_counts.txt')]
count_files = sorted(count_files)

# 初始化一个空的 DataFrame
expression_matrix = pd.DataFrame()

# 读取每个文件并将数据添加到 DataFrame 中
for count_file in count_files:
    sample_name = re.sub('_union_counts.txt', '', count_file)
    # 使用 pandas 读取文件
    df = pd.read_csv(os.path.join(count_dir, count_file), sep='\t', names=['Gene_ID', sample_name])
    # 过滤 htseq-count 的汇总行（__no_feature / __not_aligned 等）
    df = df[~df['Gene_ID'].astype(str).str.startswith('__')]
    # 只保留 Gene_ID 和 counts 列
    df = df.set_index('Gene_ID')
    
    # 合并数据
    if expression_matrix.empty:
        expression_matrix = df
    else:
        expression_matrix = pd.merge(expression_matrix, df, left_index=True, right_index=True, how='outer')

# 填充 NaN 值为 0
expression_matrix = expression_matrix.fillna(0)

# 保存表达矩阵
output_file = args.output_prefix + '_counts_matrix.csv'
expression_matrix.to_csv(output_file, index=True, sep=',')

print(f"Counts matrix saved to {output_file}")