import pandas as pd
import os
import glob
import re
import os
from collections import Counter
import argparse

# 创建命令行参数解析器
parser = argparse.ArgumentParser(description='整合各样本的STAR比对信息（默认为当前目录下所有*Log.final.out文件）: \
                                              python inf_star_map_ver0.1.py \
                                              -i input dir \
                                              -o output')

# 添加命令行参数
parser.add_argument('-i', '--input_dir', type=str, default='.', help='url to log.final.out/')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='output_file_prefix')
parser.add_argument('-c', '--count_dir', type=str, default='./', help='directory containing *_union_counts.txt files')
parser.add_argument('-g', '--gtf', type=str, default=os.environ.get('SOYBEAN_GTF'), help='gene annotation GTF file')

# 解析命令行参数
args = parser.parse_args()
if not args.gtf:
    raise SystemExit("Set SOYBEAN_GTF or pass --gtf.")

# 初始化计数器
total_info = Counter()
interest_info = Counter()

# 定义感兴趣的项
interest_ini = ['Filename', 'Sample', 'Uniquely mapped reads %', '% of reads mapped to multiple loci', '% of reads unmapped: too short', 'Number of input reads']
for ini in interest_ini:
    interest_info[ini] = {}

# 读取目录和文件
file_path = args.input_dir
file_total = os.listdir(file_path)
if not any(x.endswith('Log.final.out') for x in file_total):
    for cand in ['log_final_out'] + glob.glob('../*/log_final_out') + glob.glob('../*/*/log_final_out'):
        try:
            fl = os.listdir(cand)
        except OSError:
            continue
        if any(x.endswith('Log.final.out') for x in fl):
            file_path, file_total = cand, fl
            print('[i] auto-detected log dir:', cand)
            break
file_total = sorted(file_total)

# 打印调试信息
print(f"Reading from directory: {file_path}")
print(f"Files found: {file_total}")

# 读取对齐日志文件
for f in file_total:
    if not os.path.isdir(f) and f.endswith('Log.final.out'):
        with open(os.path.join(file_path, f)) as file:
            onefile_info = Counter()
            onefile_info['Filename'] = re.sub(r'Log\.final\.out', '', f)
            onefile_info['Sample'] = re.sub(r'Log\.final\.out$', '', f)
            for line in file:
                line = line.lstrip()
                i = line.split('|')
                if len(i) > 1:
                    onefile_info[i[0].rstrip()] = str(i[1].lstrip('\t').rstrip())
        total_info[f] = onefile_info

# 提取感兴趣的信息
for k, v in total_info.items():
    for m in interest_ini:
        # 如果不是字符串，则尝试转换为字符串
        value = str(v.get(m, '')) if not isinstance(v.get(m), str) else v.get(m)
        interest_info[m][k] = value

# 读取计数文件并计算 FPKM
count_dir = args.count_dir
gtf_file = args.gtf

# 加载基因长度
gene_lengths = {}
with open(gtf_file, 'r') as gtf:
    for line in gtf:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if parts[2] == 'gene':
            gene_id = re.search(r'gene_id "(.*?)"', parts[8]).group(1)
            start, end = map(int, parts[3:5])
            gene_lengths[gene_id] = end - start + 1

# 读取计数文件并计算 FPKM
count_files = [f for f in os.listdir(count_dir) if f.endswith('_union_counts.txt')]
count_files = sorted(count_files)

# 计算每个样本的总读段计数
total_counts = {}
for count_file in count_files:
    sample_name = re.sub('_union_counts.txt', '', count_file)
    total_counts[sample_name] = 0
    with open(os.path.join(count_dir, count_file), 'r') as file:
        for line in file:
            _, count = line.strip().split('\t')
            total_counts[sample_name] += int(count)

# 计算 FPKM
fpkm_data = {}
for count_file in count_files:
    sample_name = re.sub('_union_counts.txt', '', count_file)
    fpkm_data[sample_name] = {}
    with open(os.path.join(count_dir, count_file), 'r') as file:
        for line in file:
            gene_id, count = line.strip().split('\t')
            count = int(count)
            length = gene_lengths.get(gene_id, 1)
            fpkm = (count / length) / (total_counts[sample_name] / 1e9)
            fpkm_data[sample_name][gene_id] = fpkm

# 创建一个 DataFrame 包含 FPKM 数据
fpkm_df = pd.DataFrame.from_dict(fpkm_data)

# 将 FPKM 数据添加到 DataFrame
final_df = pd.DataFrame(interest_info)
final_df.columns = ['Filename', 'Sample', 'Uniquely mapped reads %', 'Reads mapped to multiple loci %', '% of reads unmapped: too short', 'Number of input reads']

# 强制转换列为字符串类型
final_df['Uniquely mapped reads %'] = final_df['Uniquely mapped reads %'].astype(str)
final_df['Reads mapped to multiple loci %'] = final_df['Reads mapped to multiple loci %'].astype(str)

# 添加 FPKM 列
# FPKM 单独输出到 *_fpkm.txt，MAP 表只保留比对统计

# 计算 Mapped reads %
final_df['Mapped reads %'] = final_df['Uniquely mapped reads %'].str.strip('%').astype(float) + final_df['Reads mapped to multiple loci %'].str.strip('%').astype(float)
final_df['Mapped reads %'] = final_df['Mapped reads %'].apply(lambda x: format(x, '.2f')).astype(str) + '%'

# 重新排序列
final_df = final_df[['Filename', 'Sample', 'Uniquely mapped reads %', 'Reads mapped to multiple loci %', 'Mapped reads %', '% of reads unmapped: too short', 'Number of input reads']]

# 保存 CSV 文件
output_file = args.output_prefix + '_RNA_MAP.csv'
final_df.to_csv(output_file, index=False, sep=',')

print(f"Output saved to {output_file}")

# 保存 FPKM 数据到一个单独的文本文件
fpkm_output_file = args.output_prefix + '_fpkm.txt'
fpkm_df.to_csv(fpkm_output_file, index=True, sep='\t')

print(f"FPKM data saved to {fpkm_output_file}")
