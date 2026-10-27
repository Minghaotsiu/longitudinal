# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 15:08:41 2020

@author: cyclopenta
"""

import pandas as pd
import os
import re
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description='整合各样本的STAR比对信息（默认为./log_final_out目录下所有*Log.final.out文件）: \
                                              python inf_star_map_ver0.1.py \
                                              -i input dir \
                                              -o output')

parser.add_argument('-i', '--input_dir', type=str, default='./log_final_out', help='url to /log_final_out/')
parser.add_argument('-o', '--output_prefix', type=str, required=True, help='output_file_prefix')

args = parser.parse_args()

# initialization
total_info = Counter()
interest_info = Counter()
temp_dic = Counter()
#define what you want
interest_ini = ['Filename', 'Sample', 'Uniquely mapped reads %', '% of reads mapped to multiple loci', '% of reads unmapped: too short', 'Number of input reads']
for ini in interest_ini:
    interest_info[ini] ={}

#open dir and read files    
file_path = args.input_dir
file_total = os.listdir(file_path)
file_total = sorted(file_total)

for f in file_total:
    if not os.path.isdir(f):
        file = open(file_path+'/'+f)
        onefile_info = Counter()
        onefile_info['Filename'] = re.sub('Log.final.out', '', f)
        onefile_info['Sample'] = re.sub('Log.final.out', '', f).split('_')[0]
        for line in file:
            line = line.lstrip()
            i = line.split('|')
            # prevent list index error
            if len(i) > 1:
                onefile_info[i[0].rstrip()] = i[1].lstrip('\t').rstrip()
        total_info[f] = onefile_info
    
for k,v in total_info.items():
    for m,n in interest_info.items():
        interest_info[m][k]=v[m]

final_df = pd.DataFrame(interest_info)
final_df.columns = ['Filename','Sample','Uniquely mapped reads %','Reads mapped to multiple loci %','% of reads unmapped: too short','Number of input reads']

# Ensure the columns are strings before processing
final_df['Uniquely mapped reads %'] = final_df['Uniquely mapped reads %'].astype(str)
final_df['Reads mapped to multiple loci %'] = final_df['Reads mapped to multiple loci %'].astype(str)

# Calculate Mapped reads %
final_df['Mapped reads %'] = final_df['Uniquely mapped reads %'].str.strip('%').astype(float) + final_df['Reads mapped to multiple loci %'].str.strip('%').astype(float)
final_df['Mapped reads %'] = final_df['Mapped reads %'].apply(lambda x: format(x, '.2f')).astype(str) + '%'

final_df = final_df[['Filename','Sample','Uniquely mapped reads %','Reads mapped to multiple loci %','Mapped reads %','% of reads unmapped: too short','Number of input reads']]
final_df.to_csv(args.output_prefix+'_RNA_MAP.csv', index=False, sep=',')

# Print a message indicating successful completion
print(f"Summary statistics saved to {args.output_prefix}_RNA_MAP.csv")