# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 17:24:13 2020

@author: cyclopenta
"""

import pandas as pd
import sys
import os
import re
import datetime

#define parameters
adapter_fasta = 'S1.adapters_RNA.fasta'
now_time = datetime.datetime.now()
time_str = now_time.strftime('%Y-%m-%d-%H_%M_%S')
trim_path = 'trim_fastqc_'+time_str
fastqc_out_path = 'fastqc_out_'+time_str
os.system(''' mkdir -p ./{path1}/{path2}'''.format(path1 = trim_path, path2 = fastqc_out_path))
# 同时匹配 _1.fq / _1.fastq / _R1.fq / _R1.fastq
os.system(' ls *1.fq *1.fastq > raw_file_temp.txt')

raw_file = open ('raw_file_temp.txt')

for f1 in raw_file:
    f1 = f1.rstrip()
    if f1 == '':
        continue
    # 判断原始后缀
    if f1.endswith('.fastq'):
        ext = '.fastq'
    else:
        ext = '.fq'
    # 判断是 _R1/_R2 还是 _1/_2 命名
    m = re.search(r'_R1\.(fq|fastq)$', f1)
    if m:
        # 形如 sample_R1.fastq
        f2 = re.sub(r'_R1\.(fq|fastq)$', r'_R2.\1', f1)
        r1 = re.sub(r'_R1\.(fq|fastq)$', '', f1)
        r2 = re.sub(r'_R2\.(fq|fastq)$', '', f2)
        tag1 = '_R1'
        tag2 = '_R2'
    else:
        # 形如 sample_1.fastq（用 _1 结尾定位，避免误伤文件名中间的数字）
        f2 = re.sub(r'_1\.(fq|fastq)$', r'_2.\1', f1)
        r1 = re.sub(r'_1\.(fq|fastq)$', '', f1)
        r2 = re.sub(r'_2\.(fq|fastq)$', '', f2)
        tag1 = '_1'
        tag2 = '_2'
    print('now processing: %s and %s'%(r1,r2))
    os.system('''
cutadapt -q 20 -m 20 -j 20 -b file:{adapter} -n 2 -O 6 -e 0 -o {read1}_temp{tag1}.fq {read1}{tag1}{ext};
cutadapt -q 20 -m 20 -j 20 -b file:{adapter} -n 2 -O 6 -e 0 -o {read2}_temp{tag2}.fq {read2}{tag2}{ext};
cutadapt -q 20 -m 20 -j 20 -a "A{{30}}" -a "T{{30}}" -a "G{{30}}" -n 2 -O 30 -e 0 -o {read1}_trim{tag1}.fq {read1}_temp{tag1}.fq;
cutadapt -q 20 -m 20 -j 20 -a "A{{30}}" -a "T{{30}}" -a "G{{30}}" -n 2 -O 30 -e 0 -o {read2}_trim{tag2}.fq {read2}_temp{tag2}.fq;
repair.sh in1={read1}_trim{tag1}.fq in2={read2}_trim{tag2}.fq  out1={read1}_repair{tag1}.fq out2={read2}_repair{tag2}.fq;
rm {read1}_temp{tag1}.fq {read2}_temp{tag2}.fq {read1}_trim{tag1}.fq {read2}_trim{tag2}.fq;
fastqc -o ./{path1}/{path2} -t 5 {read1}_repair{tag1}.fq;
fastqc -o ./{path1}/{path2} -t 5 {read2}_repair{tag2}.fq;'''.format(adapter = adapter_fasta, read1 = r1, read2 = r2, ext = ext, tag1 = tag1, tag2 = tag2, path1 = trim_path, path2 = fastqc_out_path))

os.system('''rm raw_file_temp.txt''')

os.system('''ls *repair_*.fq > file_repair.txt''')

out_list = open ('file_repair.txt')

for f in out_list:
    f = f.rstrip()
    os.system('''mv {out} ./{path1}'''.format(out = f,path1 = trim_path ))

os.system('''rm file_repair.txt''')

print('cut & fastqc done!!!')