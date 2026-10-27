# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 13:56:04 2020

@author: cyclopenta
"""

import sys
import os
import re
import datetime

#define parameters
source_star_path = os.environ.get("STAR_INDEX")
gtf_file_path = os.environ.get("SOYBEAN_GTF")
if not source_star_path or not gtf_file_path:
    raise SystemExit("Set STAR_INDEX and SOYBEAN_GTF before running STAR.")
now_time = datetime.datetime.now()
time_str = now_time.strftime('%Y-%m-%d-%H_%M_%S')
star_out_path = 'star_out_'+time_str
log_out_path = 'log_final_out'
os.system('''mkdir -p ./{path1}/{path2}'''.format(path1 = star_out_path, path2 = log_out_path))
# 同时匹配 _1_repair_1 / _R1_repair_R1 且 .fq/.fastq 均可
os.system(' ls *_repair_1.fq *_repair_1.fastq *_repair_R1.fq *_repair_R1.fastq > clean_file_star_temp.txt')

raw_file = open ('clean_file_star_temp.txt')

for f1 in raw_file:
    f1 = f1.rstrip()
    if f1 == '' or 'No such file' in f1:
        continue
    # 判断后缀
    if f1.endswith('.fastq'):
        ext = '.fastq'
    else:
        ext = '.fq'
    # 判断命名：R1 版 还是 1 版
    if re.search(r'_repair_R1\.(fq|fastq)$', f1):
        f2 = re.sub(r'_repair_R1\.(fq|fastq)$', r'_repair_R2.\1', f1)
        out_prefix = re.sub(r'_repair_R1\.(fq|fastq)$', '', f1)
    else:
        f2 = re.sub(r'_repair_1\.(fq|fastq)$', r'_repair_2.\1', f1)
        out_prefix = re.sub(r'_repair_1\.(fq|fastq)$', '', f1)
    r1 = re.sub(r'\.(fq|fastq)$', '', f1)
    r2 = re.sub(r'\.(fq|fastq)$', '', f2)
    print('now processing: %s and %s'%(r1,r2))
    os.system('''
STAR --runThreadN 10 --genomeDir {starpath} --outFileNamePrefix {prefix} --sjdbGTFfile {gtf} --outSAMunmapped Within --readFilesIn {read1}{ext} {read2}{ext};
mv {prefix}Log.final.out ./{path1}/{path2};
mv {prefix}Aligned.out.sam {prefix}Log.out {prefix}Log.progress.out {prefix}SJ.out.tab {prefix}_STARgenome ./{path1};
'''.format(starpath=source_star_path, read1=r1, read2=r2, ext=ext, prefix=out_prefix, gtf=gtf_file_path, path1=star_out_path, path2=log_out_path))

os.system('''rm clean_file_star_temp.txt''')

print('star map done!!!')
