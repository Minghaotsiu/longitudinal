#!/usr/bin/env python  
# -*- coding=utf-8 -*-  
  
import os  
import re  
import openpyxl  
import argparse  
  
parser = argparse.ArgumentParser(description='整合各样本的ExonIntron信息（对目录下所有*_ExonIntron.txt操作）')  
parser.add_argument('-o', '--output', type=str, required=True, help='输出文件名前缀')  
args = parser.parse_args()  
  
wb = openpyxl.Workbook()  
ws = wb.active  
ws.title = 'Sheet1'  
ws.append(['File', 'Sample', '% of bases mapped to exon', '% of bases mapped to intron', '% of bases mapped to gene', 'Exon bases', 'Gene bases', 'Total mapped bases'])  
  
data = {}  # 使用字典来存储每个文件的数据  
  
# 获取所有 *_ExonIntron.txt 文件  
files = [f for f in os.listdir('.') if f.endswith('_ExonIntron.txt')]  
  
for file_name in files:  
    sample_name = None  
    total_mapped_bases = None  
    exon_bases = None  
    p_of_bases_mapped_to_exon = None  
    gene_bases = None  
    p_of_bases_mapped_to_gene = None  
    p_of_bases_mapped_to_intron = None  
  
    with open(file_name, 'r') as file_ExonIntron:  
        for line in file_ExonIntron:  
            key, value = line.strip().split(': ', 1)  
            if key == 'Sample name':  
                sample_name = value  
            elif key == 'Total mapped bases':  
                total_mapped_bases = value  
            elif key == 'Exon bases':  
                exon_bases = value  
            elif key == '% of bases mapped to exon':  
                p_of_bases_mapped_to_exon = value  
            elif key == 'Gene bases':  
                gene_bases = value  
            elif key == '% of bases mapped to gene':  
                p_of_bases_mapped_to_gene = value  
            elif key == '% of bases mapped to intron':  
                p_of_bases_mapped_to_intron = value  
  
    if all([sample_name, total_mapped_bases, exon_bases, p_of_bases_mapped_to_exon, gene_bases, p_of_bases_mapped_to_gene, p_of_bases_mapped_to_intron]):  
        # 提取样本名的第一部分作为简化样本名  
        simplified_sample_name = sample_name.split('_')[0]  
        data[file_name] = {  
            'Sample': simplified_sample_name,  
            '% of bases mapped to exon': p_of_bases_mapped_to_exon,  
            '% of bases mapped to intron': p_of_bases_mapped_to_intron,  
            '% of bases mapped to gene': p_of_bases_mapped_to_gene,  
            'Exon bases': exon_bases,  
            'Gene bases': gene_bases,  
            'Total mapped bases': total_mapped_bases  
        }  
  
# 将数据写入 Excel  
row = 2  
for file_name, info in data.items():  
    ws.append([file_name, info['Sample'], info['% of bases mapped to exon'], info['% of bases mapped to intron'],  
               info['% of bases mapped to gene'], info['Exon bases'], info['Gene bases'], info['Total mapped bases']])  
  
wb.save(args.output + '_ExonIntron_MAP.xlsx')