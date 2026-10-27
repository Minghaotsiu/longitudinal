# ！usr/bin/env python
# -*- coding:utf-8 -*-

import re
import argparse
import os

parser = argparse.ArgumentParser(description='给出指定gene_id对应的gene_name,或将指定文件中的gene_id替换成gene_name: \
                                              python gene_id2name.py \
                                              -g gtf \
                                              -i id or file \
                                              -t type:i or f')

parser.add_argument('-g', '--gtf', type=str, default=os.environ.get('SOYBEAN_GTF'), help='path to the soybean GTF file')
parser.add_argument('-i', '--input_str', type=str, nargs='+', help='gene_id or file with gene_id in it')
parser.add_argument('-t', '--input_type', type=str, required=True, help='you should tell the type of input: i(id) or f(file)')

args = parser.parse_args()

gtf = args.gtf
if not gtf:
    raise SystemExit("Set SOYBEAN_GTF or pass --gtf.")
inp = args.input_str
typ = args.input_type

# save gene_id:gene_name to dictg={}
dictg = {}
dictg2 = {}
fileg = open(gtf, 'r')
contg = fileg.readlines()
pattern_id = r'gene_id "(.*?)";'
pattern_name = r'gene_name "(.*?)";'
row=0
for line in contg:
    row+=1
    if line.startswith('#'):
        continue
    else:
        ID = re.search(pattern_id, line)
        NAME = re.search(pattern_name, line)
        if not ((ID is None) or (NAME is None)):
            gene_id = ID.group(1)
            gene_id2 = gene_id.split('.')[0]
            gene_name = NAME.group(1)
            dictg[gene_id] = gene_name
            dictg2[gene_id2] = gene_name
fileg.close()

# gene_id to gene_name
if typ == 'i':
    for i in inp:
        if i in dictg:
            print(i+'\t'+dictg[i])
        elif i in dictg2:
            print(i+'\t'+dictg2[i])
        else:
            print('The given gene_id is not in the input gtf!')
elif typ == 'f':
    for f in inp:
        suffix = f.split('.')[-1]
        prefix = re.sub('.'+suffix, '', f)
        file_in = open(f, 'r')
        file_out = open(prefix+'.gn.'+suffix, 'w')
        cont_in = file_in.read()
        cont_out = cont_in
        for k, v in dictg.items():
            cont_out = re.sub(k, v, cont_out)
        for k2, v2 in dictg2.items():
            cont_out = re.sub(k2, v2, cont_out)
        file_out.write(cont_out)
        file_in.close()
        file_out.close()

