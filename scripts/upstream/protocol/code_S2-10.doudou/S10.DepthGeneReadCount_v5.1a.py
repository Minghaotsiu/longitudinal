#!/usr/bin/env python
# coding=utf-8

import os
import re
import numpy as np
import openpyxl
import argparse
import random
import pysam
import subprocess

# Argument parsing
parser = argparse.ArgumentParser(description='Calculate gene counts at the same sequencing depth.')
parser.add_argument('-s', '--sam', type=str, nargs='+', help='A BAM file(s) with suffix ".bam" or a list of BAM files.')
parser.add_argument('-g', '--gtf', type=str, default=os.environ.get('SOYBEAN_GTF'), help='GTF file for htseq-count.')
parser.add_argument('-mt', '--mtgene', type=str, default=os.environ.get('MT_GENE_FILE'), help='A list of mitochondrial gene IDs/names.')
parser.add_argument('-rp', '--rpgene', type=str, default=os.environ.get('RP_GENE_FILE'), help='A list of ribosomal protein gene IDs/names.')
parser.add_argument('-d', '--depth', type=str, required=True, help='Sequencing depth or "raw".')
parser.add_argument('-o', '--output', type=str, required=True, help='Output Excel file.')

args = parser.parse_args()
if not args.gtf:
    raise SystemExit("Set SOYBEAN_GTF or pass --gtf.")

# Function to perform downsampling on BAM file
def DepthGene(bam, d, depth_str):
    # Ensure the BAM file is indexed
    if not os.path.exists(bam + '.bai'):
        pysam.index(bam)
    
    # Open BAM file
    with pysam.AlignmentFile(bam, "rb") as bam_in:
        # Get number of aligned reads
        num_reads = bam_in.count(until_eof=True)
        
        # Check if the desired depth is larger than the number of reads
        if d > num_reads:
            print(f"Warning: Desired depth ({d}) is larger than the number of reads ({num_reads}). Using all reads.")
            d = num_reads
        
        # Collect all read IDs
        read_ids = [read.query_name for read in bam_in]
        
        # Reset iterator and sample read IDs
        bam_in.reset()
        sample_ids = random.sample(read_ids, k=min(d, len(read_ids)))
        
        # Create output BAM file
        out_bam = bam_in.filename.decode().rsplit('.', 1)[0] + f'_{depth_str}depth.bam'
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for read in bam_in:
                if read.query_name in sample_ids:
                    bam_out.write(read)
    return out_bam


# Function to process the output of htseq-count
def GeneReadCount(htseq_out, mtgene, rpgene):
    with open(htseq_out, 'r') as file:
        prefix = re.sub('_counts.txt', '', htseq_out)
        with open(prefix+'_GeneReadCount.txt', 'w') as file_count:
            ncount = 0
            ngene = 0
            mtcount = 0
            rpcount = 0
            first_nonzero_gene = None
            for line in file:
                if ('processed' in line) or ('__' in line):
                    continue
                else:
                    line2 = re.sub('\n', '', line)
                    l = line2.split('\t')
                    if (l[1] != '0'):
                        ncount += int(l[1])
                        ngene += 1
                        file_count.write(l[0] + '\t' + l[1] + '\n')
                        if l[0] in mtgene:
                            mtcount += int(l[1])
                        if l[0] in rpgene:
                            rpcount += int(l[1])
                        if first_nonzero_gene is None:
                            first_nonzero_gene = l[0]
                    else:
                        print(f"Warning: Zero count for {l[0]}")

            # Avoid division by zero
            if ncount == 0:
                mtprop = 0
                rpprop = 0
            else:
                mtprop = mtcount / ncount * 100
                rpprop = rpcount / ncount * 100

            file_count.write('nCount: ' + str(ncount) + '\n')
            file_count.write('nGene: ' + str(ngene) + '\n')
            file_count.write('mtCount: ' + str(mtcount) + '\n')
            file_count.write('mtProp: %.2f%%\n' % mtprop)
            file_count.write('rpCount: ' + str(rpcount) + '\n')
            file_count.write('rpProp: %.2f%%' % rpprop)

            # Additional logging
            print(f"Processed {ngene} genes with non-zero counts.")
            print(f"Total count: {ncount}")
            if first_nonzero_gene is not None:
                print(f"First non-zero count gene: {first_nonzero_gene}")

# Load gene lists
mtfile = open(args.mtgene, 'r')
mtgene = [i.rstrip() for i in mtfile.readlines()]
mtfile.close()
rpfile = open(args.rpgene, 'r')
rpgene = [i.rstrip() for i in rpfile.readlines()]
rpfile.close()

# Main script
list_bam = []
if args.sam:
    for ss in args.sam:
        if ss.endswith('.bam'):
            list_bam.append(ss)
        else:
            with open(ss, 'r') as file:
                for i in file:
                    list_bam.append(i.rstrip())
    list_bam = sorted(list(set(list_bam)))
else:
    os.system("ls *.bam > list_bam_htseqS10.tmp.txt")
    with open('list_bam_htseqS10.tmp.txt', 'r') as file:
        for i in file:
            list_bam.append(i.rstrip())
    os.system("rm list_bam_htseqS10.tmp.txt")

d = str(args.depth)
downsample = False
depth_str = 'raw'
if d == 'raw':
    pass
else:
    downsample = True
    if 'M' in d or 'm' in d:
        depth_str = re.sub('m', 'M', d)
        d = re.sub('[Mm]', '', d)
        d = int(float(d) * 1000000)
    else:
        d = int(d)
        depth_str = str(d / 1000000) + 'M'

# Initialize Excel workbook
output = args.output
wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'Sheet1'
row = 5

for bam in list_bam:
    print(f'Now processing {bam}')

    prefix = re.sub('.sorted.bam', '', bam)
    prefix = re.sub('.bam', '', prefix)
    
    if downsample:
        DepthGene_out = DepthGene(bam, d, depth_str)
        htseq_out1 = prefix + '_' + depth_str + '_union_counts.txt'
        GeneReadCount_out1 = prefix + '_' + depth_str + '_union_GeneReadCount.txt'
        output_name = output + '_' + depth_str + '_GeneReadCounts.xlsx'
    else:
        DepthGene_out = bam
        htseq_out1 = prefix + '_union_counts.txt'
        GeneReadCount_out1 = prefix + '_union_GeneReadCount.txt'
        output_name = output + '_GeneReadCounts.xlsx'

    # Execute: htseq-count
    # Check if htseq-count is available
    try:
        subprocess.run(["which", "htseq-count"], check=True)
    except subprocess.CalledProcessError:
        print("htseq-count not found. Please ensure it is installed and in your PATH.")
        exit(1)

    os.system(f"htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m union {DepthGene_out} {args.gtf} > {htseq_out1}")

    # Calculate number of genes
    result1 = subprocess.run(["wc", "-l", htseq_out1], capture_output=True, text=True).stdout
    num1 = result1.split(' ')[0]
    total_gene_num1 = int(num1) - 5

    # Execute: calculate nCount & nGene
    GeneReadCount(htseq_out1, mtgene, rpgene)

    with open(GeneReadCount_out1, 'r') as file_GeneReadCount1:
        cont1 = file_GeneReadCount1.readlines()
        ncount1 = cont1[-6].rstrip().split(': ')[1]
        ngene1 = cont1[-5].rstrip().split(': ')[1]
        mtcount1 = cont1[-4].rstrip().split(': ')[1]
        mtprop1 = cont1[-3].rstrip().split(': ')[1]
        rpcount1 = cont1[-2].rstrip().split(': ')[1]
        rpprop1 = cont1[-1].rstrip().split(': ')[1]

    ws['A' + str(row)] = prefix
    ws['B' + str(row)] = prefix.split('_')[0]
    ws['C' + str(row)] = ncount1
    ws['D' + str(row)] = ngene1
    ws['E' + str(row)] = mtcount1
    ws['F' + str(row)] = mtprop1
    ws['G' + str(row)] = rpcount1
    ws['H' + str(row)] = rpprop1

    print(f'{prefix} is finished!')
    row += 1

ws['A1'] = 'Depth reads'
ws['B1'] = depth_str
ws['A2'] = 'Total gene num'
ws['B2'] = str(total_gene_num1)
ws['C3'] = '基因总表达量'
ws['D3'] = '基因数'
ws['E3'] = '线粒体表达量'
ws['F3'] = '线粒体表达量占比'
ws['G3'] = '核糖蛋白表达量'
ws['H3'] = '核糖蛋白表达量占比'
ws['A4'] = 'File'
ws['B4'] = 'Sample'
ws['C4'] = 'nCount(union)'
ws['D4'] = 'nGene(union)'
ws['E4'] = 'mtCount(union)'
ws['F4'] = 'mtProp(union)'
ws['G4'] = 'rpCount(union)'
ws['H4'] = 'rpProp(union)'

wb.save(output_name)
