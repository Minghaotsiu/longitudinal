#!/usr/bin/env python
# coding=utf-8

import os
import re
import argparse
from multiprocessing import Pool
from tqdm import tqdm

def calculate_mapped_bases(depth_file):
    """Calculate the total number of mapped bases from a depth file."""
    total_bases = 0
    with open(depth_file, 'r', encoding='utf-8') as file_depth:
        for line in file_depth:
            _, _, depth = line.strip().split('\t')
            total_bases += int(depth)
    return total_bases

def exon_intron_analysis(sam_file, exon_bed, gene_bed, output_prefix):
    """Perform exon and intron analysis on the given SAM file."""
    # Remove the 'Aligned.out.sam' suffix to create the prefix
    prefix = re.sub(r'Aligned\.out\.sam$', '', sam_file)

    # Increase the number of threads for samtools commands
    num_threads = 10

    # Generate BAM files and sort them
    os.system(f"samtools view -@ {num_threads} -bS {sam_file} > {prefix}.bam")
    os.system(f"samtools sort -@ {num_threads} -m 32G {prefix}.bam -o {prefix}.sorted.bam")

    # Generate the total depth file
    os.system(f"samtools depth {prefix}.sorted.bam > {prefix}.bam.depth")

    # Generate exon and gene region depth files
    os.system(f"samtools depth -b {exon_bed} {prefix}.sorted.bam > {prefix}_exon.bam.depth")
    os.system(f"samtools depth -b {gene_bed} {prefix}.sorted.bam > {prefix}_gene.bam.depth")

    # Write results to a summary file
    with open(f"{prefix}_ExonIntron.txt", 'w', encoding='utf-8') as file_ExonIntron:
        file_ExonIntron.write(f"Sample name: {prefix}\n")

        total_bases = calculate_mapped_bases(f"{prefix}.bam.depth")
        file_ExonIntron.write(f"Total mapped bases: {total_bases}\n")

        exon_bases = calculate_mapped_bases(f"{prefix}_exon.bam.depth")
        p_exon = exon_bases / float(total_bases)
        file_ExonIntron.write(f"Exon bases: {exon_bases}\n")
        file_ExonIntron.write(f"% of bases mapped to exon: {p_exon:.2%}\n")

        gene_bases = calculate_mapped_bases(f"{prefix}_gene.bam.depth")
        p_gene = gene_bases / float(total_bases)
        p_intron = p_gene - p_exon
        file_ExonIntron.write(f"Gene bases: {gene_bases}\n")
        file_ExonIntron.write(f"% of bases mapped to gene: {p_gene:.2%}\n")
        file_ExonIntron.write(f"% of bases mapped to intron: {p_intron:.2%}\n")

    # Clean up intermediate files
    os.system(f"rm {prefix}.bam {prefix}.bam.depth {prefix}_exon.bam.depth {prefix}_gene.bam.depth")
    print(f"{sam_file} is finished!")

def process_sam_file(args):
    """Wrapper function to process each SAM file using the exon_intron_analysis function."""
    sam_file, exon_bed, gene_bed, output_prefix = args
    exon_intron_analysis(sam_file, exon_bed, gene_bed, output_prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate exon and intron alignment rates for all Aligned.out.sam files.")
    parser.add_argument('-e', '--exon', type=str, default=os.environ.get('SOYBEAN_EXON_BED'), help='Path to exon.bed')
    parser.add_argument('-g', '--gene', type=str, default=os.environ.get('SOYBEAN_GENE_BED'), help='Path to gene.bed')

    args = parser.parse_args()
    if not args.exon or not args.gene:
        raise SystemExit("Set SOYBEAN_EXON_BED and SOYBEAN_GENE_BED.")

    exon_bed = args.exon
    gene_bed = args.gene

    # List all Aligned.out.sam files
    os.system("ls *Aligned.out.sam > list_Aligned_out_sam_inex.tmp.txt")
    with open('list_Aligned_out_sam_inex.tmp.txt', 'r', encoding='utf-8') as file_list_sam:
        sam_files = [sam_file.strip() for sam_file in file_list_sam]

    # Define the number of processes to use (based on available cores)
    num_processes = 4  # Use 4 processes to parallelize the work

    # Prepare arguments for the pool.map method
    arguments = [(sam_file, exon_bed, gene_bed, sam_file) for sam_file in sam_files]

    # Use multiprocessing to speed up the processing with progress bar
    with Pool(processes=num_processes) as pool:
        for _ in tqdm(pool.imap_unordered(process_sam_file, arguments), total=len(arguments)):
            pass

    # Clean up temporary file
    os.system("rm list_Aligned_out_sam_inex.tmp.txt")
