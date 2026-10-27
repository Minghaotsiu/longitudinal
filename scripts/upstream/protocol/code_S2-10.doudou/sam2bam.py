# ！usr/bin/env python
# -*- coding:utf-8 -*-

import os
import re
import argparse

parser = argparse.ArgumentParser(description='将目录下的sam全部转为bam储存')
args = parser.parse_args()

os.system("ls *.sam > list_sam2bam.tmp.txt")
file_list = open('list_sam2bam.tmp.txt', 'r')
for i in file_list:
    i = re.sub('\n', '', i)
    prefix = re.sub('.sam', '', i)
    os.system("samtools view -@ 10 -bS %s > %s.bam" % (i, prefix))
file_list.close()
os.system("rm list_sam2bam.tmp.txt")
os.system("rm *.sam")

