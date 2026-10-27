import os
import re

# python raw_fastqc.py

os.makedirs("fastqc_result", exist_ok=True)
os.system("ls *.fq *.fastq > list_seq.txt")
file = open('list_seq.txt', 'r')
for i in file:
    fq = re.sub('\n', '', i)
    os.system("fastqc -o ./fastqc_result -t 5 %s" %(fq))
file.close()
os.system("rm list_seq.txt")
