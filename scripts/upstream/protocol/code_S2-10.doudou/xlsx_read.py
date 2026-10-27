#!/usr/bin/env python
#coding=utf-8

## python2 xls_read.py file.xls
## python xlsx_read.py file.xlsx

import sys
import re
import openpyxl

# 打开excel文件
data = openpyxl.load_workbook(sys.argv[1])
# 查看sheet，根据sheet名字来获取excel中的sheet
#print("sheets：" )
#print(data.sheetnames)
#name_or_index = input('Enter the sheet name or index(n or i):\n');
#if name_or_index == 'n':
#    by_name = input('Enter the sheet name:\n');
#    ws = data[by_name]
#elif name_or_index == 'i':
#    by_index = input('Enter the sheet index:\n');
#    sn = data.sheetnames[int(by_index)]
#    ws = data[sn]
by_index = 0
sn = data.sheetnames[int(by_index)]
ws = data[sn]
#else:
#    print('Input Error!')

# 某一行数据
#txt = input('Do you want to output the sheet to a txt file(y or n):\n')
txt = 'y'
if txt == 'y':
    prefix=re.sub('.xlsx','',sys.argv[1])
    file=open(prefix+'.txt','w')
    for r in ws.rows:
        line = ''
        for c in r:
            if (str(c.value)!='None'):
                line=line+str(c.value)+'\t'
        line=line.rstrip()+'\n'
        file.write(line)
    file.close()
elif txt == 'n':
    for r in ws.rows:
        line = ''
        for c in r:
            line=line+c.value+'\t'
        line=line.rstrip()+'\n'
        print(line)

