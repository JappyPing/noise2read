# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:56:24
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 20:37:50
import os

raw_dir = "../D22_D31/raw_no_adapters/"

file_names = os.listdir(raw_dir)

for file_name in file_names:
    raw = os.path.join(raw_dir, file_name)
    isomir_prfix = file_name.split(".MI_")[0]
    os.system("python IsoMiRmap.py %s --p %s" % (raw, isomir_prfix))

raw_result_dir = "../isomimap_result/raw"
if not os.path.exists(raw_result_dir):
    os.makedirs(raw_result_dir)
os.system('mv {}{} {}'.format("./", '*.txt', raw_result_dir))
os.system('mv {}{} {}'.format("./", '*.html', raw_result_dir))
os.system('mv {}{} {}'.format("./", '*.gff3', raw_result_dir))

correct_dir = "../D22_D31/corrected"

cor_file_names = os.listdir(correct_dir)

for cor_file_name in cor_file_names:
    correct = os.path.join(correct_dir, cor_file_name)
    cor_isomir_prfix = cor_file_name.split("corrected.fastq")[0]
    os.system("python IsoMiRmap.py %s --p %s" % (correct, cor_isomir_prfix))

correct_result_dir = "../isomimap_result/corrected"
if not os.path.exists(correct_result_dir):
    os.makedirs(correct_result_dir)
os.system('mv {}{} {}'.format("./", '*.txt', correct_result_dir))
os.system('mv {}{} {}'.format("./", '*.html', correct_result_dir))
os.system('mv {}{} {}'.format("./", '*.gff3', correct_result_dir))