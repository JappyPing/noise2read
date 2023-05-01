# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:56:24
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-30 14:08:09
import os

raw_dir = "../D22_D31/raw_no_adapters/"

file_names = os.listdir(raw_dir)

for file_name in file_names:
    raw = os.path.join(raw_dir, file_name)
    isomir_prfix = file_name.split(".MI_")[0]
    os.system("python IsoMiRmap.py %s --p %s" % (raw, isomir_prfix))

correct_dir = "../D22_D31/corrected/"

file_names = os.listdir(correct_dir)

for file_name in file_names:
    raw = os.path.join(correct_dir, file_name)
    isomir_prfix = file_name.split(".MI_")[0]
    os.system("python IsoMiRmap.py %s --p %s" % (raw, isomir_prfix))