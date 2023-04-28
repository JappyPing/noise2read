# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:58:34
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-02-20 16:48:57
import os

raw_dir = "./raw/"

file_names = os.listdir(raw_dir)

for file_name in file_names:
    raw = os.path.join(raw_dir, file_name)
    no_ada_name = file_name.split(".fastq.gz")[0] + "_no_adapters.fastq.gz"
    os.system("cutadapt -a TGGAATTC -m 18 --max-n 0 -o %s %s" % (no_ada_name, raw))