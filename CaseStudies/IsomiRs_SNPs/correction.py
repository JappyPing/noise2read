# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-05-03 16:55:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 17:03:11
import sys
import os

def main(raw_dir, output_dir, correct_dir):
    file_names = os.listdir(raw_dir)

    for file_name in file_names:
        raw = os.path.join(raw_dir, file_name)
        isomir_prfix = file_name.split(".MI_")[0]
        os.system("noise2read -m correction -i %s -d %s" % (raw, output_dir + "/" + isomir_prfix + "/"))
        os.system('mv {} {} {}'.format(output_dir + "/" + isomir_prfix + "/", '*corrected.fastq', correct_dir))

if __name__ == '__main__':
    raw_dir = sys.argv[1]
    output_dir = sys.argv[2]
    correct_dir = sys.argv[3]

    main(raw_dir, output_dir, correct_dir)
