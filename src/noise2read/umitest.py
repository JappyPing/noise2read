'''
@File    :   umitest.py
@Time    :   2022/10/25 23:21:30
@Author  :   Pengyao Ping
@Version :   1.0
@Contact :   Ping.Pengyao@gmail.com
@desc    : 
'''

# from collections import Counter
import collections
import math
from Bio import SeqIO
# import os
# import editdistance
# import xlsxwriter
# from noise2read.utils import parse_file_type

def parse_file_type(data_set):
    items = data_set.split(".")
    ext = items[-1]
    if ext == 'fa' or ext == 'fasta':
        f_type = 'fasta'
    elif ext == 'fq' or ext == 'fastq':
        f_type = 'fastq'
    elif ext == 'gz':
        f_type = items[-2] + '.' + ext
    return f_type

# raw_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/umi24/raw/non_umi_SRR1543964.fastq'
# ground_truth_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/umi24/true/non_umi_SRR1543964.fastq'
def data_analysis(raw_dataset, ground_truth_dataset):
    true_records = SeqIO.index(ground_truth_dataset, parse_file_type(ground_truth_dataset))
    error_records = SeqIO.index(raw_dataset, parse_file_type(raw_dataset))
    # correct_records = SeqIO.index(correct_dataset, parse_file_type(correct_dataset))
    total_reads_num = 0
    # for calculating entropy
    correct_errfree_seqs = []
    correct_err_seqs = []
    raw_errfreee_seqs = []
    raw_err_seqs = []

    true_seqs_lst = []
    raw_seqs_lst = []
    # correct_seqs_lst = []
    true_umi2seqs = {}
    raw_umi2seqs = {}
    # correct_umi2seqs = {}

    true_seq2umi = {}
    for name in error_records:
        true_seq = str(true_records[name].seq)
        raw_seq = str(error_records[name].seq)
        # correct_seq = str(correct_records[name].seq) 

        umi_base = str(true_records[name].description).split('//')[0]
        umi = umi_base.split(':')[1]
        true_umi2seqs.setdefault(umi, []).append(true_seq)
        raw_umi2seqs.setdefault(umi, []).append(raw_seq)
        # correct_umi2seqs.setdefault(umi, []).append(correct_seq)

        true_seq2umi.setdefault(true_seq, []).append(umi)

        true_seqs_lst.append(true_seq)     
        raw_seqs_lst.append(raw_seq)
        # correct_seqs_lst.append(correct_seq)

    raw_read2count = collections.Counter(raw_seqs_lst)
    true_read2count = collections.Counter(true_seqs_lst)
    # correct_read2count = collections.Counter(correct_seqs_lst)

    num = 0
    # read_num = 0
    for seq in true_seq2umi:
        unique_umi = set(true_seq2umi[seq])
        if len(unique_umi) > 1:
            num += 1
            print(unique_umi)
    print(num)

if __name__ == '__main__':

    # raw_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/umi24/raw/non_umi_SRR1543964.fastq'
    # ground_truth_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/umi24/true/non_umi_SRR1543964.fastq'

    raw_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/test12/raw/non_umi_SRR1543964.fastq'
    ground_truth_dataset = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/test12/true/non_umi_SRR1543964.fastq'

    # data_analysis(raw_dataset, ground_truth_dataset)

    umi_data = '/home/pping/Data/Repo/data/noise2read_data/umi-data/group1/umi1224/umi_SRR1543964_corrected.fastq'

    umi_records = SeqIO.index(umi_data, parse_file_type(umi_data))
    umi_seqs = []
    for name in umi_records:
        umi_seqs.append(str(umi_records[name].seq))
    
    read2counts = collections.Counter(umi_seqs)
    print(read2counts)