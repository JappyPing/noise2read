# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-07 16:47:06

from typing import Counter
import numpy as np
# import os
from sklearn.preprocessing import StandardScaler
from noise2read.encoding import EncodeScheme
from mpire import WorkerPool
from tqdm import tqdm
from noise2read.utils import *
import itertools

class Reads2Vectors():
    def __init__(self, logger, config, edit_dis):
        self.logger = logger
        self.edit_dis = edit_dis
        self.config = config

    def read2features(self, shared_objects, i):
        ES, reads_lst1, reads_lst2, other_features = shared_objects
        features = []
        # pd_fe = ES.descriptors("PairDistance", reads_lst[i])
        # features.extend(pd_fe)    
        ###########################################################################
        # features.append(read_count_lst[i]) 
        # features.append(err_pos_lst[i]) 
        # freq_pos.append([ori_seq_freq, err_pos])              
        ###########################################################################
        ft_fea1 = ES.descriptors("FourierTransform", reads_lst1[i])
        cg_fea1 = ES.descriptors("ChaosGame", reads_lst1[i])
        entropy_fea1 = ES.descriptors("Entropy", reads_lst1[i])
        fs_fea1 = ES.descriptors("FickettScore", reads_lst1[i])

        features.extend(ft_fea1)
        features.extend(cg_fea1)
        features.extend(entropy_fea1)
        features.extend(fs_fea1)
        # print(len(features))
        # ft_fea2 = ES.descriptors("FourierTransform", reads_lst2[i])
        # cg_fea2 = ES.descriptors("ChaosGame", reads_lst2[i])
        # # entropy_fea2 = ES.descriptors("Entropy", reads_lst2[i])
        # fs_fea2 = ES.descriptors("FickettScore", reads_lst2[i])

        # features.extend(ft_fea2)
        # features.extend(cg_fea2)
        # # features.extend(entropy_fea2)
        # features.extend(fs_fea2)
        # self.logger.debug(len(features))
        atomic_fea1 = ES.descriptors("atomic_number", reads_lst1[i])
        atomic_fea2 = ES.descriptors("atomic_number", reads_lst2[i])
        # atomic_fea1 = ES.descriptors("binary", reads_lst1[i])
        # atomic_fea2 = ES.descriptors("binary", reads_lst2[i])
        # print(int(other_features[i][0]))
        features.extend(atomic_fea1)
        features.extend(atomic_fea2)
        # onehot_fea = ES.descriptors("OneHot", reads_lst1[i])
        # features.extend(onehot_fea)
        features.extend(other_features[i])
        # self.logger.debug(f'FourierTransform: {len(ft_fea)}, ChaosGame: {len(cg_fea)}, Entropy: {len(entropy_fea)}, FickettScore: {len(fs_fea)}')
        return features 

    def high_all_in_one_embedding(self, genuine_df, negative_df, new_negative_df, ambiguous_df):
        self.logger.info("  Embedding genuine and high ambiguous data.")
        genuine_reads_lst1 = []
        negtive_reads_lst1 = []
        ambiguous_reads_lst1 = []

        genuine_reads_lst2 = []
        negtive_reads_lst2 = []
        ambiguous_reads_lst2 = []

        genuine_reads_features = []
        negtive_reads_features = []
        ambiguous_reads_features = [] 
        # if self.config.read_type == "DNA":
        #     error_tyes = ["A-G", "G-A", "A-T", "T-A", "A-C", "C-A", "G-T", "T-G", "G-C", "C-G", "T-C", "C-T", "T-N", "N-T", "C-N", "N-C", "A-N", "N-A", "G-N", "N-G"]
        #     base_lst = ['A', 'C', 'G', 'T', 'N']
        #     kmers = [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]
        #     # kmers3 = ['ANA', 'CNA', 'TNA', 'GNA', 'ANC', 'CNC', 'TNC', 'GNC', 'ANG', 'CNG', 'TNG', 'GNG', 'ANT', 'CNT', 'TNT', 'GNT'] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]
        #     # kmers = kmers2 + kmers3
        # elif self.config.read_type == "RNA":
        #     base_lst = ['A', 'C', 'G', 'U', 'N']
        #     kmers = [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]    
        base_lst = ['A', 'C', 'G', 'T', 'N']
        if self.config.read_type == "DNA":
            error_tyes = ["A-G", "G-A", "A-T", "T-A", "A-C", "C-A", "G-T", "T-G", "G-C", "C-G", "T-C", "C-T", "T-X", "X-T", "C-X", "X-C", "A-X", "X-A", "G-X", "X-G", "X-N", "N-X", 'A-N', 'T-N','G-N','C-N','N-A','N-T', 'N-C', 'N-G']
            kmers = ['NX', 'XN', 'XA', 'XC', 'XG', 'XT', 'AX', 'CX', 'GX', 'TX'] + [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)] + ['AXA', 'CXA', 'TXA', 'GXA', 'AXC', 'CXC', 'TXC', 'GXC', 'AXG', 'CXG', 'TXG', 'GXG', 'AXT', 'CXT', 'TXT', 'GXT', 'NXA', 'NXC', 'NXG', 'NXT', 'NXN','AXN','CXN','GXN','TXN']
            # kmers3 = ['ANA', 'CNA', 'TNA', 'GNA', 'ANC', 'CNC', 'TNC', 'GNC', 'ANG', 'CNG', 'TNG', 'GNG', 'ANT', 'CNT', 'TNT', 'GNT'] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]
            # kmers = kmers2 + kmers3
        elif self.config.read_type == "RNA":
            # error_tyes = ["A-G", "G-A", "A-U", "U-A", "A-C", "C-A", "G-U", "U-G", "G-C", "C-G", "U-C", "C-U", "U-N", "N-U", "C-N", "N-C", "A-N", "N-A", "G-N", "N-G"]
            error_tyes = ["A-G", "G-A", "A-T", "T-A", "A-C", "C-A", "G-T", "T-G", "G-C", "C-G", "T-C", "C-T", "T-X", "X-T", "C-X", "X-C", "A-X", "X-A", "G-X", "X-G", "X-N", "N-X", 'A-N', 'T-N','G-N','C-N','N-A','N-T', 'N-C', 'N-G']
            kmers = ['NX', 'XN', 'XA', 'XC', 'XG', 'XT', 'AX', 'CX', 'GX', 'TX'] + [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)] + ['AXA', 'CXA', 'TXA', 'GXA', 'AXC', 'CXC', 'TXC', 'GXC', 'AXG', 'CXG', 'TXG', 'GXG', 'AXT', 'CXT', 'TXT', 'GXT', 'NXA', 'NXC', 'NXG', 'NXU', 'NXN','AXN','CXN','GXN','UXN']

        error_tye_priors = {}
        err_tye_base_prior = 0.1
        for it in error_tyes:
            error_tye_priors[it] = err_tye_base_prior
            err_tye_base_prior += 0.01

        self.logger.debug(len(kmers))
        total_err_tyes = []
        total_err_kmers = []
        for idx, row in genuine_df.iterrows():
            total_err_tyes.append(row['ErrorTye'])
            total_err_kmers.append(row['StartErrKmer'])
            total_err_kmers.append(row['EndErrKmer'])

        for idx, row in ambiguous_df.iterrows():
            total_err_tyes.append(row['ErrorTye'])
            total_err_kmers.append(row['StartErrKmer'])
            total_err_kmers.append(row['EndErrKmer'])

        for idx, row in new_negative_df.iterrows():
            read = row['StartRead']
            total_err_tyes.append(row['ErrorTye'])
            total_err_kmers.append(row['StartErrKmer'])
            total_err_kmers.append(row['EndErrKmer'])
            
        total_err_kmers_count = len(total_err_kmers)
        total_err_tyes_count = len(total_err_tyes)
        
        err_kmers2count = Counter(total_err_kmers)
        err_tyes2count = Counter(total_err_tyes)
        self.logger.debug(len(err_kmers2count.keys()))
        kmer_keys = err_kmers2count.keys()

        not_exist_kmer = set(kmers) - set(kmer_keys)
        if not_exist_kmer:
            for mer in not_exist_kmer:
                err_kmers2count[mer] = 0

        kmers_priors = {}
        base_prior = 0.1
        for item in kmer_keys:
            kmers_priors[item] = base_prior
            base_prior += 0.01        
        self.logger.debug(err_kmers2count)
        self.logger.debug(err_tyes2count)

        for idx, row in genuine_df.iterrows():
            genuine_reads_lst1.append(row['StartRead'])
            genuine_reads_lst2.append(row['EndRead'])

            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            # self.logger.debug(row['StartRead'])
            # self.logger.debug(row['EndRead'])
            # self.logger.debug(f'{cur_kmer1}, {cur_kmer2}')
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            genuine_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
        
        for idx, row in new_negative_df.iterrows():
            negtive_reads_lst1.append(row['StartRead'])
            negtive_reads_lst2.append(row['EndRead'])
            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]

        # # isolates negative
        # for idx, row in negative_df.iterrows():
        #     read = row['StartRead']
        #     pos_reads = enumerate_ed1_seqs(read)
        #     for read2 in pos_reads:
        #         negtive_reads_lst1.append(read) 
        #         negtive_reads_lst2.append(read2)  
        #         cur_err_tye_kmers = error_type_classification(read, read2)
        #         cur_err_tye = cur_err_tye_kmers[0]
        #         cur_kmer1 = cur_err_tye_kmers[1]
        #         cur_kmer2 = cur_err_tye_kmers[2]
        #         cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
        #         cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
        #         cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
        #         negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])  

        for idx, row in ambiguous_df.iterrows():
            ambiguous_reads_lst1.append(row['StartRead'])
            ambiguous_reads_lst2.append(row['EndRead'])
            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            ambiguous_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])#, row["StartDegree"]

        genuine_fea = self.read2vec(genuine_reads_lst1, genuine_reads_lst2, genuine_reads_features)
        negative_fea = self.read2vec(negtive_reads_lst1, negtive_reads_lst2, negtive_reads_features)
        ambiguous_fea = self.read2vec(ambiguous_reads_lst1, ambiguous_reads_lst2, ambiguous_reads_features)
        del genuine_reads_lst1, genuine_reads_lst2, genuine_reads_features, negtive_reads_lst1, negtive_reads_lst2, negtive_reads_features, ambiguous_reads_lst1, ambiguous_reads_lst2, ambiguous_reads_features

        read_features = genuine_fea + negative_fea
        labels = np.array([1] * len(genuine_fea) + [0] * len(negative_fea))
        train_data = np.array(read_features)
        shape1 = (len(labels), len(read_features[0]))
        train_data.reshape(shape1)   

        ambiguous_data = np.array(ambiguous_fea)
        shape2 = (len(ambiguous_fea), len(ambiguous_fea[0]))
        ambiguous_data.reshape(shape2) 
        # scaling data
        self.logger.debug(train_data.shape)
        self.logger.debug(ambiguous_data.shape)
        train, ambiguous = self.scaler(train_data, ambiguous_data, high_flag=True)
        self.logger.debug(train[0])
        del train_data, ambiguous_data, genuine_fea, negative_fea, ambiguous_fea, read_features
        return train, labels, ambiguous

    def all_in_one_embedding(self, total_reads, genuine_df, negative_df, ambiguous_df, high_flag):
        self.logger.info("  Embedding genuine and ambiguous data.")
        genuine_reads_lst1 = []
        negtive_reads_lst1 = []
        ambiguous_reads_lst1 = []

        genuine_reads_lst2 = []
        negtive_reads_lst2 = []
        ambiguous_reads_lst2 = []

        genuine_reads_features = []
        negtive_reads_features = []
        ambiguous_reads_features = [] 

        if self.edit_dis == 1:
            base_lst = ['A', 'C', 'G', 'T', 'N']
            if self.config.read_type == "DNA":
                error_tyes = ["A-G", "G-A", "A-T", "T-A", "A-C", "C-A", "G-T", "T-G", "G-C", "C-G", "T-C", "C-T", "T-X", "X-T", "C-X", "X-C", "A-X", "X-A", "G-X", "X-G", "X-N", "N-X"]
                
                kmers = ['NX', 'XN', 'XA', 'XC', 'XG', 'XT', 'AX', 'CX', 'GX', 'TX'] + [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)] + ['AXA', 'CXA', 'TXA', 'GXA', 'AXC', 'CXC', 'TXC', 'GXC', 'AXG', 'CXG', 'TXG', 'GXG', 'AXT', 'CXT', 'TXT', 'GXT', 'NXA', 'NXC', 'NXG', 'NXT', 'NXN','AXN','CXN','GXN','TXN']
                # kmers3 = ['ANA', 'CNA', 'TNA', 'GNA', 'ANC', 'CNC', 'TNC', 'GNC', 'ANG', 'CNG', 'TNG', 'GNG', 'ANT', 'CNT', 'TNT', 'GNT'] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]
                # kmers = kmers2 + kmers3
            elif self.config.read_type == "RNA":
                # error_tyes = ["A-G", "G-A", "A-U", "U-A", "A-C", "C-A", "G-U", "U-G", "G-C", "C-G", "U-C", "C-U", "U-N", "N-U", "C-N", "N-C", "A-N", "N-A", "G-N", "N-G"]
                error_tyes = ["A-G", "G-A", "A-T", "T-A", "A-C", "C-A", "G-T", "T-G", "G-C", "C-G", "T-C", "C-T", "T-X", "X-T", "C-X", "X-C", "A-X", "X-A", "G-X", "X-G", "X-N", "N-X"]
                kmers = ['NX', 'XN', 'XA', 'XC', 'XG', 'XT', 'AX', 'CX', 'GX', 'TX'] + [''.join(i) for i in itertools.product(base_lst, repeat = 2)] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)] + ['AXA', 'CXA', 'TXA', 'GXA', 'AXC', 'CXC', 'TXC', 'GXC', 'AXG', 'CXG', 'TXG', 'GXG', 'AXT', 'CXT', 'TXT', 'GXT', 'NXC', 'NXG', 'NXU', 'NXN','AXN','CXN','GXN','UXN']

            error_tye_priors = {}
            err_tye_base_prior = 0.1
            for it in error_tyes:
                error_tye_priors[it] = err_tye_base_prior
                err_tye_base_prior += 0.01

            # base_lst = ['A', 'C', 'G', 'T']
            # kmers2 = ['AN', 'CN', 'TN', 'GN', 'NA', 'NC', 'NT', 'NG'] + [''.join(i) for i in itertools.product(base_lst, repeat = 2)]
            # kmers3 = ['ANA', 'CNA', 'TNA', 'GNA', 'ANC', 'CNC', 'TNC', 'GNC', 'ANG', 'CNG', 'TNG', 'GNG', 'ANT', 'CNT', 'TNT', 'GNT'] + [''.join(i) for i in itertools.product(base_lst, repeat = 3)]
            # kmers = kmers2 + kmers3

            self.logger.debug(f'expected kmers: {kmers}, {len(kmers)}')
            total_err_tyes = []
            total_err_kmers = []
            for idx, row in genuine_df.iterrows():
                total_err_tyes.append(row['ErrorTye'])
                total_err_kmers.append(row['StartErrKmer'])
                total_err_kmers.append(row['EndErrKmer'])

            for idx, row in ambiguous_df.iterrows():
                total_err_tyes.append(row['ErrorTye'])
                total_err_kmers.append(row['StartErrKmer'])
                total_err_kmers.append(row['EndErrKmer'])

            total_err_kmers_count = len(total_err_kmers)
            total_err_tyes_count = len(total_err_tyes)
            
            err_kmers2count = Counter(total_err_kmers)
            err_tyes2count = Counter(total_err_tyes)

            kmer_keys = err_kmers2count.keys()
            self.logger.debug(f'Real kmers: {kmer_keys}, {len(kmer_keys)}')
    
            not_exist_kmer = set(kmers) - set(kmer_keys)
            if not_exist_kmer:
                for mer in not_exist_kmer:
                    err_kmers2count[mer] = 0

            kmers_priors = {}
            base_prior = 0.1
            for item in kmer_keys:
                kmers_priors[item] = base_prior
                base_prior += 0.01        
            self.logger.debug(err_kmers2count)
            self.logger.debug(err_tyes2count)
            juge_indels = err_tyes2count.keys() & set(['N-X', 'X-N', 'X-A', 'X-C', 'X-G', 'X-T', 'A-X', 'C-X', 'G-X', 'T-X'])
            indel_num = len(juge_indels)

        for idx, row in genuine_df.iterrows():
            genuine_reads_lst1.append(row['StartRead'])
            genuine_reads_lst2.append(row['EndRead'])
            if self.edit_dis == 1:
                cur_err_tye = row['ErrorTye']
                cur_kmer1 = row['StartErrKmer']
                cur_kmer2 = row['EndErrKmer']
                # self.logger.debug(row['StartRead'])
                # self.logger.debug(row['EndRead'])
                # self.logger.debug(f'{cur_kmer1}, {cur_kmer2}')
                cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                # 
                if high_flag:
                    genuine_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
                else:
                    genuine_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
            else:
                genuine_reads_features.append([row['StartReadCount']])
            # genuine_reads_count.append(row['StartReadCount'])
        
        negative_count = 0
        neg_read2seqs = {}
        if self.edit_dis == 1:
            for idx, row in negative_df.iterrows():
                read = row['StartRead']
                if indel_num == 0:
                    pos_reads = seq2substitution(read)
                else:
                    pos_reads = enumerate_ed1_seqs(read)  
                neg_read2seqs[read] = pos_reads
                negative_count += len(pos_reads)          

        # if self.edit_dis == 2 or not self.config.over_sampling:
        # if self.edit_dis == 2 or negative_count >= 500000:
        if self.edit_dis == 2 or negative_count >= self.config.negative_sample_num:
            genuine_reads_num = len(genuine_df)
            negative_num = len(negative_df)
            if genuine_reads_num >= negative_num * 2:
                each_negative_num = int(genuine_reads_num / negative_num)
            else:
                each_negative_num = 1
        if negative_count >= self.config.negative_sample_num:
            self.logger.warning("Negative samples larger than {}, noise2read will use random sampling to generate negative samples.".format(self.config.negative_sample_num))

        for idx, row in negative_df.iterrows():
            read = row['StartRead']
            if self.edit_dis == 1:
                if negative_count < self.config.negative_sample_num:
                    pos_reads = neg_read2seqs[read]
                    for read2 in pos_reads:
                        negtive_reads_lst1.append(read) 
                        negtive_reads_lst2.append(read2)  
                        cur_err_tye_kmers = error_type_classification(read, read2)
                        cur_err_tye = cur_err_tye_kmers[0]
                        cur_kmer1 = cur_err_tye_kmers[1]
                        cur_kmer2 = cur_err_tye_kmers[2]
                        cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                        cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                        cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                        if high_flag:
                            negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                        else:
                            negtive_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
                else:
                    pos_reads = neg_read2seqs[read]
                    candidates = list(pos_reads - total_reads)
                    if len(candidates) > 0:
                        for i in range(each_negative_num):
                            if len(candidates) > 0:
                                negtive_reads_lst1.append(read) 
                                
                                select_read = random.choice(candidates) 
                                candidates.remove(select_read)
                                negtive_reads_lst2.append(select_read)  

                                cur_err_tye_kmers = error_type_classification(read, select_read)
                                cur_err_tye = cur_err_tye_kmers[0]
                                cur_kmer1 = cur_err_tye_kmers[1]
                                cur_kmer2 = cur_err_tye_kmers[2]
                                cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                                cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                                cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                                # 
                                if high_flag:
                                    negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                                else:
                                    negtive_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
            elif self.edit_dis == 2:
                # pos_reads = enumerate_ed2_seqs(read) 
                random_seqs = random_ed2_seq(read, total_reads, each_negative_num)
                for i in range(len(random_seqs)):
                    negtive_reads_lst1.append(read) 
                    negtive_reads_lst2.append(random_seqs[i])  
                    negtive_reads_features.append([row['StartReadCount']]) #, row["StartDegree"]

        '''
        for idx, row in negative_df.iterrows():
            read = row['StartRead']
            if self.edit_dis == 1:
                if self.config.over_sampling:
                    if indel_num == 0:
                        pos_reads = seq2substitution(read)
                    else:
                        pos_reads = enumerate_ed1_seqs(read)
                    # candidates = list(pos_reads - total_reads)
                    # for i in range(each_negative_num):
                    for read2 in pos_reads:
                        negtive_reads_lst1.append(read) 
                        # select_read = random.choice(candidates) 
                        # candidates.remove(select_read)
                        negtive_reads_lst2.append(read2)  

                        cur_err_tye_kmers = error_type_classification(read, read2)
                        cur_err_tye = cur_err_tye_kmers[0]
                        cur_kmer1 = cur_err_tye_kmers[1]
                        cur_kmer2 = cur_err_tye_kmers[2]
                        cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                        cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                        cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                        # 
                        if high_flag:
                            negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                        else:
                            negtive_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
                else:
                    if indel_num == 0:
                        pos_reads = seq2substitution(read)
                    else:
                        pos_reads = enumerate_ed1_seqs(read)
                    candidates = list(pos_reads - total_reads)
                    for i in range(each_negative_num):
                        negtive_reads_lst1.append(read) 
                        select_read = random.choice(candidates) 
                        candidates.remove(select_read)
                        negtive_reads_lst2.append(select_read)  

                        cur_err_tye_kmers = error_type_classification(read, select_read)
                        cur_err_tye = cur_err_tye_kmers[0]
                        cur_kmer1 = cur_err_tye_kmers[1]
                        cur_kmer2 = cur_err_tye_kmers[2]
                        cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                        cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                        cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                        # 
                        if high_flag:
                            negtive_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                        else:
                            negtive_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
            elif self.edit_dis == 2:
                # pos_reads = enumerate_ed2_seqs(read) 
                random_seqs = random_ed2_seq(read, total_reads, each_negative_num)
                for i in range(len(random_seqs)):
                    negtive_reads_lst1.append(read) 
                    negtive_reads_lst2.append(random_seqs[i])  
                    negtive_reads_features.append([row['StartReadCount']]) #, row["StartDegree"]
                    # negtive_reads_count.append(row['StartReadCount']) 
        '''
        for idx, row in ambiguous_df.iterrows():
            ambiguous_reads_lst1.append(row['StartRead'])
            ambiguous_reads_lst2.append(row['EndRead'])
            if self.edit_dis == 1:
                cur_err_tye = row['ErrorTye']
                cur_kmer1 = row['StartErrKmer']
                cur_kmer2 = row['EndErrKmer']
                cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                # 
                if high_flag:
                    ambiguous_reads_features.append([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])#, row["StartDegree"]
                else:
                    ambiguous_reads_features.append([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])#, row["StartDegree"]
            else:
                ambiguous_reads_features.append([row['StartReadCount']])

        # self.logger.debug(f'{len(negtive_reads_lst1), len(negtive_reads_lst2), len(negtive_reads_features)}')
        genuine_fea = self.read2vec(genuine_reads_lst1, genuine_reads_lst2, genuine_reads_features)
        negative_fea = self.read2vec(negtive_reads_lst1, negtive_reads_lst2, negtive_reads_features)
        ambiguous_fea = self.read2vec(ambiguous_reads_lst1, ambiguous_reads_lst2, ambiguous_reads_features)
        del genuine_reads_lst1, genuine_reads_lst2, genuine_reads_features, negtive_reads_lst1, negtive_reads_lst2, negtive_reads_features, ambiguous_reads_lst1, ambiguous_reads_lst2, ambiguous_reads_features

        read_features = genuine_fea + negative_fea
        labels = np.array([1] * len(genuine_fea) + [0] * len(negative_fea))
        train_data = np.array(read_features)
        shape1 = (len(labels), len(read_features[0]))
        train_data.reshape(shape1)   

        ambiguous_data = np.array(ambiguous_fea)
        shape2 = (len(ambiguous_fea), len(ambiguous_fea[0]))
        ambiguous_data.reshape(shape2) 
        # scaling data
        self.logger.debug(train_data.shape)
        self.logger.debug(ambiguous_data.shape)
        if self.edit_dis == 1:
            train, ambiguous = self.scaler(train_data, ambiguous_data, high_flag)
        else:
            train, ambiguous = self.scaler2(train_data, ambiguous_data)
        self.logger.debug(train[0])
        del train_data, ambiguous_data, genuine_fea, negative_fea, ambiguous_fea, read_features
        return train, labels, ambiguous

    def scaler2(self, lab_fea, ambiguous_ulab_fea):
        # pos0 = 179
        pos0 = 65
        lab_fea0 = lab_fea[:,0:pos0]

        ulab_fea0 = ambiguous_ulab_fea[:,0:pos0]
        scaler0 = StandardScaler()
        lab_scale_f0 = scaler0.fit_transform(lab_fea0)
        ulab_scale_f0 = scaler0.transform(ulab_fea0)

        # read integer encoding
        pos = (self.config.read_max_len+1) * 2
        # pos = (self.config.read_max_len+1) * 4  * 2
        lab_fea1 = lab_fea[:,pos0:pos+pos0]

        ulab_fea1 = ambiguous_ulab_fea[:,pos0:pos+pos0]
        scaler1 = StandardScaler()
        lab_scale_f1 = scaler1.fit_transform(lab_fea1)
        ulab_scale_f1 = scaler1.transform(ulab_fea1)

        # read count
        lab_fea2 = lab_fea[:,-4].reshape(-1, 1)
        ulab_fea2 = ambiguous_ulab_fea[:,-4].reshape(-1, 1)
        scaler2 = StandardScaler()
        lab_scale_f2 = scaler2.fit_transform(lab_fea2)
        ulab_scale_f2 = scaler2.transform(ulab_fea2)

        lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2))
        ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2))
        del lab_scale_f0, lab_scale_f1, lab_scale_f2, ulab_scale_f0, ulab_scale_f1, ulab_scale_f2
        return lab_scale_f, ulab_scale_f 

    def scaler(self, lab_fea, ambiguous_ulab_fea, high_flag):
        # pos0 = 179
        pos0 = 65
        lab_fea0 = lab_fea[:,0:pos0]

        ulab_fea0 = ambiguous_ulab_fea[:,0:pos0]
        scaler0 = StandardScaler()
        lab_scale_f0 = scaler0.fit_transform(lab_fea0)
        ulab_scale_f0 = scaler0.transform(ulab_fea0)

        # read integer encoding
        pos = (self.config.read_max_len+1) * 2
        # pos = self.config.read_max_len+1
        # pos = (self.config.read_max_len+1) * 4 * 2
        lab_fea1 = lab_fea[:,pos0:pos+pos0]

        ulab_fea1 = ambiguous_ulab_fea[:,pos0:pos+pos0]
        scaler1 = StandardScaler()
        lab_scale_f1 = scaler1.fit_transform(lab_fea1)
        ulab_scale_f1 = scaler1.transform(ulab_fea1)

        # error type and error kmers
        lab_fea3 = lab_fea[:,-3:]
        ulab_fea3 = ambiguous_ulab_fea[:,-3:]
        scaler3 = StandardScaler()
        lab_scale_f3 = scaler3.fit_transform(lab_fea3)
        ulab_scale_f3 = scaler3.transform(ulab_fea3)
        
        # read count
        if not high_flag:
            lab_fea2 = lab_fea[:,-4].reshape(-1, 1)
            ulab_fea2 = ambiguous_ulab_fea[:,-4].reshape(-1, 1)
            scaler2 = StandardScaler()
            lab_scale_f2 = scaler2.fit_transform(lab_fea2)
            ulab_scale_f2 = scaler2.transform(ulab_fea2)

            lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3))
            ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3))
            # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f3))
            # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f3))
            del lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3, ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3
        else:
            lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f3))
            ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f3))
            del lab_scale_f0, lab_scale_f1, lab_scale_f3, ulab_scale_f0, ulab_scale_f1, ulab_scale_f3
                
        # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3))
        # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3))
        
        return lab_scale_f, ulab_scale_f 

    def read2vec(self, reads_lst1, reads_lst2, other_features):
        read_features = []    
        ES = EncodeScheme(self.config.read_max_len, self.config.entropy_kmer, self.config.entropy_q, self.config.kmer_freq, self.config.read_type)

        shared_objects = ES, reads_lst1, reads_lst2, other_features
        # with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
        #     with tqdm(total=len(reads_lst1), desc=self.logger.info("Encoding reads")) as pbar:   
        #         for fea_lst in pool.imap(self.read2features, range(len(reads_lst1))):
        #             read_features.append(fea_lst)
        #             pbar.update()  
        with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
            for fea_lst in pool.imap(self.read2features, range(len(reads_lst1)), progress_bar=self.config.verbose):
                read_features.append(fea_lst)
        
        self.logger.debug(f'{len(read_features)}, {len(read_features[0])}')
        return read_features