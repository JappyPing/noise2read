# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-09-07 00:15:00

from typing import Counter
import numpy as np
# import os
from sklearn.preprocessing import StandardScaler
from noise2read.encoding import EncodeScheme
from mpire import WorkerPool
# from tqdm import tqdm
from noise2read.utils import *
import itertools
import pickle
from noise2read.utils import MemoryMonitor

class Reads2Vectors():
    def __init__(self, logger, config, edit_dis):
        self.logger = logger
        self.edit_dis = edit_dis
        self.config = config
        # Create an instance of the MemoryMonitor
        self.MM = MemoryMonitor(self.logger)

    def read2features(self, shared_objs, idx):
        ES, ori_features = shared_objs
        cur_feature = ori_features[idx]
        features = []
        ft_fea1 = ES.descriptors("FourierTransform", cur_feature[0])
        # cg_fea1 = ES.descriptors("ChaosGame", cur_feature[0])
        entropy_fea1 = ES.descriptors("Entropy", cur_feature[0])
        fs_fea1 = ES.descriptors("FickettScore", cur_feature[0])

        features.extend(ft_fea1)
        # features.extend(cg_fea1)
        features.extend(entropy_fea1)
        features.extend(fs_fea1)

        # atomic_fea1 = ES.descriptors("atomic_number", cur_feature[0])
        # atomic_fea2 = ES.descriptors("atomic_number", cur_feature[1])
        # features.extend(atomic_fea1)
        # features.extend(atomic_fea2)

        features.extend(cur_feature[2:])
        # self.logger.debug(f'FourierTransform: {len(ft_fea1)}, ChaosGame: {len(cg_fea1)}, Entropy: {len(entropy_fea1)}, FickettScore: {len(fs_fea1)}')
        return features 

    def all_in_one_embedding(self, total_reads, genuine_df, negative_df, ambiguous_df, high_flag):
        self.logger.info("Embedding genuine and ambiguous data.")
        self.MM.start()
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
            # for idx, row in genuine_df.iterrows():
            #     total_err_tyes.append(row['ErrorTye'])
            #     total_err_kmers.append(row['StartErrKmer'])
            #     total_err_kmers.append(row['EndErrKmer'])

            # for idx, row in ambiguous_df.iterrows():
            #     total_err_tyes.append(row['ErrorTye'])
            #     total_err_kmers.append(row['StartErrKmer'])
            #     total_err_kmers.append(row['EndErrKmer'])

            total_err_tyes.extend(genuine_df['ErrorTye'].tolist())
            total_err_tyes.extend(ambiguous_df['ErrorTye'].tolist())

            total_err_kmers.extend(genuine_df['StartErrKmer'].tolist())
            total_err_kmers.extend(ambiguous_df['StartErrKmer'].tolist())

            total_err_kmers.extend(genuine_df['EndErrKmer'].tolist())
            total_err_kmers.extend(ambiguous_df['EndErrKmer'].tolist())

            total_err_kmers_count = len(total_err_kmers)
            total_err_tyes_count = len(total_err_tyes)
            
            err_kmers2count = Counter(total_err_kmers)
            del total_err_kmers
            err_tyes2count = Counter(total_err_tyes)
            del total_err_tyes

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
        self.MM.measure()
        genuine_feature_lst = []
        for idx, row in genuine_df.iterrows():
            # genuine_feature_lst.append(row['StartRead'])
            # genuine_feature_lst.append(row['EndRead'])
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
                    # genuine_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
                    genuine_feature_lst.append((row['StartRead'], row['EndRead'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
                else:
                    # genuine_feature_lst.extend([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
                    genuine_feature_lst.append((row['StartRead'], row['EndRead'], row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
            else:
                # genuine_feature_lst.append(row['StartReadCount'])
                genuine_feature_lst.append((row['StartRead'], row['EndRead'], row['StartReadCount']))
            # genuine_reads_count.append(row['StartReadCount'])

        genuine_fea = self.read2vec(genuine_feature_lst)
        del genuine_feature_lst
        self.MM.measure()
        #################################################################################
        # print("encoding negative samples...")
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

        negative_feature_lst = []
        for idx, row in negative_df.iterrows():
            read = row['StartRead']
            if self.edit_dis == 1:
                if negative_count < self.config.negative_sample_num:
                    pos_reads = neg_read2seqs[read]
                    for read2 in pos_reads:
                        # negative_feature_lst.append(read) 
                        # negative_feature_lst.append(read2)  
                        cur_err_tye_kmers = error_type_classification(read, read2)
                        cur_err_tye = cur_err_tye_kmers[0]
                        cur_kmer1 = cur_err_tye_kmers[1]
                        cur_kmer2 = cur_err_tye_kmers[2]
                        cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                        cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                        cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                        if high_flag:
                            # negative_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                            negative_feature_lst.append((read, read2, cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
                        else:
                            # negative_feature_lst.extend([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
                            negative_feature_lst.append((read, read2, row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
                else:
                    pos_reads = neg_read2seqs[read]
                    candidates = list(pos_reads - total_reads)
                    if len(candidates) > 0:
                        for i in range(each_negative_num):
                            if len(candidates) > 0:
                                # negative_feature_lst.append(read) 
                                
                                select_read = random.choice(candidates) 
                                candidates.remove(select_read)
                                # negative_feature_lst.append(select_read)  

                                cur_err_tye_kmers = error_type_classification(read, select_read)
                                cur_err_tye = cur_err_tye_kmers[0]
                                cur_kmer1 = cur_err_tye_kmers[1]
                                cur_kmer2 = cur_err_tye_kmers[2]
                                cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                                cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                                cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                                # 
                                if high_flag:
                                    # negative_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
                                    negative_feature_lst.append((read, select_read, cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
                                else:
                                    # negative_feature_lst.extend([row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) 
                                    negative_feature_lst.append((read, select_read, row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
            elif self.edit_dis == 2:
                # pos_reads = enumerate_ed2_seqs(read) 
                random_seqs = random_ed2_seq(read, total_reads, each_negative_num)
                for i in range(len(random_seqs)):
                    # negative_feature_lst.append(read) 
                    # negative_feature_lst.append(random_seqs[i])  
                    # negative_feature_lst.append(row['StartReadCount']) #, row["StartDegree"]
                    negative_feature_lst.append((read, random_seqs[i], row['StartReadCount']))

        negative_fea = self.read2vec(negative_feature_lst)
        del negative_feature_lst
        self.MM.measure()
        ###############################################################
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
        ambiguous_feature_lst = []
        for idx, row in ambiguous_df.iterrows():
            if self.edit_dis == 1:
                cur_err_tye = row['ErrorTye']
                cur_kmer1 = row['StartErrKmer']
                cur_kmer2 = row['EndErrKmer']
                cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                # 
                if high_flag:
                    ambiguous_feature_lst.append((row['StartRead'], row['EndRead'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
                else:
                    ambiguous_feature_lst.append((row['StartRead'], row['EndRead'], row['StartReadCount'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
            else:
                ambiguous_feature_lst.append((row['StartRead'], row['EndRead'], row['StartReadCount']))

        ambiguous_fea = self.read2vec(ambiguous_feature_lst)
        del ambiguous_feature_lst
        self.MM.measure()
        ##################################################################
        read_features = genuine_fea + negative_fea
        feature_len = len(read_features[0])
        train_data = np.array(read_features, dtype=object)
        gen_len = len(genuine_fea)
        neg_len = len(negative_fea)
        del read_features, genuine_fea, negative_fea

        labels = np.array([1] * gen_len + [0] * neg_len)
        shape1 = ((gen_len + neg_len), feature_len)
        train_data.reshape(shape1)   

        ambiguous_data = np.array(ambiguous_fea, dtype=object)
        shape2 = (len(ambiguous_fea), len(ambiguous_fea[0]))
        del ambiguous_fea
        ambiguous_data.reshape(shape2) 
        # scaling data
        self.logger.debug(train_data.shape)
        self.logger.debug(ambiguous_data.shape)
        if self.edit_dis == 1:
            train, ambiguous = self.scaler(train_data, ambiguous_data, high_flag)
        else:
            train, ambiguous = self.scaler2(train_data, ambiguous_data)
        self.logger.debug(train[0])
        del train_data, ambiguous_data
        self.MM.measure()
        self.MM.stop()
        return train, labels, ambiguous

    def scaler2(self, lab_fea, ambiguous_ulab_fea):
        # pos0 = 179
        # pos0 = 65
        pos0 = 27
        lab_fea0 = lab_fea[:,0:pos0]

        ulab_fea0 = ambiguous_ulab_fea[:,0:pos0]
        scaler0 = StandardScaler()
        lab_scale_f0 = scaler0.fit_transform(lab_fea0)
        ulab_scale_f0 = scaler0.transform(ulab_fea0)

        # # read integer encoding
        # pos = (self.config.read_max_len+1) * 2
        # # pos = (self.config.read_max_len+1) * 4  * 2
        # lab_fea1 = lab_fea[:,pos0:pos+pos0]

        # ulab_fea1 = ambiguous_ulab_fea[:,pos0:pos+pos0]
        # scaler1 = StandardScaler()
        # lab_scale_f1 = scaler1.fit_transform(lab_fea1)
        # ulab_scale_f1 = scaler1.transform(ulab_fea1)

        # read count
        lab_fea2 = lab_fea[:,-4].reshape(-1, 1)
        ulab_fea2 = ambiguous_ulab_fea[:,-4].reshape(-1, 1)
        scaler2 = StandardScaler()
        lab_scale_f2 = scaler2.fit_transform(lab_fea2)
        ulab_scale_f2 = scaler2.transform(ulab_fea2)

        # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2))
        # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2))
        # del lab_scale_f0, lab_scale_f1, lab_scale_f2, ulab_scale_f0, ulab_scale_f1, ulab_scale_f2
        lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f2))
        ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f2))
        del lab_scale_f0, lab_scale_f2, ulab_scale_f0, ulab_scale_f2
        return lab_scale_f, ulab_scale_f 

    def scaler(self, lab_fea, ambiguous_ulab_fea, high_flag):
        # pos0 = 179
        # pos0 = 65
        pos0 = 27
        lab_fea0 = lab_fea[:,0:pos0]

        ulab_fea0 = ambiguous_ulab_fea[:,0:pos0]
        scaler0 = StandardScaler()
        lab_scale_f0 = scaler0.fit_transform(lab_fea0)
        ulab_scale_f0 = scaler0.transform(ulab_fea0)

        # read integer encoding
        # pos = (self.config.read_max_len+1) * 2
        # # pos = self.config.read_max_len+1
        # # pos = (self.config.read_max_len+1) * 4 * 2
        # lab_fea1 = lab_fea[:,pos0:pos+pos0]

        # ulab_fea1 = ambiguous_ulab_fea[:,pos0:pos+pos0]
        # scaler1 = StandardScaler()
        # lab_scale_f1 = scaler1.fit_transform(lab_fea1)
        # ulab_scale_f1 = scaler1.transform(ulab_fea1)

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

            # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3))
            # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3))
            # # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f3))
            # # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f3))
            # del lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3, ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3

            lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f2, lab_scale_f3))
            ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f2, ulab_scale_f3))
            # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f3))
            # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f3))
            del lab_scale_f0, lab_scale_f2, lab_scale_f3, ulab_scale_f0, ulab_scale_f2, ulab_scale_f3

        else:
            # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f3))
            # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f3))
            # del lab_scale_f0, lab_scale_f1, lab_scale_f3, ulab_scale_f0, ulab_scale_f1, ulab_scale_f3

            lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f3))
            ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f3))
            del lab_scale_f0, lab_scale_f3, ulab_scale_f0, ulab_scale_f3

        # lab_scale_f = np.hstack((lab_scale_f0, lab_scale_f1, lab_scale_f2, lab_scale_f3))
        # ulab_scale_f = np.hstack((ulab_scale_f0, ulab_scale_f1, ulab_scale_f2, ulab_scale_f3))
        
        return lab_scale_f, ulab_scale_f 

    '''
    def read2features(self, ES, ori_feature):
        # ES, reads_lst1, reads_lst2, other_features = shared_objects
        features = []
        
        ###########################################################################
        ft_fea1 = ES.descriptors("FourierTransform", ori_feature[0])
        cg_fea1 = ES.descriptors("ChaosGame", ori_feature[0])
        entropy_fea1 = ES.descriptors("Entropy", ori_feature[0])
        fs_fea1 = ES.descriptors("FickettScore", ori_feature[0])
        features.extend(ft_fea1)
        features.extend(cg_fea1)
        features.extend(entropy_fea1)
        features.extend(fs_fea1)
        atomic_fea1 = ES.descriptors("atomic_number",ori_feature[0])
        atomic_fea2 = ES.descriptors("atomic_number", ori_feature[1])
        features.extend(atomic_fea1)
        features.extend(atomic_fea2)
        
        features.extend(ori_feature[2:])
        return features 
    '''
    def read2vec(self, original_features_lst):  
        ES = EncodeScheme(self.config.read_max_len, self.config.entropy_kmer, self.config.entropy_q, self.config.kmer_freq, self.config.read_type)

        if len(original_features_lst) > (self.config.chunks_num * 10):
            chunk_size = len(original_features_lst) // self.config.chunks_num
            # remainder = len(original_features_lst) % self.config.chunks_num
            # chunks = [original_features_lst[i:i+chunk_size] for i in range(0, len(original_features_lst), chunk_size)]
            combined_data = []
            if chunk_size >= self.config.num_workers:
                chunks = [original_features_lst[i:i+chunk_size] for i in range(0, len(original_features_lst), chunk_size)]
                chunk_names = []
                # combined_data = []
                # for i in range(len(chunks)):
                i = 0
                while chunks:
                    chunk = chunks.pop(0)
                    vectors = []
                    try:
                        shared_objects = ES, chunk
                        with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
                            for item in pool.imap(self.read2features, range(len(chunk))):
                                vectors.append(item)
                        del shared_objects
                    except KeyboardInterrupt:
                        # Handle termination signal (Ctrl+C)
                        pool.terminate()  # Terminate the WorkerPool before exiting
                    except Exception:
                        # Handle other exceptions
                        pool.terminate()  # Terminate the WorkerPool before exiting
                        raise
                    del chunk        
                    # Generate the pickle file name
                    file_name = self.config.result_dir + f"chunk_{i}.pickle"
                    i += 1
                    # Write the vectors to the pickle file
                    with open(file_name, "wb") as file:
                        pickle.dump(vectors, file)
                    chunk_names.append(file_name)
                    del vectors
                del chunks
                
                for file_name in chunk_names:
                    with open(file_name, "rb") as file:
                        vectors = pickle.load(file)
                        combined_data.extend(vectors)
                        del vectors
                    os.remove(file_name)
            else:
                try:
                    shared_objects = ES, original_features_lst
                    with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
                        for item in pool.imap(self.read2features, range(len(original_features_lst))):
                            combined_data.append(item)
                    del shared_objects
                except KeyboardInterrupt:
                    # Handle termination signal (Ctrl+C)
                    pool.terminate()  # Terminate the WorkerPool before exiting
                except Exception:
                    # Handle other exceptions
                    pool.terminate()  # Terminate the WorkerPool before exiting
                    raise   
            del ES, original_features_lst         
            return combined_data
        else:
            try:
                vectors = []
                shared_objects = ES, original_features_lst
                with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
                    for item in pool.imap(self.read2features, range(len(original_features_lst))):
                        vectors.append(item)
                del shared_objects
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
            del ES, original_features_lst    
            return vectors

    
    def high_all_in_one_embedding(self, genuine_df, negative_df, new_negative_df, ambiguous_df):
        self.logger.info("Embedding genuine and high ambiguous data.")
        self.MM.start()
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
        total_err_tyes.extend(genuine_df['ErrorTye'].tolist())
        total_err_tyes.extend(ambiguous_df['ErrorTye'].tolist())
        total_err_tyes.extend(new_negative_df['ErrorTye'].tolist())

        total_err_kmers.extend(genuine_df['StartErrKmer'].tolist())
        total_err_kmers.extend(ambiguous_df['StartErrKmer'].tolist())
        total_err_kmers.extend(new_negative_df['StartErrKmer'].tolist())

        total_err_kmers.extend(genuine_df['EndErrKmer'].tolist())
        total_err_kmers.extend(ambiguous_df['EndErrKmer'].tolist())
        total_err_kmers.extend(new_negative_df['EndErrKmer'].tolist())
        self.MM.measure()
        # for idx, row in genuine_df.iterrows():
        #     total_err_tyes.append(row['ErrorTye'])
        #     total_err_kmers.append(row['StartErrKmer'])
        #     total_err_kmers.append(row['EndErrKmer'])

        # for idx, row in ambiguous_df.iterrows():
        #     total_err_tyes.append(row['ErrorTye'])
        #     total_err_kmers.append(row['StartErrKmer'])
        #     total_err_kmers.append(row['EndErrKmer'])

        # for idx, row in new_negative_df.iterrows():
        #     # read = row['StartRead']
        #     total_err_tyes.append(row['ErrorTye'])
        #     total_err_kmers.append(row['StartErrKmer'])
        #     total_err_kmers.append(row['EndErrKmer'])
            
        total_err_kmers_count = len(total_err_kmers)
        total_err_tyes_count = len(total_err_tyes)
        
        err_kmers2count = Counter(total_err_kmers)
        del total_err_kmers
        err_tyes2count = Counter(total_err_tyes)
        del total_err_tyes
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
        self.MM.measure()
        ##################################################################################
        genuine_feature_lst = []
        for idx, row in genuine_df.iterrows():
            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            # self.logger.debug(row['StartRead'])
            # self.logger.debug(row['EndRead'])
            # self.logger.debug(f'{cur_kmer1}, {cur_kmer2}')
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            # genuine_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, error_tyes[row['ErrorTye']] , row["StartDegree"], row['ErrorPosition']
            genuine_feature_lst.append((row['StartRead'], row['EndRead'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
        genuine_fea = self.read2vec(genuine_feature_lst)
        del genuine_feature_lst
        self.MM.measure()
        #################################################################################
        negative_feature_lst = []
        for idx, row in new_negative_df.iterrows():
            # negative_feature_lst.append()
            # negative_feature_lst.append()
            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            # negative_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2]) #, row["StartDegree"]
            negative_feature_lst.append((row['StartRead'], row['EndRead'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
        # isolates negative
        if self.config.iso_neg_high:
            for idx, row in negative_df.iterrows():
                read = row['StartRead']
                pos_reads = enumerate_ed1_seqs(read)
                for read2 in pos_reads:
                    # negative_feature_lst.append(read) 
                    # negative_feature_lst.append(read2)  
                    cur_err_tye_kmers = error_type_classification(read, read2)
                    cur_err_tye = cur_err_tye_kmers[0]
                    cur_kmer1 = cur_err_tye_kmers[1]
                    cur_kmer2 = cur_err_tye_kmers[2]
                    cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
                    cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
                    cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
                    # negative_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])  
                    negative_feature_lst.append((read, read2, cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
        
        negative_fea = self.read2vec(negative_feature_lst)
        del negative_feature_lst
        self.MM.measure()
        ###############################################################
        ambiguous_feature_lst = []
        for idx, row in ambiguous_df.iterrows():
            # ambiguous_feature_lst.append(row['StartRead'])
            # ambiguous_feature_lst.append(row['EndRead'])
            cur_err_tye = row['ErrorTye']
            cur_kmer1 = row['StartErrKmer']
            cur_kmer2 = row['EndErrKmer']
            cur_err_tye_val = (err_tyes2count[cur_err_tye] + error_tye_priors[cur_err_tye]) / (total_err_tyes_count + 1)
            cur_err_kmer_val1 = (err_kmers2count[cur_kmer1] + kmers_priors[cur_kmer1]) / (total_err_kmers_count + 1)
            cur_err_kmer_val2 = (err_kmers2count[cur_kmer2] + kmers_priors[cur_kmer2]) / (total_err_kmers_count + 1)
            # ambiguous_feature_lst.extend([cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2])#, row["StartDegree"]
            ambiguous_feature_lst.append((row['StartRead'], row['EndRead'], cur_err_tye_val, cur_err_kmer_val1, cur_err_kmer_val2))
        ambiguous_fea = self.read2vec(ambiguous_feature_lst)
        del ambiguous_feature_lst
        self.MM.measure()
        ##################################################################
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
        self.MM.measure()
        self.MM.stop()
        return train, labels, ambiguous

    '''
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
    '''

    '''
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
    '''

   