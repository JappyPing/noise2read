# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:04:45
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-09-09 09:58:23

import collections
from Bio import SeqIO
# import gzip
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
import copy
import random
from mpire import WorkerPool
import os
from noise2read.utils import *
from noise2read.data_generation import DataGneration
from noise2read.error_orrection import ErrorCorrection
from noise2read.data_preprocessing import DataProcessing
import numpy as np
import sys

class Simulation():
    # def __init__(self, logger, config, num_workers, f_in, error_rate, output_dir, substations, indels, seed=0):
    def __init__(self, logger, config):
        self.logger = logger
        self.config = config
        if os.path.exists(config.result_dir):
            self.logger.info("Directory '% s' already exists" % config.result_dir)
            # for f in os.listdir(self.config.result_dir):
            #     os.remove(os.path.join(self.config.result_dir, f))
        else:
            os.makedirs(config.result_dir)
            self.logger.info("Directory '% s' created" % config.result_dir)

    def enumerate_ed1_seqs(self, read):
        possible_ed1 = []
        if self.config.substations:
            possible_ed1.extend(self.seq2substitution(read))
        # if self.config.indels:
        #     possible_ed1.extend(self.seq2deletion(read))
        #     possible_ed1.extend(self.seq2insertion(read))
        return list(set(possible_ed1))

    def replace_char(self, seq, char, index):
        seq[index] = char
        return ''.join(seq)

    def sub_base(self, base):
        if base == 'A':
            return 'TCG'
        elif base == 'T':
            return 'ACG'
        elif base == 'C':
            return 'ATG'
        elif base == 'G':
            return 'ACT'
        elif base == 'N':
            return 'N'
            
    def seq2substitution(self, read):
        """
        enumerate all the substitutions for one read 
        Args:
            read (str): a sequence
        Returns:
            set: a set contains reads
        """
        editdis1_list = []
        # substitution
        raw_seq = list(read)
        n = len(raw_seq)
        for i in range(n):
            seq = copy.deepcopy(raw_seq)
            temp_base = seq[i]
            #print(temp_base)
            temp_sub_base = list(self.sub_base(temp_base))
            #print(temp_sub_base)
            for b in temp_sub_base:
                if b != 'N':
                    sub_seq = self.replace_char(seq, b, i)
                    editdis1_list.append(sub_seq)
        return set(editdis1_list)

    def seq2deletion(self, read):
        """
        enumerate all the deletions for one read 
        Args:
            read (str): a sequence
        Returns:
            set: a set contains reads
        """
        editdis1_list = []
        seq = list(read)
        n = len(seq)
        # deletion
        for i in range(n):
            del_seq = read[:i] + read[(i+1):]
            editdis1_list.append(del_seq)
        return set(editdis1_list)

    def seq2insertion(self, read):
        """
        enumerate all the insertions for one read 
        Args:
            read (str): a sequence
        Returns:
            set: a set contains reads
        """
        editdis1_list = []
        seq = list(read)
        n = len(seq)
        # insertion
        bases = ['A', 'G', 'C', 'T']
        for i in range(n+1):
            for b in bases:
                raw_seq = copy.deepcopy(seq)
                raw_seq.insert(i, b)
                editdis1_list.append(''.join(raw_seq))
        return set(editdis1_list)
        
    def read2err_read(self, read):
        possible_seqs = self.enumerate_ed1_seqs(read)
        num = len(possible_seqs)
        random.seed(self.config.sim_seed)
        random.shuffle(possible_seqs)
        idx = random.randint(0, num-1)
        return possible_seqs[idx]

    def error_correction(self, config):
        DG = DataGneration(self.logger, config)
        if config.high_ambiguous:
            isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df, high_ambiguous_df = DG.data_files(edit_dis=1)
        else:
            isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df = DG.data_files(edit_dis=1)      
        config.read_max_len = read_max_len
        ###############################################################
        EC = ErrorCorrection(self.logger, config)
        ## one model to predict
        if config.high_ambiguous:
            ##################################
            corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df)
        else:
            corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df =None) 
        if read_min_len > 30:
            genuine_df2, negative_df2, ambiguous_df2, unique_seqs2 = DG.extract_ed2_errors(corrected_file)
            correct_data = EC.all_in_one_ed2_correction(corrected_file, unique_seqs2, genuine_df2, negative_df2, ambiguous_df2)
        else:
            correct_data = corrected_file
        return correct_data        

    def simplify_error_correction(self, config):
        DG = DataGneration(self.logger, config)
        isolates_file, non_isolates_file, read_max_len, read_min_len, genuine_df, ambiguous_df = DG.simplify_data_files(config.input_file, edit_dis=1)      
        config.read_max_len = read_max_len
        ###############################################################
        EC = ErrorCorrection(self.logger, config)
        corrected_file = EC.simplify_correction(isolates_file, non_isolates_file, genuine_df, ambiguous_df)
            
        if read_min_len > config.min_read_len:
            genuine_df, ambiguous_df = DG.simplify_data_files(corrected_file, edit_dis=2) 
            correct_data = EC.simplify_2nt_correction(corrected_file, genuine_df, ambiguous_df)
        else:
            correct_data = corrected_file
        return correct_data    

    def simulation(self):
        # step 1 correct errors using noise2read
        # 
        # if os.path.exists(self.config.correct_data):
        #     corrected_data = self.config.correct_data
        # else:
            # corrected_data = self.error_correction(self.config)
        corrected_data = self.simplify_error_correction(self.config)
        # inject errors
        raw_data, true_data = self.error_injection(corrected_data)
        # raw_data, true_data = self.error_injection(self.config.input_file)
        DP = DataProcessing( 
            self.logger,
            self.config)
        umi_raw_dataset, umi_true_dataset = DP.write_mimic_umis(raw_data, true_data)
        if os.path.exists(raw_data):
            os.system("rm %s" % raw_data)
        if os.path.exists(true_data):
            os.system("rm %s" % true_data)
        if os.path.exists(corrected_data):
            os.system("rm %s" % corrected_data)
        return umi_raw_dataset, umi_true_dataset
    
    def error_injection(self, data_set):
        records_dict, file_type = parse_data_dict(data_set)
        id_lst = list(records_dict.keys())

        total_seqs_lst = []
        seq2id_dict = {}
        for item in id_lst:
            seq = str(records_dict[item].seq)
            seq2id_dict.setdefault(seq, []).append(item)
            total_seqs_lst.append(seq)
  
        read2counts = collections.Counter(total_seqs_lst)
        select_id_lst = []
        mutated_id_pool = []
        for read, count in read2counts.items():
            if count > self.config.min_freq:
                select_id_lst.extend(seq2id_dict[read])
            if count > self.config.min_read_count:
                mutated_id_pool.extend(seq2id_dict[read])

        if self.config.error_rate1 > 0 and self.config.error_rate1 < 1:
            flag1 = True
            err_1_id_lst = []
            total_id_num = len(mutated_id_pool)
            err_1base_read_num = round(total_id_num * self.config.error_rate1)

            random.shuffle(mutated_id_pool)
            err_1_id_lst.extend(random.sample(mutated_id_pool, err_1base_read_num))
      
            self.logger.info(f"Total candidature reads: {total_id_num}, 1-base-error per read error rate: {self.config.error_rate1}, 1-base-error reads number: {err_1base_read_num}")
        else:
            self.logger.error("Input invalid 1-base-error per read error rate!")
            sys.exit(1)

        # # two base errors
        if self.config.error_rate2 > 0 and self.config.error_rate2 < 1:
            flag2 = True
            remain_id_lst = list(set(mutated_id_pool) - set(err_1_id_lst))
            remain_id_num = len(remain_id_lst)
            err_2bases_read_num = round(remain_id_num * self.config.error_rate2)

            err_2_id_lst = []
            random.shuffle(remain_id_lst)
            err_2_id_lst.extend(random.sample(remain_id_lst, err_2bases_read_num))

            self.logger.info(f"Remaining candidature reads: {remain_id_num}, 2-base-errors per read error rate: {self.config.error_rate2}, 2-base-errors reads number: {err_2bases_read_num}")
        else:
            self.logger.error("Input invalid 2-base-errors per read error rate!")
            sys.exit(1)

        err_records = []
        shared_objects = records_dict, file_type
        if flag1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
                    for tmp_rec in pool.imap(self.mutate_seq, err_1_id_lst, progress_bar=self.config.verbose):
                        err_records.append(tmp_rec)
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting

            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        if flag2:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
                    for tmp_rec in pool.imap(self.mutate_2base_seq, err_2_id_lst, progress_bar=self.config.verbose):
                        err_records.append(tmp_rec)
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting

            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        del shared_objects
        err_subdataset = self.config.result_dir + "err_subdataset" + file_type
        with open(err_subdataset, "w") as handle:
            SeqIO.write(err_records, handle, file_type)

        ########################################################################################
        err_free_subdataset = self.config.result_dir + "err_free_subdataset" + file_type

        if flag1 and flag2:
            err_id_lst = err_1_id_lst + err_2_id_lst
        elif flag1 and not flag2:
            err_id_lst = err_1_id_lst
        elif not flag1 and flag2:
            err_id_lst = err_2_id_lst
        else:
            self.logger.error("No error id selected, simulation failed!")
            sys.exit(1)            
        err_free_id_lst = list(set(select_id_lst) - set(err_id_lst))
        extract_records(self.config.result_dir, err_free_id_lst, data_set, err_free_subdataset)

        base = self.config.input_file.split("/")[-1]
        file_name = base.split(file_type)[0]
        raw_data = self.config.result_dir + "simulated_raw_" + file_name + file_type
        os.system("cat %s %s > %s" % (err_subdataset, err_free_subdataset, raw_data))
        ##################################################################################
        err_free_subdataset_2 = self.config.result_dir + "err_free_subdataset2" + file_type
        extract_records(self.config.result_dir, err_id_lst, data_set, err_free_subdataset_2)
        true_data = self.config.result_dir + "simulated_true_" + file_name + file_type  
        os.system("cat %s %s > %s" % (err_free_subdataset_2, err_free_subdataset, true_data))
        if os.path.exists(err_free_subdataset_2):
            os.system("rm %s" % err_free_subdataset_2)   
        if os.path.exists(err_free_subdataset):
            os.system("rm %s" % err_free_subdataset)  
        if os.path.exists(err_subdataset):
            os.system("rm %s" % err_subdataset)          
        return raw_data, true_data
    
    def find_first_diff_char(self, string1, string2):
        min_len = min(len(string1), len(string2))

        for i in range(min_len):
            if string1[i] != string2[i]:
                return i

        # If all characters in the shorter string match the corresponding characters in the longer string
        # and the shorter string is a prefix of the longer string, return the length of the shorter string
        if len(string1) != len(string2):
            return min_len

        # If both strings are identical
        return -1

    def mutate_seq(self, shared_objects, idx):
        records_dict, file_type = shared_objects
        ori_seq = str(records_dict[idx].seq)

        possible_seqs = []
        possible_seqs.extend(self.seq2substitution(ori_seq))

        num = len(possible_seqs)
        # random.seed(idx)
        # Do not set seed here because it can produce the same position errors which is not random, or you can set different seed for different sequence
        random.shuffle(possible_seqs)
        ii = random.randint(0, num-1)
        mutation_seq = possible_seqs[ii]

        if file_type == "fastq":
            return SeqRecord(Seq(mutation_seq), id=records_dict[idx].id, description=records_dict[idx].description, letter_annotations=records_dict[idx].letter_annotations)
        elif file_type == "fasta":
            return SeqRecord(Seq(mutation_seq).seq, id=records_dict[idx].id, description=records_dict[idx].description)     

    def mutate_2base_seq(self, shared_objects, idx):
        records_dict, file_type = shared_objects
        ori_seq = str(records_dict[idx].seq)

        possible_seqs = self.enumerate_ed2_seqs(ori_seq)
        num = len(possible_seqs)
        # random.seed(idx)
        # Do not set seed here because it can produce the same position errors which is not random, or you can set different seed for different sequence
        random.shuffle(possible_seqs)
        ii = random.randint(0, num-1)
        mutation_seq = possible_seqs[ii]

        if file_type == "fastq":
            return SeqRecord(Seq(mutation_seq), id=records_dict[idx].id, description=records_dict[idx].description, letter_annotations=records_dict[idx].letter_annotations)
        elif file_type == "fasta":
            return SeqRecord(Seq(mutation_seq).seq, id=records_dict[idx].id, description=records_dict[idx].description)  


    def enumerate_ed2_seqs(self, read):
        possible_ed1 = seq2substitution(read)
        possible_ed2 = []
        for seq in possible_ed1:
            possible_ed2.extend(list(set(seq2substitution(seq)) - possible_ed1))
        return list(set(possible_ed2))

    '''
    def read2_2baseerr_read(self, read):
        possible_seqs = self.enumerate_ed2_seqs(read)
        num = len(possible_seqs)
        random.seed(self.config.sim_seed)
        random.shuffle(possible_seqs)
        idx = random.randint(0, num-1)
        return possible_seqs[idx]
        
    def extract_seq(self, shared_objects, name):
        records_dict, file_type = shared_objects
        if file_type == "fastq":
            return SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description, letter_annotations=records_dict[name].letter_annotations)
        elif file_type == "fasta":
            return SeqRecord(records_dict[name].seq, id=records_dict[name].id, description=records_dict[name].description)

    def extract_seqs(self, num_workers, name_lst, data_set, sub_dataset):
        records_dict, file_type = self.parse_data_dict(data_set)
        sub_records = []
        shared_objects = records_dict, file_type
        with WorkerPool(num_workers, shared_objects=shared_objects, start_method='fork') as pool:
            with tqdm(total=len(name_lst), desc=self.logger.info("Extract Seqs")) as pbar:   
                for tmp_rec in pool.imap(self.extract_seq, name_lst):
                    # print(tmp_rec)
                    sub_records.append(tmp_rec)
                    pbar.update() 
        # print(type(sub_records[0]), sub_records[0])
        with open(sub_dataset, "w") as handle:
            SeqIO.write(sub_records, handle, file_type)
        return 
    '''
    '''
    def error_injection(self, data_set):
        records_dict, file_type = parse_data_dict(data_set)
        id_lst = list(records_dict.keys())

        total_seqs_lst = []
        seq2id_dict = {}
        for item in id_lst:
            seq = str(records_dict[item].seq)
            seq2id_dict.setdefault(seq, []).append(item)
            total_seqs_lst.append(seq)
            # total_bases += len(seq)
            # total_reads += 1

        read2counts = collections.Counter(total_seqs_lst)
        select_id_lst = []
        total_bases = 0
        total_reads = 0
        new_read2counts = {}
        for read, count in read2counts.items():
            if count >= 5:
                select_id_lst.extend(seq2id_dict[read])
            if count >= self.config.min_read_count:
                new_read2counts[read] = count
                
                total_reads += count
                total_bases += len(read) * count

        err_reads_num = round(total_bases * self.config.error_rate1)
    
        if total_reads > err_reads_num:
            per_num =  err_reads_num / total_reads
        err_id_lst = []
        for read, count in new_read2counts.items():
            random_num = round(count * per_num)
            name_lst = seq2id_dict[read]
            random.shuffle(name_lst)
            err_id_lst.extend(random.sample(name_lst, random_num))

        # err_reads_num = round(total_bases * self.config.error_rate1 / total_reads)
        self.logger.info(f"total bases: {total_bases}, total reads: {total_reads}, error rate: {self.config.error_rate1}, calculated error reads number: {err_reads_num}, actual error reads number: {len(err_id_lst)}")

        err_subdataset = self.config.result_dir + "err_subdataset" + file_type
        err_records = []
        shared_objects = records_dict, file_type
        with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
            for tmp_rec in pool.imap(self.mutate_seq, err_id_lst, progress_bar=self.config.verbose):
                err_records.append(tmp_rec)

        # with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
        #     with tqdm(total=len(err_id_lst), desc=self.logger.info("Generate Mutated Seqs")) as pbar:   
        #         for tmp_rec in pool.imap(self.mutate_seq, err_id_lst):
        #             # print(tmp_rec)
        #             err_records.append(tmp_rec)
        #             pbar.update() 
        # print(type(sub_records[0]), sub_records[0])
        with open(err_subdataset, "w") as handle:
            SeqIO.write(err_records, handle, file_type)
        ########################################################################################
        err_free_subdataset = self.config.result_dir + "err_free_subdataset" + file_type

        err_free_id_lst = list(set(select_id_lst) - set(err_id_lst))
        extract_records(self.config.result_dir, err_free_id_lst, data_set, err_free_subdataset)

        base = self.config.input_file.split("/")[-1]
        file_name = base.split(file_type)[0]
        raw_data = self.config.result_dir + "simulated_raw_" + file_name + file_type
        os.system("cat %s %s > %s" % (err_subdataset, err_free_subdataset, raw_data))
        ##################################################################################
        err_free_subdataset_2 = self.config.result_dir + "err_free_subdataset2" + file_type
        extract_records(self.config.result_dir, err_id_lst, data_set, err_free_subdataset_2)
        true_data = self.config.result_dir + "simulated_true_" + file_name + file_type  
        os.system("cat %s %s > %s" % (err_free_subdataset_2, err_free_subdataset, true_data))     
        return raw_data, true_data
    '''

    '''
    def mutate_seq(self, shared_objects, idx):
        records_dict, file_type = shared_objects
        ori_seq = str(records_dict[idx].seq)
        mutation_seq = self.read2err_read(ori_seq)
        if file_type == "fastq":
            # dis = editdistance.eval(ori_seq, mutation_seq)
            # if dis == 0:
            #     pass
            # elif dis == 1:
            #     indel_idx = self.find_first_diff_char(ori_seq, mutation_seq)
            #     # Get the quality scores list from letter_annotations
            #     quality_scores = records_dict[idx].letter_annotations["phred_quality"]
            #     if len(ori_seq) < len(mutation_seq):
            #         # Calculate the average score
            #         min_score = min(quality_scores)
            #         # Insert the average score at the specified index
            #         quality_scores.insert(indel_idx, min_score)
            #         records_dict[idx].letter_annotations["phred_quality"] = quality_scores
            #     elif len(ori_seq) > len(mutation_seq):
            #         quality_scores.pop(indel_idx)
            #         records_dict[idx].letter_annotations["phred_quality"] = quality_scores
            # else:
            #     self.logger.error("Generate non 1-base error read.")
            #     sys.exit(1)
            return SeqRecord(Seq(mutation_seq), id=records_dict[idx].id, description=records_dict[idx].description, letter_annotations=records_dict[idx].letter_annotations)
        elif file_type == "fasta":
            return SeqRecord(Seq(mutation_seq).seq, id=records_dict[idx].id, description=records_dict[idx].description) 
    '''

    ''' 
    def error_injection(self, data_set):
        records_dict, file_type = parse_data_dict(data_set)
        id_lst = list(records_dict.keys())
        total_seqs_lst = []
        seq2id_dict = {}
        for item in id_lst:
            seq = str(records_dict[item].seq)
            seq2id_dict.setdefault(seq, []).append(item)
            total_seqs_lst.append(seq)
            # total_bases += len(seq)
            # total_reads += 1
            
        read2counts = collections.Counter(total_seqs_lst)
        select_id_lst = []
        # total_bases = 0
        total_reads = 0
        new_read2counts = {}
        for read, count in read2counts.items():
            # if count >= 5:
            if count > self.config.min_freq:
                select_id_lst.extend(seq2id_dict[read])
            if count > self.config.min_read_count:
                new_read2counts[read] = count
                total_reads += count
                # total_bases += len(read) * count
                
        # one base error
        err_1_base_num = round(total_reads * self.config.error_rate1)
        if total_reads > err_1_base_num:
            per_num =  err_1_base_num / total_reads
            err_1_id_lst = []
            for read, count in new_read2counts.items():
                random_num = round(count * per_num)
                name_lst = seq2id_dict[read]

                rng = np.random.default_rng(seed=self.config.sim_seed)
                indices = rng.permutation(len(name_lst))[:random_num]
                # print(indices)
                err_1_id_lst.extend(name_lst[i] for i in indices)
                # print(err_1_id_lst)
        else:
            self.logger.error("1 base error reads number is larger than the total reads number!")
            sys.exit(1)

        # two base errors
        if self.config.error_rate2 > 0:
            err_2_id_lst = []
            cur_total_reads = total_reads - len(err_1_id_lst)
            err_2_base_num = round(cur_total_reads * self.config.error_rate2)

            if cur_total_reads > err_2_base_num:
                per_num_2base =  err_2_base_num / (cur_total_reads * 2)
                for read, count in new_read2counts.items():
                    random_num_2base = round(count * per_num_2base)
                    cur_name_lst = list(set(seq2id_dict[read]) - set(err_1_id_lst))

                    rng = np.random.default_rng(seed=self.config.sim_seed)
                    indices = rng.permutation(len(cur_name_lst))[:random_num_2base]
                    err_2_id_lst.extend(cur_name_lst[i] for i in indices)
                    
                self.logger.info(f"Total reads: {total_reads}, 1 base error per read error rate: {self.config.error_rate1}, 1 base error reads number: {len(err_1_id_lst)}, 2 bases per read error rate: {self.config.error_rate2}, 2 bases error reads number: {len(err_2_id_lst)}")
            else:
                self.logger.error("2 bases error reads number is larger than the total reads number!")
                sys.exit(1)
        else:
            # err_reads_num = round(total_bases * self.config.error_rate / total_reads)
            self.logger.info(f"Total reads: {total_reads}, 1 base error per read error rate: {self.config.error_rate1}, 1 base error reads number: {len(err_1_id_lst)}")
        
        err_1base_records = []
        shared_objects = records_dict, file_type
        with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
            for tmp_rec in pool.imap(self.mutate_seq, err_1_id_lst, progress_bar=self.config.verbose):
                err_1base_records.append(tmp_rec)

        err_2base_records = []
        shared_objects = records_dict, file_type
        with WorkerPool(self.config.num_workers, shared_objects=shared_objects, start_method='fork') as pool:
            for tmp_rec in pool.imap(self.mutate_2base_seq, err_2_id_lst, progress_bar=self.config.verbose):
                err_2base_records.append(tmp_rec)

        err_records = err_1base_records + err_2base_records
        err_subdataset = self.config.result_dir + "err_subdataset" + file_type
        with open(err_subdataset, "w") as handle:
            SeqIO.write(err_records, handle, file_type)

        ########################################################################################
        err_free_subdataset = self.config.result_dir + "err_free_subdataset" + file_type

        err_id_lst = err_1_id_lst + err_2_id_lst
        err_free_id_lst = list(set(select_id_lst) - set(err_id_lst))
        extract_records(self.config.result_dir, err_free_id_lst, data_set, err_free_subdataset)

        base = self.config.input_file.split("/")[-1]
        file_name = base.split(file_type)[0]
        raw_data = self.config.result_dir + "simulated_raw_" + file_name + file_type
        os.system("cat %s %s > %s" % (err_subdataset, err_free_subdataset, raw_data))
        ##################################################################################
        err_free_subdataset_2 = self.config.result_dir + "err_free_subdataset2" + file_type
        extract_records(self.config.result_dir, err_id_lst, data_set, err_free_subdataset_2)
        true_data = self.config.result_dir + "simulated_true_" + file_name + file_type  
        os.system("cat %s %s > %s" % (err_free_subdataset_2, err_free_subdataset, true_data))     
        return raw_data, true_data
    '''

    '''
    total_seqs_lst = []
    seq2id_dict = {}
    for item in id_lst:
        seq = str(records_dict[item].seq)
        seq2id_dict.setdefault(seq, []).append(item)
        total_seqs_lst.append(seq)
        # total_bases += len(seq)
        # total_reads += 1

    read2counts = collections.Counter(total_seqs_lst)
    select_id_lst = []
    total_bases = 0
    total_reads = 0
    new_read2counts = {}
    for read, count in read2counts.items():
        if count >= 5:
            select_id_lst.extend(seq2id_dict[read])
        if count >= self.config.min_read_count:
            new_read2counts[read] = count
            
            total_reads += count
            total_bases += len(read) * count

    err_reads_num = round(total_bases * self.config.error_rate1)

    if total_reads > err_reads_num:
        per_num =  err_reads_num / total_reads
    err_1_id_lst = []
    for read, count in new_read2counts.items():
        random_num = round(count * per_num)
        name_lst = seq2id_dict[read]
        random.shuffle(name_lst)
        err_1_id_lst.extend(random.sample(name_lst, random_num))

    # err_reads_num = round(total_bases * self.config.error_rate1 / total_reads)
    self.logger.info(f"total bases: {total_bases}, total reads: {total_reads}, error rate: {self.config.error_rate1}, calculated error reads number: {err_reads_num}, actual error reads number: {len(err_1_id_lst)}")
    '''


# if __name__ == '__main__':
#     logger = custom_logger("simulation", debug_mode=True)
#     input_file = "/home/pping/Data/Repo/data/noise2read_data/group3/SAS-Cov-2/noise2read_result/new_no_ml/output/umi_simulated_raw_original_sars_r1_corrected_correct.fastq"
#     sim = Simulation(logger, 30, )

#     raw_data, true_data = sim.error_injection(input_file)
#     output_dir = "/home/pping/Data/Repo/data/noise2read_data/group3/SAS-Cov-2/"
#     DP = DataProcessing( 
#         logger,
#         num_workers=30,
#         output_dir=output_dir)
#     umi_raw_dataset, umi_true_dataset = DP.write_mimic_umis(raw_data, true_data)