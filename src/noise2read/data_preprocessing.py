# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-09-07 16:59:49

import collections
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from noise2read.utils import *
from noise2read.data_generation import DataGneration
from noise2read.error_orrection import ErrorCorrection
# from noise2read.data_analysis import DataAnalysis
import editdistance
from tqdm import tqdm
from Bio.Seq import Seq
import json

class DataProcessing():
    # def __init__(self, logger, num_workers, output_dir, umi_start=None, umi_end=None, non_umi_start=None):
    def __init__(self, logger, config):
        self.logger = logger
        self.config = config
        
        if os.path.exists(self.config.result_dir):
            self.logger.info("Directory '%s' already exists" % self.config.result_dir)
        else:
            os.makedirs(self.config.result_dir)
        if os.path.exists(self.config.result_dir + 'raw/'):
            self.logger.info(f"Directory {self.config.result_dir}raw/ already exists")
        else:
            os.makedirs(self.config.result_dir + 'raw/')       
        if os.path.exists(self.config.result_dir + 'true/'):
            self.logger.info(f"Directory {self.config.result_dir}true/ already exists")
        else:
            os.makedirs(self.config.result_dir + 'true/') 
               
    def extract_umis(self, original_data, umi_start, umi_end, non_umi_start):
        record_iterator, ff_type = parse_data(original_data) 
        self.logger.info("Read Data")
        umi_records = []
        non_umi_records = []

        for item in tqdm(record_iterator):
            seq = str(item.seq)
            if 'N' not in seq:
                umi_qual = {}
                non_umi_qual = {}
                if ff_type == 'fastq' or ff_type == 'fq' or ff_type == 'fastq.gz' or ff_type == 'fq.gz':
                    val = item.letter_annotations['phred_quality']            
                umi_desc = "length=" + str(umi_end-umi_start)
                non_umi_desc = "length=" + str(len(seq) - non_umi_start)
                if ff_type == "fastq": 
                    umi_qual['phred_quality'] = val[umi_start:umi_end]
                    non_umi_qual['phred_quality'] = val[non_umi_start:]
                    umi_rec = SeqRecord(Seq(seq[umi_start:umi_end]), id=item.id, description=umi_desc, letter_annotations=umi_qual)   
                    non_umi_rec = SeqRecord(Seq(seq[non_umi_start:]), id=item.id, description=non_umi_desc, letter_annotations=non_umi_qual)   
                elif ff_type == "fasta":
                    umi_rec = SeqRecord(Seq(seq[umi_start:umi_end]), id=item.id, description=umi_desc)        
                    non_umi_rec = SeqRecord(Seq(seq[non_umi_start:]), id=item.id, description=non_umi_desc)  
                umi_records.append(umi_rec)
                non_umi_records.append(non_umi_rec)        
        self.logger.info("Rewrite umi and non-umi Data")
        raw_base_name = original_data.split('/')[-1]
        umi_raw_dataset = self.config.result_dir + "umi_" + raw_base_name
        non_umi_raw_dataset = self.config.result_dir + "umi_in_name_" + raw_base_name
        with open(umi_raw_dataset, "w") as handle:
            SeqIO.write(umi_records, handle, ff_type)        
        with open(non_umi_raw_dataset, "w") as handle:
            SeqIO.write(non_umi_records, handle, ff_type) 
        self.logger.info("Split umi and non-umis completed.") 
        return umi_raw_dataset, non_umi_raw_dataset

    def real_umi_data(self, original_data):
        # step1: extract umis to a sperate fastq file
        umi_raw_dataset, non_umi_raw_dataset = self.extract_umis(original_data, self.config.umi_start, self.config.umi_end, self.config.non_umi_start)
        # step2: correct umi errors
        # self.config.input_file = umi_raw_dataset
        # DG = DataGneration(self.logger, self.config)

        # genuine_df = DG.extract_umi_genuine_errs(umi_raw_dataset)
        # # ##############################################################
        # EC = ErrorCorrection(self.logger, self.config)
        # corrected_file = EC.umi_correction(umi_raw_dataset, genuine_df)

        # base_name = umi_raw_dataset.split("/")[-1].split(".fastq")[0]
        # step3: use umi and real dataset to generate raw and true dataset
        # raw_file, true_file = self.raw_true_umi(corrected_file, non_umi_raw_dataset)
        raw_file, true_file = self.raw_true_umi(umi_raw_dataset, non_umi_raw_dataset)
        return raw_file, true_file

    def raw_true_umi(self, correct_umi_f, non_umi_f):
        umi2id = {}
        umi_record_iterator, umi_f_type = parse_data(correct_umi_f)
        for item in umi_record_iterator:
            umi = str(item.seq)
            umi2id.setdefault(umi, []).append(str(item.id))

        id2read, ori_file_type = parse_data_dict(non_umi_f)
        raw_records = []
        true_records = []
        for umi in umi2id:
            ids = umi2id[umi]
            # n = 0
            tmp_raw_rec = []
            tmp_true_rec = []
            if len(ids) >= self.config.group_read_number:
                cur_seqs = []
                for id in ids:
                    cur_seqs.append(str(id2read[id].seq))
                seq2count = collections.Counter(cur_seqs)
                max_key = max(seq2count, key=seq2count.get)    
                total_count = sum(list(seq2count.values()))
                max_key_count = seq2count[max_key]
                if max_key_count >= (total_count-max_key_count):     
                    for i in ids:
                        cur_seq = str(id2read[i].seq)
                        dis = editdistance.eval(max_key, cur_seq)
                        desc = "umi:" + str(umi) + "//" + str(id2read[i].description)
                        if dis == 0:
                            # n += 1
                            if ori_file_type == "fastq":
                                rec = SeqRecord(Seq(cur_seq), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                            else:
                                rec = SeqRecord(Seq(cur_seq), id=i, description=desc)
                            # raw_records.append(rec)
                            # true_records.append(rec)
                            tmp_raw_rec.append(rec)
                            tmp_true_rec.append(rec)
                        elif dis <= self.config.read_edit_dif:
                        # else:
                            # n += 1
                            if ori_file_type == "fastq":
                                raw_rec = SeqRecord(Seq(cur_seq), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                                true_rec = SeqRecord(Seq(max_key), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                            else:
                                raw_rec = SeqRecord(Seq(cur_seq), id=i, description=desc)
                                true_rec = SeqRecord(Seq(max_key), id=i, description=desc)
                            # raw_records.append(raw_rec)
                            # true_records.append(true_rec)
                            tmp_raw_rec.append(raw_rec)
                            tmp_true_rec.append(true_rec)
                    # if n >= 10:
                    raw_records.extend(tmp_raw_rec)
                    true_records.extend(tmp_true_rec)                    
        self.logger.info("Write raw and true data to file")
        raw_file = self.config.result_dir + "raw/" + non_umi_f.split('/')[-1]
        true_file = self.config.result_dir + "true/" + non_umi_f.split('/')[-1]
        with open(raw_file, "w") as handle:
            SeqIO.write(raw_records, handle, ori_file_type)        
        with open(true_file, "w") as handle:
            SeqIO.write(true_records, handle, ori_file_type) 
        self.logger.info("Real umi data processing completed.")                 

        return raw_file, true_file

    def raw_true_umi_1(self, correct_umi_f, non_umi_f):
        umi2id = {}
        umi_record_iterator, umi_f_type = parse_data(correct_umi_f)
        umi_n = 0
        for item in umi_record_iterator:
            umi = str(item.seq)
            if 'N' in umi:
                pass
            else:
                umi_n += 1
                umi2id.setdefault(umi, []).append(str(item.id))
        print("test")
        print(umi_n)
        id2read, ori_file_type = parse_data_dict(non_umi_f)
        # print(list(id2read.keys()))
        raw_records = []
        true_records = []
        n = 0
        nn = 0
        for umi in umi2id:
            seq_ids = umi2id[umi]
            tmp_raw_rec = []
            tmp_true_rec = []
            if len(seq_ids) >= 4:
                n += 1
                # print(ids[0])
                # print("#")
                # print(str(id2read[ids[0]].seq))
                cur_seqs = []
                for seq_id in seq_ids:
                    cur_seqs.append(str(id2read[seq_id].seq))
                seq2count = collections.Counter(cur_seqs)
                max_key = max(seq2count, key=seq2count.get)    
                
                total_count = sum(list(seq2count.values()))
                max_key_count = seq2count[max_key]
                print(max_key_count)
                sorted_lst = sorted(list(seq2count.values()), reverse=True)

                # if max_key_count >= (total_count-max_key_count):  
                # if max_key_count >= 0.6 * total_count:  
                if len(sorted_lst) >= 2:
                    if sorted_lst[1] < max_key_count:# and max_key_count >= (total_count-max_key_count):   
                    # if max_key_count >= (total_count-max_key_count):
                    # if max_key_count >= 4:
                        for read_id in seq_ids:
                            cur_seq = str(id2read[read_id].seq)
                            dis = editdistance.eval(max_key, cur_seq)
                            desc = "umi:" + str(umi) + "//" + str(id2read[read_id].description)
                            if dis == 0:
                                nn += 1
                                if ori_file_type == "fastq":
                                    rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                else:
                                    rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc)
                                # raw_records.append(rec)
                                # true_records.append(rec)
                                tmp_raw_rec.append(rec)
                                tmp_true_rec.append(rec)
                            # elif dis == 1 or dis == 2:
                            # elif dis <= 10: #and seq2count[cur_seq] < max_key_count
                            else:
                                nn += 1
                                if ori_file_type == "fastq":
                                    raw_rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                    true_rec = SeqRecord(Seq(max_key), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                else:
                                    raw_rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc)
                                    true_rec = SeqRecord(Seq(max_key), id=read_id, description=desc)
                                # raw_records.append(raw_rec)
                                # true_records.append(true_rec)
                                tmp_raw_rec.append(raw_rec)
                                tmp_true_rec.append(true_rec)
                            if dis > 8:
                                print(seq2count.values())
                            # if n >= 10:
                            # raw_records.extend(tmp_raw_rec)
                            # true_records.extend(tmp_true_rec)
                elif len(sorted_lst) == 1:
                    for i in seq_ids:
                        cur_seq = str(id2read[i].seq)
                        dis = editdistance.eval(max_key, cur_seq)
                        desc = "umi:" + str(umi) + "//" + str(id2read[i].description)
                        if dis == 0:
                            n += 1
                            if ori_file_type == "fastq":
                                rec = SeqRecord(Seq(cur_seq), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                            else:
                                rec = SeqRecord(Seq(cur_seq), id=i, description=desc)
                            # raw_records.append(rec)
                            # true_records.append(rec)
                            tmp_raw_rec.append(rec)
                            tmp_true_rec.append(rec)
                        # elif dis == 1 or dis == 2:
                        # elif dis <= 4:
                        else:
                            n += 1
                            if ori_file_type == "fastq":
                                raw_rec = SeqRecord(Seq(cur_seq), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                                true_rec = SeqRecord(Seq(max_key), id=i, description=desc, letter_annotations=id2read[i].letter_annotations)
                            else:
                                raw_rec = SeqRecord(Seq(cur_seq), id=i, description=desc)
                                true_rec = SeqRecord(Seq(max_key), id=i, description=desc)
                            # raw_records.append(raw_rec)
                            # true_records.append(true_rec)
                            tmp_raw_rec.append(raw_rec)
                            tmp_true_rec.append(true_rec)
                        # if dis > 8:
                        #     print(seq2count.values())
                    # if n >= 10:
            if tmp_raw_rec and tmp_true_rec:
                raw_records.extend(tmp_raw_rec)
                true_records.extend(tmp_true_rec)                                        
        print(n)
        print(nn)
        self.logger.info("Write raw and true data to file")
        raw_file = self.config.result_dir + "raw/" + non_umi_f.split('/')[-1]
        true_file = self.config.result_dir + "true/" + non_umi_f.split('/')[-1]
        with open(raw_file, "w") as handle:
            SeqIO.write(raw_records, handle, ori_file_type)        
        with open(true_file, "w") as handle:
            SeqIO.write(true_records, handle, ori_file_type) 
        self.logger.info("Real umi data processing completed.")                 

        return raw_file, true_file

    def umi_cluster_analysis(self, correct_umi_f, non_umi_f):
        umi2id = {}
        umi_record_iterator, umi_f_type = parse_data(correct_umi_f)
        umi_n = 0
        for item in umi_record_iterator:
            umi = str(item.seq)
            if 'N' in umi:
                pass
            else:
                umi_n += 1
                umi2id.setdefault(umi, []).append(str(item.id))
        print("test")
        print(umi_n)
        id2read, ori_file_type = parse_data_dict(non_umi_f)
        # print(list(id2read.keys()))
        raw_records = []
        true_records = []
        n = 0
        nn = 0
        nnn = 0
        high_seq_umi_number = 0

        umi2high_read2lowread_edit_dis = {}
        umi2high_read2lowread_edit_dis2 = {}
        umi2high_read2lowread_edit_dis3 = {}
        
        read_editDis2count = {}

        for umi in umi2id:
            seq_ids = umi2id[umi]
            # tmp_raw_rec = []
            # tmp_true_rec = []
            high_read2lowread_edit_dis = {}
            high_read2lowread_edit_dis2 = {}
            if len(seq_ids) >= self.config.group_read_number:
                n += 1
                # print(ids[0])
                # print("#")
                # print(str(id2read[ids[0]].seq))
                cur_seqs = []
                for seq_id in seq_ids:
                    cur_seqs.append(str(id2read[seq_id].seq))
                seq2count = collections.Counter(cur_seqs)
                # sorted_lst = sorted(list(seq2count.values()), reverse=True)
                if len(seq2count) >= 2:
                    # seq2count_above = {key: value for key, value in seq2count.items() if value > 4}
                    high_seq2count = {}
                    low_seq2count = {}
                    
                    for seq, count in seq2count.items():
                        if count > self.config.high_freq_thre:
                            high_seq2count[seq] = count
                        else:
                            low_seq2count[seq] = count
                    if len(high_seq2count)==0:
                        high_seq_umi_number += 1
                    # print(high_seq2count)
                    if len(high_seq2count) > 0:
                        
                        high_seqs = list(high_seq2count.keys())

                        for read_id in seq_ids:
                            cur_seq = str(id2read[read_id].seq)
                            desc = "umi:" + str(umi) + "//" + str(id2read[read_id].description)
                            if cur_seq in high_seqs:
                                nn += 1
                                if ori_file_type == "fastq":
                                    rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                else:
                                    rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc)
                                raw_records.append(rec)
                                true_records.append(rec)
                                del rec
                                # tmp_raw_rec.append(rec)
                                # tmp_true_rec.append(rec) 
                            else:  
                                sorted_seq2count_above = dict(sorted(high_seq2count.items(), key=lambda item: item[1]), reverse=True)  
                                # print(sorted_seq2count_above)                       
                                for high_read, high_count in sorted_seq2count_above.items():
                                    dis = editdistance.eval(high_read, cur_seq)
                                    if dis <= self.config.read_edit_dif:
                                    # if dis <= 4:
                                        nn += 1
                                        if ori_file_type == "fastq":
                                            raw_rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                            true_rec = SeqRecord(Seq(high_read), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                                        else:
                                            raw_rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc)
                                            true_rec = SeqRecord(Seq(high_read), id=read_id, description=desc)
                                        raw_records.append(raw_rec)
                                        true_records.append(true_rec)
                                        del raw_rec, true_rec
                                        # tmp_raw_rec.append(raw_rec)
                                        # tmp_true_rec.append(true_rec) 
                                        break 
                        # for high_read, high_count in high_seq2count.items():
                        for low_seq, low_seq_count in low_seq2count.items():
                            high_seq, edit_dis = self.get_smallest_edit_dis_seq(high_seqs, low_seq)

                            high_read2lowread_edit_dis.setdefault((high_seq, high_seq2count[high_seq]), []).append([low_seq, low_seq_count, edit_dis]) 
                            high_read2lowread_edit_dis2.setdefault((high_seq, high_seq2count[high_seq]), []).append(edit_dis)
                        new_high_read2lowread_edit_dis2 = {}

                        for key, val in high_read2lowread_edit_dis2.items():
                            edit_dis2count = collections.Counter(val)
                            new_high_read2lowread_edit_dis2[key] = edit_dis2count
                            for k, value in edit_dis2count.items():
                                if k in read_editDis2count:
                                    read_editDis2count[k] += value
                                else:
                                    read_editDis2count[k] = value                            
                            
                        new_high_read2lowread_edit_dis = {}
                        for key, val in high_read2lowread_edit_dis.items():
                            temp_dict = {}
                            for item in val:
                                edit_dis = item[2]
                                if edit_dis in temp_dict:
                                    temp_dict[edit_dis] += item[1]
                                else:
                                    temp_dict[edit_dis] = item[1]
                            new_high_read2lowread_edit_dis[key] = temp_dict

                        if len(high_seqs) >= 2:
                            edit_dis_high = self.pairwise_edit_distances(high_seqs)
                            new_high_read2lowread_edit_dis["edit_dis_between_high_freq_read"] = edit_dis_high

                    else:
                        nnn += sum(list(seq2count.values()))
                elif len(seq2count) == 1:
                    for read_id in seq_ids:
                        cur_seq = str(id2read[read_id].seq)
                        desc = "umi:" + str(umi) + "//" + str(id2read[read_id].description)
                        nn += 1
                        if ori_file_type == "fastq":
                            rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc, letter_annotations=id2read[read_id].letter_annotations)
                        else:
                            rec = SeqRecord(Seq(cur_seq), id=read_id, description=desc)
                        # tmp_raw_rec.append(rec)
                        # tmp_true_rec.append(rec)
                        raw_records.append(rec)
                        true_records.append(rec)
                        del rec

            if len(high_read2lowread_edit_dis)>0:
                data_str_keys = {str(key): value for key, value in high_read2lowread_edit_dis.items()}
                
                umi2high_read2lowread_edit_dis[umi] = data_str_keys

                data_str_keys2 = {str(key): value for key, value in new_high_read2lowread_edit_dis2.items()}
                umi2high_read2lowread_edit_dis2[umi] = data_str_keys2

                data_str_keys3 = {str(key): value for key, value in new_high_read2lowread_edit_dis.items()}
                umi2high_read2lowread_edit_dis3[umi] = data_str_keys3



            # if tmp_raw_rec and tmp_true_rec:
            #     raw_records.extend(tmp_raw_rec)
            #     true_records.extend(tmp_true_rec) 
        new_read_editDis2count = {str(key): value for key, value in read_editDis2count.items()}
        with open(self.config.result_dir + non_umi_f.split('/')[-1] + "1.json", "w") as json_file:
            json.dump(umi2high_read2lowread_edit_dis, json_file, indent=4)
        with open(self.config.result_dir + non_umi_f.split('/')[-1] + "2.json", "w") as json_file:
            json.dump(umi2high_read2lowread_edit_dis2, json_file, indent=4)
        with open(self.config.result_dir + non_umi_f.split('/')[-1] + "3.json", "w") as json_file:
            json.dump(umi2high_read2lowread_edit_dis3, json_file, indent=4)
            
        with open(self.config.result_dir + non_umi_f.split('/')[-1] + "edit_dis2number.json", "w") as json_file:
            json.dump(new_read_editDis2count, json_file, indent=4)

        print(n)
        print(nn)
        print(nnn)
        print(high_seq_umi_number)
        self.logger.info("Write raw and true data to file")
        raw_file = self.config.result_dir + "raw/" + non_umi_f.split('/')[-1]
        true_file = self.config.result_dir + "true/" + non_umi_f.split('/')[-1]
        with open(raw_file, "w") as handle:
            SeqIO.write(raw_records, handle, ori_file_type)        
        with open(true_file, "w") as handle:
            SeqIO.write(true_records, handle, ori_file_type) 
        self.logger.info("Real umi data processing completed.")                 

        return raw_file, true_file

    def pairwise_edit_distances(self, strings):
        n = len(strings)
        distances = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                distance = editdistance.eval(strings[i], strings[j])
                distances[i][j] = distance
                distances[j][i] = distance

        return distances

    def get_smallest_edit_dis_seq(self, seq_lst, read):
        smallest_edit_dis = 100000
        target_read = ""
        for seq in seq_lst:
            dis = editdistance.eval(seq, read)
            if dis < smallest_edit_dis:
                smallest_edit_dis = dis
                target_read = seq
        return target_read, smallest_edit_dis


    def write_mimic_umis(self, raw_dataset, true_dataset):
        true_records_dict, true_file_type = parse_data_dict(true_dataset)

        true_seqs_lst = []
        # raw_seqs_lst = []
        for name in true_records_dict:
            true_seq = str(true_records_dict[name].seq)
            # raw_seq = str(error_records[name].seq)
            true_seqs_lst.append(true_seq)     
            # raw_seqs_lst.append(raw_seq)
        # raw_read2count = collections.Counter(raw_seqs_lst)
        true_read2count = collections.Counter(true_seqs_lst)
        # write the same read the same id in true dataset, these ids also write to raw dataset
        true_read2id = {}
        i = 1
        for read in true_read2count:
            true_read2id[read] = i
            i += 1
        true_records = []
        raw_records = []
        raw_records_dict, raw_file_type = parse_data_dict(raw_dataset)
        self.logger.info("Read Data")
        for name in tqdm(true_records_dict):
            true_seq = str(true_records_dict[name].seq)
            true_des = "umi:" + str(true_read2id[true_seq]) + "//" + str(true_records_dict[name].description)
            raw_des = "umi:" + str(true_read2id[true_seq]) + "//" + str(raw_records_dict[name].description)
            if true_file_type == "fastq": 
                true_rec = SeqRecord(true_records_dict[name].seq, id=true_records_dict[name].id, description=true_des, letter_annotations=true_records_dict[name].letter_annotations)   
                raw_rec = SeqRecord(raw_records_dict[name].seq, id=raw_records_dict[name].id, description=raw_des, letter_annotations=raw_records_dict[name].letter_annotations)   
            elif true_file_type == "fasta":
                true_rec = SeqRecord(true_records_dict[name].seq, id=true_records_dict[name].id, description=true_des)        
                raw_rec = SeqRecord(raw_records_dict[name].seq, id=raw_records_dict[name].id, description=raw_des)  
            raw_records.append(raw_rec)
            true_records.append(true_rec)
        raw_base_name = raw_dataset.split('/')[-1]
        true_base_name = true_dataset.split('/')[-1]
        umi_raw_dataset = self.config.result_dir + "raw/" + "umi_" + raw_base_name
        umi_true_dataset = self.config.result_dir + "true/" + "umi_"  + true_base_name
        self.logger.info("Rewrite Data")
        with open(umi_raw_dataset, "w") as handle:
            SeqIO.write(raw_records, handle, raw_file_type)        
        with open(umi_true_dataset, "w") as handle:
            SeqIO.write(true_records, handle, true_file_type) 
        self.logger.info("Mimic_umi preprocessing completed.") 
        return umi_raw_dataset, umi_true_dataset
    
    def extract_umis2fasta(self):
        record_iterator, ff_type = parse_data(self.config.input_file) 
        self.logger.info("Read Data")
        umi_records = []
        # non_umi_records = []

        for item in tqdm(record_iterator):
            # seq = str(item.seq)
            seq_id = str(item.id)
            seq_des = str(item.description)
            # print(seq_des)
            umi = seq_des.split(self.config.separator1)[self.config.separator1_idx].split(self.config.separator2)[self.config.separator2_idx]

            umi_rec = SeqRecord(Seq(umi), id=seq_id, description=seq_des)        
            umi_records.append(umi_rec)
                   
        self.logger.info("Rewrite umi and non-umi Data")
        raw_base_name = self.config.input_file.split('/')[-1].split('.')[0]
        umi_raw_dataset = self.config.result_dir + "umi_" + raw_base_name + '.fasta'
        with open(umi_raw_dataset, "w") as handle:
            SeqIO.write(umi_records, handle, 'fasta')        
        self.logger.info("Extract umi to fasta completed.") 
        return umi_raw_dataset, self.config.input_file
    
    def real_umi_in_name_data(self):
        # step1: extract umis to a sperate fasta file
        # umi_raw_dataset, non_umi_raw_dataset = self.extract_umis(original_data, self.config.umi_start, self.config.umi_end, self.config.non_umi_start)
        umi_raw_dataset, non_umi_raw_dataset = self.extract_umis2fasta()
        # step2: correct umi errors
        # self.config.input_file = umi_raw_dataset
        DG = DataGneration(self.logger, self.config)

        genuine_df = DG.extract_umi_genuine_errs(umi_raw_dataset)
        # ##############################################################
        EC = ErrorCorrection(self.logger, self.config)

        corrected_file = EC.umi_correction(umi_raw_dataset, genuine_df)

        raw_file, true_file = self.raw_true_umi(corrected_file, non_umi_raw_dataset)
        
        return raw_file, true_file

    def umi2groundtruth(self):
        raw_file, true_file = self.raw_true_umi(self.config.umi_file, self.config.input_file)
        return raw_file, true_file