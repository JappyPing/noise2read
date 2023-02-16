# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-02-16 11:11:05

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

class DataProcessing():
    def __init__(self, logger, num_workers, output_dir, umi_start=None, umi_end=None, non_umi_start=None):
        self.logger = logger
        self.num_workers = num_workers
        self.output_dir = output_dir
        self.umi_start = umi_start
        self.umi_end = umi_end
        self.non_umi_start = non_umi_start
        if os.path.exists(self.output_dir):
            self.logger.info("Directory '% s' already exists" % self.output_dir)
        else:
            os.makedirs(self.output_dir)
        if os.path.exists(self.output_dir + 'raw/'):
            self.logger.info("Directory '% s' already exists" % self.output_dir + 'raw/')
        else:
            os.makedirs(self.output_dir + 'raw/')       
        if os.path.exists(self.output_dir + 'true/'):
            self.logger.info("Directory '% s' already exists" % self.output_dir + 'true/')
        else:
            os.makedirs(self.output_dir + 'true/') 
               
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
        umi_raw_dataset = self.output_dir + "umi_" + raw_base_name
        non_umi_raw_dataset = self.output_dir + "non_umi_" + raw_base_name
        with open(umi_raw_dataset, "w") as handle:
            SeqIO.write(umi_records, handle, ff_type)        
        with open(non_umi_raw_dataset, "w") as handle:
            SeqIO.write(non_umi_records, handle, ff_type) 
        self.logger.info("Split umi and non-umis completed.") 
        return umi_raw_dataset, non_umi_raw_dataset

    def real_umi_data(self, original_data):
        # step1: extract umis to a sperate fastq file
        umi_raw_dataset, non_umi_raw_dataset = self.extract_umis(original_data, self.umi_start, self.umi_end, self.non_umi_start)
        # step2: correct umi errors
        DG = DataGneration(
            self.logger,
            self.num_workers,
            umi_raw_dataset, 
            self.output_dir,
            high_freq_thre = 5,
            max_error_freq = 4,
            save_graph = False,
            graph_visualization = False,
            drawing_graph_num = False,
            high_ambiguous=False, 
            verbose=True
            )
        genuine_df = DG.extract_umi_genuine_errs()
        # ##############################################################
        EC = ErrorCorrection(
            self.logger,
            self.num_workers,
            umi_raw_dataset,
            self.output_dir,
            read_max_len = self.umi_end - self.umi_start,
            entropy_kmer=3, 
            entropy_q=2, 
            kmer_freq=3,
            read_type="DNA",
            iso_change_detail=False,
            min_iters=100)
        corrected_file = EC.umi_correction(umi_raw_dataset, genuine_df)

        # base_name = umi_raw_dataset.split("/")[-1].split(".fastq")[0]
        # step3: use umi and real dataset to generate raw and true dataset
        raw_file, true_file = self.raw_true_umi(corrected_file, non_umi_raw_dataset)
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
            if len(ids) >= 10:
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
                        elif dis == 1 or dis == 2:
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
        raw_file = self.output_dir + "raw/" + non_umi_f.split('/')[-1]
        true_file = self.output_dir + "true/" + non_umi_f.split('/')[-1]
        with open(raw_file, "w") as handle:
            SeqIO.write(raw_records, handle, ori_file_type)        
        with open(true_file, "w") as handle:
            SeqIO.write(true_records, handle, ori_file_type) 
        self.logger.info("Real umi data processing completed.")                 

        return raw_file, true_file

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
        umi_raw_dataset = self.output_dir + "raw/" + "umi_" + raw_base_name
        umi_true_dataset = self.output_dir + "true/" + "umi_"  + true_base_name
        self.logger.info("Rewrite Data")
        with open(umi_raw_dataset, "w") as handle:
            SeqIO.write(raw_records, handle, raw_file_type)        
        with open(umi_true_dataset, "w") as handle:
            SeqIO.write(true_records, handle, true_file_type) 
        self.logger.info("Mimic_umi preprocessing completed.") 
        return umi_raw_dataset, umi_true_dataset