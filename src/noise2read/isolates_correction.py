# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-02-16 11:11:32

from Bio import SeqIO
import os
import editdistance
import operator
from collections import Counter
from noise2read.utils import *
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from noise2read.utils import MemoryMonitor

class IsolatesErrorCorrection():
    def __init__(self, 
            logger, 
            num_workers,
            isolates, 
            correct_non_isolates,
            output_dir,
            iso_change_detail,
            min_iters):
        self.logger = logger
        self.num_workers = num_workers
        self.isolates = isolates
        self.correct_non_isolates = correct_non_isolates
        self.output_dir = output_dir
        bases = self.isolates.split('/')[-1]
        self.base = bases.split('.')
        self.iso_change_detail = iso_change_detail
        self.min_iters = min_iters
        self.MM = MemoryMonitor(logger)
        self.MM.start()

    def bcool_correct_isolates(self):
        # # bcool correction
        self.logger.info("Correcting Isolated nodes using bcool.")
        bcool_dir = self.output_dir + "bcool/input/"
        if not Path(bcool_dir).is_dir():
            os.makedirs(bcool_dir)
        iso_file_tye = parse_file_type(self.isolates)
        if iso_file_tye == "fastq":
            read_fa = bcool_dir + "input.fasta"
            # self.covert_fq2fa(self.isolates, read_fa)
            SeqIO.convert(self.isolates, "fastq", read_fa, "fasta-2line")
            # os.system("seqtk seq -a %s > %s" % (self.isolates, read_fa))
            os.system("bcool -u %s -t %s -o %s" % (read_fa, self.num_workers, bcool_dir))
        else:
            os.system("bcool -u %s -t %s -o %s" % (self.isolates, self.num_workers, bcool_dir))
        bcool_isolates = self.output_dir + "bcool/reads_corrected.fa"
        frequency_file = self.output_dir + self.base[0] + '_frequency.txt'
        self.MM.measure()
        corrected_isolates = self.select_correction(self.correct_non_isolates, self.isolates, bcool_isolates, self.output_dir + self.base[0], frequency_file)
        # os.system("rm %s" % bcool_isolates)
        self.logger.info("Isolated nodes Correction finished!")
        return corrected_isolates

    def freqency_seq_extraction(self, f_original, f_correct, frequency_file):
        correct_seq2name_dict = {}
        original_name2seq_dict = {}

        ori_record_iterator, _ = parse_data(f_original)
        cor_record_iterator, _ = parse_data(f_correct)
        for ori_rec, cor_rec in zip(ori_record_iterator, cor_record_iterator):
            ori_seq = str(ori_rec.seq)
            cor_seq = str(cor_rec.seq)

            ori_id = str(ori_rec.id)
            cor_id = str(cor_rec.id)

            original_name2seq_dict.setdefault(ori_id, []).append(ori_seq)
            correct_seq2name_dict.setdefault(cor_seq, []).append(cor_id)

        seq_freq = {}
        for seq, name_lst in correct_seq2name_dict.items():
            seq_freq[seq] = len(name_lst)

        seq_freq_order = dict(sorted(seq_freq.items(), key=operator.itemgetter(1), reverse=True))

        with open(frequency_file, 'a') as f:        
            f.write('Top frequency: ')
            f.write('\n')
            f.write(str(list(seq_freq_order.values())))
            f.write('\n')
            for seq, freq in seq_freq_order.items():
                f.write('##########################################')
                f.write('\n')
                f.write(str(freq))
                f.write('\n')
                f.write(seq)
                f.write('\n')
                f.write(str(correct_seq2name_dict[seq]))
                f.write('\n')
                for name in correct_seq2name_dict[seq]:
                    # print(name)
                    f.write(name)
                    f.write('\n')
                    f.write(str(original_name2seq_dict[name]))
                    f.write('\n')
                f.write('##########################################')
        del ori_record_iterator, cor_record_iterator, seq_freq_order, seq_freq, correct_seq2name_dict, original_name2seq_dict
        return

    # def bcool_seqs(self, num_workers, keep_bcool_name_lst, f_bcool, f_original, keep_correct_f):
    def bcool_seqs(self, keep_bcool_name_lst, f_bcool, f_original, keep_correct_f):
        bcool_records_dict, _ = parse_data_dict(f_bcool)
        ori_records_dict, file_type = parse_data_dict(f_original)
        sub_records = []
        for name in keep_bcool_name_lst:
            tmp_rec = SeqRecord(bcool_records_dict[name].seq, id=ori_records_dict[name].id, description=ori_records_dict[name].description, letter_annotations=ori_records_dict[name].letter_annotations)
            sub_records.append(tmp_rec)
        del bcool_records_dict, ori_records_dict
        # shared_objects = copy.deepcopy(bcool_records_dict, ori_records_dict)
        # workerpool and shared_objects use too large memory
        # name_len = len(keep_bcool_name_lst)
        # with WorkerPool(num_workers, shared_objects=shared_objects, start_method='fork') as pool:
        #     with tqdm(total=name_len, desc=self.logger.info("Extract Seqs"), miniters=int(name_len/self.min_iters)) as pbar:   
        #         for tmp_rec in pool.imap(self.bcool_seq, keep_bcool_name_lst):
        #             sub_records.append(tmp_rec)
        #             pbar.update() 
        with open(keep_correct_f, "w") as handle:
            SeqIO.write(sub_records, handle, file_type)  
        # del shared_objects, keep_bcool_name_lst, sub_records
        del keep_bcool_name_lst, sub_records
        return

    def select_correction(self, f_corrected_file, f_original, f_bcool, prefix, frequency_file):
        correct_non_isolates_seqs_lst = []
        correct_record_iterator, f_correct_type = parse_data(f_corrected_file)
        for rec in correct_record_iterator:
            correct_non_isolates_seqs_lst.append(str(rec.seq))

        ori_seqs_lst = []
        cor_seqs_lst = []
        ori_seq2name_dict = {}
        bcool_seq2name_dict = {}
        ori_record_iterator, ori_file_type = parse_data(f_original)
        bcool_record_iterator, bcool_file_type = parse_data(f_bcool)
        for ori_rec, cor_rec in zip(ori_record_iterator, bcool_record_iterator):
            ori_seq = str(ori_rec.seq)
            cor_seq = str(cor_rec.seq)
            ori_name = str(ori_rec.id)
            cor_name = str(cor_rec.id)
            ori_seqs_lst.append(ori_seq)
            cor_seqs_lst.append(cor_seq)
            
            ori_seq2name_dict.setdefault(ori_seq,[]).append(ori_name)
            bcool_seq2name_dict.setdefault(cor_seq,[]).append(cor_name)

        del ori_record_iterator, bcool_record_iterator, correct_record_iterator
        self.MM.measure()
        ori_seqs_set = set(ori_seqs_lst)
        cor_seqs_set = set(cor_seqs_lst)

        total_seqs_set = ori_seqs_set | set(correct_non_isolates_seqs_lst)
        del correct_non_isolates_seqs_lst

        self.logger.debug("No. of sequences(isolates) before bcool correction: {}".format(len(ori_seqs_set)))
        if len(ori_seqs_set) == len(ori_seqs_lst):
            self.logger.debug("Each the isolates' frequency equals to one before bcool correction.")
        else:
            self.logger.debug("Before bcool correction, not all the isolates' frequency equals to one.")
        
        self.logger.debug("No. of unique sequences (isolates) after bcool correction: {}".format(len(cor_seqs_set)))

        intersect = ori_seqs_set & cor_seqs_set
        self.logger.debug("No. of the sequences' intersection before and after bcool correction: {}".format(len(intersect)))

        union = cor_seqs_set | ori_seqs_set
        self.logger.debug("No. of the sequences' union for isolates before and after bcool correction: {}".format(len(union)))

        new_sequences = cor_seqs_set - ori_seqs_set
        self.logger.debug("No. of newly generated sequences for isolates before and after bcool correction: {}".format(len(new_sequences)))

        total_intersect = total_seqs_set & cor_seqs_set
        self.logger.debug("No. of the sequences' intersection between corrected isolates and corrected non-isolates and not corrected isolates: {}".format(len(total_intersect)))

        total_union = cor_seqs_set | total_seqs_set
        self.logger.debug("No. of the sequences' union for corrected isolates with corrected non-isolates and not corrected isolates: {}".format(len(total_union)))
        
        select_new_sequences = new_sequences & total_seqs_set
        self.logger.debug("No. of newly generated sequences for corrected isolates existing in corrected non-isolates or not corrected isolates: {}".format(len(select_new_sequences)))

        del total_seqs_set

        cor_seqs_lst_dict = Counter(cor_seqs_lst)
        
        seqs_not_1 = []
        for key, v in cor_seqs_lst_dict.items():
            if v != 1:
                seqs_not_1.append(key)
        seqs_not_1_set = set(seqs_not_1)

        intersect_sequences_set = seqs_not_1_set & total_intersect
        new_sequences_set = seqs_not_1_set & select_new_sequences

        del cor_seqs_lst_dict
        del union
        del new_sequences
        del intersect
        del ori_seqs_lst
        del ori_seqs_set
        del cor_seqs_lst
        del cor_seqs_set
        self.MM.measure()
        ############################################################################
        if intersect_sequences_set:
            original_records, _ = parse_data_index(f_original)
            keep_inter_bcool_name_lst = []

            for seq in intersect_sequences_set:
                bcool_names = bcool_seq2name_dict[seq]
                temp_bcool_name_lst = []
                for name in bcool_names:
                    ori_seq = str(original_records[name].seq)
                    if editdistance.eval(seq, ori_seq) <= 1:
                        temp_bcool_name_lst.append(name)                        
                if len(temp_bcool_name_lst) >= 2:
                    keep_inter_bcool_name_lst.extend(temp_bcool_name_lst)
                del temp_bcool_name_lst, bcool_names

            # save these frquency changing reads before and after correction
            len1 = len(keep_inter_bcool_name_lst)
            if len1 and self.iso_change_detail:
                original_bcool_inter_only_fastq_file = prefix + '_original_bcool_inter_only.' + ori_file_type
                extract_records(self.output_dir, keep_inter_bcool_name_lst, f_original, original_bcool_inter_only_fastq_file)

                bcool_inter_only_fastq_file = prefix + '_bcool_inter_only.' + bcool_file_type 
                extract_records(self.output_dir, keep_inter_bcool_name_lst, f_bcool, bcool_inter_only_fastq_file)
                with open(frequency_file, 'w') as f:
                    f.write('===========================================')
                    f.write('\n')
                    f.write('Intersection sequences: ')
                    f.write('\n')
                    f.write('===========================================')
                self.freqency_seq_extraction(original_bcool_inter_only_fastq_file, bcool_inter_only_fastq_file, frequency_file)
                os.system("rm %s" % original_bcool_inter_only_fastq_file)
                os.system("rm %s" % bcool_inter_only_fastq_file)
        self.MM.measure()
        ######################################################################
        if new_sequences_set:
            keep_new_seq_bcool_name_lst = []
            for new_seq in new_sequences_set:
                bcool_names = bcool_seq2name_dict[seq]
                temp_bcool_name_lst = []
                for name in bcool_names:
                    ori_seq = str(original_records[name].seq)
                    if editdistance.eval(new_seq, ori_seq) <= 3:
                        temp_bcool_name_lst.append(name)
                if len(temp_bcool_name_lst) >= 2:
                    keep_new_seq_bcool_name_lst.extend(temp_bcool_name_lst)
                del bcool_names, temp_bcool_name_lst

            len2 = len(keep_new_seq_bcool_name_lst)
            if len2 and self.iso_change_detail:
                # save these frquency changing reads before and after correction
                original_bcool_new_seq_only_fastq_file = prefix + '_original_bcool_new_seq_only.' + ori_file_type
                extract_records(self.output_dir, keep_new_seq_bcool_name_lst, f_original, original_bcool_new_seq_only_fastq_file)

                bcool_new_seq_only_fastq_file = prefix + '_bcool_new_seq_only.' + bcool_file_type
                extract_records(self.output_dir, keep_new_seq_bcool_name_lst, f_bcool, bcool_new_seq_only_fastq_file)
                
                with open(frequency_file, 'a') as f:
                    f.write('===========================================')
                    f.write('\n')
                    f.write('Newly generated sequences: ')
                    f.write('\n')
                    f.write('===========================================')
                self.freqency_seq_extraction(original_bcool_new_seq_only_fastq_file, bcool_new_seq_only_fastq_file, frequency_file)
                os.system("rm %s" % original_bcool_new_seq_only_fastq_file)
                os.system("rm %s" % bcool_new_seq_only_fastq_file)
        self.MM.measure()
        ##################################################################################################################
        corrected_isolates = self.output_dir + 'corrected_isolates.' + ori_file_type
        keep_correct_f = prefix + '_keep_corrected_records.' + ori_file_type 
        no_change_fastq_file = prefix + 'no_change_isolates.' + ori_file_type
        
        bcool_records, _ = parse_data_index(f_bcool)
        total_name_lst = list(bcool_records)

        if new_sequences_set and intersect_sequences_set:
            keep_bcool_name_lst = list(set(keep_new_seq_bcool_name_lst).union(set(keep_inter_bcool_name_lst)))
        elif new_sequences_set:
            keep_bcool_name_lst = keep_new_seq_bcool_name_lst
        elif intersect_sequences_set:
            keep_bcool_name_lst = keep_inter_bcool_name_lst
        else:
            keep_bcool_name_lst = []
        
        del intersect_sequences_set, new_sequences_set, bcool_records

        if len(keep_bcool_name_lst):
            if bcool_file_type != ori_file_type:
                # self.bcool_seqs(self.num_workers, keep_bcool_name_lst, f_bcool, f_original, keep_correct_f) 
                self.bcool_seqs(keep_bcool_name_lst, f_bcool, f_original, keep_correct_f)  
            else:    
                extract_records(self.output_dir, keep_bcool_name_lst, f_bcool, keep_correct_f)

            rest_name_lst = list(set(total_name_lst) - set(keep_bcool_name_lst))
            if len(rest_name_lst):
                extract_records(self.output_dir, rest_name_lst, f_original, no_change_fastq_file)

                self.logger.debug("rest_name_lst:{}".format(len(rest_name_lst)))  
        corrected_isolates = self.output_dir + 'corrected_isolates.' + ori_file_type
        if os.path.exists(keep_correct_f):
            os.system("cat %s %s > %s" % (keep_correct_f, no_change_fastq_file, corrected_isolates))
            os.system("rm %s" % keep_correct_f)
            os.system("rm %s" % no_change_fastq_file)          
        else:   
            return self.isolates
        self.MM.measure()
        self.MM.stop()
        return corrected_isolates
'''
    def covert_fq2fa(self, fin, fout):
        # records_dict, _ = parse_data_dict(fin)
        id2seq_desc = {}
        records, _ = parse_data_index(fin)
        name_lst = list(records)
        for name in name_lst:
            seq = str(records[name].seq)
            desc = str(records[name].description)
            id2seq_desc[name] = [seq, desc]
        del records
        fa_records = []
        name_len = len(name_lst)
        with WorkerPool(self.num_workers, shared_objects=id2seq_desc, start_method='fork') as pool:
            with tqdm(total=name_len, desc=self.logger.info("Convert fq to fa"), miniters=int(name_len/self.min_iters)) as pbar:   
                for tmp_rec in pool.imap(self.fq2fa, name_lst):
                    fa_records.append(tmp_rec)
                    pbar.update() 
        with open(fout, "w") as handle:
            SeqIO.write(fa_records, handle, 'fasta')
        del id2seq_desc, name_lst, fa_records, name_len
        return

    def fq2fa(self, shared_objects, name):
        return SeqRecord(shared_objects[name][0], id=name, description=shared_objects[name][1])
    def bcool_seq(self, shared_objects, name):
        bcool_records_dict, ori_records_dict = shared_objects
        return SeqRecord(bcool_records_dict[name].seq, id=ori_records_dict[name].id, description=ori_records_dict[name].description, letter_annotations=ori_records_dict[name].letter_annotations)
'''
