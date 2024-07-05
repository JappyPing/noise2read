import os
from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import copy
from noise2read.utils import *
import sys
from collections import Counter

class UMIReadErrorCorrection():
    """
    A class contains umi and read error correction related functions
    """
    def __init__(self, logger, config):
        """
        initialize the UMIReadErrorCorrection class

        Args:
            logger (class): customized logging
            config (class): parameters below setting using configparser
        """
        self.logger = logger
        self.config = config

    def umi_read_error_correction(self, umi_correct_data, read_correct_data):
        # load umi data
        umi_record_iterator, umi_file_type = parse_data(umi_correct_data)
        umis2id_dict = {}
        id2umi_dict = {}
        total_umis = []
        for item in umi_record_iterator:
            umi = str(item.seq)
            # ll = len(seq)
            # seq_lens_set.add(ll)
            total_umis.append(umi)
            umis2id_dict.setdefault(umi, []).append(str(item.id))  
            id2umi_dict.setdefault(str(item.id), []).append(umi)     
        # load read data 
        read_record_iterator, read_file_type = parse_data(read_correct_data)
        reads2id_dict = {}
        id2read_dict = {}
        total_reads = []
        for item in read_record_iterator:
            read = str(item.seq)
            total_reads.append(read)
            reads2id_dict.setdefault(read, []).append(str(item.id))
            id2read_dict[str(item.id)] = read 
        no_change_id_lst = []
        changed_id2read = {}
        for umi, ids in umis2id_dict.items():
            cur_reads_lst = []
            for id in ids:
                cur_reads_lst.append(id2read_dict[id])
            cur_read2count = Counter(cur_reads_lst)
            num = len(cur_read2count)
            if num == 1:
                no_change_id_lst.extend(ids)
            else:
                filtered_dict = {k: v for k, v in cur_read2count.items() if v > self.config.high_freq_thre}

                if filtered_dict:
                    sorted_dict = dict(sorted(filtered_dict.items(), key=lambda item: item[1], reverse=True))
                    # top_reads = [item[0] for item in sorted_dict[:self.config.top_count]]
                    sorted_items = list(sorted_dict.items())
                    top_reads = [item[0] for item in sorted_items[:self.config.top_count]]
                    for id in ids:
                        cur_read = id2read_dict[id]
                        if cur_read in top_reads:
                            no_change_id_lst.append(id)
                        else:
                            closest_read = None
                            dis = float('inf')
                            for top_read in top_reads:
                                cur_dis = editdistance.eval(cur_read, top_read)

                                if cur_dis < dis:
                                    dis = cur_dis
                                    closest_read = top_read
                                if dis <= self.config.max_dis:
                                    changed_id2read[id] = closest_read
                                else:
                                    no_change_id_lst.append(id)
                else:
                    no_change_id_lst.extend(ids)
        # write corrected results
        record_dict, ori_file_type = parse_data_dict(read_correct_data)
        corrected_read_records = []
        for cur_id in no_change_id_lst:
            if read_file_type == 'fastq' or read_file_type == 'fq' or read_file_type == 'fastq.gz' or read_file_type == 'fq.gz':
                tmp_rec = SeqRecord(Seq(id2read_dict[cur_id]), id=cur_id, description=record_dict[cur_id].description, letter_annotations=record_dict[cur_id].letter_annotations)  
            else:
                tmp_rec = SeqRecord(Seq(id2read_dict[cur_id]), id=cur_id, description=record_dict[cur_id].description)  
            corrected_read_records.append(tmp_rec)

        base = read_correct_data.split('/')[-1].split('_corrected')[0]

        if ".gz" in read_file_type:
            out_file_tye = read_file_type.split(".gz")[0]
        else:
            out_file_tye = read_file_type

        final_corrected_read_file = self.config.result_dir + base + ".corrected.reads." + out_file_tye
        with open(final_corrected_read_file, "w") as handle:
            SeqIO.write(corrected_read_records, handle, read_file_type)

        if os.path.exists(read_correct_data):
            os.system("rm %s" % read_correct_data)

        return final_corrected_read_file


                
