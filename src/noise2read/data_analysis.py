# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-01-30 09:35:18
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-08-29 10:14:41

from collections import Counter
import collections
import math
from Bio import SeqIO
import os
import editdistance
import xlsxwriter
from noise2read.utils import *
from numpy.linalg import norm
from numpy import array
from numpy.fft import fft,rfft
# import numpy as np
# import antropy as ant
import cmath
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from random import shuffle
from tqdm import tqdm
from mpire import WorkerPool

class DataAnalysis():
    """
    A class to evaluate an error correction method.
    """
    def __init__(self, logger, config):
        """
        initialize the DataAnalysis class

        Args:
            logger (class): customized logging
            config (class): parameters setting using configparser
        """
        self.logger = logger
        self.config = config
        self.logger.info("-------------------------------------------------------------")
        self.logger.info(f'Raw dataset: {config.input_file}')
        self.logger.info(f'Correct dataset: {config.correct_data}')
        self.logger.info(f'Ground truth dataset: {config.ground_truth_data}')
        if os.path.exists(self.config.result_dir):
            self.logger.info("Directory '% s' already exists" % self.config.result_dir)
            # for f in os.listdir(self.config.result_dir):
            #     os.remove(os.path.join(self.config.result_dir, f))
        else:
            os.makedirs(self.config.result_dir)
            self.logger.info("Directory '% s' created" % self.config.result_dir)
        bases = config.correct_data.split('/')[-1]
        self.prefix = bases.split('.' + parse_file_type(config.correct_data))[0]
        self.result_file = os.path.join(config.result_dir, self.prefix + '.xlsx')
        self.workbook_flie = xlsxwriter.Workbook(self.result_file)
        
    def percentage(self, part, whole):
        """
        percentage calculation

        Args:
            part (int): _description_
            whole (int): _description_

        Returns:
            str: percentage, if whole equals 0, return 'Not exist in raw'
        """
        if whole != 0:
            Percentage = 100 * float(part)/float(whole)
            return str(Percentage) + "%"
        else:
            return str('Not exist in raw')
        
    def ranking(self, num, num_lst):
        """
        get read count's ranking

        Args:
            num (int): read count
            num_lst (list): a list of read counts

        Returns:
            int: ranking
        """
        for n in num_lst:
            if num == n:
                return num_lst.index(n) + 1

    def evaluation(self):
        """
        Invoke the function to run the whole evaluation process, and write the result to files. 
        """
        # read the input using SeqIO
        raw_record_iterator, raw_file_tye = parse_data(self.config.input_file)
        correct_record_iterator, correct_file_tye = parse_data(self.config.correct_data)
        raw_seqs = []
        correct_seqs = []

        self.raw_len2seqs_dict = {}
        self.raw_len_lst = []
        total_reads_num = 0
        
        raw_record_dict = {}
        correct_record_dict = {}
        id_lst = []
        for raw_item, correct_item in zip(raw_record_iterator, correct_record_iterator):
            raw_id = str(raw_item.id)
            raw_seq = str(raw_item.seq)

            cor_id = str(correct_item.id)
            cor_seq = str(correct_item.seq)
            
            raw_seqs.append(raw_seq)
            correct_seqs.append(cor_seq)
            
            ll = len(raw_seq)
            self.raw_len_lst.append(ll)
            self.raw_len2seqs_dict.setdefault(ll, []).append(raw_seq)
            total_reads_num += 1
            
            raw_record_dict[raw_id] = raw_seq
            correct_record_dict[cor_id] = cor_seq
            
            id_lst.append(raw_id)
        del raw_record_iterator, correct_record_iterator
        corrected_reads_num = 0

        for item in tqdm(id_lst):
            ori_read = raw_record_dict[item]
            cor_read = correct_record_dict[item]
            if str(cor_read) != str(ori_read):
                corrected_reads_num += 1
                
        raw_read2count = Counter(raw_seqs)
        del raw_seqs
        correct_read2count = Counter(correct_seqs)
        del correct_seqs

        raw_read2count_val_lst = list(raw_read2count.values())
        raw_read2count_val_lst.sort(reverse=True)

        correct_read2count_val_lst = list(correct_read2count.values())
        correct_read2count_val_lst.sort(reverse=True)

        # raw_counts_lst = raw_read2count_val_lst[ : self.config.top_n * 10]
        correct_counts_lst = correct_read2count_val_lst[ : self.config.top_n * 10]
        
        # self.workbook_flie = xlsxwriter.Workbook(self.result_file)
        worksheet = self.workbook_flie.add_worksheet('Read Count')
        # Use the worksheet object to write
        # data via the write() method.
        worksheet.write('A1', 'Read')
        worksheet.write('B1', 'Count after Correction')
        worksheet.write('C1', 'Count before Correction')
        worksheet.write('D1', 'Change in percentage')
        worksheet.write('E1', 'Rank after Correction')
        worksheet.write('F1', 'Rank before Correction')

        row = 1
        col = 0
        for read, count in correct_read2count.most_common(self.config.top_n):
            worksheet.write(row, col, read)
            worksheet.write(row, col + 1, count)
            
            raw_read_count = raw_read2count[read]
            worksheet.write(row, col + 2, raw_read_count)
            abundance_percentage = self.percentage(count - raw_read_count, raw_read_count)
            worksheet.write(row, col + 3, abundance_percentage)

            correct_rank = self.ranking(count, correct_counts_lst)
            worksheet.write(row, col + 4, correct_rank)
            raw_rank = self.ranking(raw_read_count, raw_read2count_val_lst)
            worksheet.write(row, col + 5, raw_rank)
            row += 1

        self.logger.debug(f'Ground truth: {self.config.ground_truth_data}')
        if self.config.ground_truth_data:
            self.evaluation_with_groundtruth()
        self.logger.info("corrected {} reads out of {} ({:.6f}) reads".format(corrected_reads_num, total_reads_num, corrected_reads_num/total_reads_num))
        raw_unique_num = len(raw_read2count)
        correct_unique_num = len(correct_read2count)
        self.logger.info("Unique reads decreased by {} from {} to {} ({:.6f}) reads".format(raw_unique_num - correct_unique_num, raw_unique_num, correct_unique_num, (raw_unique_num - correct_unique_num)/raw_unique_num))
        # no ground truth entropy
        self.evaluation_without_groundtruth(raw_read2count, correct_read2count, total_reads_num)
        del raw_read2count, correct_read2count, total_reads_num
        return

    # rely on read frequency instead of sequecing id
    def evaluation_with_groundtruth(self):
        """
        Evaluating an error correction method using raw, true and corrected data based on read frequency instead of sequecing id
        """
        read_level = {'tp':0, 'fp':0, 'tn':0, 'fn':0}
        base_level = {'tp':0, 'fp':0, 'tn':0, 'fn':0}
        
        true_records, true_f_tye = parse_data_index(self.config.ground_truth_data)
        error_records, err_f_tye = parse_data_index(self.config.input_file)
        correct_records, cor_f_tye = parse_data_index(self.config.correct_data)

        total_reads_num = 0
        # for calculating entropy
        correct_errfree_seqs = []
        correct_err_seqs = []
        raw_errfreee_seqs = []
        raw_err_seqs = []

        true_seqs_lst = []
        raw_seqs_lst = []
        correct_seqs_lst = []
        true_umi2seqs = {}
        raw_umi2seqs = {}
        correct_umi2seqs = {}

        true_seq2umi = {}
        for name in error_records:
            true_seq = str(true_records[name].seq)
            raw_seq = str(error_records[name].seq)
            correct_seq = str(correct_records[name].seq) 

            umi_base = str(true_records[name].description).split('//')[0]
            umi = umi_base.split(':')[1]
            true_umi2seqs.setdefault(umi, []).append(true_seq)
            raw_umi2seqs.setdefault(umi, []).append(raw_seq)
            correct_umi2seqs.setdefault(umi, []).append(correct_seq)

            true_seq2umi.setdefault(true_seq, []).append(umi)

            true_seqs_lst.append(true_seq)     
            raw_seqs_lst.append(raw_seq)
            correct_seqs_lst.append(correct_seq)
        del true_records, correct_records, error_records
        raw_read2count = collections.Counter(raw_seqs_lst)
        # true_read2count = collections.Counter(true_seqs_lst)
        correct_read2count = collections.Counter(correct_seqs_lst)

        num = 0
        read_num = 0
        for seq in true_seq2umi:
            unique_umi = set(true_seq2umi[seq])
            if len(unique_umi) > 1:
                num += 1
                read_num += correct_read2count[seq]
                # print(unique_umi, )
        print(num, read_num)
############################################################################################
        # merge same sequences with different umis
        true_umis2seqs = {}
        raw_umis2seqs = {}
        correct_umis2seqs = {}
        for seq in true_seq2umi:
            unique_umi = list(set(true_seq2umi[seq]))
            umi_num = len(unique_umi)
            umi = tuple(unique_umi)
            if umi_num == 1:
                tmp_umi = unique_umi[0]
                true_umis2seqs[umi] = true_umi2seqs[tmp_umi]
                raw_umis2seqs[umi] = raw_umi2seqs[tmp_umi]
                correct_umis2seqs[umi] = correct_umi2seqs[tmp_umi]   
            else:
                for tmp_umi in unique_umi:
                    true_umis2seqs.setdefault(umi, []).extend(true_umi2seqs[tmp_umi])   
                    raw_umis2seqs.setdefault(umi, []).extend(raw_umi2seqs[tmp_umi])
                    correct_umis2seqs.setdefault(umi, []).extend(correct_umi2seqs[tmp_umi])      
############################################################################################
        fn_lst = []
        positive_fp_lst = []
        # ii = 0
        original_high_ambiguous = []

        for umi in true_umis2seqs:
            true_seqs = true_umis2seqs[umi]
            raw_seqs = raw_umis2seqs[umi]
            correct_seqs = correct_umis2seqs[umi]

            true_seqs_len = len(true_seqs)
            total_reads_num += true_seqs_len
            # raw base level
            total_positive_bases = 0
            total_negative_bases = 0
            # raw read level
            total_positive_reads = 0
            total_negative_reads = 0

            # correction base level
            cor_positive_bases = 0
            cor_negative_bases = 0
            # correction read level
            cor_positive_reads = 0
            cor_negative_reads = 0
            
            if len(set(true_seqs)) != 1:
                self.logger.exception("UMI(s) contain(s) more than one true sequences.")
            t_seq = true_seqs[0]
            r_seqs = list(set(raw_seqs))
            r_seq2counts = collections.Counter(raw_seqs)
            if len(r_seqs) > 1:
                for r_seq in r_seqs:
                    t_r_dis = editdistance.eval(t_seq, r_seq)
                    # if t_r_dis == 1 or t_r_dis == 2:
                    if t_r_dis == 1:
                        r_seq_count = raw_read2count[r_seq]
                        t_seq_count = raw_read2count[t_seq]
                        if r_seq_count > 5 and t_seq_count > 5:
                            # original_high_ambiguous.append([umi, t_r_dis, t_seq, t_seq_count, r_seq, r_seq_count, r_seq2counts[r_seq]])
                            line = [str(umi), t_r_dis, t_seq, t_seq_count, r_seq, r_seq_count, r_seq2counts[r_seq]]
                            if line not in original_high_ambiguous:
                                original_high_ambiguous.append(line)
                # ii += 1

            for i in range(true_seqs_len):
                raw_seq = raw_seqs[i]
                true_seq = true_seqs[i]
                seq_len = len(true_seq)
                true_raw_dis = editdistance.eval(true_seq, raw_seq)
                # raw base level
                total_positive_bases += true_raw_dis
                total_negative_bases += seq_len - true_raw_dis
                # raw read level
                if true_raw_dis == 0:
                    raw_errfreee_seqs.append(raw_seq)
                    total_negative_reads += 1
                else:
                    raw_err_seqs.append(raw_seq)
                    total_positive_reads += 1
                # after correction
                correct_seq = correct_seqs[i]
                true_cor_dis = editdistance.eval(true_seq, correct_seq)
                # correction base level
                cor_positive_bases += true_cor_dis
                cor_negative_bases += seq_len - true_cor_dis
                # correction read level
                if true_cor_dis == 0:
                    correct_errfree_seqs.append(correct_seq)
                    cor_negative_reads += 1
                else: 
                    correct_err_seqs.append(correct_seq)
                    cor_positive_reads += 1
            # note: for any postive read, if it is modified correctly, tp will increase one. if it is modified incorrectly, it will keep as a fn rather than a fp. Similarly, for any negative read, 
            # base level
            if total_positive_bases >= cor_positive_bases:
                base_level['tp'] += total_positive_bases - cor_positive_bases
                base_level['fn'] += cor_positive_bases  
            if total_negative_bases >= cor_negative_bases:  
                base_level['tn'] += cor_negative_bases
                base_level['fp'] += total_negative_bases - cor_negative_bases

            if total_positive_bases < cor_positive_bases:
                base_level['tp'] += 0
                base_level['fn'] += total_positive_bases
                # self.logger.warning(f'Base-level Evaluation: introduced additional {cor_positive_bases - total_positive_bases} positive bases after correction for UMI {str(umi)}')  
            if total_negative_bases < cor_negative_bases:
                base_level['tn'] += total_negative_bases
                base_level['fp'] += 0    
                # self.logger.warning(f'Base-level Evaluation: introduced additional {cor_negative_bases - total_negative_bases} negative bases after correction for UMI {str(umi)}')   

            # read level
            if total_positive_reads >= cor_positive_reads:
                read_level['tp'] += total_positive_reads - cor_positive_reads
                read_level['fn'] += cor_positive_reads  
            if total_negative_reads >= cor_negative_reads:
                read_level['tn'] += cor_negative_reads
                read_level['fp'] += total_negative_reads - cor_negative_reads  
            if total_positive_reads < cor_positive_reads:
                read_level['tp'] += 0
                read_level['fn'] += total_positive_reads
                # self.logger.warning(f'Read-level Evaluation: introduced additional {cor_positive_reads - total_positive_reads} positive reads after correction for UMI {str(umi)}') 
            if total_negative_reads < cor_negative_reads:
                read_level['tn'] += total_negative_reads
                read_level['fp'] += 0
                # self.logger.warning(f'Read-level Evaluation: introduced additional {cor_negative_reads - total_negative_reads} negative reads after correction for UMI {str(umi)}') 

            ## fn
            raw_positive_seqs = set(raw_seqs) - set(true_seqs)
            untouched_positive_seqs = set(correct_seqs) - set(true_seqs)

            if len(raw_positive_seqs) > 0:
                if len(untouched_positive_seqs) > 0:
                    # true_group = [t_seq, true_read2count[t_seq], raw_read2count[t_seq], correct_read2count[t_seq]]
                    for seq in untouched_positive_seqs:
                        fn_lst.append([str(umi), t_seq, raw_read2count[t_seq], seq, raw_read2count[seq], correct_read2count[t_seq], correct_read2count[seq]]) 
        del raw_read2count, correct_read2count                                        
        self.evaluation_metircs(base_level, 'Base level')
        self.evaluation_metircs(read_level, 'Read level') 
        rawset_entropy = self.set_entropy(len(raw_err_seqs), len(raw_errfreee_seqs), total_reads_num)
        correctset_entropy = self.set_entropy(len(correct_err_seqs), len(correct_errfree_seqs), total_reads_num)
        self.save_entropy('Purity entropy', rawset_entropy, correctset_entropy)
        del correct_err_seqs, correct_errfree_seqs, raw_errfreee_seqs, raw_err_seqs
        return

    def evaluation_without_groundtruth(self, raw_read2count, correct_read2count, total_reads_num):
        """
        Evaluating an error correction method using raw and corrected data (before and after) base on defined read counts entropy

        Args:
            raw_read2count (dict): A dictionary to save reads (keys) and its counts (values) for raw dataset
            correct_read2count (dict): A dictionary to save reads (keys) and its counts (values) for corrected dataset of raw dataset
            total_reads_num (int): the number of the total reads
        """
        # rawset_entropy_noTruth = self.read_counts_entropy(raw_read2count, total_reads_num)
        # correctset_entropy_noTruth = self.read_counts_entropy(correct_read2count, total_reads_num)
        # self.save_entropy('ReadCountEntropy', rawset_entropy_noTruth, correctset_entropy_noTruth) 

        # entropy = self.noise2signal_entropy(raw_read2count, correct_read2count, total_reads_num)
        entropy = self.noise2read_entropy(raw_read2count, correct_read2count, total_reads_num)
        # self.save_entropy('Entropy H', entropy[0], entropy[1]) 
        self.logger.info("{}: raw dataset entropy: {}, correct dataset entropy: {}".format('Entropy H', entropy[0], entropy[1]))
        self.logger.info("Information Gain: {}".format(entropy[0] - entropy[1]))
        
        worksheet3 = self.workbook_flie.add_worksheet('Non-frequent Entropy')
        worksheet3.write('A1', 'H')
        worksheet3.write('A2', entropy[0])
        worksheet3.write('B1', "H'")
        worksheet3.write('B2', entropy[1]) 
        worksheet3.write('C1', '\u0394 H')
        worksheet3.write('C2', entropy[0] - entropy[1]) 

        self.workbook_flie.close()
        # total_variation = self.total_variation_differnce(raw_read2count, correct_read2count, total_reads_num)
        # print("Total variation distance: {}".format(total_variation))
        # relative_entropyq2p, entropyp2q = self.total_variation_differnce(raw_read2count, correct_read2count, total_reads_num)
        # entropy = self.spectral_entropy(raw_read2count, correct_read2count, total_reads_num)
        # raw_entropy = self.spectral_entropy(raw_read2count, raw_read2count, total_reads_num)
        
        return
        
    def evaluation_metircs(self, confusion_dict, sheet_name):
        """
        Evaluation metrics

        Args:
            confusion_dict (dict): A dictionary consists the values of TP, TN, FP and FN
            sheet_name (str): A string to indicate a sheet name of a csv file
        """
        tp = confusion_dict['tp']
        tn = confusion_dict['tn']
        fp = confusion_dict['fp']
        fn = confusion_dict['fn']
        if tp + tn + fp + fn != 0: 
            accuracy = (tp + tn) / (tp + tn + fp + fn)
        else:
            accuracy = "None"
        if tp + fp != 0:
            precision = tp / (tp + fp)
        else:
            precision = "None"
        if tp + fn != 0:
            recall = tp / (tp + fn)
            gain = (tp - fp) / (tp + fn)
        else:
            recall = "None"
            gain = "None"
        if fp + tn != 0:
            fall_out = fp / (fp + tn)
        else:
            fall_out = "None"
        # self.logger.info(sheet_name)
        self.logger.info("{}: Accuracy: {}, Precision: {}, Recall: {}, Positive Gain: {}, Fall-out: {}".format(sheet_name, accuracy, precision, recall, gain, fall_out))
        worksheet3 = self.workbook_flie.add_worksheet(sheet_name)
        worksheet3.write('A1', 'TP')
        worksheet3.write('A2', tp)
        worksheet3.write('B1', 'FP')
        worksheet3.write('B2', fp)
        worksheet3.write('C1', 'FN')
        worksheet3.write('C2', fn)
        worksheet3.write('D1', 'TN')
        worksheet3.write('D2', tn)

        worksheet3.write('E1', 'Accuracy')
        worksheet3.write('E2', accuracy)
        worksheet3.write('F1', 'Precision')
        worksheet3.write('F2', precision)
        worksheet3.write('G1', 'Recall')
        worksheet3.write('G2', recall)
        worksheet3.write('H1', 'Positive Gain')
        worksheet3.write('H2', gain)
        worksheet3.write('I1', 'Fall-out')
        worksheet3.write('I2', fall_out)
        return

    def entropy(self, seqs_num, total_num):
        """
        part of the dataset entropy to measure the reads impurity of error-contained and error-free
        Args:
            seqs_num (int): the number of error-contained or error-free reads
            total_num (int): the number of total reads

        Returns:
            float: part of the dataset entropy
        """
        p = seqs_num / total_num
        if p == 1 or p == 0:
            entropy_val = 0
        else:
            entropy_val = -(p * math.log2(p))
        return entropy_val

    def gain2heatmap(self, correct2raw_diff):
        
        x_len = len(correct2raw_diff)
        x_len_sqrt = round(math.sqrt(x_len))
        xx = x_len_sqrt * x_len_sqrt
        val = x_len - xx
        if xx == x_len:
            m = x_len_sqrt
            n = x_len_sqrt
            x = np.array((correct2raw_diff))
            x_res=x.reshape(m, n)
 
        elif val < x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 1
            remain = (m - val) * [0]
            x = np.array((correct2raw_diff + remain))
            x_res = x.reshape(m, n)

        elif val > x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 2
            remain = (val - m) * [0]
            x = np.array((correct2raw_diff + remain))
            x_res=x.reshape(m, n)

        # thre_max = 0
        # thre_max = 1000000
        fig, ax = plt.subplots(figsize=(6, 6))
        # cmap = mpl.colors.ListedColormap(['whitesmoke', 'cornflowerblue', 'blueviolet', 'fuchsia', 'yellow', 'lawngreen', 'orange', 'darkorange', 'black', 'forestgreen', 'deepskyblue', 'yellow', 'orange'])
        # cmap = mpl.colors.ListedColormap([ '#8dd3c7', '#ffffb3', '#bebada', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 'deepskyblue', '#bc80bd', '#ccebc5', '#ffed6f', 'whitesmoke'])
        # name_lst = []
        # i = 0
        # for name, _ in mpl.colors.cnames.items():
        #     if name not in ['black', 'maroon', 'darked', 'darkgreen', 'darkslategray', 'darkslategrey', 'midnightblue', 'navy', 'darkblue', 'mediumblue']:
        #         name_lst.append(name)
        #     if i == 149:
        #         break
        #     i += 1
        # print(name_lst)
        # name_lst = ['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'mediumaquamarine', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple', 'rebeccapurple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey',  'whitesmoke', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#bc80bd', '#ccebc5', '#ffed6f']
        # cmap = ListedColormap(name_lst)
        # cmap = mpl.colors.ListedColormap(['whitesmoke', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3'])                                  
        # cmap.set_over('yellow')
        # cmap.set_under('red')
        # masked_array = np.ma.array (x_res, mask=np.isnan(x_res))
        # cur_cmap = mpl.cm.get_cmap(cmap)
        cur_cmap = mpl.cm.get_cmap('tab20c')#
        # cur_cmap.set_over('red')
        cur_cmap.set_under('red')
        cur_cmap.set_bad(color='red')
        image = ax.imshow(x_res, cmap=cur_cmap) #interpolation='none', , vmax=thre_max
        # image = ax.pcolor(x_res, cmap=cmap, vmin=thre_min, antialiased=True)
        cbar = fig.colorbar(image, extend='min', shrink=0.8)
        cbar.cmap.set_under('red')
        # cbar.cmap.set_over('red')
        fig.tight_layout()
        fig.savefig(os.path.join(self.config.result_dir, self.prefix + '_information_gain.png'), transparent=True)

        return

    def entropy_item(self, total_freq, freq):
        if freq > 0:
            p = freq / total_freq
            result = - p * math.log2(p)
        return result

    def noise2read_entropy(self, raw_read2count, correct_read2count, total_num):
        raw_unique_reads = set(raw_read2count.keys())
        frequent_reads = set([k for k, v in raw_read2count.items() if v > self.config.high_freq_thre])

        # raw dataset
        non_frequent_raw_reads = raw_unique_reads - frequent_reads
        raw_entropy_items = []
        for read in non_frequent_raw_reads:
            raw_entropy_items.append(raw_read2count[read])
        del non_frequent_raw_reads
        
        # correct dateset
        correct_unique_reads = set(correct_read2count.keys())
        non_frequent_correct_reads = correct_unique_reads - frequent_reads
        correct_entropy_items = []
        for read in non_frequent_correct_reads:
            correct_entropy_items.append(correct_read2count[read])
        del non_frequent_correct_reads

        # new reads
        new_reads = correct_unique_reads - raw_unique_reads
        new_reads_num = len(new_reads)
        self.logger.info("Wrongly introduced {} new reads".format(new_reads_num))
        del new_reads
        
        # information gain
        raw_kept_counts = []
        correct_kept_counts = []
        kept_reads = correct_unique_reads & raw_unique_reads
        for read in kept_reads:
            correct_kept_counts.append(correct_read2count[read])
            raw_kept_counts.append(raw_read2count[read])
        del kept_reads, correct_read2count

        raw_removed_reads = raw_unique_reads - correct_unique_reads
        raw_removed_items = []
        for read in raw_removed_reads:
            raw_removed_items.append(raw_read2count[read])

        del raw_removed_reads, raw_unique_reads, correct_unique_reads, raw_read2count 
        ###################################################################
        # raw entropy
        raw_nonFre_reads_total_num = sum(raw_entropy_items)
        try:
            with WorkerPool(self.config.num_workers, shared_objects=raw_nonFre_reads_total_num, start_method='fork') as pool:
                raw_entropy_lst = pool.map(self.entropy_item, raw_entropy_items)
            raw_entropy = sum(raw_entropy_lst) 
        except KeyboardInterrupt:
            # Handle termination signal (Ctrl+C)
            pool.terminate()  # Terminate the WorkerPool before exiting
        except Exception:
            # Handle other exceptions
            pool.terminate()  # Terminate the WorkerPool before exiting
            raise
        del raw_entropy_items, raw_entropy_lst, raw_nonFre_reads_total_num
        # correct entropy
        correct_nonFre_reads_total_num = sum(correct_entropy_items)
        try:
            with WorkerPool(self.config.num_workers, shared_objects=correct_nonFre_reads_total_num, start_method='fork') as pool:
                correct_entropy_lst = pool.map(self.entropy_item, correct_entropy_items) 
            correct_entropy = sum(correct_entropy_lst)
        except KeyboardInterrupt:
            # Handle termination signal (Ctrl+C)
            pool.terminate()  # Terminate the WorkerPool before exiting
        except Exception:
            # Handle other exceptions
            pool.terminate()  # Terminate the WorkerPool before exiting
            raise
        del correct_entropy_items, correct_entropy_lst, correct_nonFre_reads_total_num
        ##################################################################################
        #information gain (\delta I) heatmap
        try:
            with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                raw_kept_entropy_lst = pool.map(self.entropy_item, raw_kept_counts)
        except KeyboardInterrupt:
            # Handle termination signal (Ctrl+C)
            pool.terminate()  # Terminate the WorkerPool before exiting
        except Exception:
            # Handle other exceptions
            pool.terminate()  # Terminate the WorkerPool before exiting
            raise
        del raw_kept_counts
        ########
        try:
            with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                correct_kept_entropy_lst = pool.map(self.entropy_item, correct_kept_counts) 
        except KeyboardInterrupt:
            # Handle termination signal (Ctrl+C)
            pool.terminate()  # Terminate the WorkerPool before exiting
        except Exception:
            # Handle other exceptions
            pool.terminate()  # Terminate the WorkerPool before exiting
            raise
        del correct_kept_counts
        ######
        try:
            with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                raw_removed_entropy_items_lst = pool.map(self.entropy_item, raw_removed_items)
        except KeyboardInterrupt:
            # Handle termination signal (Ctrl+C)
            pool.terminate()  # Terminate the WorkerPool before exiting
        except Exception:
            # Handle other exceptions
            pool.terminate()  # Terminate the WorkerPool before exiting
            raise
        del raw_removed_items
        #######
        # for i in raw_removed_entropy_items_lst:
        #     if i <=0:
        #         print('Warning')
        entropy_item_lst = []
        for i, j in zip(raw_kept_entropy_lst, correct_kept_entropy_lst):
            entropy_item_lst.append(i - j)
        entropy_item_lst.extend(raw_removed_entropy_items_lst)
        if new_reads_num > 0:
            entropy_item_lst.extend([np.nan] * new_reads_num)
        shuffle(entropy_item_lst)
        self.gain2heatmap(entropy_item_lst)
        return [raw_entropy, correct_entropy]

    def set_entropy(self, err_seqs_num, errfree_seqs_num, total_num):
        """
        dataset entropy to measure the reads impurity of error-contained and error-free

        Args:
            err_seqs_num (int): the number of error-conatined reads
            errfree_seqs_num (int): the number of error-free reads
            total_num (int): the number of total reads

        Returns:
            float: dataset entropy
        """
        # total_seqs = err_seqs + errfree_seqs
        # unique_len = len(set(total_seqs))
        err_entropy = self.entropy(err_seqs_num, total_num)
        errfree_entropy = self.entropy(errfree_seqs_num, total_num)
        return errfree_entropy + err_entropy

    def save_entropy(self, sheet_name, raw_set_entropy, correct_entropy):
        """
        save entropy result to a csv sheet

        Args:
            sheet_name (str): A string to indicate a sheet name of a csv file
            raw_set_entropy (float): raw dataset entropy
            correct_entropy (float): corrected dataset entropy
        """
        self.logger.info("{}: raw dataset entropy: {}, correct dataset entropy: {}".format(sheet_name, raw_set_entropy, correct_entropy))
        worksheet3 = self.workbook_flie.add_worksheet(sheet_name)
        worksheet3.write('A1', 'E')
        worksheet3.write('A2', raw_set_entropy)
        worksheet3.write('B1', "E'")
        worksheet3.write('B2', correct_entropy) 
        worksheet3.write('C1', '\u0394 E')
        worksheet3.write('C2', raw_set_entropy - correct_entropy) 
        return   

    #############################################################################################################
    # Warning: the following functions have been deprecated and may contain bugs. If you want to use, you must check every row carefully.

    '''
    def noise2signal_entropy(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = set(list(raw_read2count.keys()))
        correct_unique_reads = set(list(correct_read2count.keys()))
        total_unique_reads = raw_unique_reads.union(correct_unique_reads)
        id2read = {}
        idx = 0
        for read in total_unique_reads:
            id2read[idx] = read
            idx += 1
    
        raw_only_reads = raw_unique_reads - correct_unique_reads
        correct_only_reads = correct_unique_reads - raw_unique_reads
        raw_correct_inter = correct_unique_reads & raw_unique_reads

        for i in len(total_unique_reads):
            read = id2read[i]
            if read in 

        entropy = 0
        n = len(raw_unique_reads)
        correct_in_counts = []
        out_counts = []        
        for read in correct_unique_reads:
            correct_count = correct_read2count[read]
            # p = correct_count / total_num  
            if read in raw_unique_reads:
                correct_in_counts.append(correct_count)
            else:
                out_counts.append(correct_count)

        correct_in_counts.sort(reverse=True)
        out_counts.sort(reverse=True)

        shared_vars = total_num, n, 1
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            in_entropy = pool.map(self.entropy_item, in_counts)
        if len(out_counts) > 0:
            shared_vars = total_num, n, self.config.delta
            with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
                out_entropy = pool.map(self.entropy_item, out_counts)
            entropy = sum(in_entropy) + sum(out_entropy)
            shared_vars = total_num, n, 1
            with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
                out_entropy2 = pool.map(self.entropy_item, out_counts)
        else:
            entropy = sum(in_entropy)

        # for read in correct_unique_reads:
        #     correct_count = correct_read2count[read]
        #     p = correct_count / total_num  
        #     if read in raw_unique_reads:
        #         entropy += -(p * math.log2(p) / math.log2(n))
        #     else:
        #         entropy += -10*(p * math.log2(p) / math.log2(n))
        # print(entropy)
        # visualisation
        if raw_read2count == correct_read2count:
            if len(out_counts) > 0:
                raw_lst = in_entropy + len(out_entropy2) * [0]
            else:
                raw_lst = in_entropy
            raw_heatmap = os.path.join(self.config.result_dir, 'raw_entropy_heatmap.png')
            self.entropy_visualisation(raw_lst, )
        else:
            correct_lst = in_entropy + out_entropy2
            correct_heatmap = os.path.join(self.config.result_dir, 'correct_entropy_heatmap.png')
            self.entropy_visualisation(correct_lst, correct_heatmap)
        return entropy
    '''

    '''
    def evaluation_with_groundtruth0(self):
        read_level = {'tp':0, 'fp':0, 'tn':0, 'fn':0}
        base_level = {'tp':0, 'fp':0, 'tn':0, 'fn':0}
        true_records = SeqIO.index(self.config.ground_truth_data, parse_file_type(self.config.ground_truth_data))
        error_records = SeqIO.index(self.config.input_file, parse_file_type(self.config.input_file))
        correct_records = SeqIO.index(self.config.correct_data, parse_file_type(self.config.input_file))
        total_reads_num = 0
        # for calculating entropy
        correct_errfree_seqs = []
        correct_err_seqs = []
        raw_errfreee_seqs = []
        raw_err_seqs = []

        true_seqs_lst = []
        raw_seqs_lst = []
        correct_seqs_lst = []
        for name in error_records:
            true_seq = str(true_records[name].seq)
            raw_seq = str(error_records[name].seq)
            correct_seq = str(correct_records[name].seq)   
            true_seqs_lst.append(true_seq)     
            raw_seqs_lst.append(raw_seq)
            correct_seqs_lst.append(correct_seq)
        raw_read2count = collections.Counter(raw_seqs_lst)
        true_read2count = collections.Counter(true_seqs_lst)
        correct_read2count = collections.Counter(correct_seqs_lst)

        fn_lst = []

        for name in error_records:
            total_reads_num += 1
            true_seq = str(true_records[name].seq)
            raw_seq = str(error_records[name].seq)
            correct_seq = str(correct_records[name].seq)

            true_raw_dis = editdistance.eval(true_seq, raw_seq)
            true_cor_dis = editdistance.eval(true_seq, correct_seq)
            raw_cor_dis = editdistance.eval(raw_seq, correct_seq)
            #
            if true_cor_dis == 0:
                correct_errfree_seqs.append(correct_seq)
            else:
                correct_err_seqs.append(correct_seq)

            # actual condition with no errors
            # negative conditions
            if true_raw_dis == 0:
                raw_errfreee_seqs.append(raw_seq)
                if true_cor_dis == 0:
                    base_level['tn'] += len(true_seq)
                    read_level['tn'] += 1
                else:
                    base_level['tn'] += len(true_seq) - true_cor_dis
                    base_level['fp'] += true_cor_dis
                    read_level['fp'] += 1
            # true_raw_dis != 0 which means actual condition contain errors 
            else:
                raw_err_seqs.append(raw_seq)
                if true_cor_dis == 0:
                    base_level['tp'] += true_raw_dis
                    base_level['tn'] += len(true_seq) - true_raw_dis
                    read_level['tp'] += 1
                else:
                    read_level['fn'] += 1
                    fn_lst.append([raw_seq, raw_read2count[raw_seq], true_seq, true_read2count[true_seq], correct_seq, correct_read2count[correct_seq]])
                    if raw_cor_dis == 0:
                        base_level['tn'] += len(true_seq) - true_raw_dis
                        base_level['fn'] += true_raw_dis # true_raw_dis == true_cor_dis
                    else:
                        # measure the distance change with ture sequence before and after correction
                        gap = true_raw_dis - true_cor_dis
                        if gap > 0: # distance decrease after correction
                            base_level['tp'] += gap
                            base_level['tn'] += len(true_seq) - true_raw_dis
                            base_level['fn'] += true_cor_dis
                        elif gap < 0:
                            base_level['tn'] += len(true_seq) - true_cor_dis
                            base_level['fp'] += abs(gap)
                            base_level['fn'] += true_raw_dis
                        else:
                            base_level['tn'] += len(true_seq) - true_cor_dis
                            base_level['fn'] += true_cor_dis
        worksheet = self.workbook_flie.add_worksheet('fn')
        for row_num, data in enumerate(fn_lst):
            worksheet.write_row(row_num, 0, data)
        self.evaluation_metircs(base_level, 'Base level')
        self.evaluation_metircs(read_level, 'Read level') 
        rawset_entropy = self.set_entropy(raw_err_seqs, raw_errfreee_seqs, total_reads_num)
        correctset_entropy = self.set_entropy(correct_err_seqs, correct_errfree_seqs, total_reads_num)
        self.save_entropy('Dataset Entropy', rawset_entropy, correctset_entropy)
        # self.workbook_flie.close()

        return
    '''        
    '''
    # def entropy_visualisation(self, x_lst, out_heatmap, v_min, v_max):
    def entropy_visualisation(self, raw_lst, correct_lst):
        val_unique = list(set(raw_lst + correct_lst))

        norm_raw_lst = []
        for raw_v in raw_lst:
            if raw_v == 0:
                cur_bin_num = 0
            else:
                inter_val = int(raw_v / 5)
                cur_bin_num = 0
                for i in range(inter_val):
                    cur_bin_num += i+1
                    if cur_bin_num > inter_val:
                        cur_bin_num -= 1
                        break
            norm_raw_lst.append(cur_bin_num)

        norm_correct_lst = []
        for correct_v in correct_lst:
            if correct_v == 0:
                cur_bin_num = 0
            else:
                inter_val = int(correct_v / 5)
                cur_bin_num = 0
                for i in range(inter_val):
                    cur_bin_num += i+1
                    if cur_bin_num > inter_val:
                        cur_bin_num -= 1
                        break
            norm_correct_lst.append(cur_bin_num)

        # color_num = len(norm_raw_lst)
        
        x_len = len(norm_raw_lst)
        x_len_sqrt = round(math.sqrt(x_len))
        xx = x_len_sqrt * x_len_sqrt
        val = x_len - xx
        if xx == x_len:
            m = x_len_sqrt
            n = x_len_sqrt
            raw_x = np.array((norm_raw_lst))
            raw_x_res=raw_x.reshape(m, n)

            correct_x = np.array((norm_correct_lst))
            correct_x_res=correct_x.reshape(m, n)
            
        elif val < x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 1
            remain = (m - val) * [0]
            raw_x = np.array((norm_raw_lst + remain))
            raw_x_res = raw_x.reshape(m, n)

            correct_x = np.array((norm_correct_lst + remain))
            correct_x_res = correct_x.reshape(m, n)

        elif val > x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 2
            remain = (val - m) * [0]
            raw_x = np.array((norm_raw_lst + remain))
            raw_x_res=raw_x.reshape(m, n)

            correct_x = np.array((norm_correct_lst + remain))
            correct_x_res=correct_x.reshape(m, n)

        # x_len = len(correct_lst)
        # x_len_sqrt = round(math.sqrt(x_len))
        # xx = x_len_sqrt * x_len_sqrt
        # val = x_len - xx
        # if xx == x_len:
        #     m = x_len_sqrt
        #     n = x_len_sqrt
        #     correct_x = np.array((correct_lst))
        #     correct_x_res=correct_x.reshape(m, n)
            
        # elif val < x_len_sqrt:
        #     m = x_len_sqrt
        #     n = x_len_sqrt + 1
        #     remain = (m - val) * [0]
        #     correct_x = np.array((correct_lst + remain))
        #     correct_x_res = correct_x.reshape(m, n)

        # elif val > x_len_sqrt:
        #     m = x_len_sqrt
        #     n = x_len_sqrt + 2
        #     remain = (val - m) * [0]
        #     correct_x = np.array((correct_lst + remain))
        #     correct_x_res=correct_x.reshape(m, n)


        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(32,16)) #

        fig.subplots_adjust(wspace=0.01)
        # unique_val = set(x_lst)
        # color_num = len(unique_val)
        # print(color_num)
        # cmap = sns.color_palette(n_colors = color_num, linecolor='white')#"deep",
        # cmap = sns.color_palette("Paired")
        # cmap = LinearSegmentedColormap.from_list(N=256)

        # color_num = len(set(norm_raw_lst))
        max_v1 = max(norm_raw_lst)
        max_v2 = max(norm_correct_lst)
        color_num = max(max_v1, max_v2) + 1
        name_lst = []
        i = 0
        for name, _ in mpl.colors.cnames.items():
            name_lst.append(name)
            if i == color_num-1:
                break
        cmap = ListedColormap(name_lst)

        v_min = min(val_unique)
        v_max = max(val_unique)
        sns.heatmap(raw_x_res, square=True, ax=ax1, vmin=v_min, vmax=v_max)
        # sns.heatmap(raw_x_res, square=True, ax=ax1, cmap=cmap, cbar=True)
        # fig.colorbar(ax1.collections[0], ax=ax1, location="left", use_gridspec=False, pad=0.2) #, pad=0.2 ticks=val_unique, 

        # color_num2 = len(set(correct_lst))
        # name_lst2 = []
        # i = 0
        # for name2, _ in mpl.colors.cnames.items():
        #     name_lst2.append(name2)
        #     if i == color_num2-1:
        #         break
        # cmap2 = ListedColormap(name_lst2)
        # sns.heatmap(correct_x_res, square=True, cmap=cmap, ax=ax2, cbar=True)
        sns.heatmap(correct_x_res, square=True, ax=ax2, vmin=v_min, vmax=v_max)
        # fig.colorbar(ax2.collections[0], ax=ax2, location="right", use_gridspec=False, pad=0.2) #, pad=0.2  ticks=val_unique,
        # colorbar = ax1.collections[0].colorbar
        # colorbar.set_ticks(val_unique)
        # colorbar.set_ticklabels(val_unique)

        # plt.tight_layout()
        plt.savefig(os.path.join(self.config.result_dir, 'entropy_heatmap.png')) 
        plt.close()
        return
    def noise2signal_entropy0(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = list(raw_read2count.keys())
        correct_unique_reads = list(correct_read2count.keys())

        new_reads = list(set(correct_unique_reads) - set(raw_unique_reads))
        raw_unique_reads.sort()
        new_reads.sort()

        raw_in_counts = []
        correct_in_counts = []
        correct_out_counts = []
        for read in raw_unique_reads:
            raw_in_counts.append(raw_read2count[read])
            if read in correct_unique_reads:
                correct_in_counts.append(correct_read2count[read])
            else:
                correct_in_counts.append(0) # noise of identical reads removed by marked as 0
        if len(new_reads) > 0:
            for read in new_reads:
                correct_out_counts.append(correct_read2count[read])
            raw_out_counts = [0] * len(new_reads) # noise of new sequence generated by marked as 0

        # raw entropy
        n = len(raw_unique_reads)
        shared_vars = total_num, n, 1
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            raw_entropy_lst = pool.map(self.entropy_item, raw_in_counts)
        # correct entropy
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            correct_in_entropy_lst = pool.map(self.entropy_item, correct_in_counts) 

        raw_entropy = sum(raw_entropy_lst)
        new_seq_len = len(correct_out_counts) 
        if new_seq_len > 0:
            shared_vars = total_num, n, self.config.delta
            with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
                correct_out_entropy_lst = pool.map(self.entropy_item, correct_out_counts)

            correct_entropy = sum(correct_out_entropy_lst)
            for v in correct_in_entropy_lst:
                if v != -2:
                    correct_entropy += v
        else:
            correct_entropy = 0
            for v in correct_in_entropy_lst:
                if v != -2:
                    correct_entropy += v

        if new_seq_len > 0:
            # raw
            # raw_lst = copy.deepcopy(raw_entropy_lst + raw_out_counts)
            raw_lst = raw_in_counts + raw_out_counts
            # raw_lst = raw_in_counts
            # correct_lst = copy.deepcopy(correct_in_entropy_lst + correct_out_entropy_lst)
            correct_lst = correct_in_counts + correct_out_counts
        else:
            # raw_lst = copy.deepcopy(raw_entropy_lst)
            # correct_lst = copy.deepcopy(correct_in_entropy_lst)
            raw_lst = raw_in_counts
            correct_lst = correct_in_counts   
        self.entropy_visualisation(raw_lst, correct_lst)
        # raw_heatmap = os.path.join(self.config.result_dir, 'raw_entropy_heatmap.png')

        # self.entropy_visualisation(raw_lst, raw_heatmap, len(set(raw_lst)))
        # # correct
        # correct_heatmap = os.path.join(self.config.result_dir, 'correct_entropy_heatmap.png')
        # # self.entropy_visualisation(correct_lst, correct_heatmap,v_min, v_max)
        # # self.entropy_visualisation(correct_lst, correct_heatmap, color_num)
        # self.entropy_visualisation(correct_lst, correct_heatmap, len(set(correct_lst)))
        return [raw_entropy, correct_entropy]
    '''
    '''
    def noise2read_entropy(self, raw_read2count, correct_read2count):

        raw_unique_reads = set(raw_read2count.keys())
        correct_unique_reads = set(correct_read2count.keys())

        frequent_reads = set([k for k, v in raw_read2count.items() if v > self.config.high_freq_thre])

        # raw dataset
        non_frequent_raw_reads = raw_unique_reads - frequent_reads
        raw_entropy_items = []
        
        for read in non_frequent_raw_reads:
            raw_entropy_items.append(raw_read2count[read])

        raw_nonFre_reads_num = sum(raw_entropy_items)
        # raw entropy
        n = len(non_frequent_raw_reads)
        shared_vars1 = raw_nonFre_reads_num, n, 1
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars1, start_method='fork') as pool:
            raw_entropy_lst = pool.map(self.entropy_item, raw_entropy_items)
        raw_entropy = sum(raw_entropy_lst) 

        # correct dateset
        non_frequent_correct_reads = correct_unique_reads - frequent_reads

        correct_entropy_items = []

        for read in non_frequent_correct_reads:
            correct_entropy_items.append(correct_read2count[read])

        # correct entropy
        correct_nonFre_reads_num = len(non_frequent_correct_reads)
        shared_vars2 = correct_nonFre_reads_num, n, 1

        with WorkerPool(self.config.num_workers, shared_objects=shared_vars2, start_method='fork') as pool:
            correct_entropy_lst = pool.map(self.entropy_item, correct_entropy_items) 
        correct_entropy = sum(correct_entropy_lst)
        ##################################################################################
        new_reads = correct_unique_reads - raw_unique_reads
        new_reads_num = len(new_reads)
        self.logger.info("Wrongly introduced {} new reads".format(new_reads_num))

        # information gain visualisation
        correct2raw_diff = []
        common_non_frequent_reads = non_frequent_raw_reads & non_frequent_correct_reads

        common_non_frequent_raw_items = []
        common_non_frequent_correct_items = []
        for read in common_non_frequent_reads:
            common_non_frequent_raw_items.append(raw_read2count[read])
            common_non_frequent_correct_items.append(correct_read2count[read])

        with WorkerPool(self.config.num_workers, shared_objects=shared_vars1, start_method='fork') as pool:
            raw_common_lst = pool.map(self.entropy_item, common_non_frequent_raw_items) 

        with WorkerPool(self.config.num_workers, shared_objects=shared_vars2, start_method='fork') as pool:
            correct_common_lst = pool.map(self.entropy_item, common_non_frequent_correct_items) 

        if new_reads_num > 0:
            correct2raw_diff.extend([np.nan] * new_reads_num)

        ###### removed reads
        # removed_reads = non_frequent_raw_reads - non_frequent_correct_reads
        # if len(removed_reads) > 0:
        #     correct2raw_diff.extend([1000000] * len(removed_reads))

        for i, j in zip(raw_common_lst, correct_common_lst):
            correct2raw_diff.append(i - j)

        shuffle(correct2raw_diff)
        self.gain2heatmap(correct2raw_diff)

        return [raw_entropy, correct_entropy]
    '''

    '''
    def fft_spectrum(self, freq):
        return freq * freq

    def spectral_entropy(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = list(raw_read2count.keys())
        correct_unique_reads = list(correct_read2count.keys())
        n = len(raw_unique_reads)

        in_counts = []
        out_counts = []
        counts = []
        for read in correct_unique_reads:
            correct_count = correct_read2count[read]
            p = correct_count / total_num  
            if read in raw_unique_reads:
                in_counts.append(correct_count)
            else:
                out_counts.append(correct_count)
            counts.append(correct_count)
        # print(ant.spectral_entropy(in_counts, sf=100, method='welch', normalize=True))
        # print(ant.spectral_entropy(out_counts, sf=100, method='welch', normalize=True))
        # print(ant.spectral_entropy(counts, sf=100, method='welch', normalize=True))
        with WorkerPool(self.config.num_workers) as pool:
            in_fft_vals = pool.map(self.fft_spectrum, rfft(in_counts).tolist())
        if len(out_counts) > 0:
            with WorkerPool(self.config.num_workers) as pool:
                out_fft_vals = pool.map(self.fft_spectrum, rfft(out_counts).tolist())

            val_sum = sum(in_fft_vals) + sum(out_fft_vals)
        else:
            val_sum = sum(in_fft_vals)

        shared_vars = val_sum, n, 1
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            in_entropy = pool.map(self.spectrum2entropy, in_fft_vals)
        if len(out_counts) > 0:
            shared_vars = val_sum, n, 10
            with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
                out_entropy = pool.map(self.spectrum2entropy, out_fft_vals)
            entropy = sum(in_entropy) + sum(out_entropy)
        else:
            entropy = sum(in_entropy)
        
        print(entropy)
        # return 

    def spectrum2entropy(self, shared_vars, freq):
        total_freq, n, delta = shared_vars
        p = freq / total_freq
        result = - delta * (p * cmath.log(p,2)/ cmath.log(n,2))
        return result

    def gain2heatmap0(self, correct2raw_diff):
        
        x_len = len(correct2raw_diff)
        x_len_sqrt = round(math.sqrt(x_len))
        xx = x_len_sqrt * x_len_sqrt
        val = x_len - xx
        if xx == x_len:
            m = x_len_sqrt
            n = x_len_sqrt
            x = np.array((correct2raw_diff))
            x_res=x.reshape(m, n)
 
        elif val < x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 1
            remain = (m - val) * [0]
            x = np.array((correct2raw_diff + remain))
            x_res = x.reshape(m, n)

        elif val > x_len_sqrt:
            m = x_len_sqrt
            n = x_len_sqrt + 2
            remain = (val - m) * [0]
            x = np.array((correct2raw_diff + remain))
            x_res=x.reshape(m, n)

        # Specify the indices of the two types of bad data
        outliers1 = np.where(x_res == -1000000)
        # outliers2 = np.where(x_res == 1000000)

        # Plot the heatmap
        sns.heatmap(x_res, cmap="Set3")

        # Plot the first type of outliers in red
        plt.scatter(outliers1[0], outliers1[1], color="red", s=200)

        # Plot the second type of outliers in blue
        # plt.scatter(outliers2[0], outliers2[1], color="blue", s=200)
        plt.savefig(os.path.join(self.config.result_dir, self.prefix + '_information_gain.png'), transparent=True)

        # # Plot the first type of bad data using a red color
        # fig, ax = plt.subplots(figsize=(6, 6))
        # x_res[bad_indices_1] = np.nan

        # cur_cmap = mpl.cm.get_cmap('Set3')
        # cur_cmap.set_under('red')
        # cur_cmap.set_bad(color='red')
        # heatmap1 = ax.imshow(x_res, cmap=cur_cmap)

        # # Plot the second type of bad data using a green color
        # x_res[bad_indices_2] = np.nan
        # cur_cmap = mpl.cm.get_cmap('Set3')
        # cur_cmap.set_over('black')
        # cur_cmap.set_bad(color='Black')
        # heatmap2 = ax.imshow(x_res, cmap=cur_cmap)

        # # image = ax.imshow(x_res, cmap=cur_cmap) #interpolation='none', , vmax=thre_max
        # # # image = ax.pcolor(x_res, cmap=cmap, vmin=thre_min, antialiased=True)

        # cbar = fig.colorbar(heatmap1, extend='both', shrink=0.8)
        # cbar.cmap.set_under('red')
        # cbar.cmap.set_over('black')
        # fig.tight_layout()
        # fig.savefig(os.path.join(self.config.result_dir, self.prefix + '_information_gain.png'), transparent=True)

        return

    '''
    '''
    def total_variation_differnce(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = list(raw_read2count.keys())
        correct_unique_reads = list(correct_read2count.keys())
        total_unique_reads = set(raw_unique_reads + correct_unique_reads)
        variation = 0
        dif_vals = []
        for read in total_unique_reads:
            if read in raw_unique_reads:
                raw_count = raw_read2count[read]
            else:
                raw_count = 0
            if read in correct_unique_reads:
                correct_count = correct_read2count[read]
            else:
                correct_count = 0
            q = raw_count / total_num
            p = correct_count / total_num   
            dif_vals.append(p-q) 
        return norm(array(dif_vals))

    def relative_entropy(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = list(raw_read2count.keys())
        correct_unique_reads = list(correct_read2count.keys())
        total_unique_reads = raw_unique_reads + correct_unique_reads
        relative_entropy = 0
        relative_entropy_p2q = 0
        for read in total_unique_reads:
            if read in raw_unique_reads:
                raw_count = raw_read2count[read]
            else:
                raw_count = 0
            if read in correct_unique_reads:
                correct_count = correct_read2count[read]
            else:
                correct_count = 0
            q = raw_count / total_num
            p = correct_count / total_num   
            if q != 0:
                relative_entropy += p*math.log2(p/q)   
            if p != 0: 
                relative_entropy_p2q += q*math.log2(q/p) 
        print(relative_entropy, relative_entropy_p2q)
        return 

    '''

    '''
    def read_counts_entropy(self, read_count, total_num):
        """
        read count entropy calculation

        Args:
            read_count (dict): A dictionary to save reads (keys) and its counts (values) for dataset
            total_num (int): the number of total reads

        Returns:
            float: read count entropy 
        """
        entropy_val = 0
        # n = len(read_count.values())
        for count in read_count.values():
            p = count / total_num
            entropy_val += -(p * math.log2(p))
            # entropy_val += -math.log2(p)
            # entropy_val += -(p * math.log2(p) / math.log2(n))
            # entropy_val += -(p * math.log10(p))
        return entropy_val

    def noise2signal_entropy(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = list(raw_read2count.keys())
        correct_unique_reads = list(correct_read2count.keys())

        new_reads = list(set(correct_unique_reads) - set(raw_unique_reads))
        raw_unique_reads.sort()
        new_reads.sort()

        raw_counts = []
        raw_kept_counts = []
        correct_kept_counts = []
        
        raw_entropy_items = []

        correct2raw_diff = []

        correct_new_counts = []
        for read in raw_unique_reads:
            raw_counts.append(raw_read2count[read])
            if read in correct_unique_reads:
                correct_kept_counts.append(correct_read2count[read])
                raw_kept_counts.append(raw_read2count[read])
                correct2raw_diff.append( correct_read2count[read] - raw_read2count[read])
            else:
                # correct_removed_counts.append(0) # noise of identical reads removed by marked as 0
                # raw_removed_counts.append(raw_read2count[read])
                correct2raw_diff.append( 0 - raw_read2count[read])
                raw_entropy_items.append(raw_read2count[read])

        new_reads_num = len(new_reads)
        self.logger.info("Wrongly introduced {} new reads".format(new_reads_num))
        if new_reads_num > 0:
            entropy_new_items = []
            for read in new_reads:
                correct_new_counts.append(correct_read2count[read])
                correct2raw_diff.append(np.nan)
                entropy_new_items.append(np.nan)

        # raw entropy
        n = len(raw_unique_reads)
        shared_vars = total_num, n, 1
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            raw_entropy_lst = pool.map(self.entropy_item, raw_counts)
        raw_entropy = sum(raw_entropy_lst) 

        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            raw_kept_entropy_lst = pool.map(self.entropy_item, raw_kept_counts)
        # correct entropy
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            correct_kept_entropy_lst = pool.map(self.entropy_item, correct_kept_counts) 

        if new_reads_num > 0:
            shared_vars = total_num, n, self.config.delta
            with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
                correct_new_entropy_lst = pool.map(self.entropy_item, correct_new_counts)
            correct_entropy = sum(correct_new_entropy_lst) + sum(correct_kept_entropy_lst)
        else:
            correct_entropy = sum(correct_kept_entropy_lst)
            # correct_new_entropy_lst = []

        # shuffle(correct2raw_diff)
        # self.gain2heatmap(correct2raw_diff)
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars, start_method='fork') as pool:
            raw_entropy_items_lst = pool.map(self.entropy_item, raw_entropy_items)
        for i in raw_entropy_items_lst:
            if i <=0:
                print('no')
        entropy_item_lst = []
        for i, j in zip(raw_kept_entropy_lst, correct_kept_entropy_lst):
            entropy_item_lst.append(i - j)
        entropy_item_lst.extend(raw_entropy_items_lst)
        if new_reads_num > 0:
            entropy_item_lst.extend(entropy_new_items)
        shuffle(entropy_item_lst)
        self.gain2heatmap(entropy_item_lst)

        return [raw_entropy, correct_entropy]

    def entropy_item(self, shared_vars, freq):
        total_freq, n, delta = shared_vars
        if freq > 0:
            p = freq / total_freq
            result = - delta * (p * math.log2(p)/ math.log2(n))
        return result


    def entropy_item(self, shared_vars, freq):
        total_freq, n = shared_vars
        if freq > 0:
            p = freq / total_freq
            # result = - p * math.log2(p) / math.log2(n)
            result = - p * math.log2(p)
        return result

    def noise2read_entropy(self, raw_read2count, correct_read2count, total_num):

        raw_unique_reads = set(raw_read2count.keys())
        correct_unique_reads = set(correct_read2count.keys())

        frequent_reads = set([k for k, v in raw_read2count.items() if v > self.config.high_freq_thre])

        # raw dataset
        non_frequent_raw_reads = raw_unique_reads - frequent_reads
        raw_entropy_items = []
        for read in non_frequent_raw_reads:
            raw_entropy_items.append(raw_read2count[read])

        raw_nonFre_reads_total_num = sum(raw_entropy_items)
        # raw entropy
        non_frequent_raw_unique_num = len(non_frequent_raw_reads)

        shared_vars1 = raw_nonFre_reads_total_num, non_frequent_raw_unique_num
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars1, start_method='fork') as pool:
            raw_entropy_lst = pool.map(self.entropy_item, raw_entropy_items)
        raw_entropy = sum(raw_entropy_lst) 

        # correct dateset
        non_frequent_correct_reads = correct_unique_reads - frequent_reads
        correct_entropy_items = []
        for read in non_frequent_correct_reads:
            correct_entropy_items.append(correct_read2count[read])

        # correct entropy
        correct_nonFre_reads_total_num = sum(correct_entropy_items)

        # correct_unique_num = len(non_frequent_correct_reads)
        shared_vars2 = correct_nonFre_reads_total_num, non_frequent_raw_unique_num

        with WorkerPool(self.config.num_workers, shared_objects=shared_vars2, start_method='fork') as pool:
            correct_entropy_lst = pool.map(self.entropy_item, correct_entropy_items) 
        correct_entropy = sum(correct_entropy_lst)
        ##################################################################################
        #information gain (\delta I) heatmap
        new_reads = correct_unique_reads - raw_unique_reads
        new_reads_num = len(new_reads)
        self.logger.info("Wrongly introduced {} new reads".format(new_reads_num))

        raw_kept_counts = []
        correct_kept_counts = []
        kept_reads = correct_unique_reads & raw_unique_reads
        for read in kept_reads:
            correct_kept_counts.append(correct_read2count[read])
            raw_kept_counts.append(raw_read2count[read])

        # unique_kept_reads_num = len(kept_reads)
        raw_unique_num = len(raw_unique_reads)
        shared_vars3 = total_num, raw_unique_num
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars3, start_method='fork') as pool:
            raw_kept_entropy_lst = pool.map(self.entropy_item, raw_kept_counts)

        correct_kept_total_num = sum(correct_kept_counts)
        # shared_vars4 = correct_kept_total_num, unique_kept_reads_num
        shared_vars4 = total_num, raw_unique_num
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars4, start_method='fork') as pool:
            correct_kept_entropy_lst = pool.map(self.entropy_item, correct_kept_counts) 

        raw_removed_reads = raw_unique_reads - correct_unique_reads
        raw_removed_items = []
        for read in raw_removed_reads:
            raw_removed_items.append(raw_read2count[read])
        # removed_total_num = sum(raw_removed_items)
        # unique_removed_num = len(raw_removed_reads)
        shared_vars5 = total_num, raw_unique_num
        with WorkerPool(self.config.num_workers, shared_objects=shared_vars5, start_method='fork') as pool:
            raw_removed_entropy_items_lst = pool.map(self.entropy_item, raw_removed_items)
        for i in raw_removed_entropy_items_lst:
            if i <=0:
                print('Warning')
        entropy_item_lst = []
        for i, j in zip(raw_kept_entropy_lst, correct_kept_entropy_lst):
            entropy_item_lst.append(i - j)
        entropy_item_lst.extend(raw_removed_entropy_items_lst)
        if new_reads_num > 0:
            entropy_item_lst.extend([np.nan] * new_reads_num)
        shuffle(entropy_item_lst)
        self.gain2heatmap(entropy_item_lst)
        return [raw_entropy, correct_entropy]

    '''