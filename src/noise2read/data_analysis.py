# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-01-30 09:35:18
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-08-29 10:14:41

from collections import Counter
import collections
import math
import os
import editdistance
import xlsxwriter
from noise2read.utils import *
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from random import shuffle
from mpire import WorkerPool
from matplotlib.ticker import MaxNLocator
# import gc

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

        for item in id_lst:
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
        #gc.collect()
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
        #gc.collect()
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
        #gc.collect()
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
        #gc.collect()
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
        #gc.collect()
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
        fig, ax = plt.subplots(figsize=(8, 8))
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
        cbar = fig.colorbar(image, extend='min', shrink=0.7)
        cbar.cmap.set_under('red')
        # cbar.cmap.set_over('red')
        # cbar.set_ticks(np.linspace(x_res.min(), x_res.max(), 8))
        cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=8))
        cbar.ax.tick_params(labelsize=18)
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=18)
        fig.tight_layout()
        fig.savefig(os.path.join(self.config.result_dir, self.prefix + '_information_gain.png'), transparent=True)
        #gc.collect()
        return

    def entropy_item(self, total_freq, freq):
        if freq > 0:
            p = freq / total_freq
            result = - p * math.log2(p)
        #gc.collect()
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

        if self.config.reads_chunks_num == 1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=raw_nonFre_reads_total_num, start_method='fork') as pool:
                    raw_entropy_lst = pool.map(self.entropy_item, raw_entropy_items)
                raw_entropy = sum(raw_entropy_lst) 
                #gc.collect()
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        else:
            chunk_size = len(raw_entropy_items) // self.config.reads_chunks_num
            groups = [raw_entropy_items[i:i+chunk_size] for i in range(0, len(raw_entropy_items), chunk_size)]
            for group in groups:
                try:
                    with WorkerPool(self.config.num_workers, shared_objects=raw_nonFre_reads_total_num, start_method='fork') as pool:
                        raw_entropy_lst = pool.map(self.entropy_item, group)
                    raw_entropy = sum(raw_entropy_lst) 
                    #gc.collect()
                except KeyboardInterrupt:
                    # Handle termination signal (Ctrl+C)
                    pool.terminate()  # Terminate the WorkerPool before exiting
                except Exception:
                    # Handle other exceptions
                    pool.terminate()  # Terminate the WorkerPool before exiting
                    raise
        del raw_entropy_items, raw_entropy_lst, raw_nonFre_reads_total_num
        #############################################################
        # correct entropy
        correct_nonFre_reads_total_num = sum(correct_entropy_items)

        if self.config.reads_chunks_num == 1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=correct_nonFre_reads_total_num, start_method='fork') as pool:
                    correct_entropy_lst = pool.map(self.entropy_item, correct_entropy_items) 
                correct_entropy = sum(correct_entropy_lst)
                #gc.collect()
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        else:
            chunk_size = len(correct_entropy_items) // self.config.reads_chunks_num
            groups = [correct_entropy_items[i:i+chunk_size] for i in range(0, len(correct_entropy_items), chunk_size)]
            for group in groups:
                try:
                    with WorkerPool(self.config.num_workers, shared_objects=correct_nonFre_reads_total_num, start_method='fork') as pool:
                        correct_entropy_lst = pool.map(self.entropy_item, group) 
                    correct_entropy = sum(correct_entropy_lst)
                    #gc.collect()
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
        if self.config.reads_chunks_num == 1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                    raw_kept_entropy_lst = pool.map(self.entropy_item, raw_kept_counts)
                #gc.collect()
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        else:
            chunk_size = len(raw_kept_counts) // self.config.reads_chunks_num
            groups = [raw_kept_counts[i:i+chunk_size] for i in range(0, len(raw_kept_counts), chunk_size)]
            for group in groups:
                try:
                    with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                        raw_kept_entropy_lst = pool.map(self.entropy_item, group)
                    #gc.collect()
                except KeyboardInterrupt:
                    # Handle termination signal (Ctrl+C)
                    pool.terminate()  # Terminate the WorkerPool before exiting
                except Exception:
                    # Handle other exceptions
                    pool.terminate()  # Terminate the WorkerPool before exiting
                    raise

        del raw_kept_counts
        ########################################################

        if self.config.reads_chunks_num == 1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                    correct_kept_entropy_lst = pool.map(self.entropy_item, correct_kept_counts) 
                #gc.collect()
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        else:
            chunk_size = len(correct_kept_counts) // self.config.reads_chunks_num
            groups = [correct_kept_counts[i:i+chunk_size] for i in range(0, len(correct_kept_counts), chunk_size)]
            for group in groups:
                try:
                    with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                        correct_kept_entropy_lst = pool.map(self.entropy_item, group) 
                    #gc.collect()
                except KeyboardInterrupt:
                    # Handle termination signal (Ctrl+C)
                    pool.terminate()  # Terminate the WorkerPool before exiting
                except Exception:
                    # Handle other exceptions
                    pool.terminate()  # Terminate the WorkerPool before exiting
                    raise
        del correct_kept_counts
        ############################################################################

        if self.config.reads_chunks_num == 1:
            try:
                with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                    raw_removed_entropy_items_lst = pool.map(self.entropy_item, raw_removed_items)
                #gc.collect()
            except KeyboardInterrupt:
                # Handle termination signal (Ctrl+C)
                pool.terminate()  # Terminate the WorkerPool before exiting
            except Exception:
                # Handle other exceptions
                pool.terminate()  # Terminate the WorkerPool before exiting
                raise
        else:
            chunk_size = len(raw_removed_items) // self.config.reads_chunks_num
            groups = [raw_removed_items[i:i+chunk_size] for i in range(0, len(raw_removed_items), chunk_size)]
            for group in groups:
                try:
                    with WorkerPool(self.config.num_workers, shared_objects=total_num, start_method='fork') as pool:
                        raw_removed_entropy_items_lst = pool.map(self.entropy_item, group)
                    #gc.collect()
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
        #gc.collect()
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
        #gc.collect()
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
        #gc.collect()
        return   