'''
Author: Pengyao PING
Date: 2022-05-27 14:08:11
LastEditors: Pengyao PING
LastEditTime: 2022-08-04 00:26:06
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
import unittest
from unittest import result
from noise2read.data_generation import DataGneration

class TestDataGeneration(unittest.TestCase):
    def __init__(self, methodName: str = ...) -> None:
        super().__init__(methodName)
        input_filename = '/data/pping/Repo/srec/data/test/test.fastq'
        file_type = 'fastq'
        output_dir = '/data/pping/Repo/srec/data/test'
        self.DG = DataGneration(input_filename, file_type, output_dir, max_error_freq = 2)
    def test_sub_base(self):
        self.assertEqual('TCG', self.DG.sub_base('A'))  
        self.assertEqual('ACG', self.DG.sub_base('T'))
        self.assertEqual('ATG', self.DG.sub_base('C'))
        self.assertEqual('ACT', self.DG.sub_base('G'))   
        self.assertEqual('N', self.DG.sub_base('N'))    

    def test_seq2deletion(self):
        read = 'ACGTGA'
        del_list = self.DG.seq2deletion(read)
        self.assertEqual(del_list, set(["CGTGA", 'AGTGA', 'ACTGA', 'ACGGA', 'ACGTA','ACGTG']))

    def test_seq2insertion(self):
        read = 'ACGT'
        read_set = self.DG.seq2insertion(read)
        exp_result = set(['AACGT', 'GACGT', 'CACGT', 'TACGT', 'AACGT', 'AGCGT','ACCGT', 'ATCGT', 'ACAGT', 'ACGGT', 'ACCGT', 'ACTGT', 'ACGAT', 'ACGGT', 'ACGCT', 'ACGTT', 'ACGTA', 'ACGTG', 'ACGTC', 'ACGTT'])
        self.assertEqual(read_set, exp_result)     

    def test_seq2substitution(self):
        read1 = 'ACGT'
        del_list1 = self.DG.seq2substitution(read1)
        exp_result1 = set(['CCGT', 'GCGT', 'TCGT', 'AAGT', 'AGGT', 'ATGT', 'ACAT', 'ACCT', 'ACTT', 'ACGA', 'ACGG', 'ACGC'])
        self.assertEqual(del_list1, exp_result1)     
        read2 = 'ACGN'
        del_list2 = self.DG.seq2substitution(read2)
        exp_result2 = set(['CCGN', 'GCGN', 'TCGN', 'AAGN', 'AGGN', 'ATGN', 'ACAN', 'ACCN', 'ACTN'])
        self.assertEqual(del_list2, exp_result2)    

    def test_ErrorTypeClassification(self):
        line1 = ['ACGTA', 5, 'ACGTAA', 1]
        real_result1 = self.DG.ErrorTypeClassification(line1)
        exp_result1 = ['ACGTA', 5, 'N-A', 5, 'ACGTAA', 1]
        self.assertEqual(real_result1, exp_result1) 

        line1 = ['ACGTAG', 5, 'ACGTAA', 1]
        real_result1 = self.DG.ErrorTypeClassification(line1)
        exp_result1 = ['ACGTAG', 5, 'G-A', 5, 'ACGTAA', 1]
        self.assertEqual(real_result1, exp_result1) 

        line1 = ['TACGTA', 5, 'ACGTA', 1]
        real_result1 = self.DG.ErrorTypeClassification(line1)
        exp_result1 = ['TACGTA', 5, 'T-N', 0, 'ACGTA', 1]
        self.assertEqual(real_result1, exp_result1) 

        # line1 = ['TACG', 5, 'ACGTA', 1]
        # # real_result1 = self.DG.ErrorTypeClassification(line1)
        # # raise IOError("The editdistance of two reads in the input list must equal to one!")
        # self.assertRaises(IOError, self.DG.ErrorTypeClassification(line1)) 

    def test_EditDis1Seqs(self):
        read1 = 'ACGTA'
        del_list1 = self.DG.EditDis1Seqs(read1, 'A-G')
        exp_result1 = set(['GCGTA', 'ACGTG'])
        self.assertEqual(del_list1, exp_result1)     
        read2 = 'ACGT'
        del_list2 = self.DG.EditDis1Seqs(read2, 'N-A')
        exp_result2 = set(['AACGT', 'ACAGT', 'ACGAT', 'ACGTA'])
        self.assertEqual(del_list2, exp_result2) 

        read2 = 'ACGTT'
        del_list2 = self.DG.EditDis1Seqs(read2, 'T-N')
        exp_result2 = set(['ACGT', 'ACGT'])
        self.assertEqual(del_list2, exp_result2) 

if __name__ == '__main__':
    unittest.main()