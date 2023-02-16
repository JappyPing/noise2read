'''
Author: Pengyao PING
Date: 2022-06-03 01:36:38
LastEditors: Pengyao PING
LastEditTime: 2022-08-04 00:26:03
Email: Pengyao.Ping@student.uts.edu.au
Description: 
'''
import unittest
from noise2read.reads2vectors import Reads2Vectors

class TestReads2Vectors(unittest.TestCase):
    def __init__(self):
        data_dir = '/data/pping/Repo/srec/data/test'
        lab_csv = '/data/pping/Repo/srec/data/test/test.csv'
        unlab_csv = '/data/pping/Repo/srec/data/test/test.csv'
        k_lst = '1, 2, 3'
        self.RV = Reads2Vectors(data_dir, lab_csv, unlab_csv, k_lst)
    
    def test_doubling_kmer_contribution(self):
        read = 'ACGTACTGAACT'
        k = 3
        kmers1 = ['ACG', 'ACG']
        res = self.RV.doubling_kmer_contribution(read, k, kmers1)
        self.assertEqual(res, 0)

        kmers2 = ['ACT', 'ACT']
        res2 = self.RV.doubling_kmer_contribution(read, k, kmers2)
        self.assertEqual(res2, 0)


if __name__ == '__main__':
    unittest.main()