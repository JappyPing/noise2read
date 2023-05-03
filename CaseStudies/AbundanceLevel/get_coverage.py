# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-19 18:51:21
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 15:10:56
import re
import sys
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
# import scipy.interpolate
# import statistics
import matplotlib.pyplot as plt

dict_tonu = {'A': 1, 'C': 2, 'T': 3, 'G': 4, 'N': 5, '-': 6}
dict_tole = dict(zip(dict_tonu.values(), dict_tonu.keys()))

def read_sam(R_file,freq=False):
	read_list = []
	with open(R_file,"r") as rf:
		for line in rf:
			if len(line.strip().split(" ")) == 5:
				r = pd.read_csv(R_file, delimiter=' ', names=['ID', 'strand', 'sta_p', 'sam_q', 'cigar'], encoding='unicode_escape')
			else:
				assert "sam file format error, should contain 5 or 6 columns of data. please re run find_sub.sh"
			break
	read_number = r.shape[0]
	# print("r is", r)
	for i in range(read_number):
		read = [str(r["ID"].loc[i]), int(r["strand"].loc[i]), int(r["sta_p"].loc[i]), str(r["sam_q"].loc[i]),
						  str(r["cigar"].loc[i])]

		read_list.append(read)
	return read_list

def narrow_reads(ref, narrowed_read, out_dir, brute_force=True):
	# global narrowed_read,  half_real_reads, half_real_ID
	
	print(len(narrowed_read), " 100% M reads in narrowed_extract.sam")
	true_total_match = []
	cigar_match = []
	for ri,rl in enumerate(narrowed_read):
		index = rl[2] - 1
		if ref[index:index + len(rl[3])] == rl[3]:
			true_total_match.append(rl)
		elif re.search('^[0-9]+m$', rl[4])  is not  None:
			cigar_match.append(rl)

	print(len(cigar_match)," reads have total match cigar string")
	print(len(true_total_match), " truely matched reads in real_narrowed_extract.sam")
	
	with open(out_dir + "real_narrowed_extract.sam", "w+") as nf1:
		for line in true_total_match:
			nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")
			
	# with open(out_dir + "narrowed_extract.sam", "w+") as nf1:
	# 	for line in cigar_match:
	# 		nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")

	# narrowed_read = true_total_match
	# mated = []
	# for rl in narrowed_read:
	# 	if rl[0] in ID_count.keys():
	# 		ID_count[rl[0]] += 1
	# 	else:
	# 		ID_count[rl[0]] = 1
			
	# for mate_rl in narrowed_read:
	# 	if ID_count[mate_rl[0]] == 2:
	# 		mated.append(mate_rl)

	# narrowed_read = true_total_match
	mated = []
	ID_count = {}
	for rl in true_total_match:
		if rl[0] in ID_count.keys():
			ID_count[rl[0]] += 1
		else:
			ID_count[rl[0]] = 1
			
	for mate_rl in true_total_match:
		if ID_count[mate_rl[0]] == 2:
			mated.append(mate_rl)
		elif ID_count[mate_rl[0]] > 2 :
			print("error")

	# narrowed_read = mated
	
	print(len(mated), " reads in paired_real_narrowed_extract.sam")
	print(len(ID_count))
	
	# with open(out_dir + "paired_real_narrowed_extract.sam", "w+") as nf1:
	# 	for line in narrowed_read:
	# 		nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")

	with open(out_dir + "paired_real_narrowed_extract.sam", "w+") as nf1:
		for line in mated:
			nf1.write(line[0] + " " + str(line[1]) + " " + str(line[2]) + " " + line[3] + " " + line[4] + "\n")

	if len(mated) == 0:
		print("no read satisfies condition, exiting")
		exit(-2)
		
	return cigar_match, true_total_match, mated

def matrix_from_readlist(readlist):
	insertion_reads = {}
	row_l = []
	col_l = []
	val_l = []
	narrowed = []
	maxposition = 0  # debug variable
	included_i = 0

	for i in range(len(readlist)):
		iread = readlist[i]
		index = iread[2] - 1
		sam_q = iread[3]
		tmp_length = len(iread[3])
		# configuring sam_q_num for matrix
		sam_q_num = []
		for j in range(0,len(sam_q)):
			sam_q_num.append(dict_tonu[sam_q[j]])
		val_l.extend(sam_q_num)
		row_tmp = [int(included_i)] * tmp_length
		row_l.extend(row_tmp)
		col_tmp = [n for n in range(index, index + tmp_length)]
		col_l.extend(col_tmp)
		if (index + len(sam_q)) > maxposition:
			maxposition = index + tmp_length
			maxindex = index
		included_i += 1
	target_matrix =	sp.coo_matrix((val_l, (row_l,col_l))).tocsc()  # matrix
	if len(val_l) > 0:
		csc = target_matrix
	else:
		csc = sp.coo_matrix((0,0)).tocsc()

	return target_matrix


def calc_coverage(ref_file, sam_file):

	readlist = read_sam(sam_file)
	# matrix = matrix_from_readlist(readlist)

	ref = ""
	with open(ref_file, 'r') as refg:
		for line in refg:
			if ">" not in line:
				ref += line.strip()

	narrowed, real_narrowed, paired_real_narrowed = narrow_reads(ref, readlist, "", True)

	# real_narrowed_matrix = matrix_from_readlist(real_narrowed)
	paired_real_narrowed_matrix = matrix_from_readlist(paired_real_narrowed)
	# rn_cvg_list = []
	prn_cvg_list = []

	# for i in range(0, len(ref)):
	# 	if i < real_narrowed_matrix.shape[1]:
	# 		tmp = np.squeeze(real_narrowed_matrix.getcol(i).toarray())
	# 		tmp_count = np.bincount(tmp)[1:]
	# 	else:
	# 		tmp_count = [0]
	# 	rn_cvg_list.append(sum(tmp_count))

	for i in range(0, len(ref)):
		if i < paired_real_narrowed_matrix.shape[1]:
			tmp = np.squeeze(paired_real_narrowed_matrix.getcol(i).toarray())
			tmp_count = np.bincount(tmp)[1:]
		else:
			tmp_count = [0]
		prn_cvg_list.append(sum(tmp_count))

	# return rn_cvg_list, prn_cvg_list
	return prn_cvg_list

def write_coverage_output(ref_file, sam_file):
	prn_cvg_list = calc_coverage(ref_file, sam_file)
	with open("prn_cvg.txt", "w+") as rncf:
		for i1, v1 in enumerate(prn_cvg_list):
			rncf.write(str(i1 + 1) + ": " + str(v1) + ", ")
	return "prn_cvg.txt"

if __name__ == '__main__':
	ref_file = sys.argv[1]
	sam_file = sys.argv[2]

	cov_file = write_coverage_output(ref_file, sam_file)

