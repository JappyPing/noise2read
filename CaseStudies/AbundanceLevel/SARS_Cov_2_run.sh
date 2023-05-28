# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-27 23:04:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-27 13:03:07

virus_name="SARS_Cov_2"                                                 # the name of the virus
ref="../SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta"                               # the reference with full path
raw_r1="../SARS_Cov_2/raw/D18_SRR11092062_reduced_r1.fastq"              # original read 1
raw_r2="../SARS_Cov_2/raw/D19_SRR11092062_reduced_r2.fastq"              # original read 2
correct_r1="../SARS_Cov_2/corrected/D18_SRR11092062_reduced_r1_corrected_corrected.fastq"    # corrected read 1
correct_r2="../SARS_Cov_2/corrected/D19_SRR11092062_reduced_r2_corrected_corrected.fastq"    # corrected read 2

mkdir -p ./result/SARS_Cov_2/raw
./get_coverage.sh -r ${ref} -1 ${raw_r1} -2 ${raw_r2} -o ./result/SARS_Cov_2/raw

mkdir -p ./result/SARS_Cov_2/correct
./get_coverage.sh -r ${ref} -1 ${correct_r1} -2 ${correct_r2} -o ./result/SARS_Cov_2/correct

raw_cvg="./result/SARS_Cov_2/raw/prn_cvg.txt"
correct_cvg="./result/SARS_Cov_2/correct/prn_cvg.txt"

python ./SARS_Cov_2_draw.py ${virus_name} ${raw_cvg} ${correct_cvg}