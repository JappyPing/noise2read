# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-27 23:04:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-03 01:25:11

virus_name="SARS_Cov_2"                                                 # the name of the virus
ref="./D18_D21/SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta"                               # the reference with full path
raw_r1="./D18_D21/SARS_Cov_2/raw/D18_SRR11092062_reduced_r1.fastq"              # original read 1
raw_r2="./D18_D21/SARS_Cov_2/raw/D19_SRR11092062_reduced_r2.fastq"              # original read 2
correct_r1="./D18_D21/SARS_Cov_2/corrected/D18_SRR11092062_reduced_r1.fastq"    # corrected read 1
correct_r2="./D18_D21/SARS_Cov_2/corrected/D19_SRR11092062_reduced_r2.fastq"    # corrected read 2

./get_coverage.sh -r ${ref} -1 ${raw_r1} -2 ${raw_r2} -o ./raw/

./get_coverage.sh -r ${ref} -1 ${correct_r1} -2 ${correct_r2} -o ./correct/

raw_cvg="./raw/prn_cvg.txt"
correct_cvg="./correct/prn_cvg.txt"

python "${DIR}"/draw.py ${virus_name} ${raw_cvg} ${correct_cvg}