# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-27 23:04:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-28 09:51:35

virus_name="Monkeypox"                                                 # the name of the virus
ref="./D18_D21/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta"                               # the reference with full path
raw_r1="./D18_D21/Monkeypox/raw/D20_SRR22085311_1.fastq"              # original read 1
raw_r2="./D18_D21/Monkeypox/raw/D21_SRR22085311_2.fastq"              # original read 2
correct_r1="./D18_D21/Monkeypox/corrected/D20_SRR22085311_1.fastq"    # corrected read 1
correct_r2="./D18_D21/Monkeypox/corrected/D21_SRR22085311_2.fastq"    # corrected read 2

./get_coverage.sh -r ${ref} -1 ${raw_r1} -2 ${raw_r2} -o ./raw/

./get_coverage.sh -r ${ref} -1 ${correct_r1} -2 ${correct_r2} -o ./correct/

raw_cvg="./raw/prn_cvg.txt"
correct_cvg="./correct/prn_cvg.txt"

python "${DIR}"/draw.py ${virus_name} ${raw_cvg} ${correct_cvg}