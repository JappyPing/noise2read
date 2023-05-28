# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-04-27 23:04:36
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-27 21:10:06

virus_name="Monkeypox"                                                 # the name of the virus
ref="../Monkeypox/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta"                               # the reference with full path
raw_r1="../Monkeypox/raw/D20_SRR22085311_1.fastq"              # original read 1
raw_r2="../Monkeypox/raw/D21_SRR22085311_2.fastq"              # original read 2
correct_r1="../Monkeypox/corrected/D20_SRR22085311_1_corrected_corrected.fastq"    # corrected read 1
correct_r2="../Monkeypox/corrected/D21_SRR22085311_2_corrected_corrected.fastq"    # corrected read 2

mkdir -p ./result/Monkeypox/raw
./get_coverage.sh -r ${ref} -1 ${raw_r1} -2 ${raw_r2} -o ./result/Monkeypox/raw

mkdir -p ./result/Monkeypox/correct
./get_coverage.sh -r ${ref} -1 ${correct_r1} -2 ${correct_r2} -o ./result/Monkeypox/correct

raw_cvg="./result/Monkeypox/raw/prn_cvg.txt"
correct_cvg="./result/Monkeypox/correct/prn_cvg.txt"

python ./monkey_draw.py ${virus_name} ${raw_cvg} ${correct_cvg}