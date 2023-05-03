AbundanceLevel
--------------

The instructions for the case studies of the abundance changes for the reference genome of SARS_Cov_2 and Monkeypox.

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies/AbundanceLevel

* Download the reference and original and corrected datasets from `D18_D21 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD18%5FD21&view=0>`_

and move the datasets to the folder \*/noise2read/CaseStudies/AbundanceLevel

* For the abundance analysis of SARS_Cov_2

  * Get the base coverage by perfectly matching the raw reads to the SARS-Cov-2 genome
    
  .. code-block:: console  

    ./get_coverage.sh -r ./D18_D21/SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta -1 ./D18_D21/SARS_Cov_2/raw/D18_SRR11092062_reduced_r1.fastq -2 ./D18_D21/SARS_Cov_2/raw/D19_SRR11092062_reduced_r2.fastq -o ./result/SARS_Cov_2/raw/
    python ./get_coverage.py ./D18_D21/SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta ./result/SARS_Cov_2/raw/paired_real_narrowed_extract.sam 
    mv prn_cvg.txt  ./result/SARS_Cov_2/raw/prn_cvg.txt

  * Get the base coverage by perfectly matching the corrected reads to the SARS-Cov-2 genome

  .. code-block:: console  

    ./get_coverage.sh -r ./D18_D21/SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta -1 ./D18_D21/SARS_Cov_2/corrected/D18_SRR11092062_reduced_r1.fastq -2 ./D18_D21/SARS_Cov_2/corrected/D19_SRR11092062_reduced_r2.fastq -o ./result/SARS_Cov_2/correct/
    python ./get_coverage.py ./D18_D21/SARS_Cov_2/ref/sars_cov_ref_MN996528.1.fasta ./result/SARS_Cov_2/correct/paired_real_narrowed_extract.sam
    mv prn_cvg.txt  ./result/SARS_Cov_2/correct/prn_cvg.txt

  * Draw the base coverage results before and after correction

  .. code-block:: console  

      python ./draw.py SARS_Cov_2 ./result/SARS_Cov_2/raw/prn_cvg.txt ./result/SARS_Cov_2/correct/prn_cvg.txt

* For the abundance analysis of Monkeypox

  * Get the base coverage by perfectly matching the raw reads to the SARS-Cov-2 genome
    
  .. code-block:: console  

    ./get_coverage.sh -r ./D18_D21/Monkeypox/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta -1 ./D18_D21/Monkeypox/raw/SRR22085311_1.fastq -2 ./D18_D21/Monkeypox/raw/SRR22085311_2.fastq -o ./result/Monkeypox/raw/
    python ./get_coverage.py ./D18_D21/Monkeypox/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta ./result/Monkeypox/raw/paired_real_narrowed_extract.sam 
    mv prn_cvg.txt  ./result/Monkeypox/raw/prn_cvg.txt

  * Get the base coverage by perfectly matching the corrected reads to the SARS-Cov-2 genome

  .. code-block:: console  

    ./get_coverage.sh -r ./D18_D21/Monkeypox/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta -1 ./D18_D21/Monkeypox/corrected/SRR22085311_1.fastq -2 ./D18_D21/Monkeypox/corrected/SRR22085311_2.fastq -o ./result/Monkeypox/correct/
    python ./get_coverage.py ./D18_D21/Monkeypox/ref/GCA_025947495.1_ASM2594749v1_genomic.fasta ./result/Monkeypox/correct/paired_real_narrowed_extract.sam 
    mv prn_cvg.txt  ./result/Monkeypox/correct/prn_cvg.txt

  * Draw the base coverage results before and after correction

  .. code-block:: console  

      python ./draw.py Monkeypox ./result/Monkeypox/raw/prn_cvg.txt ./result/Monkeypox/correct/prn_cvg.txt