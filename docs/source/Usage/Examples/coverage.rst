Real_UMI
--------

Take SARS_Cov_2 data as an example

* Download the reference from 
`D18_D21 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD18%5FD21&view=0>`_

Or NCBI

.. code-block:: console

    wget https://www.ncbi.nlm.nih.gov/nuccore/MN996528.1?report=fasta

* Download original and corrected datasets
`D18_D21 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD18%5FD21&view=0>`_

* Get base coverage

  * get the codes
    
  .. code-block:: console

      git clone https://github.com/Jappy0/noise2read
      cd get_coverage

  * Setting the input info in "run.sh"

  .. code-block:: console

      virus_name="SARS_Cov_2"                                                 # the name of the virus
      ref="./ref/sars_cov_ref_MN996528.1.fasta"                               # the reference with full path
      raw_r1="./SARS_Cov_2/raw/D18_SRR11092062_reduced_r1.fastq"              # original read 1
      raw_r2="./SARS_Cov_2/raw/D19_SRR11092062_reduced_r2.fastq"              # original read 2
      correct_r1="./SARS_Cov_2/corrected/D18_SRR11092062_reduced_r1.fastq"    # corrected read 1
      correct_r2="./SARS_Cov_2/corrected/D19_SRR11092062_reduced_r2.fastq"    # corrected read 2

  * then run run.sh

  .. code-block:: console

      ./run.sh
  