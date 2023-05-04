ABEs_CBEs
---------

Correcting outcome sequence of ABEs and CBEs by noise2read

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies
    mkdir ABEs_CBEs
    cd ABEs_CBEs

* Download datasets `D32_D33 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EmjKFVI9QklJrR8Xe0YJP1kBEq8F_SPeUa-Xwx98JQZRNw>`_.

* Using noise2read to correct the datasets. The running time of each experiment is about 13 minutes using 26 cores and GPU for training.

.. code-block:: console

    noise2read -m correction -i ./D32_D33/raw/D32_ABE_outcome_seqs.fasta -a False -d ./ABE/
    noise2read -m correction -i ./D32_D33/raw/D33_CBE_outcome_seqs.fasta -a False -d ./CBE/