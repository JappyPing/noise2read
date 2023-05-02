ABEs_CBEs
---------

Outcome sequence analysis for ABEs and CBEs.

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies
    mkdir ABEs_CBEs
    cd ABEs_CBEs

* Download datasets `D32_D33 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EmjKFVI9QklJrR8Xe0YJP1kBEq8F_SPeUa-Xwx98JQZRNw`_.

* Using noise2read to correct the datasets

.. code-block:: console

    noise2read -m correction -i ./D32_D33/raw/D32_ABE_outcome_seqs.fasta -a False -d ./ABE/
    noise2read -m correction -i ./D32_D33/raw/D33_CBE_outcome_seqs.fasta -a False -d ./CBE/