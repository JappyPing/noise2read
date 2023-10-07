ABEs_CBEs
---------

Correcting outcome sequence of ABEs and CBEs by noise2read

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies
    mkdir ABEs_CBEs
    cd ABEs_CBEs

* Download datasets under data of `D32_D33 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EokIIeQd2nFHjlpurzDaBywB7Smy6Sm0dBR86GIJt0PSdg>`_.

* Using noise2read to correct the datasets. The running time of each experiment is about 13 minutes using 26 cores and GPU for training.

.. code-block:: console

    noise2read -m correction -i ./data/ABE_outcome_seqs.fasta -a False -d ./ABE/
    noise2read -m correction -i ./data/CBE_outcome_seqs.fasta -a False -d ./CBE/