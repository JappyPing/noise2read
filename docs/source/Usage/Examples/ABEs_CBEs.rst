ABEs_CBEs
---------

Outcome sequence analysis for ABEs and CBEs.

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd CaseStudies/ABEs_CBEs

* Download datasets `D32_D33 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD32%5FD33&view=0>`_.

* Using noise2read to correct the datasets

.. code-block:: console

    noise2read -m correction -i ./D32_D33/raw/D32_ABE_outcome_seqs.fasta -a False -d ./ABE/
    noise2read -m correction -i ./D32_D33/raw/D33_CBE_outcome_seqs.fasta -a False -d ./CBE/