AbundanceLevel
--------------

The instructions for the case studies of the abundance changes for the reference genome of SARS_Cov_2 and Monkeypox.

* Clone the codes

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd CaseStudies/AbundanceLevel

* Download the reference and original and corrected datasets from `D18_D21 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD18%5FD21&view=0>`_

    and move the datasets to the folder \*/CaseStudies/AbundanceLevel

* For the abundance analysis of SARS_Cov_2

.. code-block:: console

    chmod +x SARS_Cov_2_run.sh
    ./SARS_Cov_2_run.sh

* For the abundance analysis of Monkeypox

.. code-block:: console

    chmod +x Monkeypox_run.sh
    ./Monkeypox_run.sh