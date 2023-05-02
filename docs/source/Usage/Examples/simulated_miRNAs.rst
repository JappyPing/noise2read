simulated_miRNAs
----------------

These examples implement the results for correcting simulated miRNAs data with mimic UMIs (`D14 and D16 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EjBTpjExiShHg0kO72fVpzABn_Krd0K61xdLlK5_03JB5A?e=5GXsg8>`_) by nois2read.

* nois2read installation
   
Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

* Clone the codes with datasets in github

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/Examples/simulated_miRNAs

* correcting D14

  * with high ambiguous errors correction

  .. code-block:: console

      nois2read -m correction -c ../../config/simulated_miRNA.ini -i ./raw/D14_umi_miRNA_mix.fa.fastq -t ./true/D14_umi_miRNA_mix.fa.fastq -a True

  * without high ambiguous errors correction

  .. code-block:: console

      nois2read -m correction -c ../../config/simulated_miRNA.ini -i ./raw/D14_umi_miRNA_mix.fa.fastq -t ./true/D14_umi_miRNA_mix.fa.fastq -a False

* correcting D16

  * with high ambiguous errors correction

  .. code-block:: console

      nois2read -m correction -c ../../config/simulated_miRNA.ini -i ./raw/D16_umi_miRNA_mix.fa.fastq -t ./true/D16_umi_miRNA_mix.fa.fastq -a True

  * without high ambiguous errors correction

  .. code-block:: console

      nois2read -m correction -c ../../config/simulated_miRNA.ini -i ./raw/D16_umi_miRNA_mix.fa.fastq -t ./true/D16_umi_miRNA_mix.fa.fastq -a False