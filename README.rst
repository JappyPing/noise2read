.. _noise2read-documentation:

.. image:: ./logo/logo.svg
   :align: center
   :target: https://noise2read.readthedocs.io/en/latest/

.. image:: https://readthedocs.org/projects/noise2read/badge/?version=latest
    :target: https://noise2read.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

Turn 'noise' to signal: accurately rectify millions of erroneous short reads through graph learning on edit distances
=====================================================================================================================

`noise2read <https://noise2read.readthedocs.io/en/latest/>`__, originated in a computable rule translated from PCR erring mechanism that: a rare read is erroneous if it has a neighboring read of high abundance, turns erroneous reads into their original state without bringing up any non-existing sequences into the short read set(<300bp) including DNA and RNA sequencing (DNA/RNA-seq), small RNA, unique molecular identifiers (UMI) and amplicon sequencing data.

Click `noise2read <https://noise2read.readthedocs.io/en/latest/>`__ to jump to its documentation
================================================================================================


Examples
========

These examples implement the results for correcting simulated miRNAs data with mimic UMIs (`D14 and D16 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EjBTpjExiShHg0kO72fVpzABn_Krd0K61xdLlK5_03JB5A?e=5GXsg8>`_) by nois2read.

* nois2read installation
   
Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

* Clone the codes with datasets in this repository

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