simulated_miRNAs
----------------

Examples for correcting simulated miRNAs data with mimic UMIs by noise2read

Take data sets `D14 and D16 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EjBTpjExiShHg0kO72fVpzABn_Krd0K61xdLlK5_03JB5A?e=5GXsg8>`_) as examples.

* noise2read installation
   
Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

* Clone the codes with datasets in github

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/Examples/simulated_miRNAs

* Reproduce the evaluation results for D14 and D16 from raw, true and corrected datasets

.. code-block:: console

    noise2read -m evaluation -i ./simulated_miRNAs/raw/D14_umi_miRNA_mix.fa -t ./simulated_miRNAs/true/D14_umi_miRNA_mix.fa -r ./simulated_miRNAs/correct/D14_umi_miRNA_mix.fasta -d ./result
    noise2read -m evaluation -i ./simulated_miRNAs/raw/D16_umi_miRNA_subs.fa -t ./simulated_miRNAs/true/D16_umi_miRNA_subs.fa -r ./simulated_miRNAs/correct/D16_umi_miRNA_subs.fasta -d ./result

* correcting D14

  * with high ambiguous errors correction and using GPU for training (running about 4 mins with 26 cores and GPU)

  .. code-block:: console

      noise2read -m correction -c ../../config/D14.ini -a True -g gpu_hist

  * without high ambiguous errors correction and using CPU (default) for training (running about 4 mins with 26 cores)

  .. code-block:: console

      noise2read -m correction -c ../../config/D14.ini -a False

.. note:: 

    The latest noise2read  runs fast but produces slightly different corrected result from these under Examples/simulated_miRNAs/correct

* correcting D16

  * with high ambiguous errors correction and using GPU for training (running about 3 mins with 26 cores and GPU)

  .. code-block:: console

      noise2read -m correction -c ../../config/D16.ini -a True -g gpu_hist

  * without high ambiguous errors correction and using CPU (default) for training (running about 3 mins with 26 cores)

  .. code-block:: console

      noise2read -m correction -c ../../config/D16.ini -a False

.. note:: 

    The latest noise2read runs fast but produces slightly different corrected result from these under Examples/simulated_miRNAs/correct