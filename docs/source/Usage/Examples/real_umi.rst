Real_UMI
--------

Take the dataset `SRR1543694 <https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR1543694&display=data-access>`_ as an example to generate real sequencing data with UMI-based ground truth.

* Download dataset

.. code-block:: console

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1543694/SRR1543694
    
* Open the dataset using fastq-dump

.. code-block:: console

    fastq-dump --split-3 SRR1543694

* Configuration

  * Download it from `SRR1543694.ini <https://github.com/Jappy0/noise2read/blob/master/configs/D1_D8/real_umi/SRR1543964.ini>`_

  **Or**

  * create it by yourself and copy the following contents

  .. code-block:: console

        [Paths]
        result_dir = ./

        [SourceInputData]
        input_file = /path_to_data/SRR1543964.fastq

        [General]
        num_workers = -1

        [GraphSetup]
        high_freq_thre = 4
        max_error_freq = 4

        [RealUMI]
        umi_in_read = True
        umi_start = 0
        umi_end = 12
        non_umi_start = 24
        group_read_number = 10
        read_edit_dif = 2

* Run
 
.. code-block:: console

    noise2read -m real_umi -c SRR1543694.ini
