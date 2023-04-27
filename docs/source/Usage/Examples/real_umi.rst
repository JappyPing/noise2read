Real_UMI
--------

Take the dataset `SRR1543694 <https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR1543694&display=data-access>`_ as an example.

1. Download dataset

.. code-block:: console

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1543694/SRR1543694
    
Open the dataset using fastq-dump

.. code-block:: console

    fastq-dump --split-3 SRR1543694

2. Generate real sequencing data with UMI-based ground truth
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* Configuration file (config.ini)
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Download it from `real_umi_config.ini <>`

Or create it by yourself and copy the following contents

.. code-block:: console

    [Paths]
    ResultDir = "./result/" # set output directory

    [SourceInputData]
    input_file = ./SRR1543694.fastq # set your data to be corrected

    [RealUMI]
    umi_start = 0
    umi_end = 12
    non_umi_start = 24

* Run
 
.. code-block:: console

    nois2read -m real_umi -i config.ini
