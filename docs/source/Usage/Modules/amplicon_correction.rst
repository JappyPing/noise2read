Amplicon_correction
-------------------

Correction errors for Amplicon sequencing data

1. Run with ini configuration
   
.. code-block:: console

    nois2read -m amplicon_correction -c config.ini

Or

2. Run the commands only 

* Evaluating without ground truth data

.. code-block:: console

    nois2read -m amplicon_correction -i path/to/raw.fastq -d ./results/ -p 60

* Evaluating with ground truth data

.. code-block:: console

    nois2read -m amplicon_correction -i path/to/raw.fastq -t path/to/true.fastq -r path/to/corrected.fastq -d ./results/ 

* Training with GPU (default CPU)
  
.. code-block:: console

    nois2read -m amplicon_correction -i path/to/raw.fastq -d ./results/ -g gpu_hist