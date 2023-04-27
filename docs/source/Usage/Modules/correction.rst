Correction
----------
1. Run with ini configuration
   
.. code-block:: console

    nois2read -m correction -c config.ini

Or

2. Run the commands only 

without ground truth data

.. code-block:: console

    nois2read -m correction -i path/to/raw.fastq -d ./results/ -p 60 -g gpu_hist

with ground truth data

.. code-block:: console

    nois2read -m evaluation -i path/to/raw.fastq -t path/to/true.fastq -r path/to/corrected.fastq -d ./results/ 