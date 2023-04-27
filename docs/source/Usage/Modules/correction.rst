Correction
----------
#. Run with ini configuration
   
.. code-block:: console

    nois2read -m correction -c config.ini

Or

#. Run the commands only 

* Evaluating without ground truth data

.. code-block:: console

    nois2read -m correction -i path/to/raw.fastq -d ./results/ -p 60

* Evaluating with ground truth data

.. code-block:: console

    nois2read -m correction -i path/to/raw.fastq -t path/to/true.fastq -r path/to/corrected.fastq -d ./results/ 

* Training with GPU (default CPU)
  
.. code-block:: console

    nois2read -m correction -i path/to/raw.fastq -d ./results/ -g gpu_hist