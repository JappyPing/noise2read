Simulation
----------
Take the dataset `SRR12060401 <https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR12060401&display=data-access>`_ as an example to generate simulated data with mimic UMIs.

* Download dataset

.. code-block:: console

    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12060401/SRR12060401
    
* Open the dataset using fastq-dump

.. code-block:: console

    fastq-dump --split-3 SRR12060401

* Configuration

  * Download `D9_simi.ini <https://github.com/Jappy0/noise2read/blob/master/configs/D9_D13/simulation/D9_simi.ini>`_.

  **Or**

  * Create a file and copy below
    
.. code-block:: console

        [Paths]
        result_dir = ./D9_sim_test/

        [SourceInputData]
        input_file = ./data/SRR12060401.fastq

        [General]
        num_workers = 60
        chunks_num = 100

        [GraphSetup]
        high_freq_thre = 4
        max_error_freq = 4

        [Simulation]
        min_freq = 4
        min_read_count = 30
        substations = True
        indels = False
        error_rate1 = 0.09
        error_rate2 = 0.02
        sim_random_state = 42

* Run
  
.. code-block:: console

    noise2read -m simulation -c D9_simi.ini