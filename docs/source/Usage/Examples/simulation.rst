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

  * Download `simulation.ini <https://raw.githubusercontent.com/Jappy0/noise2read/master/examples/simulation.ini>`_.

  **Or**

  * Create a file and copy below
    
  .. code-block:: console

      [Paths]
      ResultDir = "./result/" # set output directory

      [SourceInputData]
      input_file = ./SRR12060401.fastq # set your data to be corrected

      [General]
      num_workers = -1 # if num_workers = -1 or 0, noise2read will use all the available cpus 
      verbose = True 
      min_iters = 100
      iso_change_detail = True
      top_n = 100

      [GraphSetup]
      high_freq_thre = 5
      max_error_freq = 4
      save_graph = False
      graph_visualization = False
      drawing_graph_num = 50

      [EmbeddingSetup]
      entropy_kmer = 3
      entropy_q = 2
      kmer_freq = 3
      read_type = DNA

      [ModelTuningSetup]
      n_trials = 1
      n_estimators = 10 
      test_size = 0.1 # default        
      random_state = 32 # default  
      tree_method = 'auto'
      learning_rate_min = 1e-3 # default     
      learning_rate_max = 1e-1 # default 
      max_depth_min = 3 # default     
      max_depth_max = 15 # default     
      max_depth_step = 1 # default 
      num_boost_round_min = 200 # default     
      num_boost_round_max = 300 # default     
      num_boost_round_step = 10 # default 
      subsample_min = 0.8 # default     
      subsample_max = 1 # default     
      colsample_bytree_min = 0.8 # default     
      colsample_bytree_max = 1 # default     
      verbose_eval = True
      # xgboostclassifier seed
      seed = 32 # default 
      # optuna best trial accuracy
      best_accuracy = 0.75

      [Simulation]
      substations = True
      indels = False
      error_rate = 0.001

* Run
  
.. code-block:: console

    noise2read -m simulation -c simulation.ini