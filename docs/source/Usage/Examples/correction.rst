Correction
----------

Take dataset "D1_umi_SRR1543964.fastq" as an example, if you want to run the other datasets, change the dataset name in the configuration.

1. noise2read installation
<<<<<<<<<<<<<<<<<<<<<<<<<   

Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

2. Download datasets
<<<<<<<<<<<<<<<<<<<<

Download UMI-based ground truth datasets `raw <https://studentutsedu-my.sharepoint.com/:u:/g/personal/pengyao_ping_student_uts_edu_au/EZnprFyUT2xPgeIsgpZBam8BFyuxfnLwnquLx1ek7bCOIA?e=7G8z3S>`_ and `true <https://studentutsedu-my.sharepoint.com/:u:/g/personal/pengyao_ping_student_uts_edu_au/EVzmag9mPHhAl7WU4wdVcnQBgO1s-PHxR0AYvh59WMhcAg?e=xmPrKc>`_ from `D1_D8 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/ElxypUHIIqtDuyeQmmlZtQMBIzOa2YzFsMsqr7E6h0rVhQ?e=nWvTOh>`_

.. code-block:: console

  cd your_working_directory
  mkdir data & cd data
  mkdir raw & mkdir true
  cd your_working_directory

.. note:: 

  Please note that the datasets raw and true have the same file name, move them in the different folders raw and true

1. Configuration
<<<<<<<<<<<<<<<<

* Download it using

.. code-block:: console

  wget https://raw.githubusercontent.com/Jappy0/noise2read/master/examples/D1_D8_config.ini

**Or** 

* create a new file by yourself and copy the following contents

.. code-block:: console

    [Paths]
    ResultDir = "./result/" # set output directory

    [SourceInputData]
    input_file = ./data/D1_D8/raw/D1_umi_SRR1543964.fastq # set your data to be corrected
    ground_truth_data = ./data/D1_D8/true/D1_umi_SRR1543964.fastq # only set when you have groundtruth data, otherwise comment it

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

    [AmbiguousSetup]
    ambiguous_error_node_degree = 4
    high_ambiguous = False 
    # high ambiguous predict probability difference
    proba_deviation = 0.6  

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

4. Run
<<<<<<
    
.. code-block:: console

    noise2read -m correction -c D1_D8_config.ini
