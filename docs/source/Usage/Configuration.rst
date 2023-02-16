Configuration
-------------

noise2read is a command-line interface (CLI) based tool to eliminate PCR and sequencing errors for short reads. It utilises CLI means and INI file for configuring the parameters. 

.. 1. Command line Options
.. <<<<<<<<<<<<<<<<<<<<<<<

1. Two Indispensable CLI setting
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

1.1 Module selection
====================
  
  Using noise2read, you must select one module from module list ["correction", "amplicon_correction", "mimic_umi", "real_umi", "umi_correction", "simulation", "evaluation"] first.

.. code-block:: console

    -m | --module module_name

1.2 Setting INI file or input dataset
=====================================

* INI file

.. code-block:: console

    -c | --config config.ini

* Input Read dataset

.. code-block:: console

    -i | --input data.fastq

1.2 Optional CLI setting
<<<<<<<<<<<<<<<<<<<<<<<<

You can configure some parameters using CLI mode with/without INI file configuration. INI file configuration can set all the parameters except for module selection. With INI file configuration, the following parameters settings in the INI file will be invalid when setting them using CLI mode. The following parameters settings

.. code-block:: console


.. code-block:: console


.. code-block:: console


.. code-block:: console

.. code-block:: console


1. INI file configuration
<<<<<<<<<<<<<<<<<<<<<<<<< 

The following contents is the default setting for an INI file configuration.

.. code-block:: console

    [Paths]
    ResultDir = "./result/"

    [SourceInputData]
    input_file = path/to/data.fastq
    ground_truth_data = path/to/data.fastq # optional

    [General]
    num_workers = -1
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

    [RealUMI]
    umi_start = 0
    umi_end = 12
    non_umi_start = 24

    [Amplicon]
    amplicon_low_freq = 50
    amplicon_high_freq = 1500
    amplicon_threshold_proba = 0.9
    amplicon_error_node_degree = 4

    [Simulation]
    substations = True
    indels = False
    error_rate = 0.001

