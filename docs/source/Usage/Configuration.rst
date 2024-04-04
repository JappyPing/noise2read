Configuration
-------------

noise2read is a command-line interface (CLI) based tool to eliminate PCR and sequencing errors for short reads. It utilises CLI mode and INI file for configuring the parameters. 

1. Indispensable CLI setting
<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* Module selection

Using noise2read, we must select the module name from ["correction", "amplicon_correction", "mimic_umi", "real_umi", "umi_correction", "simulation", "evaluation"] first.

.. code-block:: console

  -m | --module module_name

* Setting configuration file or input dataset

   * configuration

   .. code-block:: console

       -c | --config config.ini

   * Input Read dataset

   .. code-block:: console

       -i | --input data.fastq

2. Optional CLI setting
<<<<<<<<<<<<<<<<<<<<<<<

You can set some parameters using CLI mode with/without INI file configuration. INI file configuration can set all the parameters except for module selection. The following parameters settings in the INI file will be invalid when setting them using CLI mode.

.. code-block:: console

    -u | --umi_file umi.fastq

.. code-block:: console

    -t | --true ground_truth.data.fastq

.. code-block:: console

    -r | --rectification corrected.data.fastq

.. code-block:: console

    -p | --parallel num_of_cpu_core

.. code-block:: console

    -a | --high_ambiguous True/False

.. code-block:: console

    -g | --tree_method gpu_hist/auto

.. code-block:: console

    -d | --directory */output_dir/

3. INI file configuration 
<<<<<<<<<<<<<<<<<<<<<<<<< 

The following are the default setting for all parameters of noise2read using an INI file configuration. There is no need to set these if you only want to use the default parameters which have been embedded in the program.

.. code-block:: console

    [Paths]
    ResultDir = "./result/"
    ; set output directory

    [SourceInputData]
    input_file = path/to/data.fastq 
    ; set your data to be corrected
    ground_truth_data = None
    ; only set when you have groundtruth data, otherwise comment it

    [General]
    num_workers = -1
    ; if num_workers = -1 or 0, noise2read will use all the available cpus 
    chunks_num = 100
    reads_chunks_num = 1
    verbose = True 
    iso_change_detail = False
    top_n = 100
    negative_sample_num = 300000
    min_read_len = 30

    [GraphSetup]
    high_freq_thre = 4
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
    high_ambiguous = True 
    proba_deviation = 0.95 
    iso_neg_high = True  

    [ModelTuningSetup]
    n_trials = 20
    n_estimators = 400 
    test_size = 0.1        
    random_state = 42  
    tree_method = auto
    learning_rate_min = 1e-3     
    learning_rate_max = 1e-1 
    max_depth_min = 3     
    max_depth_max = 15     
    max_depth_step = 1 
    subsample_min = 0.8     
    subsample_max = 1     
    colsample_bytree_min = 0.8     
    colsample_bytree_max = 1     
    verbose_eval = False
    xgboost_seed = 42 

    [RealUMI]
    umi_in_read = False
    umi_start = 0
    umi_end = 12
    non_umi_start = 24
    group_read_number = 10
    read_edit_dif = 2
    separator1 = '_'
    separator1_idx = 2
    separator2 = ' '
    separator2_idx = 0

    [Amplicon]
    amplicon_low_freq = 50
    amplicon_high_freq = 1500
    amplicon_threshold_proba = 0.85

    [Simulation]
    min_freq = 4
    min_read_count = 30
    error_rate1 = 0.09
    error_rate2 = 0.02

