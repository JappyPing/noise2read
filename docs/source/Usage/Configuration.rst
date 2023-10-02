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

The following contents are the default setting for noise2read using an INI file configuration. There is no need to set these if you only want to use the default parameters which have been embedded in the program.

.. code-block:: console

    [Paths]
    ResultDir = "./result/"
    ; set output directory

    [SourceInputData]
    input_file = path/to/data.fastq 
    ; set your data to be corrected
    ground_truth_data = path/to/data.fastq
    ; only set when you have groundtruth data, otherwise comment it

    [General]
    num_workers = -1
    ; if num_workers = -1 or 0, noise2read will use all the available cpus 
    verbose = True 
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
    proba_deviation = 0.6  

    [ModelTuningSetup]
    n_trials = 1
    n_estimators = 10 
    test_size = 0.1        
    random_state = 32  
    tree_method = 'auto'
    learning_rate_min = 1e-3     
    learning_rate_max = 1e-1 
    max_depth_min = 3     
    max_depth_max = 15     
    max_depth_step = 1 
    subsample_min = 0.8     
    subsample_max = 1     
    colsample_bytree_min = 0.8     
    colsample_bytree_max = 1     
    verbose_eval = True
    seed = 32 
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

