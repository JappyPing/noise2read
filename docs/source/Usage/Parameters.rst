Parameters
----------

Noise2read is a command-line interface (CLI) based tool to eliminate PCR and sequencing errors for short reads. It utilises CLI mode and INI file for configuring the parameters. Noise2read was mainly developed to correct short-read sequencing data, but it also provides several other modules. Therefore, to run noise2read, we must specify the module name from ["simplify_correction", "correction", "amplicon_correction", "mimic_umi", "real_umi", "umi_correction", "simulation", "evaluation"] first. Then, we set the relevant parameters required by the each module, otherwise noise2read will use the default parameters.

1. Guidance on setting noise2read parameters:
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Noise2read has many parameters. But, most of these parameters do not necessarily need to change for different datasets. 

* We strongly suggest the users to set the ``tree_method`` as "gpu_hist" which means using GPU for training the model. If you do not have GPU resources, please use the simplified version of noise2read to do error correction because using CPU training for large datasets may require several days.

* For the large datasets, we do not suggest use a big number of multiprocessing processes (``num_workers``) for noise2Read to run, as we have observed that those situations could suddenly consume a significant amount of memory, and the program ran out of memory and terminated. 

    - If noise2read runs out of memory and terminates during searching 1nt- or 2nt-edit-distance edges, please increase the parameter ``reads_chunks_num`` (default 1) and decrease ``num_workers`` to try again. Beware, a bigger ``reads_chunks_num`` may slow down the noise2read. 
    
    - If noise2read runs out of memory and terminates during training, please increase the ``chunks_num`` (default 100) and decrease ``num_workers`` to try again.

2. The parameters required by different modules are summarised as follows:
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* simplify_correction

``result_dir``, ``input_file``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``, ``iso_change_detail``, ``top_n``, ``min_read_len``

``high_freq_thre``, ``max_error_freq``, ``save_graph``, ``graph_visualization``, ``drawing_graph_num``

* correction

All the parameters of "Paths", "SourceInputData", "General", "GraphSetup", "EmbeddingSetup", "AmbiguousSetup" and "ModelTuningSetup"

.. ``result_dir``, ``input_file``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``, ``iso_change_detail``, ``top_n``, ``min_read_len``, ``negative_sample_num``

.. ``high_freq_thre``, ``max_error_freq``, ``save_graph``, ``graph_visualization``, ``drawing_graph_num``
* amplicon_correction

All the parameters of "Paths", "SourceInputData", "General", "GraphSetup", "EmbeddingSetup", "AmbiguousSetup", "ModelTuningSetup" and "Amplicon".

* umi_correction

``result_dir``, ``input_file``, ``ground_truth_data``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``, ``top_n``

``high_freq_thre``, ``max_error_freq``, ``save_graph``, ``graph_visualization``, ``drawing_graph_num``

* real_umi

``result_dir``, ``input_file``, ``num_workers``

All the parameters of "real_umi" 

* mimic_umi

``result_dir``, ``input_file``, ``ground_truth_data``

* simulation

All the parameters required by "simplify_correction" and those parameters of "simulation"

* evaluation

``result_dir``, ``input_file``, ``ground_truth_data``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``

3. Parameters Description
<<<<<<<<<<<<<<<<<<<<<<<<<

*****
Paths
*****

* ``result_dir`` [default= ``working_dirctory/result/``]

    - the directory where all the outputs saved.

***************
SourceInputData
***************

* ``input_file``

    - the path to the input data required by noise2read

* ``ground_truth_data`` [default= ``None`` ]

    - the path to the ground truth data required by noise2read evaluation module

*******
General
*******

* ``num_workers`` [default= ``-1`` ]

    - num_workers is the number of worker processes to use. If num_workers is -1 then the number returned by os.cpu_count() is used.

* ``chunks_num`` [default= ``100`` ]

    - chunks_num is the number of worker processes to use. If num_workers is -1 then the number returned by os.cpu_count() is used.

* ``reads_chunks_num`` [default= ``1`` ]

    - reads_chunks_num is used to divide list of reads into chunks when searching 1nt- or 2nt-edit-distance edges using multiprocessing. It also used by the multiprocessing module during evaluation process. Beware, if 

* ``verbose`` [default= ``False`` ]

    - If true, noise2read will save the extracted training instances such as genuine, ambiguous errors and negative reads as csv to ``result_dir``.     

* ``iso_change_detail`` [default= ``False`` ]

    - If true, noise2read will save the frquency changing reads of isolated nodes before and after correction as .txt file to ``result_dir``.     

* ``top_n`` [default= ``10`` ]

    - During the evaluation process, noise2read saves the frequency changes of the top ``top_n`` sequence frequencies into a sheet of the output .excel file to ``result_dir``.   

* ``min_read_len`` [default= ``30`` ]

    - The threshold of the sequence's minimum length to determine whether to perform 2nt-edit-distance-based error correction.  

* ``negative_sample_num`` [default= ``300000`` ]

    - When the number of negative samples larger than preseting threshold ``negative_sample_num``, noise2read will downsample negative samples for training. 

**********
GraphSetup
**********
* ``high_freq_thre`` [default= ``4`` ]

    - The threshold of τ to determine whether a read is high-frequency or lwo-frequency. 

* ``max_error_freq`` [default= ``4`` ]

    - A read is considered as an error read when its frequency is smaller than the predefined maximum frequency threshold ``max_error_freq``. 

* ``save_graph`` [default= ``False`` ]

    - If true, noise2read will save the construted graph as 'graph.gexf' to ``result_dir``.    

* ``graph_visualization`` [default= ``False`` ]

    - If true, noise2read will visualize the connected subgraphs as ".svg" and save them to ``result_dir``. 

* ``drawing_graph_num`` [default= ``50`` ]

    - The number of the connected subgraphs to be drawed.

**************
EmbeddingSetup
**************
* ``entropy_kmer`` [default= ``3`` ]

    - The kmer size used to calculate Shannon's and Tsallis's entropy during feature extraction. 

* ``entropy_q`` [default= ``2`` ]

    - The real number q used to calculate Tsallis's entropy during feature extraction.

* ``kmer_freq`` [default= ``3`` ]

    - The kmer frequency used to calculate Shannon's and Tsallis's entropy during feature extraction.    

* ``read_type`` [default= ``DNA`` ]

    - The read type of the sequencing data required to be corrected. Set to ``RNA`` when the nitrogenous base Uracil(U) instead of  Thymine (T) appears in reads of RNA sequencing data. 

**************
AmbiguousSetup
**************
* ``high_ambiguous`` [default= ``True`` ]

    - If Ture, noise2read will correct the potential ambiguous errors between high-frequency reads. 

* ``proba_deviation`` [default= ``0.95`` ]

    - The mutation observed in high-frequency reads exhibits a bidirectional nature.Therefore, we only consider the prediction result with a higher probability when the bidirectional predictions match. In other words, if the absolute difference between the probabilities of the two-way predictions is less than a specific value ``proba_deviation``, we discard the prediction; otherwise, we choose the prediction having a higher probability.

* ``iso_neg_high`` [default= ``False`` ]

    - If True, the high frequency isolated nodes aslso included as negative samples for high ambiguous prediction. This will rquire quite a lot computational resources (memory) for embeeding and model training.  

****************
ModelTuningSetup
****************
* ``n_trials`` [default= ``20`` ]

    - An Optuna trial is a process of evaluating an objective function. ``n_trials`` refers to the number of the trials for optimizing the best model.

* ``n_estimators`` [default= ``400`` ]

    - Number of boosting rounds.

* ``test_size`` [default= ``0.1`` ]

    - ``test_size`` represents the proportion of the dataset to serve as independent test for evaluating the models.

* ``random_state`` [default= ``42`` ]

    - Controls the shuffling applied to the data before applying the sklearn.model_selection.train_test_split.

    - The seed used by the random number generator to control the randomization of the algorithm of performing over-sampling using SMOTE.

* ``tree_method`` [default= ``auto`` ]

    - The tree construction algorithm used in XGBoost. See description in XGBoost documentation.

    - Choices: auto, exact, approx, hist, gpu_hist, this is a combination of commonly used updaters. For other updaters like refresh, set the parameter updater directly.

        - auto: Use heuristic to choose the fastest method.

        - exact: Exact greedy algorithm. Enumerates all split candidates.

        - approx: Approximate greedy algorithm using quantile sketch and gradient histogram.

        - hist: Faster histogram optimized approximate greedy algorithm.

        - gpu_hist: GPU implementation of hist algorithm.

* ``learning_rate_min`` [default= ``1e-3`` ]

    - The minimum learning rate of the setted learning rate intervel. Optuna will choose the learning rate from the predifined intervel to optimize a best XGBoost model. The learning rate is a step size shrinkage used in update to prevents overfitting. 
    - range: (0,1]

* ``learning_rate_max`` [default= ``1e-1`` ]

    - The maximum learning rate of the setted learning rate intervel. ``learning_rate_max`` > ``learning_rate_min``.
    - range: (0,1] 

* ``max_depth_min`` [default= ``3`` ]

    - The minimum of the setted maximum depth of a tree. Optuna will choose the maximum depth from the predifined intervel to optimize a best XGBoost model. XGBoost aggressively consumes memory when training a deep tree.

    - range: [0,∞]

* ``max_depth_max`` [default= ``15`` ]

    - The maximum of the setted maximum depth of a tree. ``max_depth_max` > ``max_depth_min``.

    - range: [0,∞]

* ``max_depth_step`` [default= ``1`` ]

    - The step size for choosing max_depth of tree from the intervel [max_depth_min, max_depth_max].

* ``subsample_min`` [default= ``0.8`` ]

    - The minimum of the subsample ratio of the training instances. Optuna will choose the subsample ratio from the predifined intervel to optimize a best XGBoost model.

    - range: (0,1]

* ``subsample_max`` [default= ``1`` ]

    - The minimum of the subsample ratio of the training instances. ``subsample_max`` >  ``subsample_min``.

    - range: (0,1]

* ``colsample_bytree_min`` [default= ``0.8`` ]

    - The minimum of the subsample ratio of columns when constructing each tree. Optuna will choose the subsample ratio from the predifined intervel to optimize a best XGBoost model. Subsampling occurs once for every tree constructed.

    - range: (0,1]

* ``colsample_bytree_max`` [default= ``1`` ]

    - The maximum of the subsample ratio of columns when constructing each tree. ``colsample_bytree_max`` > ``colsample_bytree_min``.

    - range: (0,1]

* ``verbose_eval`` [default= ``False`` ]

    -  If verbose and an evaluation set is used, writes the evaluation metric measured on the validation set to stderr.

* ``xgboost_seed`` [default= ``42`` ]

    - Random number seed.

* ``optuna_seed`` [default= ``42`` ]

    - Seed for random number generator used in optuna.samplers.TPESampler. 

********
real umi
********

* ``umi_in_read`` [default= ``False`` ]

    - If true indicates that the UMI sequences are contained in the reads.

        * ``umi_start`` [default= ``0`` ]

            - When ``umi_in_read`` is true. ``umi_start`` represents the start position of the UMIs in the reads.

        * ``umi_end`` [default= ``12`` ]

            - When ``umi_in_read`` is true. ``umi_end`` represents the end position of the UMIs in the reads.

        * ``non_umi_start`` [default= ``24`` ]

            - ``non_umi_start`` represents the start position of the sequenced target fragments which does not include the other sequence such as barcode and UMIs.

        * ``group_read_number`` [default= ``10`` ]

            - The minimum number of reads in an UMI cluster to be selected for constructing UMI-based ground truth data set. 

        * ``read_edit_dif`` [default= ``2`` ]

            - The edit difference between each low-frequency read and high-frequency read in a UMI culster. If the edit distance <= ``read_edit_dif``, the low-frequency read will be retained for constructing UMI-based ground truth data.

    - If false indicates that the UMI sequences are contained in the sequence description. Then we may use two customized separators and indices to split the description and extract the UMIs.

        * ``separator1`` [default= ``_`` ]

            - The first separtor to split the sequence description. 

        * ``separator1_idx`` [default= ``2`` ]

            - The first index to get the string containing the UMI sequence from the splited string list. 

        * ``separator2`` [default= ``_`` ]

            - The second separtor to split the splited string containing the UMI sequence. 

        * ``separator2_idx`` [default= ``0`` ]

            - The second index to get the UMI sequence from the splited string list. 

* ``read_edit_dif`` [default= ``2`` ]

    - The edit difference between each low-frequency read and high-frequency read in a UMI culster. If the edit distance <= ``read_edit_dif``, the low-frequency read will be retained for constructing UMI-based ground truth data.

********
Amplicon
********

* ``amplicon_low_freq`` [default= ``50`` ]

    - The threshold to indicate a read is a low-frequency when its frequency <= ``amplicon_low_freq`` for the additional amplicon sequencing correction.

* ``amplicon_high_freq`` [default= ``1500`` ]

    - The threshold to indicate a read is a high-frequency when its frequency >= ``amplicon_low_freq`` for the additional amplicon sequencing correction.

* ``amplicon_threshold_proba`` [default= ``0.85`` ]

    - The probability threshold to determine whether potential amplicon errors mutated from its neighbouring high-frequency reads. If the prediceted probability >= ``amplicon_threshold_proba``, then noise2read retain this prediction, otherwise discard.

**********
simulation
**********

* ``min_freq`` [default= ``5`` ]

    - The predetermined threshold to filtered out low-frequency reads after correction by noise2read simplify_correction to eliminate noise for simulation. 

* ``min_read_count`` [default= ``30`` ]

    - The minimum counts of reads to select reads for constituting an error-prone subset.  Then 1 or 2 errors are randomly injecting induced within these error-prone reads according to the predefined error rates per read.

* ``error_rate1`` [default= ``0.09`` ]

    - The 1nt-based-error rate per read.

* ``error_rate2`` [default= ``0.02`` ]

    - The 2nt-based-error rate per read.

********************
Required CLI setting
********************

* Module selection

Using noise2read, you must select the module name from ["correction", "amplicon_correction", "mimic_umi", "real_umi", "umi_correction", "simulation", "evaluation"] first.

.. code-block:: console

  -m | --module module_name

* Setting configuration file or input dataset

   * configuration

   .. code-block:: console

       -c | --config config.ini

   * Input Read dataset

   .. code-block:: console

       -i | --input data.fastq

********************
Optional CLI setting
********************

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