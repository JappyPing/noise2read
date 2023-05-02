Evaluation
----------

Correcting NGS sequencing datasets using noise2read, it will automatically produce the evaluation results. If you want to use this module to evaluate your own algorithms, 

* Using the following configuration to set your input

.. code-block:: console

    [Paths]
    ResultDir = "./result/" # set output directory

    [SourceInputData]
    input_file = path/to/data.fastq # set your data to be corrected
    # ground_truth_data = path/to/data.fastq # only set when you have groundtruth data, otherwise comment it

    [General]
    num_workers = -1 # if num_workers = -1 or 0, noise2read will use all the available cpus 
    verbose = True 
    min_iters = 100
    top_n = 100

and run 

.. code-block:: console

    noise2read -m evaluation -c config.ini


* Using the commands only 

  * without ground truth data

  .. code-block:: console

      noise2read -m evaluation -i path/to/raw.fastq -r path/to/corrected.fastq -d ./results/ 

  * with ground truth data

  .. code-block:: console

      noise2read -m evaluation -i path/to/raw.fastq -t path/to/true.fastq -r path/to/corrected.fastq -d ./results/ 
