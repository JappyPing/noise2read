Correction
----------
1. Download data set
<<<<<<<<<<<<<<<<<<<<



2. Configuration file (config.ini)
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Download it from **

Or create it by yourself and copy the following contents

.. code-block:: console

    [Paths]
    ResultDir = "./result/" # set output directory

    [SourceInputData]
    input_file = path/to/data.fastq # set your data to be corrected
    # ground_truth_data = path/to/data.fastq # only set when you have groundtruth data, otherwise comment it

    [General]
    num_workers = -1 # if num_workers = -1 or 0, nois2read will use all the available cpus 
    verbose = True 
    min_iters = 100
    iso_change_detail = True
    top_n = 100



1. correction
<<<<<<<<<<<<<
.. code-block:: console

    nois2read -m correction -i config.ini



