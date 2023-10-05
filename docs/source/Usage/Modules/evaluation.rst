Evaluation
----------

Correcting NGS sequencing datasets using noise2read, it will automatically produce the evaluation results. If you want to use this module to evaluate your own algorithms, 

* Parameters required by "evaluation"

``result_dir``, ``input_file``, ``ground_truth_data``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``

* Run with INI configuration

.. code-block:: console

    noise2read -m evaluation -c config.ini


* Using the commands only 

  * without ground truth data

  .. code-block:: console

      noise2read -m evaluation -i path/to/raw.fastq -r path/to/corrected.fastq -d ./results/ 

  * with ground truth data

  .. code-block:: console

      noise2read -m evaluation -i path/to/raw.fastq -t path/to/true.fastq -r path/to/corrected.fastq -d ./results/ 
