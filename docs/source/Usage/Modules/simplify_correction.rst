simplify_correction
-------------------


* Parameters required by "simplify_correction"

``result_dir``, ``input_file``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``, ``iso_change_detail``, ``top_n``, ``min_read_len``

``high_freq_thre``, ``max_error_freq``, ``save_graph``, ``graph_visualization``, ``drawing_graph_num``


* Run with ini configuration.
   
.. code-block:: console

    noise2read -m simplify_correction -c config.ini

**Or**

* Run the commands only 

  * Evaluating without ground truth data

  .. code-block:: console

      noise2read -m simplify_correction -i path/to/raw.fastq -d ./results/

  * Evaluating with ground truth data

  .. code-block:: console

      noise2read -m simplify_correction -i path/to/raw.fastq -t path/to/true.fastq -d ./results/
