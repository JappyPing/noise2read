Simulation
----------

This module generates simulated datasets with mimic UMIs.

* Parameters required by "simulation"

``result_dir``, ``input_file``, ``num_workers``, ``chunks_num``, ``reads_chunks_num``, ``iso_change_detail``, ``top_n``, ``min_read_len``

``high_freq_thre``, ``max_error_freq``, ``save_graph``, ``graph_visualization``, ``drawing_graph_num``

See parameter descriptions at `Parameters <https://noise2read.readthedocs.io/en/latest/Usage/Parameters.html>`_

* Usage

.. code-block:: console

    noise2read -m simulation -i path/to/raw.fastq -d ./results/ 
