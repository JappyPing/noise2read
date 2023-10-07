Real_UMI
--------

This module generates UMI-based ground truth datasets using UMI-based sequencing data to evaluate computational error correction algorithms.

* Parameters required by "real_umi"

``result_dir``, ``input_file``, ``num_workers``

See parameter descriptions at `Parameters <https://noise2read.readthedocs.io/en/latest/Usage/Parameters.html>`_

* Usage

  * Run with INI configuration
    
  .. code-block:: console

      noise2read -m real_umi -c config.ini
