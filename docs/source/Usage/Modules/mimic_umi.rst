Mimic_UMI
---------

This module generates mimic UMIs for simulated datasets in order to evaluate using UMI tags rather sequencing IDs by the instruments.

* Parameters required by "mimic_umi"

``result_dir``, ``input_file``, ``ground_truth_data``

See parameter descriptions at `Parameters <https://noise2read.readthedocs.io/en/latest/Usage/Parameters.html>`_

* Usage

.. code-block:: console

    noise2read -m mimi_umi -i path/to/raw.fastq -t path/to/true.fastq -d ./results/ 