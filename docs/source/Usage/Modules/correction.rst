Correction
----------

* Parameters required by "correction"

All the parameters of "Paths", "SourceInputData", "General", "GraphSetup", "EmbeddingSetup", "AmbiguousSetup" and "ModelTuningSetup"

* Run with ini configuration
   
.. code-block:: console

    noise2read -m correction -c config.ini

**Or**

* Run the commands only 

  * Evaluating without ground truth data

  .. code-block:: console

      noise2read -m correction -i path/to/raw.fastq -d ./results/

  * Evaluating with ground truth data

  .. code-block:: console

      noise2read -m correction -i path/to/raw.fastq -t path/to/true.fastq -d ./results/ 

  * Training with GPU (default CPU)
    
  .. code-block:: console

      noise2read -m correction -i path/to/raw.fastq -d ./results/ -g gpu_hist