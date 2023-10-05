Amplicon_correction
-------------------

Correction errors for Amplicon sequencing data

* Parameters required by "amplicon_correction"

All the parameters of "Paths", "SourceInputData", "General", "GraphSetup", "EmbeddingSetup", "AmbiguousSetup", "ModelTuningSetup" and "Amplicon".

* Run with ini configuration
   
.. code-block:: console

    noise2read -m amplicon_correction -c config.ini

**Or**

* Run the commands only 

  * Evaluating without ground truth data

  .. code-block:: console

      noise2read -m amplicon_correction -i path/to/raw.fastq -d ./results/

  * Evaluating with ground truth data

  .. code-block:: console

      noise2read -m amplicon_correction -i path/to/raw.fastq -t path/to/true.fastq -d ./results/

  * Training with GPU (default CPU)
    
  .. code-block:: console

      noise2read -m amplicon_correction -i path/to/raw.fastq -d ./results/ -g gpu_hist