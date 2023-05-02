Real_UMI
--------

This module generates UMI-based ground truth datasets using UMI-based sequencing data to evaluate computational error correction algorithms.

* Usage

  * INI configuration
    
  .. code-block:: console

      [Paths]
      ResultDir = "./result/"
      ; set output directory

      [SourceInputData]
      input_file = path/to/data.fastq
      ; set your data to be corrected

      [RealUMI]
      umi_start = 0 
      ; the UMI start position in a sequence
      umi_end = 12 
      ; the UMI end position in a sequence
      non_umi_start = 24 
      ; the start position of the read

  * Run
    
  .. code-block:: console

      noise2read -m real_umi -c config.ini
