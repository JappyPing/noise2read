QuickStart
----------

============
Easy Install
============

`noise2read <https://pypi.org/project/noise2read/>`_ is currently distributed on `PyPI <https://pypi.org/project/noise2read/>`_ and `bioconda <https://anaconda.org/bioconda/noise2read>`_, quick install noise2read in the build environment.

* Install via bioconda

.. code-block:: console

    conda install -c bioconda noise2read

.. code-block:: console
  
   conda install py-xgboost-gpu
   
.. Note:: 
  
  Currently, noise2read at bioconda requires installing py-xgboost-gpu separately. I will include this dependency in later release.

* Install via pip

.. code-block:: console

    pip install noise2read

and then install bioconda distributed packages of seqtk and bcool.

.. code-block:: console

    conda install -c bioconda seqtk bcool

Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console

    conda install -c conda-forge pygraphviz

========
Examples
========

.. #. A simplified version of noise2read which excludes machine learning instead uses heuristics for error correction for short reads set with default parameters.

..    .. code-block:: console

..        noise2read -m simplify_correction -i *.fa/fasta/fastq/fq -d output_directory

#. General correction for short reads set with default parameters.
   
   * Training with CPU
     
   .. code-block:: console

       noise2read -m correction -i *.fa/fasta/fastq/fq -a True -d output_directory

   * Training with GPU

   .. code-block:: console

       noise2read -m correction -i *.fa/fasta/fastq/fq -a True -g gpu_hist -d output_directory

#. Correcting amplicon sequencing data with default parameters

    * Training with CPU
    
    .. code-block:: console

        noise2read -m amplicon_correction -i *.fa/fasta/fastq/fq -a True -d output_directory

    * Training with GPU
    
    .. code-block:: console
        
        noise2read -m amplicon_correction -i *.fa/fasta/fastq/fq -a True -g gpu_hist -d output_directory

.. Note:: 
  
  We strongly recommend utilizing GPU for model training and prediction, especially for large data sets, rather than using a CPU. If a GPU resource is available. 
..   ; otherwise, using the simplified version of noise2read (simplify_correction) is better.