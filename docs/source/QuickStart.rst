QuickStart
----------

============
Easy Install
============

`noise2read <https://pypi.org/project/noise2read/>`_ is currently distributed on `PyPI <https://pypi.org/project/noise2read/>`_, quick install noise2read in the build environment.

.. code-block:: console

    pip install noise2read

and then install three conda distributed packages of seqtk, bcool and pygraphviz, respectively.

.. code-block:: console

    conda install -c bioconda seqtk

.. code-block:: console

    conda install -c bioconda bcool

Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console

    conda install -c conda-forge pygraphviz

========
Examples
========
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
