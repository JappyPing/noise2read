============
Installation
============

Create virtual environments to run noise2read.

1. Creating Virtual Environment
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* using conda (recommended because I developed noise2read using conda environment )

.. code-block:: console

  conda create -n noise2read_env python=3.11
  conda activate noise2read_env

**Or**

* using python venv

.. code-block:: console

  mkdir noise2read
  cd noise2read
  python -m venv noise2read_env
  source noise2read_env/bin/activate

.. Note:: 
  
  make sure you have installed and currently using python 3

2. Installing noise2read
<<<<<<<<<<<<<<<<<<<<<<<<

* Using bioconda

.. code-block:: console

   conda install -c bioconda noise2read


* Using pip
  
   .. code-block:: console

      pip install noise2read

   * dependencies available at pip will be installed automatically

   * install dependencies only available at bioconda

      .. code-block:: console

         conda install -c bioconda seqtk bcool

* source code clone and Installation 

.. code-block:: console

   git clone https://github.com/Jappy0/noise2read.git
   cd noise2read
   pip install -e .

* Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console

   conda install -c conda-forge pygraphviz