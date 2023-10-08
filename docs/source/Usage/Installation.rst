============
Installation
============

Create virtual environments to run noise2read.

1. Create a Virtual Environment
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* using conda (recommended because I developed noise2read using conda environment)

.. code-block:: console

  conda create -n noise2read_env python=3.8

.. code-block:: console 

  conda activate noise2read_env

**Or**

* using python venv

.. code-block:: console

  mkdir noise2read

.. code-block:: console 

  cd noise2read

.. code-block:: console 

  python -m venv noise2read_env

.. code-block:: console 

  source noise2read_env/bin/activate

.. Note:: 
  
  make sure you have installed and currently using python >=3.8

2. Install noise2read
<<<<<<<<<<<<<<<<<<<<<

.. * Install via bioconda

.. .. code-block:: console

..    conda install -c bioconda noise2read

* Install via pip
  
   .. code-block:: console

      pip install noise2read

   * dependencies available at pip will be installed automatically

   * install dependencies only available at bioconda

      .. code-block:: console

         conda install -c bioconda seqtk bcool

* source code clone and Installation 

.. code-block:: console

   git clone https://github.com/Jappy0/noise2read.git

.. code-block:: console 

   cd noise2read

.. code-block:: console 

   pip install -e .

* Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console
  
   conda install -c conda-forge pygraphviz