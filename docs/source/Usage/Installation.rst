============
Installation
============

Create virtual environments to run noise2read.

1. Creating Virtual Environment
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

* using conda (recommended)

.. code-block:: console

  (base) $ conda create -n noise2read_env python=3.8

**Or**

* using python venv

.. code-block:: console

  $ mkdir noise2read
  $ cd noise2read
  $ python -m venv noise2read_env
  $ source noise2read_env/bin/activate

.. Note:: 
  
  make sure you have installed and currently using python 3

2. Installing dependency
------------------------

.. code-block:: console

   (noise2read_env) $ conda install -c bioconda seqtk
   (noise2read_env) $ conda install -c bioconda bcool

* Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console

   conda install -c conda-forge pygraphviz

3. Installing noise2read
------------------------

* Using pip
  
.. code-block:: console

   (noise2read_env) $ pip install noise2read

.. code-block:: console

   pip dependencies of
      biopython == 1.79
      xgboost == 1.6.1
      Xlsxwriter == 3.0.3
      tqdm == 4.64.0
      scikit-learn == 1.1.1
      networkx == 2.8.5
      pandas == 1.4.3
      optuna == 2.10.1
      matplotlib == 3.5.2
      mpire >= 2.5.0
      editdistance == 0.6.0
      imbalanced-learn == 0.9.1
      seaborn >= 0.12.1

will be installed automatically

* source code clone and Installation 
.. code-block:: console

   (noise2read_env)$ git clone https://github.com/Jappy0/noise2read.git
   (noise2read_env)$ cd noise2read
   (noise2read_env)$ pip install -e .

4. Bioconda version
<<<<<<<<<<<<<<<<<<<

Bioconda channel-based noise2read will be released after paper published.

5. Singularity version
<<<<<<<<<<<<<<<<<<<<<<

noise2read.simg will be released after paper published.