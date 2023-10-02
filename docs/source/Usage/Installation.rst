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

2. Installing dependency
------------------------

.. code-block:: console

   conda install -c bioconda seqtk
   conda install -c bioconda bcool

* Optional to install pygraphviz if you need the visualised read graph.

.. code-block:: console

   conda install -c conda-forge pygraphviz

3. Installing noise2read
------------------------

* Using pip
  
.. code-block:: console

   pip install noise2read

.. code-block:: console

   pip dependencies of
      biopython == 1.79
      xgboost == 1.6.1
      Xlsxwriter == 3.0.3
      tqdm == 4.66.1
      scikit-learn == 1.3
      networkx == 2.8.5
      pandas == 2.1.0
      optuna >= 3.1.1
      matplotlib == 3.5.2
      mpire >= 2.8.0
      editdistance == 0.6.2
      imbalanced-learn == 0.9.1
      seaborn >= 0.12.1
      psutil == 5.9.5

will be installed automatically

* source code clone and Installation 
.. code-block:: console

   git clone https://github.com/Jappy0/noise2read.git
   cd noise2read
   pip install -e .

4. Bioconda version
-------------------

Bioconda channel-based noise2read will be released after paper published.

5. Singularity version
----------------------

noise2read.simg will be released after paper published.