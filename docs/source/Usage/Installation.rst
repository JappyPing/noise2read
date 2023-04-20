============
Installation
============

To use noise2read, please create anaconda virtual environment (because current noise2read depend on two Conda packages and several pip packages), and then install noise2read using pip. Besides, when installing noise2read, it will automatically install all the pip dependencies. 

1. pip and conda
<<<<<<<<<<<<<<<<

1.1 Creating Conda Environment
------------------------------
.. code-block:: console

   (base) $ conda create -n noise2read python=3.8

1.2 Installing conda distributed dependency
-------------------------------------------

.. code-block:: console

   (noise2read) $ conda install -c bioconda seqtk
   (noise2read) $ conda install -c bioconda bcool

Optional to install pygraphviz if you don't need the visualised read graph.

.. code-block:: console

   conda install -c conda-forge pygraphviz

1.3 Installing noise2read
--------------------

.. code-block:: console

   (noise2read) $ pip install noise2read

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

will be installed automatically

2. source code clone
<<<<<<<<<<<<<<<<<<<<
2.1 Creating Virtual Environment 
--------------------------------
conda
-----
.. code-block:: console

   (base) $ conda create -n noise2read python=3.8

pip
---
.. code-block:: console

   (base) $ conda create -n noise2read python=3.8

-- code-block:: console

   git clone https://github.com/Jappy0/noise2read.git


1. Bioconda version
<<<<<<<<<<<<<<<<<<<

Bioconda channel-based noise2read version will be released after published.