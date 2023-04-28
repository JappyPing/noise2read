isomiRs_SNPs
------------

* Clone the codes optional if you have already cloned

.. code-block:: console  

    git clone https://github.com/Jappy0/noise2read
    cd CaseStudies/IsomiRs_SNPs      

* Download datasets `D22_D31 <>`_ and move it to the directory

* Then run

.. code-block:: console

    cd isoMiRmap
    python run.py

* If you want to remove adapters yourself using `cutadapt <>`_,

  * Make sure you have installed cutadapt

  * Then run 
    
  .. code-block:: console
    
      cd CaseStudies/IsomiRs_SNPs  
      python cutadapters.py

