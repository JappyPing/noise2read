isomiRs_SNPs
------------

* Clone the codes optional if you have already cloned

.. code-block:: console  

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies/IsomiRs_SNPs      

* Download datasets `D22_D31 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/Eo2YQ1tynupLs_CxDGh-rBYBAXw77X2tt-2fZM4Jxv1n6Q>`_ and move it to the directory of \*CaseStudies/IsomiRs_SNPs

* Run correction

.. code-block:: console

    python correction.py ./data/no_adapters ./result/ ./corrected

* Then run

.. code-block:: console

    cd isoMiRmap
    python run.py
    python analysis.py ../isomimap_result/raw ../isomimap_result/corrected ../isomimap_result

* If you want to remove adapters yourself using `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_,

  * Make sure you have installed cutadapt

  * Then run 
    
  .. code-block:: console
    
      cd CaseStudies/IsomiRs_SNPs  
      python cutadapters.py

