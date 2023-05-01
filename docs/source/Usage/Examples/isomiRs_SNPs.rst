isomiRs_SNPs
------------

* Clone the codes optional if you have already cloned

.. code-block:: console  

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/CaseStudies/IsomiRs_SNPs      

* Download datasets `D22_D31 <https://studentutsedu-my.sharepoint.com/personal/pengyao_ping_student_uts_edu_au/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpengyao%5Fping%5Fstudent%5Futs%5Fedu%5Fau%2FDocuments%2Fnoise2read%5Fdata%2FD22%5FD31&view=0>`_ and move it to the directory of \*CaseStudies/IsomiRs_SNPs

* Then run

.. code-block:: console

    cd isoMiRmap
    python run.py

* If you want to remove adapters yourself using `cutadapt <https://cutadapt.readthedocs.io/en/stable/>`_,

  * Make sure you have installed cutadapt

  * Then run 
    
  .. code-block:: console
    
      cd CaseStudies/IsomiRs_SNPs  
      python cutadapters.py

