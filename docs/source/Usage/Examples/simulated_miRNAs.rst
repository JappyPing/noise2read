simulated_miRNAs
----------------

This example implements the results for correcting simulated miRNAs data with mimic UMIs (`D14_D17 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EjBTpjExiShHg0kO72fVpzABn_Krd0K61xdLlK5_03JB5A?e=5GXsg8>`_) by nois2read.

* nois2read installation
   
Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

* Clone the codes with datasets in github

.. code-block:: console

    (noise2read_env)$ git clone https://github.com/Jappy0/noise2read
    (noise2read_env)$ cd noise2read/Examples/simulated_miRNAs

* run

* with high ambiguous errors correction

.. code-block:: console

    (noise2read_env)$ nois2read -m correction -i ../../config/simulated_miRNA.ini -a True

* without high ambiguous errors correction

.. code-block:: console

    (noise2read_env)$ nois2read -m correction -i ../../config/simulated_miRNA.ini -a False