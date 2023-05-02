simulated_miRNAs
----------------

These examples implement the results for correcting simulated miRNAs data with mimic UMIs (`D14 and D16 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/EjBTpjExiShHg0kO72fVpzABn_Krd0K61xdLlK5_03JB5A?e=5GXsg8>`_) by noise2read.

* noise2read installation
   
Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

* Clone the codes with datasets in github

.. code-block:: console

    git clone https://github.com/Jappy0/noise2read
    cd noise2read/Examples/simulated_miRNAs
* correcting D14

  * with high ambiguous errors correction and using GPU for training

  .. code-block:: console

      noise2read -m correction -c ../../config/D14.ini -a True -g gpu_hist

  * without high ambiguous errors correction and using CPU (default) for training

  .. code-block:: console

      noise2read -m correction -c ../../config/D14.ini -a False

* correcting D16

  * with high ambiguous errors correction and using GPU for training

  .. code-block:: console

      noise2read -m correction -c ../../config/D16.ini -a True -g gpu_hist

  * without high ambiguous errors correction and using CPU (default) for training

  .. code-block:: console

      noise2read -m correction -c ../../config/D16.ini -a False