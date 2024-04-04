Correction
----------

Take datasets D1-D8 as examples

1. noise2read installation  

  Please refer to `QuickStart <https://noise2read.readthedocs.io/en/latest/QuickStart.html>`_ or `Installation <https://noise2read.readthedocs.io/en/latest/Usage/Installation.html>`_.

2. Download datasets

  Download UMI-based ground truth datasets from `D1_D8 <https://studentutsedu-my.sharepoint.com/:f:/g/personal/pengyao_ping_student_uts_edu_au/Eu08Ycnf-mNOqvo9_kNesccBIekAmqNTd_ck2692R36GhQ?e=prmqsb>`_ 

3. Configuration

   * Download config files from D1_D8 of `configs <https://github.com/Jappy0/noise2read/tree/master/configs>`_

  
4. Run correction on D1

  .. code-block:: console

      noise2read -m correction -c ./configs/correction/D1_D8/D1.ini

5. Run amplicon correction on D1

   .. code-block:: console

      noise2read -m amplicon_correction -c ./configs/D1_D8/amplicon_correction/D1.ini

.. 6. Run simplify_correction on D1

..    .. code-block:: console

..       noise2read -m simplify_correction -c ./configs/D1_D8/simplify_correction/D1.ini      