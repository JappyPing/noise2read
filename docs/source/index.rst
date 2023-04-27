.. image:: ../../logo/logo.svg
   :align: center

========================
noise2read Documentation
========================

.. image:: https://readthedocs.org/projects/noise2read/badge/?version=latest
    :target: https://noise2read.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

**noise2read**, originated in a computable rule translated from PCR erring mechanism that: a rare read is erroneous if it has a neighboring read of high abundance, turns erroneous reads into their original state without bringing up any non-existing sequences into the short read set(<300bp).

.. note::

   The easy-usable and automatic tuning of the classifiers' parameters facilitates wide-range explorations, but we note that noise2read may yield a slightly different result at different trials, even setting the same seeds.

.. toctree::
   .. :caption: Overview
   :maxdepth: 2
   :hidden:

   QuickStart

.. toctree::
   :caption: Usage
   :maxdepth: 3
   :hidden:

   Usage/Installation
   Usage/Configuration
   Usage/Modules/Index
   Usage/Examples/Index

.. toctree::
   :caption: About
   :maxdepth: 2
   :hidden:

   About/Issues
   About/Releases
   About/License
   About/Citation

.. .. toctree::
..    :caption: Applications
..    :maxdepth: 2
..    :hidden:

..    Applications/BaseCoverage
..    Applications/RNAQuantification
..    Applications/IsomiRs
..    Applications/BaseEditors

.. tidelift-referral-banner::