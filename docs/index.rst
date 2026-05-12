.. Euclid Surface Photometry Pipeline documentation master file

Euclid Surface Photometry Pipeline
===================================

.. image:: figures/Euclid_Pipeline.png
   :alt: Overview of the five-module pipeline applied to galaxy EUCL J040743.76−465230.1
   :width: 100%
   :align: center

|

**Version:** 0.1.2  |  **Paper:** Sánchez-Alarcón et al. (2026, A&A) — *Euclid Quick Data Release (Q1): Disc breaks highlight the structural evolution of galaxies through cosmic time*

Overview
--------

This package implements an automatic surface brightness profile analysis pipeline
designed to extract and characterise the radial light distribution of disc galaxies
observed by the *Euclid* space mission.
In *Sánchez-Alarcón et al. 2026* the pipeline was applied to **8 748 disc galaxies** selected from the *Euclid* Quick
Release 1 (Q1) dataset, covering 63.1 deg² across the three *Euclid* Deep Fields (EDF-N,
EDF-S, EDF-Fornax).

The pipeline produces, for each galaxy:

* **Module 1**: Multi-band image **cutouts** centred on the source and combined masks that remove
  contaminating neighbours while preserving the low surface brightness outskirts of
  the target.
* **Module 2**: Radial **position-angle (PA) and ellipticity (ε) profiles** measured via elliptical
  isophote fitting.
* **Module 3**:  Per-band **surface brightness profiles** calibrated to a common sky background measured
  adaptively at the edge of each galaxy.
* **Module 4**: Integrated **photometric and structural parameters** — asymptotic magnitudes, effective
  radii, isophotal radii, concentration indices, and more.
* **Module 5**: Automatic **disc break classification** into Type I (pure exponential), Type II
  (down-bending), Type III (up-bending), and composite (e.g. Type II+III) profiles,
  together with best-fit structural parameters for each disc component.

The pipeline is parallelised with Python's ``multiprocessing.Pool`` and is driven by a
YAML configuration file. It is also implemented on the scientific platform of ESA, `Datalabs`_ see :doc:`usage` for details.

.. _Datalabs: https://datalabs.esa.int/

.. toctree::
   :maxdepth: 2
   :caption: Contents

   science
   usage
   pipeline/index
   mer
   api/index
   changelog

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
