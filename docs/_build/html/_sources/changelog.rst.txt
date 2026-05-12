.. _changelog:

Changelog
=========

Version 0.1.2
-------------
* Added ``sb_threshold`` parameter to M5 to apply a rest-frame surface
  brightness cut before the change-point analysis.
* Parallel execution via ``multiprocessing.Pool`` in
  ``euclid_profile_pipeline_pool.py``.
* Extended master profile in M3 to the maximum :math:`R_\mathrm{max}` across
  all bands.

Version 0.1.1
-------------
* Background estimation in M3 upgraded to the ``background_estimation_euclid``
  adaptive method; sigma-clipped fallback retained.
* Cosine-similarity filter added to M1 mask generation.
* MTObjects deblending integrated into M1.

Version 0.1.0
-------------
* Initial release.
* Five-module pipeline: cutouts, isophote fitting, multi-band profiles,
  photometry, and disc break classification.
* Applied to Euclid Q1 dataset (8 748 disc galaxies).
