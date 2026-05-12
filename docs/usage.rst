.. _usage:

Usage
=====

Requirements
------------

Install the required Python packages:

.. code-block:: bash

   pip install numpy pandas astropy scipy scikit-image photutils tqdm pyyaml fabada

The pipeline also depends on the ``astropipe`` package and external tools
**NoiseChisel**/**Segment** (part of `GNU Astronomy Utilities <https://www.gnu.org/software/gnuastro/>`_)
and **MTObjects** for source detection within each cutout.

Configuration file
------------------

All pipeline options are controlled by a YAML configuration file.
A minimal example (``euclid_config.yaml``):

.. code-block:: yaml

   absolute_paths:
     data_path:   /path/to/MER/mosaics      # Root directory of Euclid MER tiles
     output_path: /path/to/output           # Results are written here

   catalogs:
     objects: /path/to/master_catalog.csv   # MER object catalogue with all columns
     ids:     /path/to/object_ids.csv       # CSV file with an 'object_id' column

   parameters:
     n_radius:     5       # Cutout size as a multiple of the Kron radius (default 5)
     sb_threshold: 27.5    # Rest-frame SB limit applied in M5 (mag arcsec⁻²)
     max_size:     1000    # Hard cap on cutout size (pixels, default 1000)

   pipeline_steps:
     M1: true   # Cutouts and masks
     M2: true   # PA and ellipticity profiles
     M3: true   # Sky background and surface brightness profiles
     M4: true   # Derived photometric parameters
     M5: true   # Disc break classification

   logging:
     # Minimum severity to emit: DEBUG | INFO | WARNING | ERROR | CRITICAL
     level: INFO
     # Path for the rotating log file; set to null to disable file logging
     file: pipeline.log
     # strftime-compatible format string
     format: "%(asctime)s [%(levelname)-8s] %(processName)s - %(message)s"
     datefmt: "%Y-%m-%d %H:%M:%S"

Setting any pipeline step to ``false`` instructs the pipeline to read
previously saved results from the output directory instead of recomputing them.

Logging levels
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Level
     - Emitted messages
   * - ``DEBUG``
     - Per-galaxy module entry/exit, file paths, radial-bin counts, break
       classifications, full tracebacks on errors.
   * - ``INFO``
     - Pipeline start/finish, galaxy completion status (M1–M5 flags), final
       result file path and execution time.  **Recommended for production runs.**
   * - ``WARNING``
     - Non-fatal issues that do not cause a module to fail (e.g. fallback to
       sigma-clipped sky estimate).
   * - ``ERROR``
     - Per-galaxy pipeline errors (exception message only at this level).
   * - ``CRITICAL``
     - Fatal errors that stop the entire pipeline.

File logging uses a :class:`logging.handlers.RotatingFileHandler` (10 MB
per file, 5 backup files) so log files never grow unbounded.  Set ``file`` to
``null`` to log to the console only.

Running the pipeline
--------------------

Single process
^^^^^^^^^^^^^^

.. code-block:: python

   from euclid_profile_pipeline_pool import euclid_pipeline, read_yaml
   import pandas as pd

   config  = read_yaml('euclid_config.yaml')
   table   = pd.read_csv(config['catalogs']['ids'])
   dataset = [(oid, config) for oid in table['object_id']]

   for data in dataset:
       result = euclid_pipeline(data)
       print(result)

Parallel execution
^^^^^^^^^^^^^^^^^^

Use the command-line entry point, specifying the number of parallel workers
with ``-n``:

.. code-block:: bash

   python euclid_profile_pipeline_pool.py euclid_config.yaml -n 8

The script prints progress using ``tqdm`` and writes a timestamped results CSV
once all galaxies have been processed:

.. code-block:: text

   Running Euclid Pipeline v0.1.2
   with config file: euclid_config.yaml, 8 processors, and for 8748 galaxies.
   Starting at: 2026-05-11_10:00
   100%|████████████████| 8748/8748 [18:00<00:00, ...]
   Execution time:  0 days  18h:2m:17s
   Pipeline completed successfully. Check: profile_results_2026-05-11_10:00.csv

Output files
------------

For each galaxy (``<object_id>``), the following files are written inside
``<output_path>/source_<object_id>/``:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - File
     - Description
   * - ``EUC_SEGMID_<id>_VIS_cutout.fits``
     - VIS-band image cutout (background restored)
   * - ``EUC_SEGMID_<id>_<band>_cutout.fits``
     - One cutout per band (VIS, Y_E, J_E, H_E + external optical)
   * - ``EUC_SEGMAP_<id>_segmap.fits``
     - Cutout of the MER segmentation map
   * - ``EUC_SEGMID_<id>_mask.fits``
     - Combined contamination mask
   * - ``EUC_SEGMID_<id>_master-profile.csv``
     - VIS isophote-fit PA/ε profile (ECSV format)
   * - ``EUC_SEGMID_<id>_profiles.csv``
     - Multi-band surface brightness profiles (ECSV format)
   * - ``EUC_SEGMID_<id>_photometry.csv``
     - Integrated photometry and structural parameters (CSV format)
   * - ``EUC_SEGMID_<id>_breaks.csv``
     - Disc break classification and component parameters (CSV format)
   * - ``profile_results_<timestamp>.csv``
     - Batch summary: pipeline step status per galaxy (root directory)
