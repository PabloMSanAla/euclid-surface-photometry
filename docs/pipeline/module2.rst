.. _module2:

Module 2 – Position Angle and Ellipticity Profiles
====================================================

.. autofunction:: euclid_pipelines.euclid_module2

Scientific description
-----------------------

Module 2 measures the radial position-angle (PA) and ellipticity (:math:`\varepsilon`)
profiles of the galaxy in the VIS band.  These profiles define the master elliptical
geometry used to extract surface brightness profiles in all bands (Module 3).

Isophote fitting
^^^^^^^^^^^^^^^^^

The measurement uses the **AutoProf** implementation of the
:cite:t:`Jedrzejewski` harmonic-expansion method, specifically the
``Isophote_Fit_FFT_Robust`` function (:cite:t:`Autoprof`).  Elliptical
apertures are fitted to the VIS image at a series of semi-major axis (SMA)
positions.

Key algorithmic choices:

* **Logarithmic radial sampling**: consecutive apertures are spaced by a
  constant logarithmic increment such that :math:`r_{i+1} = 1.03\,r_i`.  This
  gives fine sampling in the high-S/N inner regions while naturally increasing
  the aperture width in the faint outer disc, where photon noise is dominant.

* **Regularisation**: the loss function includes regularisation terms that
  penalise hard jumps between adjacent isophotes with a scale factor
  ``ap_regularize_scale = 1.5``.  This produces smooth, physically meaningful
  PA and :math:`\varepsilon` profiles.

* **Fit limit**: isophote fitting stops when the residual reaches 3 % of the
  noise level (``ap_fit_limit = 0.03``).  The noise is estimated as the
  sigma-clipped standard deviation of all non-masked pixels in the cutout.

* **Mask application**: the mask produced by M1 is applied during each
  isophotal fitting step to prevent contaminating sources from biasing the
  geometric parameters.

* **Radial growth parameter**: the profile is grown at a rate of 3 % per step
  (``growth = 0.03``).

The result is a table of PA, :math:`\varepsilon`, and their uncertainties as a
function of SMA, forming the **master geometry profile** for the source.

.. figure:: ../figures/Euclid_Pipeline.pdf
   :alt: M2 panel of the pipeline illustration
   :width: 50%
   :align: center

   *M2 panel: resulting PA (upper left) and :math:`\varepsilon` (lower left)
   profiles for EUCL J040743.76−465230.1, and the corresponding elliptical
   apertures overlaid on the colour image (right, masked pixels in red).*

Output file
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - File
     - Description
   * - ``EUC_SEGMID_<id>_master-profile.csv``
     - ECSV table with columns ``fit R``, ``fit pa``, ``fit pa_err``,
       ``fit ellip``, ``fit ellip_err``, and Fourier mode amplitudes
       (``fit Fmode A2``, ``fit Fmode Phi2``, ``fit Fmode A4``,
       ``fit Fmode Phi4``) as a function of SMA.  Metadata includes the sky
       background, the centre coordinates, and the global PA and :math:`\varepsilon`.

.. note::

   The ``result`` table returned by this function is passed directly as the
   ``master`` argument to :func:`~euclid_pipelines.euclid_module3`, carrying
   both the geometry columns and the metadata dictionary.
