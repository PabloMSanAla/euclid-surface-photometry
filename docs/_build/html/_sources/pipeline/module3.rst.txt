.. _module3:

Module 3 – Sky Background and Surface Brightness Profiles
===========================================================

.. autofunction:: euclid_pipelines.euclid_module3

Scientific description
-----------------------

Module 3 performs two tasks in sequence:

1. Determine the **sky background level** and the **maximum profile radius**
   :math:`R_\mathrm{max}` for each band.
2. Extract **surface brightness profiles** in all available bands using the
   elliptical geometry from M2.

Sky background and profile extension
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sky background cannot be determined from a fixed aperture because the
relevant radial extent varies from galaxy to galaxy.  The following adaptive
method is used.

**Step 1 — Fixed elliptical intensity profile.**
The outer PA and ellipticity are fixed to the median of the outermost five
isophotes from the M2 profile.  An intensity profile :math:`I_b(r)` is
measured in elliptical annuli out to the edge of each cutout.

**Step 2 — Locating the background transition.**
For a pure exponential disc, the ratio :math:`I_b / I^\prime_b`, where
:math:`I^\prime_b = \Delta I_b / \Delta r`, is constant.  When the profile
reaches the sky background a plateau forms, :math:`I^\prime_b\rightarrow 0`,
causing :math:`I_b/I^\prime_b \rightarrow \infty`.  The transition is
identified as the first radius at which :math:`I_b/I^\prime_b` exceeds the
mean of the profile by more than :math:`3\,\sigma`.  The derivative is
estimated by a linear fit to the five adjacent profile points to suppress local
noise, and the intensity profile :math:`I_b` is pre-smoothed with a
Savitzky–Golay filter (window = 5, polynomial degree = 1; :cite:t:`SG-Filter`
via :cite:t:`scipy`) to remove high-frequency artefacts.

This radius is set as :math:`R_\mathrm{max}` for that band.

**Step 3 — Background value and uncertainty.**
The sky value for each band is the sigma-clipped median of the unmasked pixels
in an elliptical annulus of width :math:`0.1\,R_\mathrm{max}` centred at
:math:`R_\mathrm{max}`.  Square sky boxes placed at the same distance provide
the background uncertainty (standard deviation of box medians normalised by
:math:`\sqrt{N_\mathrm{boxes}}`).

The maximum :math:`R_\mathrm{max}` across all bands defines a common outer
radius out to which the master profile is optionally extended via logarithmic
interpolation before profile extraction.

.. figure:: ../figures/Euclid_Pipeline.pdf
   :alt: M3 panel – sky background estimation
   :width: 50%
   :align: center

   *M3 panels for EUCL J040743.76−465230.1.  Top left: VIS cutout visualised
   around the sky level, with the elliptical and square apertures used to
   measure the background.  Top right: intensity profile (grey dots = all
   pixels; red = sigma-clipped mean; green = Savitzky–Golay smoothed) and
   :math:`I_b/I_b^\prime` profile below.  Bottom left: final SB profiles in
   all* Euclid *bands.  Bottom right:* Euclid *colour profiles.*

Surface brightness profile extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once :math:`R_\mathrm{max}` and the sky backgrounds are known, surface
brightness profiles are extracted for every available band using the
radially varying PA and :math:`\varepsilon` from the M2 master profile.

For each elliptical annulus:

* The **surface brightness** is computed as the sigma-clipped average intensity
  of the unmasked pixels in the annulus.
* **Flux** is also measured; flux lost to masked pixels is replaced by
  assigning the annulus mean surface brightness to those pixels before
  integration.
* Formal uncertainties are propagated from the per-pixel noise.

In addition to the surface brightness and flux, the following geometric
quantities from M2 are included in the output table: PA, :math:`\varepsilon`,
Fourier mode amplitudes (A2, :math:`\Phi_2`, A4, :math:`\Phi_4`), pixel
centre coordinates, and number of contributing pixels per annulus.

Output file
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - File
     - Description
   * - ``EUC_SEGMID_<id>_profiles.csv``
     - ECSV table.  Columns include ``sma`` (semi-major axis in pixels),
       ``pa``, ``eps``, morphological Fourier mode coefficients, and for each
       band ``<b>``: ``int_<b>``, ``intstd_<b>``, ``flux_<b>``, ``fluxstd_<b>``.
       Table metadata contains the sky background (``bkg_<b>``), background
       uncertainty (``bkgstd_<b>``), background radius (``bkgrad_<b>``),
       image zero-point (``zp_<b>``), and pixel scale (``pixscale_<b>``).

.. note::

   If the primary background estimation (``background_estimation_euclid``) fails
   (e.g., when the galaxy fills the cutout), the fallback is the global
   sigma-clipped statistics of the masked image.
