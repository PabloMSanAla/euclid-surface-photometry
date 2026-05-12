.. _module1:

Module 1 – Cutouts and Masks
==============================

.. autofunction:: euclid_pipelines.euclid_module1

Scientific description
-----------------------

The first step of the pipeline constructs multi-band cutout images centred on
each galaxy and generates a contamination mask capable of preserving the very
faint outskirts of the target.

Cutout creation
^^^^^^^^^^^^^^^

For each galaxy a square cutout of side :math:`5\,R_\mathrm{Kron}` is extracted
from every available band in the MER mosaics:

* **VIS** (:math:`I_\mathrm{E}`, 0.1 arcsec pixel⁻¹)
* **NISP near-infrared**: :math:`Y_\mathrm{E}`, :math:`J_\mathrm{E}`, :math:`H_\mathrm{E}` (0.3 arcsec pixel⁻¹, resampled to 0.1 by MER)
* **External optical** (UNIONS *ugriz* in EDF-N; DES *griz* in EDF-S/F)

The large cutout size ensures that enough sky background is available for a
reliable background estimation in M3 and for accurate masking in the outer
regions.

Background model restoration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MER pipeline subtracts a large-scale background model before delivering
the science mosaics (``EUC_MER_BGSUB-MOSAIC``).  This subtraction can
over-remove flux at the positions of large or bright galaxies (:cite:t:`Urbano25`).
To mitigate this effect, the MER background model (``EUC_MER_BGMOD``) is added
back to each cutout before any further analysis.

Centre refinement
^^^^^^^^^^^^^^^^^

After constructing the VIS cutout, the galaxy centre is refined using the
**Hill Climb** algorithm of **AutoProf** (:cite:t:`Autoprof`).  The algorithm
applies iterative FFT-based gradient detection followed by parabolic fitting
and Nelder-Mead optimisation to converge on the global flux peak, ignoring
local features such as bars or HII regions.  The updated (RA, Dec) coordinates
are propagated to all subsequent modules.

Mask generation
^^^^^^^^^^^^^^^

Reliable masks are critical for measuring surface brightness profiles to low
surface brightness levels.  A multi-step hierarchical masking strategy is used:

1. **MER segmentation map** (``EUC_MER_FINAL-SEGMAP``): provides an initial
   catalogue of detected sources from **SourceXtractor++**.

2. **NoiseChisel** (:cite:t:`noisechisel`): detects all pixels with signal down
   to very low S/N on the VIS image.  A detection map is produced that serves
   as the target map for watershed filling.

3. **Watershed growing**: the MER segments are grown to fill the NoiseChisel
   detection map using the ``scikit-image.segmentation.watershed`` algorithm
   (:cite:t:`scikit-image`) applied to an LRGB image constructed from
   :math:`J_\mathrm{E}`, :math:`H_\mathrm{E}`, :math:`Y_\mathrm{E}` (colour
   channels) and :math:`I_\mathrm{E}` (luminance), stretched with the
   ``asinh`` transform.  The colour information improves deblending performance
   compared to using a single band.

4. **MTObjects deblending** (:cite:t:`mtobjects`): run on the VIS image with
   ``move_factor = 0.3`` to identify sub-structure within the galaxy footprint.

5. **Cosine similarity filter**: distinguishes foreground/background sources
   from real galaxy substructure using the multi-band colour vector.

   For each pixel at position :math:`(i,j)`, the cosine similarity to the
   mean galaxy colour :math:`\vec{g} = (g_{I_\mathrm{E}}, g_{Y_\mathrm{E}},
   g_{J_\mathrm{E}}, g_{H_\mathrm{E}})` is:

   .. math::

      \vec{S}_C = \frac{\vec{C} \cdot \vec{g}}{\|\vec{C}\|\,\|\vec{g}\|}
               = \frac{\sum_b C_b\,g_b}
                      {\sqrt{\sum_b C_b^2}\,\sqrt{\sum_b g_b^2}} \,,

   where :math:`b` runs over the four *Euclid* photometric bands and
   :math:`\vec{C}` is the per-pixel colour vector.  Sources whose median
   cosine similarity differs from the galaxy's own value by more than
   :math:`3\,\sigma_{\vec{g}}` are flagged as external contaminants and masked.

6. **Colour outlier masking**: pixels within the galaxy footprint that deviate
   by more than :math:`5\,\sigma_{\vec{g}}` in colour space (typically
   high-redshift point sources not resolved by MTObjects) are morphologically
   eroded and dilated, and regions larger than FWHM² pixels are masked.
   The central 20 % of each galaxy segment is always left unmasked.

.. figure:: ../figures/EUC_SEGMID_-612927984481008950_mask_segment.pdf
   :alt: Mask comparison for EUCL J040510.27−444836.3
   :width: 95%
   :align: center

   *Mask comparison for EUCL J040510.27−444836.3.  Top row (left to right):
   LRGB image, cosine similarity map (:math:`\vec{S}_C`), MER segmentation map,
   and watershed output.  Bottom row (left to right): NoiseChisel+Segment,
   MTObjects, cosine similarity regions, and the final combined mask.*

Output files
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - File
     - Description
   * - ``EUC_SEGMID_<id>_VIS_cutout.fits``
     - VIS cutout with MER background added (BGMOD=1 in header)
   * - ``EUC_SEGMID_<id>_<band>_cutout.fits``
     - One cutout per band
   * - ``EUC_SEGMAP_<id>_segmap.fits``
     - MER segmentation map cutout
   * - ``EUC_SEGMID_<id>_mask.fits``
     - Final contamination mask (0 = valid, non-zero = masked)

.. note::

   The function returns updated (RA, Dec) coordinates that reflect the refined
   centre and should be passed to all downstream modules.
