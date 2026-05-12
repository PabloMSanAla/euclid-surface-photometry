.. _module4:

Module 4 – Derived Photometric Parameters
==========================================

.. autofunction:: euclid_pipelines.euclid_module4

Scientific description
-----------------------

Module 4 converts the surface brightness profiles from M3 into a set of
physically meaningful integrated and structural parameters for each photometric
band.  The methodology follows the **S⁴G** (Spitzer Survey of Stellar Structure
in Galaxies; :cite:t:`Seth10`; :cite:t:`Munoz-Mateos13`) and **CS⁴G** 
(:cite:t:`Sanchez-Alarcon25`) pipelines.

Asymptotic magnitude
^^^^^^^^^^^^^^^^^^^^^

Total magnitudes are derived from the **curve of growth** (CoG), i.e. the
cumulative flux as a function of aperture radius.  Because the CoG does not
always reach a true plateau (faint flux continues to accumulate), the
asymptotic magnitude is obtained by fitting a straight line to the final,
settling portion of the CoG and extrapolating to its y-intercept
(:cite:t:`Watkins22`).

.. figure:: ../figures/Euclid_Pipeline_M4.png
   :alt: M4 panel – curve of growth
   :width: 85%
   :align: center

   *Left M4 panel: curve of growth for EUCL J040743.76−465230.1 used to derive
   the asymptotic magnitude in the* :math:`I_\mathrm{E}` *band.  Right M4
   panel: comparison between our asymptotic magnitudes and the MER catalogue
   magnitudes; our method recovers on average 0.09 mag more flux.*

Structural radii
^^^^^^^^^^^^^^^^^

The following radii are measured for each band:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Parameter
     - Definition
   * - :math:`R_\mathrm{eff}` (**reff**)
     - Half-light (effective) radius enclosing 50 % of the asymptotic flux.
   * - :math:`R_{25.5}`, :math:`R_{26.5}`, :math:`R_{28}` (**r25p5**, **r26p5**, **r28**)
     - Isophotal radii at surface brightness levels of 25.5, 26.5, and
       28 mag arcsec⁻²; also measured with inclination correction, redshift
       dimming, and both corrections combined.
   * - :math:`R_\mathrm{petro}` (**rpetro**)
     - Petrosian radius.

At each isophotal radius the covariant PA and :math:`\varepsilon` are also
recorded.

Surface brightness at the effective radius
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mean surface brightness within :math:`R_\mathrm{eff}` is:

.. math::

   \langle \mu \rangle_\mathrm{eff} = m_\mathrm{tot} + 2.5\log_{10}(2\pi R_\mathrm{eff}^2) \,,

where :math:`m_\mathrm{tot}` is the asymptotic magnitude and all quantities
are in consistent units.

Concentration indices
^^^^^^^^^^^^^^^^^^^^^

Two non-parametric morphology proxies are measured:

* :math:`C_{82} = 5\log_{10}(R_{80}/R_{20})` — ratio of the radii enclosing
  80 % and 20 % of the total flux.
* :math:`C_{31} = R_{75}/R_{25}` — ratio of the radii enclosing 75 % and 25 %.

Inclination and redshift corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Inclination-corrected radii are derived using the outer axis ratio
:math:`(b/a) = 1 - \varepsilon_\mathrm{out}`, where :math:`\varepsilon_\mathrm{out}`
is the median ellipticity at radii beyond :math:`R_{25.5}`.  The inclination
correction applied to surface brightness is:

.. math::

   \Delta \mu_\mathrm{incl} = -2.5\log_{10}(1 - \varepsilon_\mathrm{out}) \,.

A redshift dimming correction :math:`\Delta \mu_\mathrm{dim}` is also applied
when the source redshift :math:`z` is known.  Corrected isophotal radii are
provided for three cases: inclination only, dimming only, and both.

Output file
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - File
     - Description
   * - ``EUC_SEGMID_<id>_photometry.csv``
     - One-row CSV with ``object_id`` and, for each band ``<b>``\ :
       ``mag_<b>``, ``reff_<b>``, ``sb_eff_<b>``, ``c82_<b>``, ``c31_<b>``,
       ``rpetro_<b>``, ``r25p5_<b>``, ``r26p5_<b>``, ``r28_<b>`` (plus
       ``pa`` and ``eps`` at each isophotal level and their corrected variants),
       and ``axisRatio_<b>``.
