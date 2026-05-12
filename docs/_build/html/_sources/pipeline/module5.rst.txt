.. _module5:

Module 5 – Disc Break Classification
======================================

.. autofunction:: euclid_pipelines.euclid_module5

Scientific description
-----------------------

Module 5 classifies each galaxy's surface brightness profile into one of
several disc break types and fits parametric models to each disc component.
The method is adapted from the change-point analysis of
:cite:t:`Watkins19` and :cite:t:`Sanchez-Alarcon23`, extended with automated
piecewise fitting and Bayesian model selection.

Profile preparation
^^^^^^^^^^^^^^^^^^^^

The VIS intensity profile from M3 is converted to surface brightness
magnitudes:

.. math::

   \mu_{I_\mathrm{E}} = \mathrm{ZP}_{I_\mathrm{E}} - 2.5\log_{10}(I_{I_\mathrm{E}})
                       + 5\log_{10}(\mathrm{scale}) + c \,,

where ``scale`` is the pixel scale in arcsec pixel⁻¹ and the correction
:math:`c = c_1 + c_2` includes:

**Redshift surface brightness dimming** (:cite:t:`Ribeiro2016`):

.. math::

   c_1 = -7.5\log_{10}(1+z) \,.

**Inclination correction** (following :cite:t:`Sanchez-Alarcon25`):

.. math::

   c_2 = -2.5\log_{10}(1-\varepsilon) \,,

where :math:`\varepsilon` is the median ellipticity at :math:`\mu > 25.5`
mag arcsec⁻².

A **surface brightness threshold** of :math:`\mu_\mathrm{thresh} = 27.5`
mag arcsec⁻² in the galaxy rest frame (configurable via ``sb_threshold``) is
applied to ensure a homogeneous depth across the redshift range
:math:`0 < z < 1`.  This level corresponds approximately to the *Euclid*
sensitivity limit of ~29.5 mag arcsec⁻² observed at :math:`z\approx 1`.

The profile is additionally smoothed with a Savitzky–Golay filter (width
adapted to 4 % of :math:`R_\mathrm{max}` with a minimum of 5 points,
polynomial order 1), preserving the inner region (:math:`< 1.5` arcsec) in its
unsmoothed form.

Change-point analysis
^^^^^^^^^^^^^^^^^^^^^^

The derivative of the smoothed profile, :math:`\mu^\prime(r)`, is computed via
local linear regression over the five adjacent profile points.  The
weighted mean derivative :math:`\bar{\mu}^\prime` is defined using arc-length
weights to account for the logarithmic radial spacing:

.. math::

   \bar{\mu}^\prime = \frac{\sum_j \mu^\prime(r_j)\,\Delta r_j}{\sum_j \Delta r_j} \,.

The **cumulative sum** (CS) of :math:`(\mu^\prime - \bar{\mu}^\prime)` is:

.. math::

   \mathrm{CS}_0 = 0, \quad
   \mathrm{CS}_N = \sum_{j=1}^N \left[\mu^\prime(j) - \bar{\mu}^\prime\right] \,.

Peaks and minima of CS are detected with ``scipy.signal.find_peaks``
(minimum prominence 0.15 mag arcsec⁻³) to locate candidate break positions.
:math:`N_\mathrm{peaks}` breaks divide the profile into
:math:`N_\mathrm{peaks} + 1` regions.

To isolate the outer disc, a second CS is computed after masking the inner
region (defined by the innermost peak), using only the outermost portion of
:math:`\mu^\prime` to determine the weighted mean derivative.

Piecewise model fitting
^^^^^^^^^^^^^^^^^^^^^^^^

**Inner region (Sérsic model):**
The inner region — from the centre to 110 % of the first CS peak — is always
fitted with a Sérsic profile in mag arcsec⁻² units:

.. math::

   \mu(r) = \mu_\mathrm{e} + \frac{2.5}{\ln 10}\,b_n
            \left[\left(\frac{r}{R_\mathrm{e}}\right)^{1/n} - 1\right] \,,

where :math:`b_n` is the Sérsic shape constant, and the free parameters are
:math:`(\mu_\mathrm{e},\, R_\mathrm{e},\, n)`.

**Outer regions (broken or single exponential):**
For each outer region, two models are considered:

*Broken exponential* (:cite:t:`2016A&A...587A..70S`):

.. math::

   \mu(r) =
   \left(\mu_\mathrm{in} + \frac{2.5}{\ln 10}\,\frac{r}{h_\mathrm{in}}\right)(1-W_b)
   + \left(\mu_\mathrm{out}+\frac{2.5}{\ln 10}\,\frac{r}{h_\mathrm{out}}\right)W_b \,,

with the softening weight:

.. math::

   W_b = \frac{\pi/2 + \arctan\!\left(\frac{r - r_b}{\beta}\right)}{\pi} \,,

and continuity condition at :math:`r_b`:

.. math::

   \mu_\mathrm{in} + \frac{2.5}{\ln 10}\,\frac{r_b}{h_\mathrm{in}}
   = \mu_\mathrm{out} + \frac{2.5}{\ln 10}\,\frac{r_b}{h_\mathrm{out}} \,.

The sharpness parameter is fixed at :math:`\beta = 10^{-8}` arcsec.  The loss
function used to fit the broken exponential adds two regularisation terms to
the :math:`\chi^2` to prevent degenerate solutions (identical inner and outer
scale lengths) and breaks too close to the region boundaries:

.. math::

   \ell = \chi^2
          + \frac{\alpha}{(h_\mathrm{in}-h_\mathrm{out})^2+\delta}
          + \left[\frac{\gamma}{1+e^{3(r_b-r_\mathrm{min})}}
                 +\frac{\gamma}{1+e^{2(r_\mathrm{max}-r_b)}}\right] \,,

with :math:`\alpha=10^{-2}`, :math:`\delta=10^{-4}`, and :math:`\gamma = N`
(number of data points in the region).

*Single exponential*: the first term of the broken-exponential function with
:math:`W_b = 1` everywhere.

**Model selection per region:**
The broken exponential is adopted only if both conditions are met:

a) :math:`|h_\mathrm{in} - h_\mathrm{out}|/\langle h \rangle > 0.1` —
   the two scale lengths differ by at least 10 %  (:cite:t:`Wang18`;
   :cite:t:`Yu25`).
b) :math:`\mathrm{BIC}_\mathrm{broken} < \mathrm{BIC}_\mathrm{single}`, where

   .. math::

      \mathrm{BIC} = N\ln\!\left(\chi^2/N\right) + \ln(N)\,N_\mathrm{varys} \,.

Global model and final classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The per-region components are assembled into a **combined model** (at most one
Sérsic + up to four disc components, with continuity enforced between
consecutive discs).  The combined model is re-fitted to the full profile with
all free parameters initialised from the per-region fits.

The combined model is accepted over a single global Sérsic fit only if its
BIC is lower.  Profiles requiring more than four disc components are flagged as
**Null**.

The resulting scale lengths :math:`(h_1, h_2, h_3)` — ordered from the galaxy
centre outward — determine the break type:

.. list-table:: Break type classification
   :header-rows: 1
   :widths: 15 10 75

   * - Type
     - # discs
     - Criterion
   * - 0 (Sérsic only)
     - 0
     - Single Sérsic profile; no exponential disc component
   * - **I**
     - 1
     - Single disc; no break
   * - **II**
     - 2
     - :math:`h_1 > h_2` — down-bending (truncated)
   * - **III**
     - 2
     - :math:`h_1 < h_2` — up-bending (anti-truncated)
   * - **II+II**
     - 3
     - :math:`h_1 > h_2` and :math:`h_2 > h_3`
   * - **III+III**
     - 3
     - :math:`h_1 < h_2` and :math:`h_2 < h_3`
   * - **II+III**
     - 3
     - :math:`h_1 > h_2` and :math:`h_2 < h_3`
   * - **III+II**
     - 3
     - :math:`h_1 < h_2` and :math:`h_2 > h_3`
   * - Null
     - >3
     - More than four components required

Estimated classification accuracy is **~70 %** based on visual inspection of
490 randomly selected galaxies from the reliable sample
(:math:`\chi^2_\nu < 0.75`, :math:`\chi^2_\mathrm{max} < 7.4`).

.. figure:: ../figures/Euclid_Pipeline_M5.png
   :alt: M5 panel – change-point analysis and profile fit
   :width: 100%
   :align: center

   *M5 panels for EUCL J040743.76−465230.1.  Top left: :math:`\mu_{I_\mathrm{E}}`
   profile and its derivative :math:`\mu^\prime`.  Bottom left: cumulative sum
   CS with the single detected peak (after masking the inner region) used to
   define two disc regions.  Bottom right: best-fit broken exponential (Type II)
   overlaid on the observed profile.*

.. figure:: ../figures/average_profiles_VIS.png
   :alt: Stacked profiles by disc break type
   :width: 90%
   :align: center

   *Stacked normalised SB profiles for Type I, II, III, and II+III galaxies.
   Profiles are scaled to the break radius and the surface brightness at the
   break.  For Type I, the scale length and central SB are used.*

Output file
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - Column group
     - Description
   * - ``sersic_*``
     - Global Sérsic fit parameters (:math:`\mu_e`, :math:`R_e`, :math:`n`) and χ²
   * - ``sersic_inner_*``
     - Inner-region Sérsic fit parameters and radial limits
   * - ``disc[1|2]_*``
     - Per-region exponential/broken-exponential fit results (``broken_flag``,
       :math:`\mu_0`, :math:`h_0`, :math:`r_\mathrm{break}`, :math:`h_1`,
       BIC, BIC ratio)
   * - ``final_*``
     - Combined model parameters: :math:`\mu_e`, :math:`R_e`, :math:`n`,
       :math:`\mu_0`, :math:`h_1 \ldots h_4`,
       :math:`r_\mathrm{break,1} \ldots r_\mathrm{break,3}`,
       number of disc components, χ², BIC, BIC ratio
   * - ``final_classification``
     - Integer break type (−1 = not classified; corresponds to types in the
       table above)
   * - ``sb_threshold``, ``rad_cut``, ``sb_cut``
     - Applied SB threshold and corresponding profile truncation radius/SB
