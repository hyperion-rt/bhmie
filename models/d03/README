About
=====

This directory contains code to compute the dust absorption and scattering
properties for the models in Draine (2003b,c). Note that these are not
perfectly reproduced since there is no special treatment for small grains
here.

The published properties are taken from:

    http://www.astro.princeton.edu/~draine/dust/dustmix.html

The description from the website reads:

    Draine (2003a) discussed the available evidence for dust extinction per
    unit H column in regions with different extinction curves. Based on this
    discussion, we have renormalized the size distributions of Weingartner &
    Draine (2001). For R_V=3.1, we have reduced the grain abundance per H by a
    factor 0.93 For R_V=4.0, we have taken the "A" model of Weingartner and
    Draine (2001) and increased the grain abundance per H by a factor
    0.93*1.27=1.18. For R_V=5.5, we have taken the R_V=5.5 "A" model of
    Weingartner & Draine (2001) and increased the grain abundance per H by a
    factor 0.93*1.52=1.42

    Dielectric functions were recently revised by Draine (2003b,c), including
    realistic structure at X-ray absorption edges of C, O, Mg, Si, and Fe.

    As described in Draine (2003a) and Draine (2003b), extinction, absorption,
    albedo, <cos(theta)>, and <cos^2(theta)> have been calculated for
    wavelengths from 1 cm (30 GHz) to 1 Angstrom (12.4 keV), for selected
    mixtures of carbonaceous grains and amorphous silicate grains:

    * Milky Way, R_V = 3.1: Weingartner & Draine (2001) Milky Way size
      distribution for R_V=3.1 with C/H = b_C = 60 ppm in log-normal size
      dists, but renormalized by a factor 0.93 (now has C/H= 60*0.93 = 55.8
      ppm in log-normal size distributions). This grain model is considered to
      be appropriate for the typical diffuse HI cloud in the Milky Way.

    * Milky Way, R_V = 4.0: Weingartner & Draine (2001) Milky Way size
      distribution "A" for R_V=4.0 with C/H = b_C = 40 ppm in log-normal size
      dists, renormalized by a factor 1.18 (now has C/H = 40*1.18 ppm = 47.2
      ppm in log-normal size distributions). Milky Way, R_V = 5.5: Weingartner
      & Draine (2001)

    * Milky Way size distribution "A" for R_V=5.5 with C/H = b_C = 30 ppm in
      log-normal size dists, renormalized by a factor 1.42 (now has C/H =
      30*1.42 = 42.6 ppm in log-normal size distribution).

These three models are the ones included in this directory (*.orig) and the
*.in files are parameter files to reproduce them with the bhmie code.

References
==========

  - Draine, B.T. 2003a, ARA&AA 41, 241-289
  - Draine, B.T. 2003b, ApJ, 598, 1017-1025
  - Draine, B.T. 2003c, ApJ, 598, 1026-1037
  - Weingartner, J.C., & Draine, B.T. 2001, ApJ, 548, 296-309

Running
=======

To generate the dust properties:

    bhmie d03_3.1_6.0_A.in
    bhmie d03_4.0_4.0_A.in
    bhmie d03_5.5_3.0_A.in

To compare to the original properties by B. Draine:

    python compare.py
