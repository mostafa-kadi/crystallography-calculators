# Crystallography & XRD Calculators

Nine browser-based calculators: d-spacing (cubic and general triclinic), lattice parameter (cubic), Bragg's law, XRD peak position, general unit-cell volume, interplanar angle (cubic), reciprocal lattice parameter (cubic), and Scherrer crystallite size.

## Scientific conventions
- Reciprocal lengths use the crystallographic convention without 2π.
- Scherrer calculator accepts observed and instrumental FWHM in degrees of 2θ and applies a Gaussian correction; output is nm.
- Interplanar angle formula assumes a cubic lattice.
- Unit-cell volume supports general α, β, γ cell angles.

Open `index.html` directly or deploy the repository root as a static site. No backend or external dependencies.

Version 2 adds general triclinic d-spacing and Gaussian instrumental broadening correction for Scherrer analysis. See VALIDATION.md for reference calculations and assumptions.
