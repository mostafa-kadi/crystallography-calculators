# Crystallography & XRD Calculators

Nine browser-based calculators: d-spacing (cubic), lattice parameter (cubic), Bragg's law, XRD peak position, general unit-cell volume, interplanar angle (cubic), reciprocal lattice (cubic), and Scherrer crystallite size.

## Scientific conventions
- Reciprocal lengths use the crystallographic convention without 2π.
- Scherrer FWHM must be instrument-corrected and entered in degrees of 2θ; output is nm.
- Interplanar angle formula assumes a cubic lattice.
- Unit-cell volume supports general α, β, γ cell angles.

Open `index.html` directly or deploy the repository root as a static site. No backend or external dependencies.

Version 2 development branch: general triclinic d-spacing and profile-aware instrumental broadening correction are under development. These advanced calculations will be documented with explicit assumptions and test cases before merging.
