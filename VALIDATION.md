# Validation plan

## Analytic checks to run before production use

- Cubic synthetic dataset: a = 4.000 Å; reflections (100), (110), (111), (200), with d = a / sqrt(h²+k²+l²). Fit must recover a = 4.000 Å.
- Tetragonal synthetic dataset: a = b = 3.900 Å, c = 4.100 Å; use independent hkl rows that separate basal and c-axis terms.
- Orthorhombic synthetic dataset: a = 4.0 Å, b = 5.0 Å, c = 6.0 Å; use (100), (010), (001) plus redundant mixed reflections.
- Hexagonal synthetic dataset: a = 3.0 Å, c = 5.0 Å; use (100), (002), (101) and additional independent reflections.
- Rhombohedral, monoclinic and triclinic datasets should be generated from known direct-cell metrics and checked for parameter recovery, positive definiteness and residuals.
- Density check: BaTiO3 formula mass with known Z and V should match ZM/(NA Vcell).
- Bulk/relative density check: m = 1.000 g, V = 0.200 cm³ gives 5.000 g cm⁻³; relative density is 100 × 5 / rho_theory.

## Important

These are the validation cases to execute; they are not a claim that every case has already passed. Browser/device testing and comparison with a trusted crystallographic reference remain required.
