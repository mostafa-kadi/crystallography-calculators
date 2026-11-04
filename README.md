# Ceramic Density Analysis

An integrated, browser-based workflow for ceramic density analysis:

1. Enter indexed XRD reflections as h, k, l and d-spacing (Å).
2. Fit lattice parameters for a selected crystal system by linear least squares in reciprocal-metric space.
3. Calculate unit-cell volume and theoretical density from chemical formula and Z.
4. Calculate bulk density from sample mass and measured volume.
5. Calculate relative density automatically.

## Supported crystal systems

Cubic, tetragonal, orthorhombic, hexagonal, rhombohedral (rhombohedral axes), monoclinic (unique b), and triclinic.

## Scientific assumptions

- Least-squares fitting minimizes residuals in 1/d². Enter correctly indexed, independent reflections and inspect the residual table.
- The fit is a metric least-squares estimate, not a full-profile Rietveld refinement. It does not correct zero shift, specimen displacement, wavelength uncertainty, peak overlap or anisotropic strain.
- Formula parser supports common element symbols, decimal occupancies and parentheses. Unsupported elements/formula syntax are rejected.
- Theoretical density assumes a single phase, full occupancy, correct formula and correct Z.
- Bulk density is mass divided by user-entered external volume. Relative density = 100 × bulk density / theoretical density.
- A relative density above 100% is flagged for review; the app does not silently clip results.

Open index.html directly or deploy the repository root as a static site. No backend or external dependencies.
