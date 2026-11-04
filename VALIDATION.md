# Scientific validation notes

## Independent numerical reference checks

| Calculator | Reference inputs | Expected result |
|---|---|---:|
| d-Spacing | Cubic a = 4 Å; (110) | 2.828427 Å |
| Lattice parameter | Cubic d = 2 Å; (110) | 2.828427 Å |
| Bragg's law | d = 2 Å; λ = 1.5406 Å; n = 1 | θ = 22.653053° |
| XRD peak position | d = 2 Å; λ = 1.5406 Å | 2θ = 45.306106° |
| Unit-cell volume | a = b = c = 4 Å; α = β = γ = 90° | 64 Å³ |
| Interplanar angle (cubic) | (100) and (110) | 45° |
| Reciprocal parameter (cubic) | a = 4 Å | a* = 0.25 Å⁻¹ |
| Scherrer size | K = 0.9; λ = 1.5406 Å; βobs = 0.25°; βinst = 0.15°; θ = 20° | 42.270679 nm |
| General triclinic d-spacing | a = b = 4 Å; c = 5 Å; α = 90°; β = 100°; γ = 90°; (110) | 2.806695 Å |

## Scope and assumptions

- d-spacing, lattice parameter, interplanar angle and reciprocal parameter tools currently use cubic-lattice relations.
- The unit-cell volume tool supports general cell lengths and angles.
- Reciprocal parameter uses the crystallographic convention without 2π. It is a cubic reciprocal-lattice parameter, not a full reciprocal-basis matrix.
- Scherrer β is the instrument-corrected FWHM of a peak in 2θ, entered in degrees. The calculation converts β to radians and reports a coherent-domain size, not necessarily particle size.
- Bragg order n must be a positive integer; Miller indices must be integers and not all zero.

## Validation status

The default examples above were independently recalculated and matched to the application's numerical outputs. Input rejection checks cover invalid lattice lengths, zero Miller indices, non-integer indices, non-integer Bragg order and invalid unit-cell angles.

This is a numerical and formula-level check, not yet a browser/device compatibility test or a comparison against a certified diffraction reference dataset.