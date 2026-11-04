# Validation results

## Synthetic lattice-parameter regression

The reciprocal-metric fitting equations were checked against synthetic d-spacings generated from known direct-cell metrics. Each dataset includes at least the number of independent reflections required by its model.

| Crystal system | Reference cell (Å; °) | Maximum absolute parameter error |
|---|---|---:|
| Cubic | a = b = c = 4.000; α = β = γ = 90 | 0 |
| Tetragonal | a = b = 3.900, c = 4.100; α = β = γ = 90 | < 1×10⁻¹⁵ |
| Orthorhombic | a = 4.000, b = 5.000, c = 6.000; α = β = γ = 90 | < 3×10⁻¹⁵ |
| Hexagonal | a = b = 3.000, c = 5.000; α = β = 90, γ = 120 | < 2×10⁻¹⁴ |
| Rhombohedral (rhombohedral axes) | a = b = c = 4.000; α = β = γ = 75 | < 2×10⁻¹⁴ |
| Monoclinic (unique b) | a = 4.000, b = 5.000, c = 6.000; α = γ = 90, β = 105 | < 3×10⁻¹⁵ |
| Triclinic | a = 4.000, b = 5.000, c = 6.000; α = 75, β = 82, γ = 70 | < 6×10⁻¹⁴ |

Maximum calculated-versus-generated d-spacing differences were below 4×10⁻¹⁵ Å in these exact synthetic cases (floating-point numerical precision).

**Hexagonal correction:** the reciprocal metric is reconstructed with diagonal terms 4p/3 and off-diagonal terms +2p/3, where p = 1/a². This recovers γ = 120° for the conventional hexagonal cell. The previous reconstruction incorrectly returned γ = 60° and an incorrect a parameter.

## Formula and density checks

- Formula parser checks: BaTi0.96Nb0.04O3 = 234.99256 g mol⁻¹; Ca3(PO4)2 = 310.17400 g mol⁻¹; Al2(SO4)3 = 342.13200 g mol⁻¹; Ba(Zr0.8Y0.2)O2.9 = 274.48450 g mol⁻¹.
- Unsupported element (XxO2) and unmatched closing parenthesis were rejected.
- Density equation check: for BaTi0.96Nb0.04O3, Z = 1 and Vcell = 64.000 Å³, calculated theoretical density = 6.09710 g cm⁻³.
- Bulk-density check: 1.000 g / 0.200 cm³ = 5.000 g cm⁻³.

## Scope of verification

The regression checks above were run independently against the mathematical model and parser. They are not a browser/device end-to-end test, and they do not establish accuracy for experimental XRD data. Experimental use still requires correctly indexed reflections, appropriate crystal-system selection, adequate independent reflections, and comparison with a trusted crystallographic refinement/reference. The method is a least-squares fit to d-spacing data, not a Rietveld refinement.
