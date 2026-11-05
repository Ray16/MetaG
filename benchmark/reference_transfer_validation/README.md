# Charged-reference chemical proton-transfer controls

The chemical reactions are acetic acid + ethanolamine -> acetate +
ethanolammonium, and 2-mercaptoethanol + ethanolamine ->
2-mercaptoethanolate + ethanolammonium. Both conserve elements, charge and
molecular count. NIST-evaluated ionization constants give chemical standard
energies -27.07 and +1.44 kJ/mol at 298.15 K and zero ionic strength in water.
The common molality-to-molarity change cancels in each two-to-two exchange.
The source is Goldberg, Kishore and Lennen (2002), Tables 7.2, 7.28 and 7.43:
[DOI 10.1063/1.1416902](https://doi.org/10.1063/1.1416902).
The full primary tables were read via the web tool; direct binary retrieval
was refused. The curated constants and locators are supplied, with no claim
of local source-PDF ownership. Their uncertainty codes are subjective
assessments, not standard deviations or declared confidence intervals.

All 48 planned optimizations finished: six chemical states, four fixed ETKDG
starts and both gas and effective-water optimization. Thirty geometries remain
under the frozen full-H RMSD/energy clustering rule. Water forces use
E(UMA) + E(GFN2/ALPB) - E(GFN2/gas), dense 5810-point SASA quadrature and
xTB accuracy 0.001. Gas thermochemistry uses exact UMA curvature; water
thermochemistry includes independently reconstructed solvent curvature at two
finite-difference steps. All retained states pass the declared force/mode gates.
Two acetate states have soft negative modes below 50 cm^-1: this is an approximate
mode policy, not proof of strictly positive curvature or conformational completeness.

Eight fixed arms cross geometry, electronic level (UMA/parentDFT) and solvent
(ALPB/SMD). No pKa correction or fitted offset enters their raw predictions.
The frozen results are unfavorable: UMA/ALPB MAE is 27.46 kJ/mol on gas geometries
and 31.21 on water geometries; the respective UMA/SMD values are 39.02 and 36.65.
The source compilation and related constants were previously consulted during
development, and underlying solvation-parameter exposure is not fully indexed.
This is a chemical component control, not independent metabolic validation.
Scoring the same equations as transformed metabolite pools would make them
identities; quantity_scope_audit.json proves that distinction algebraically.

Independent QR/native replay reconstructs 3876 retained solvent executions.
Three additional charged CPU fixtures reproduce the GPU calculation within
0.00030 kJ/mol. This rules out that tested implementation discrepancy as an
explanation of the large errors, without validating the physical models.
component_diagnosis.json records a post-outcome base-cancelling contrast; it
does not add another independent observation or change the two-reference MAEs.

Replay from the MetaG repository, using NumPy and ASE without CUDA, UMA assets,
PySCF or xTB:

```bash
python benchmark/reference_transfer.py
```

For this folder alone use `python reference_transfer.py .`. The archive retains
all initial scientific inputs, frozen code/analysis, optimization traces,
native solvent records, Hessians, SCF outputs, CPU fixtures and ownership records.
Original absolute producer paths are provenance, not portable launch commands.
The independent record-only replay does not authenticate unavailable model assets
or recompute quantum calculations. No production reference policy is adopted,
and no primary independent biochemical reference or ModelSEED prediction is added.
