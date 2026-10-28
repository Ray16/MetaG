# Benchmark: Mg-bounded openTECR reference

`reference_mg_bounded.json` holds the experimental ΔrG'° labels (298.15 K reference conditions, pH 7,
I = 0, no Mg²⁺ binding) for the reactions of the packaged panel
(`src/metag/data/reactions_opentecr_std.json`). They were derived from the openTECR re-curation of TECRDB
by a Legendre transform of every measurement to the target conditions with the eQuilibrator dissociation
model. A measurement with unknown free Mg²⁺ is admitted only when its label changes by at most 2 kJ/mol
over a physically bounded pMg interval. The rule was fixed before any MetaG residual was computed
([PREREGISTRATION.md](PREREGISTRATION.md)). Thresholds 1 and 4 kJ/mol are reported as sensitivities.

```bash
python benchmark/score.py                              # MetaG 0.3.0 panel predictions
python benchmark/score.py metag_out/metag_all.json     # a fresh `metag score` run
```

`metag_0.3.0_predictions.json` holds the 0.3.0 panel predictions (320 estimates for 364 inputs). With the
primary policy they give MAE 8.42 kJ/mol on 239 of 274 labels; the zero predictor gives 11.68 on the
same labels. They combine the 0.2.0.dev24 outcomes for reactions the 0.3.0 corrections do not touch
with fresh 0.3.0 calculations for the 24 that they change; see [docs/validation.md](../docs/validation.md).
`metag_dev24_predictions.json` keeps the dev24 outcomes (MAE 8.44) for reference.

Caveats, as reported in the manuscript:

- Discrete method choices were made while looking at TECRDB residuals, so these scores are
  development-exposed rather than independent validation. No parameter fitted to reaction free energies
  enters an estimate.
- Labels are medians over measurements at their native temperatures (273–348 K). Many reactions have no
  measurement at 298.15 K.
- Ionic strength is imputed as 0.25 M where it was not reported.
- Fitted comparators (component contribution, dGbyG, dGPredictor) have to be compared by regime: in-sample,
  out-of-fold, or compound held out. Their predictions are not redistributed here.

Data provenance: openTECR (Wagner et al.) and the original TECRDB measurements. Cite them when you use
these labels.
