# Apatite U-Pb runners

## Detrital samples: primary workflow

For detrital apatite, do not fit every grain in a sample to one common-Pb
isochron as the primary result. The grains need not share one age or one
initial/common-Pb composition. Use the grain-by-grain detrital launcher:

```matlab
result = ladd_apatite_upb_detrital;
```

Choose the folder containing the iolite CSV files. The launcher refreshes the
normalized input, applies IsoplotR's age-dependent Stacey-Kramers common-Pb
correction independently to each grain, and tests one through four
error-aware Gaussian components for each sample. It writes:

- `apatite_upb_detrital_distribution_summary.csv`: corrected unweighted mean,
  median, spread, correction shift, and high-level distribution class.
- `apatite_upb_detrital_population_modes.csv`: candidate mode ages,
  approximate uncertainties, proportions, intrinsic dispersion, and role.
- `apatite_upb_detrital_corrected_grains.csv`: corrected grain ages,
  analytical uncertainty, common-Pb estimate, and mixture membership.
- `apatite_upb_detrital_model_selection.csv`: BIC results for every candidate
  component count.
- one before/after and fitted-distribution PNG for each sample.

The mean and median describe the corrected detrital distribution; they are not
weighted-mean crystallisation ages. Mixture components are exploratory and are
labelled separately as candidate tight modes, broad components, or minor tails.
No grain is discarded merely for having an inconvenient age.

Set `RefreshNormalizedInput` to `false` only when intentionally reusing the
existing normalized input:

```matlab
result = ladd_apatite_upb_detrital("/path/to/ApPb", ...
    "RefreshNormalizedInput", false, ...
    "MaxComponents", 4);
```

## Cogenetic samples: common-Pb isochron workflow

This runner batches cogenetic apatite analyses by sample and uses the same
IsoplotR/Ludwig semitotal-Pb/U maximum-likelihood regression as the IsoplotR
site. It estimates the sample-specific common `207Pb/206Pb` composition and
the isochron age together. It also writes diagnostic common-Pb-corrected
ratios and per-analysis dates.

## Quick start

For a folder containing one CSV per sample, use the batch launcher:

```matlab
result = ladd_apatite_upb_batch;
```

Choose the folder when prompted. The batch runner leaves the raw CSVs
unchanged, treats `2SE(prop)` as absolute propagated 2SE by default, writes a
single sample summary, and creates one Tera-Wasserburg plot per sample. If the
output folder contains `apatite_upb_column_overrides.csv`, those documented
column corrections are reused automatically.

For one file at a time:

In MATLAB:

```matlab
result = ladd_apatite_upb_isochron;
```

Choose the iolite CSV or Excel export when prompted. The runner attempts to
identify the sample, analysis, `207Pb/235U`, `206Pb/238U`, uncertainty, and
error-correlation columns. Review `apatite_upb_column_mapping.csv` in the
output folder before interpreting the ages.

If the automatic mapping is ambiguous, name the columns explicitly:

```matlab
result = ladd_apatite_upb_isochron("my_iolite_export.csv", ...
    "SampleColumn", "Sample", ...
    "AnalysisColumn", "Analysis", ...
    "Pb207U235Column", "Final Pb207/U235 mean", ...
    "Pb207U235ErrorColumn", "Final Pb207/U235 2SE", ...
    "Pb206U238Column", "Final Pb206/U238 mean", ...
    "Pb206U238ErrorColumn", "Final Pb206/U238 2SE", ...
    "RhoColumn", "Error correlation", ...
    "InputErrors", "2se_abs");
```

The accepted uncertainty declarations are `1se_abs`, `2se_abs`, `1se_pct`,
and `2se_pct`. The normalized input passed to IsoplotR is always 1SE absolute.

## Outputs

- `apatite_upb_isochron_summary.csv`: sample ages, internal and external
  uncertainties, common-Pb composition, MSWD, p-value, and fit status.
- `apatite_upb_commonPb_corrected_analyses.csv`: corrected ratios and
  per-analysis diagnostic dates.
- one numbered Tera-Wasserburg isochron PNG per sample.
- `apatite_upb_column_mapping.csv`: the exact source-to-analysis mapping.
- `apatite_upb_excluded_or_invalid.csv`: rows not used in a fit.
- `apatite_upb_methods.txt`: calculation and uncertainty provenance.

No analysis is automatically rejected because of a high MSWD or an old age.
If the model-1 p-value is below 0.05, the summary reports the IsoplotR-style
MSWD-expanded uncertainty and a model-3 overdispersion sensitivity result.

## Disequilibrium

The default is no U-series disequilibrium correction:

```matlab
result = ladd_apatite_upb_isochron("my_iolite_export.csv", ...
    "DisequilibriumMode", "none");
```

A fixed initial `230Th/238U` activity ratio can be tested explicitly:

```matlab
result = ladd_apatite_upb_isochron("my_iolite_export.csv", ...
    "DisequilibriumMode", "fixed_th230_u238", ...
    "InitialTh230U238", 0.5);
```

Do not describe `0.5` as a magma `Th/U` correction when the input only contains
the two U-Pb ratios. A defensible magma-Th/U correction requires the additional
Th-bearing measurements used by IsoplotR formats 7 or 8. For most non-Quaternary
apatite datasets, run the common-Pb model first and treat disequilibrium as a
documented sensitivity test unless the required Th data and petrologic basis
are available.

## Interpretation checks

- Fit only analyses expected to be cogenetic and to share one common-Pb
  composition.
- Preserve the ratio error correlation (`rho`) from iolite when available.
  The runner can use zero when it is absent, but flags this assumption.
- A useful isochron needs real spread in U/Pb or common-Pb fraction. A tight
  cluster can return a numerical fit with a poorly constrained intercept.
- Treat `EXCESS_SCATTER`, implausible common-Pb intercepts, or unstable model-3
  results as evidence to inspect zoning, inclusions, mixed age domains, Pb
  loss, and sample grouping—not as permission to delete points automatically.
- Per-analysis corrected uncertainties are analytical diagnostics. The sample
  isochron age and its shared uncertainty are the primary result.

## R dependency

The MATLAB runner uses `Rscript` and a project-local IsoplotR installation in
`.r-lib`. If the project is moved to another computer, install it from this
folder with:

```sh
mkdir -p .r-lib
Rscript -e "install.packages('IsoplotR',repos='https://cloud.r-project.org',lib='.r-lib')"
```
