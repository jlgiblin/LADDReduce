# LADDReduce user guide

## 1. Purpose and scope

LADDReduce reduces laser-ablation double-dating measurements for apatite and
zircon. The supported workflow combines:

1. an ordinary blank-corrected helium measurement from the He pit;
2. an outer U-Th-(Sm) measurement;
3. measured He and U-Th pit volumes; and
4. user-declared calibration and reference-material information.

The final calculation reports sub-grain (U-Th-Sm)/He dates for apatite or
(U-Th)/He dates for zircon. No alpha-ejection correction is applied because
the calculation represents nested ablation volumes rather than a conventional
whole grain.

## 2. Requirements and setup

The package requires MATLAB. Parent-isotope reduction uses Statistics and
Machine Learning Toolbox functions.

Add the source folder at the beginning of a MATLAB session:

~~~matlab
addpath('/path/to/LADDReduce/src')
~~~

All input uncertainties described below are absolute 1-sigma values. Convert
other uncertainty conventions before running LADDReduce and retain the
conversion in the analytical record.

## 3. Grain identifiers

GrainID connects the helium, pit-volume, metadata, parent-isotope, and age
tables. Use one stable identifier for each analyzed grain in every file.

Matching is case-insensitive and ignores common punctuation differences, but
do not rely on that tolerance to repair inconsistent identifiers. Inspect the
final UThMatchCount, UThMatchDecision, and flags columns.

## 4. Helium reduction

Call ladd_reduce_helium with four input files and an explicit Mineral value.

~~~matlab
he = ladd_reduce_helium( ...
    'helium_input.csv', ...
    'he_pit_volumes.csv', ...
    'air_calibration.csv', ...
    'helium_metadata.csv', ...
    'Mineral', 'apatite', ...
    'OutputFile', 'helium_reduced.csv');
~~~

### Helium measurement table

Required columns:

| Column | Meaning |
|---|---|
| RunID | Measurement or acquisition identifier |
| SampleName | Name used to join the sample-type table |
| GrainID | Stable grain identifier |
| He4_cps | Ordinary blank-corrected helium signal in counts per second |
| He4_1SD | Absolute 1-sigma uncertainty on He4_cps |

He4_cps must already contain the ordinary correction produced by the
acquisition workflow.

### Helium metadata table

The table maps SampleName to RunScript:

| RunScript | Meaning |
|---:|---|
| 1 | Air-calibration measurement |
| 2 | Blank measurement |
| 3 | Unknown sample |
| 4 | Mineral reference material |

Rows absent from this table remain unclassified and are not silently treated
as unknowns.

### Air-calibration table

Required columns are AirID, FourHeAir, and FourHeAir1SD. FourHeAir is the
calibration factor in atoms per cps; FourHeAir1SD is its absolute 1-sigma
uncertainty. Unknowns and mineral standards use the nearest valid air row in
measurement order.

### He pit-volume table

Required columns are GrainID, PitVol_um3, and PV1SD_um3. The volume uncertainty
is absolute 1 sigma. Mineral density converts the measured volume to mass.

### Principal helium outputs

The output retains the measured signal and reports pit volume, mass, assigned
air calibration, He atoms, He atoms per gram, propagated 1-sigma uncertainty,
blank diagnostics, standard diagnostics, and flags.

## 5. Parent-isotope inputs

Both mineral workflows require:

- a folder of raw-style parent-isotope time-series CSV files;
- one metadata CSV;
- a He pit-volume CSV; and
- a U-Th pit-volume CSV or an explicitly declared session-average U-Th pit
  volume and its 1-sigma uncertainty.

The time-series files need a time channel and the declared isotope channels.
The public reducers use median signal-plateau CPS minus median pre-ablation
background CPS. Row-sum integration is not available.

### Parent metadata

Required or strongly recommended columns:

| Column | Meaning |
|---|---|
| file | Exact time-series filename |
| type | NIST612, MineralStd, or Unknown |
| stdname | Exact reference-material name; blank for unknowns |
| grainid | Stable GrainID |
| known_u_ppm | Declared U concentration for calibration rows |
| known_th_ppm | Declared Th concentration for calibration rows |
| known_sm_ppm | Declared Sm concentration for apatite calibration rows |

Rows with type NIST612 must represent NIST612 glass. The zircon reducer stops
on inconsistent type and stdname combinations instead of relabeling them.

Reference concentrations are user supplied. LADDReduce does not provide
laboratory defaults. Supply `ReferenceLookupFile` to either public parent
reducer to fill blank `known_*` metadata fields from the lookup table. Existing
nonblank values are never overwritten. Exact normalized standard names are
preferred; ambiguous partial matches stop for correction.

### Nested-pit volumes

For unknowns with both volumes, the parent-isotope volume is:

~~~text
analyzed parent volume = outer U-Th pit volume - inner He pit volume
~~~

The difference must be positive. Bridge-standard rows use their full outer
U-Th pit volume. If a session-average outer volume is used, the same declared
value and uncertainty are applied consistently and recorded in the output.
Missing, nonpositive, or mismatched volumes are retained with review flags
rather than being repaired silently.

## 6. Apatite parent reduction

~~~matlab
parents = ladd_reduce_apatite( ...
    'parent_timeseries', ...
    'parent_metadata.csv', ...
    'he_pit_volumes.csv', ...
    'uth_pit_volumes.csv', ...
    'BridgeStandardName', 'ReferenceMaterial', ...
    'ReferenceLookupFile', 'reference_material_lookup.csv', ...
    'SmReferenceBasis', 'total', ...
    'AnchorMode', 'median', ...
    'OutputFile', 'parents_reduced.csv');
~~~

BridgeStandardName must match metadata stdname. SmReferenceBasis must be total
or 147isotope and must agree with known_sm_ppm.

Apatite AnchorMode choices:

| Mode | Behavior |
|---|---|
| nearest | Uses the nearest following retained bridge replicate; the terminal block uses the nearest retained bridge |
| median | Uses one session-wide median bridge factor for each parent element |
| nist_following | Uses the next following NIST612 analysis, with terminal fallback |
| nist_interpolated | Retains the time-interpolated NIST612 calibration |

The default is nearest. Selection should reflect the reference material,
acquisition design, and intended calibration strategy.

The fourth positional input can instead be one positive numeric session-average
U-Th pit volume. In that case, `UthAverage1SD` must also be supplied. Measured
per-analysis volumes are preferable when available.

NIST612 is not mathematically required for the `nearest` or `median` bridge
calibration. If no genuine NIST612 rows are present, set `AllowBridgeOnly` to
true; at least three replicates of the named bridge material are required and
the output records that no independent NIST comparison was available. The
`nist_following` and `nist_interpolated` modes do require actual NIST612 rows.

## 7. Zircon parent reduction

~~~matlab
parents = ladd_reduce_zircon( ...
    'parent_timeseries', ...
    'parent_metadata.csv', ...
    'he_pit_volumes.csv', ...
    'uth_pit_volumes.csv', ...
    'BridgeStandardName', 'ReferenceMaterial', ...
    'ReferenceLookupFile', 'reference_material_lookup.csv', ...
    'AnchorMode', 'median', ...
    'OutputFile', 'parents_reduced.csv');
~~~

Zircon AnchorMode choices:

| Mode | Behavior |
|---|---|
| following | Uses the next following retained bridge replicate; the terminal block uses the nearest retained bridge |
| median | Uses one session-wide median bridge factor for U and Th |
| single | Uses the explicitly named SingleBridgeID |

The default is following.

The fourth positional input can instead be one positive numeric session-average
U-Th pit volume. In that case, UthAverage1SD must also be supplied explicitly.
Per-analysis measured volumes are preferable when available.

A zircon session without genuine NIST612 rows stops by default. `AllowBridgeOnly`
must be set explicitly to use at least three named bridge rows as the sole
parent calibration. That output clearly records that no independent NIST check
was available.

## 8. Age calculation

~~~matlab
ages = ladd_calculate_ages( ...
    'helium_reduced.csv', ...
    'parents_reduced.csv', ...
    'Mineral', 'apatite', ...
    'OutputFile', 'ladd_ages.csv');
~~~

The calculation matches GrainID, solves the radioactive-production equation
iteratively, and propagates helium and parent uncertainties. Important outputs
include:

| Column | Meaning |
|---|---|
| Age_Ma | Calculated age |
| Age_1SD_Ma | Absolute 1-sigma age uncertainty |
| Age_2SD_Ma | Two times the 1-sigma result |
| Age_1SDpct | Relative 1-sigma age uncertainty |
| converged | Whether the age solver met its tolerance |
| UThMatchCount | Number of parent rows matching the normalized GrainID |
| UThMatchDecision | Unique, automatically resolved duplicate, unresolved duplicate, or no match |
| ManualExclude | Explicit user exclusion; row retained with age fields blank |
| flags | Combined helium, parent, matching, and age flags |

Duplicate matches are not resolved by taking the first row. Automatic
resolution occurs only when exactly one duplicate has physically usable parent
measurements; otherwise the result requires review.

## 9. Review flags

Flags provide information; they are not automatic geological rejection rules.

| Flag | Review meaning |
|---|---|
| HIGH_BLANK | Session median blank is large relative to the grain signal |
| NO_PITVOL | Required pit volume was not matched |
| NO_AIRSTD | No usable air calibration was assigned |
| ZERO_HE | Helium signal is nonpositive |
| NOISY_HE | Helium signal uncertainty exceeds the current review criterion |
| STD_REVIEW | Mineral-standard result is far from the session median |
| STD_HAMPEL_REVIEW | Parent bridge replicate was statistically unusual but retained |
| STD_EXCLUDE or STD_CAL_EXCLUDED_MANUAL | User explicitly named a documented bad standard |
| POSSIBLE_SKIP | Grain lies near an explicitly excluded He standard in run order |
| NO_HEPV | U-Th volume exists but He volume is missing; full U-Th volume was used |
| NO_UTHPV | U-Th pit volume is missing |
| BAD_PV | Outer U-Th volume is not greater than the inner He volume |
| NO_PARENT_SIGNAL | No positive usable parent signal |
| SM_WEAK | Sm signal is below the review threshold |
| SM_UNSTABLE | Sm relative signal uncertainty is high |
| SM_OUTLIER | Sm concentration exceeds the plausibility review threshold |
| SM_NIST612_FALLBACK | Apatite Sm used the recorded NIST fallback path |
| SM_NOT_CALCULATED | No usable Sm calibration was available |
| BRIDGE_NIST_MISMATCH_REVIEW | Bridge result differs materially from the independent NIST result |
| BRIDGE_NEAREST_FALLBACK | No following bridge existed; the nearest retained bridge was used |
| NO_CONVERGE | Age solver did not meet the requested tolerance |
| NEG_AGE | Calculated age is negative |
| NO_UTH_MATCH | No unique usable parent row was matched |
| LOW_PARENT | Parent production is too low for a stable result |
| INVALID_PARENT | Parent input is nonfinite or physically invalid |
| MANUAL_EXCLUDE | User explicitly withheld the GrainID; the row remains auditable |

Review the numerical provenance columns beside each flag. A flagged row can be
valid, while an unflagged row can still be unsuitable for a particular study.

## 10. Exclusions

Statistical review flags do not automatically remove calibration rows or
grains. Explicit exclusions are available through ExcludeStandardIDs,
ExcludeNistFiles, ExcludeBridgeIDs, and ExcludeGrainIDs. Use them only for an
independently documented issue such as a misfire, misidentified material, or
confirmed sequence problem.

Source rows remain in the output with their exclusion status. Keep the
exclusion list and rationale with the run record.

## 11. Synthetic example and tests

The complete synthetic walkthrough is in examples/synthetic:

~~~matlab
addpath('/path/to/LADDReduce/examples/synthetic')
example = run_synthetic_examples;
~~~

It creates generic raw-style inputs, runs both mineral workflows, and verifies
target ages of 50 and 65 Ma for apatite and 80 and 100 Ma for zircon. The data
are mathematical fixtures, not measurements or a recommended acquisition
design.

Run the package tests with:

~~~matlab
results = runtests('/path/to/LADDReduce/tests');
assert(all([results.Passed]));
~~~

## 12. Minimum reporting checklist

For reproducible use, record:

- LADDReduce release or commit;
- MATLAB release and relevant toolbox versions;
- mineral and isotope channels;
- confirmation that He4_cps was the ordinary blank-corrected signal;
- air-calibration source and uncertainty;
- reference-material identities, concentrations, Sm basis, and citations;
- apatite or zircon anchor mode;
- measured pit-volume source and uncertainties, or the explicitly declared
  session-average volume;
- all explicit exclusions and their independent rationale; and
- the review flags retained in the reported dataset.
