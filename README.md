# LADDReduce

LADDReduce is a MATLAB toolkit for reducing laser-ablation double-dating
(LADD) measurements. It combines helium measurements, U-Th-(Sm)
measurements, and nested-pit volumes to calculate sub-grain (U-Th[-Sm])/He
dates and 1-sigma uncertainties for apatite or zircon.

The public workflow and numerical core have been validated with synthetic
tests and archived apatite and zircon production datasets.

## Input assumption

The public workflow assumes that the mass-spectrometer export already
contains the ordinary blank-corrected helium signal produced by the
acquisition software.

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox for parent-isotope reduction

The package has been checked with MATLAB R2025b.

See `USER_GUIDE.md` for input definitions, calibration modes, output columns,
review flags, exclusions, troubleshooting context, and a reporting checklist.

## Installation

Download or clone the package, then add `src` to the MATLAB path:

```matlab
addpath('/path/to/LADDReduce/src')
```

## Workflow

1. Reduce the ordinary blank-corrected helium signal:

```matlab
he = ladd_reduce_helium( ...
    'helium_input.csv', ...
    'he_pit_volumes.csv', ...
    'air_calibration.csv', ...
    'helium_metadata.csv', ...
    'Mineral', 'apatite', ...
    'OutputFile', 'helium_reduced.csv');
```

2. Reduce parent isotopes from the time-series files and nested-pit volumes.
   The standard lookup table supplies the user's reference-material values;
   no mineral reference material or concentration is built into the code.
   The lookup may also supply absolute 1-sigma uncertainties for those
   concentrations, which are propagated into the final age uncertainty:

```matlab
parents = ladd_reduce_apatite( ...
    'parent_timeseries', ...
    'parent_metadata.csv', ...
    'he_pit_volumes.csv', ...
    'uth_pit_volumes.csv', ...
    'BridgeStandardName', 'ReferenceMaterial', ...
    'ReferenceLookupFile', 'reference_material_lookup.csv', ...
    'SmReferenceBasis', 'total', ...
    'OutputFile', 'parents_reduced.csv');
```

For zircon, call `ladd_reduce_zircon`. For either mineral, the fourth input can
be a per-analysis U-Th pit-volume CSV or a positive, explicitly declared
session-average volume; the average also requires `UthAverage1SD`.

3. Match the two reductions by `GrainID` and calculate dates:

```matlab
ages = ladd_calculate_ages( ...
    'helium_reduced.csv', ...
    'parents_reduced.csv', ...
    'Mineral', 'apatite', ...
    'OutputFile', 'ladd_ages.csv');
```

Input column templates are in `examples/templates`. A complete generated
apatite-and-zircon walkthrough is in `examples/synthetic`:

```matlab
addpath('/path/to/LADDReduce/examples/synthetic')
example = run_synthetic_examples;
```

## Parent-reduction choices

The primary measured-volume workflow subtracts each grain's inner He-pit
volume from its measured outer U-Th-pit volume. A session-average outer-pit
volume remains available as a documented alternative when individual outer
measurements are unavailable. Every output records which volume mode and
source were used.

Bridge calibration is available for both minerals. `BridgeStandardName` is
the exact user-defined `stdname` in the metadata and lookup table; no mineral
reference-material identity is hardcoded. NIST612 is optional for bridge
calibration: when present it supplies an independent comparison, and when
absent the user explicitly selects `AllowBridgeOnly`, using at least three
replicates of the named bridge material. Apatite's direct-NIST sensitivity
modes require actual NIST612 analyses.

Apatite reduction measures U, Th, and Sm and requires the declared Sm
reference basis. Zircon reduction measures U and Th without an Sm term. The
age-calculation step then joins the He and parent outputs by `GrainID` and uses
their atoms-per-gram values and 1-sigma uncertainties.

Reference-material concentration uncertainties are supplied through the optional
`known_u_1sd_ppm`, `known_th_1sd_ppm`, and `known_sm_1sd_ppm` columns. Historical
inputs without these columns remain valid, but the output records
`NOT_SUPPLIED_ASSUMED_ZERO` because that run does not contain the complete
reference-composition uncertainty contribution described by the LADD method.

## Review and exclusion policy

Review flags report calibration, matching, or analytical conditions without
automatically removing an analysis. A statistical flag alone is not treated
as evidence that a grain is invalid. Explicit exclusions are supported only
through the corresponding named options and should be tied to an independently
documented problem.

Parent time-series integration is always reported as median signal-plateau CPS
minus median pre-ablation background CPS. Historical row-sum integration is
not available in the public package.

## Package layout

- `src/`: four user-facing entry points and their calculation engines
- `examples/templates/`: generic CSV templates
- `examples/synthetic/`: runnable end-to-end apatite and zircon examples
- `tests/`: self-contained synthetic regression and safeguard tests

## Run the tests

```matlab
results = runtests('/path/to/LADDReduce/tests');
assert(all([results.Passed]));
```

The tests verify ordinary blank-corrected helium reduction, guard the CPS
integration unit, exercise the documented calibration and pit-volume options,
and run both mineral workflows through final age calculation.
