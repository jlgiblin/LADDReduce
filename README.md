# LADDReduce

LADDReduce is a MATLAB toolkit for reducing laser-ablation double-dating
(LADD) measurements. It combines helium measurements, U-Th-(Sm)
measurements, and nested-pit volumes to calculate sub-grain (U-Th[-Sm])/He
dates and 1-sigma uncertainties for apatite or zircon.

This is a release-candidate package. Its numerical core and safeguards are
being tested before the first public release.

## Input assumption

The public workflow assumes that the mass-spectrometer export already
contains the ordinary blank-corrected helium signal produced by the
acquisition software.

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox for parent-isotope reduction
- R and IsoplotR only if the optional apatite U-Pb module is used

The release candidate has been checked with MATLAB R2025b.

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
    'sample_types.csv', ...
    'Mineral', 'apatite', ...
    'OutputFile', 'helium_reduced.csv');
```

2. Reduce parent isotopes from the time-series files and measured nested-pit
   volumes. Reference-material identity and the Sm reference basis must be
   declared explicitly:

```matlab
parents = ladd_reduce_apatite( ...
    'parent_timeseries', ...
    'parent_metadata.csv', ...
    'he_pit_volumes.csv', ...
    'uth_pit_volumes.csv', ...
    'BridgeStandardName', 'ReferenceMaterial', ...
    'SmReferenceBasis', 'total', ...
    'OutputFile', 'parents_reduced.csv');
```

For zircon, call `ladd_reduce_zircon`. Its fourth input can be either a
per-analysis U-Th pit-volume CSV or a positive, explicitly declared session
average in cubic micrometres.

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
- `optional_apatite_upb/`: optional apatite U-Pb/common-Pb tools

## Run the tests

```matlab
results = runtests('/path/to/LADDReduce/tests');
assert(all([results.Passed]));
```

The current tests verify ordinary blank-corrected helium reduction, guard the
CPS integration unit, and run both mineral workflows through final age
calculation. The cleaned reducer also reproduces three preserved accepted
fixture calculations spanning apatite and two zircon instrument exports;
those private fixtures are not included in the public package.
