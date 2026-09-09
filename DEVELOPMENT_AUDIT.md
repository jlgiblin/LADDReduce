# LADDReduce public-package audit

## Package decision

The original `LADDReduce` working directory is a mixed research workspace, not
a release-ready repository. It contains reusable reducers, session-specific
runners, private sample identifiers, raw and reduced data, working documents,
diagnostics, historical sensitivity code, and unrelated project outputs.
Publishing that directory directly would expose unrelated material and make the
supported workflow difficult to identify.

The release is therefore being assembled in a separate folder. No original
code or result has been removed or overwritten.

## Public core selected for cleanup

- Ordinary blank-corrected helium reduction
- Apatite U-Th-Sm reduction from measured nested-pit volumes
- Zircon U-Th reduction from measured or explicitly declared nested-pit
  volumes
- Grain matching and (U-Th[-Sm])/He date calculation
- Metadata enrichment using user-supplied reference-material values
- Optional apatite U-Pb/common-Pb tools, kept in a separate module

Four short public entry points now wrap the more detailed reduction engines:

- `ladd_reduce_helium`
- `ladd_reduce_apatite`
- `ladd_reduce_zircon`
- `ladd_calculate_ages`

These wrappers require users to declare mineral type, bridge-standard identity,
and Sm reference basis where applicable.

## Material intentionally excluded

- Lab session runners and session-specific metadata repairs
- Research sample identifiers, private paths, raw data, and working documents
- Instrument-specific empirical standard values presented as universal
  defaults
- Pecube, DetritalChronFilter, and MultichronFitTSF working outputs

## Validation baseline

- A self-contained synthetic test guards median signal-minus-background CPS.
- The cleaned reducer reproduced the accepted CPS values for apatite and for
  both instrument generations represented by three private zircon/apatite
  fixtures.
- A self-contained helium test verifies ordinary blank-corrected reduction and
  validates the public output schema.
- A self-contained end-to-end generator runs generic apatite and zircon
  datasets through all four public entry points and recovers their declared
  target ages.
- The existing zircon policy test suite is stale: two tests point to metadata
  files that have moved, and one expected rejection no longer occurs under the
  current implementation. This test suite cannot be presented publicly as
  passing until it is repaired or replaced with isolated synthetic fixtures.
- MATLAB's code analyzer found no syntax errors in the public entry points or
  engines. Remaining messages in the public engines are advisory style and
  performance notes, not calculation failures.

## Remaining release blockers

1. Repair or replace the stale zircon metadata-policy tests.
2. Complete the detailed user manual.
3. Perform a clean-folder installation test.
4. Initialize and review the public Git repository.
