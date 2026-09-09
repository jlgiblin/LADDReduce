# Synthetic end-to-end examples

`run_synthetic_examples.m` creates complete, generic apatite and zircon input
sets and runs them through helium reduction, parent-isotope reduction, grain
matching, and age calculation.

The generated data are mathematical fixtures, not measurements and not a
recommended acquisition design. Their purpose is to show file relationships,
exercise the public interface, and provide known target ages for verification.

From MATLAB:

```matlab
addpath('/path/to/LADDReduce/examples/synthetic')
result = run_synthetic_examples;
```

By default, the inputs and outputs are written to a new timestamped folder in
MATLAB's temporary directory. To retain them somewhere specific, supply a new
folder path that does not already exist:

```matlab
result = run_synthetic_examples('/path/to/new/example_output');
```

Each mineral folder contains:

- ordinary blank-corrected helium input;
- air calibration and helium metadata tables;
- He and U-Th pit-volume tables;
- raw-style parent-isotope time-series CSVs;
- parent metadata and a reference-material lookup table;
- reduced helium and parent tables; and
- the final LADD age table.

The apatite example targets 50 and 65 Ma. The zircon example targets 80 and
100 Ma. The script stops with an error if either workflow fails to recover its
targets within 0.02 Ma.
