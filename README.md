# Polynomial lemniscate verification

This private repository contains the computation code and certificate data
accompanying a manuscript on the maximal arclength of polynomial lemniscates.
The analytic proof is supplied separately. There are no TeX or PDF files here,
and no document archives or precompiled executables are inputs to the replay.

The files were extracted from the sealed verification supplement dated
8 November 2026. The low-degree runner was adapted to read the extracted
`prior/quartic/` and `prior/quintic/` directories directly.

## Reproduce the checks

Use a 64-bit Python 3.12 interpreter, with assertions enabled. The quartic
and the second quintic reconstruction require:

```text
numpy==2.3.5
sympy==1.14.0
python-flint==0.9.0
```

Install these dependencies, enter the repository, and run:

```text
python -B verification/replay_all_degrees.py --workers 2
```

Windows needs the 64-bit .NET Framework C# compiler. On Unix the runners use
`mcs` and `mono`. Set `CSC` to a compiler path when necessary. A separate
python-flint installation can be selected with `EHP_PYTHON_FLINT`.
Once dependencies are installed, verification is offline. The programs
rebuild their counters and tables; allow several GiB of free memory and
sufficient runtime for the complete parameter covers.

The `--part earlier`, `--part high`, and `--part extension` options replay
the quartic/quintic, degree-at-least-eight, and sextic/septic parts separately.
The extension also retains an optional independent cubic computation. The
revised manuscript instead cites Dahlke's direct cubic proof.

## What each directory verifies

* `prior/quartic/`: the complete compact parameter domain for quartics,
  including both collisions. The two integer root-counting implementations
  share analytic parameter enclosures but reconstruct source controls by
  different methods. The certificate has 20,809 positive leaves and 646
  excluded leaves.
* `prior/quintic/`: both complete cubic-discriminant phase domains, the
  critical-point configuration cover, and the diameter comparisons.
* `verification/`: analytic endpoint inequalities for the unbounded degree
  range, finite interpolation and coefficient bounds, real-root strip
  estimates, and transfer from lines to circles.
* `verification/cubic_terminal_6/` and `cubic_terminal_7/`: rational goals,
  complete subdivision trees, exact integer counters, and a separate
  factored evaluator for the sextic and septic comparisons.

All acceptance tests concern entire closed parameter boxes or explicit
scalar inequalities. Floating output summaries are not acceptance tests.
The exact counters are compiled with checked integer arithmetic. Sine
enclosures are reconstructed using rational bounds and integer square
roots. The alternate floating discriminant evaluator uses a separately
checked uniform error allowance.

Some legacy C# progress lines print a `proved` count for an older constant.
The terminal-tree Python verifier independently evaluates the actual
rational goal for every recorded leaf. The Python acceptance assertions,
and the final empty unresolved set, determine success.

## Scope

Executing the code does not verify the analytic lemmas connecting these
finite inequalities with lemniscate length. Those lemmas must be checked
in the accompanying manuscript. The source-provenance audit and final
replay summary will be added when the current review is complete.

No repository license or authorship information is inferred from the
manuscript's mathematical claims. This repository remains private.
