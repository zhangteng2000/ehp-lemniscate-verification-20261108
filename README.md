# Polynomial lemniscate verification

This private repository contains the computation code and certificate data
accompanying a manuscript on the maximal arclength of polynomial lemniscates.
The analytic proof is supplied separately. There are no TeX or PDF files here,
and no document archives or precompiled executables are inputs to the replay.

The original code was extracted from the sealed supplement dated 8 November
2026. The subsequent proof simplification adds a one-parameter quartic line
check, signed short-arc estimates, and shorter finite circle comparisons.
The archived low-degree
sources remain in `prior/` for comparison.

## Reproduce the checks

Use a 64-bit Python 3.12 interpreter, with assertions enabled. The complete
current replay uses NumPy for unbounded-integer array operations in one
quintic reconstruction and SymPy for the symbolic cross-checks:

```text
numpy==2.3.5
sympy==1.14.0
```

Install these dependencies, enter the repository, and run:

```text
python -B verification/replay_all_degrees.py --workers 2
```

Windows needs the 64-bit .NET Framework C# compiler. On Unix the runners use
`mcs` and `mono`. Set `CSC` to a compiler path when necessary; the new quartic line runner
also accepts `EHP_CSC`. A separate
python-flint installation can be selected with `EHP_PYTHON_FLINT` only for
the optional archived quartic replay. That optional replay also needs
`python-flint==0.9.0` from `requirements-archival.txt`; the current proof
does not use Arb or python-flint.
The complete replay has been run on Windows; the Unix compiler paths are
provided but have not been exercised in this audit. The current replay entry points
reject `python -O` and `python -OO` explicitly, because disabling assertions
would disable their acceptance tests.
Once dependencies are installed, verification is offline. The programs
rebuild their counters and tables; allow several GiB of free memory and
sufficient runtime for the complete parameter covers.
The `.gitattributes` file preserves all supplement files byte for byte,
so cloning with automatic line-ending conversion does not alter the
certificate inputs or their recorded hashes.

The `--part earlier`, `--part high`, and `--part extension` options replay
the current quartic/quintic, degree-at-least-eight, and degree-six/seven
parts separately. The `--part archival-quartic` option replays the superseded
two-parameter quartic proof. The revised manuscript cites Dahlke's direct
cubic proof, so the older cubic computation is not part of the main replay.

## What each directory verifies

* `verification/audit_cubic_literature.py`: exact rational audit of
  the elementary constants in the imported direct cubic proof, together
  with its Laurent-polynomial identities and final strict sign. This
  supplements primary-text review of the branch coverage and integrals.
* `verification/inverse_beta_moment_audit.py`: the rational endpoints
  in the inverse-beta moment lemma, which replaces the three-event
  probability estimate for recovery of high coefficients. The centered
  moment recurrence is proved analytically; the extra raw-moment identity
  tests are not used as a sampling argument for unbounded degrees.
* `verification/uniform_circle_transfer.py`: the rational arithmetic in
  the uniform analytic circle implication for every degree at least 28.
  This check has no degree loop and no diameter subdivision.
* `verification/quartic_line/`: the current quartic proof. Both generators
  rebuild three closed radical charts, 128 intervals each; complex
  conjugation covers the fourth branch. The two
  checked integer circuits count positive discriminant rectangles separately.
  Independent rational calculations verify the line-to-circle constants
  and the quantitative gap coefficient 1/32. No NumPy, SymPy, or flint is
  used in this quartic check.
* `verification/stationary_circle_transfer.py` and
  `independent_stationary_circle_transfer.py`: the current common circle
  comparison in degrees five through eight. Stationarity gives an extra
  orthogonal Laurent norm, and a direct affine-flow calculation gives a
  second constant coefficient. The primary bisection covers have
  4, 10, 10, and 6 closed intervals; the independent trisection covers
  have 4, 12, 12, and 6. The programs reconstruct all interval endpoints,
  scalar constants, and the short-arc curvature bounds independently.
* `verification/audit_stationary_identity.py` and
  `audit_stationary_direct.py`: supplemental algebra checks. The second
  reconstructs the polynomial division flow and inverse germs directly
  over exact dual rationals. These finite identity checks support the
  general analytic derivation; they do not replace it.
* `verification/audit_second_coefficient_extraction.py`: an optional
  supplemental check for the shorter second-coefficient proof. It reconstructs
  the polynomial flow directly over dual rationals and checks 36 exact
  identities in twelve cases. Run it with
  `python -B verification/audit_second_coefficient_extraction.py`.
  The universal identity is proved analytically in the manuscript; these
  examples are not a proof for arbitrary degree and are not inputs to the
  scalar certificate acceptance tests.
* `verification/audit_capacity_majorants_exact.py`: the explicit rational
  algebra replacing the two adaptive unit-ball majorant covers.
* `verification/single_coefficient_transfer.py` and
  `independent_single_coefficient_transfer.py`: one complete exterior
  coefficient sum, rather than a finite row truncation, proves the circle
  comparisons in every degree from nine through 27. The two programs use
  different constant enclosures and different interval partitions.
* `verification/audit_octic_majorant_exact.py` and
  `audit_nonic_majorant_exact.py`: elementary rational inequalities for
  the nonlinear inverse coefficients in degrees eight and nine. The
  degree-nine bound is 15/16; its proof shares the quadratic area budget
  with the cubic and quartic terms.
* `verification/audit_inverse_coefficients.py` and
  `audit_common_identities.py`: symbolic checks of the inverse logarithm,
  resultant signs, nonlinear coefficients, cubic discriminant, and row
  identities. They supplement the analytic derivations rather than
  replace them.
* `prior/quartic/`: the archived complete compact parameter domain for quartics,
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

Several older drivers and modules retain superseded circle comparisons
as historical checks, including `quintic_circle_simplified.py`,
`octic_circle_simplified.py`, `verify_arc_capacity_transfer.py`, and their
companions. The current circle arguments use the stationary two-coefficient
comparison in degrees 5--8, the complete one-coefficient sum in degrees
9--27, and the uniform analytic implication from degree 28 onward. Their
validity does not depend on the older diameter partitions.

Some legacy C# progress lines print a `proved` count for an older constant.
The terminal-tree Python verifier independently evaluates the actual
rational goal for every recorded leaf. The Python acceptance assertions,
and the final empty unresolved set, determine success.

## Scope

Executing the code does not verify the analytic lemmas connecting these
finite inequalities with lemniscate length. Those lemmas must be checked
in the accompanying manuscript. The audit directory records the source-provenance review and the final
replay summary. Historical alternative checks are identified above.

No repository license or authorship information is inferred from the
manuscript's mathematical claims. This repository remains private.
