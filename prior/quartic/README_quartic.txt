Verification supplement for the quartic Erdős–Herzog–Piranian lemniscate inequality
7 November 2026

The proof is the section “The quartic extremal problem” of the accompanying article.
The certificate establishes a strict inequality for every quartic whose critical
values have modulus one and include two distinct points. The proof in the article
then treats arbitrary monic quartics and all equality cases.

Reproduction
Extract the quartic directory to a writable folder. Use 64-bit Python 3.12 with
numpy==2.3.5, sympy==1.14.0 and python-flint==0.9.0. For example, run

    python -m pip install -r requirements.txt
    python verify_quartic.py --workers 3

Windows can use the 64-bit .NET Framework 4 C# compiler. On Unix the runner
calls the Mono C# compiler (mcs) and mono. Set CSC to a compiler path if needed.
For an existing python-flint installation in a separate library directory, set
EHP_PYTHON_FLINT to that directory. The proof uses Python assertions. Use the
normal interpreter with PYTHONOPTIMIZE unset. No network access is needed once
the stated packages and the compiler are installed. Allow 3 GiB of free RAM for
one verification at a time. Runtime depends on the machine. The runner generates
the source control tables, counters and reports in the extracted directory.

Certificate and checks
quartic_atlas_certificate.json stores all 224 tiles, their binary subdivisions,
and the primary integer results. Digits 0 and 1 split the first coordinate at its
midpoint. Digits 2 and 3 split the second. Every retained node branches in both
directions of its chosen coordinate. P denotes a positive lower length deficit.
O denotes a rectangle disjoint from the admissible lens. The replay follows every
edge, checks the initial tiling, and verifies every terminal assertion.
The certificate SHA-256 is
964241c49502d2d7c26dad063cca2edf679fb1d0ac64a61379fe0fc70b1059d9.

slice_data.py forms the initial Bernstein controls by rational polynomial
multiplication. independent_controls.py reconstructs them by tensor interpolation
of unfactored values on an integer grid. SliceCounter.cs uses sign-parity losses,
and CheckSlices.cs uses a two-state recurrence and a narrower symmetric control
interval. They independently bound univariate slice root counts. Both are
compiled with overflow checks. chord_refined.py and chord_geometry.py supply the
shared analytic parameter enclosures. They use Arb with directed ball bounds.
The small amount of binary floating-point in the code chooses a positive power
of two for rescaling or search priorities. The acceptance comparisons use
integer/rational inequalities or outward interval bounds.
audit_atlas.py verifies the algebraic identities, reference beta integral bounds,
discriminant-budget constants used in the sums, binary domain coverage, and all
root-count sums and comparison margins. The analytic proofs of the projection
budget, target-circle tail bound and extremizer reduction are in the article.
The two replays use the same analytic parameter formulas and Arb. They are
separate implementations of the control-table and real-root counting stages.

Expected output
The 224 tiles yield 20,809 positive rectangles and 646 excluded rectangles.
Of the positive rectangles 16,598 use the entire circle and 4,211 use a
truncation with an analytic tail bound. The common denominator of the reported
minimum margins in the table coordinate is 7,902,739,824,640. The 104-bit primary
numerator is 54,996,863. The 150-bit second numerator is 899,883,551.
The primary and secondary share the same finite partition. Both imply the
coefficient 2^(-18) in the quantitative quartic proposition.
audit_primary.json and audit_independent.json contain the original complete
replays. The numerical time fields depend on the machine.

