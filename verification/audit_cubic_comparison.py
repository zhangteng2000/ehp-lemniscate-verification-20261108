"""Exact scalar and Laurent-polynomial checks for the cubic comparison.

The manuscript proves the existence, ray coverage and integral inequalities
analytically. These checks are supplementary, not degree or angle sampling.
Dahlke's cubic reduction: https://doi.org/10.5281/zenodo.20318611 .
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

if not __debug__:
    raise RuntimeError("Run without -O or -OO.")


def add(*polys):
    h = {}
    for p in polys:
        for j, c in p.items():
            h[j] = h.get(j, 0) + c
    return {j: c for j, c in h.items() if c}


def scale(k, p):
    return {j: k*c for j, c in p.items() if k*c}


def mul(p, q):
    h = {}
    for i, a in p.items():
        for j, b in q.items():
            h[i+j] = h.get(i+j, 0) + a*b
    return add(h)


def deriv(p):
    return {j-1: j*c for j, c in p.items() if j*c}


def verify():
    one = {0: 1}
    u = {1: 1}
    plus, minus = add(one, u), add(one, scale(-1, u))
    # Laurent identity for the Joukowski branches and their ray parameter.
    J = {1: Q(1, 2), -1: Q(1, 2)}
    T = add(scale(4, mul(J, mul(J, J))), scale(-3, J))
    assert T == {3: Q(1, 2), -3: Q(1, 2)}
    num, den = {3: 1, 0: -1}, {3: 1, 0: 1}
    assert mul(add(T, scale(-1, one)), mul(den, den)) == mul(add(T, one), mul(num, num))
    # Bound the radial envelope, with u=rho^2. Both sides before squaring
    # are nonnegative on u>=1.
    major = {1: 1, 0: Q(-1, 2), -1: Q(1, 2)}
    assert add(mul(major, major), {2: -1, 1: 1, 0: -1}) == {0: Q(1, 4), -1: Q(-1, 2), -2: Q(1, 4)}
    # f'/f=N/D for f^3=(1-u)(1+u)^2; differentiate a second time.
    G = mul(minus, mul(plus, plus))
    N, D = {0: 1, 1: -3}, {0: 3, 2: -3}
    assert mul(D, deriv(G)) == scale(3, mul(N, G))
    assert add(mul(deriv(N), D), scale(-1, mul(N, deriv(D))), mul(N, N)) == {0: -8}
    # Integrating 3f'' <= -(8/3)(1+u/3) twice gives the cubic numerator.
    major = {0: 3, 1: 1, 2: -Q(4, 3), 3: -Q(4, 27)}
    assert deriv(deriv(major)) == {0: -Q(8, 3), 1: -Q(8, 9)}
    assert add(major, minus) == add(scale(Q(4, 3), mul(plus, minus)),
                                     {0: Q(8, 3), 3: -Q(4, 27)})
    # Integral of t^4(1-t)^4/(1+t^2) is 22/7 - pi.
    quotient = {6: 1, 5: -4, 4: 5, 2: -4, 0: 4}
    assert add(mul({0: 1, 2: 1}, quotient), {0: -4}) == {
        4+k: (-1)**k*comb(4, k) for k in range(5)}
    pi_upper = sum(Q(v, j+1) for j, v in quotient.items())
    assert pi_upper == Q(22, 7)
    sqrt2_upper = Q(99, 70)
    assert sqrt2_upper**2 > 2
    upper = pi_upper*(2+4*sqrt2_upper)/3 - 4/(81*sqrt2_upper)
    assert 8-upper == Q(25832, 1964655) > 0
    # Remaining bracket after arctan(lambda)<lambda; v=lambda^2.
    bracket = add(scale(-9, plus), scale(Q(6, 7), mul(plus, u)), plus, {0: 8})
    assert bracket == {1: -Q(50, 7), 2: Q(6, 7)}
    assert bracket[1]+bracket[2] < 0
    receipt = {"uniform_small_ray_upper": str(upper), "margin_to_8": str(8-upper),
               "final_bracket_upper_coefficients": ["0", "-50/7", "6/7"],
               "parameter_samples": 0, "scope": "Exact identities and scalar bounds in the analytic cubic proof."}
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(receipt, indent=2)+"\n", encoding="utf-8", newline="\n")
    print("Cubic: Laurent identities and both uniform inequalities passed.")
    return receipt


if __name__ == "__main__":
    verify()
