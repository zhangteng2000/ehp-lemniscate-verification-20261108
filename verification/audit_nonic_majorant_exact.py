"""Exact proof constants for |R_9| < 9/10 on the area-theorem ball."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json

if not __debug__:
    raise RuntimeError("Run without -O or -OO.")


def verify():
    # Expand sqrt(1-Z)-1-log(1-Z)/2 where Z=sum j*a_j*t^(j+1).
    # Derive each coefficient using the binomial recurrence.
    binom = Q(1)
    terms = {}
    for k in range(1, 5):
        binom *= (Q(1, 2)-k+1)/k
        weight = (-1)**k*binom + Q(1, 2*k)
        for js in product(range(1, 7), repeat=k):
            if sum(js)+k != 9 or not weight:
                continue
            coef = weight
            for j in js:
                coef *= j
            key = tuple(sorted(js))
            terms[key] = terms.get(key, 0)+coef
    assert terms == {(1,6):Q(3,2), (2,5):Q(5,2), (3,4):Q(3),
                     (1,1,4):Q(5,4), (1,2,3):Q(15,4), (2,2,2):Q(5,6),
                     (1,1,1,2):Q(11,16)}
    # y_j = sqrt(j)|a_j|. Pair high and low indices by Young's inequality.
    a, b = Q(11,40), Q(23,40)
    assert 4*a*b > Q(10,16) >= Q(6,16)
    assert b*b > Q(12,64)
    assert Q(5,24)**2*2 < Q(3,10)**2
    assert Q(11,64)**2*2 < Q(3,10)**2
    assert a+Q(3,10) == b
    # u(6s-5u) = 9s^2/5 - 5(u-3s/5)^2,
    # and 4/27-s^2(1-s) = (s-2/3)^2(s+1/3).
    assert (5*Q(3,5)**2, 10*Q(3,5)) == (Q(9,5), 6)
    y = Q(2,3)
    assert (y*y/3, y*y-2*y/3, Q(1,3)-2*y) == (Q(4,27), 0, -1)
    square = Q(5,8)**2*Q(9,5)*Q(4,27)
    assert square == Q(15,144) < Q(13,40)**2
    upper = b+Q(13,40)
    assert upper == Q(9,10)
    out = {"R9_upper":str(upper), "signed_arc_upper_constant":str(upper/2),
           "quadratic_and_diagonal_upper":str(b), "remaining_square_upper":str(square),
           "parameter_subdivisions":0}
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(out, indent=2)+"\n", encoding="utf-8", newline="\n")
    print("Nonic: |R9| < 9/10; signed short-arc constant 9/20.")
    return out


if __name__ == "__main__":
    verify()
