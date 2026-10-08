import mpmath as mp
from bolzakit import arith
from bolzakit.exact import group as G


def test_golden_geodesic_primes():
    from sympy import primerange
    hits = [p for p in primerange(3, 200) if arith.golden_geodesic_is_systole(p)]
    assert hits == [3, 7, 47]


def test_golden_rejects_nonprime():
    assert not arith.golden_geodesic_is_systole(9) and not arith.golden_geodesic_is_systole(2)


def test_golden_adversarial():
    for n in (9, 15, 49, 4, 1, 0):
        assert not arith.golden_geodesic_is_systole(n)
    assert arith.golden_geodesic_is_systole(2207)
    assert not arith.golden_geodesic_is_systole(2203)


def test_lucas():
    assert [arith.lucas(n) for n in range(8)] == [2, 1, 3, 4, 7, 11, 18, 29]
    assert arith.lucas(8) == 47 and arith.lucas(16) == 2207


def test_salem():
    assert arith.bolza_quartic() == [1, -4, -2, -4, 1]
    s = arith.salem_number(50)
    with mp.workdps(50):   # reference must be evaluated at 50 digits too (default dps=15 would cap agreement at 1e-16)
        ref = 1 + mp.sqrt(2) + mp.sqrt(2 + 2 * mp.sqrt(2))
    assert mp.almosteq(s, ref, 1e-45)
    assert abs(float(mp.log(s)) - G.systole() / 2) < 1e-14
    roots = mp.polyroots([1, -4, -2, -4, 1], extraprec=60)
    assert sum(1 for r in roots if abs(abs(r) - 1) < 1e-30) == 2
    with mp.workdps(50):
        assert abs(mp.polyval([1, -4, -2, -4, 1], s)) < mp.mpf(10) ** -40
    with mp.workdps(50):
        r50 = mp.polyroots([1, -4, -2, -4, 1], maxsteps=200, extraprec=200)
        assert mp.almosteq(s, max(r.real for r in r50 if abs(r.imag) < mp.mpf(10) ** -30), mp.mpf(10) ** -40)
        assert sum(1 for r in r50 if abs(r) >= s * (1 - mp.mpf(10) ** -40)) == 1
    import sympy
    x = sympy.symbols("x")
    assert sympy.Poly(x**4 - 4*x**3 - 2*x**2 - 4*x + 1).is_irreducible
