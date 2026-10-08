import math
import random
import numpy as np
from bolzakit.exact import zeta8 as Z, group as G

S2 = math.sqrt(2)


# ------------------------------------------------------------------ Z[zeta8]

def test_zeta8_arithmetic():
    z = Z.ZETA
    assert Z.mul(Z.mul(Z.mul(z, z), z), z) == Z.neg(Z.ONE)          # zeta^4 = -1
    assert Z.mul(Z.SQRT2, Z.SQRT2) == (2, 0, 0, 0)
    assert Z.mul(Z.conj(z), z) == Z.ONE
    for _ in range(200):
        x = tuple(random.randint(-9, 9) for _ in range(4))
        y = tuple(random.randint(-9, 9) for _ in range(4))
        assert abs(Z.to_complex(Z.mul(x, y)) - Z.to_complex(x) * Z.to_complex(y)) < 1e-9
        assert abs(Z.to_complex(Z.conj(x)) - Z.to_complex(x).conjugate()) < 1e-12
        p, q = random.randint(-50, 50), random.randint(-50, 50)
        s = Z.zs_sign(p, q)
        v = p + q * S2
        assert s == (0 if v == 0 else (1 if v > 0 else -1))


# ------------------------------------------------------------------ group

def test_relation_exact():
    assert G.word(G.RELATOR) == G.IDENTITY          # exact equality in Z[zeta8]^2, not 1e-13
    # the PPU "abABcdCD" adjacent-pairing relation must FAIL for these generators
    assert G.word([0, 1, 4, 5, 2, 3, 6, 7]) != G.IDENTITY


def test_generators():
    for k, g in enumerate(G.GENS):
        assert G.det(g) == (1, 0)
        assert G.half_trace(g) == (1, 1)               # 1 + sqrt2
        assert G.mul(g, G.GENS[G.inv_index(k)]) == G.IDENTITY
        assert abs(G.disp0(g) - G.L_SYS) < 1e-12
        assert abs(G.act(g, 0) - math.tanh(G.L_SYS / 2) * np.exp(1j * k * math.pi / 4)) < 1e-12
    assert abs(G.L_SYS - 3.0571418390) < 1e-9



def test_mult_matrices():
    rnd = random.Random(1)
    for _ in range(100):
        g = G.word([rnd.randrange(8) for _ in range(5)])
        h = G.word([rnd.randrange(8) for _ in range(4)])
        assert G.elt(G.row(g) @ G.right_mult_matrix(h)) == G.mul(g, h)
        assert G.elt(G.row(g) @ G.left_mult_matrix(h)) == G.mul(h, g)
        assert np.allclose(G.to_matrix(G.mul(g, h)), G.to_matrix(g) @ G.to_matrix(h))


def test_spheres_and_torsion_free():
    S = G.spheres_by_word_length(6)
    sizes = [len(s) for s in S]
    assert sizes == [1, 8, 56, 392, 2736, 19096, 133288]
    # 2736 = 8*7^3 - 8: exactly the 8 coincidences of length-4 words coming from the length-8 relator
    X = np.concatenate(S[1:])
    assert len(G._unique_rows(np.concatenate(S))) == sum(sizes)       # spheres disjoint (parity well defined)
    assert len(G.rows_classify_nonhyperbolic(X)) == 0                  # no elliptic/parabolic, no -I, in |w| <= 6
    assert not np.any(np.all(X == G.row(G.MINUS_IDENTITY), axis=1))


def test_halftrace_integrality_certificate():
    """Re A in Z[sqrt2] for ALL of Gamma: the residues (A,B) mod 2 reachable from the identity under right
    multiplication by generators form a finite closed set on which a1 = a3 (mod 2).  Exact proof by closure."""
    Rm = [R % 2 for R in G.R_GENS]
    start = tuple(G.row(G.IDENTITY) % 2)
    seen = {start}
    front = [start]
    while front:
        nf = []
        for s in front:
            for R in Rm:
                t = tuple((np.array(s) @ R) % 2)
                if t not in seen:
                    seen.add(t)
                    nf.append(t)
        front = nf
    assert all((s[1] - s[3]) % 2 == 0 for s in seen)
    assert len(seen) < 256


def test_ball_completeness_vs_word_enumeration():
    """Every element of word length <= 4 with displacement <= 5 is found by ball_by_displacement(5)."""
    B = G.ball_by_displacement(5.0)
    keys = {tuple(r) for r in B.tolist()}
    S = G.spheres_by_word_length(4)
    for s in S:
        for r in s:
            if G.rows_cosh_disp0(r[None])[0] <= math.cosh(5.0):
                assert tuple(r.tolist()) in keys
    # orbit counting vs hyperbolic area: #{d <= R} ~ (cosh R - 1)/2 * (2pi/A) * 2 = (cosh R - 1)/2
    B10 = G.ball_by_displacement(10.0)
    assert abs(len(B10) / ((math.cosh(10) - 1) / 2) - 1) < 0.05


def test_halftraces_small():
    S = G.spheres_by_word_length(5)
    X = np.concatenate(S[1:])
    P, Q = G.rows_trace(X)
    hts = sorted({(abs(int(p)) // 2, abs(int(q)) // 2) if p + q * S2 > 0 else (-int(p) // 2, -int(q) // 2)
                  for p, q in zip(P, Q)}, key=lambda t: t[0] + t[1] * S2)
    first = [(1, 1), (3, 2), (5, 3), (5, 4), (7, 5)]
    assert hts[:5] == first


# ------------------------------------------------------------------ length spectrum

def test_systole_value_and_label():
    from bolzakit.provenance import info
    assert abs(G.systole() - 2 * math.acosh(1 + math.sqrt(2))) < 1e-15
    assert abs(G.systole() - 3.0571418389619960) < 1e-13
    assert info(G.systole).label == "EXACT"


def test_relator_is_identity():
    assert G.word(G.RELATOR) in (G.IDENTITY, G.MINUS_IDENTITY)


def test_every_generator_det_one():
    assert all(G.det(g) == (1, 0) for g in G.GENS)


def test_elements_counts():
    assert [len(G.elements(L)) for L in (0, 1, 2)] == [1, 9, 65]       # 1, 1+8, 1+8+56
    els = G.elements(2)
    assert len(set(els)) == 65 and G.IDENTITY in els
