"""Independent checker for bolzakit-cert/1 certificates (python-flint, mpmath and the stdlib only).

check_file(path, deep=False) -> (ok, message).  A certificate is VALID iff
  * it is strict JSON with the bolzakit-cert/1 schema and its sha256 matches the canonical core;
  * its claims have exactly the expected (lhs, op, rhs) for its kind, and every stored `holds` equals the truth value
    re-derived here from the certificate's inputs and data (a claim that cannot be re-derived is INVALID);
  * every claim holds (and the stored `certified` flag says so).
What is re-derived per kind:
  rr_upper            M and mu M - K positive definite (own LDL^T); degree 2: tr(M^-1 K) <= the cluster-sum rhs.
  count               lp.verify: exact Sturm sign conditions, Arb re-evaluation of (I + G_P - F(0)) / min F vs target.
                      Mode B is VALID only conditionally on its hypothesis lambda_1 >= a.
  rank                own Bareiss determinant of the 6 x 6 minor J excludes 0; the 5 constraints recomputed from the
                      16 Klein balls enclose 0; the stored stacked determinant / kernel residual enclosures are finite
                      (and the determinant excludes 0).
  cluster_first_order assembly.b_ee: B + E recomputed from the transfer terms and the lam1 / gap inputs, < 0.
  kyfan_second_order  assembly.g2: g2 rebuilt from data["raw"], < 0.
"""
from __future__ import annotations

import random

from flint import arb

from . import assembly, deep as deep_mod, hermitian, lp
from ._cert import SCHEMA, CheckError, arb_q, ball, core_hash, frac, load, lt, matrix, precision


# ---------------------------------------------------------------- per-kind re-derivations: [(lhs, op, rhs, truth)]
def _rr(raw):
    inputs, data = raw["inputs"], raw["data"]
    k, degree = inputs["k"], inputs["degree"]
    if not isinstance(k, int) or k < 1 or degree not in (1, 2):
        raise CheckError("bad rr_upper inputs k / degree")
    mu_s = data["mu"]
    with precision(128):
        K = hermitian.check_symmetric(matrix(data["K"], k + 1))
        M = hermitian.check_symmetric(matrix(data["M"], k + 1))
        m_pd = hermitian.is_pd(M)
        if mu_s == "inf":
            ok_mu = False
        else:
            mu = arb_q(frac(mu_s))
            ok_mu = m_pd and hermitian.is_pd([[mu * M[i][j] - K[i][j] for j in range(k + 1)]
                                               for i in range(k + 1)])
        claims = [(f"lambda_{k}", "<=", mu_s, ok_mu)]
        if degree == 2:
            # Ky Fan on the (k+1)-dim trial space (lambda_0 = 0): lambda_1 + .. + lambda_k <= tr(M^-1 K)
            rhs = raw["claims"][1]["rhs"] if len(raw["claims"]) > 1 else None
            ok_cs = False
            if isinstance(rhs, str) and rhs != "inf":
                tr = hermitian.trace_minv_k(K, M)
                ok_cs = m_pd and lt(tr, arb_q(frac(rhs)))
            claims.append(("+".join(f"lambda_{i}" for i in range(1, k + 1)), "<=", rhs, ok_cs))
    return claims, None


def _count(raw):
    return lp.verify(raw["data"], raw["inputs"], raw["dependencies"])


def _constraints(klein):
    """The 5 octagon constraints (opposite sides equal, angle sum 2 pi) at Klein balls (16 values)."""
    B = lambda x, y: x[0] * y[0] - x[1] * y[1] - x[2] * y[2]  # noqa: E731
    H = []
    for k in range(8):
        x, y = klein[2 * k], klein[2 * k + 1]
        t = 1 / (1 - x * x - y * y).sqrt()
        H.append([t, t * x, t * y])
    side = [B(H[k], H[(k + 1) % 8]).acosh() for k in range(8)]
    total = arb(0)
    for k in range(8):
        a, p, q = H[k], H[(k - 1) % 8], H[(k + 1) % 8]
        t1 = [p[i] - B(a, p) * a[i] for i in range(3)]      # tangents at a toward p and q (spacelike: B < 0)
        t2 = [q[i] - B(a, q) * a[i] for i in range(3)]
        total += (-B(t1, t2) / (B(t1, t1) * B(t2, t2)).sqrt()).acos()   # cos = <t1,t2> / |t1||t2|, <,> = -B
    return [side[k + 4] - side[k] for k in range(4)] + [total - 2 * arb.pi()]


def _rank(raw):
    inputs, data = raw["inputs"], raw["data"]
    fin0 = lambda v: v.is_finite() and v.contains(0)  # noqa: E731
    with precision(128):
        cols = data["cols"]
        ok_det = False
        if isinstance(cols, list) and len(cols) == 6 and len(set(cols)) == 6 and all(
                isinstance(c, int) and 0 <= c < 11 for c in cols):
            d = hermitian.bareiss_det(matrix(data["J"], 6))
            ok_det = bool(d.is_finite() and not d.contains(0))
        sc, sd = data["stacked_cols"], data["stacked_det"]
        ok_st = (isinstance(sc, list) and len(sc) == 11 and len(set(sc)) == 11
                 and all(isinstance(c, int) and 0 <= c < 16 for c in sc)
                 and isinstance(sd, str) and ball(sd).is_finite() and not ball(sd).contains(0))
        kr = data["kernel_residual_upper"]
        ok_ker = isinstance(kr, str) and ball(kr).is_finite()
        klein = [ball(x) for x in inputs["klein"]]
        if len(klein) != 16:
            raise CheckError("inputs.klein must hold 16 balls")
        cons = _constraints(klein)
        stored = data["constraints"]
        ok_con = (all(fin0(c) for c in cons) and isinstance(stored, list) and len(stored) == 5
                  and all(fin0(ball(x)) for x in stored))
    return [("det J[:, cols]", "!=0", "0", ok_det), ("det [Jc; Jf][:, stacked_cols]", "!=0", "0", ok_st),
            ("Jc . Tb (all entries)", "encloses", "0", ok_ker), ("constraints(x)", "encloses", "0", ok_con)], None


def _first_order(raw):
    s = assembly.b_ee(raw["data"], raw["inputs"])
    return [("b_ee_upper_over_m", "<", "0", lt(s, 0))], None


def _kyfan(raw):
    g = assembly.g2(raw["data"], raw["inputs"])
    return [("g2_upper", "<", "0", lt(g, 0))], None


KINDS = {"rr_upper": _rr, "count": _count, "rank": _rank, "cluster_first_order": _first_order,
         "kyfan_second_order": _kyfan}


def _check_sample_selection(raw):
    """The element_sample indices are the deterministic selection of the producer's documented rule."""
    es = raw["data"].get("element_sample")
    if not es:
        return
    n_el, items = es["n_elements"], es["items"]
    core = dict(raw, data={k: v for k, v in raw["data"].items() if k != "element_sample"})
    seed = int(core_hash(core)[:8], 16)
    want = sorted(random.Random(seed).sample(range(n_el), max(1, min(200, n_el // 100))))
    if [it["index"] for it in items] != want:
        raise CheckError("element_sample indices are not the documented deterministic selection")


def _check(raw, deep):
    for key in ("schema", "kind", "inputs", "claims", "data", "dependencies", "sha256", "certified"):
        if key not in raw:
            raise CheckError(f"missing field {key!r}")
    if raw["schema"] != SCHEMA:
        raise CheckError(f"unknown schema {raw['schema']!r}")
    if core_hash(raw) != raw["sha256"]:
        return False, "INVALID: sha256 mismatch"
    kind = raw["kind"]
    if kind not in KINDS:
        raise CheckError(f"unknown kind {kind!r}")
    stored = raw["claims"]
    if not isinstance(stored, list) or not stored:
        raise CheckError("no claims")
    try:
        derived, note = KINDS[kind](raw)
    except (KeyError, IndexError, TypeError, ValueError, ZeroDivisionError) as e:
        raise CheckError(f"cannot re-derive the claims ({type(e).__name__}: {e})") from e
    if len(derived) != len(stored):
        raise CheckError(f"expected {len(derived)} claims for kind {kind}, found {len(stored)}")
    for (lhs, op, rhs, truth), c in zip(derived, stored):
        if (c.get("lhs"), c.get("op"), c.get("rhs")) != (lhs, op, rhs):
            raise CheckError(f"claim {c.get('lhs')!r} {c.get('op')} {c.get('rhs')!r} is not the expected "
                             f"{lhs!r} {op} {rhs!r}")
        if c.get("holds") is not truth:
            raise CheckError(f"claim {lhs} {op} {rhs}: stored holds={c.get('holds')!r}, re-derived {truth}")
    if raw["certified"] is not all(t for *_, t in derived):
        raise CheckError("stored 'certified' flag disagrees with the claims")
    failed = [f"{lhs} {op} {rhs}" for lhs, op, rhs, t in derived if not t]
    if failed:
        return False, "INVALID: not certified (re-derivation agrees): " + "; ".join(failed) + " does not hold"
    msg = "VALID" if note is None else f"VALID (conditional on {note}); assumes: {note}"
    if deep:
        if kind == "rr_upper":
            _check_sample_selection(raw)
            ok, dmsg = deep_mod.check_elements(raw["data"], raw["inputs"]["degree"])
            if not ok:
                return False, "INVALID: " + dmsg
            msg += " " + dmsg
        else:
            msg += " (deep: consistency only; nothing to re-integrate for this kind)"
    return True, msg


def check_file(path, deep=False) -> tuple:
    """(ok, message): message is 'VALID...' or 'INVALID: <reason>'."""
    try:
        return _check(load(path), deep)
    except CheckError as e:
        return False, f"INVALID: {e}"
    except (OSError, ValueError) as e:
        return False, f"INVALID: unreadable certificate ({type(e).__name__}: {e})"
