"""Independent checker for bolzakit-cert/1 certificates (python-flint, mpmath and the stdlib only).

check_file(path, deep=False) -> (ok, message).  A certificate is VALID iff
  * it is strict JSON with the bolzakit-cert/1 schema and its sha256 matches the canonical core;
  * its claims have exactly the expected (lhs, op, rhs) for its kind, and every stored `holds` equals the truth value
    re-derived here from the certificate's inputs and data (a claim that cannot be re-derived is INVALID);
  * every claim holds (and the stored `certified` flag says so).
Messages: "VALID", "VALID given stored Gram data" (the claim is re-derived, but from Gram / element-sum data that
the checker cannot recompute), "VALID (conditional on lambda_1 >= a)" (mode-B count), or "INVALID: <reason>".

Trust boundary per kind (what is re-derived / what is trusted):
  rr_upper            re-derived: M and mu M - K positive definite (own LDL^T, Sylvester fallback); degree 2:
                      tr(M^-1 K) <= the cluster-sum rhs.  Trusted: the stored Gram balls K, M (element integration;
                      --deep re-integrates a 1% element sample, consistency only).
  count               re-derived from c_coeffs, s, the window and closed-form lengths only: exact Sturm sign conditions,
                      Arb (I + G_P - F(0)) / min F vs target.  Trusted: FBP trace formula and geodesic data
                      (dependencies D2, D3); mode B is conditional on its hypothesis lambda_1 >= a (H_a).
  rank                re-derived from the 16 Klein balls and the words only: own dual-number Jacobians of the
                      constraints and traces, Tb, the residual Jc Tb, det of the 6 x 6 minor and of the stacked minor.
                      Stored J / dets / residual are cross-checked only.  Trusted: nothing beyond Arb.
  cluster_first_order re-derived from data["raw"] (element sums for the trial vector) and the lam1 / gap inputs:
                      mu, B, E and B + E < 0.  Trusted: the raw sums (Gram data) and the lam1 / gap dependencies.
  kyfan_second_order  re-derived: g2 rebuilt from data["raw"] < 0.  Trusted: the raw Gram matrices / weighted sums and
                      the lam1 / gap dependencies.
"""
from __future__ import annotations

import random

from flint import arb

from . import assembly, deep as deep_mod, hermitian, lp, rank
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


def _rank(raw):
    with precision(128):
        return rank.derive(raw["inputs"], raw["data"], ball), None


def _first_order(raw):
    s = assembly.b_ee(raw["data"], raw["inputs"])
    return [("b_ee_upper_over_m", "<", "0", lt(s, 0))], None


def _kyfan(raw):
    g = assembly.g2(raw["data"], raw["inputs"])
    return [("g2_upper", "<", "0", lt(g, 0))], None


GRAM_KINDS = ("rr_upper", "cluster_first_order", "kyfan_second_order")
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
    if note is not None:
        msg = f"VALID (conditional on {note})"
    elif kind in GRAM_KINDS:
        msg = "VALID given stored Gram data"
    else:
        msg = "VALID"
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
    except Exception as e:  # noqa: BLE001  malformed input never aborts a multi-file run
        return False, f"INVALID: malformed certificate ({type(e).__name__}: {e})"
