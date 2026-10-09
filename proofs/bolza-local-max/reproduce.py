#!/usr/bin/env python
"""Reproduce the certified computations behind "the Bolza surface is a strict local maximum of lambda_1".

    reproduce.py --quick            check certificates/*.json with bolzakit_check, then re-run the n <= 24 analogues
    reproduce.py --full [--out DIR] regenerate the golden certificates (about 30-45 min on one core) into DIR
                                    (default: certificates/ next to this file); --only NAME... restricts the set

Both modes print a markdown table: claim | certified value | published value | match.  "match" means the value
agrees with the published one to the stated tolerance AND the certificate's claims hold (or, for the one
non-certifying n = 24 analogue, that they do not hold, as published).  The exit code is 0 iff every row matches
(and, in --quick, no certificate is INVALID).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
CERTS = HERE / "certificates"

# Published values (spec §7.1; item-23 certificates and LB_RESULTS.md; the two Lehmann-Goerisch rows pin bolzakit's
# own values, which use the outward decimal rounding 3.8388757010 of the LG lower bound).  name -> (claim, value, tolerance, relative).
# rr rows compare the in-process bisection value (SpectrumResult.upper_inprocess), which is what item 23 published.
PUBLISHED = {
    "lp_A_3.3": ("LP count in (0, 3.3] < 1 (mode A)", 0.912612, 1e-6, False),
    "lp_B_3.3_4.0": ("LP count in [3.3, 4.0] < 4 (mode B, given lambda_1 >= 3.3)", 3.869661, 1e-6, False),
    "rr_p1_n192": ("lambda_3 <= (P1 Rayleigh-Ritz, n=192)", 3.839184804702046, 1e-9, False),
    "rr_p2_n96": ("lambda_3 <= (P2 Rayleigh-Ritz, n=96)", 3.8388872829261946, 1e-12, False),
    "first_order_FA_n192_su13": ("first order F_A, n=192: b(e,e)/m < 0 (SU13 lambda_1 ball)",
                                 -3.0875694043934345, 1e-6, False),
    "first_order_FA_n192_lg": ("first order F_A, n=192: b(e,e)/m < 0 (Lehmann-Goerisch lambda_1 ball)",
                               -3.0383475150201082, 1e-6, False),   # pinned: this API (item 23: -3.0383476)
    "kyfan_FB_n96_su13": ("Ky Fan F_B, P2 n=96: D2 Lambda(v,v) < 0 (SU13 lambda_1 ball)",
                          -12.850476214196528, 1e-6, False),
    "kyfan_FB_n96_lg": ("Ky Fan F_B, P2 n=96: D2 Lambda(v,v) < 0 (Lehmann-Goerisch lambda_1 ball)",
                        -9.341148662380872, 1e-6, False),
    "rank": ("rank d(pi) = 6: det of the chosen 6x6 minor != 0", 7872428136359854080.0, 1e-6, True),
}
# n <= 24 analogues for --quick (published: item-23 cert_rr_bolza_12_24, cert_rr_p2_bolza_24_48, cert_first_order_24_48,
# cert_sign_24_48).  The n=24 first-order test does NOT certify (b(e,e)/m > 0); that is the published outcome.
QUICK = {
    "rr_p1_n24": ("lambda_3 <= (P1, n=24)", 3.857985420933427, 1e-5, False),
    "rr_p2_n24": ("lambda_3 <= (P2, n=24)", 3.838893606899935, 1e-9, False),
    "first_order_FA_n24_su13": ("first order F_A, n=24 (does not certify, as published)", 21.700387418284663, 1e-6,
                                False),
    "kyfan_FB_n24_su13": ("Ky Fan F_B, P2 n=24: D2 Lambda(v,v) < 0 (SU13 ball)", -9.713754495605823, 1e-6, False),
}
EXPECT_UNCERTIFIED = {"first_order_FA_n24_su13"}
# The stored rr certificate attests mu = the bisection value rounded up to a stored ball (SpectrumResult.upper);
# it sits at most this far above the in-process value that item 23 published.
STORED_RR_SLACK = 1e-9


def _inputs():
    from bolzakit.certified import GapInput, Lam1Input
    # SU13 as item 23's scripts used it: arb('3.8388872588', '1.0000001e-6') (the published SU13 numbers use this).
    su13 = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001",
                     source="SU13 as used in item-23 scripts: arb('3.8388872588','1.0000001e-6')")
    lg = Lam1Input("3.8388757010", "3.8388872830",
                   "Lehmann–Goerisch lower bound (math-pursuits item 23 LB_RESULTS) + P2 upper bound")
    gap = GapInput(lower="4.0", source="trace-formula LP count (certified)")
    return {"su13": su13, "lg": lg}, gap


def compute(name):
    """(value, certificate) for a golden or quick-analogue name, through the bolzakit API."""
    from bolzakit.certified import Surface, certify_cluster, certify_rank, certify_spectrum, count_eigenvalues
    lam, gap = _inputs()
    if name.startswith("lp_"):
        window = (0, 3.3) if name == "lp_A_3.3" else (3.3, 4.0)
        r = count_eigenvalues("bolza", window=window)
        return r.bound, r.certificate
    if name.startswith("rr_"):
        degree, n = int(name[4]), int(name.split("_n")[1])
        r = certify_spectrum(Surface.bolza(), k=3, degree=degree, n=n)
        return r.upper_inprocess[-1], r.certificate
    if name.startswith("first_order_FA_"):
        n, route = name.split("_n")[1].split("_")
        r = certify_cluster("F_A", degree=1, n=int(n), lam1=lam[route], gap=gap)
        return r.b_ee_upper, r.certificates[0]
    if name.startswith("kyfan_FB_"):
        n, route = name.split("_n")[1].split("_")
        r = certify_cluster("F_B", degree=2, n=int(n), lam1=lam[route], gap=gap, second_order=True)
        return r.kyfan_upper, r.certificates[0]
    if name == "rank":
        r = certify_rank()
        return float(_arb(r.det).mid()), r.certificate
    raise KeyError(name)


def _arb(s):
    from flint import arb
    return arb(s)


def stored_value(cert: dict) -> float:
    """The certified value a stored certificate attests (read from its data, not recomputed)."""
    k, d = cert["kind"], cert["data"]
    if k == "count":
        return float(d["bound_upper"])
    if k == "rr_upper":
        return float(_arb(d["mu"]).upper())
    if k == "cluster_first_order":
        return float(_arb(d["b_ee_upper_over_m"]).upper())
    if k == "kyfan_second_order":
        return float(_arb(d["g2_enclosure"]).upper())
    if k == "rank":
        from flint import arb_mat
        return float(arb_mat([[_arb(x) for x in row] for row in d["J"]]).det().mid())
    raise KeyError(k)


def _matches(value, published, tol, relative):
    if value is None:
        return False
    err = abs(value - published)
    return err <= tol * abs(published) if relative else err <= tol


def _row(claim, value, published, ok):
    return f"| {claim} | {value!r} | {published!r} | {'match' if ok else 'MISMATCH'} |"


def _header():
    return ["| claim | certified value | published value | match |", "|---|---|---|---|"]


def full(out: pathlib.Path, only=None) -> bool:
    out.mkdir(parents=True, exist_ok=True)
    rows, all_ok = _header(), True
    for name, (claim, pub, tol, rel) in PUBLISHED.items():
        if only and name not in only:
            continue
        t0 = time.time()
        value, cert = compute(name)
        cert.save(out / f"{name}.json")
        ok = cert.certified and _matches(value, pub, tol, rel)
        all_ok &= ok
        print(f"[{name}] {time.time() - t0:.0f} s, certified={cert.certified}, value={value!r}", flush=True)
        rows.append(_row(claim, value, pub, ok))
    print("\n".join(rows))
    return all_ok


def quick() -> bool:
    files = sorted(CERTS.glob("*.json"))
    missing = [n for n in PUBLISHED if not (CERTS / f"{n}.json").exists()]
    r = subprocess.run([sys.executable, "-m", "bolzakit_check", *map(str, files)], capture_output=True, text=True)
    print(r.stdout.replace(str(HERE) + "/", ""), end="")
    print(f"bolzakit_check exit code {r.returncode} (0 = all VALID; 2 = VALID, some only conditionally; "
          f"1 = some INVALID)\n")
    all_ok = r.returncode in (0, 2) and not missing and bool(files)
    if missing:
        print("missing golden certificates: " + ", ".join(missing))
    rows = _header()
    for name, (claim, pub, tol, rel) in PUBLISHED.items():
        p = CERTS / f"{name}.json"
        if not p.exists():
            continue
        cert = json.loads(p.read_text(encoding="utf-8"))
        value = stored_value(cert)
        if cert["kind"] == "rr_upper":    # attested mu: at or above the published bisection value, within the slack
            ok = pub - tol <= value <= pub + tol + STORED_RR_SLACK
        else:
            ok = _matches(value, pub, tol, rel)
        ok = ok and cert["certified"]
        all_ok &= ok
        rows.append(_row(claim + " [stored]", value, pub, ok))
    for name, (claim, pub, tol, rel) in QUICK.items():
        value, cert = compute(name)
        ok = (cert.certified != (name in EXPECT_UNCERTIFIED)) and _matches(value, pub, tol, rel)
        all_ok &= ok
        rows.append(_row(claim, value, pub, ok))
    print("\n".join(rows))
    return all_ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true", help="check stored certificates + n <= 24 analogues (default)")
    mode.add_argument("--full", action="store_true", help="regenerate the golden certificates")
    ap.add_argument("--out", type=pathlib.Path, default=CERTS, help="output directory for --full")
    ap.add_argument("--only", nargs="+", choices=list(PUBLISHED), help="with --full: regenerate only these")
    args = ap.parse_args(argv)
    t0 = time.time()
    ok = full(args.out, set(args.only or ())) if args.full else quick()
    print(f"\n{'ALL MATCH' if ok else 'SOME ROWS DO NOT MATCH'} ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
