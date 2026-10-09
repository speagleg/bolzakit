#!/usr/bin/env python
"""Reproduce the certified computations behind "the Bolza surface is a strict local maximum of lambda_1".

    reproduce.py --quick            check certificates/*.json with bolzakit_check, then re-run the n <= 24 analogues
    reproduce.py --full [--out DIR] regenerate the golden certificates (about 35-50 min on one core) into DIR
                                    (default: certificates/ next to this file); --only NAME... restricts the set

Both modes print a markdown table: claim | certified value | published value | provenance | match.
- certified value: the certificate-attested bound.
- published value: the reference value; provenance says whether it is published in an item-23 file or is a bolzakit
  regression pin.
- match: the value agrees with the reference to the stated tolerance (API rows compare the in-process value; stored
  rows compare the attested bound one-sidedly), AND the certificate's claims hold (for the one non-certifying n = 24
  analogue: that they do not hold).
The exit code is 0 iff every row matches (and, in --quick, no certificate is INVALID).
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

# Reference values.  name -> (claim, value, tolerance, relative, provenance).  provenance says where the value comes
# from: an item-23 file (math-pursuits items/23-bolza-lambda1-extremality/work/, SRC commit 53bc014), or a bolzakit
# regression pin (this implementation, fixed ARPACK start vector) with the ruling that sanctions it.
# The table shows certificate-attested values; rr, first-order, Ky Fan and LG rows from the API are matched on the
# in-process values (SpectrumResult.upper_inprocess, ClusterResult.*_upper_inprocess,
# LowerBoundResult.lam1_lower_inprocess), which is what item 23 published.
I23 = "item 23 "
PUBLISHED = {
    "lp_A_3.3": ("LP count in (0, 3.3] < 1 (mode A)", 0.912612, 1e-6, False, I23 + "lp_A_3.3.json"),
    "lp_B_3.3_4.0": ("LP count in [3.3, 4.0] < 4 (mode B, given lambda_1 >= 3.3)", 3.869661, 1e-6, False,
                     I23 + "lp_B_3.3_4.0.json"),
    "rr_p1_n192": ("lambda_3 <= (P1 Rayleigh-Ritz, n=192)", 3.839184804702046, 1e-9, False,
                   I23 + "cert_rr_bolza_48_96_192.json"),
    "rr_p2_n96": ("lambda_3 <= (P2 Rayleigh-Ritz, n=96)", 3.8388872829261946, 1e-12, False,
                  I23 + "cert_rr_p2_bolza_24_48_96.json"),
    "first_order_FA_n192_su13": ("first order F_A, n=192: b(e,e)/m < 0 (SU13 lambda_1 ball)",
                                 -3.0875694043934345, 1e-6, False, I23 + "cert_first_order_96_192.json"),
    "first_order_FA_n192_lg": ("first order F_A, n=192: b(e,e)/m < 0 (Lehmann-Goerisch lambda_1 ball)",
                               -3.0383476242423058, 1e-6, False, I23 + "lb_cert_first_order_96_192.json"),
    "kyfan_FB_n96_su13": ("Ky Fan F_B, P2 n=96: D2 Lambda(v,v) < 0 (SU13 lambda_1 ball)",
                          -12.850476214196528, 1e-6, False, I23 + "cert_sign_48_96.json"),
    # item 23's lb_cert_sign_96.json publishes -9.341154560 from a lambda ball starting at 3.8388757010263 (above
    # the decimal 3.8388757010 used here, below the certified LG bound); for exactly [3.8388757010, 3.8388872830]
    # the value is this pin.
    "kyfan_FB_n96_lg": ("Ky Fan F_B, P2 n=96: D2 Lambda(v,v) < 0 (Lehmann-Goerisch lambda_1 ball)",
                        -9.341148662380872, 1e-6, False, "bolzakit pin (v0.2 ruling, Task 9 LB golden)"),
    "rank": ("rank d(pi) = 6: det of the chosen 6x6 minor != 0", 7872428136359854080.0, 1e-6, True,
             I23 + "cert_dpi_rank.json"),
    # item 23 publishes float(lb.lower()) rounded to nearest (review_lb_1 F2: ~2.5e-16 above the Arb endpoint); the
    # in-process value here is rounded down, and the stored row is the attested (rounded-down) enclosure end.
    "lg_lower_n48": ("lambda_1 >= (Lehmann-Goerisch, P2 + RT2, n=48, rho=5, given lambda_4 >= 4.0)",
                     3.8388757010527246, 1e-12, False, I23 + "lb_bolza_48.json"),
}
# n <= 24 analogues for --quick.  The n=24 first-order test does NOT certify (b(e,e)/m > 0).  Its value is a bolzakit
# regression pin: no item-23 file has it (cert_first_order_24_48.json predates the energy-weighted transfer bound).
QUICK = {
    "rr_p1_n24": ("lambda_3 <= (P1, n=24)", 3.857985420933427, 1e-5, False, I23 + "cert_rr_bolza_12_24.json"),
    "rr_p2_n24": ("lambda_3 <= (P2, n=24)", 3.838893606899935, 1e-9, False, I23 + "cert_rr_p2_bolza_24_48_96.json"),
    "first_order_FA_n24_su13": ("first order F_A, n=24 (does not certify)", 21.700387418284663, 1e-6, False,
                                "bolzakit pin (v0.2 Task 8 regression pin)"),
    "kyfan_FB_n24_su13": ("Ky Fan F_B, P2 n=24: D2 Lambda(v,v) < 0 (SU13 ball)", -9.713754495605823, 1e-6, False,
                          I23 + "cert_sign_24_48.json"),
}
EXPECT_UNCERTIFIED = {"first_order_FA_n24_su13"}
# A stored certificate attests an upper bound rounded up from the in-process value: mu for rr (to a stored ball), the
# 7-digit outward-rounded radius of the stored enclosure for first-order / Ky Fan.  Stored rows are therefore matched
# one-sidedly: reference - tol <= stored <= reference + tol + slack.
STORED_SLACK = {"rr_upper": 1e-9, "cluster_first_order": 1e-6, "kyfan_second_order": 1e-6}
# A stored lower bound (lg_lower) attests the enclosure's lower end rounded down: matched one-sidedly the other way,
# reference - tol - slack <= stored <= reference + tol.
STORED_SLACK_LOWER = {"lg_lower": 1e-15}


def _inputs():
    from bolzakit.certified import GapInput, Lam1Input
    # SU13 as item 23's scripts used it: arb('3.8388872588', '1.0000001e-6') (the published SU13 numbers use this).
    # dep_id follows the source (PROOF.md §Dependencies): D1 = SU13; D15 = the Lehmann-Goerisch theorem, whose
    # instance here is the lg_lower_n48 certificate.
    su13 = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001",
                     source="SU13 as used in item-23 scripts: arb('3.8388872588','1.0000001e-6')", dep_id="D1")
    lg = Lam1Input("3.8388757010", "3.8388872830",
                   "lower end: Lehmann–Goerisch bound, certificate lg_lower_n48.json (D15; math-pursuits item 23 "
                   "LB_RESULTS); upper end: P2 Rayleigh-Ritz bound, certificate rr_p2_n96.json (min-max, D9)",
                   dep_id="D15")
    gap = GapInput(lower="4.0", source="trace-formula LP count (certified)")
    return {"su13": su13, "lg": lg}, gap


def compute(name):
    """(attested, pinned, certificate) for a golden or quick-analogue name, through the bolzakit API.

    attested: the certificate-attested bound (shown in the table); pinned: the value compared with the published
    one (for rr / first-order / Ky Fan the in-process value item 23 published; otherwise the same as attested)."""
    from bolzakit.certified import Surface, certify_cluster, certify_rank, certify_spectrum, count_eigenvalues
    lam, gap = _inputs()
    if name.startswith("lp_"):
        window = (0, 3.3) if name == "lp_A_3.3" else (3.3, 4.0)
        r = count_eigenvalues("bolza", window=window)
        return r.bound, r.bound, r.certificate
    if name.startswith("rr_"):
        degree, n = int(name[4]), int(name.split("_n")[1])
        r = certify_spectrum(Surface.bolza(), k=3, degree=degree, n=n)
        return r.upper[-1], r.upper_inprocess[-1], r.certificate
    if name.startswith("first_order_FA_"):
        n, route = name.split("_n")[1].split("_")
        r = certify_cluster("F_A", degree=1, n=int(n), lam1=lam[route], gap=gap)
        return r.b_ee_upper, r.b_ee_upper_inprocess, r.certificates[0]
    if name.startswith("kyfan_FB_"):
        n, route = name.split("_n")[1].split("_")
        r = certify_cluster("F_B", degree=2, n=int(n), lam1=lam[route], gap=gap, second_order=True)
        return r.kyfan_upper, r.kyfan_upper_inprocess, r.certificates[0]
    if name.startswith("lg_lower_n"):
        from bolzakit.certified.lower_bounds import certify_lam1_lower
        r = certify_lam1_lower(n=int(name.split("_n")[1]), rho=5, gap=gap)
        return r.lam1_lower, r.lam1_lower_inprocess, r.certificate
    if name == "rank":
        r = certify_rank()
        det = float(_arb(r.det).mid())
        return det, det, r.certificate
    raise KeyError(name)


def _arb(s):
    from flint import arb
    return arb(s)


def stored_value(cert: dict) -> float:
    """The certified value a stored certificate attests (read from its data, not recomputed): upper bounds are rounded
    up from the stored rational / ball (exact arithmetic), lower bounds rounded down."""
    from bolzakit.certified._arb import attested_upper
    k, d = cert["kind"], cert["data"]
    if k == "count":
        return attested_upper(d["bound"])
    if k == "rr_upper":
        return attested_upper(d["mu"])
    if k == "cluster_first_order":
        return attested_upper(d["b_ee_upper_over_m"])
    if k == "kyfan_second_order":
        return attested_upper(d["g2_enclosure"])
    if k == "lg_lower":
        from bolzakit.certified.lower_bounds import _attested_lower
        return _attested_lower(d["lam1_lower_enclosure"])
    if k == "rank":
        from flint import arb_mat
        return float(arb_mat([[_arb(x) for x in row] for row in d["J"]]).det().mid())
    raise KeyError(k)


def _matches(value, published, tol, relative):
    if value is None:
        return False
    err = abs(value - published)
    return err <= tol * abs(published) if relative else err <= tol


def _row(claim, value, published, provenance, ok):
    return f"| {claim} | {value!r} | {published!r} | {provenance} | {'match' if ok else 'MISMATCH'} |"


def _header():
    return ["| claim | certified value | published value | provenance | match |", "|---|---|---|---|---|"]


def full(out: pathlib.Path, only=None) -> bool:
    out.mkdir(parents=True, exist_ok=True)
    rows, all_ok = _header(), True
    for name, (claim, pub, tol, rel, src) in PUBLISHED.items():
        if only and name not in only:
            continue
        t0 = time.time()
        value, pinned, cert = compute(name)
        cert.save(out / f"{name}.json")
        ok = cert.certified and _matches(pinned, pub, tol, rel)
        all_ok &= ok
        print(f"[{name}] {time.time() - t0:.0f} s, certified={cert.certified}, attested={value!r}, "
              f"pinned={pinned!r}", flush=True)
        rows.append(_row(claim, value, pub, src, ok))
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
    for name, (claim, pub, tol, rel, src) in PUBLISHED.items():
        p = CERTS / f"{name}.json"
        if not p.exists():
            continue
        cert = json.loads(p.read_text(encoding="utf-8"))
        value = stored_value(cert)
        if cert["kind"] in STORED_SLACK:    # attested upper bound: at or above the reference, within the slack
            ok = pub - tol <= value <= pub + tol + STORED_SLACK[cert["kind"]]
        elif cert["kind"] in STORED_SLACK_LOWER:    # attested lower bound: at or below the reference
            ok = value is not None and pub - tol - STORED_SLACK_LOWER[cert["kind"]] <= value <= pub + tol
        else:
            ok = _matches(value, pub, tol, rel)
        ok = ok and cert["certified"]
        all_ok &= ok
        rows.append(_row(claim + " [stored]", value, pub, src, ok))
    for name, (claim, pub, tol, rel, src) in QUICK.items():
        value, pinned, cert = compute(name)
        ok = (cert.certified != (name in EXPECT_UNCERTIFIED)) and _matches(pinned, pub, tol, rel)
        all_ok &= ok
        rows.append(_row(claim, value, pub, src, ok))
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
    if args.only and not args.full:
        ap.error("--only applies to --full only")
    t0 = time.time()
    ok = full(args.out, set(args.only or ())) if args.full else quick()
    print(f"\n{'ALL MATCH' if ok else 'SOME ROWS DO NOT MATCH'} ({time.time() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
