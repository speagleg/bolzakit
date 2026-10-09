"""python -m bolzakit_check FILE... [--deep]: prints 'FILE: VALID' / 'FILE: INVALID: <reason>'; exit 0 iff all VALID."""
import argparse
import sys

from .check import check_file


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m bolzakit_check", description=__doc__)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--deep", action="store_true", help="also re-integrate sampled elements (consistency only)")
    args = ap.parse_args(argv)
    all_ok = True
    for f in args.files:
        ok, msg = check_file(f, deep=args.deep)
        print(f"{f}: {msg}")
        all_ok = all_ok and ok
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
